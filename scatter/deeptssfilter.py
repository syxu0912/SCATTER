"""CNN-BiLSTM DeepTSSFilter inference utilities."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import numpy as np
import torch
from torch import nn
from torch.nn import functional as F


_DNA = {"A": (1, 0, 0, 0), "T": (0, 1, 0, 0), "C": (0, 0, 1, 0), "G": (0, 0, 0, 1)}
_COMPLEMENT = str.maketrans("ATCGNatcgn", "TAGCNtagcn")


class DeepTSSFilter(nn.Module):
    """The frozen CNN-BiLSTM architecture used by SCATTER inference."""

    def __init__(self) -> None:
        super().__init__()
        self.conv1 = nn.Conv1d(4, 128, kernel_size=12)
        self.conv2 = nn.Conv1d(128, 64, kernel_size=6)
        self.pool = nn.MaxPool1d(kernel_size=4)
        self.dropout = nn.Dropout(0.3)
        self.lstm = nn.LSTM(64, 128, batch_first=True, bidirectional=True)
        self.fc1 = nn.Linear(256, 64)
        self.fc2 = nn.Linear(64, 2)

    def forward(self, value: torch.Tensor) -> torch.Tensor:
        value = F.relu(self.conv1(value))
        value = F.relu(self.conv2(value))
        value = self.pool(value)
        value = self.dropout(value)
        value = value.transpose(1, 2)
        value, _ = self.lstm(value)
        value = value[:, -1, :]
        value = F.relu(self.fc1(value))
        return self.fc2(value)


class DeepTSSFilterDeepPASS(nn.Module):
    """DeepPASS-style CNN-BiLSTM architecture for a new model version.

    This keeps the sequence input and convolutional/recurrent stack used by
    the legacy model, while matching the paper's 1024-unit dense layer and
    placing dropout after that layer.  It is intentionally a separate class
    so existing v1 weights remain loadable and reproducible.
    """

    def __init__(self) -> None:
        super().__init__()
        self.conv1 = nn.Conv1d(4, 128, kernel_size=12)
        self.conv2 = nn.Conv1d(128, 64, kernel_size=6)
        self.pool = nn.MaxPool1d(kernel_size=4)
        self.lstm = nn.LSTM(64, 128, batch_first=True, bidirectional=True)
        self.fc1 = nn.Linear(256, 1024)
        self.dropout = nn.Dropout(0.3)
        self.fc2 = nn.Linear(1024, 2)

    def forward(self, value: torch.Tensor) -> torch.Tensor:
        value = F.relu(self.conv1(value))
        value = F.relu(self.conv2(value))
        value = self.pool(value)
        value = value.transpose(1, 2)
        value, _ = self.lstm(value)
        value = value[:, -1, :]
        value = F.relu(self.fc1(value))
        value = self.dropout(value)
        return self.fc2(value)


class OptimizedDeepTSSFilter(nn.Module):
    """Configurable single-model CNN-BiLSTM used by architecture v2 packages."""

    def __init__(
        self,
        *,
        conv_mode: str = "sequential",
        readout: str = "hidden",
        dense_units: int = 512,
        dropout: float = 0.3,
        layer_norm: bool = False,
        conv1_channels: int = 128,
        conv2_channels: int = 64,
        kernel1: int = 12,
        kernel2: int = 6,
        pool_size: int = 4,
        lstm_hidden: int = 128,
    ) -> None:
        super().__init__()
        self.conv_mode = conv_mode
        self.readout = readout
        self.dense_units = dense_units
        self.dropout_rate = dropout
        self.use_layer_norm = layer_norm
        self.conv1_channels = conv1_channels
        self.conv2_channels = conv2_channels
        self.kernel1 = kernel1
        self.kernel2 = kernel2
        self.pool_size = pool_size
        self.lstm_hidden = lstm_hidden
        if conv_mode == "sequential":
            self.conv1 = nn.Conv1d(4, conv1_channels, kernel_size=kernel1)
            self.conv2 = nn.Conv1d(conv1_channels, conv2_channels, kernel_size=kernel2)
            lstm_input = conv2_channels
        elif conv_mode == "multiscale":
            self.conv7 = nn.Conv1d(4, 32, kernel_size=7, padding=3)
            self.conv15 = nn.Conv1d(4, 32, kernel_size=15, padding=7)
            self.conv31 = nn.Conv1d(4, 32, kernel_size=31, padding=15)
            self.conv_fuse = nn.Conv1d(96, 64, kernel_size=5, padding=2)
            lstm_input = 64
        else:
            raise ValueError(f"Unknown conv mode: {conv_mode}")
        self.pool = nn.MaxPool1d(kernel_size=pool_size)
        self.lstm = nn.LSTM(lstm_input, lstm_hidden, batch_first=True, bidirectional=True)
        representation_units = 2 * lstm_hidden if readout in {"legacy_last", "hidden"} else 4 * lstm_hidden
        if readout not in {"legacy_last", "hidden", "hidden_center"}:
            raise ValueError(f"Unknown readout: {readout}")
        self.norm = nn.LayerNorm(representation_units) if layer_norm else nn.Identity()
        self.fc1 = nn.Linear(representation_units, dense_units)
        self.dropout = nn.Dropout(dropout)
        self.fc2 = nn.Linear(dense_units, 2)

    def architecture_config(self) -> dict[str, object]:
        return {
            "conv_mode": self.conv_mode,
            "readout": self.readout,
            "dense_units": self.dense_units,
            "dropout": self.dropout_rate,
            "layer_norm": self.use_layer_norm,
            "conv1_channels": self.conv1_channels,
            "conv2_channels": self.conv2_channels,
            "kernel1": self.kernel1,
            "kernel2": self.kernel2,
            "pool_size": self.pool_size,
            "lstm_hidden": self.lstm_hidden,
        }

    def forward(self, value: torch.Tensor) -> torch.Tensor:
        if self.conv_mode == "sequential":
            value = F.relu(self.conv1(value))
            value = F.relu(self.conv2(value))
        else:
            value = torch.cat(
                [F.relu(self.conv7(value)), F.relu(self.conv15(value)), F.relu(self.conv31(value))],
                dim=1,
            )
            value = F.relu(self.conv_fuse(value))
        value = self.pool(value).transpose(1, 2)
        output, (hidden, _cell) = self.lstm(value)
        if self.readout == "legacy_last":
            representation = output[:, -1, :]
        else:
            global_state = torch.cat([hidden[-2], hidden[-1]], dim=1)
            if self.readout == "hidden_center":
                representation = torch.cat([global_state, output[:, output.shape[1] // 2, :]], dim=1)
            else:
                representation = global_state
        representation = self.norm(representation)
        return self.fc2(self.dropout(F.relu(self.fc1(representation))))


class DeepTSSFilterProbabilityEnsemble(nn.Module):
    """Weighted probability ensemble loaded from a small JSON manifest."""

    def __init__(self, models: list[nn.Module], weights: list[float]) -> None:
        super().__init__()
        if not models or len(models) != len(weights):
            raise ValueError("An ensemble requires equally sized, non-empty model and weight lists")
        if any(weight < 0 for weight in weights):
            raise ValueError("Ensemble weights must be non-negative")
        weight_sum = float(sum(weights))
        if not np.isclose(weight_sum, 1.0):
            raise ValueError(f"Ensemble weights must sum to 1; observed {weight_sum}")
        self.models = nn.ModuleList(models)
        self.register_buffer("weights", torch.tensor(weights, dtype=torch.float32))

    def forward(self, value: torch.Tensor) -> torch.Tensor:
        probabilities = torch.stack(
            [torch.softmax(model(value), dim=1) for model in self.models], dim=0
        )
        mean_probability = torch.sum(
            probabilities * self.weights[:, None, None], dim=0
        ).clamp_min(torch.finfo(probabilities.dtype).tiny)
        # score_sequences applies softmax to logits.  Log probabilities are
        # valid logits and recover the weighted probabilities exactly.
        return torch.log(mean_probability)


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def one_hot(sequence: str) -> np.ndarray:
    return np.asarray([_DNA.get(base, (0, 0, 0, 0)) for base in sequence.upper()], dtype=np.float32).T


def load_model(weights: str | Path, device: str = "cpu") -> nn.Module:
    path = Path(weights).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"DeepTSSFilter weights do not exist: {path}")
    if path.suffix.lower() == ".json":
        manifest = json.loads(path.read_text())
        if manifest.get("type") != "deeptssfilter_probability_ensemble":
            raise ValueError(f"Unsupported DeepTSSFilter manifest type in {path}")
        members = manifest.get("members")
        if not isinstance(members, list) or not members:
            raise ValueError(f"DeepTSSFilter ensemble has no members: {path}")
        member_models = []
        member_weights = []
        for member in members:
            member_path = Path(member["model"])
            if not member_path.is_absolute():
                member_path = path.parent / member_path
            member_models.append(load_model(member_path, device=device))
            member_weights.append(float(member["weight"]))
        model = DeepTSSFilterProbabilityEnsemble(member_models, member_weights)
        lengths = {getattr(member, 'sequence_length_bp', 200) for member in member_models}
        if len(lengths) != 1:
            raise ValueError('Ensemble members must have the same input length')
        model.sequence_length_bp = lengths.pop()
        model.to(device)
        model.eval()
        return model
    try:
        state = torch.load(path, map_location=device, weights_only=True)
    except TypeError:  # pragma: no cover - compatibility with older torch
        state = torch.load(path, map_location=device)
    if isinstance(state, dict) and state.get("format_version") == 2 and state.get("architecture") == "optimized_cnn_bilstm":
        model = OptimizedDeepTSSFilter(**state["config"])
        model.load_state_dict(state["state_dict"], strict=True)
        model.sequence_length_bp = int(state.get('sequence_length_bp', 200))
        model.to(device)
        model.eval()
        return model
    dense_units = int(state["fc1.weight"].shape[0])
    if dense_units == 1024:
        model = DeepTSSFilterDeepPASS()
    elif dense_units == 64:
        model = DeepTSSFilter()
    else:
        raise ValueError(f"Unsupported DeepTSSFilter fc1 width in {path}: {dense_units}")
    model.load_state_dict(state, strict=True)
    model.sequence_length_bp = 200
    model.to(device)
    model.eval()
    return model


def score_sequences(model: nn.Module, sequences: Iterable[str], device: str = "cpu", batch_size: int = 256) -> np.ndarray:
    prepared = list(sequences)
    if not prepared:
        return np.empty(0, dtype=np.float32)
    lengths = {len(sequence) for sequence in prepared}
    expected = int(getattr(model, 'sequence_length_bp', 200))
    if lengths != {expected}:
        raise ValueError(f"DeepTSSFilter requires {expected} bp sequences; observed lengths: {sorted(lengths)}")
    if batch_size < 1:
        raise ValueError('batch_size must be positive')
    if any(set(seq.upper()) - set('ATCG') for seq in prepared):
        raise ValueError('Model scoring requires unambiguous A/T/C/G sequences')
    scores: list[np.ndarray] = []
    with torch.no_grad():
        for offset in range(0, len(prepared), batch_size):
            encoded = np.stack([one_hot(sequence) for sequence in prepared[offset : offset + batch_size]])
            logits = model(torch.from_numpy(encoded).to(device))
            scores.append(torch.softmax(logits, dim=1)[:, 1].cpu().numpy())
    return np.concatenate(scores)
