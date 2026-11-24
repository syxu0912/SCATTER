"""BAM-to-single-cell-TSS production entry point."""
from __future__ import annotations
from dataclasses import asdict
import hashlib
import html
import json
from pathlib import Path
import time

import pysam
import yaml

from . import __version__
from .calling import call_clusters
from .clustering import SplitConfig
from .io import read_barcodes
from .model import bundled_model
from .preflight import inspect_bam, inspect_fasta, inspect_gtf, resolve_index_path, validate_contig_overlap
from .reference import build_reference


def input_profile(bam, requested, evidence_mate='auto', tso_tag='ZT'):
    counts = dict(paired=0, single=0, read1=0, read2=0, selected=0, tag=0)
    records = []
    with pysam.AlignmentFile(str(bam), 'rb') as handle:
        if handle.header.to_dict().get('HD', {}).get('SO') != 'coordinate':
            raise ValueError('Input BAM must be coordinate-sorted')
        for read in handle.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary: continue
            counts['paired' if read.is_paired else 'single'] += 1
            counts['read1'] += int(read.is_read1); counts['read2'] += int(read.is_read2)
            records.append(read)
            if len(records) >= 20000: break
    if not records: raise ValueError('BAM has no primary mapped alignments')
    if counts['paired'] and counts['single']:
        raise ValueError('BAM mixes paired and single alignment flags; provide one consistent layout')
    layout = 'paired' if counts['paired'] else 'single'
    profile = requested
    if profile == 'auto':
        profile = 'r2only' if layout == 'single' or not counts['read1'] else 'paired'
    if profile == 'paired' and layout != 'paired':
        raise ValueError('--layout paired requires paired-end BAM flags; use r2only or auto for single R2 BAM')
    mate = evidence_mate
    if mate == 'auto': mate = 'read1' if profile == 'paired' else ('read2' if layout == 'paired' else 'single')
    if (layout == 'single' and mate != 'single') or (layout == 'paired' and mate == 'single'):
        raise ValueError('Selected evidence mate does not match BAM flags')
    for read in records:
        selected = (read.is_read1 if mate == 'read1' else read.is_read2 if mate == 'read2' else not read.is_read1 and not read.is_read2)
        if selected:
            counts['selected'] += 1; counts['tag'] += int(read.has_tag(tso_tag))
    return dict(profile=profile, alignment_layout=layout, evidence_mate=mate,
                endpoint_rule='physical-5prime' if profile == 'paired' else 'transcript-5prime', counts=counts)


def run(args):
    tick = time.monotonic()
    if not 0 <= args.threshold <= 1: raise ValueError('Threshold must be between 0 and 1')
    if args.threads < 1 or args.batch_size < 1: raise ValueError('Threads and batch size must be positive')
    for name in ('bam', 'fasta', 'gtf'):
        if not Path(getattr(args, name)).is_file(): raise FileNotFoundError(f'{name}: {getattr(args, name)}')
    model, model_metadata = bundled_model()
    if args.model: model, model_metadata = Path(args.model).resolve(), {'custom_model':True}
    if not model.is_file(): raise FileNotFoundError(model)
    profile = input_profile(args.bam, args.layout, args.evidence_mate, args.tso_tag)
    source = args.tso_source
    if args.mode == 'off': source = 'none'
    elif source == 'auto':
        if profile['counts']['selected'] and profile['counts']['tag'] == profile['counts']['selected']:
            source = 'tag'
        elif profile['profile'] == 'r2only':
            raise ValueError('R2-only tso/mix requires a transferred TSO tag on each evidence record; use --mode off or supply --tso-source tag --tso-tag ZT')
        else: source = 'evidence'
    out = Path(args.outdir).resolve()
    if out.exists() and any(out.iterdir()): raise FileExistsError(f'Output directory is not empty: {out}')
    out.mkdir(parents=True, exist_ok=True)
    metadata = out/'metadata'; metadata.mkdir()
    status = metadata/'run_status.json'
    status.write_text(json.dumps({'status':'RUNNING','version':__version__})+'\n')
    try:
        reference = out/'reference'; reference.mkdir()
        reference_stats = build_reference(args.gtf, reference)
        if not reference_stats['genes'] or not reference_stats['transcripts']:
            raise ValueError('Annotation must contain usable gene/transcript records')
        reference_tss = reference/'reference_tss.tsv'
        fasta = Path(args.fasta).resolve()
        if not Path(str(fasta)+'.fai').is_file():
            local = reference/('genome.fa.gz' if fasta.suffix=='.gz' else 'genome.fa')
            local.symlink_to(fasta); pysam.faidx(str(local)); fasta = local
        fasta_summary = inspect_fasta(fasta)
        gtf_summary = inspect_gtf(args.gtf)
        contigs = validate_contig_overlap(args.bam, fasta_summary, gtf_summary)
        with pysam.AlignmentFile(args.bam, 'rb') as bam_header:
            if not set(bam_header.references) & set(fasta_summary.references):
                raise ValueError('BAM and FASTA have no identical contig names; check chr prefixes and reference assembly')
        try: index = resolve_index_path(args.bam, args.bam_index)
        except Exception:
            if args.bam_index: raise
            index = metadata/'input.bam.bai'; pysam.index(str(Path(args.bam).resolve()), str(index))
        whitelist = None
        if args.barcodes:
            cells = read_barcodes(args.barcodes)
            whitelist = metadata/'barcodes.tsv'; whitelist.write_text(''.join(c+'\n' for c in cells))
        audit = inspect_bam(args.bam, index, None, args.cb_tag, args.umi_tag,
            str(whitelist) if whitelist else None, args.protocol_check_records,
            alignment_layout=profile['alignment_layout'], evidence_mate=profile['evidence_mate'],
            reference_tss_path=reference_tss, library_strand='gene', tso_source=source,
            tso_tag=args.tso_tag, min_mapq=args.min_mapq)
        config = {k:v for k,v in vars(args).items() if k!='handler'}
        config.update(resolved_profile=profile, resolved_tso_source=source,
                      resolved_cb_tag=audit.selected_cb_tag, resolved_umi_tag=audit.selected_umi_tag)
        (metadata/'run_config.yaml').write_text(yaml.safe_dump(config,sort_keys=True))
        (metadata/'preflight_report.json').write_text(json.dumps({'status':'PASS','bam':asdict(audit),'contigs':contigs},indent=2)+'\n')
        from .cli import fingerprint
        fingerprints={key:fingerprint(getattr(args,key),'fast') for key in ('bam','fasta','gtf')}
        if args.barcodes: fingerprints['barcodes']=fingerprint(args.barcodes,'full')
        (metadata/'input_checksums.json').write_text(json.dumps(fingerprints,indent=2)+'\n')
        result = call_clusters(bam_path=args.bam,bam_index=str(index),fasta_path=str(fasta),
            reference_tss_path=str(reference_tss),model_path=str(model),outdir=str(out),
            alignment_layout=profile['alignment_layout'],evidence_mate=profile['evidence_mate'],
            endpoint_rule=profile['endpoint_rule'],assay_type='5prime',library_strand='gene',
            tso_source=source,tso_tag=args.tso_tag,cb_tag=audit.selected_cb_tag,umi_tag=audit.selected_umi_tag,
            barcode_whitelist=str(whitelist) if whitelist else None,min_mapq=args.min_mapq,
            min_gene_reads=args.min_gene_reads,min_cluster_reads=args.min_cluster_reads,
            threshold=args.threshold,tso_mode=args.mode,tso_sequence=args.tso_sequence,
            evidence_technical_prefix_bp=args.technical_prefix_bp,contigs=set(args.contigs.split(',')) if args.contigs else None,
            max_records=args.max_records,split_config=SplitConfig(max_width_bp=args.cluster_max_width_bp),
            device=args.device,batch_size=args.batch_size,threads=args.threads,
            sample_id=args.sample_id or Path(args.bam).stem)
        result.update(status='PASS',scatter_version=__version__,sample_id=args.sample_id or Path(args.bam).stem,
            model=model_metadata,diagnostic_subset=args.max_records is not None or args.contigs is not None,
            elapsed_seconds=time.monotonic()-tick)
        (metadata/'call_summary.json').write_text(json.dumps(result,indent=2)+'\n')
        (metadata/'call_config.yaml').write_text(yaml.safe_dump({'barcode_whitelist':str(whitelist) if whitelist else None}))
        (out/'run_summary.json').write_text(json.dumps(result,indent=2)+'\n')
        status.write_text(json.dumps({'status':'PASS','elapsed_seconds':result['elapsed_seconds']})+'\n')
        write_report(out,result)
        (out/'COMPLETE').write_text('SCATTER '+__version__+' BAM-to-cell-TSS pipeline completed.\n')
        print(json.dumps(result,indent=2))
        return 0
    except Exception as exc:
        status.write_text(json.dumps({'status':'FAILED','error':str(exc)})+'\n')
        raise


def write_report(out,result):
    import anndata as ad
    matrix=ad.read_h5ad(out/'count/scTSS_count_all.h5ad')
    qc=dict(cells=matrix.n_obs,retained_TSS=matrix.n_vars,total_TSS_UMIs=int(matrix.X.sum()),
            zero_TSS_cells=int((matrix.obs.total_counts==0).sum()),
            raw_clusters=result['raw_clusters'],model_passing_clusters=result['high_confidence_clusters'])
    (out/'metadata/matrix_qc.json').write_text(json.dumps(qc,indent=2)+'\n')
    rows=''.join(f'<tr><th>{html.escape(k)}</th><td>{v:,}</td></tr>' for k,v in qc.items())
    document='''<!doctype html><html lang="en"><meta charset="utf-8"><title>SCATTER report</title>
<style>body{font:16px system-ui;max-width:1000px;margin:50px auto;padding:0 24px;color:#182431}table{border-collapse:collapse}td,th{padding:10px 24px;border-bottom:1px solid #ddd;text-align:left}a{color:#2c7fb8}</style>
<h1>SCATTER single-cell TSS results</h1>'''
    document+=f'<p>Sample: {html.escape(str(result["sample_id"]))}; mode: {html.escape(result["tso_mode"])}; model threshold: {result["model_threshold"]}</p><table>{rows}</table>'
    document+='<p>Counts are corrected UMI molecules. Zero-count requested cells are retained. TSS model support is a prediction, not experimental proof.</p>'
    document+='<ul><li><a href="count/scTSS_count_all.h5ad">All cell-by-TSS counts (AnnData)</a></li><li><a href="ctss/all_ctss.h5ad">Cell-by-CTSS counts</a></li><li><a href="count/matrix/">Matrix Market files</a></li><li><a href="run_summary.json">Full run metadata</a></li></ul></html>'
    (out/'report.html').write_text(document)


def doctor(args):
    import importlib.metadata
    import torch
    model,metadata=bundled_model()
    print(json.dumps({'scatter':__version__,'default_model':metadata,'model_available':model.is_file(),
        'CUDA_available':torch.cuda.is_available(),
        'dependencies':{p:importlib.metadata.version(p) for p in ('numpy','pandas','scipy','anndata','pysam','torch')}},indent=2))
    return 0


def demo(args):
    from importlib.resources import files
    from .cli import main
    data=files('scatter').joinpath('data')
    return main(['run','--bam',str(data.joinpath(args.layout+'.bam')),
        '--fasta',str(data.joinpath('genome.fa')),'--gtf',str(data.joinpath('genes.gtf')),
        '--barcodes',str(data.joinpath('barcodes.tsv')),'--mode',args.mode,
        '--layout',args.layout,'--device','cpu','--threads','2','--outdir',args.outdir])
