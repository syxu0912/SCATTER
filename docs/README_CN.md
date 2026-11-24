# SCATTER：从BAM到单细胞TSS结果

本软件整合已验证的paired/R2-only取坐标、TSO/OFF/MIX聚簇、最新Shift模型过滤和逐细胞UMI计数。默认模型为 **hg38、1000 bp、σ=10、N1/N2负例的Shift模型**，默认保留分数≥0.9的簇。

## 安装与自检

建议Python≥3.10，使用独立虚拟环境，先安装适合本机CPU/CUDA的PyTorch，再安装本软件：

```bash
python -m pip install .
scatter doctor
scatter demo --outdir demo --layout paired --mode mix
```

模型权重已随软件打包，安装后无需另外下载模型。`doctor`会核对模型SHA256并显示依赖与GPU可用情况。

## 输入

1. 坐标排序的BAM，可复用已有BAI/CSI索引；未提供索引时在输出目录创建。
2. 与该BAM比对版本一致的FASTA和GTF；GTF可为gzip格式。默认模型按hg38训练。
3. 实际有效细胞的barcode列表，例如Cell Ranger输出的`filtered_feature_bc_matrix/barcodes.tsv.gz`，顺序和`-1`后缀原样保留。

BAM应包含校正后的`CB`、`UB`和唯一基因分配`GX`。`--cb-tag`、`--umi-tag`允许显式选择标签。使用原始`CR/UR`需要提供白名单；本软件不重新做cell calling或barcode纠错。确需计入所有合格BAM barcode时，必须显式指定`--all-barcodes`。

## 一条命令运行

paired示例：

```bash
scatter run --bam sample.bam --fasta genome.fa --gtf genes.gtf.gz \
  --barcodes barcodes.tsv.gz --layout paired --mode mix \
  --outdir results/sample --sample-id sample --threads 4 --device auto
```

R2-only示例：

```bash
scatter run --bam r2only.bam --fasta genome.fa --gtf genes.gtf \
  --barcodes barcodes.tsv.gz --layout r2only --mode off \
  --outdir results/r2only
```

三种模式的含义：

- `tso`：只使用TSO阳性的证据形成候选，不要求候选靠近参考TSS。
- `off`：使用全部合格证据，不进行TSO选择；当前沿用实验流程，候选需在该基因参考TSS的1000 bp以内。
- `mix`：TSO与全读段路线独立聚簇，允许TSO重叠证据救回OFF候选，再合并重叠簇；UMI计数不把两条路线重复相加。

这三个模式决定“用哪些证据找簇”，不是CamoTSS中TC/CTSS输出模式的对应替换：**这里默认同时输出TSS簇与CTSS两个层面的结果**。

默认每条路线每基因至少50条合格reads、每原始簇至少50条reads，均可通过`--min-gene-reads`和`--min-cluster-reads`调整。两者是reads阈值；最终矩阵是UMI纠错后的分子数。原始簇按现有分谷规则限制宽度≤200 bp。

## paired与R2-only方向

- paired默认使用R1物理5′端，要求本项目使用的长R1协议；预检的主导R1长度至少100 bp。
- R2-only记录基因/RNA链。RNA正链取R2比对区间左端，RNA负链取右端，与本项目修正后的原始簇保持一致；不按片段长度虚构缺失R1的位置。
- `--layout auto`可根据BAM flags判断，但正式运行建议写明布局。带paired flags但只想选择R2的输入也可以指定`--layout r2only --evidence-mate read2`。

TSO来源自动选择：若证据记录均有`ZT`标签，优先使用；否则paired使用物理5′soft clip检测。保留了技术前缀的R1可用`--technical-prefix-bp`说明前缀长度。若比对软件已经去掉TSO、BAM里也没有转移标签，软件无法从BAM恢复丢失的TSO证据。

**R2-only运行tso/mix需要提前转移的TSO标签**。缺少标签时直接报错，不会悄悄切换成off；可直接运行off。

## 核心输出

| 文件 | 内容 |
|---|---|
| `count/scTSS_count_all.h5ad` | 细胞×全部保留TSS簇，稀疏整数UMI矩阵 |
| `count/scTSS_count_two.h5ad` | 保留至少两个TSS簇的基因对应矩阵 |
| `ctss/all_ctss.h5ad` | 细胞×CTSS矩阵 |
| `ctss/all_ctss_two.h5ad` | 至少两个CTSS特征的基因对应矩阵 |
| `count/matrix/`、`ctss/matrix/` | Seurat/通用Matrix Market格式，特征×细胞 |
| `clusters/` | 各路线原始簇与模型保留簇，含模型分数 |
| `metadata/`、`run_summary.json` | 参数、输入指纹、模型版本、链与端点规则、计数QC |
| `report.html` | 可直接打开的运行报告 |
| `work/ctss_evidence.sqlite` | 重建逐细胞UMI计数所需的证据库 |

零TSS计数的白名单细胞也会保留，以便和原基因表达矩阵对齐。所有矩阵使用一致的细胞顺序。`.X`和`.layers['counts']`保存整数计数；`.layers['usage']`保存已观察到的细胞内比例：TSS按基因归一化，CTSS按所属TSS簇归一化。它不是插补值，也不是差异使用检验。

`.obs`包含细胞计数QC与sample_id；`.var`包含gene_id、坐标、链、模型分数等；`.uns['scatter']`保存软件、模型、阈值和坐标约定。详见[输出定义](OUTPUTS.md)。

```python
import anndata as ad
a = ad.read_h5ad('results/sample/count/scTSS_count_all.h5ad')
print(a.shape)
print(a.var[['gene_id', 'summit_1', 'deeptssfilter_score']].head())
usage = a.layers['usage']
```

## 使用边界

软件基于参考DNA打分，相同基因组位点的分数不会因细胞类型而改变；单细胞特异信息来自每个细胞的UMI计数。模型保留表示序列模型支持的候选，不是实验确认的真TSS。

当前0 bp、四种外部数据任一支持的规则用于研究中的事后证据评价；运行软件不要求这些数据，也不把外部证据标签作为模型保留条件。未来如另外注释外部证据，应明确细胞系与匹配规则。

不同样本单独运行的`SCATTER_...` ID是各自局部编号，不能直接把它们作为共同位点合并。跨样本差异TSS分析需要先确定共同位点目录并一致重计数。当前版本提供可供Scanpy/Seurat及后续TSS使用分析的计数和注释，不自动宣称完成差异TSS统计推断。

使用新/空输出目录；软件拒绝覆盖已有结果。大型BAM需预留SQLite证据库与矩阵输出空间。`--max-records`、`--contigs`只用于小规模验证，结果中会明确标记为诊断子集。
