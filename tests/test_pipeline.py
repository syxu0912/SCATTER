"""Real bundled-model end-to-end tests: BAM -> molecule counts -> AnnData/MTX."""
from importlib.resources import files
import gzip
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pytest
from scipy.io import mmread

from scatter.cli import main, make_parser
from scatter.model import bundled_model


@pytest.mark.parametrize('layout', ['paired', 'r2only'])
@pytest.mark.parametrize('mode,features,molecules', [('off',3,7), ('mix',3,7), ('tso',2,5)])
def test_real_bam_to_single_cell_results(tmp_path, layout, mode, features, molecules):
    out = tmp_path/'results'
    data = Path(str(files('scatter').joinpath('data')))
    assert not (data/'genome.fa.fai').exists()
    assert not (data/(layout+'.bam.bai')).exists()
    assert main(['demo','--layout',layout,'--mode',mode,'--outdir',str(out)]) == 0
    assert (out/'COMPLETE').is_file()
    assert json.loads((out/'metadata/run_status.json').read_text())['status']=='PASS'
    a = ad.read_h5ad(out/'count/scTSS_count_all.h5ad')
    c = ad.read_h5ad(out/'ctss/all_ctss.h5ad')
    assert a.shape == c.shape == (3,features)
    assert a.obs_names.tolist() == c.obs_names.tolist() == ['cellB-1','cellA-1','cellC-1']
    assert int(a.X.sum()) == int(c.X.sum()) == molecules
    assert a.obs.loc['cellC-1','total_counts'] == 0
    assert int(a.layers['counts'].sum()) == molecules
    assert np.allclose(a.var.TSS_start, a.var.start_1-1)
    assert (a.var.deeptssfilter_score >= .9).all()
    assert set(a.obs.sample_id)=={layout}
    assert a.uns['scatter']['model_sequence_length_bp']==1000
    assert a.uns['scatter']['model_threshold']==.9
    assert a.var.Strand.tolist().count('-') == 1
    summary = json.loads((out/'run_summary.json').read_text())
    assert summary['model_inference']['sequence_length_bp']==1000
    assert summary['endpoint_rule']==('physical-5prime' if layout=='paired' else 'transcript-5prime')
    matrix = mmread(out/'count/matrix/matrix.mtx.gz').tocsr()
    assert matrix.shape==(features,3)
    assert np.array_equal(matrix.toarray(),a.X.T.toarray())
    with gzip.open(out/'count/matrix/barcodes.tsv.gz','rt') as f:
        assert f.read().splitlines()==a.obs_names.tolist()
    if mode != 'tso':
        gene = a[:,a.var.gene_id=='G1']
        assert np.allclose(np.asarray(gene.layers['usage'].sum(axis=1)).ravel(), [1,1,0])
        assert ad.read_h5ad(out/'count/scTSS_count_two.h5ad').shape==(3,2)
    assert not (data/'genome.fa.fai').exists()
    assert not (data/(layout+'.bam.bai')).exists()
    assert main(['demo','--layout',layout,'--mode',mode,'--outdir',str(out)]) == 2


def test_model_is_the_selected_shift_checkpoint():
    path, meta = bundled_model()
    assert path.is_file() and meta['sequence_length_bp']==1000 and meta['training_sigma_bp']==10
    assert meta['sha256']=='58410fdd8a913a95f3075dd5c35741b5418d20effab8338d7f329c077c4ce80b'


def test_primary_run_requires_an_explicit_cell_selection():
    parser=make_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(['run','--bam','x','--fasta','x','--gtf','x','--mode','off','--outdir','x'])


def test_r2only_missing_tso_evidence_is_rejected_before_writing(tmp_path):
    from scatter.pipeline import input_profile
    import pysam
    source=Path(str(files('scatter').joinpath('data/r2only.bam')))
    dest=tmp_path/'untagged.bam'
    with pysam.AlignmentFile(source,'rb') as src,pysam.AlignmentFile(dest,'wb',header=src.header) as out:
        for r in src:
            r.set_tag('ZT',None);out.write(r)
    data=source.parent
    output=tmp_path/'out'
    assert main(['run','--bam',str(dest),'--fasta',str(data/'genome.fa'),'--gtf',str(data/'genes.gtf'),
                 '--barcodes',str(data/'barcodes.tsv'),'--mode','mix','--layout','r2only','--outdir',str(output)])==2
    assert not output.exists()


def test_gzip_gtf_and_called_cell_list(tmp_path):
    data=Path(str(files('scatter').joinpath('data')))
    gtf=tmp_path/'genes.gtf.gz';barcodes=tmp_path/'cells.tsv.gz'
    with gzip.open(gtf,'wt') as f:f.write((data/'genes.gtf').read_text())
    with gzip.open(barcodes,'wt') as f:f.write('barcode\n'+(data/'barcodes.tsv').read_text())
    out=tmp_path/'out'
    assert main(['run','--bam',str(data/'paired.bam'),'--fasta',str(data/'genome.fa'),
        '--gtf',str(gtf),'--barcodes',str(barcodes),'--layout','auto','--mode','off',
        '--sample-id','gzip-test','--device','cpu','--threads','1','--outdir',str(out)])==0
    a=ad.read_h5ad(out/'count/scTSS_count_all.h5ad')
    assert a.shape==(3,3) and int(a.X.sum())==7
    assert set(a.obs.sample_id)=={'gzip-test'}


def test_length_aware_scoring_invalid_windows_and_minus_strand(tmp_path):
    import pysam
    from scatter.model import score_loci
    source=Path(str(files('scatter').joinpath('data/genome.fa')))
    text=''.join(source.read_text().splitlines()[1:])
    text=text[:3500]+'N'+text[3501:]
    fasta=tmp_path/'ref.fa';fasta.write_text('>chr1\n'+text+'\n');pysam.faidx(str(fasta))
    loci=[('chr1',1501,'+'),('chr1',8501,'-'),('chr1',3501,'+'),('chr1',200,'+'),('unknown',1000,'+')]
    rows=[dict(chrom=c,summit_1=p,strand=s) for c,p,s in loci]
    scores,status,meta=score_loci(rows,fasta,bundled_model()[0],device='cpu',batch_size=2,threads=1)
    assert scores[loci[0]]>.99 and np.isclose(scores[loci[0]],scores[loci[1]])
    assert status[loci[2]]=='ambiguous_bases' and scores[loci[2]] is None
    assert status[loci[3]]=='out_of_bounds' and scores[loci[3]] is None
    assert status[loci[4]]=='missing_contig' and scores[loci[4]] is None
    assert meta['sequence_length_bp']==1000 and meta['unscorable_loci']==3
