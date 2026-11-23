import csv
import importlib.util
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[1]

def module(filename):
    spec=importlib.util.spec_from_file_location(filename,ROOT/'bin'/filename)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m

def test_complete_candidate_integrity(tmp_path):
    m=module('mm_candidate_integrity.py')
    p=tmp_path/'in.tsv';out=tmp_path/'qc.tsv'
    p.write_text('fragment_id\thi\tnh\tpaired\tr1_chrom\tr2_chrom\n'+
                 'f1\t1\t2\t0\tchr1\t.\n'+
                 'f1\t2\t2\t0\tchr2\t.\n'+
                 'f2\t1\t3\t0\tchr1\t.\n')
    result=m.validate(str(p),str(out))
    assert result['complete']==1 and result['incomplete']==1

def test_reverse_homology_coordinates(tmp_path):
    m=module('mm_homology_index.py')
    paf=tmp_path/'self.paf';out=tmp_path/'homology.tsv'
    paf.write_text('chr1\t100\t10\t18\t-\tchr2\t100\t20\t28\t8\t8\t60\tcg:Z:8M\n')
    assert m.build(str(paf),str(out),1,.9)==2
    rows=list(csv.DictReader(open(out),delimiter='\t'))
    assert rows[0]['source_start0']=='10' and rows[0]['strand']=='-'

def test_quality_likelihood_ranks_correct_origin(tmp_path):
    m=module('mm_quality_scores.py')
    fasta=tmp_path/'r.fa';fasta.write_bytes(b'>chr1\nACGTACGT\n>chr2\nACATACGT\n')
    # offset for chr1=6; chr2=21
    (tmp_path/'r.fa.fai').write_text('chr1\t8\t6\t8\t9\nchr2\t8\t21\t8\t9\n')
    fa=m.Faidx(str(fasta))
    def rec(ch):return ['f','0',ch,'1','0','4M','*','0','0','ACGT','IIII']
    p1,b1=m.score_record(rec('chr1'),fa)
    p2,b2=m.score_record(rec('chr2'),fa)
    assert p1>p2 and b1==b2==4
    fa.close()

def test_em_conserves_fragment_mass_and_uses_uniques():
    m=module('mm_locus_em.py')
    def groups():
        for f in ('f1','f2'):
            yield f, [({'candidate_id':f+'|HI=1'},'A',0.),({'candidate_id':f+'|HI=2'},'B',0.)], 'COMPLETE'
    names,theta,unique,frags,cands,source,resp,rejected,iters,converged=m.em(groups(),{'A':10},alpha=1.01)
    assert names==['A','B']
    assert theta[0]>.5
    assert np.allclose(resp.reshape(2,2).sum(axis=1),1)
    assert not rejected

def test_homology_projection_requires_full_span(tmp_path):
    m=module('mm_homology_propose.py')
    p=tmp_path/'h.tsv'; c=tmp_path/'cand.tsv'; out=tmp_path/'proposals.tsv'
    p.write_text('source_chrom\tsource_start0\tsource_end0\ttarget_chrom\ttarget_start0\ttarget_end0\tstrand\n'+
                 'chr1\t10\t30\tchr2\t50\t70\t-\n')
    c.write_text('fragment_id\tcandidate_id\tr1_chrom\tr1_pos1\tr1_cigar\n'+
                 'f1\tf1|HI=1\tchr1\t11\t5M\n'+
                 'f2\tf2|HI=1\tchr1\t28\t5M\n')
    m.propose(str(p),str(c),str(out))
    x=list(csv.DictReader(open(out),delimiter='\t'))
    assert len(x)==1 and x[0]['target_start0']=='65' and x[0]['status']=='UNVALIDATED_PROPOSAL'
