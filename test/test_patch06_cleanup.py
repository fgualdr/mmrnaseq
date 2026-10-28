#!/usr/bin/env python3
from pathlib import Path
import json

root = Path(__file__).resolve().parents[1]

# Schema / manifest remain parseable and on the validated runtime line.
schema = json.loads((root / 'nextflow_schema.json').read_text())
config = (root / 'nextflow.config').read_text()
assert "nf-boost@0.7.0" in config
assert "nextflowVersion = '!>=25.10.4'" in config
assert "name            = 'fgualdr/mmrnaseq'" in config
assert "version         = '4.0.0'" in config

# Unsupported/dead branches must not remain in active runtime configuration.
runtime_files = [
    root / 'nextflow.config',
    root / 'nextflow_schema.json',
    root / 'main.nf',
    root / 'lib/WorkflowRnaseq.groovy',
    root / 'workflows/rnaseq.nf',
    root / 'subworkflows/local/prepare_genome.nf',
    root / 'conf/modules.config',
    root / 'conf/test.config',
    root / 'conf/test_full.config',
]
for token in [
    'star_rsem', 'hisat2_index', 'rsem_index', 'hisat2_build_memory',
    'stringtie_ignore_gtf', 'skip_stringtie', 'pseudo_aligner',
    'em_locus_merge_distance', 'em_hotspot_distance', 'label_overlap',
    'em_iter', 'em_seed', 'BAM_EM_PREP', 'BAM_EM_OUT',
]:
    for f in runtime_files:
        assert token not in f.read_text(), (token, f)

# Deleted implementation files must really be gone.
deleted = [
    'bin/clean_work_files.sh',
    'bin/em_algorithm_bedpe_sm_bis.py',
    'bin/prepare-for-rsem.py',
    'modules/local/bam_em.nf',
    'modules/local/bam_em_out.nf',
    'modules/local/bam_em_prep.nf',
    'modules/local/bam_fragment.nf',
    'modules/local/bedtools_genomecov.nf',
    'modules/local/clean_utils.nf',
    'modules/local/rsem_merge_counts.nf',
    'modules/local/umitools_prepareforrsem.nf',
    'modules/nf-core/hisat2/align/main.nf',
    'modules/nf-core/hisat2/build/main.nf',
    'modules/nf-core/hisat2/extractsplicesites/main.nf',
    'modules/nf-core/rsem/calculateexpression/main.nf',
    'modules/nf-core/rsem/preparereference/main.nf',
    'modules/nf-core/stringtie/stringtie/main.nf',
    'subworkflows/local/quantify_rsem.nf',
    'subworkflows/nf-core/fastq_align_hisat2/main.nf',
]
for f in deleted:
    assert not (root / f).exists(), f

# RSEM had been abused solely to create transcript FASTA; gffread now owns that job.
assert (root / 'modules/local/gtf_to_transcript_fasta.nf').exists()
prep = (root / 'subworkflows/local/prepare_genome.nf').read_text()
assert 'GTF_TO_TRANSCRIPT_FASTA' in prep
assert 'RSEM_' not in prep
assert 'rsem_index' not in prep

# The renamed UMI helper must be actually imported and used.
wf = (root / 'workflows/rnaseq.nf').read_text()
assert "from '../modules/local/umitools_prepareforsalmon'" in wf
assert 'UMITOOLS_PREPAREFORSALMON (' in wf
assert (root / 'bin/prepare-for-salmon.py').exists()

# PATCH03/PATCH05 BAM semantics and QC diagnostics must survive cleanup unchanged in meaning.
sub = (root / 'subworkflows/local/bam_filter_em.nf').read_text()
for token in [
    'bam_quantitative = ch_quantitative_bam',
    'bam_highconf = ch_highconf_bam',
    'bam_full = ch_full_bam',
    'mm_density_cv',
    'mm_qc =',
    'MM_QC_MERGE',
]:
    assert token in sub, token
for token in [
    'BAM_FILTER_EM.out.bam_quantitative',
    'BAM_FILTER_EM.out.bam_highconf',
    'BAM_FILTER_EM.out.bam_full',
    'BAM_SORT_STATS_SAMTOOLS_MM_HIGHCONF',
    'BAM_SORT_STATS_SAMTOOLS_MM_FULL',
]:
    assert token in wf, token
mods = (root / 'conf/modules.config').read_text()
for token in [
    '.mm.quantitative.sorted', '.mm.highconf.sorted', '.mm.full.sorted',
    'MM_QC_MERGE', 'MM_BAYES_RESOLVER', 'MM_BAM_REALIZE'
]:
    assert token in mods, token

# Do not let coordinate duplicate marking act on arbitrary STAR primary placements.
assert 'Coordinate Picard duplicate marking before multimapper filtering/resolution is disabled' in wf
assert 'This applies to both --skip_em false and --skip_em true' in wf
assert 'BAM_FILTER_EM(ch_bam_genome_bai, PREPARE_GENOME.out.chrom_sizes)' in wf
assert (root / 'modules/local/mm_skip_resolution.nf').exists()
assert 'ch_quantitative_bam = BAM_FILTER.out.bam' not in sub

# Remove accidental debug channel materialization.
assert '.view(' not in wf

# The MultiQC positional signature must no longer contain removed aligner slots.
mult = (root / 'modules/local/multiqc.nf').read_text()
assert "path ('hisat2/*')" not in mult
assert "path ('rsem/*')" not in mult
assert 'pseudoaligner' not in mult
assert 'ch_hisat2_multiqc' not in wf
assert 'ch_rsem_multiqc' not in wf
assert 'ch_pseudoaligner_' not in wf
assert 'ch_star_multiqc.collect' in wf and 'ch_salmon_multiqc.collect' in wf

# Resource configs should contain one, not duplicate, high-memory override.
for f in [root / 'conf/ieo.config', root / 'conf/ieo.slurm.config']:
    assert f.read_text().count("withLabel: 'process_high_memory'") == 1, f

# modules.json should no longer advertise removed modules/subworkflow.
modules = json.loads((root / 'modules.json').read_text())
blob = json.dumps(modules).lower()
for token in ['hisat2/align', 'hisat2/build', 'hisat2/extractsplicesites',
              'rsem/calculateexpression', 'rsem/preparereference',
              'stringtie/stringtie', 'fastq_align_hisat2']:
    assert token not in blob, token

print('PASS PATCH06 workflow cleanup rebased on PATCH05 Bayesian QC')
