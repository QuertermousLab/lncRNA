# *de novo* lncRNA discovery pipeline
-----
This is a shell/awk re-written [PLAR](http://www.weizmann.ac.il/Biological_Regulation/IgorUlitsky/PLAR) script for lncRNA discovery. 

> Requirements:
> bash, awk, [bedtools](https://bedtools.readthedocs.io) (>=2.26), [CPC2](http://cpc2.gao-lab.org/download.php) (beta or 1.x, needs Python/Biopython), [HMMER3](http://hmmer.org) ([with Pfam-A dataset](https://www.ebi.ac.uk/interpro/download/Pfam/))
>
> All of them are available from bioconda: `conda install -c bioconda bedtools cpc2 hmmer`. bbmap and tophat (gtf_to_fasta) are no longer needed.

[![996.icu](https://img.shields.io/badge/link-996.icu-red.svg)](https://996.icu) [![LICENSE](https://img.shields.io/badge/license-Anti%20996-blue.svg)](https://github.com/996icu/996.ICU/blob/master/LICENSE)

----

### Input

This script uses original output GTF file of [stringtie](https://ccb.jhu.edu/software/stringtie/). Attributes are read by name: `transcript_id` (unique), `FPKM` (transcript lines) and `reference_id` (set by stringtie for transcripts matching the [GENCODE](https://www.gencodegenes.org/human/release_29lift37.html) or NONCODE guide GTF).

### How it works

The filtering principle and threshold are based on [PLAR](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4576741/) with several steps simplification and improvement.

#### **For known lncRNAs:**

1. Separate reference lncRNA from input by matching the unique ID of  known lncRNAs of the reference.
> [GENCODE  release29 (GRCh37)](https://www.gencodegenes.org/human/release_29lift37.html) was used by default, another transcript ID list can be given with `-l`.
>
> !NOTE: Use the SAME version of ref GTF annotation for transcript assembly.

2. Filter the known lncRNA by FPKM (default>0.1, `-r`).

#### **For *de novo* lncRNAs:**

1. Get all *de novo* transcripts by removing the transcripts which have been assigned with reference ID.
2. Sort the rest of transcripts into single- or multi- exon groups.
3. Filter 1 - FPKM and length: transcript length (sum of exons) has to be >=200nt and FPKM has to be >1 for single-exon or >0.1 for multi-exon to be kept. Transcripts on non-`chr` contigs are dropped.
4. Filter 2 - repeat/gap regions: Any transcripts which have at least 1 exon overlaps with a genomic repeat sequence (same strand) or gap region (any strand) for >=50% fraction will be removed.
5. Filter 3 - CDS: Any transcripts which overlap with protein coding sequences on the same strand will be removed.
6. Filter 4 - distance to genes: For the single-exon transcripts which have distance < 2000bp to protein coding gene and <500 for multi-exon transcripts will be removed.
7. Filter 5 - CPC2 and HMMER: Non-coding potential and protein coding potential of the transcripts will be calculated by *CPC2* and *hmmsearch* (3-frame translation against Pfam-A, E-value<0.001). Transcripts passed both will be kept and could be consider as potential *de novo* lncRNAs.

> Thresholds are set in the configuration block at the top of `lncRNA.sh`.
>
> Transcripts without strand (`.`, e.g. unstranded single-exon transcripts) are compared with features on both strands in filters 2-4 and their coding potential is tested on both strands.

Filters 2-4 run in parallel, CPC2 runs in parallel chunks and HMMER runs once for single- and multi-exon transcripts together.

Then the script combines known and *de novo* lncRNAs together to a final.gtf.

### Options

help message can be shown by `lncRNA.sh -h`

```shell
Usage: lncRNA.sh <options> -g|-n <transcripts.gtf>

### INPUT: GTF file (stringtie output) ###
### bedtools/CPC2/HMMER(Pfam-A) required ###

Options:
    -g  GENCODE reference was used for transcript assembly
    -n  NONCODE reference was used for transcript assembly
    -s  single-exon FPKM cutoff (default: 1)
    -m  multi-exon FPKM cutoff (default: 0.1)
    -r  referenced transcript FPKM cutoff (default: 0.1)
    -l  known lncRNA transcript ID list, plain or .gz
        (default: bed/gencode.v29lift37.long_noncoding_RNAs.list.gz for -g,
         bed/NONCODE.list.gz for -n)
    -o  output directory (default: ./)
    -t  threads for CPC2/HMMER (default: number of CPUs, max 16)
    -G  genome FASTA (default: $LNC_GENOME)
    -C  CPC2.py command (default: $LNC_CPC2 or CPC2.py in PATH)
    -P  Pfam-A.hmm (default: $LNC_PFAM)
    -k  keep intermediate files
    -v  print version
    -h  print this help message
```

Genome, CPC2 and Pfam paths can also be set once in the configuration block of the script or with the environment variables `LNC_GENOME`, `LNC_CPC2`, `LNC_PFAM` (and `LNC_BED_DIR` for the annotation BED files). The genome has to be the one used for read mapping (GRCh37/hg19 with `chr` names for the bundled annotation).

### Example

**1.** Reads mapping and transcripts assembly. [trans_assemble.sh](https://github.com/zhaoshuoxp/Pipelines-Wrappers#trans_assemblesh) can be used. i.e.

```shell
wget https://raw.githubusercontent.com/zhaoshuoxp/Pipelines-Wrappers/master/trans_assemble.sh
chmod 755 trans_assemble.sh
./trans_assemble.sh test_R1.fastq.gz test_R2.fastq.gz rf
```

The output test.gtf can be used for lncRNA discovery.

**2.** lncRNA filtering:

```shell
git clone https://github.com/zhaoshuoxp/lncRNA
cd lncRNA
chmod 755 lncRNA.sh
./lncRNA.sh -g -t 16 -G GRCh37.p13.genome.fa -P Pfam-A.hmm test.gtf
```

###  Output

All results will be stored in the output directory (`-o`, default ./). Log will be printed during running.

* test_final.gtf: combined final transcripts in GTF format.
* test_known_lncRNA_f1.gtf: FPKM>1 filtered reference lncRNA transcripts.
* test_multi_exon_f5.gtf: all filter passed multi-exon *de novo* lncRNA transcripts.
* test_single_exon_f5.gtf: all filter passed single-exon *de novo* lncRNA transcripts.
* non_add: for NONCODE reference transcripts, the pipeline will put them into *de novo* filters and output the filtered transcripts in this subdirectory.
* test_lncRNA_status.tsv: one line per input transcript with its type, length, FPKM, source (known/denovo/reference) and the filter which removed it (`pass` for final lncRNAs).

Further transcript deduplication could be performed if you merged multiple GTFs before or after running this pipeline:

```shell
wget https://raw.githubusercontent.com/zhaoshuoxp/Pipelines-Wrappers/master/GTF_rmdup.sh
chmod 755 GTF_rmdup.sh
./GTF_rmdup.sh final.gtf final_uniq.gtf
```

> NOTE:[UCSC Genome Browser utility](http://hgdownload.soe.ucsc.edu/admin/exe/) gtfToGenePred and genePredToBed are required.



### Changes in v2.0

Bug fixes (they change results, see below):

* Minus-strand transcripts were tested for coding potential on the wrong strand: tophat `gtf_to_fasta` does not reverse-complement them, so CPC2 and HMMER saw the antisense sequence. Sequences are now extracted with `bedtools getfasta -s -split`.
* Gap regions in `hg19.repeat_gap.bed` have no strand and never matched with `bedtools -s`, so filter 2 ignored them.
* Unstranded transcripts (`.`) never matched anything with `-s` in filters 2-3, and got distance `-1` from `closestBed` in filter 4, so all of them were removed there. They are now compared with both strands.
* BED intervals were built from 1-based GTF starts, and transcript length was computed as `end - start` (1 nt short).
* The BED files were never gunzipped (`for i in files:` loop), the script failed silently on any error and always reported "Run succeed", multi-exon length cutoff used the single-exon variable, FPKM/IDs were read from fixed `"`-separated column positions (broken when e.g. `ref_gene_name` is missing), temporary files in `./` made parallel runs in one directory overwrite each other.

Performance (chr21/chr22 test set, 8142 transcripts, 16 threads): 28.7 min -> 1.4 min.

* HMMER: `hmmsearch -Z <number of Pfam models>` (identical hits and E-values to `hmmscan`) on peptide chunks in parallel, instead of `hmmscan` (re-reads the whole Pfam database for every query). This is ~20x faster.
* CPC2 in parallel chunks, single- and multi-exon transcripts in one pass, filters 2-4 in parallel, GTF parsed once.
* bbmap (Java) and tophat are no longer needed: translation is done in awk (validated against Biopython).

On the test set, known lncRNAs were identical; 95 de novo lncRNAs were gained and 118 lost, all explained by the strand/gap fixes above. `<prefix>_lncRNA_status.tsv` records which filter removed each transcript.

----

Author [@zhaoshuoxp](https://github.com/zhaoshuoxp)  

Nov 8 2019  