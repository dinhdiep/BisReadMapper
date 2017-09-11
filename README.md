# Outline for pipeline

1) use genomePrep.pl to generate the in-silico bisulfite converted references, C->T for Watson, and G->A for Crick strands. Note that bisulfite conversion makes the two strands non-complementary.

2) use bwa to create index files from the reference sequence for alignment. NOTE: Both strands (*.bis.CT and *.bis.GA) can be concatenated into one file and only one index needs to be created so long as the aligner can support larger index files.

3) use samtools to generate the *.fai file from the reference sequence.

4) generate the list_files

```
<sample id> <dir> <read1.fq | read1.fq,read2.fq | *.sam> <phred> <clonal method> <adaptor r1> <adaptor r2>
```

5) generate the list_paths file

```
# Required 
scripts_dir="/home/ddiep/scripts/BisReadMapper/src"
samtools="/home/ddiep/softwares/samtools-0.1.18/samtools"
bamUtils="bam"
cpg_list="/media/3TB_AG/bisHg19/hg19_lambda.cpg.positions.txt"
ref_fai="/media/3TB_AG/bisHg19/hg19.fa.fai"
ref_fa="/media/3TB_AG/bisHg19/hg19.fa"

template_fwd="/media/3TB_AG/bisHg19/hg19_lambda.bis.CT"
template_rev="/media/3TB_AG/bisHg19/hg19_lambda.bis.GA"
mapper="/home/ddiep/softwares/bwa-0.7.5a/bwa"
trimgalore="/home/ddiep/softwares/trim_galore_latest/trim_galore"

# for snp calling
samtools_snp="/home/ddiep/softwares/samtools-0.1.8/samtools"
snp_list="NA"
ref_dbsnp="/media/3TB_AG/bisHg19/snp134_snv.txt"

# for targeted sequencing
target_bed="NA"

```

6) run MasterBisReadMapper.pl with list_sam_files and list_paths as inputs.

# Sample usage

```
Usage: genomePrep.pl genome.fa[.gz] context=[cg/all] convert=[yes/no]
convert: yes = convert genome, no = don't convert genome
context: CG = CG only, ALL = all C context
```
This program will generate the CpG position files when the context option is set to CG or ALL


```
Usage: MasterBisReadMapper.pl -i </path/to/list_sam_file> -s </path/to/list_path_file> [options] &> log
Required:
	-i <list_files>    : <list_sam_files> is a table (tab separated values) of all the fastq and sam files to be processes[Required]
                                 <sample id> <dir> <read1.fq | read1.fq,read2.fq | *.sam> <phred> <clonal method> <adaptor r1> <adaptor r2> <TRIM mode>
                                 Note: umi in sam files are indicated in the reads name field, separated by '_' and is only part with no [0-9]
                                 Note: umi are first X number of bases in read 1, use UMI:X for clonal method, where X is length of UMI.
        -s <list_paths>        : <list_paths> is a file that lists all of the required paths[Required]
Optional settings:
	-c [yes/no]            : split processed files by chromosomes or not[Default no]
        -v [yes/no]            : indicate whether to call SNPs or not[Default no]
        -b [yes/no]            : indicate whether to generate BED format file for methylation frequencies[Default no]
        -p [yes/no]            : keep pileup yes or no[Default no]
        -m [yes/no]            : Map only, do not make methylFreq [Default no]
        -d <mindepth>          : minimum depth for BED and frCorr[Default 5]
===================================================================

```
