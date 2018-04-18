# Outline for pipeline

To start, use genomePrep.pl to generate the in-silico bisulfite converted references, C->T for Watson, and G->A for Crick strands. Note that bisulfite conversion makes the two strands non-complementary. Next, use '''bwa''' to create the index files from the reference sequences (i. e. genome.fa file) for alignment. Then, use '''samtools''' to quickly generate the genome index file from the reference sequences. After the previous three preparation steps are finished, we need to generate the master paths file (i. e. list_of_paths.txt) similar to the one below:

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

Generate a table with the following columns for each sample to be processed. Each sample can have more than one row for more than one fastq file. Multiple samples can be processed together in one table for batch jobs but multiple samples can be processed simultaneously within separate tables and in separate jobs. See the usage below for how to start a mapping job.

```
<sample id> <dir> <read1.fq | read1.fq,read2.fq | *.sam> <phred> <clonal method> <adaptor r1> <adaptor r2>
```

# Sample usage

The following program will generate the bisulfite reference files.
```
Usage: genomePrep.pl genome.fa[.gz] context=[cg/all] convert=[yes/no]
convert: yes = convert genome, no = don't convert genome
context: CG = CG only, ALL = all C context
```

The following program will perform the alignment steps and generate the methylation frequency files.
```
Usage: MasterBisReadMapper.pl -i </path/to/input.table> -s </path/to/list_of_paths_file> [options] &> log
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

```
