#!/usr/bin/perl -w

# Contact: Dinh Diep
# Version 1.5

use strict;
use warnings;
use Getopt::Std;

my @reads = ();
my ($encodedFqName1, $encodedFqName2);
my $align_mode = 'S';
my $cpu = 2;
my $fqName = "Sequence1";
my $trim3=0;
my $trim5=0;
my $qualtrim=0;
my $nomap = 0;
my $noencode = 0;
my $keep_bam = 1;
my $TMP_DIR;
my $qual_base = 64;
my $maxMismatches=5;
my $score_min= "L,-0.6,-0.6";
my $gem_allow="-m 0.04 -e 0.04 -s 0 --unique-mapping";
my $min_identity = 0.95;

my ($template_fwd, $template_rev, $template_idx);
my ($soap2_exe, $bowtie2_exe, $bwa_exe, $last_dir, $gem_dir) = (0,0,0,0,0);

my %rcTable;
$rcTable{'A'}='T';
$rcTable{'T'}='A';
$rcTable{'G'}='C';
$rcTable{'C'}='G';
$rcTable{'N'}='N';
$rcTable{'R'}='Y';
$rcTable{'Y'}='R';
$rcTable{'M'}='K';
$rcTable{'K'}='M';
$rcTable{'S'}='S';
$rcTable{'W'}='W';

my %chrSizes;
my %countReads;
my $missed_guesses = 0;
my $total_reads = 0;
my $total_bases = 0;
my $total_mbases = 0;

&main;
exit;

sub main(){
        $template_idx = $ARGV[0];
        die("Need to have faidx file: $template_idx\n") if(!getChromSizes());
        my $start_time = time;
        my $map_file = $ARGV[1];
        $fqName = $map_file;
        $fqName =~ s/.combined.sorted//g;
        my $end_time = time;
        my $time_taken = $end_time - $start_time;
        print "Finished mapping and converting reads in $time_taken\n";
        sortsam($map_file);
        undef %chrSizes;
        my $total_mreads = 0;
        foreach my $chr(keys %countReads){
                print $chr, "\t", $countReads{$chr}, "\n";
                $total_mreads += $countReads{$chr};
        }
        print "Total sequences\t$total_reads\n";
        print "Total bases\t$total_bases\n";
        print "Total mapped bases\t$total_mbases\n";
        print "Total wrong guesses\t", $missed_guesses, "\n";
        undef %countReads;
}

sub processhit{
	my $array = shift;
	my @f = @{$array};
	my ($id, $orig_seq, $guess) = split /\|/, $f[0];
	$f[0] = $id;
	$f[9] = $orig_seq;
	my @qual = split "", $f[10];
	if($f[2] =~ s/_Watson//){
		if($f[1] & 16){
			return if($guess ne "R");
			# make everything on Watson maps to forward
			$f[9] = revComp($orig_seq);
			$f[10] = join("", reverse(@qual));
			#$f[0] = $f[0].":R";
		}else{
			return if($guess ne "F");
			#$f[0] = $f[0].":F";
		}
		$f[1] = 0;
	}elsif($f[2] =~ s/_Crick//){
		if($f[1] & 16){
			return if($guess ne "F");
			#$f[0] = $f[0].":R";
			# reverse on Crick is forward on Watson, make it reverse on Watson
			$f[9] = revComp($orig_seq);
			$f[10] = join("", reverse(@qual));
		}else{
			return if($guess ne "R");
			#$f[0] = $f[0].":F";
		}
		$f[1] = 16;
	}
	###-----------------BEGIN deal with CIGAR---------------------###
	my @CIGARS;
		# I - insertion to the reference
		# D - deletion to the reference
		# N - skipped region from the reference
		# S - soft clip on the read
		# H - hard clip on the read
		# P - padding - silent deletion from the padded reference sequence
	my @values = split /(\d+)/, $f[5];
	my $hard_clip_beg = 0;
	my $hard_clip_end = 0;
	for(my $i = 1; $i<scalar(@values)-1; $i=$i+2){
		if($values[$i+1] eq "H"){ #hide hard clippings
			$hard_clip_beg = $values[$i] if($i == 1);
			$hard_clip_end = $values[$i] if($i > 1);
			next;
		}
		push(@CIGARS, $values[$i].$values[$i+1]);
	}
	# match orig-seq to CIGARS
	my $clipped_s = substr($f[9], $hard_clip_beg, length($f[9])-$hard_clip_beg-$hard_clip_end);
	$f[5] = join("", @CIGARS);
	$f[9] = $clipped_s;
	###-----------------END deal with CIGAR---------------------###

	###----------------BEGIN deal with alt refs (with SNPs) ----###
	if(!$chrSizes{$f[2]}){
		my ($chr, $offset) = split "," , $f[2];
		$f[2] = $chr;
		$f[3] = $f[3] + $offset;
	}
	###----------------END deal with alt refs ------------------###

	$total_mbases+=length($f[9]);
	$countReads{$f[2]}++;
	return @f[0...10];
}

sub sortsam{
	my $map_file = shift;
	my $combined_sorted_map_file = $fqName.".combined.sorted";
	
	open(SAM_IN, "$combined_sorted_map_file") || die("Error in opening $combined_sorted_map_file.");

	my $unsorted = $fqName . ".sam";
	open ( SAM_OUT , ">$unsorted") || die("Error writing to file $unsorted\n");
	my $last_line = 0;
	my @last_fields;
	my $last_cnt = 0;
	while(my $line =  <SAM_IN>){
		my @fields = split(/\t/, $line);
		next if(scalar(@fields) < 5 or $fields[0] =~ /^@/);
		if($fields[5]){ # CIGAR
			next if($fields[5] eq "*");
			my @cigar_str = split(/(\d+)/, $fields[5]);
			my $identity = 0;
			for(my $i = 0; $i < scalar(@cigar_str); $i++){
				if($cigar_str[$i] eq "M"){
					$identity+= $cigar_str[$i-1];
				}
			}
			next if($identity/length($fields[9]) < $min_identity); # skip read if less than $min_identity
		}
		if(!$last_line){
			$last_line = $line;
			@last_fields = @fields;
			$last_cnt = 0;
			next;
		}
		if($fields[0] eq $last_fields[0]){
			my $score = $fields[13];
			$score =~ s/\D//g;
			my $last_score = $last_fields[13];
			$last_score =~ s/\D//g;
			if($score > $last_score){
				#2nd line is a better hit
				$last_line = $line;
				@last_fields = @fields;
				$last_cnt = 0;
			}elsif($score == $last_score){
				#check if one is from alternate reference (alternate references are not listed in fai file)
				#two equivalent good hits, increment last_cnt
				my $last_ref = $last_fields[2];
				my $ref = $fields[2];
				$last_ref =~ s/_Crick//;
				$last_ref =~ s/_Watson//;
				$ref =~ s/_Crick//;
				$ref =~ s/_Watson//;
				$last_cnt++ if($chrSizes{$ref} and $chrSizes{$last_ref});
			}else{
				#1st line is a better hit, do nothing
			}
		}else{	
			if($last_cnt ne 0){ # not a unique best hit
				undef $last_line;
				undef @last_fields;
				next;
			}
			if(my @hit = processhit(\@last_fields)){
				print SAM_OUT join("\t", @hit), "\n";
			}else{
				print "Wrong guess for ", $last_fields[0], "\n";
				$missed_guesses++;
			}
			$last_line = $line;
			@last_fields = @fields;
		}
	}
	#print the last line
	if($last_line and $last_cnt eq 0){
		if(my @hit = processhit(\@last_fields)){
			print SAM_OUT join("\t", @hit), "\n";
		}else{
			print "Wrong guess for ", $last_fields[0], "\n";
			$missed_guesses++;
		}
	}
	close(SAM_IN);
	close(SAM_OUT);
}

sub revComp{
	my $seq = shift;
	my $rcSeq='';
	for(my $i=0; $i<length($seq); $i++){
		$rcSeq = $rcTable{uc(substr($seq,$i,1))} . $rcSeq;
	}
	return $rcSeq;
}

sub guess_strand{
	my $seq = shift;
	my %baseCounts;
	$baseCounts{'A'}=0.001;
	$baseCounts{'T'}=0.001;
	$baseCounts{'G'}=0.001;
	$baseCounts{'C'}=0.001;
	while(my $base = chop($seq)){
		$baseCounts{$base}++;
	}
	if($baseCounts{'T'}/$baseCounts{'C'} > $baseCounts{'A'}/$baseCounts{'G'}) {
		return "F";
	}else{
		return "R";
	}
}

sub getChromSizes{
	return 0 if($template_idx !~ m/fai/);
	open(GENOME_INDEX, "$template_idx") || return 0;
	while(my $line = <GENOME_INDEX>){
		chomp($line);
		my @f = split "\t", $line;
		my $cur_chr = $f[0];
		$cur_chr =~ s/_Watson//;
		$cur_chr =~ s/_Crick//;
		$f[0] = $cur_chr;
		$chrSizes{$cur_chr} = $f[1];
	}
	close(GENOME_INDEX);
	return 1;
}

