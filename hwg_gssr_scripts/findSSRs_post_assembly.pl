#!/usr/bin/env perl
################################################################################
# Author: Meg Staton & Stephen Ficklin
# Date: Jan 23rd, 2014
# Version: 1
#
# DESCRIPTION
# -----------
# This script identifies simple sequence repeats (SSRs) and calls primers from
# the sequences in a fastq formatted file.
#
# Dependencies:
# ------------
# Perl must have access to the packages:
# Getopt::Long
# Bio::SeqIO
# Excel::Writer::XLSX
# All are available from CPAN.
#
# Also path to the primer3 executable and primer3 config files must be specified
# in the global variables section of the script.
#
# Usage:
# -----
# Usage: findSSRs.pl <arguments>
#
# The list of arguments includes:
#
# -f|--fasta_file <fasta_file>
# Required.  The file of the sequences to be searched.
#
# -m|--masked_file <masked_fasta_file>
# Required.  A soft-masked version of the fasta file (soft masked means low
# complexity sequences are in lower case bases.)
#
# Output:
# ------
# Eight output files are produced:
#
# <input-file-name>.ssr.fasta
# A fasta file with sequences with a single identified SSR.
#
# <input-file-name>.ssr_multi_seqs.fasta
# A fasta with sequences with more than one identified SSR.
#
# <input-file-name>.ssr_stats.txt
# A text file of statistics about the SSRs discovered.
#
# <input-file-name>.ssr_report.txt
# A tab-delimited file with each SSR.  The columns are sequence name,
# motif, number of repeats, start position and end position.
#
# <input-file-name>.ssr_report.xlsx
# A excel file with SSR results and stats
#
# <input-file-name>.di_primer_report.txt
# A tab-delimited file with sequences with a 2-bp SSR motif.  Columns are
# sequence name, motif, start position, end position, left primer,
# right primer, left primer Tm, right primer Tm, amplicon size, full 
# sequence, masked sequence
#
# <input-file-name>.tri_primer_report.txt
# A tab-delimited file with sequences with a 3-bp SSR motif.  Columns are
# sequence name, motif, start position, end position, left primer,
# right primer, left primer Tm, right primer Tm, amplicon size, full 
# sequence, masked sequence
#
# <input-file-name>.tetra_primer_report.txt
# A tab-delimited file with sequences with a 4-bp SSR motif.  Columns are
# sequence name, motif, start position, end position, left primer,
# right primer, left primer Tm, right primer Tm, amplicon size, full
# sequence, masked sequence
#
#
# Details:
# -------
# By default the script finds:
# 2 bp motifs repeated from 8 to 40 times,
# 3 bp motifs repeated from 7 to 30 times,
# 4 bp motifs repeated from 6 to 20 times,
#
# The script only reports SSRs that are not within 15 bases of either
# end of the sequence, in order to allow for primer design.
#
# These parameters may be changed in the "GLOBAL PARAMETERS" part of
# the script.
#


use strict;

#-------------------------------------------------------------------------------
#  DEPENDENCIES
#-------------------------------------------------------------------------------

use Getopt::Long;
use Bio::SeqIO;
use Excel::Writer::XLSX;

#-------------------------------------------------------------------------------
# GLOBAL PARAMETERS
#-------------------------------------------------------------------------------

#--------------
# REPEAT IDENTIFICATION PARAMETERS
# Specify Motif Frequency
# Motifs that occur less frequently than indicated below will be ignored.
# A 0 indicates that this motif length will be ignored.

our $MIN_REPS_2bp = 8;
our $MIN_REPS_3bp = 7;
our $MIN_REPS_4bp = 6;

our $MAX_REPS_2bp = 40;
our $MAX_REPS_3bp = 30;
our $MAX_REPS_4bp = 20;

# SSRs at the beginning or end of a sequence prevents proper primers design.
# This is how close we will allow an SSR to be to the ends of the sequence.
our $LENGTH_FROM_END = 15;

#------------
# PRIMER PARAMETERS

my $PRIMER3 = "/lustre/projects/staton/software/primer3-2.3.6/src/primer3_core";
my $PRIMER3_CONFIG = "/lustre/projects/staton/software/primer3-2.3.6/src/primer3_config/";

my $PRIMER_OPT_SIZE="20";  # default 20
my $PRIMER_MIN_SIZE="18";  # default 18
my $PRIMER_MAX_SIZE="25";  # default 27

my $PRIMER_NUM_NS_ACCEPTED = "0";  # default 0

my $PRIMER_PRODUCT_SIZE_RANGE = "100-200";

my $PRIMER_OPT_TM = "60.0";
my $PRIMER_MIN_TM = "55.0";
my $PRIMER_MAX_TM = "65.0";

my $PRIMER_MIN_GC = "40";
my $PRIMER_MAX_GC = "60";

my $PRIMER_MAX_POLY_X = "3";
my $PRIMER_GC_CLAMP = "2";

#-------------------------------------------------------------------------------
# GLOBAL HASHES
#-------------------------------------------------------------------------------
# This makes life much easier than passing a bunch of hash refs all over the place.
#
#-------------------------------------------------------------------------------
# Data structures:

## CONTIG_SSR_STARTS structure:
## key: contig_name
##  value: array of starts of SSRs in that contig
my %CONTIG_SSR_STARTS = ();

## SSR_STATS structure:
## key: ssr_id
##  value -> keys: MOTIF START END MOTIF_LENGTH NO_REPEATS
my %SSR_STATS = ();
my $SEQ_COUNT = 0;

# Set up the Motif specifications, based on the chosen motif types:.
my @MOTIF_SPECS;
push(@MOTIF_SPECS,[2, $MIN_REPS_2bp, $MAX_REPS_2bp, 'dinucleotides']);
push(@MOTIF_SPECS,[3, $MIN_REPS_3bp, $MAX_REPS_3bp, 'trinucleotides']);
push(@MOTIF_SPECS,[4, $MIN_REPS_4bp, $MAX_REPS_4bp, 'tetranucleotides']);

#-------------------------------------------------------------------------------
# Generating statistics
my %MOTIFS = ('|AT|TA|'                   => 0,
              '|AG|GA|CT|TC|'             => 0,
              '|AC|CA|TG|GT|'             => 0,
              '|GC|CG|'                   => 0,

              '|AAT|ATA|TAA|ATT|TTA|TAT|' => 0,
              '|AAG|AGA|GAA|CTT|TTC|TCT|' => 0,
              '|AAC|ACA|CAA|GTT|TTG|TGT|' => 0,

              '|CCA|CAC|CCA|TGG|GTG|TGG|' => 0,
              '|GGC|GCG|CGG|GCC|CCG|CGC|' => 0,
              '|AGG|GAG|GGA|CCT|CTC|TCC|' => 0,

              '|ATG|TGA|GAT|CAT|ATC|TCA|' => 0,
              '|AGT|GTA|TAG|ACT|CTA|TAC|' => 0,
              '|AGC|GCA|CAG|GCT|CTG|TGC|' => 0,
              '|ACG|CGA|GAC|CGT|GTC|TCG|' => 0);

my %MOTIFLEN = ('2'  => 0,
                '3'  => 0,
                '4'  => 0);

my %MOTIFLEN_w_PRIMERS = ('2'  => 0,
                          '3'  => 0,
                          '4'  => 0);

my $SSR_COUNT = 0;
my $SSR_w_PRIMER_COUNT = 0;

#-------------------------------------------------------------------------------
#  EXECUTE
#-------------------------------------------------------------------------------
main();
#-------------------------------------------------------------------------------
#  PUBLIC SUBROUTINES
#-------------------------------------------------------------------------------
# Function Name:  main()

sub main{
    my $fasta_file;
    my $masked_file;
    my $project;

    my $p3_input;
    my $p3_output;

    my $ssr_out;  # the tab-delimited output file
    my $ssr_xlsx;
    my $fasta_out;
    my $fasta_out_multi;
    my $stats_out;

    my $di_primer_out;
    my $tri_primer_out;
    my $tetra_primer_out;

    Getopt::Long::Configure ('bundling');
    GetOptions('f|fasta_file=s'  => \$fasta_file,
               'm|masked_file=s' => \$masked_file,
               'p|project=s'     => \$project);

    ## Check that all required parameters have been included
    if(!$fasta_file){ print "A fasta file is required.\n"; _printUsage(); exit;}
    if(!$masked_file){ print "A masked file is required.\n"; _printUsage(); exit;}

    ## Check that fasta file exists
    if(! -e $fasta_file) { print "Fasta file $fasta_file does not exist\n"; exit; }
    if(! -e $masked_file) { print "Masked file $masked_file does not exist\n"; exit; }

    $p3_input  = "$fasta_file.p3in.txt";
    $p3_output = "$fasta_file.p3out.txt";

    $ssr_out        = "$fasta_file.ssr_report.txt";
    $ssr_xlsx       = "$fasta_file.ssr_report.xlsx";
    $fasta_out      = "$fasta_file.ssr_filtered.fasta";
    $fasta_out_multi = "$fasta_file.ssr_multi_seqs.fasta";
    $stats_out      = "$fasta_file.ssr_stats.txt";

    $di_primer_out = "$fasta_file.di_primer_report.txt";
    $tri_primer_out = "$fasta_file.tri_primer_report.txt";
    $tetra_primer_out = "$fasta_file.tetra_primer_report.txt";

	##---------------------------------------------------------------
    print "finding SSRs...\n";
    open(PI3, ">$p3_input") or die $!;
    my $p3_input_fh = *PI3;
    process_file($fasta_file, $masked_file, $ssr_out, $p3_input_fh);
    close PI3;

    print "running primer3...\n";
    print "$PRIMER3 < $p3_input > $p3_output\n";
    my $status = system("$PRIMER3 < $p3_input > $p3_output");

    print "creating Excel workbook...";
    my ($workbook,$formats) = createExcelWorkbook($ssr_xlsx);
    print "done.\n";

    print "parsing primer3...";
    parseP3_output($p3_output);
    print "done.\n";

	##---------------------------------------------------------------
	## Producing output - Excel and statistics

	# open filehandles
    open (DI, ">$di_primer_out");
    open (TRI, ">$tri_primer_out");
    open (TETRA, ">$tetra_primer_out");
    open (FASTAOUT, ">$fasta_out");
    open (FASTAMULTI, ">$fasta_out_multi");
    my $di_fh = *DI;
    my $tri_fh = *TRI;
    my $tetra_fh = *TETRA;
    my $fastaout_fh = *FASTAOUT;
    my $fastamulti_fh = *FASTAMULTI;
	#initiate_workbooks($di_fh, $tri_fh, $tetra_fh, $workbook, $formats, $project, $fastaout_fh, $fastamulti_fh);
    print "done.\n";
    close DI;
    close TRI;
    close TETRA;
    close FASTAOUT;
    close FASTAMULTI;

	#print "stats...\n";
	#my $worksheet_stats = printStats($stats_out, $workbook, $formats, $project);
	#
	#$worksheet_stats->activate();
	#$worksheet_stats->select();
	#$workbook->close();

    print "done.\n";
}

###############################################################

sub process_file{
    my $fasta_file  = $_[0]; # file name
    my $masked_file = $_[1]; # file name
    my $ssr_out     = $_[2]; # file name
    my $p3_input_fh = $_[3]; # file name

    my $seqio;
    my $seqioM;

    my $seqobj;
    my $seqobjM;

    my $seqname;
    my $seqnameM;

    my $seqstr;
    my $seqstrM;

    open (OUT, ">".$ssr_out) || die "ERROR cannot open $_\n";
    my $out_fh = *OUT;

    $seqio = Bio::SeqIO->new('-format' => 'fasta', -file => $fasta_file);
    $seqioM = Bio::SeqIO->new('-format' => 'fasta', -file => $masked_file);

    # Get seq obj from io stream
    while($seqobj = $seqio->next_seq){
        $seqobjM = $seqioM->next_seq;

        $SEQ_COUNT++;

        $seqname = $seqobj->id;        # get actual sequence as a string
        $seqnameM = $seqobjM->id;        # get actual sequence as a string

        $seqstr = $seqobj->seq();        # get actual sequence as a string
        $seqstrM = $seqobjM->seq();        # get actual sequence as a string

        if($seqname ne $seqnameM){
            die "masked sequence $seqnameM not in same order as regular sequence $seqname\n";
        }

        process_seq($seqname, $seqstr, $seqstrM, $out_fh, $p3_input_fh);
    }

}

###############################################################

sub process_seq{
    my $contig_name = shift;
    my $seq         = shift;
    my $seq_masked  = shift;
    my $out_fh      = shift;
    my $p3_input_fh = shift;


    my %seen;     # used to keep track of start positions we've already seen
    my $index;    # used to iterate through the 2D motif specs array

    ## LOOPB
    # iterate through the motif specifications
    for $index (0 .. (scalar @MOTIF_SPECS - 1)) {

        # holds the motif size (1,2,3,4,5 or 6)
        my $motifLength = $MOTIF_SPECS[$index][0];
        # holds the minimum number of repeats
        my $min_number_of_repeats = $MOTIF_SPECS[$index][1]-1;
        # holds the maximum number of repeats
        my $max_number_of_repeats = $MOTIF_SPECS[$index][2]-1;
        # the regular expression for looking for SSRs
        my $regex = "(([gatc]{$motifLength})\\2{$min_number_of_repeats,$max_number_of_repeats})";

        # run through the sequence and check for this motif spec
        while ($seq =~ /$regex/ig) {
            # Get the ssr and motif that were found
            my $ssr = $1;
            my $motif = lc $2;
			my $start_index = $-[0];
			my $end_index = $+[0];

			if(quality_check_ssr($contig_name, $ssr, $motif, $start_index, $end_index, $seq)){
				process_ssr($contig_name, $ssr, $motif, $start_index, $end_index, $seq, $seq_masked, $out_fh, $p3_input_fh);

			}

		}
	}
}



sub quality_check_ssr{
	my $contig_name = shift;
	my $ssr = shift;
 	my $motif = shift;
 	my $start_index = shift;
 	my $end_index = shift;
	my $seq = shift;

	##-------------------------------------
	## CHECKS to see if this is a good ssr
	my $flag_same_base = 0;
	my $flag_already_seen = 0;
	my $flag_too_close_to_end = 0;

	## Check #1
	## ignore SSRs that are the same base repeated
	if ($ssr !~ /^g+$/i &&
	    $ssr !~ /^a+$/i &&
	    $ssr !~ /^c+$/i &&
	    $ssr !~ /^t+$/i ) {
		$flag_same_base = 1;
	}

	# Check #2
	# Make sure this isn't an already called SSR in disguise
	# (a dinucleotide repeat posing as a tetranucleotide repeat, for instance)
	if (!exists $SSR_STATS{$contig_name."_ssr".$start_index} &&
	    !exists $SSR_STATS{$contig_name."_ssr".($start_index-1)} &&
	    !exists $SSR_STATS{$contig_name."_ssr".($start_index-2)} &&
	    !exists $SSR_STATS{$contig_name."_ssr".($start_index+1)} &&
	    !exists $SSR_STATS{$contig_name."_ssr".($start_index+2)}
	) {
		$flag_already_seen = 1;
	}

	# Check #3
	# Distance from end
    my $seqLen = length $seq;
    if($start_index >= $LENGTH_FROM_END && $end_index <= ($seqLen-$LENGTH_FROM_END)){
		$flag_too_close_to_end = 1;
	}

	if($flag_same_base && $flag_already_seen && $flag_too_close_to_end){
		return 1;
	}
	else{
		return 0;
	}

}



sub process_ssr{
	my $contig_name = shift;
	my $ssr = shift;
 	my $motif = shift;
 	my $start_index = shift;
 	my $end_index = shift;
	my $seq = shift;
	my $seq_masked = shift;
    my $out_fh      = shift;
    my $p3_input_fh = shift;

	##-------------------------------------
	## generate a few more stats and variables
	my $motif_len = length $motif;
	my $num_of_repeats = ($end_index-$start_index)/$motif_len;
    my $ssr_id = $contig_name."_ssr".$start_index;

	#print "$contig_name $ssr $motif\n";
	#print "Start: $start_index\n";
	#print "End: $end_index\n";

	##-------------------------------------
	## store in data structures

    $SSR_COUNT++;
    $MOTIFLEN{$motif_len}++;

    if(exists $CONTIG_SSR_STARTS{$contig_name}){
        push @{ $CONTIG_SSR_STARTS{$contig_name} }, $start_index;
    }
    else{
       $CONTIG_SSR_STARTS{$contig_name} = [$start_index];
    }

	$SSR_STATS{$ssr_id}{MOTIF}        = $motif;
	$SSR_STATS{$ssr_id}{START}        = $start_index;
	$SSR_STATS{$ssr_id}{END}          = $end_index;
	$SSR_STATS{$ssr_id}{MOTIF_LENGTH} = $motif_len;
	$SSR_STATS{$ssr_id}{NO_REPEATS}   = $num_of_repeats;
	$SSR_STATS{$ssr_id}{SEQ}          = $seq;
	$SSR_STATS{$ssr_id}{SEQM}         = $seq_masked;

	foreach my $group (keys %MOTIFS) {
		my $motifUC = uc($motif);
		if($group =~ /\|$motifUC\|/){
			# If this group contains this motif
			#print "Incrementing $group for $motif\n";
			$MOTIFS{$group}++;
		}
	}

	##-------------------------------------
	# print to text file of all repeats found and primer 3 input file

	printf $out_fh ("$contig_name\t$motif\t$num_of_repeats\t$start_index\t$end_index\n");
	_addToPrimer3InputFile($p3_input_fh, $contig_name, $start_index, $end_index, $seq_masked);


}

###############################################################
sub _addToPrimer3InputFile{
    my $p3_input_fh = shift;
    my $contig_name = shift;
    my $ssrStart    = shift;
    my $ssrEnd      = shift;
    my $seq         = shift;

	# change from soft mask to hard mask
	$seq =~ s/[actg]/N/g;

    my $len         = $ssrEnd-$ssrStart;

    printf $p3_input_fh ("SEQUENCE_ID=$contig_name\_ssr$ssrStart\n");
    printf $p3_input_fh ("SEQUENCE_TEMPLATE=$seq\n");
    printf $p3_input_fh ("SEQUENCE_TARGET=$ssrStart,$len\n");
    printf $p3_input_fh ("PRIMER_TASK=generic\n");
    printf $p3_input_fh ("PRIMER_PICK_LEFT_PRIMER=1\n");
    printf $p3_input_fh ("PRIMER_PICK_INTERNAL_OLIGO=0\n");
    printf $p3_input_fh ("PRIMER_PICK_RIGHT_PRIMER=1\n");
    printf $p3_input_fh ("PRIMER_OPT_SIZE=$PRIMER_OPT_SIZE\n");
    printf $p3_input_fh ("PRIMER_MIN_SIZE=$PRIMER_MIN_SIZE\n");
    printf $p3_input_fh ("PRIMER_MAX_SIZE=$PRIMER_MAX_SIZE\n");
    printf $p3_input_fh ("PRIMER_NUM_NS_ACCEPTED=$PRIMER_NUM_NS_ACCEPTED\n");
    printf $p3_input_fh ("PRIMER_PRODUCT_SIZE_RANGE=$PRIMER_PRODUCT_SIZE_RANGE\n");
    printf $p3_input_fh ("PRIMER_OPT_TM=$PRIMER_OPT_TM\n");
    printf $p3_input_fh ("PRIMER_MIN_TM=$PRIMER_MIN_TM\n");
    printf $p3_input_fh ("PRIMER_MAX_TM=$PRIMER_MAX_TM\n");
    printf $p3_input_fh ("PRIMER_MIN_GC=$PRIMER_MIN_GC\n");
    printf $p3_input_fh ("PRIMER_MAX_GC=$PRIMER_MAX_GC\n");
    printf $p3_input_fh ("PRIMER_MAX_POLY_X=$PRIMER_MAX_POLY_X\n");
    printf $p3_input_fh ("PRIMER_GC_CLAMP=$PRIMER_GC_CLAMP\n");
    printf $p3_input_fh ("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=$PRIMER3_CONFIG\n");
    printf $p3_input_fh ("=\n");

}

###############################################################
sub parseP3_output{
    my $p3_output = $_[0]; # file name

	# We are going to keep track of a weird phenomenon only seen in one
	# project - the generation of identical forward and reverse primers. The
	# sequences from this project were overlapping paired ends that were joined.
	# Apparently something went wrong and weird sequences were obtained, all of
	# which yield the identical primers.
	# This is not reported in the final stats, just as part of the standard output.
    my $identical_primer_cnt = 0;

    # The primers output file separates information about different sequences
    # with an equal sign on a single line. So, we want to set the file line
    # delimiter (for looping on the input file below) to a single equal sign
    # followed by a line feed.  This way were guranteed to have all the primer
    # information together per line
    local $/ = "=\n";

    open (P3O, $p3_output) || die "could not open $_\n";

    # Read in all of the lines of the input file
	my $primer_record;
    while ($primer_record = <P3O>) {
        my $start = "";
        my $seq_id = "";
        my $ssr_id = "";
        my $forward = "";
        my $reverse = "";
        my $product_size = "";
        my $left_tm = "";
        my $right_tm = "";

		if ($primer_record =~ /SEQUENCE_ID=(\S+)/) {
		    $ssr_id = $1;
		}
		# get the primary primers only
		if ($primer_record =~ /PRIMER_LEFT_0_SEQUENCE=(\S+)/) {
		    $forward = $1;
		}
		if ($primer_record =~ /PRIMER_RIGHT_0_SEQUENCE=(\S+)/) {
		    $reverse = $1;
		}
		if ($primer_record =~ /PRIMER_LEFT_0_TM=(\S+)/) {
		    $left_tm = $1;
		}
		if ($primer_record =~ /PRIMER_RIGHT_0_TM=(\S+)/) {
		    $right_tm = $1;
		}
		if ($primer_record =~ /PRIMER_PAIR_0_PRODUCT_SIZE=(\S+)/) {
		    $product_size = $1;
		}
		print "record\n";
		print "ssr_id $ssr_id\n";
		print "forward $forward\n";
		print "reverse $reverse\n";
		print "left_tm $left_tm\n";
		print "right_tm $right_tm\n";
		print "product_size $product_size\n\n";

		if(length $forward > 1){
			if($forward eq $reverse){
				print "FLAG: problem with $ssr_id\n";
				$identical_primer_cnt++;
			}
			else{
				$SSR_STATS{$ssr_id}{FORWARD} = $forward;
				$SSR_STATS{$ssr_id}{REVERSE} = $reverse;
				$SSR_STATS{$ssr_id}{PRODUCT_SIZE} = $product_size;
				$SSR_STATS{$ssr_id}{LEFT_TM} = $left_tm;
				$SSR_STATS{$ssr_id}{RIGHT_TM} = $right_tm;

				my $motif      = $SSR_STATS{$ssr_id}{MOTIF};
				my $ssrStart   = $SSR_STATS{$ssr_id}{START};
				my $ssrEnd     = $SSR_STATS{$ssr_id}{END};
				my $seq        = $SSR_STATS{$ssr_id}{SEQ};
				my $seq_masked = $SSR_STATS{$ssr_id}{SEQM};

				$ssr_id =~ /(\S+)_ssr\d+/;
				my $contig = $1;
			}
		}
	}

    print "total identical primers: $identical_primer_cnt\n";
}

#    my $di_fh = $_[1]; # file name
#    my $tri_fh = $_[2]; # file name
#    my $tetra_fh = $_[3]; # file name
#
#    my $workbook = $_[4]; # file name
#    my $formats  = $_[5]; # file name
#    my $project  = $_[6]; # file name
#
#    my $fastaout_fh = $_[7]; # file name
#    my $fastamulti_fh = $_[8]; # file name
#
#    my $start;
#    my $seq_id;
#    my $ssr_id;
#    my $forward;
#    my $reverse;
#    my $product_size;
#    my $left_tm;
#    my $right_tm;
#
#    my $di_worksheet = $workbook->add_worksheet("Di");
#    $di_worksheet->set_column('A:A', 60, $formats->{text});
#    $di_worksheet->set_column('F:G', 30, $formats->{text});
#    #$di_worksheet->set_column('J:J', 100, $formats->{text});
#    $di_worksheet->write('A1', "Dinucleotide Repeats for $project", $formats->{header});
#    $di_worksheet->write('A2', 'Sequence Name', $formats->{header});
#        $di_worksheet->write('B2', 'Motif', $formats->{header});
#        $di_worksheet->write('C2', '# Repeats', $formats->{header});
#        $di_worksheet->write('D2', 'Start', $formats->{header});
#        $di_worksheet->write('E2', 'End', $formats->{header});
#        $di_worksheet->write('F2', 'Forward Primer', $formats->{header});
#        $di_worksheet->write('G2', 'Reverse Primer', $formats->{header});
#        $di_worksheet->write('H2', 'Forward Tm', $formats->{header});
#        $di_worksheet->write('I2', 'Reverse Tm', $formats->{header});
#        $di_worksheet->write('J2', 'Fragment Size', $formats->{header});
#        #$di_worksheet->write('J2', 'Sequence', $formats->{header});
#    my $di_index = 3;
#
#    my $tri_worksheet = $workbook->add_worksheet("Tri");
#    $tri_worksheet->set_column('A:A', 60, $formats->{text});
#    $tri_worksheet->set_column('F:G', 30, $formats->{text});
#    #$tri_worksheet->set_column('J:J', 100, $formats->{text});
#    $tri_worksheet->write('A1', "Trinucleotide Repeats for $project", $formats->{header});
#    $tri_worksheet->write('A2', 'Sequence Name', $formats->{header});
#        $tri_worksheet->write('B2', 'Motif', $formats->{header});
#        $tri_worksheet->write('C2', '# Repeats', $formats->{header});
#        $tri_worksheet->write('D2', 'Start', $formats->{header});
#        $tri_worksheet->write('E2', 'End', $formats->{header});
#        $tri_worksheet->write('F2', 'Forward Primer', $formats->{header});
#        $tri_worksheet->write('G2', 'Reverse Primer', $formats->{header});
#        $tri_worksheet->write('H2', 'Forward Tm', $formats->{header});
#        $tri_worksheet->write('I2', 'Reverse Tm', $formats->{header});
#        $tri_worksheet->write('J2', 'Fragment Size', $formats->{header});
#        #$tri_worksheet->write('J2', 'Sequence', $formats->{header});
#    my $tri_index = 3;
#
#    my $tetra_worksheet = $workbook->add_worksheet("Tetra");
#    $tetra_worksheet->set_column('A:A', 60, $formats->{text});
#    $tetra_worksheet->set_column('F:G', 30, $formats->{text});
#    #$tetra_worksheet->set_column('J:J', 100, $formats->{text});
#    $tetra_worksheet->write('A1', "Tetranucleotide Repeats for $project", $formats->{header});
#    $tetra_worksheet->write('A2', 'Sequence Name', $formats->{header});
#        $tetra_worksheet->write('B2', 'Motif', $formats->{header});
#        $tetra_worksheet->write('C2', '# Repeats', $formats->{header});
#        $tetra_worksheet->write('D2', 'Start', $formats->{header});
#        $tetra_worksheet->write('E2', 'End', $formats->{header});
#        $tetra_worksheet->write('F2', 'Forward Primer', $formats->{header});
#        $tetra_worksheet->write('G2', 'Reverse Primer', $formats->{header});
#        $tetra_worksheet->write('H2', 'Forward Tm', $formats->{header});
#        $tetra_worksheet->write('I2', 'Reverse Tm', $formats->{header});
#        $tetra_worksheet->write('J2', 'Fragment Size', $formats->{header});
#        #$tetra_worksheet->write('J2', 'Sequence', $formats->{header});
#    my $tetra_index = 3;
#
#                    my $multi_flag = 0;
#                    ## skip contigs with more than one ssr
#                    if(scalar @{ $CONTIG_SSR_STARTS{$contig}} == 1){
#                        print $fastaout_fh ">$contig $motif.$ssrStart-$ssrEnd\n$seq\n";
#                        #print "\t$forward\n";
#                        $SSR_w_PRIMER_COUNT++;
#                        if(length $motif == 2){
#                            print $di_fh join("\t", $contig, $motif, $ssrStart, $ssrEnd, $forward, $reverse, $left_tm, $right_tm, $product_size, $seq, $seq_masked);
#                            print $di_fh "\n";
#                            my $tmp = $MOTIFLEN_w_PRIMERS{2};
#                            $tmp++;
#                            $MOTIFLEN_w_PRIMERS{2} = $tmp;
#                            my $cnt = ($ssrEnd-$ssrStart+1)/2;
#
#                            $di_worksheet->write("A$di_index", $contig, $formats->{text});
#                                $di_worksheet->write("B$di_index", $motif, $formats->{text});
#                                $di_worksheet->write("C$di_index", $cnt, $formats->{text});
#                                $di_worksheet->write("D$di_index", $ssrStart, $formats->{text});
#                                $di_worksheet->write("E$di_index", $ssrEnd, $formats->{text});
#                                $di_worksheet->write("F$di_index", $forward, $formats->{text});
#                                $di_worksheet->write("G$di_index", $reverse, $formats->{text});
#                                $di_worksheet->write("H$di_index", $left_tm, $formats->{text});
#                                $di_worksheet->write("I$di_index", $right_tm, $formats->{text});
#                                $di_worksheet->write("J$di_index", $product_size, $formats->{text});
#                                #$di_worksheet->write("J$di_index", $seq, $formats->{text});
#                            $di_index++;
#
#                            # Increment motif count
#                            foreach my $group (keys %MOTIFS) {
#                                my $motifUC = uc($motif);
#                                if($group =~ /\|$motifUC\|/){
#                                    # If this group contains this motif
#                                    my $tmp = $MOTIFS{$group}++;
#                                    $tmp++;
#                                    $MOTIFS{$group} = $tmp;
#                                }
#                            }# end foreach $group
#                        }
#                        elsif(length $motif == 3){
#                            print $tri_fh join("\t", $contig, $motif, $ssrStart, $ssrEnd, $forward, $reverse, $left_tm, $right_tm, $product_size, $seq, $seq_masked);
#                            print $tri_fh "\n";
#                            my $tmp = $MOTIFLEN_w_PRIMERS{3};
#                            $tmp++;
#                            $MOTIFLEN_w_PRIMERS{3} = $tmp;
#
#                            my $cnt = ($ssrEnd-$ssrStart+1)/3;
#                            $tri_worksheet->write("A$tri_index", $contig, $formats->{text});
#                                $tri_worksheet->write("B$tri_index", $motif, $formats->{text});
#                                $tri_worksheet->write("C$tri_index", $cnt, $formats->{text});
#                                $tri_worksheet->write("D$tri_index", $ssrStart, $formats->{text});
#                                $tri_worksheet->write("E$tri_index", $ssrEnd, $formats->{text});
#                                $tri_worksheet->write("F$tri_index", $forward, $formats->{text});
#                                $tri_worksheet->write("G$tri_index", $reverse, $formats->{text});
#                                $tri_worksheet->write("H$tri_index", $left_tm, $formats->{text});
#                                $tri_worksheet->write("I$tri_index", $right_tm, $formats->{text});
#                                $tri_worksheet->write("J$tri_index", $product_size, $formats->{text});
#                                #$tri_worksheet->write("J$tri_index", $seq, $formats->{text});
#                            $tri_index++;
#                        }
#                        elsif(length $motif == 4){
#							_printLineToWorksheet();
#                        }
#                    }
#                    else{
#                        print $fastamulti_fh ">$contig\n$seq\n";
#                    }
#                }
#            }
#        }
#    } # end while <INPUT>
#
#    close P3O;
#
#    return;
#}
#
#################################################################
#sub _printLineToWorksheet{
#	my $fh = shift;
#	my $index = shift;
#	my $worksheet = shift;
#
#	my $contig = shift;
#	my $motif = shift;
#	my $ssrStart = shift;
#	my $ssrEnd = shift;
#	my $forward = shift;
#	my $reverse = shift;
#	my $left_tm = shift;
#	my $right_tm = shift;
#	my $product_size = shift;
#	my $seq = shift;
#	my $seq_masked = shift;
#
#	print $fh join("\t", $contig, $motif, $ssrStart, $ssrEnd, $forward, $reverse, $left_tm, $right_tm, $product_size, $seq, $seq_masked);
#	print $fh "\n";
#	my $tmp = $MOTIFLEN_w_PRIMERS{4};
#	$tmp++;
#	$MOTIFLEN_w_PRIMERS{4} = $tmp;
#
#	my $cnt = ($ssrEnd-$ssrStart+1)/4;
#	$worksheet->write("A$index", $contig, $formats->{text});
#		$worksheet->write("B$index", $motif, $formats->{text});
#		$worksheet->write("C$index", $cnt, $formats->{text});
#		$worksheet->write("D$index", $ssrStart, $formats->{text});
#		$worksheet->write("E$index", $ssrEnd, $formats->{text});
#		$worksheet->write("F$index", $forward, $formats->{text});
#		$worksheet->write("G$index", $reverse, $formats->{text});
#		$worksheet->write("H$index", $left_tm, $formats->{text});
#		$worksheet->write("I$index", $right_tm, $formats->{text});
#		$worksheet->write("J$index", $product_size, $formats->{text});
#	$index++;
#
#}
#

################################################################
sub printStats{
    my $stats_out = $_[0]; # file name
    my $workbook  = $_[1]; # file name
    my $formats   = $_[2]; # file name
    my $project   = $_[3]; # file name


    ##--------------------------------------------------------------------
    ## calculate some info

    my $time = scalar localtime;                   # Get the current time

    ## count number of seqs with single or multiple SSRs
    my $SINGLE_SSR_COUNT_w_primers = 0;
    my $SINGLE_SSR_COUNT = 0;
    my $MULTI_SSR_COUNT = 0;

    foreach my $contig_name (keys %CONTIG_SSR_STARTS){
        my $starts = scalar @{ $CONTIG_SSR_STARTS{$contig_name} };
        if($starts == 1){
            $SINGLE_SSR_COUNT++;
            my @starts = @{ $CONTIG_SSR_STARTS{$contig_name} };
            my $start = $starts[0];
            my $ssr_id = $contig_name."_ssr".$start;
            if(exists $SSR_STATS{$ssr_id} && $SSR_STATS{$ssr_id}{FORWARD} =~ /\S/){
                $SINGLE_SSR_COUNT_w_primers++;
            }
        }
        elsif($starts > 1){ $MULTI_SSR_COUNT++; }
    }

    my $SSR_COUNT_w_primers = 0;
    foreach my $ssr_id (keys %SSR_STATS){
        if($SSR_STATS{$ssr_id}{FORWARD} =~ /\S/){
            $SSR_COUNT_w_primers++;
        }
    }

    ##--------------------------------------------------------------------
    ## print text file
    open (OUTS, ">".$stats_out) || die "ERROR cannot open $stats_out\n";

    print OUTS 'SSR Summary Report\n';
    print OUTS "Analsis of $SEQ_COUNT sequences\n";
    print OUTS "$time\n";
    print OUTS "Number of SSRs identified\t$SSR_COUNT\n";
    print OUTS "Number of sequences with 1 SSR: $SINGLE_SSR_COUNT\n";
    print OUTS "Number of sequences with more than one SSR: $MULTI_SSR_COUNT\n";
    print OUTS "\n";
    print OUTS "Base Pairs in Motif\tMin # Reps\tMax # Reps\n";
    print OUTS "--------------------------------------\n";
    print OUTS "2 (Dinucleotides)\t$MIN_REPS_2bp\t$MAX_REPS_2bp\n";
    print OUTS "3 (Trinucleotides)\t$MIN_REPS_3bp\t$MAX_REPS_3bp\n";
    print OUTS "4 (Tetranucleotides)\t$MIN_REPS_4bp\t$MAX_REPS_4bp\n";
    print OUTS "\n";
    print OUTS "Motif Patterns\tNumber of SSRs Found\n";
    print OUTS "--------------------------------------\n";
    my $group;
    foreach $group (sort {length $a <=> length $b} keys %MOTIFS){
        $group =~ s/^|//;
        $group =~ s/|$//;
        print OUTS "$group\t$MOTIFS{$group}\n";
    }
    print OUTS "\n";
    print OUTS "Motif Pattern Length\tNumber of SSRs Found\n";
    print OUTS "--------------------------------------\n";

    foreach $group (sort keys %MOTIFLEN){
        print OUTS "$group\t$MOTIFLEN{$group}\n";
    }

    print OUTS "SSRS with PRIMERS\n";
    print OUTS "Number of SSRs identified with successful primer design: $SSR_COUNT_w_primers\n";
    print OUTS "Number of sequences with 1 SSR and successful primer design: $SINGLE_SSR_COUNT_w_primers\n";
    print OUTS "Motif Pattern Length\tNumber of SSRs Found\n";
    print OUTS "--------------------------------------\n";

    foreach $group (sort keys %MOTIFLEN_w_PRIMERS){
        print OUTS "$group\t$MOTIFLEN_w_PRIMERS{$group}\n";
    }


    close OUTS;

    ##--------------------------------------------------------------------
    ## print excel file
    my $worksheet = $workbook->add_worksheet("Summary");

    $worksheet->set_column('A:A', 75, $formats->{text});
    $worksheet->set_column('B:B', 30, $formats->{text});

    $worksheet->write('A1',"SSR Summary Report for $project", $formats->{header});
    $worksheet->write('A2',"Analsis of $SEQ_COUNT sequences", $formats->{text});
    $worksheet->write('A3',"$time", $formats->{text});

    $worksheet->write('A4',"Number of SSRs identified", $formats->{text});
        $worksheet->write('B4',"$SSR_COUNT", $formats->{text});
    $worksheet->write('A5',"Number of sequences with 1 SSR", $formats->{text});
        $worksheet->write('B5',"$SINGLE_SSR_COUNT", $formats->{text});
    $worksheet->write('A6',"Number of sequences with more than one SSR", $formats->{text});
        $worksheet->write('B6',"$MULTI_SSR_COUNT", $formats->{text});

    $worksheet->write('A8','Base Pairs in Motif', $formats->{header});
        $worksheet->write('B8','Min # Reps', $formats->{header});
        $worksheet->write('C8','Max # Reps', $formats->{header});
    $worksheet->write('A9','2 (Dinucleotides)', $formats->{text});
        $worksheet->write('B9',"$MIN_REPS_2bp", $formats->{text});
        $worksheet->write('C9',"$MAX_REPS_2bp", $formats->{text});
    $worksheet->write('A10','3 (Trinucleotides)', $formats->{text});
        $worksheet->write('B10',"$MIN_REPS_3bp", $formats->{text});
        $worksheet->write('C10',"$MAX_REPS_3bp", $formats->{text});
    $worksheet->write('A11','4 (Tetranucleotides)', $formats->{text});
        $worksheet->write('B11',"$MIN_REPS_4bp", $formats->{text});
        $worksheet->write('C11',"$MAX_REPS_4bp", $formats->{text});

    $worksheet->write('A13','Motif Patterns', $formats->{header});
        $worksheet->write('B13','Number of SSRs Found', $formats->{header});
    my $group;
    my $i = 13;
    foreach $group (sort {length $a <=> length $b} keys %MOTIFS){
        $group =~ s/^|//;
        $group =~ s/|$//;
        $i++;
        $worksheet->write("A$i", $group, $formats->{text});
        $worksheet->write("B$i", $MOTIFS{$group}, $formats->{text});
    }

    $i++;
    $i++;
    $worksheet->write("A$i",'Motif Pattern Length', $formats->{header});
        $worksheet->write("B$i",'Number of SSRs Found', $formats->{header});
    foreach $group (sort keys %MOTIFLEN){
        $i++;
        $worksheet->write("A$i", "$group bp", $formats->{text});
        $worksheet->write("B$i", $MOTIFLEN{$group}, $formats->{text});
    }

    $i++;
    $i++;
    $worksheet->write("A$i",'SSRs with Primers', $formats->{header});
    $i++;
    $worksheet->write("A$i", "Number of SSRs identified with successful primer design", $formats->{text});
        $worksheet->write("B$i", $SSR_COUNT_w_primers, $formats->{text});
    $i++;
    $worksheet->write("A$i", "Number of sequences with 1 SSR and successful primer design", $formats->{text});
        $worksheet->write("B$i", $SINGLE_SSR_COUNT_w_primers, $formats->{text});

    $i++;
    $i++;
    $worksheet->write("A$i",'Motif Pattern Length (Only SSRs with Primers)', $formats->{header});
        $worksheet->write("B$i",'Number of SSRs Found', $formats->{header});
    foreach $group (sort keys %MOTIFLEN_w_PRIMERS){
        $i++;
        $worksheet->write("A$i", "$group bp", $formats->{text});
        $worksheet->write("B$i", $MOTIFLEN_w_PRIMERS{$group}, $formats->{text});
    }

    close OUTS;
    return $worksheet;

}

###############################################################

sub createExcelWorkbook{

    my $ssr_xlsx = $_[0];

    my $workbook;        # the excel workbook
    my %formats;
    my %header;
    my %text;
    my %bigheader;
    my %highlight;


    # Create an excel workbook
    $workbook = Excel::Writer::XLSX->new("$ssr_xlsx");

    # Setup the four formats that will be necessary for the excel spreadsheet
    %header = (font         => 'Calibri',
                size         => 12,
                bold         => 1,
                color        => 'black',
                align        => 'left',
                text_wrap    => 1);

    %text = (font         => 'Calibri',
                size         => 12,
                color        => 'black',
                align        => 'left',
                text_wrap    => 1);

    #add the formats to the workbook
    $formats{header} = $workbook->add_format(%header);
    $formats{text} = $workbook->add_format(%text);

    return ($workbook,\%formats);

}


###############################################################
sub _printUsage {
    print "Usage: $0.pl <arguments>";
    print qq(
    The list of arguments includes:

    -f|--fasta_file <fasta_file>
        Required.  The file of the sequences to be searched. 

    -m|--masked_file <masked_fasta_file>
        Required.  A soft-masked version of the fasta file (soft masked means low
        complexity sequences are in lower case bases.)

    -p|--project "project name"
        Optional.  A project name for use in the Excel output.

    );
    print "\n";
    return;
}


1;
