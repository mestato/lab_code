#!/usr/bin/env perl
################################################################################
# Author: Meg Staton & Stephen Ficklin
# Date: Jan 14th, 2014
# Version: 1
#
# DESCRIPTION
# -----------
# This script identifies simple sequence repeats (SSRs) from
# the sequences in a fastq formatted file.
#
# Dependencies:
# ------------
# Perl must have access to the package Getopt::Long, available from 
# CPAN.
# 
# Usage:
# -----
# Usage: findSSRs.pl <arguments>
#
# The list of arguments includes:
# 
# -f|--fastq_file <fastq_file>
# Required.  The file of the sequences to be searched.  
# 
# -a|--fasta_format 
# Optional.  The two output sequence files will be in fasta format 
# instead of fastq format.
#
# Output:
# ------
# Four output files are produced:
#
# <input-file-name>.ssr.fastq (or fasta)
# A fastq (or fasta file if -a is specified) with sequences with a 
# single identified SSR.
#
# <input-file-name>.multissr.fastq (or fasta)
# A fastq (or fasta file if -a is specified) with sequences with more 
# than one identified SSR.
#
# <input-file-name>.ssr_stats.txt
# A text file of statistics about the SSRs discovered.
#
# <input-file-name>.ssr_report.txt
# A tab-delimited file with each SSR.  The columns are sequence name, 
# motif, number of repeats, start position and end position.
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

use Getopt::Long;

#-------------------------------------------------------------------------------
# GLOBAL PARAMETERS
#-------------------------------------------------------------------------------

# REPEAT IDENTIFICATION PARAMETERS
# Specify Motif Frequency  

# Motifs that occur less frequently than indicated below will be ignored.
our $MIN_REPS_2bp = 8;
our $MIN_REPS_3bp = 7;
our $MIN_REPS_4bp = 6;

# Motifs that occur more frequently than indicated below will be ignored.
our $MAX_REPS_2bp = 40;
our $MAX_REPS_3bp = 30;
our $MAX_REPS_4bp = 20;

# SSRs at the beginning or end of a sequence prevent flanking primer design.
# This is how close we will allow an SSR to be to the ends of the sequence.
our $LENGTH_FROM_END = 15;

#-------------------------------------------------------------------------------
# GLOBAL HASHES
#-------------------------------------------------------------------------------

# For counting the number of times each motif is found
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

# For counting the number of times each motif length is found
my %MOTIFLEN = ('2'  => 0,
                '3'  => 0,
                '4'  => 0);
                
    
# Set up the Motif specifications, based on the chosen motif types:.
my @MOTIF_SPECS;
push(@MOTIF_SPECS,[2, $MIN_REPS_2bp, $MAX_REPS_2bp, 'dinucleotides']);
push(@MOTIF_SPECS,[3, $MIN_REPS_3bp, $MAX_REPS_3bp, 'trinucleotides']);
push(@MOTIF_SPECS,[4, $MIN_REPS_4bp, $MAX_REPS_4bp, 'tetranucleotides']);
 
my $SSR_COUNT = 0;  	# total number of SSRs found in sequences
my %STATS = ();			# stores information about each SSR (motif, motif length, etc)
my %SEQ_SSR_LOC = ();	# stores the start positions of SSRs found in each sequence
my $SEQ_COUNT = 0;		# counts the number of sequences seen


#-------------------------------------------------------------------------------
#  EXECUTE
#-------------------------------------------------------------------------------
main();
#-------------------------------------------------------------------------------
#  PUBLIC SUBROUTINES
#-------------------------------------------------------------------------------
# Function Name:  main()

sub main{
    my $fastq_file;		# holds the command line parameter for fastq file name
    my $fasta_flag = 0;		# whether to print fastq or fasta files, set on command line

    Getopt::Long::Configure ('bundling');
    GetOptions('f|fastq_file=s'  => \$fastq_file,
               'a|fasta_format'  => \$fasta_flag);

    ## Check that all required parameters have been included
    if(!$fastq_file){ print "A fastq file is required.\n"; _printUsage(); exit;}

    ## Check that fastq file exists
    if(! -e $fastq_file) { print "Fasta file $fastq_file does not exist\n"; exit; }


    print "finding SSRs...\n";
    my $ssr_out         = "$fastq_file.ssr_report.txt"; 	# the tab-delimited output file
    getContigHash($fastq_file, $ssr_out);

    print "stats...\n";
    my $stats_out       = "$fastq_file.ssr_stats.txt";  	# file of statistics
    printStats($stats_out);
     
    if($fasta_flag){
    	my $fasta_out       = "$fastq_file.ssr.fasta";  	# fasta file of sequences with 1 SSR
    	my $fasta_out_multi = "$fastq_file.multissr.fasta"; # fasta file of sequences with >1 SSR
    	print "printing fasta file...";
    	printFasta($fastq_file, $fasta_out, $fasta_out_multi);
    	print "done.\n";
    }
    else{
    	my $fastq_out       = "$fastq_file.ssr.fastq";  	# fastq file of sequences with 1 SSR
    	my $fastq_out_multi = "$fastq_file.multissr.fastq"; # fastq file of sequences with >1 SSR
    	print "printing fastq file...";
    	printFastq($fastq_file, $fastq_out, $fastq_out_multi);
    	print "done.\n";
    }
     
}

###############################################################

sub getContigHash{
    my $fastq_file  = $_[0]; # file name
    my $ssr_out     = $_[1]; # file name

    open (OUT, ">".$ssr_out) || die "ERROR cannot open $_\n";
    my $out_fh = *OUT;

	open (IN, $fastq_file) || die "Can't open $fastq_file\n";
	my @aux = undef;
	my ($seqname, $seqstr, $qual);
	while (($seqname, $seqstr, $qual) = readfq(\*IN, \@aux)) {
        $SEQ_COUNT++;
        process_seq($seqname, $seqstr, $out_fh);
    }

}

###############################################################

sub process_seq{
    my $contig_name = shift;
    my $seq = shift;
    my $out_fh = shift;


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
        my $regex = "(([gatc]{$motifLength})\\2{$min_number_of_repeats,})";

        ## LOOPC
        # run through the sequence and check for this motif spec 
        while ($seq =~ /$regex/ig) {
            # Get the ssr and motif that were found
            my $ssr = $1;
            my $motif = lc $2;

            ## initialize the repeats varaible to (motifLength - 1)
            ## or 1 if motif length - 1 equals 0;
            my $repeats = $motifLength - 1;
            $repeats = 1 if(!$repeats);

            my $noRepeats = ((length $ssr)/$motifLength);
            my $ssrEnd = pos $seq;
            my $ssrStart = ($ssrEnd - (length $ssr) + 1);

            ## LOOP D - check to make sure the motif isn't the same letter over and over again AND we don't have too many repeats

            if (!($motif =~ /([gatc])\1{$repeats}/) && $noRepeats <= $max_number_of_repeats) {

                ## LOOPE
                # Only store the information if we have never
                # seen this starting position before.
                if (!exists $seen{$contig_name."_ssr".$ssrStart}) {

                    # LOOP F
                    # Check to make sure the SSR is not within $LENGTH_FROM_END bp of either end of the sequence
                    my $seqLen = length $seq;
                    if($ssrStart >= $LENGTH_FROM_END && $ssrEnd <= ($seqLen-$LENGTH_FROM_END)){

                        ## FOUND A SSR TO REPORT
                        my $ssr_id = $contig_name."_ssr".$ssrStart;
                        $seen{$ssr_id} = 1;
    
                        if(exists $SEQ_SSR_LOC{$contig_name}){
                            push @{ $SEQ_SSR_LOC{$contig_name} }, $ssrStart;
                        }
                        else{   
                           $SEQ_SSR_LOC{$contig_name} = [$ssrStart];
                        }
    
                        $SSR_COUNT++;
                        printf $out_fh ("$contig_name\t$motif\t$noRepeats\t$ssrStart\t$ssrEnd\n");
                        # Store STATS
                        $STATS{$ssr_id}{MOTIF}        = $motif;
                        $STATS{$ssr_id}{START}        = $ssrStart;
                        $STATS{$ssr_id}{END}          = $ssrEnd;
                        $STATS{$ssr_id}{MOTIF_LENGTH} = $motifLength;
                        $STATS{$ssr_id}{NO_REPEATS}   = $noRepeats;
    
                        my $motiflen = length $motif;
                        my $tmp = $MOTIFLEN{$motiflen};
                        $tmp++;
                        $MOTIFLEN{$motiflen} = $tmp;
    
                        # Increment motif count
                        foreach my $group (keys %MOTIFS) {
                            my $motifUC = uc($motif);
                            if($group =~ /\|$motifUC\|/){
                                my $tmp = $MOTIFS{$group}++;
                                $tmp++;
                                $MOTIFS{$group} = $tmp;
                            }
                        }# end foreach $group
                    } ## LOOP F

                ## LOOPE
               }

            # LOOPD
            }

        ## LOOPC
        } # end while (sequence =~ /$regex/ig);

    ## LOOPB
    } # end for $index (0 .. (scalar @{$MOTIF_SPECS} - 1)) 

}

################################################################
sub readfq {
        my ($fh, $aux) = @_;
        @$aux = [undef, 0] if (!(@$aux));
        return if ($aux->[1]);
        if (!($aux->[0])) {
                while (<$fh>) {
                        chomp;
                        if (substr($_, 0, 1) eq '>' || substr($_, 0, 1) eq '@') {
                                $aux->[0] = $_;
                                last;
                        }
                }
                if (!($aux->[0])) {
                        $aux->[1] = 1;
                        return;
                }
        }
        my $name = /^.(\S+)/? $1 : '';
        my $seq = '';
        my $c;
        $aux->[0] = undef;
        while (<$fh>) {
                chomp;
                $c = substr($_, 0, 1);
                last if ($c eq '>' || $c eq '@' || $c eq '+');
                $seq .= $_;
        }
        $aux->[0] = $_;
        $aux->[1] = 1 if (!($aux->[0]));
        return ($name, $seq) if ($c ne '+');
        my $qual = '';
        while (<$fh>) {
                chomp;
                $qual .= $_;
                if (length($qual) >= length($seq)) {
                        $aux->[0] = undef;
                        return ($name, $seq, $qual);
                }
        }
        $aux->[1] = 1;
        return ($name, $seq);
}


################################################################
sub printStats{
    my $stats_out = $_[0]; # file name

    open (OUTS, ">".$stats_out) || die "ERROR cannot open $stats_out\n";

    my $time = scalar localtime;                   # Get the current time

    ## count number of seqs with single or multiple SSRs
    my $SINGLE_SSR_COUNT = 0;
    my $MULTI_SSR_COUNT = 0;

    foreach my $contig_name (keys %SEQ_SSR_LOC){
        my $starts = scalar @{ $SEQ_SSR_LOC{$contig_name} };
        if($starts == 1){ $SINGLE_SSR_COUNT++;}
        elsif($starts > 1){ $MULTI_SSR_COUNT++; }
    }
    
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

    close OUTS;

}

sub printFasta{
    my $fastq_file      = $_[0]; # file name
    my $fasta_out       = $_[1]; # file name
    my $fasta_out_multi = $_[2]; # file name

    my $line;
    my $name;

    open (IN, $fastq_file) || die "ERROR cannot open $_\n";
    open (OUTF, ">".$fasta_out) || die "ERROR cannot open $_\n";
    open (OUTFM, ">".$fasta_out_multi) || die "ERROR cannot open $_\n";

    while($line = <IN>){
        if($line =~ /^@(\S+)/){
            $name = $1;
            if(exists $SEQ_SSR_LOC{$name}){
                my $ssr_count = scalar @{ $SEQ_SSR_LOC{$name} };
                if($ssr_count == 1){
                    print OUTF ">$name\n";
                    my $seq = <IN>;
		    print OUTF "$seq\n";
                    $line = <IN>;
                    $line = <IN>;
                }
                else{
                    print OUTFM ">$name\n";
                    my $seq = <IN>;
		    print OUTFM "$seq\n";
                    $line = <IN>;
                    $line = <IN>;
                }
            }
            else{
                $line = <IN>;
                $line = <IN>;
                $line = <IN>;
            }
        }
    } 
        
    close IN;
    close OUTF;
    close OUTFM;

}

sub printFastq{
    my $fastq_file      = $_[0]; # file name
    my $fastq_out       = $_[1]; # file name
    my $fastq_out_multi = $_[2]; # file name

    my $line;
    my $name;

    open (IN, $fastq_file) || die "ERROR cannot open $_\n";
    open (OUTF, ">".$fastq_out) || die "ERROR cannot open $_\n";
    open (OUTFM, ">".$fastq_out_multi) || die "ERROR cannot open $_\n";

    while($line = <IN>){
        if($line =~ /^@(\S+)/){
            $name = $1;
            if(exists $SEQ_SSR_LOC{$name}){
                my $ssr_count = scalar @{ $SEQ_SSR_LOC{$name} };
                if($ssr_count == 1){
                    print OUTF $line;
                    $line = <IN>;
                    print OUTF $line;
                    $line = <IN>;
                    print OUTF $line;
                    $line = <IN>;
                    print OUTF $line;
                }
                else{
                    print OUTFM $line;
                    $line = <IN>;
                    print OUTFM $line;
                    $line = <IN>;
                    print OUTFM $line;
                    $line = <IN>;
                    print OUTFM $line;
                }
            }
            else{
                $line = <IN>;
                $line = <IN>;
                $line = <IN>;
            }
        }
    } 
        
    close IN;
    close OUTF;
    close OUTFM;

}

sub _printUsage {
    print "Usage: $0 <arguments>";
    print qq(
    The list of arguments includes:

    -f|--fastq_file <fastq_file>
        Required.  The file of the sequences to be searched.  

    -a|--fasta_format 
	Optional.  The two output sequence files will be in fasta format
        instead of fastq format.

    );
    print "\n";
    return;
}


1;
