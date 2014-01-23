#!/usr/bin/env perl
################################################################################
# Author: Meg Staton 
# Date: Jan 23rd, 2014
# Version: 1
#
# DESCRIPTION
# -----------
# This script identifies forward and reverse reads with the same beginning bases.
#
# Dependencies:
# ------------
# Perl must have access to the package Getopt::Long, available from 
# CPAN.
#
# Usage:
# -----
# Usage: find_dup_reads.pl -f <forward fastq file> -r <reverse fastq file>
#
# The list of arguments includes:
# 
# -f|--forward_fastq_file <fastq_file>
# Required.  The file of the forward direction sequences to be searched.  
# 
# -r|--reverse_fastq_file <fastq_file>
# Required.  The file of the reverse direction sequences to be searched.  
# 
# Output:
# ------
# Four output files are produced, in the same directory as the input files:
#
# <forward_fastq_file>.uniq
# A fastq file with forward sequences that do not share duplicated beginning
# with the corresponding reverse read.
#
# <forward_fastq_file>.dup
# A fastq file with forward sequences that do share duplicated beginning
# with the corresponding reverse read.
#
# <reverse_fastq_file>.uniq
# A fastq file with reverse sequences that do not share duplicated beginning
# with the corresponding forward read.
#
# <reverse_fastq_file>.dup
# A fastq file with reverse sequences that do share duplicated beginning
# with the corresponding forward read.
#
# Details: 
# -------
# This script compares from base 0 to base 20 of the forward read against
# base 0 to base 20 of the reverse read.  In order to detect duplications
# where a base has been miscalled, it also compares subsequences from base
# 10 to base 30 and from base 20 to base 40.  If any of these subsequnces
# from both the forward and the reverse read is identical it kicks the reads
# out into the duplicate file output.  All other reads go into the unique
# file output.
#
# The length of the subsequence samples, currently set at 20, can be modified
# with the $LEN variable below.
#
# WARNING: This script puts the entire forward file into memory in a hash, so
# use with caution on large files.
#
# Needed improvements and caveat emptors: 
# - if the duplicated region is shifted by even one base in either
#   direction, this script will not detect it.
#
# - memory requirements are ridiculous, a better strategy is needed
#
# - the output files are put whereever the input files are because the path
#   is not parsed off, could be fixed with File::Basename
# 
#
#---------------------------------------------------------------------------------
use strict;

#-------------------------------------------------------------------------------
#  DEPENDENCIES

use Getopt::Long;

#-------------------------------------------------------------------------------
# GLOBAL PARAMETERS
#-------------------------------------------------------------------------------
my $LEN = 20;

my $for_file; 
my $rev_file;

Getopt::Long::Configure ('bundling');
GetOptions('f|forward=s'  => \$for_file,
           'r|reverse=s'  => \$rev_file);

## Check that all required parameters have been included
if(!$for_file){ print "Forward fastq file is required.\n"; _printUsage(); exit;}
if(!$rev_file){ print "Reverse fastq file is required.\n"; _printUsage(); exit;}

## Check that fastq file exists
if(! -e $for_file) { print "Forward fastq file $for_file does not exist\n"; exit; }
if(! -e $rev_file) { print "Reverse fastq file $rev_file does not exist\n"; exit; }


my %seq;
my $line;
open IN,  $for_file;
while($line = <IN>){
    $line =~ /\@(\S+)/;
    my $name = $1;
	# store the first of the four fastq lines
    $seq{$name}{FIRST} = $line;
    my $seq = <IN>;

	# store the second of the four fastq lines
    $seq{$name}{SEQ} = $seq;

	# store the third of the four fastq lines
    $seq{$name}{THIRD} = <IN>;

	# store the fourth of the four fastq lines
    $seq{$name}{FOURTH} = <IN>;
}
close IN;


my $count = 0;
my $dups = 0;

open IN,  $rev_file;
open OUT_1_OK, ">".$for_file.".uniq";
open OUT_2_OK, ">".$rev_file.".uniq";
open OUT_1_DUP, ">".$for_file.".dup";
open OUT_2_DUP, ">".$rev_file.".dup";

while($line = <IN>){
    $count++;
    $line =~ /\@(\S+)/;
    my $name = $1;

	# get the rest of the lines associated with this fastq record
    my $seq = <IN>;
    my $third = <IN>;
    my $fourth = <IN>;
	
	# get a substring of the sequence characters starting from base 0
	# length determined by global var $LEN
    my $rev_begin = substr($seq, 0, $LEN);

	# get a substring of the sequence characters starting from base 10
	# length determined by global var $LEN
    my $rev_begin2 = substr($seq, 10, $LEN);

	# get a substring of the sequence characters starting from base 20
	# length determined by global var $LEN
    my $rev_begin3 = substr($seq, 20, $LEN);

    if(!exists $seq{$name}){
        print "Can't find $name forward seq\n";
    }
    else{
		# get the same three substrings from the forward read
        my $for_begin = substr($seq{$name}{SEQ}, 0, $LEN);
        my $for_begin2 = substr($seq{$name}{SEQ}, 10, $LEN);
        my $for_begin3 = substr($seq{$name}{SEQ}, 20, $LEN);

		## if any of the substrings between the forward and reverse sequences match,
		## its a duplicate that needs to be filtered.
        if(($for_begin eq $rev_begin) || ($for_begin2 eq $rev_begin2) || ($for_begin3 eq $rev_begin3)){
            ## found duplicated bases at begining of r1 and r2
            $dups++;
            print OUT_1_DUP $seq{$name}{FIRST};
            print OUT_1_DUP $seq{$name}{SEQ};
            print OUT_1_DUP $seq{$name}{THIRD};
            print OUT_1_DUP $seq{$name}{FOURTH};

            print OUT_2_DUP $line;
            print OUT_2_DUP $seq;
            print OUT_2_DUP $third;
            print OUT_2_DUP $fourth;

        }
        else{
            ## r1 and r2 look uniq
            print OUT_1_OK $seq{$name}{FIRST};
            print OUT_1_OK $seq{$name}{SEQ};
            print OUT_1_OK $seq{$name}{THIRD};
            print OUT_1_OK $seq{$name}{FOURTH};

            print OUT_2_OK $line;
            print OUT_2_OK $seq;
            print OUT_2_OK $third;
            print OUT_2_OK $fourth;

        }
    }
}
close IN;


my $pct = ($dups/$count)*100;
my $uniq2 = $count - $dups;
print "total pairs examined: $count\n";
print "pairs with dup beginning: $dups\n";
print "pairs with uniq beginning: $uniq2\n";
print "percent duplicates: $pct %\n";


################################################################################
sub _printUsage {
    print qq(
Usage: $0 -f <forward fastq file> -r <reverse fastq file>

The list of arguments includes:
 
-f|--forward
Required.  The file of the forward direction sequences to be searched.  
 
-r|--reverse
Required.  The file of the reverse direction sequences to be searched.  
 
WARNING: This script puts the entire forward fastq file into memory, so
use with extreme caution on large files.

    );
    print "\n";
    return;
}
