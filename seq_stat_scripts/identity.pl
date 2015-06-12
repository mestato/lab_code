#!/data/apps/perl/5.16.2/bin/perl

use strict;
use Bio::SeqIO;
BEGIN { $ENV{CLUSTALDIR} = '/lustre/projects/staton/software/clustalw-2.1-linux-x86_64-libcppstatic/' }
use Bio::Tools::Run::Alignment::Clustalw;
use Bio::SimpleAlign;
use Bio::AlignIO;
 
# Check for input
if($ARGV[0] eq ""){
	die "\n\t Usage : perl <thisScript.pl> <file.fasta>\n\n";
}

# Capture input in variable, initialize variables
my $fasta = $ARGV[0];
my %sequences;
my $seqio;
my @seq_array;

# load fasta file into BioPerl seqio object
$seqio = Bio::SeqIO->new(-file => $fasta, -format => "fasta");
while(my$seqobj = $seqio->next_seq) {
	my $id  = $seqobj->display_id;
	my $seq = $seqobj->seq;
	$sequences{$id} = $seqobj;
	# add seqio object to an array
	push(@seq_array, $sequences{$id});
}

# run clustalw
my $factory = Bio::Tools::Run::Alignment::Clustalw->new(-matrix => 'BLOSUM', -quiet => '1');
my $ktuple = 3;
$factory->ktuple($ktuple);  # change the parameter before executing
my $seq_array_ref = \@seq_array;
my $aln = $factory->align($seq_array_ref);

# print results
print $seq_array[0]->id();
print "\t";
print $seq_array[1]->id();
print "\t";
print $aln->percentage_identity();
print "\n";
