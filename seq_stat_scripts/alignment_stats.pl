#!/data/apps/perl/5.16.2/bin/perl

use strict;
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;
 
# Check for input
if($ARGV[0] eq ""){
	die "\n\t Usage : perl <thisScript.pl> <file.fasta>\n\n";
}

# Capture input in variable, initialize variables
my $fasta = $ARGV[0];
my %sequences;
my $seqio;

# load fasta file into BioPerl seqio object
$seqio = Bio::SeqIO->new(-file => $fasta, -format => "fasta");
while(my$seqobj = $seqio->next_seq) {
	my $id  = $seqobj->display_id;
	my $seq = $seqobj->seq;
	$sequences{$id} = $seqobj;
}

#my $factory = Bio::Tools::Run::StandAloneBlast->new(-outfile => 'bl2seq.out');
#my $bl2seq_report = $factory->bl2seq($sequences{"Liriodendron_tulipifera_10132014_comp2_c0_seq1"}, $sequences{"Liriodendron_tulipifera_10132014_comp0_c0_seq1"});
