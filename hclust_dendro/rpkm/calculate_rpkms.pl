#!/usr/bin/perl
use strict;

my %lens;
get_transcript_lengths();

my @files = `ls ./raw_counts/*`;
foreach my $f (@files){
	my $total_reads = get_total_mapped_reads($f);
	#print "$total_reads total_reads\n";

	$f =~ /.*\/(\S+)$/;
	my $file_base = $1;
	$file_base =~ s/.txt//;
	#print $file_base."\n";

	get_fpkms($f, $file_base, $total_reads);
}

#-------------------------------------------
sub get_transcript_lengths{
	open IN, "./Fraxinus_pennsylvanica_120313_transcripts.fasta";
	while(<IN>){
		if(/>(\S+).*len=(\d+)/){
			$lens{$1} = $2;
			#print "$1 $2\n";
		}
	}
	close IN;
}

#-------------------------------------------
sub get_total_mapped_reads{
	my $file = shift;
	my $total = 0;

	open IN, $file;
	while(<IN>){
		/(\d+)$/;
		$total += $1;
	}
	close IN;
	return $total;
}

#-------------------------------------------
sub get_fpkms{
	my $file = shift;
	my $file_base = shift;
	my $total_reads = shift;

	print "$file_base $total_reads\n";
	open IN, $file;
	open OUT, ">$file_base.rpkm";
	while(<IN>){
		/^(\S+)\s+(\d+)$/;
		my $transcript = $1;
		my $count = $2;
		my $transcript_len = $lens{$transcript};
		my $fpkm = ($count*1000000000)/($transcript_len*$total_reads);
		print OUT "$transcript\t$fpkm\n";
	}
	close IN;
	close OUT;
}

