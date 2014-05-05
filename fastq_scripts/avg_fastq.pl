#!/usr/bin/perl
use strict;
my $file = shift;
open IN, $file || die "File not found -3 \n";

my $count = 0;
my $len = 0;
my $total = 0;
my $smallest = 1000000000000;
my $largest = 0;
my $line;
my $total_ns;
my $total_gcs;

while($line = <IN>){
    $count++;
    my $seq = <IN>;
    chomp $seq;
    $line = <IN>;
    $line = <IN>;

    my $len = length $seq;
    $total+=$len;
    if($len < $smallest){ $smallest = $len;}
    if($len > $largest){ $largest = $len;}

    my $n_cnt = () = $seq =~ /[Nn]/g;
    $total_ns += $n_cnt;

    my $gc_cnt = () = $seq =~ /[GCgc]/g;
    $total_gcs += $gc_cnt;

}

my $average = $total/$count;
print "The file is $file\n";
print "The number of contigs is $count\n";
print "The total length of the contigs is $total\n";
print "The average length of the contigs is $average\n"; 
print "The smallest contig is $smallest\n";
print "The largest contig is $largest\n";

my $gc_pct = $total_gcs/$total;
my $n_pct = $total_ns/$total;
print "Total Ns: $total_ns\n";
print "Percent Ns: $n_pct\n";
print "GC percent: $gc_pct\n";

close IN;
