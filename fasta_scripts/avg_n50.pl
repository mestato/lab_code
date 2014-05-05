#!/usr/bin/perl
use strict;
my $file = shift;
open IN, $file || die "File not found -3 \n";

my $count = 0;
my $len = 0;
my $total = 0;
my @x;
my $largest = 0;
my $smallest = 1000000000000;

while(<IN>){
 
 	if(/^>(.*)/){
        $count++;
        if($len>0){
            $total+=$len;
            push @x,$len;
            if($len > $largest){
                $largest = $len;
            }
            if($len < $smallest){
                $smallest = $len;
            }
        }
        $len=0;
    }
    else{
        s/\s//g;
        $len+=length($_);
    }
}
if ($len>0){
    $total+=$len;
    push @x,$len;
    if($len > $largest){
        $largest = $len;
    }
    if($len < $smallest){
        $smallest = $len;
    }
}

my $average = $total/$count;
print "The file is $file\n";
print "The number of contigs is $count\n";
print "The total length of the contigs is $total\n";
print "The average length of the contigs is $average\n"; 
print "The largest contig is $largest\n";
print "The smallest contig is $smallest\n";

@x=sort{$b<=>$a} @x; 
my ($count,$half)=(0,0);
for (my $j=0;$j<@x;$j++){
    $count+=$x[$j];
    if (($count>=$total/2)&&($half==0)){
        print "N50: $x[$j]\n";
        $half=$x[$j]
    }elsif ($count>=$total*0.9){
        print "N90: $x[$j]\n";
        exit;
    }
}
 

close IN;
