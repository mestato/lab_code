#!/usr/local/bin/perl


$file = shift;


 open IN, $file || die "File not found -3 \n";
 
 $num = 0;
 $total = 0;
 $sequence = "";
 
while($lineIn = <IN>){
 
 	if($lineIn =~ /^>(.*)/){
		##This is the beginning of a sequence
		
		if(length $sequence > 0){
			$total += length $sequence;
			$num += 1;
		}
		
		$sequence = "";
 	}
 	
 	else{
 		##This is part of a sequence

 		$sequence .= $lineIn;
 		chomp($sequence);
 	}
}
if(length $sequence > 0){
		$total += length $sequence;
		$num += 1;
}

$average = $total/$num;
print "The file is $file\n";
print "The number of contigs is $num\n";
print "The total length of the contigs is $total\n";
print "The average length of the contigs is $average\n"; 
close IN;
