use strict;

my @files = `ls *rpkm`;
foreach my $f (@files){
	my $total = 0;
	open IN, $f;
	while(<IN>){
		/(\S+)\t([0-9.]+)$/;
		#print "$1 $2\n";
		my $cnt = $2;
		if($cnt > 0.1){ $total++; }
	}
	close IN;

	chomp $f;
	print "$f $total\n";

}
