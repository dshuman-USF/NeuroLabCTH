#!/usr/bin/perl

if (`file $ARGV[0]` =~ /ascii|text/i) {
    @t = <>;
    print STDERR "csv";
}
else {
    @t = `xls2csv.pl $ARGV[0]`;
    print STDERR "xls";
}

for (@t) {
    s/\"//g;
    if (/^([A-Z])\d+[a-z]?,\d+,([-\.0-9]+),([-\.0-9]+),([-\.0-9]+)/) {
	($ltr, $ap,$rl,$dp) = ($1,$2,$3,$4);

	if ($ltr eq 'P') {
	    $ap += 12.5;
	    $dp += 2.00;
	}	    

	$x =  $rl;
	$y = -$dp;
	$z = -$ap;
	
	push @cells, [$x, $y, $z];

    }
}

printf "object 1 class array type float rank 1 shape 3 items %d data follows\n", scalar @cells;

for (@cells) {
    print "@$_\n";
}

printf ("object 2 class constantarray type float rank 0 items %d data follows\n1\n", scalar @cells);

print "object \"cells\" class field\n",
    "component \"positions\" value 1\n",
    "component \"data\" value 2\n",
    "end\n";
