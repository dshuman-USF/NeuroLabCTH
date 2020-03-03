#!/usr/bin/perl

use vars qw /$scale/;
use File::Basename;
use Cwd;

for $n (0..$#ARGV) {if ($ARGV[$n] eq 'noref'   ) {splice @ARGV, $n, 1; $noref    = 1; last;}}
for $n (0..$#ARGV) {if ($ARGV[$n] eq 'keep_oob') {splice @ARGV, $n, 1; $keep_oob = 1; last;}}
for $n (0..$#ARGV) {if ($ARGV[$n] eq 'triangle') {splice @ARGV, $n, 1; $triangle = 1; last;}}

if (@ARGV >= 3) {
    $blue  = pop @ARGV;
    $green = pop @ARGV;
    $red   = pop @ARGV;
    $rgb = 1;
}
$which = $ARGV[0];
$which =~ s/\.csv$//;
$which =~ s/\.xls$//;

#if ($which) {
#    for (`cat colors`) {
#        s/[\r\n]*$//;
#        ($red, $green, $blue, $name) = split " ",$_,4;
#        last if $name eq $which;
#    }
#    $rgb = 1 if $name eq $which;
#    if ($which eq $ARGV[0] && !-e $which) {
#        splice @ARGV, 0, 1;
#    }
#}

# we only use this for one thing, adjusting the stereotaxic cell coords,
# so assume other scripts are on same path as this script
$ourpath = Cwd::abs_path($0);
$srcpath = File::Basename::dirname($ourpath);
$_  = $srcpath . '/cth_data.pl';
do $_ || die "$_: $!";

$plate_count = @level - 1;

sub distance
{
    my ($p0, $p1) = @_;
    my ($x0, $y0, $z0) = @{$points[$p0]};
    my ($x1, $y1, $z1) = @{$points[$p1]};
    return (sqrt (($x1 - $x0)**2 + ($y1 - $y0)**2 + ($z1 - $z0)**2));
}

$point_count = 0;
for ($pn = 1; $pn <= $plate_count; $pn++) {
    $point_count += @{$pt[$pn][1]} / 2;
}

@points = ();
for ($pn = 1; $pn <= $plate_count; $pn++) {
    $d = sqrt (($scale0[$pn]{x} - $scale1[$pn]{x})**2 + ($scale0[$pn]{y} - $scale1[$pn]{y})**2);
    $d_per_mm = $d / $scale[$pn];
    $pcnt = @{$pt[$pn][1]} / 2;
    for ($n = 0; $n < $pcnt; $n++) {
	my $live_over_berman = 12.5 / 11.61;
	my $x = $pt[$pn][1][2*$n] / $d_per_mm * $live_over_berman;
	my $y = ($pt[$pn][1][2*$n+1] + 742.5) / $d_per_mm * $live_over_berman;
	my $z = -$level[$pn] * $live_over_berman;
	my $a = -16.7 / 180 * 3.141592653589793;
	my $yr = cos ($a) * $y - sin ($a) * $z;
	my $zr = sin ($a) * $y + cos ($a) * $z;
	$y = $yr;
	$z = $zr;
	push @points, [$x, $y, $z];
    }
}

$start1 = 0;
@triangles = ();
for ($pn = 1; $pn < $plate_count; $pn++) {
    $start0 = $start1;
    $pcnt0 = @{$pt[$pn][1]} / 2;
    $start1 = $start0 + $pcnt0;
    next if $pcnt0 == 0;
    $pcnt1 = @{$pt[$pn+1][1]} / 2;
    last if $pcnt1 == 0;
    $first_p0 = $p0 = 0;
    $min = distance ($start0 + p0, $start1);
    $n_at_min = 0;
    for ($n = 1; $n < $pcnt1; $n++) {
	my $d = distance ($start0 + $p0, $start1 + $n);
	if ($d < $min) {
	    $min = $d;
	    $n_at_min = $n;
	}
    }
    $first_p1 = $p1 = $n_at_min;
    $p0cnt = $p1cnt = 0;
    while (1) {
	my $d0 = distance ($start0 + ($p0 + 1) % $pcnt0, $start1 + $p1);
	my $d1 = distance ($start0 + $p0, $start1 + ($p1 + 1) % $pcnt1);

	if ($d0 < $d1) {
	    push @triangles, [$start0 + $p0, $start1 + $p1, $start0 + ($p0 + 1) % $pcnt0];
	    $p0 = ($p0 + 1) % $pcnt0;
	    $p0cnt++;
	}
	else {
	    push @triangles, [$start0 + $p0, $start1 + $p1, $start1 + ($p1 + 1) % $pcnt1];
	    $p1 = ($p1 + 1) % $pcnt1;
	    $p1cnt++;
	}
	if ($p0cnt > 0 && $p0 == $first_p0) {
	    do {
		push @triangles, [$start0 + $p0, $start1 + $p1, $start1 + ($p1 + 1) % $pcnt1];
		$p1 = ($p1 + 1) % $pcnt1;
	    } while ($p1 != $first_p1);
	    last;
	}
	if ($p1cnt > 0 && $p1 == $first_p1) {
	    do {
		push @triangles, [$start0 + $p0, $start1 + $p1, $start0 + ($p0 + 1) % $pcnt0];
		$p0 = ($p0 + 1) % $pcnt0;
	    } while ($p0 != $first_p0);
	    last;
	}
    }
}
if ($triangle) {
    printf "#define TRIANGLE_COUNT %d\n", scalar @triangles;
    print <<EOT;

typedef struct
{
    double x;
    double y;
    double z;
} Point;

typedef struct
{
    Point a;
    Point b;
    Point c;
} Triangle;

EOT
    print "static Triangle triangle[TRIANGLE_COUNT] = {\n";
    for (@triangles) {
        my ($a, $b, $c) = @$_;
        my ($xa, $ya, $za) = @{$points[$a]};
        my ($xb, $yb, $zb) = @{$points[$b]};
        my ($xc, $yc, $zc) = @{$points[$c]};
        print "{{$xa, $ya, $za}, {$xb, $yb, $zb}, {$xc, $yc, $zc}},\n";
    }
    print "};\n";
    exit;
}

sub print_triangle
{
    my $tref = shift;
    my ($a, $b, $c) = @$tref;
    my ($xa, $ya, $za) = @{$points[$a]};
    my ($xb, $yb, $zb) = @{$points[$b]};
    my ($xc, $yc, $zc) = @{$points[$c]};
    printf STDERR "%6.2f\t%6.2f\t%6.2f\n", $xa, $ya, $za;
    printf STDERR "%6.2f\t%6.2f\t%6.2f\n", $xb, $yb, $zb;
    printf STDERR "%6.2f\t%6.2f\t%6.2f\n", $xc, $yc, $zc;
}

# for (@triangles) {
#     my ($a, $b, $c) = @$_;
#     my ($xa, $ya, $za) = @{$points[$a]};
#     my ($xb, $yb, $zb) = @{$points[$b]};
#     my ($xc, $yc, $zc) = @{$points[$c]};
#     printf ("%6.2f\t%6.2f\t%6.2f\n", $xa, $ya, $za);
#     printf ("%6.2f\t%6.2f\t%6.2f\n", $xb, $yb, $zb);
#     printf ("%6.2f\t%6.2f\t%6.2f\n", $xc, $yc, $zc);
#     printf ("\n");
# #    printf ("%6.2f %6.2f %6.2f, %6.2f %6.2f %6.2f, %6.2f %6.2f %6.2f\n", $xa, $ya, $za, $xb, $yb, $zb, $xc, $yc, $zc);
# #    print "$xa $ya $za\n$xb $yb $zb\n$xc $yc $zc\n\n";
# }
# exit;

sub sign
{
    my ($p1, $p2, $x, $z) = @_;
    my ($x1, $y1, $z1) = @{$points[$p1]};
    my ($x2, $y2, $z2) = @{$points[$p2]};
    my $t = ($x2 - $x1) * ($z - $z1) - ($x - $x1) * ($z2 - $z1);
#    printf "point count: %d\n", scalar @points;
#    print "p1, p2: $p1 $p2\n";
#    print "my $t = ($x2 - $x1) * ($z - $z1) - ($x - $x1) * ($z2 - $z1);\n";
#    print "t: $t\n";
    return $t < 0 ? -1 : ($t > 0 ? 1 : 0);
}

sub in_triangle
{
    my ($tref, $x, $z) = @_;
    my ($p0, $p1, $p2) = @$tref;
    my $sign0 = sign ($p0, $p1, $x, $z);
    my $sign1 = sign ($p1, $p2, $x, $z);
    my $sign2 = sign ($p2, $p0, $x, $z);
#    print "signs: $sign0 $sign1 $sign2\n";
    return  (($sign0 == $sign1 && $sign1 == $sign2)
	      || ($sign0 == 0 && $sign1 == $sign2)
	      || ($sign1 == 0 && $sign0 == $sign2)
	      || ($sign2 == 0 && $sign0 == $sign1))
}

sub y_at_x_z
{
    my ($tref, $x, $z) = @_;
    my ($p1, $p2, $p3) = @$tref;
    my ($x1, $y1, $z1) = @{$points[$p1]};
    my ($x2, $y2, $z2) = @{$points[$p2]};
    my ($x3, $y3, $z3) = @{$points[$p3]};

    $D = ($x1*($z3-$z2)-$x2*$z3+$x3*$z2+($x2-$x3)*$z1);
    $A = ($y1*($z3-$z2)-$y2*$z3+$y3*$z2+($y2-$y3)*$z1) / $D;
    $B = ($x1*($y3-$y2)-$x2*$y3+$x3*$y2+($x2-$x3)*$y1) / $D;
    $C = -($x1*($y3*$z2-$y2*$z3)+$y1*($x2*$z3-$x3*$z2)+($x3*$y2-$x2*$y3)*$z1) / $D;

    if ($D == 0) { return undef; }
    $y = $A * $x + $B * $z + $C;
#    print "$x $z is in ($x1, $y1, $z1) ($x2, $y2, $z2) ($x3, $y3, $z3) at $y\n";
    return $y;
}

sub surface
{
    my @yvals = ();
    for (@triangles) {
	if (in_triangle ($_, @_)) {
#	    print STDERR "in triangle\n";
	    my $y = y_at_x_z ($_, @_);
	    if (defined $y) {
		push @yvals, $y;
	    }
	}
    }
    push @yvals, ($yvals[0] > -2 ? -7.35 : 1.86) if @yvals == 1;
#    print STDERR "surface(@_) = @yvals\n";
    return undef if @yvals == 0;
    @yvals = sort {$a <=> $b} @yvals;
    return [$yvals[$#yvals], $yvals[0]];
}

if (@ARGV == 0 || `file "$ARGV[0]"` =~ /ascii|text/i) {@txt = <>;}
else {
   $cvtcmd = "$srcpath" . "/cth_xls2csv.pl " . "$ARGV[0]";
   @txt = `$cvtcmd`;
}

%map = (
        CName          => "name",
        "Cell Name"    => "name",
        Cell_Name      => "name",
        Cell           => "name",
        AP             => "ap",
        "AP coord"     => "ap",
        AP_coord       => "ap",
        RL             => "rl",
        "RL coord"     => "rl",
        RL_coord       => "rl",
        Depth          => "dp",
        "Depth coord"  => "dp",
        Depth_coord    => "dp",
        "UNIT NAME"    => "name",
        "MERGED CH. #" => "mchan",
        "A/P"          => "ap",
        "R/L"          => "rl",
        DEPTH          => "dp"
        );


$_ = shift @txt;
$orig = $_;
s/\s+$//;

$cellfile = 1;
if (/^EXPERIMENT DATE:,,,.*,,EXPERIMENT:,,,.*,RECORDING:,,,/) {
    $cellfile = 0;
    $_ = shift @txt; #blank line
    $_ = shift @txt; #AA/STA line
    $_ = shift @txt; #header
}
print STDERR $cellfile ? "cell file\n" : "info file\n";

@vars = split /,/;

for (@vars) {
#    print STDERR "map $_ to $map{$_}\n" if $map{$_};
    $_  = $map{$_} if $map{$_};
}

for (@vars) {$n{$_}++;}

if ($n{'ap'} + $n{'rl'} + $n{'dp'} + $n{'name'} == 0) { 
    my $n = 0;
    my $p = "^-?1?\\d\\.\\d\\d\$";
    $orig = $_;
    for (@vars) {
        if (/^[a-z][a-z]?[0-9]+[a-z]?$/i) {
            print STDERR "name: $_\n";
            $_ = "name";
        }
        elsif (/$p/ && $vars[$n + 1] =~ /$p/ && $vars[$n + 2] =~ /$p/ && $vars[$n + 3] !~ /$p/) {
            print STDERR "ap: $_\n";
            $_ = 'ap';
        }
        elsif ($vars[$n - 1] eq 'ap') {
            print STDERR "rl: $_\n";
            $_ = 'rl';
        }
        elsif ($vars[$n - 1] eq 'rl') {
            print STDERR "dp: $_\n";
            $_ = 'dp';
        }
        $n++;
    }
    %n = ();
    for (@vars) {$n{$_}++;}
    unshift @txt, $orig;
}

unless ($n{'ap'} == 1 && $n{'rl'} == 1 && $n{'dp'} == 1 && $n{'name'} == 1) {
    die "can't parse header";
}

@okvar = (qw/name mchan ap rl dp ref dchan r g b/); # dale mod
@okvar{@okvar} = @okvar;
for (@vars) {$_ = "" unless $okvar{$_};}

for (@txt) {
    s/[\r\n]*$//;
    s/\.,/,/g;

    {
	my ($ltr);
	undef $ref;
	undef $dchan;
        my @vals = split /,/;
        for $n (0..$#vars) {
            if ($vars[$n] =~ /\S/) {
                $ {$vars[$n]} = $vals[$n] if $vars[$n] ne "";
            }
        }
	if (($ap !~ /\d/) || ($rl !~ /\d/) || ($dp !~ /\d/) || ($name !~ /[a-z]/i)) {
	    print STDERR "$_: \"$name\", \"$ap\", \"$rl\", \"$dp\"\n";
	    next;
	}

        if ($rgb) {
            $r = $red;
            $g = $green;
            $b = $blue;
        }
        else {
            $r = 1 unless defined $r;
            $g = 1 unless defined $g;
            $b = 1 unless defined $b;
        }

	if ($name =~ /^([pq])/i) {
            $ltr = $1;
	    $ap += 12.5;
	}	    

	my $x =  $rl;
	my $y = -$dp;
	my $z = -$ap;

        my $xz = "$x, $z";
        $surface{$xz} = surface ($x, $z) unless $surface{$xz};
        my ($top, $bottom) = @{$surface{$xz}};

	if (!defined $ref || !defined $dchan || $ref == $dchan) {
	    $ref = $xz unless defined $ref;
	    if (!$dy{$ref} || $name eq $ltr . '0') {
		if (!defined ($dy{$ref} = $top)) {
                    printf STDERR "no depth reference for cell at $ap,$rl,$dp, ref = $ref\n";
                    $dy{$ref} = 0;
                }
	    }
	}
        # reference (pseudo) cells end in a non-digit followed by a 0
        # unless the user said "noref" and it is a ref cell, use the cell
        unless ($noref && $name =~ /\D0$/) {
            push @ref, $ref;
            if ($name =~ /\D0$/) { # reference cell
                push @tcells, [$x, 0, $z]; # ignore depth
            }
            else {
                push @tcells, [$x, $y, $z,$mchan];   # dale mod
            }
            push @tcolors, [$r, $g, $b];
            push @twhere, "$ap,$rl,$dp,$mchan";      # dale mod
        }
	if (!defined $dy{$ref}) {
	    print STDERR "dy is not defined for ref: $ref\n";
	    die;
	}
    }
	
    printf STDERR ("%d\n", $txtidx) if ++$txtidx % 100 == 0;
}
@ref == @tcells || die;
for $n (0..$#tcells) {
    if (!defined $dy{$ref[$n]}) {
	print STDERR "dy is not defined for $ref[$n]\n";
	die;
    }
    $tcells[$n][1] += $dy{$ref[$n]};
    my ($x, $y, $z) = @{$tcells[$n]};
    my ($top, $bottom) = @{$surface{"$x, $z"}};
    if ((defined $top && $y < $top && $y > $bottom) || $keep_oob) {
        push @cells, $tcells[$n];
        push @colors, $tcolors[$n];
    } else {
        print STDERR "dropping out of bounds cell at $twhere[$n]\n";
    }
}

printf "object 1 class array type float rank 1 shape 3 items %d data follows\n", scalar @cells;

for (@cells) {
    print "@$_\n";
}

printf ("object 2 class constantarray type float rank 0 items %d data follows\n1\n", scalar @cells);

printf "object 3 class array type float rank 1 shape 3 items %d data follows\n", scalar @cells;

for (@colors) {
    print "@$_\n";
}
printf ("attribute \"dep\" string \"positions\"\n");

print "object \"cells\" class field\n",
    "component \"positions\" value 1\n",
    "component \"data\" value 2\n",
    "component \"colors\" value 3\n",
    "end\n";
