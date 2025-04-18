#! /usr/bin/perl

# 30 nov 2020 mk@mat.ethz.ch

sub USAGE { print<<EOF; 
\nperl $0 [-unfolded] Z1+initconfig.dat > myfile.dump
OR
perl $0 [-unfolded] Z1+SP.dat > myfile.dump
OR
perl $0 [-unfolded] Z1+PPA.dat > myfile.dump

This script converts a Z1+ generated .dat configuration or trajectory file to 
a lammps dump-formatted configuration or trajectory file (myfile.dump). If the option 
-unfolded is entered, the dump-file contains unfolded coordinates (xu yu zu),
otherwise it contains folded (wrapped) coordinates (x,y,z). 

The dump-file contains two bead types: type 1 (interior bead), type 2 (terminal bead).
EOF
exit;
};

sub round { my $TMP="@_"; $TMP+=0; if ($TMP eq 0) { } else { $TMP=($TMP/abs($TMP))*int(abs($TMP)+0.5); }; $TMP; };

foreach $arg (@ARGV) {
    if ($arg eq "-unfolded") { 
        $unfolded=1; 
    } else {
        $file=$arg; if (-s "$file") { } else { print "missing file $file\n"; exit; }; 
    };
};
if (!$file) { USAGE; };


sub strip { chomp $_[0]; $_[0]=~s/^\s+//g; $_[0]=~s/\s+$//; $_[0]=~s/\s+/ /g; $_[0]; }; 

$frame=0;
open(A,"<$file"); 
while (!eof(A)) { 
 $frame+=1; 
 $chains=<A>+0; if ($chains eq 0) { print "This is not a dat-formatted file created by Z1+\n"; exit; }; 
 $line=<A>; $line=strip($line); ($boxx,$boxy,$boxz)=split(/ /,$line); 
 print "ITEM: TIMESTEP\n$frame\n"; 
 $atoms = 0;
 foreach $mol (1 .. $chains) { 
    $N[$mol]=<A>+0; 
    foreach $b (1 .. $N[$mol]) { 
        $atoms+=1; 
        $line=<A>; chomp $line; strip($line); ($x[$atoms],$y[$atoms],$z[$atoms],$rest)=split(/ /,$line);
    };
 };
 print "ITEM: NUMBER OF ATOMS\n$atoms\n";
 print "ITEM: BOX BOUNDS pp pp pp\n";
 $boxxh=$boxx/2; 
 $boxyh=$boxy/2;
 $boxzh=$boxz/2;
 print "-$boxxh $boxxh\n-$boxyh $boxyh\n-$boxzh $boxzh\n";
 if ($unfolded) { 
    print "ITEM: ATOMS id mol type xu yu zu\n";
 } else {
     print "ITEM: ATOMS id mol type x y z\n";
 };
 $id=0;
 foreach $mol (1 .. $chains) { 
  $id+=1;
  if ($unfolded) { 
    print "$id $mol 2 $x[$id] $y[$id] $z[$id]\n";
  } else { 
    $xf=$x[$id]-$boxx*round($x[$id]/$boxx); 
    $yf=$y[$id]-$boxy*round($y[$id]/$boxy);
    $zf=$z[$id]-$boxz*round($z[$id]/$boxz);
    print "$id $mol 2 $xf $yf $zf\n";
  };
  foreach $j (2 .. $N[$mol]-1) { 
    $id+=1;
    if ($unfolded) {
        print "$id $mol 1 $x[$id] $y[$id] $z[$id]\n";
    } else {
        $xf=$x[$id]-$boxx*round($x[$id]/$boxx);
        $yf=$y[$id]-$boxy*round($y[$id]/$boxy);
        $zf=$z[$id]-$boxz*round($z[$id]/$boxz);
        print "$id $mol 1 $xf $yf $zf\n";
    }; 
  };
  $id+=1;
  if ($unfolded) {   
    print "$id $mol 2 $x[$id] $y[$id] $z[$id]\n";
  } else {
        $xf=$x[$id]-$boxx*round($x[$id]/$boxx);
        $yf=$y[$id]-$boxy*round($y[$id]/$boxy);
        $zf=$z[$id]-$boxz*round($z[$id]/$boxz);
        print "$id $mol 2 $xf $yf $zf\n";
  };
 };
}; 
