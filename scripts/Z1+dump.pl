#! /usr/bin/perl

# 30 nov 2020 mk@mat.ethz.ch

sub USAGE { print<<EOF; 
\nperl $0 [-unfolded] <Z1-formatted-file> 

This script converts a Z1-formatted configuration or trajectory file to 
a lammps dump-formatted configuration or trajectory file. If the option 
-unfolded is entered, the dump-file contains unfolded coordinates (xu yu zu),
otherwise it contains folded (wrapped) coordinates (x,y,z). 

The dump-file contains two bead types: type 1 (interior beads), type 2 (terminal beads).
EOF
exit;
};

sub round { $_[0]="@_"; $_[0]=$_[0]+0; if ($_[0] eq 0) { } else { $_[0]=($_[0]/abs($_[0]))*int(abs($_[0])+0.5); }; $_[0]; }

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
 $chains=<A>+0; if ($chains eq 0) { print "This is not a Z1-formatted file\n"; exit; }; 
 $line=<A>; $line=strip($line); ($boxx,$boxy,$boxz)=split(/ /,$line); 
 print "ITEM: TIMESTEP\n$frame\n"; 
 $line=<A>; $line=strip($line); @fmt=split(/\*/,$line); 
 if ($#fmt eq 1) { foreach $i (1 .. $chains) { $N[$i-1]=$fmt[1]; }; } else { @N=split(/ /,$line); }; 
 $atoms=0; foreach $i (1 .. $chains) { $atoms+=$N[$i-1]; }; 
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
  $line=<A>; $line=strip($line); ($x,$y,$z)=split(/ /,$line);
  $id+=1;
  if ($unfolded) { 
    print "$id $mol 2 $x $y $z\n";
  } else { 
    $xf=$x-$boxx*round($x/$boxx); 
    $yf=$y-$boxy*round($y/$boxy);
    $zf=$z-$boxz*round($z/$boxz);
    print "$id $mol 2 $xf $yf $zf\n";
  };
  foreach $j (2 .. $N[$mol-1]-1) { 
    $line=<A>; $line=strip($line); ($x,$y,$z)=split(/ /,$line); 
    $id+=1;
    if ($unfolded) {
        print "$id $mol 1 $x $y $z\n";
    } else {
        $xf=$x-$boxx*round($x/$boxx);
        $yf=$y-$boxy*round($y/$boxy);
        $zf=$z-$boxz*round($z/$boxz);
        print "$id $mol 1 $xf $yf $zf\n";
    }; 
  };
  $line=<A>; $line=strip($line); ($x,$y,$z)=split(/ /,$line);
  $id+=1;
  if ($unfolded) {   
    print "$id $mol 2 $x $y $z\n";
  } else {
    $xf=$x-$boxx*round($x/$boxx);
        $xf=$x-$boxx*round($x/$boxx);
        $yf=$y-$boxy*round($y/$boxy);
        $zf=$z-$boxz*round($z/$boxz);
        print "$id $mol $xf $yf $zf\n";
  };
 };
}; 
