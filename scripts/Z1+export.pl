#! /usr/bin/perl

$config = "Z1+initconfig.dat";
$SP     = "Z1+SP.dat";

sub strip      { chomp $_[0]; $_[0]=~s/^\s+//g; $_[0]=~s/\s+$//; $_[0]=~s/\s+/ /g; $_[0]; }; 
sub greenprint { print "\033[1;32;40m$_[0]\033[1;37;40m\n"; };
sub redprint   { print "\033[1;31;40m$_[0]\033[1;37;40m\n"; };


sub USAGE { 
greenprint("\nperl $0 -all");
print "OR\n";
greenprint("perl $0 -data -all"); 
print "OR\n";
greenprint("perl $0 -from=<frame-no> -to=<frame-no>");
print "OR\n";
greenprint("perl $0 -data -from=<frame-no> -to=<frame-no>");
print<<EOF;
\nThis script can be called immediately after Z1+ finished. 
It creates a dump trajectory file containing both original system (bead types 1,2,3) 
and primitive path (types 4,5,6), where original and primitive chain belong to the same molecule id. 
Using -data, the script creates one data file for each of the selected frames. 
Select frames via -from=.. -to=... 
EOF
exit;
};

if ($#ARGV eq -1) { USAGE; }; 

$dump=1; # default
$from=1; # default
$to=1;   # default
foreach $arg (@ARGV) { ($field,$value)=split(/=/,$arg);
    if ($field eq "-data")      { $data=1; $dump=0;  
    } elsif ($field eq "-dump") { $dump=1; $data=0; 
    } elsif ($field eq "-from") { $from=$value; 
    } elsif ($field eq "-to")   { $to=$value; 
    } elsif ($field eq "-all")  { $from=1; $to=1e7; 
    };     
};

sub create_merged { 
$frame = 0;

open(C,"<$config"); 
open(S,"<$SP");
while (!eof(C)) { 
$frame += 1; 
$chains=<C>+0; # print "frame $frame with $chains";
$line=<S>+0;   # print " == $line chains ";
$line=<C>; $line=strip($line); ($boxx,$boxy,$boxz)=split(/ /,$line); # print "in box $boxx $boxy $boxz ";
$xlo=-$boxx/2;
$xhi= $boxx/2;
$ylo=-$boxy/2;
$yhi= $boxy/2;
$zlo=-$boxz/2;
$zhi= $boxz/2;
$line=<S>;
if ($dump eq 1) { $ZEROS = ""; } else { $ZEROS=" 0 0 0"; }; 
foreach $mol (1 .. $chains) {
    $N=<C>+0;
    # ---
    $line=<C>; $line=strip($line); ($x,$y,$z,$rest)=split(/ /,$line);
    $type = 1;
    $id += 1; 
    $ATOM .= "$id $mol $type $x $y $z$ZEROS\n"; 
    # ---
    foreach $i (2 .. $N-1) {
        $line=<C>; $line=strip($line); ($x,$y,$z,$rest)=split(/ /,$line); 
        $type = 2;  
        $btype = 1; 
        $id0 = $id;
        $id += 1;
        $bid += 1;
        $ATOM .= "$id $mol $type $x $y $z$ZEROS\n";
        $BOND .= "$bid $btype $id0 $id\n";
    };
    # ---
    $line=<C>; $line=strip($line); ($x,$y,$z,$rest)=split(/ /,$line);
    $type = 3;
    $btype = 1;
    $id0 = $id;
    $id += 1;
    $bid += 1;
    $ATOM .= "$id $mol $type $x $y $z$ZEROS\n";
    $BOND .= "$bid $btype $id0 $id\n";
    # ---
    $n=<S>+0;
    # ---
    $line=<S>; $line=strip($line); ($X,$Y,$Z,$rest)=split(/ /,$line);
    $type = 4; 
    $ATOM .= "$id $mol $type $X $Y $Z$ZEROS\n";
    # ---
    foreach $i (2 .. $n-1) {
        $line=<S>; $line=strip($line); ($X,$Y,$Z,$rest)=split(/ /,$line);
        $type = 5;
        $btype = 2; 
        $id0 = $id; 
        $id += 1;
        $bid += 1;
        $ATOM .= "$id $mol $type $X $Y $Z$ZEROS\n";
        $BOND .= "$bid $btype $id0 $id\n";
    };
    # ---
    $line=<S>; $line=strip($line); ($X,$Y,$Z,$rest)=split(/ /,$line);
    $type = 6;
    $btype = 2;
    $id0 = $id;
    $id += 1;
    $bid += 1;
    $ATOM .= "$id $mol $type $X $Y $Z$ZEROS\n";
    $BOND .= "$bid $btype $id0 $id\n";
};

# --- trailer
$line=<C>; $line=strip($line);
$line=<C>; $line=strip($line);

chomp $ATOM;
chomp $BOND;

# -----------------
if (($frame>=$from)&&($frame<=$to)) { 
# -----------------
if ($data) { 
print "created Z1+merged-$frame.data with $id merged atoms and $bid merged bonds\n"; 
open(O,">Z1+merged-$frame.data"); 
print O<<EOF;
LAMMPS data file via Z1+export.pl (mk\@mat.ethz.ch)

$id atoms
6 atom types
$bid bonds
2 bond types

$xlo $xhi xlo xhi
$ylo $yhi ylo yhi
$zlo $zhi zlo zhi

Atoms # angle

$ATOM

Bonds 

$BOND
EOF
close(O);
}; 
# -----------------
if ($dump) {
if ($frame eq $from) { `rm -f Z1+merged.dump`; print "creating Z1+merged.dump\n"; }; 
print "added frame $frame to Z1+merged.dump with $id merged atoms and $bid merged bonds\n";
open(O,">>Z1+merged.dump"); 
print O<<EOF;
ITEM: TIMESTEP
$frame
ITEM: NUMBER OF ATOMS
$id
ITEM: BOX BOUNDS pp pp pp
$xlo $xhi
$ylo $yhi
$zlo $zhi
ITEM: ATOMS id mol type x y z
$ATOM
EOF
close(O);
};
# -----------------
}; 
# -----------------

}; # while frame

close(C);
close(S);
}; 

create_merged;
