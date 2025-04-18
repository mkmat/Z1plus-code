#! /usr/bin/perl

sub USAGE { print<<EOF; 

NAME 
    Z1+SP-to-data.pl

SYNOPSIS 
    perl $0 [-i] [-s] [-si] 

DESCRIPTION 
    This script converts the initial configuration and the shortest path configuration
    obtained by Z1+ to lammps-data formatted files so that they can be visualized using vmd,
    ovito, or other software that recognizes lammps-data format.

    Atom types:
        1   first atom of a chain in the initial configuration
        2   non-terminal atom of a chain in the initial configuration
        3   last atom of a chain in the initial configuration

        4   first node of a primitive path in the shortest path configuration
        5   entangled node of a primitive path in the shortest path configuration
        6   last node of a primitive path in the shortest path configuration

    Bond types:
        1   bond in the initial configuration
        2   bond in the primitive path


    (c) 20 april 2024 Martin Kroger mk\@mat.ethz.ch addon for the Z1+ package,
    hosted at https://github.com/mkmat/Z1plus-code

OPTIONS

    -i   
        converts existing Z1+initconfig.dat to Z1+initconfig.data
    -s  
        converts existing Z1+SP.dat to Z1+SP.data
    -si 
        merges initial and shortest path configuration in a single file Z1+merge.data
        If the original system carries M molecules, the molecule number of the 
        shortest path corresponding to molecule X is X+M. The merged file thus contains
        2*M molecules, each (X,X+M) pair of molecules shares the positions of their terminal 
        atoms (nodes). If the -si option is chosen, and if Z1+ had been called with the +
        or -+SP option, the Z1+SP.dat file contains information about atom IDs corresponding
        to interior kinks of the shortest path. The atom IDs are mentioned in the 
        Z1+merge.data file after entangled node IDs in the Atoms section (before Bonds section).
        
EOF
exit;
};

$IN = "Z1+initconfig.dat"; 
$SP = "Z1+SP.dat"; 

$atomtypes = 6;
$bondtypes = 2; 

$option_si = 0;
if ($#ARGV eq -1) { USAGE; }; 
foreach $arg (@ARGV) { 
    if      ($arg eq "-i")  { $option_i=1; 
    } elsif ($arg eq "-s")  { $option_s=1; 
    } elsif ($arg eq "-si") { $option_i=1; $option_s=1; $option_si=1; 
    } else { USAGE; 
    }; 
};

sub strip { chomp $_[0]; $_[0]=~s/^\s+//g; $_[0]=~s/\s+$//; $_[0]=~s/\s+/ /g; $_[0]; }; 
sub round { if ($_[0] eq 0) { } else { $_[0]/=abs($_[0]); $_[0]*=int(abs($_[0])+0.5); }; $_[0]; };

# ------------------------
if ($option_i) { 
# ------------------------
if (-s "$IN") { } else { print "missing $IN\n"; exit; };
open(IN,"<$IN"); {
    $id   = 0;
    $bonds = 0;
    @mol  = ();
    @type = ();
    @x = ();
    @y = ();
    @z = ();
    @b1 = ();
    @b2 = ();
    $chains=<IN>+0;
    print "$chains chains .. ";
    $line=<IN>; $line=strip($line); ($boxx,$boxy,$boxz)=split(/ /,$line);
    foreach $c (1 .. $chains) {
        $N[$c] = <IN>+0;
        # print "chain $c with $N[$c] atoms [id=$id]\n";
        foreach $j (1 .. $N[$c]) { 
            $id+=1; 
            $line=<IN>; $line=strip($line); ($xu[$id],$yu[$id],$zu[$id])=split(/ /,$line);
            $mol[$id]  = $c; 
            if ($j eq 1) { 
                $type[$id] = 1;  
                $firstID_of_chain[$c]=$id; 
            } elsif ($j eq $N[$c]) {
                $type[$id] = 3;
                $bonds+=1; 
                $b1[$bonds] = $id-1;
                $b2[$bonds] = $id;
            } else {
                $type[$id] = 2; 
                $bonds+=1; 
                $b1[$bonds] = $id-1;
                $b2[$bonds] = $id;
            };
        }; 
    };
};
close(IN);
$atoms = $id; 
print "$atoms atoms .. ";
print "$bonds bonds\n";
$xlo=-$boxx/2; $xhi=-$xlo; 
$ylo=-$boxy/2; $yhi=-$ylo;
$zlo=-$boxz/2; $zhi=-$zlo;
foreach $id (1 .. $atoms) {
    $ix[$id] = round($xu[$id]/$boxx); 
    $x[$id]  = $xu[$id] - $boxx*$ix[$id];
    $iy[$id] = round($yu[$id]/$boxy);
    $y[$id]  = $yu[$id] - $boxy*$iy[$id];
    $iz[$id] = round($zu[$id]/$boxz);
    $z[$id]  = $zu[$id] - $boxz*$iz[$id];
};

if ($option_si eq 0) { 
$INdata = $IN; $INdata.="a";
open(DATA,">$INdata"); 
print DATA<<EOF;
lammps data file created by $0 (software at https://github.com/mkmat/Z1plus-code)

$atoms atoms
$bonds bonds
$atomtypes atom types
$bondtypes bond types

$xlo $xhi xlo xhi
$ylo $yhi ylo yhi
$zlo $zhi zlo zhi

Atoms

EOF
foreach $id (1 .. $atoms) {
    print DATA "$id $mol[$id] $type[$id] $x[$id] $y[$id] $z[$id] $ix[$id] $iy[$id] $iz[$id]\n";
};
print DATA<<EOF;

Bonds

EOF
foreach $bid (1 .. $bonds) {
    print DATA "$bid 1 $b1[$bid] $b2[$bid]\n";
};
close(DATA);
print "created: $INdata\n";
}; 
# ------------------------
};
# ------------------------

# ------------------------
if ($option_s) { 
# ------------------------
if (-s "$SP") { } else { print "missing $SP\n"; exit; };
open(IN,"<$SP"); {
    $id   = 0;
    $SPbonds = 0;
    @SPmol  = ();
    @SPtype = ();
    @SPx = ();
    @SPy = ();
    @SPz = ();
    @SPb1 = ();
    @SPb2 = ();
    $SPchains=<IN>+0;
    print "$SPchains SP chains .. ";
    $line=<IN>; $line=strip($line); ($boxx,$boxy,$boxz)=split(/ /,$line);
    foreach $c (1 .. $SPchains) {
        $SPN[$c] = <IN>+0;
        foreach $j (1 .. $SPN[$c]) {
            $id+=1;
            $line=<IN>; $line=strip($line); ($SPxu[$id],$SPyu[$id],$SPzu[$id],$SPpos,$entangled,$entmol,$entbead)=split(/ /,$line);
            $SPmol[$id]  = $c;
            if ($j eq 1) {
                $SPtype[$id] = 4;
            } elsif ($j eq $N[$c]) {
                $SPtype[$id] = 6;
                $SPbonds+=1;
                $SPb1[$SPbonds] = $id-1;
                $SPb2[$SPbonds] = $id;
            } else {
                $SPtype[$id] = 5;
                $SPbonds+=1;
                $SPb1[$SPbonds] = $id-1;
                $SPb2[$SPbonds] = $id;
                if ($entmol) { $entangled_with_ID[$id] = $firstID_of_chain[$entmol] + $entbead-1; 
                               $entangled_with_ID[$id] = " # $entangled_with_ID[$id]";
                };
            };
        };
    };
};
close(IN);
$SPatoms = $id;
print "$SPatoms SP atoms .. ";
print "$SPbonds SP bonds\n";
$xlo=-$boxx/2; $xhi=-$xlo;
$ylo=-$boxy/2; $yhi=-$ylo;
$zlo=-$boxz/2; $zhi=-$zlo;
foreach $id (1 .. $SPatoms) {
    $SPix[$id] = round($SPxu[$id]/$boxx);
    $SPx[$id]  = $SPxu[$id] - $boxx*$SPix[$id];
    $SPiy[$id] = round($SPyu[$id]/$boxy);
    $SPy[$id]  = $SPyu[$id] - $boxy*$SPiy[$id];
    $SPiz[$id] = round($SPzu[$id]/$boxz);
    $SPz[$id]  = $SPzu[$id] - $boxz*$SPiz[$id];
};
if ($option_si eq 0) { 
$SPdata=$SP; $SPdata.="a";
open(DATA,">$SPdata");
print DATA<<EOF;
lammps data file created by $0 (software at https://github.com/mkmat/Z1plus-code)

$SPatoms atoms
$SPbonds bonds
$atomtypes atom types
$bondtypes bond types

$xlo $xhi xlo xhi
$ylo $yhi ylo yhi
$zlo $zhi zlo zhi

Atoms

EOF
foreach $id (1 .. $SPatoms) {
    print DATA "$id $SPmol[$id] $SPtype[$id] $SPx[$id] $SPy[$id] $SPz[$id] $SPix[$id] $SPiy[$id] $SPiz[$id]\n";
};
print DATA<<EOF;

Bonds

EOF
foreach $bid (1 .. $SPbonds) {
    print DATA "$bid 1 $SPb1[$bid] $SPb2[$bid]\n";
};
close(DATA);
print "created: $SPdata\n";
}; 
# ------------------------
};
# ------------------------


# ------------------------
if ($option_si) { 
# ------------------------
$atoms_both = $atoms + $SPatoms; 
$bonds_both = $bonds + $SPbonds; 

open(DATA,">Z1+merge.data");
print DATA<<EOF;
lammps data file created by $0 (software at https://github.com/mkmat/Z1plus-code)

$atoms_both atoms
$bonds_both bonds
$atomtypes atom types
$bondtypes bond types

$xlo $xhi xlo xhi
$ylo $yhi ylo yhi
$zlo $zhi zlo zhi

Atoms

EOF
foreach $id (1 .. $atoms) {
    print DATA "$id $mol[$id] $type[$id] $x[$id] $y[$id] $z[$id] $ix[$id] $iy[$id] $iz[$id]\n";
};
foreach $SPid (1 .. $SPatoms) {
    $id = $atoms+$SPid; 
    $MOL = $SPmol[$SPid]+$chains; 
    print DATA "$id $MOL $SPtype[$SPid] $SPx[$SPid] $SPy[$SPid] $SPz[$SPid] $SPix[$SPid] $SPiy[$SPid] $SPiz[$SPid] $entangled_with_ID[$SPid]\n";
};
print DATA<<EOF;

Bonds

EOF
foreach $bid (1 .. $bonds) {
    print DATA "$bid 1 $b1[$bid] $b2[$bid]\n";
};
foreach $SPbid (1 .. $SPbonds) {
    $bid = $bonds+$SPbid; 
    $B1 = $SPb1[$SPbid]+$atoms;
    $B2 = $SPb2[$SPbid]+$atoms;
    print DATA "$bid 2 $B1 $B2\n";
};
close(DATA);
print "created: Z1+merge.data\n";
# ------------------------
}; 
# ------------------------
