#! /usr/bin/perl

# (c) mk@mat.ethz.ch  19 oct 2023 Martin Kroger (ETH Zurich)
# added -ignore-types 29 feb 2024 

sub USAGE { print<<EOF;
use as:\n 
perl $0 <lammps-data-file> [<lammps-dump-trajectory>] [-ignore-types=<type,type,..>]

(c) mk@mat.ethz.ch 19 oct 2023

This software operates on a lammps data file and returns config.Z1, 
a Z1-formatted file that can be used with Z1+. This script identifies the
linear backbones of atomistically detailed chains from the data-file. 
If a dump-trajectory file is not provided, the linear backbone is 
saved in config.Z1. If both a data-file and dump-trajectory-file are provided, 
this script produces a file config.Z1 that contains the trajectory.
The created config.Z1 can be analyzed via: perl ./Z1+ config.Z1

For user's convenience the script furthermore creates a file 
backbone-info.txt. For each chain molecule (chain number, number of atoms) 
it provides a list of consecutively bonded original atom ids that form
the linear chain. 

Optionally, the script takes the argument -ignore-types=<type,type,..>
such as -ignore-types=4 or -ignore-types=2,4. In that case the bead type is 
completely ignored and won't appear in the resulting files. Any corresponding 
bonds to other bead types are ignored as well.

This script, without the -ignore-types option can also be used to
create a Z1-formatted trajectory file from a lammps data and unsorted
lammps dump file.
EOF
exit;
}; 

sub strip { chomp $_[0]; $_[0]=~s/^\s+//g; $_[0]=~s/\s+$//; $_[0]=~s/\s+/ /g; $_[0]; };

if ($#ARGV eq -1) { USAGE; }; 
$datafile=$ARGV[0]; 
if (-s "$datafile") { } else { print "missing file $datafile\n"; exit; }; 
if ($#ARGV > 0) { 
    $dumpfile=$ARGV[1]; 
    if (-s "$dumpfile") { } else { print "missing file $dumpfile\n"; exit; };
}; 
@ignoretype=(); 
foreach $arg (@ARGV) { ($field,$value)=split(/=/,$arg);
    if ($field eq "-ignore-types") { @ignoretype=split(/,/,$value); }; 
};
foreach $j (0 .. $#ignoretype) { print "ignoring type $ignoretype[$j] atoms\n"; }; 

# read data file, give away new ids = row number
open(D,"<$datafile"); $mols=0; $linearatoms=0;
while (!eof(D)) {
    $line=<D>; $line=strip($line); 
    if ($line=~/atoms/) { ($original_atoms,$rest)=split(/ /,$line); };
    if ($line=~/bonds/) { ($original_bonds,$rest)=split(/ /,$line); };  
    if ($line=~/xlo/)   { ($xlo,$xhi,$rest)=split(/ /,$line); };
    if ($line=~/ylo/)   { ($ylo,$yhi,$rest)=split(/ /,$line); };
    if ($line=~/zlo/)   { ($zlo,$zhi,$rest)=split(/ /,$line); };
    if ($line=~/^Atoms/) { $line=<D>;
        $j = 0; 
        foreach $jj (1 .. $original_atoms) { $line=<D>; $line=strip($line);
            $j += 1; 
            @tmp=split(/ /,$line); 
            if ($#tmp eq 8) {
                ($id,$mol[$j],$type[$j],$x[$j],$y[$j],$z[$j],$rest)=split(/ /,$line); 
            } elsif ($#tmp eq 9) {
                ($id,$mol[$j],$type[$j],$q[$j],$x[$j],$y[$j],$z[$j],$rest)=split(/ /,$line);
            } elsif ($#tmp eq 5) {
                ($id,$mol[$j],$type[$j],$x[$j],$y[$j],$z[$j],$rest)=split(/ /,$line);
            } elsif ($#tmp eq 6) {
                ($id,$mol[$j],$type[$j],$q[$j],$x[$j],$y[$j],$z[$j],$rest)=split(/ /,$line);
            } else {
                print "data format problem, we expect either id mol type x y z or id mol type q x y z\n"; exit; 
            };
            $ignore = 0; 
            foreach $k (0 .. $#ignoretype) { if ($type[$j] eq $ignoretype[$k]) { $ignore=1; }; };
            if ($ignore eq 0) {
                $row[$id] = $j;         # -> newly assigned ids
                $anchor[$mol[$j]]=$j; 
                if ($mol[$j]>$mols) { $mols=$mol[$j]; }; 
            } else {    
                $row[$id] = 0; 
                $j -= 1; 
            }; 
        };
        $atoms = $j;    # overwriting atoms
    }; 
    if ($line=~/^Bonds/) { $line=<D>;
        $j = 0; 
        foreach $jj (1 .. $original_bonds) { $line=<D>; $line=strip($line);
            $j += 1; 
            ($bid,$btype,$b1[$j],$b2[$j])=split(/ /,$line);
            $ignore = 0;
            foreach $k (0 .. $#ignoretype) { 
                if ($type[$b1[$j]] eq $ignoretype[$k]) { $ignore=1; };
                if ($type[$b2[$j]] eq $ignoretype[$k]) { $ignore=1; };
            };
            if ($ignore eq 0) { 
                $b1[$j]=$row[$b1[$j]];
                $b2[$j]=$row[$b2[$j]];
            } else {
                $j -= 1; 
            };
        };
        $bonds = $j;    # overwriting bonds
    };
};
if (!$original_atoms) { print "This file $datafile does not seem to be a lammps data file\n"; exit; }; 
if ($#ignoretype>-1) {
    print "$atoms/$original_atoms atoms, $bonds/$original_bonds bonds after erasing atom types, largest mol id is $mols\n";
} else { 
    print "$atoms atoms, $bonds bonds, largest mol id is $mols\n";
}; 
print "box $xlo $xhi $ylo $yhi $zlo $zhi\n";

# create connectivity information
foreach $j (1 .. $bonds) { 
    $conn[$b1[$j]].="$b2[$j] ";
    $conn[$b2[$j]].="$b1[$j] "
};
foreach $j (1 .. $atoms) { $conn[$j]=~s/ $//; }; 

# report
$maxfunc=0; $minfunc=100; 
foreach $j (1 .. $atoms) {
    @tmp=split(/ /,$conn[$j]); 
    if ($#tmp+1>$maxfunc) { $maxfunc=$#tmp+1; }; 
    if ($#tmp+1<$minfunc) { $minfunc=$#tmp+1; }; 
};
print "functionalities range between $minfunc and $maxfunc\n"; 
if (($minfunc eq 1)&&($maxfunc eq 2)) { print "You may try using Z1+ with the -h option (strip H-bonds), but I continue creating a config.Z1 file.\n"; }; 
if (($minfunc eq 2)&&($maxfunc eq 2)) { print "all molecules are linear molecules, there is nothing to do.\n"; exit; };

# loop over molecules, find longest path
@N=();
foreach $m (1 .. $mols) {
    @member=();
    foreach $j (1 .. $atoms) {
        if ($mol[$j] eq $m) { 
            $member[$#member+1]=$j; 
        }; 
    }; 
    $members=$#member+1;
    if ($members>0) { 
        print "mol $m has $members members\n";

        # find one end of the molecule
        foreach $k (@member) { $chemdist[$k]=-1; }; 
        $id = $anchor[$m];
        $dist = 0; 
        $chemdist[$id] = $dist;  
        $found = 1; 
        while ($found eq 1) { 
            $found = 0;
            $dist += 1; 
            foreach $k (@member) { 
                if ($chemdist[$k] eq $dist-1) {
                    @neighbors=split(/ /,$conn[$k]); 
                    foreach $n (@neighbors) { 
                        if ($chemdist[$n] eq -1) { 
                            $chemdist[$n]=$dist; $found=1; $chainbeg=$n; 
                        }; 
                    };     
                };
            };
        }; 
        # report
        # print "molecule $m:\n"; foreach $k (@member) { print "member-id $k type $type[$k] has chemdist $chemdist[$k]\n"; };    

        # find the other end
        foreach $k (@member) { $chemdist[$k]=-1; };
        $id = $chainbeg; 
        $dist = 0;
        $chemdist[$id] = $dist;
        $found = 1;
        while ($found eq 1) {
            $found = 0;
            $dist += 1;
            foreach $k (@member) {
                if ($chemdist[$k] eq $dist-1) {
                    @neighbors=split(/ /,$conn[$k]);
                    foreach $n (@neighbors) {
                        if ($chemdist[$n] eq -1) { $chemdist[$n]=$dist; $found=1; $chainend=$n; };
                    };
                };
            };
        }; 
        $dist-=1; 
        # report
        # print "molecule $m:\n"; foreach $k (@member) { print "member-id $k type $type[$k] has chemdist $chemdist[$k]\n"; };

        print "longest path with $dist == $chemdist[$chainend] bonds from id $chainbeg [type $type[$chainbeg]] to $chainend [type $type[$chainend]]\n"; 
    
        # now follow the reverse path to construct the linear chain
        @list = ();
        $list[0] = $chainend; 
        # print "id $chainend at distance $chemdist[$chainend]\n";
        while ($dist > 0) { 
            $dist -= 1; 
            @neighbors = split(/ /,$conn[$list[$#list]]);
            foreach $n (@neighbors) {
                if ($chemdist[$n] eq $dist) {
                    $list[$#list+1] = $n; 
                };
            }; 
        }; 
        $n=$#list+1; 
        $N[$#N+1]=$n;
    
        # report
        # print "linear molecule $m: @list\n";

        $TRANS[$#TRANS+1]="$m $n\n";
        foreach $j (0 .. $#list) { $k=$j+1; 
            $XYZ                   .= "$x[$list[$j]] $y[$list[$j]] $z[$list[$j]]\n";    
            $TRANS[$#TRANS+1]       = "$list[$j]\n";
            $LINEAR_mol[$list[$j]]  = $m;
            $LINEAR_id[$list[$j]]   = $linearatoms+$j+1;
        }; 
        $linearatoms += $n;
    }; 

};

print "$linearatoms linear atoms\n"; 

# save Z1-formatted data file
$molecules = $#N+1; 
$boxx = $xhi-$xlo;
$boxy = $yhi-$ylo;
$boxz = $zhi-$zlo;
$Z1="$molecules\n$boxx $boxy $boxz\n";
foreach $j (0 .. $#N) { $Z1.="$N[$j] "; }; $Z1.="\n$XYZ";
open(Z1,">config.Z1"); print Z1 $Z1; close(Z1); 
print "created config.Z1\n";

# dump translation table for visualization purposes
open(D,">backbone-info.txt"); print D @TRANS; close(D);
print "created backbone-info.txt\n";

if (-s "$dumpfile") {
    open(D,"<$dumpfile"); 
    open(Z1,">config.Z1"); 
    while (!eof(D)) {
        $line=<D>; $line=strip($line); 
        if      ($line =~ /NUMBER OF ATOMS/) { $line=<D>; $line=strip($line);
            $dumpatoms=$line+0; 
            if ($dumpatoms eq $original_atoms) { } else { print "conflicting data and dump files [$original_atoms versus $dumpatoms atoms]\n"; exit; }; 
        } elsif ($line =~ /ITEM: ATOMS id/) {
            @XYZ=();
            @tmp=split(/ /,$line);   
            $col=-1;
            foreach $k (0 .. $#tmp) { 
                if (($tmp[$k] eq "x")||($tmp[$k] eq "xu")) { $col=$k-2; };
            };
            if ($col eq -1) { print "format error in dumpfile $dumpfile (must contain x or xu)\n"; };
            foreach $j (1 .. $original_atoms) {
                $line=<D>; $line=strip($line);
                @tmp   = split(/ /,$line);
                $id    = $row[$tmp[0]]; 
                if ($id) { 
                    $xyz   = "$tmp[$col] $tmp[$col+1] $tmp[$col+2]";
                    if ($LINEAR_mol[$id]) { $XYZ[$LINEAR_id[$id]] = "$xyz\n"; }; 
                }; 
            };
            print Z1 "$molecules\n$boxx $boxy $boxz\n";
            foreach $j (0 .. $#N)   { print Z1 "$N[$j] "; }; print Z1 "\n";
            print Z1 @XYZ; 
            # print Z1 "0\n";     # separator
        } elsif ($line =~ /TIMESTEP/) {
            $line=<D>+0; print "[$dumpfile] processing time step $line\n"; 
        };
    };  
    close(D);
    close(Z1);
    print "created config.Z1 (trajectory file)\n";
};
