#! /susr/bin/perl

sub USAGE { print<<EOF;

perl $0 -data=<data-file> [-dump=<dump-file>] [-angle]

[helper file for Z1+ software at https://github.com/mkmat/Z1plus-code]

Converts a lammps data file that does not contain proper information
about molecules into a proper data file that contains this information.
If a dump-file is provided in addition, also the dump-file is corrected
using the bond information from the data file. The new data and dump
files have the names of the original files, followed by "-corrected".

By default, this script assumes that the lammps atom style is "full".
In that case the Atoms section contains: id mol type q x y z .. 
If the atom style is "angle" or "bond", add the -angle option.
If you use a different atom style, edit this script at the position marked by "STYLE".

Script written by (c) 2024 Martin Kroger mk\@mat.ethz.ch
EOF
exit;
};

if ($#ARGV eq -1) { USAGE; }; 

$full=1; # the default
foreach $arg (@ARGV) { ($field,$value)=split(/=/,$arg);
    if ($field eq "-data") { $data=$value; 
    } elsif ($field eq "-dump") { $dump=$value; 
    } elsif ($field eq "-angle") { $full=0; 
    } else { USAGE; 
    };
}; 

if (-s "$data") { } else { print "____________\n\nmissing data file\n____________\n"; `sleep 1`; USAGE; }; 
if ($dump) { if (-s "$dump") { } else { print "missing dump file [$dump]\n"; `sleep 1`; USAGE; }; }; 

sub strip { chomp $_[0]; $_[0]=~s/^\s+//g; $_[0]=~s/\s+$//; $_[0]=~s/\s+/ /g; $_[0]; };

open(DATA,"<$data"); $max_Id=0;
while (!eof(DATA)) { 
    $line=<DATA>; $line=strip($line); @values=split(/ /,$line);
    if      ($line =~ /atoms/) { $atoms=$values[0]; 
    } elsif ($line =~ /bonds/) { $bonds=$values[0]; print "data file has $atoms atoms, $bonds bonds, "; 
    } elsif ($line =~ /xlo/)   { $xlo=$values[0]; $xhi=$values[1]; 
    } elsif ($line =~ /ylo/)   { $ylo=$values[0]; $yhi=$values[1];
    } elsif ($line =~ /zlo/)   { $zlo=$values[0]; $zhi=$values[1];
    } elsif ($line =~ /^Atoms/) { $line=<DATA>; 
        foreach $i (1 .. $atoms) {
            $line=<DATA>; $line=strip($line); @values=split(/ /,$line);
            $id         = $values[0]; if ($id>$max_ID) { $max_ID=$id; }; 
            $nomol      = $values[1];   # ignore this entry
            $type[$id]  = $values[2];
            if ($full) { 
                $q[$id] = $values[3];
                $x[$id] = $values[4];   # atom style full [STYLE]
                $y[$id] = $values[5];
                $z[$id] = $values[6];
            } else {
                $x[$id] = $values[3];
                $y[$id] = $values[4];   # atom style bond or angle [STYLE]
                $z[$id] = $values[5];
            }; 
        }; 
    } elsif ($line =~ /^Bonds/) { $line=<DATA>;
        foreach $i (1 .. $bonds) {
            $line=<DATA>; $line=strip($line); @values=split(/ /,$line);
            $bid          = $values[0];
            $btype[$bid]  = $values[1]; 
            $b1[$bid]     = $values[2];
            $b2[$bid]     = $values[3];
        };
    }; 
};
close(DATA);
print "max ID $max_ID, ";

$chains = 0; 
foreach $id (1 .. $atoms) { $mol[$id]=0; }; 
foreach $bid (1 .. $bonds) {
    if      ($mol[$b1[$bid]]) {
        $mol[$b2[$bid]] = $mol[$b1[$bid]];  
    } elsif ($mol[$b2[$bid]]) { 
        $mol[$b1[$bid]] = $mol[$b2[$bid]];
    } else { 
        $chains+=1; 
        $mol[$b1[$bid]] = $chains;
        $mol[$b2[$bid]] = $chains; 
    }    
};
print "$chains chains.\n";

open(DATAC,">$data-corrected");
open(DATA,"<$data"); $max_Id=0;
while (!eof(DATA)) {
    $line=<DATA>; $line=strip($line); @values=split(/ /,$line);
    if      ($line =~ /^Atoms/) { print DATAC "@values\n\n"; $line=<DATA>;
        foreach $i (1 .. $atoms) {
            $line=<DATA>; $line=strip($line); @values=split(/ /,$line);
            $id     = $values[0]; 
            if ($full) { 
                print DATAC "$id $mol[$id] $type[$id] $q[$id] $x[$id] $y[$id] $z[$id]\n";
            } else {
                print DATAC "$id $mol[$id] $type[$id] $x[$id] $y[$id] $z[$id]\n";
            }; 
        };
    } else { 
        print DATAC "@values\n";              
    };
};
close(DATA);
close(DATAC); 
print "created $data-corrected\n";

if ($dump) { 
    open(DUMPC,">$dump-corrected"); 
    open(DUMP,"<$dump");
    $frames=0;
    while (!eof(DUMP)) {
        $line=<DUMP>; $line=strip($line); @values=split(/ /,$line);
        if      ($line=~/ITEM: NUMBER OF ATOMS/) { $dumpatoms=<DUMP>+0; print DUMPC "@values\n$dumpatoms\n"; 
            $frames+=1; 
            if ($atoms eq $dumpatoms) { } else { print "incompatible data and dump files. $dumpatoms atoms in dump file\n"; exit; }; 
        } elsif ($line=~/ITEM: ATOMS/) { 
            $col_id = -1; 
            $col_x = -1; 
            foreach $j (2 .. $#values) {
                if      ($values[$j] eq "id")   { $col_id=$j-2;        if ($frames eq 1) { print "id found in column $col_id\n"; };
                } elsif ($values[$j] eq "mol")  { $col_mol=$j-2;  
                } elsif ($values[$j] eq "type") { $col_type=$j-2; 
                } elsif ($values[$j] eq "x")    { $col_x=$j-2; $X="x"; if ($frames eq 1) { print "x found in column $col_x\n"; };
                } elsif ($values[$j] eq "y")    { $col_y=$j-2; $Y="y";
                } elsif ($values[$j] eq "z")    { $col_z=$j-2; $Z="z";
                } elsif ($values[$j] eq "xu")   { $col_x=$j-2; $X="xu"; if ($frames eq 1) { print "xu found in column $col_x\n"; };
                } elsif ($values[$j] eq "yu")   { $col_y=$j-2; $Y="yu";
                } elsif ($values[$j] eq "zu")   { $col_z=$j-2; $Z="zu";
                }; 
            };
            if ($col_id eq -1)   { print "missing id information in dump file\n"; exit; }; 
            if ($col_x eq -1)    { print "missing x or xu information in dump file\n"; exit; };
            print DUMPC "ITEM: ATOMS id mol type $X $Y $Z\n";
            foreach $j (1 .. $atoms) { 
                $line=<DUMP>; $line=strip($line); @values=split(/ /,$line);                             
                $id = $values[$col_id];
                print DUMPC "$id $mol[$id] $type[$id] $values[$col_x] $values[$col_y] $values[$col_z]\n";
            };
        } else {
            print DUMPC "@values\n";
        };
    };
    close(DUMP);
    close(DUMPC);
    print "created $dump-corrected with $frames frames\n";
};
