#! /usr/bin/perl

# (c) 01 June 2022 mk@mat.ethz.ch (added -ignore_types)
# this software should not be re-distributed

# 22 feb 2023 added CREATE_ID_CONVERSION_TABLE
# 13 nov 2024 added xy in USE_COORDINATES_FROM_DATA_FILE

sub USAGE { if ($ERROR) { $ADD="**** ERROR ****"; }; print<<EOF;

_______________________ make_Z1+_on_lammps ___________________________

$ADD $ERROR 

usage: perl ./Z1+import-lammps.pl
    -dump=<dumpfile> 
    -data=<datafile>
    [-from=<snapshot no>]
    [-to=<snapshot no>]
    [-each=<number>]
    [-ignore_H]
    [-ignore_dumbbells]
    [-ignore_types=<type>[,<type>]]
    [-Z1]
    [-self-entanglements]
    [-branched]
    [-xml=<xml-file>]
    [-out-dump=<dumpfile>]   format id mol x y z
    [-verbose]

Example:
perl ./Z1+import-lammps.pl -dump=example.dump -data=example.data -ignore_H -from=1 -to=2

____________________________________________________________________

This tool usually requires 2 files, both conveniently created from within lammps: 

1) <dumpfile> (trajectory of snapshots) is generated via the LAMMPS command
dump ID group-ID custom N dumpfile id ... x y z ...     OR 
dump ID group-ID custom N dumpfile id ... xu yu zu ...  OR
dump ID group-ID custom N dumpfile id ... xs ys zs ...
Example: 
dump mydump all custom 100 file.dump id mol type xu yu zu

2) <datafile> (containing connectivity information) is generated via the LAMMPS command
write_data <datafile> 
Example:
write_data file.data

_____________________________________________________________________

If only a <dumpfile> (snapshot or trajectory) is available, it must have
been saved using: dump_modify sort id
and every chain should have a different mol number, adjacent atoms within 
chains must have adjacent ids. Then call this tool without the -data option.
_____________________________________________________________________

If only a <datafile> is available, the analysis is performed on the coordinates 
contained in the <datafile>. Then call this tool without the -dump option.

_____________________________________________________________________
EOF
exit;
}; 

sub EX { if (-s "$_[0]") { } else { print "file $_ does not exist.\n"; USAGE; }; };
sub INT { return $_[0]+0; }; 
sub STRIP_WITH_COMMENTS { $_[0]=~s/\t/ /g; $_[0]=~s/\r/ /g; $_[0]=~s/^\s+//g; $_[0]=~s/\s+$//; $_[0]=~s/\s+/ /g; return $_[0]; }; 
sub STRIP { my @var=split(/\#/,$_[0]); return STRIP_WITH_COMMENTS($var[0]); };
sub min { if ($_[0]>$_[1]) { return($_[1]); } else { return($_[0]); }; };
sub max { if ($_[0]>$_[1]) { return($_[0]); } else { return($_[1]); }; };

sub PARSE {
    $REQ=0; $ignoreH=0; @ignoreType=();
    $from=1; 
    $each=1; 
    $to=10**10;
    $self="";
    $verbose=0;  
    foreach $i (0 .. $#ARGV) {
        if      ($ARGV[$i]=~/^-dump=/)    { $dumpfile=$ARGV[$i]; $dumpfile=~s/-dump=//; EX($dumpfile); $REQ+=1; $REQDUMP=1;  
        } elsif ($ARGV[$i]=~/^-data=/)    { $datafile=$ARGV[$i]; $datafile=~s/-data=//; EX($datafile); $REQ+=1; $REQDATA=1;  
        } elsif ($ARGV[$i]=~/^-from=/)    { ($arg,$from)=split(/=/,$ARGV[$i]); $from=INT($from); 
        } elsif ($ARGV[$i]=~/^-to=/)      { ($arg,$to)=split(/=/,$ARGV[$i]); $to=INT($to); 
        } elsif ($ARGV[$i]=~/^-each=/)    { ($arg,$each)=split(/=/,$ARGV[$i]); $each=INT($each); 
        } elsif ($ARGV[$i]=~/^-ignore_H/) { $ignoreH=1; 
        } elsif ($ARGV[$i]=~/^-ignore_dumbbells/) { $ignore_dumbbells=1; 
        } elsif ($ARGV[$i]=~/^-ignore_types/) { ($arg,$ignoretypes)=split(/=/,$ARGV[$i]); @ignoreType=split(/,/,$ignoretypes); 
        } elsif ($ARGV[$i]=~/^-Z1/)       { $Z1=1; 
        } elsif ($ARGV[$i]=~/^-verbose/)  { $verbose=1;  
        } elsif ($ARGV[$i]=~/^-self-entanglements/) { $self="self";  
        } elsif ($ARGV[$i]=~/^-xml/)      { $xml=1; $xmlfile=$ARGV[$i]; $xmlfile=~s/-xml=//; EX($xmlfile); $REQ=2; 
        } elsif ($ARGV[$i]=~/^-branched/) { $accept_branched=1; 
        } elsif ($ARGV[$i]=~/^-out-dump/) { ($arg,$outdumpfile)=split(/=/,$ARGV[$i]); 
        } elsif ($ARGV[$i] eq "")         { 
        } else { $ERROR="option >$ARGV[$i]< does not exist.\n"; USAGE; };
    };
    if ($REQ eq 2) { } else { 
        if ($REQDUMP eq 1) { } else { 
            if ($REQDATA eq 1) { } else { 
                $ERROR="missing dump AND/OR data file\n"; USAGE; 
            }; 
        }; 
    };
};

sub CREATE_ID_CONVERSION_TABLE {
    foreach $ic (1 .. $chains) { $beginID[$ic] = $ATOMS; }; 
    foreach $id (0 .. $maxid) {
        if ($newid[$id] < $beginID[$newmol[$id]]) { $beginID[$newmol[$id]] = $newid[$id]; }; 
    }; 
    open(CT,">table-lammps-ids-to-Z1-id-mol-bead.txt");
    foreach $id (0 .. $maxid) {
        if ($newid[$id]) {
            $beadno = $newid[$id] - $beginID[$newmol[$id]] + 1; 
            print CT "$id $mol[$id] $newid[$id] $newmol[$id] $beadno\n";    # lammps atom id, lammps mol id, Z1 id, Z1 mol id, Z1 bead-in-chain
        };
    };
    close(CT);
};

sub OPEN_DATA_FILE  { open(DATA,"<$datafile"); };
sub CLOSE_DATA_FILE { close(DATA); CREATE_ID_CONVERSION_TABLE; }; 

sub DEBUG {
    print "--------------- DEBUG ----------------\n";
    print "maxid=$maxid atoms=$atoms chains=$chains\n";
    foreach $trymol (1 .. $chains) {
        foreach $id (0 .. $maxid) {
            if ($newid[$id]) {
                if ($newmol[$id] eq $trymol) {
                    print "id=$id newid=$newid[$id] mass=$mass[$id] type=$type[$id]=$TYPE[$newid[$id]] ";
                    print "newmol=$newmol[$id]=$MOL[$newid[$id]] N=$N[$newmol[$id]] conn=$conn[$id]\n";
                };
            };
        };
    };
};

sub IGNORE_DUMBBELLS {
    if (!$ignore_dumbbells) { return; }; 
    if ($verbose) { print "ignoring dumbbells ..\n"; };
    @newN=();
    $newATOMS=0;
    $newchains=0;
    @newx=(); 
    @newy=(); 
    @newz=(); 
    $j=0;
    foreach $i (1 .. $chains) {
        if ($N[$i] eq 2) {
            $j+=$N[$i]; 
            if ($verbose) { print "ignore dumbbell chain $i\n"; };
        } else {
            $newchains+=1; 
            $newN[$newchains]=$N[$i];     
            foreach $k (1 .. $N[$i]) { 
                $j+=1; 
                $newATOMS+=1; 
                $newx[$newATOMS]=$x[$j]; 
                $newy[$newATOMS]=$y[$j];
                $newz[$newATOMS]=$z[$j];
            }; 
        };
    };
    $ATOMS=$newATOMS;
    $chains=$newchains;
    @N=@newN;
    @x=@newx; @newx=();
    @y=@newy; @newy=();
    @z=@newz; @newz=();
    $Ns=""; foreach $i (1 .. $chains) { $Ns.="$N[$i] "; }; $Ns=~s/ $//;
};

sub OPTIONAL_WRITE_TRUECHAINS_TO_DUMP {		# added 15 oct 2019
    if ($outdumpfile) { } else { return; }; 
    open(OUTDUMP,">$outdumpfile"); 
    print OUTDUMP "ITEM: TIMESTEP\n0\nITEM: ATOMS\n$ATOMS\n";
    if ($xy) { 
        print OUTDUMP "ITEM: BOX BOUNDS xy\n";
        print OUTDUMP "$xlo $xhi $xy\n";
        print OUTDUMP "$ylo $yhi 0\n";
        print OUTDUMP "$zlo $zhi 0\n";
    } else {
        print OUTDUMP "ITEM: BOX BOUNDS\n";
        print OUTDUMP "$xlo $xhi\n";
        print OUTDUMP "$ylo $yhi\n";
        print OUTDUMP "$zlo $zhi\n";
    };
    print OUTDUMP "ITEM: ATOMS id mol x y z\n";
    foreach $i (1 .. $ATOMS) {
        print OUTDUMP "$i $MOL[$i] $x[$i] $y[$i] $z[$i]\n";
    };
    close(OUTDUMP);  
};

sub APPEND_TO_Z1_FORMATTED {
    IGNORE_DUMBBELLS;
    if ($verbose) { print "$ATOMS atoms, $chains chains, box sizes $boxx $boxy $boxz\n"; }; 
    if ($verbose) { print "Ns: $Ns\n"; }; 
    if ($verbose) { print "creating Z1-formatted\n"; };
    open(Z,">>config.Z1");
    print Z<<EOF;
$chains
$boxx $boxy $boxz
$Ns
EOF
    foreach $i (1 .. $ATOMS) { print Z "$x[$i] $y[$i] $z[$i]\n"; };
    if ($xz eq 0) { } else { print "ERROR: BOX xz <> 0"; };
    if ($yz eq 0) { } else { print "ERROR: BOX yz <> 0"; };
    if ($Z1 eq 1) {
        if ($xy) { print Z "$xy\n"; };
    } else {
        print Z "-2\n";
        print Z "$xy\n";
        print Z "$timestep\n";
    };
    close(Z);
};

sub USE_COORDINATES_FROM_DATA_FILE {
 $found=0;
 $timestep=0; 
 $xy=0;
 $xz=0;
 $yz=0;
 $boxx=$xhi-$xlo;
 $boxy=$yhi-$ylo;
 $boxz=$zhi-$zlo;
 foreach $c (1 .. $chains) { $N[$c]=0; }; 
 until ($found eq 1) {
  $line=<DATA>; $line=STRIP($line);
  if ($line=~/ atoms/) { ($origatoms,$key)=split(/ /,$line); };
  if ($line=~/Atoms/) { $found=1; }; 
  if ($line=~/ xy xz yz/) { ($xy,$xz,$yz)=split(/ /,$line); $xy+=0; $xz+=0; $yz+=0; };  # added 13 nov 2024
 };
 $line=<DATA>;
 foreach $i (1 .. $origatoms) {
  $line=<DATA>; $line=STRIP($line); @tmp=split(/ /,$line);
  if ($#tmp eq 5) {
   ($id,$mymol,$mytype,$myx,$myy,$myz)=split(/ /,$line);
  } elsif ($#tmp eq 6) {
   ($id,$mymol,$mytype,$mycharge,$myx,$myy,$myz)=split(/ /,$line);
  } elsif ($#tmp eq 8) {
   ($id,$mymol,$mytype,$myx,$myy,$myz,$myix,$myiy,$myiz)=split(/ /,$line);
  } elsif ($#tmp eq 9) {
   ($id,$mymol,$mytype,$mycharge,$myx,$myy,$myz,$myix,$myiy,$myiz)=split(/ /,$line);
  } else {
   print "unrecognized data format ($#tmp+1 columns)\n"; exit;
  };  
  if ($newid[$id]) {
   $ID=$newid[$id];
   $MOL[$ID]=$newmol[$id];
   $TYPE[$ID]=$mytype;
   $x[$ID]=$myx;
   $y[$ID]=$myy; 
   $z[$ID]=$myz; 
   $N[$MOL[$ID]]+=1;
  }; 
 }; 
 $ATOMS=0; 
 $Ns=""; foreach $i (1 .. $chains) { $Ns.="$N[$i] "; $ATOMS+=$N[$i]; }; $Ns=~s/ $//;
 APPEND_TO_Z1_FORMATTED; 
};

sub INSPECT_DATA_FILE {
    $found=0; 
    until ($found eq 1) { 
        $line=<DATA>; $line=STRIP($line); 
        if ($line=~/ atoms/) { ($atoms,$key)=split(/ /,$line); };
        if ($line=~/ bonds/) { ($bonds,$key)=split(/ /,$line); };
        if ($line=~/ atom types/) { ($atomtypes,$key)=split(/ /,$line); };  
        if ($line=~/ bond types/) { ($bondtypes,$key)=split(/ /,$line); };
        if ($line=~/ xlo/)   { ($xlo,$xhi,$key)=split(/ /,$line); };
        if ($line=~/ ylo/)   { ($ylo,$yhi,$key)=split(/ /,$line); };
        if ($line=~/ zlo/)   { ($zlo,$zhi,$key)=split(/ /,$line); };
        if ($line=~/ xy xz yz/) { ($xy,$xz,$yz)=split(/ /,$line); $xy+=0; $xz+=0; $yz+=0; };
        if ($line=~/^Masses/) { $found=1; $line=<DATA>; };
    };
    if ($verbose) { print "$atoms atoms\n"; };
    if ($verbose) { print "$bonds bonds\n"; };
    foreach $at (1 .. $atomtypes) { 
        $line=<DATA>; $line=STRIP($line); ($aid,$masstype[$at])=split(/ /,$line); 
        $masstype[$at]=int(0.5+$masstype[$at]); 
        if ($verbose) { print "type $at mass $masstype[$at]\n"; };
    };
    $found=0; until ($found eq 1) { $line=<DATA>; if ($line=~/Atoms/) { $found=1; }; }; $line=<DATA>;
    $maxid=0;
    foreach $i (1 .. $atoms) {
        $line=<DATA>; $line=STRIP($line); @tmp=split(/ /,$line); 
        if ($#tmp eq 5) {
            ($id,$mymol,$mytype,$myx,$myy,$myz)=split(/ /,$line); 
        } elsif ($#tmp eq 6) {
            ($id,$mymol,$mytype,$mycharge,$myx,$myy,$myz)=split(/ /,$line);
        } elsif ($#tmp eq 8) {
            ($id,$mymol,$mytype,$myx,$myy,$myz,$myix,$myiy,$myiz)=split(/ /,$line);
        } elsif ($#tmp eq 9) {
            ($id,$mymol,$mytype,$mycharge,$myx,$myy,$myz,$myix,$myiy,$myiz)=split(/ /,$line);
        } else {
            print "unrecognized data format ($#tmp+1 columns)\n"; exit; 
        }; 
        $q[$id]=$mycharge;
        $mol[$id]=$mymol;
        $type[$id]=$mytype; 
        $mass[$id]=$masstype[$mytype];
         if (($ignoreH eq 1)&(abs($mass[$id]-1) < 0.01)) { $active[$id]=0; } else { $active[$id]=1; }; 
        if ($mytype ~~ @ignoreType) { $active[$id]=0; }; 
        if ($id>$maxid) { $maxid=$id; }; 
    }; 
    $found=0; until ($found eq 1) { $line=<DATA>; if ($line=~/Bonds/) { $found=1; }; }; $line=<DATA>;
    foreach $i (1 .. $bonds) { 
        $line=<DATA>; $line=STRIP($line); 
        ($bid,$btype,$b1,$b2)=split(/ /,$line); 
        if (($active[$b1])&($active[$b2])) { 
            $conn[$b1].="$b2 ";
            $conn[$b2].="$b1 ";
        }; 
    }; 
    foreach $i (1 .. $maxid) { $conn[$i]=~s/ $//; }; 
    foreach $i (1 .. $maxid) { @tmp=split(/ /,$conn[$i]); $conns[$i]=$#tmp+1; };  
    if ($verbose) { foreach $i (1 .. $maxid) { if ($conns[$i] eq 0) { print "disconnected atom [id=$i] ignored\n"; }; }; }; 
    $deleted_bonds=0;
    if ($accept_branched eq 1) {
        foreach $i (1 .. $maxid) { 
            if ($conns[$i] > 2) { 
                $to_be_deleted = $conns[$i]-2; 
                print "branched: id $i has $conns[$i] connections $conn[$i] - cut - "; 
                @tmp1=split(/ /,$conn[$i]); $conn[$i]=""; foreach $itmp1 (0 .. $#tmp1-$to_be_deleted) { $conn[$i].="$tmp1[$itmp1] "; }; 
                $conns[$i]-=$to_be_deleted; 
                foreach $to_be (1 .. $to_be_deleted) { 
                    $ii = $tmp1[$#tmp1+1-$to_be]; 
                    @tmp2=split(/ /,$conn[$ii]); $conn[$ii]=""; foreach $itmp2 (0 .. $#tmp2) { if ($tmp2[$itmp2] eq $i) { } else { $conn[$ii].="$tmp2[$itmp2] "; }; }; 
                    $conns[$ii]-=1; 
                    print "@conn[$i] ($conns[$i] connections)\n";
                    $deleted_bonds+=1;
                }; 
            }; 
        }; 
        print "-branched option: $deleted_bonds deleted bonds\n"; 
    } else { 
        foreach $i (1 .. $maxid) { if ($conns[$i]>2) { $ERROR =<<EOF;
    ******* Branched structure! 
	1) If your system carries H atoms, restart with the additional option: -ignore_H
	2) additional info: conn[$i]=$conn[$i]
	3) convert your chain into a linear chain, ie erase atoms unnecessary for Z1 analysis
EOF
                USAGE; 
            }; 
        }; 
    };

    # find chains and erase holes in ids and mols: create newid[..] and newmol[..]
    $chains=0;
    foreach $i (1 .. $maxid) { $conn[$i]=~s/ $//;
        @tmp=split(/ /,$conn[$i]);
        if ($#tmp eq 0) { 
            @tmp2=split(/ /,$conn[$tmp[0]]); 		# added 9 oct 2019
            if ($#tmp2 eq 0) {				# added 9 oct 2019
                if ($verbose) { print "pair $tmp[0] $tmp2[0] is a dimer\n"; }; 
                if ($ignore_dumbbells) { } else {
                    $chains+=1; $endids[$#endids+1]=$i;
                }; 
            } elsif ($#tmp2 eq -1) {
                print "single atom [id=$id] not part of any chain detected and ignored. BUG-E2. contact mk.\n"; die; 
            } else { 
                $chains+=1; $endids[$#endids+1]=$i; 
            }; 
        };
    };
    if (($#endids %2) eq 1) { } else { print "odd number of terminal atoms. BUG.\n"; USAGE; }; 
    $chains/=2;
    if ($verbose) { print "$chains chains\n"; };
    $chain=0; $novelid=0;
    foreach $endid (@endids) {
        if ($deactivated[$endid]) { } else {
            $chain+=1; $bead=1;
            $id=$endid; $id=~s/ $//;
            $OUT1="$chain";
            $novelid+=1; 
            $newid[$id]=$novelid;
            $newmol[$id]=$chain; 
            @ids=split(/ /,$conn[$id]);
            while ($#ids eq 0) {
                $bead+=1;
                $novelid+=1; 
                $newid[$ids[0]]=$novelid;
                $newmol[$ids[0]]=$chain;
                @tmp2=split(/ /,$conn[$ids[0]]); 
                if ($tmp2[0] eq $id) { 		# added 9 oct 2019
                    $conn[$ids[0]]=$tmp2[1]; 
                } elsif ($tmp2[1] eq $id) {
                    $conn[$ids[0]]=$tmp2[0]; 
                } else { 
                    die 'BUG-E1. contact mk'; 
                };
                $id=$ids[0];
                @ids=split(/ /,$conn[$id]);
            };
            # if ($bead eq 1) { print "chain with single bead. BUG\n"; die; };
            if ($deactivated[$id]) { $ERROR="$id was already deactivated! STOP. contact mk.\n"; USAGE; };
            $deactivated[$endid]=1; 		
            $deactivated[$id]=1;
        };
    };
    $atoms=$novelid; 
    $bonds=$atoms-$chains;
    if ($chain eq $chains) { } else { $ERROR = "format error [$chain] [$chains] [$#endids]"; USAGE; }; 
    if ($ignoreH eq 1) { 
        if ($verbose) { print "after ignoring H atoms:\n"; };
        if ($verbose) { print "$atoms atoms left\n"; };
        if ($verbose) { print "$bonds bonds left\n"; };
    }; 
}; 

sub OPEN_DUMP_FILE { open(DUMP,"<$dumpfile"); $snapshot=0; };
sub CLOSE_DUMP_FILE { close(DUMP); CREATE_ID_CONVERSION_TABLE; }; 

sub OPEN_Z1_FORMATTED {
 `rm -f config.Z1`; 
};

sub HEADER { print "timestep chains N Ree Lpp Z app bpp Lpp2 NeCK NeMK NeCC NeMC\n"; };

sub READ_SNAPSHOT_USING_INFO_FROM_DATA {
 $snapshot+=1; 
 # print "snapshot $snapshot\n";
 $found=0;
 until ($found eq 1) { $line=<DUMP>; $line=STRIP($line);
  if ($line=~/ITEM: TIMESTEP/) { $timestep=<DUMP>; chomp $timestep; $timestep+=0; }; 
  if ($line=~/ITEM: NUMBER OF ATOMS/) { 
   $dumpatoms=<DUMP>; chomp $dumpatoms; $dumpatoms+=0; 
   $line=<DUMP>; $line=STRIP($line);
   if ($line=~/ITEM: BOX BOUNDS xy/) { 
    $line=<DUMP>; $line=STRIP($line); ($xlobound,$xhibound,$xy)=split(/ /,$line);
    $line=<DUMP>; $line=STRIP($line); ($ylobound,$yhibound,$xz)=split(/ /,$line);
    $line=<DUMP>; $line=STRIP($line); ($zlobound,$zhibound,$yz)=split(/ /,$line);
    $xlo=$xlobound+min(0.0,$xy);
    $xhi=$xhibound-max(0,$xy);
   } elsif ($line=~/ITEM: BOX BOUNDS/) {
     $line=<DUMP>; $line=STRIP($line); ($xlo,$xhi)=split(/ /,$line); $xy=0; 
     $line=<DUMP>; $line=STRIP($line); ($ylo,$yhi)=split(/ /,$line); $xz=0; 
     $line=<DUMP>; $line=STRIP($line); ($zlo,$zhi)=split(/ /,$line); $yz=0; 
   } else {
    $ERROR="unrecognized BOX BOUNDS format\n"; USAGE; 
   };
   $boxx=$xhi-$xlo;
   $boxy=$yhi-$ylo;
   $boxz=$zhi-$zlo;
   $line=<DUMP>; $line=STRIP($line);
   if ($line=~/ITEM: ATOMS/) {
    $REQ=0; @tmp=split(/ /,$line); 
    foreach $i (0 .. $#tmp-2) {
     if ($verbose) { print "column $i: $tmp[$i+2]\n"; };
     if ("$tmp[$i+2]" eq "id")  { $REQ+=1; $col_id=$i; }; 
     if ("$tmp[$i+2]" eq "mol") { $REQ+=1; $col_mol=$i; };
     if ("$tmp[$i+2]" eq "type"){ $REQ+=1; $col_type=$i; }; 
     if ("$tmp[$i+2]" eq "xs")  { $REQ+=1; $col_x=$i; $scalex=$boxx; };
     if ("$tmp[$i+2]" eq "ys")  { $REQ+=1; $col_y=$i; $scaley=$boxy; };
     if ("$tmp[$i+2]" eq "zs")  { $REQ+=1; $col_z=$i; $scalez=$boxz; };
     if ("$tmp[$i+2]" eq "xu")  { $REQ+=1; $col_x=$i; $scalex=1.0; };
     if ("$tmp[$i+2]" eq "yu")  { $REQ+=1; $col_y=$i; $scaley=1.0; };
     if ("$tmp[$i+2]" eq "zu")  { $REQ+=1; $col_z=$i; $scalez=1.0; };
     if ("$tmp[$i+2]" eq "x")   { $REQ+=1; $col_x=$i; $scalex=1.0; };
     if ("$tmp[$i+2]" eq "y")   { $REQ+=1; $col_y=$i; $scaley=1.0; };
     if ("$tmp[$i+2]" eq "z")   { $REQ+=1; $col_z=$i; $scalez=1.0; };
    };
    if ($verbose) { print "id in column $col_id\n"; };
    if ($verbose) { print "x in column $col_x\n"; };
    if ($verbose) { print "y in column $col_y\n"; };
    if ($verbose) { print "z in column $col_z\n"; };
    $found=1; $ATOMS=0; 
    foreach $j (1 .. $chains) { $N[$j]=0; }; 
    foreach $j (1 .. $dumpatoms) {
     $line=<DUMP>; $line=STRIP($line);
     @tmp=split(/ /,$line); 
     $id=$tmp[$col_id]; 
     if ($newid[$id]) { 
      $ID=$newid[$id];
      $MOL[$ID]=$newmol[$id]; 
      $TYPE[$ID]=$type[$id];
      $x[$ID]=$tmp[$col_x]*$scalex;
      $y[$ID]=$tmp[$col_y]*$scaley;
      $z[$ID]=$tmp[$col_z]*$scalez; 
      $N[$MOL[$ID]]+=1; 
      if ($ID>$ATOMS) { $ATOMS=$ID; }; 
     }; 
    };
    if ($verbose) { print "$ATOMS atoms\n"; };
    $Ns=""; foreach $i (1 .. $chains) { $Ns.="$N[$i] "; }; $Ns=~s/ $//;
    if (($snapshot>=$from)&($snapshot<=$to)&((($snapshot-$from)%$each)eq 0)) { APPEND_TO_Z1_FORMATTED; };
   } else { 
    $ERROR="missing ITEM: ATOMS line\n"; USAGE;
   };
  }; 
 };
 return(0); 
};

sub READ_SNAPSHOT_WITHOUT_INFO_FROM_DATA {
 $snapshot+=1; 
 $found=0;
 until ($found eq 1) { $line=<DUMP>; $line=STRIP($line);
  if ($line=~/ITEM: TIMESTEP/) { $timestep=<DUMP>; chomp $timestep; $timestep+=0; };
  if ($line=~/ITEM: NUMBER OF ATOMS/) {
   $dumpatoms=<DUMP>; chomp $dumpatoms; $dumpatoms+=0;
   $line=<DUMP>; $line=STRIP($line);
   if ($line=~/ITEM: BOX BOUNDS xy/) {
    $line=<DUMP>; $line=STRIP($line); ($xlobound,$xhibound,$xy)=split(/ /,$line);
    $line=<DUMP>; $line=STRIP($line); ($ylobound,$yhibound,$xz)=split(/ /,$line);
    $line=<DUMP>; $line=STRIP($line); ($zlobound,$zhibound,$yz)=split(/ /,$line);
    $xlo=$xlobound+min(0.0,$xy);
    $xhi=$xhibound-max(0,$xy);
    $ylo=$ylobound;
    $yhi=$yhibound;
    $zlo=$zlobound;
    $zhi=$zhibound;
   } elsif ($line=~/ITEM: BOX BOUNDS/) {
     $line=<DUMP>; $line=STRIP($line); ($xlo,$xhi)=split(/ /,$line); $xy=0;
     $line=<DUMP>; $line=STRIP($line); ($ylo,$yhi)=split(/ /,$line); $xz=0;
     $line=<DUMP>; $line=STRIP($line); ($zlo,$zhi)=split(/ /,$line); $yz=0;
   } else {
    $ERROR="unrecognized BOX BOUNDS format\n"; USAGE;
   };
   $boxx=$xhi-$xlo;
   $boxy=$yhi-$ylo;
   $boxz=$zhi-$zlo;
   $line=<DUMP>; $line=STRIP($line);
   if ($line=~/ITEM: ATOMS/) {
    $REQ=0; @tmp=split(/ /,$line);
    foreach $i (0 .. $#tmp-2) {
     if ($verbose) { print "column $i: $tmp[$i+2]\n"; };
     if ("$tmp[$i+2]" eq "id")  { $REQ+=1; $col_id=$i; };
     if ("$tmp[$i+2]" eq "mol") { $REQ+=1; $col_mol=$i; };
     if ("$tmp[$i+2]" eq "type"){ $REQ+=1; $col_type=$i; };
     if ("$tmp[$i+2]" eq "xs")  { $REQ+=1; $col_x=$i; $scalex=$boxx; };
     if ("$tmp[$i+2]" eq "ys")  { $REQ+=1; $col_y=$i; $scaley=$boxy; };
     if ("$tmp[$i+2]" eq "zs")  { $REQ+=1; $col_z=$i; $scalez=$boxz; };
     if ("$tmp[$i+2]" eq "xu")  { $REQ+=1; $col_x=$i; $scalex=1.0; };
     if ("$tmp[$i+2]" eq "yu")  { $REQ+=1; $col_y=$i; $scaley=1.0; };
     if ("$tmp[$i+2]" eq "zu")  { $REQ+=1; $col_z=$i; $scalez=1.0; };
     if ("$tmp[$i+2]" eq "x")   { $REQ+=1; $col_x=$i; $scalex=1.0; };
     if ("$tmp[$i+2]" eq "y")   { $REQ+=1; $col_y=$i; $scaley=1.0; };
     if ("$tmp[$i+2]" eq "z")   { $REQ+=1; $col_z=$i; $scalez=1.0; };
    };
    if ($verbose) { print "id in column $col_id\n"; };
    if ($verbose) { print "x in column $col_x\n"; };
    if ($verbose) { print "y in column $col_y\n"; };
    if ($verbose) { print "z in column $col_z\n"; };
    $found=1; $ATOMS=$dumpatoms; $chains=0; $lastmol=0; @N=(); 
    # foreach $j (1 .. $chains) { $N[$j]=0; };
    foreach $j (1 .. $dumpatoms) {
     $line=<DUMP>; $line=STRIP($line);
     @tmp=split(/ /,$line);
     $id=$tmp[$col_id];
     if ($id eq $j) { } else { print "ERROR: atoms not saved using dump_modify sort id\n"; exit; }; 
     $mol=$tmp[$col_mol]; 
     if ($mol<$lastmol)   { print "ERROR: chains not saved with increasing mol number\n"; exit; }; 
     if ($mol>$lastmol+1) { print "ERROR: chains not saved with increasing mol number\n"; exit; }; 
     $lastmol=$mol;
     if ($N[$mol]) { $N[$mol]+=1; } else { $N[$mol]=1; $chains+=1; }; 
     $ID=$j; 
     $MOL[$ID]=$mol;
     $TYPE[$ID]=$tmp[$col_type]; 
     $x[$ID]=$tmp[$col_x]*$scalex;
     $y[$ID]=$tmp[$col_y]*$scaley;
     $z[$ID]=$tmp[$col_z]*$scalez;
    };
    if ($verbose) { print "$ATOMS atoms\n"; };
    $Ns=""; foreach $i (1 .. $chains) { $Ns.="$N[$i] "; }; $Ns=~s/ $//;
    if (($snapshot>=$from)&($snapshot<=$to)&((($snapshot-$from)%$each)eq 0)) { APPEND_TO_Z1_FORMATTED; };
   } else {
    $ERROR="missing ITEM: ATOMS line\n"; USAGE;
   };
  };
 };
 return(0);
};

sub HANDLE_XML_FILE {   # added 4 july 2022
    open(XML,"<$xmlfile"); 
    while (!eof(XML)) {
        $line=<XML>; $line=STRIP($line); 
        if ($line=~/<box /) {
            $line=~s/\"/ /g; 
            @tmp=split(/ /,$line); 
            foreach $k (0 .. $#tmp-1) {
                if ($tmp[$k] eq "lx=") { $boxx=$tmp[$k+1]; }; 
                if ($tmp[$k] eq "ly=") { $boxy=$tmp[$k+1]; };
                if ($tmp[$k] eq "lz=") { $boxz=$tmp[$k+1]; };
            }; 
        }; 
        if ($line=~/<position /) { 
            $line=~s/\"/ /g; 
            @tmp=split(/ /,$line);
            foreach $k (0 .. $#tmp-1) {
                if ($tmp[$k] eq "num=") { $atoms=$tmp[$k+1]; }; 
            }; 
            foreach $k (0 .. $atoms-1) {
                $line=<XML>; $line=STRIP($line);
                ($x[$k],$y[$k],$z[$k])=split(/ /,$line); 
            }; 
        }; 
        if ($line=~/<bond /) {
            @conn=(); 
            $line=~s/\"/ /g;
            @tmp=split(/ /,$line);
            foreach $k (0 .. $#tmp-1) {
                if ($tmp[$k] eq "num=") { $bonds=$tmp[$k+1]; };
            };
            foreach $k (0 .. $atoms-1) { $conns[$k]=0; }; 
            foreach $k (1 .. $bonds) {
                $line=<XML>; $line=STRIP($line);
                ($btype,$b1[$k],$b2[$k])=split(/ /,$line);
                $conn[$b1[$k]].=" $b2[$k]"; $conns[$b1[$k]]+=1; if ($conns[$b1[$k]]>2) { die 'branched xml not implemented'; }; 
                $conn[$b2[$k]].=" $b1[$k]"; $conns[$b2[$k]]+=1; if ($conns[$b2[$k]]>2) { die 'branched xml not implemented'; };
            };
            foreach $k (0 .. $atoms-1) { $conn[$k]=~s/^ //; }; 
            $chains=0; 
            $Z1part2="";
            foreach $k (0 .. $atoms-1) {
                @tmp=split(/ /,$conn[$k]); 
                if ($#tmp eq 0) { 
                    $chains+=1;  
                    $N[$chains]=1; 
                    $Z1part2.="$x[$k] $y[$k] $z[$k]\n";
                    $last = $k; 
                    $next = $tmp[0]; 
                    while ($next>-1) { 
                        $N[$chains]+=1; 
                        $Z1part2.="$x[$next] $y[$next] $z[$next]\n";
                        @tmp=split(/ /,$conn[$next]);
                        if ($#tmp eq 0) { 
                            $last=$next;
                            $next=-1; 
                            $conn[$last] = ""; 
                        } elsif ($tmp[0] eq $last) { 
                            $last=$next;
                            $next=$tmp[1]; 
                        } elsif ($tmp[1] eq $last) {
                            $last=$next;
                            $next=$tmp[0]; 
                        } else {
                            print "[k=$k] [tmp=@tmp] [$#tmp] [$tmp[0]] [$tmp[1]]\n";
                            die 'ERROR non-linear';
                        };    
                    };  
                }; 
            }; 
            $Z1part1="$chains\n$boxx $boxy $boxz\n"; foreach $k (1 .. $chains) { $Z1part1.="$N[$k] "; }; 
            open(Z1,">config.Z1"); 
            print Z1 "$Z1part1\n$Z1part2"; 
            close(Z1);
        }; 
    }; 
    close(XML);
    print "box $boxx $boxy $boxz with $atoms atoms and $bonds bonds\n"; 
}; 

PARSE; 
if ($verbose) { print "REQDATA $REQDATA REQDUMP $REQDUMP verbose $verbose\n"; }; 
if (($REQDATA eq 1)&&($REQDUMP eq 1)) { 
 OPEN_DATA_FILE;
 INSPECT_DATA_FILE; 
 CLOSE_DATA_FILE; 
 HEADER; 
 OPEN_DUMP_FILE;
 OPEN_Z1_FORMATTED; 
 until (eof(DUMP)||($snapshot>$to)) { READ_SNAPSHOT_USING_INFO_FROM_DATA; }; 
 CLOSE_DUMP_FILE; 
} elsif ($REQDUMP eq 1) {
 HEADER;
 OPEN_DUMP_FILE;
 OPEN_Z1_FORMATTED;
 until ((eof(DUMP))||($snapshot>$to)) { READ_SNAPSHOT_WITHOUT_INFO_FROM_DATA; };
 CLOSE_DUMP_FILE;
} elsif ($REQDATA eq 1) { 
 OPEN_DATA_FILE;		
 INSPECT_DATA_FILE;		
 CLOSE_DATA_FILE;		
 OPEN_DATA_FILE;                
 OPEN_Z1_FORMATTED;		
 USE_COORDINATES_FROM_DATA_FILE; 
 CLOSE_DATA_FILE; 
 OPTIONAL_WRITE_TRUECHAINS_TO_DUMP;
} elsif ($xml eq 1) { 
 HANDLE_XML_FILE; 
};
