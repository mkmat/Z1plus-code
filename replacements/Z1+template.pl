#! INTERPRETER

# (c) 01 june 2022 mk@mat.ethz.ch 
#     11 july 2024 allows for blanks in directory and file names

# Copyright 2022 Martin Kroger
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#   http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


# installation directory
$installation_directory = INSTALLATION-DIRECTORY.use.Z1+install.pl

# working directory
$working_directory = `pwd`; chomp $working_directory; 

$perl = "PERL";
$user = "USER";

sub greenprint { print "\033[1;32;40m$_[0]\033[1;37;40m\n"; }; 
sub redprint   { print "\033[1;31;40m$_[0]\033[1;37;40m\n"; }; 
sub show_directories {
    print "working directory:      $working_directory\n";
    print "installation directory: $installation_directory\n";
};

sub check_directories { 
    $code = "$installation_directory/Z1+.ex";
    if ("$working_directory" eq "$installation_directory") { 
        redprint("do not start Z1+ from within the installation_directory"); exit;
    };
};

sub load_options { @OPTIONS=(); @FILEFORMATS=(); @EXAMPLES=(); 
    $OPTIONS[$#OPTIONS+1] = "-h or -help";
    $OPTIONS[$#OPTIONS+1] = "-c or -clean";
    $OPTIONS[$#OPTIONS+1] = "-l or -log";
    $OPTIONS[$#OPTIONS+1] = "-s or -stats";
    $OPTIONS[$#OPTIONS+1] = "-0 or -selfZ";
    $OPTIONS[$#OPTIONS+1] = "+  or -SP+";
    $OPTIONS[$#OPTIONS+1] = "-t or -PPA";
    $OPTIONS[$#OPTIONS+1] = "-p or -PPA+";
    $OPTIONS[$#OPTIONS+1] = "The following options can be used to pick snapshots from a trajectory";
    $OPTIONS[$#OPTIONS+1] = "-from=snapshot-no";
    $OPTIONS[$#OPTIONS+1] = "-to=snapshot-no"; 
    $OPTIONS[$#OPTIONS+1] = "The following options are available starting from lammps data files containing bond information";
    $OPTIONS[$#OPTIONS+1] = "-ignore_H";
    $OPTIONS[$#OPTIONS+1] = "-ignore_type1=type-no";
    $OPTIONS[$#OPTIONS+1] = "-ignore_type2=type-no";
    if ($user eq "mk") { 
        $OPTIONS[$#OPTIONS+1] = "-d or -debug";
        $OPTIONS[$#OPTIONS+1] = "-m or -movie";
        $OPTIONS[$#OPTIONS+1] = "-v or -visualize";
    };
    $FILEFORMATS[$#FILEFORMATS+1] = "Z1-formatted";
    $FILEFORMATS[$#FILEFORMATS+1] = "lammps.data";
    $FILEFORMATS[$#FILEFORMATS+1] = "lammps.dump";
    $FILEFORMATS[$#FILEFORMATS+1] = "hoomd-blue-galamost.xml"; 
    $EXAMPLES[$#EXAMPLES+1] = "config.Z1";
    $EXAMPLES[$#EXAMPLES+1] = "-clean";
    $EXAMPLES[$#EXAMPLES+1] = "config.data";
    $EXAMPLES[$#EXAMPLES+1] = "-self config.data";
    $EXAMPLES[$#EXAMPLES+1] = "-from=3 -to=10 -ignore_H config.data";
};

sub LOGO {
    greenprint("    _____  ___  ");
    greenprint("   /__  / /__ |    _    ");
    greenprint("     / /    | |  _| |_     (c) 2022 mk@mat.ethz.ch");
    greenprint("    / /_    | | |_   _| ");
    greenprint("   /____/   |_|   |_|   ");
    greenprint("                        ");
};

sub USAGE {
    greenprint("\nusage:");
    print "Z1+ [options] configuration-file\n";
    greenprint("\nlist of options:");
    foreach $option (@OPTIONS) { 
        print "$option\n";
    }; 
    greenprint("\nlist of allowed configuration file formats:\n");
    foreach $fileformat (@FILEFORMATS) { 
        print "$fileformat\n";
    };
    greenprint("\nexamples:");
    foreach $example (@EXAMPLES) { 
        print "Z1+ $example\n";
    };
    exit;
}; 

sub create_config_Z1 { 
    if (-s "$configfile") { } else { redprint("$configfile does not exist"); exit; }; 
    open(fp,"<$configfile"); 
    foreach $i (0 .. 20) { 
        $line = <fp>; chomp $line; 
        if ($line=~/atom types/) { 
            greenprint("\nconverting lammps.data -> config.Z1");
            print "$perl \"$installation_directory/Z1+import-lammps.pl\" -data=\"$configfile\" $arg_from $arg_to $arg_noH $arg_ignore_type1 $arg_ignore_type2\n"; 
            `$perl "$installation_directory/Z1+import-lammps.pl" -data="$configfile" $arg_from $arg_to $arg_noH $arg_ignore_type1 $arg_ignore_type2`; 
            return;
        }; 
        if ($line=~/ITEM: TIMESTEP/) {
            greenprint("\nconverting lammps.dump -> config.Z1");
            `$perl "$installation_directory/Z1+import-lammps.pl" -dump="$configfile" $arg_from $arg_to`; 
            return;
        }; 
        if ($line=~/galamost/) { 
            greenprint("\nconverting hoomd-blue-galamost.xml -> config.Z1"); 
            `$perl "$installation_directory/Z1+import-lammps.pl" -xml="$configfile"`;  
            return;
        }; 
    }; 
    if ($configfile eq "config.Z1") { } else { 
        `cp "$working_directory/$configfile" config.Z1`; 
        greenprint("\nassuming Z1-formatted configuration file");
    }; 
    if (-s "config.Z1") { } else { 
        redprint("Your configuration could not be processed, config.Z1 was not created.");
        USAGE;
    };
};


sub READ_ARGUMENTS { 
    if ($#ARGV eq -1) { USAGE; }; 
    # defaults
    $True       = "True"; 
    $False      = "False"; 
    $configfile = "config.Z1";
    $cleanfiles = $False;
    $basic_SP   = $True;
    $visualize  = $False;
    $movie      = $False;
    $selfZ      = $False;
    $debug      = $False;
    $stats      = $False;
    $log        = $False;
    $PPA        = $False;
    $CPPA       = $False;
    $arg_from   = "";
    $arg_to     = "";
    $arg_noH    = "";
    $arg_ignore_type1 = "";
    $arg_ignore_type2 = "";
    $SMDP       = "Z1+SP.dat";

    foreach $arg (@ARGV) { ($field,$value)=split(/=/,$arg); 
        
        if (($arg eq "-help")||($argv eq "-h")) { USAGE; 
        } elsif (($arg eq "-SP+")      ||($arg eq "+"))  { $basic_SP=$False; 
        } elsif (($arg eq "-clean")    ||($arg eq "-c")) { $cleanfiles=$True; 
        } elsif (($arg eq "-visualize")||($arg eq "-v")) { $visualize=$True;
        } elsif (($arg eq "-movie")    ||($arg eq "-m")) { $movie=$True; 
        } elsif (($arg eq "-selfZ")    ||($arg eq "-0")) { $selfZ=$True; 
        } elsif (($arg eq "-debug")    ||($arg eq "-d")) { $debug=$True; 
        } elsif (($arg eq "-log")      ||($arg eq "-l")) { $log=$True; 
        } elsif (($arg eq "-stats")    ||($arg eq "-s")) { $stats=$True; 
        } elsif (($arg eq "-PPA+")     ||($arg eq "-p")) { $PPA=$True;      $SMDP="PPA+.dat"; 
        } elsif (($arg eq "-PPA")      ||($arg eq "-t")) { $CPPA=$True;     $SMDP="PPA.dat"; 
        } elsif ($field eq "-from")     { $arg_from="$field=$value"; 
        } elsif ($field eq "-to")       { $arg_to="$field=$value"; 
        } elsif ($field eq "-ignore_H") { $arg_noH="-ignore_H"; 
        } elsif ($field eq "-ignore_type1") { $arg_ignore_type1="-ignore_type1";
        } elsif ($field eq "-ignore_type2") { $arg_ignore_type2="-ignore_type2";
        } else { $configfile=$arg;
        }; 
    }; 
    create_config_Z1; 
}; 

sub handle_clean_before_start { @DEL=();
    $DEL[$#DEL+1] = "Z1+SP.dat";
    $DEL[$#DEL+1] = "log.Z1";
    $DEL[$#DEL+1] = "Lpp_values.dat";
    $DEL[$#DEL+1] = "N_values.dat";
    $DEL[$#DEL+1] = "Ree_values.dat";
    $DEL[$#DEL+1] = "Z1+NODES.dat";
    $DEL[$#DEL+1] = "NODES-unfolded.dat";
    $DEL[$#DEL+1] = "Z1+initconfig.dat";
    $DEL[$#DEL+1] = "$SMDP"; 
    if (($PPA eq $False)&&($CPPA eq $False)) { $DEL[$#DEL+1] = "Z_values.dat"; };  
    foreach $filen (@DEL) { 
        if (-s "$filen") { print "[Z1+] removing $filen\n"; `rm -f "$filen"`; }; 
    }; 
    if ($PPA eq $True) { 
        foreach $filen ("PPA+summary.dat","PPA+summary.html") { print "[Z1+] removing $filen\n"; `rm -f "$filen"`; };  
    } elsif ($CPPA eq $True) {
        foreach $filen ("PPA-summary.dat","PPA-summary.html") { print "[Z1+] removing $filen\n"; `rm -f "$filen"`; };
    } else { 
        foreach $filen ("Z1+summary.dat","Z1+summary.html")   { print "[Z1+] removing $filen\n"; `rm -f "$filen"`; };
    }; 
    
    print "[Z1+] finished cleaning\n";
};

sub handle_clean_after_start { `rm -f Z1+parameters`; };

sub handle_clean { 
    handle_clean_before_start;
    handle_clean_after_start;
    if ("$configfile" eq "config.Z1") { } else { print "[Z1+] removing config.Z1\n"; `rm -f config.Z1`; }; 
    print "[Z1+] finished cleaning\n";
}; 

sub create_parfile { 
    if (($calc_PPA eq $True)&&($calc_CPPA eq $True)) { die 'both -p and -t options simultaneously is forbidden'; }; 
    open(f,">Z1+parameters");
    print f "\&parameters\n";
    print f "stats = .$stats.\n";
    print f "logfile = .$log.\n";
    print f "lmax_factor = 1.0\n";
    print f "thickness = 0.002\n";
    print f "self_entanglement = .$selfZ.\n";
    print f "basic_SP = .$basic_SP.\n";
    print f "calc_PPA = .$PPA.\n";
    print f "calc_CPPA = .$CPPA.\n"; 
    if ($user eq "mk") { 
        print f "debug = .$debug.\n";
        print f "user = '$user'\n";
        print f "visualize = .$visualize.\n";
        print f "movie = .$movie.\n";
    }; 
    print f "/\n";
    close(f);
    greenprint("\nZ1+parameters contains:");
    print `cat Z1+parameters`; print "\n";
}; 

sub check_if_files_exist { 
    if (-s "$code") { } else { redprint("STOP. $code is missing"); exit; }; 
}; 

# ------------------------ main -----------------------

LOGO;
show_directories;
load_options; 
check_directories;
check_if_files_exist;
READ_ARGUMENTS;

if ($cleanfiles eq $True) { handle_clean; exit; }; 

create_parfile;
handle_clean_before_start;
greenprint("\n"+code+" launched ..");
if ($log eq $True) { 
 `"$code"`; 
} else {
 open my $cmd_mk, "$code |"; while (<$cmd_mk>) { print "$_"; };
};
handle_clean_after_start;

if (($debug eq $True)&&(-s "_error")) { 
    if (-s "_error") { 
        greenprint("Errors");
        print `cat _error`; 
        open(fp,">_error"); 
    };
}; 

if (-s "Z1+SP.dat") { } else { if ((-s "PPA.dat")||(-s "PPA+.dat")) { } else { redprint("CRASHED. contact mk"); exit; }; }; 

`$perl "$installation_directory/Z1+rearrange.pl" $SMDP`; 
open(f,"<.Z1+cpu+secs");
$cpuseconds = <f>; $cpuseconds+=0;
close(f);
print "Z1+ finished after $cpuseconds seconds\n";
