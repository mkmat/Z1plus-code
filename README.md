# Z1+ code

The Z1+ code creating the shortest multiple disconnected path for the analysis of entanglements in macromolecular systems is available for download [here at mendeley](https://data.mendeley.com/datasets/m425t6xtwr/1)
The related publication describing all features is available for free [here at Comput. Phys. Commun.](https://www.sciencedirect.com/science/article/pii/S0010465522002867?via%3Dihub)

Here we collect questions, answers, and additional scripts that may be useful for Z1+ users. 

## How to extract linear backbones from fully atomistic LAMMPS models 

Question: I am simulating atomistically detailed PMMA via LAMMPS, and have saved both a LAMMPS data file and a LAMMPS dump trajectory. Z1+ crashes as the LAMMPS files carry branched structures (chemical formula below). How to make it work? 

<img src="images/PMMA-chemical-formula.jpg" width="25%">

Answer: I created a script that automatically recognizes and extracts the linear backbones from your LAMMPS data file, and saves the linear conformation as Z1-formatted file config.Z1. The Z1+ code can then be directly applied to config.Z1. If you have both a LAMMPS data file and LAMMPS dump trajectory for the same system, this script creates a Z1-formatted trajectory file. The script is available for download in the scripts folder. Call it via

    perl ./extract-backbone.pl

or 

    perl ./extract-backbone.pl <lammps-data-file>

or

    perl ./extract-backbone.pl <lammps-data-file> <lammps-dump-trajectory>

## How to extract linear backbones from fully atomistic LAMMPS models, if the atomistic model contains non-polymers in addition?

Question posed by Jingqi Zhang in Feb 2024. I have LAMMPS data and dump-trajectories (id mol type xu yu zu) for a system that contains branched polymers as well as individual C60 beads (bead type 4). How to convert the dump-trajectory file to a Z1-trajectory file that contains only the linear backbones of the polymers? Such Z1-trajectory file can be analyzed directly using the Z1+ code, while the LAMMPS dump-trajectory file produces errors. 

Answer: The extract-backbone.pl script had been extended to contain a -ignore-types=<type1,type2,..> option. Call it via 


    perl ./extract-backbone.pl <lammps-data-file> <lammps-dump-trajectory> -ignore-types=4

If your system has more than a single atom type that need to be ignored, such as types 2,4, and 10, use -ignore-types=2,4,10.

## How to produce a Z1-formatted trajectory file from an unsorted LAMMPS dump-trajectory?

A LAMMPS dump-file does not contain information about bonds. Only if the dump-file had been generated using the dump_modify sort id option, and if your bead id is bonded to the adjacent bead id, Z1+ can recognize the chains. The LAMMPS data file, on the other hand, contains bond information. With a LAMMPS data and unsorted LAMMPS dump-trajectory at hand, you can use the extract-backbone.pl script to create a Z1-formatted trajectory file, as the extract-backbone.pl script retrieves information about the bonds from the data-file, and uses it to sort the trajectory file. 

    perl ./extract-backbone.pl <lammps-data-file> <unsorted-lammps-dump-trajectory>

## How to convert a Z1-formatted configuration or trajectory file to LAMMPS-dump-formatted file? 

    perl ./Z1+dump [-unfolded] <Z1-formatted-file>  

creates a LAMMPS-dump file or LAMMPS-dump trajectory file. If the option -unfolded is given, the dump-file contains unfolded coordinates (xu yu zu), otherwise it contains folded (wrapped) coordinates (x,y,z). The dump-file contains two bead types: type 1 (interior beads), type 2 (terminal beads). The script is available for download in the scripts folder. 

## How to merge shortest path and original configuration file into a single data or dump trajectory?

Call

    perl ./Z1+export.pl

to see the options. It creates data or dump files or trajectories for selected (or all) snapshots and assigns bead types 1,2,3 for the original chains, and bead types 4,5,6 for the shortest path. 

## Are there benchmark configurations to test my own implementation of Z1+?

Yes, some of the benchmark configurations treated in the publication are available from the benchmark-configurations directory. 

## The -PPA (and -PPA+) option produces no useful result

This happens if the system is not of standard Kremer-Grest type, with a maximum bond length of 1.5. The PPA option has been added using classical PPA parameters, and can therefore only be applied if the system respects the constraint. We did not invent new PPA parameters to allow for a comparison with classical PPA results, and because results depend on the choice of parameters. If your system has a bond length that exceeds the maximum allowed value 1.5 is seen in this line:

        PPA+ init max bondl (all)     1.54848

If you still want to use the PPA or PPA+ options, you have to scale your box sizes and particle coordinates in your configuration file.

## Z1+ crashes because the the lammps data and/or dump files created by vmd or other software seem to be corrupt.

Z1+ crashes, because the data and/or dump files may not contain the molecule IDs. To heal this problem we offer a script that corrects data and dump files using just one command, and saves the new files with "-corrected" appended to their original names. Call

        perl convert_vmd_data_to_proper_data.pl

to see the description. A typical call is 

        perl convert_vmd_data_to_proper_data.pl -data=MyLammps.data -dump=MyLammps.dump 

## How to cite the Z1+ code?

    M. Kröger, J. D. Dietz, R. S. Hoy and C. Luap,
    The Z1+package: Shortest multiple disconnected path for the analysis of entanglements in macromolecular systems,
    Comput. Phys. Commun. 283 (2023) 108567. DOI 10.1016/j.cpc.2022.108567

or if you are using bibtex:

    @article{Z1+,
     author = {M. Kr\"oger and J. D. Dietz and R. S. Hoy and C. Luap},
     title = {The Z1+package: Shortest multiple disconnected path for the analysis of entanglements in macromolecular systems},
     journal = {Comput. Phys. Commun.},
     volume = {283},
     pages = {108567},
     year = {2023},
     doi = {10.1016/j.cpc.2022.108567}
    }
