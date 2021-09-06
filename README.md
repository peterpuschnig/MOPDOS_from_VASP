# MOPDOS_from_VASP
This repository contains a fortran code which computes the molecular orbital projected density of states (MOPDOS) based on VASP calculations (WAVECAR-files). The MOPDOS is calculated by projecting the Kohn-Sham orbitals of a "full" system onto the Kohn-Sham orbitals of a fragment of the full system, the "molecular" system. A detailed defintion of the MOPDOS can be found in the following publication:

[1] D. Lüftner, P. Hurdax, G. Koller, M. G. Ramsey, S. Weiß, X. Yang, S. Soubatch, F.S. Tautz, V. Feyer, A. Gottwald and P. Puschnig,
Understanding the photoemission distribution of strongly interacting two-dimensional overlayers,
Phys. Rev. B 96, 125402 (2017).
https://doi.org/10.1103/PhysRevB.96.125402

Please cite publication [1] when using results from this code in your publications.

## Authors
- Peter Puschnig (peter.puschnig@uni-graz.at)
- Daniel Lüftner (daniel.lueftner@hotmail.com)

## Quick-Start

Compilation of program:

    cd src
    gfortran mopdos_spinpol.f90
    mv a.out ../bin/mopdos_spinpol.x

Usage:

    cd ../example/
    ../bin/mopdos_spinpol.x

The program runs for less than a second for this example data (benzene adsorbed on graphene) and produces the text outputfile `MOPDOS.out` containing the MOPDOS data and the log-file `MOPDOS.log`. 
Compare your output with the sample outputfile provided in the `example/output/` folder which also contains the file `MOPDOS.png` displaying the data containd in `MOPDOS.out` as well as the python script `plot_mopdos.py` which reads in the data from `MOPDOS.out` and creates the plot `MOPDOS.png`. 

## Example-data

### VASP data

The `example/VASP` folder contains the results of two VASP caclulations. The subfolder `full` contains the VASP calculation for the full system, here, it is a benzene ring adsorbed on a graphene layer. The subfolder `molecule` contains the VASP calculation for the molecular system, here, it is the benzene ring. Note that it is absolutely necessary that the unit cell dimensions and the plane wave cut-off must be identical in the full and molecular system. Also, the atomic coordintes of the molecular fragment in 'molecule' and 'full' must be identical-

### Computation of the MOPDOS

In order to compute the MOPDOS go to the respective folder,

    cd ../example
    
If necessary, edit the text file `MOPDOS.in` which provides the input data for running the MOPDOS program:
    "VASP/full/WAVECAR"         ! name of WAVECAR file of full calculation
    "VASP/molecule/WAVECAR"     ! name of WAVECAR file of molecule calculation
    4                           ! number of k-points in IBZKPT
    11 15                       ! band range for analysis in molecule 
    1 82                        ! band range for analysis in full system
    -1.5894                     ! Fermi-level of full system (energy of HOMO) in eV
    0.1                         ! sigma value for Gauss broadening in MOPDOS curves
    -15.0    5.00    0.02       ! range of energies for analysis (1 - 2) and energy step (3) [eV]  

Make sure that the VASP output file `IBZKPT` is also located in the same folder as the `MOPDOS.in`. Execute the program  

    ../bin/mopdos_spinpol.x   
    
The resulting text output file `MOPDOS.out` contains the computed MOPDOS curves. You can use the python script 
`plot_mopdos.py` to read in and plot the MOPDOS curves.


## Generating VASP files for your own systems



-----
September 6th, 2021
