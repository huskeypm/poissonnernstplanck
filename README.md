# README #

Fenics-based poisson-nernst-planck solver


5/5/17
PKH copied over files from /home/AD/bsu233/labscripts/poissonnernstplanck

12/23/2017
Bin Sun added necessary files to reproduce data in the PNP paper

### What is this repository for? ###

* Quick summary
Code in support of nanoporous paper

Two main directories are "Geometries" and "fenics_scripts".

"Geometries" contains all 3D .geo files used for fenics calculations.

"fenics_scripts" contains fenics scripts:

  1) KClonlyPNP.py is used to calculate the KCl ionic conductance in a 3D single nanopre.
     The corresponding .geo file is "3D_pore.geo" in "Geometries" directory.
  2) CaCl2onlyPNP.py is used to calculate the CaCl2 ionic conductance in a 3D channel.
     The corresponding .geo file is "3D_channel.geo" in "Geometries" directory.
  3) CF_buffer.py is used to calculate the CF permeability in the hexgonal unit cell with buffers
     The corresponding .geo file is "hexgonal.geo" in "Geometries" directory
* Version

### How do I get set up? ###

** Summary of set up

Compile .geo files

gmsh -3 3D_pore.geo
 
dolfin-convert 3D_pore.msh 3D_pore.xml 

** Sample run 

mpirun -np 16 python KCl_CaCl_PNP.py 100 -30

"100" means the [KCl] = 100mM

"-30" means the electric potential at the nanopore wall is -30mV

** Data analysis in Paraview

After each run, the flux density, concentration and electric potential data are stored in corresponding ".pvd" files,these files can be opened in Paraview.

All surface integrations in Paraview are accomplished fist by an "Slice" operationa and then by the "Integrate Variables" filter.


** WARNING: 
* Mesh files are hard-coded into the scripts

* Note also that the lengths/etc attributes of the mesh are also hardcoded in the scripts

** Computer vision/matched filering techniques 
* are available on bitbucket: https://bitbucket.org/pkhlab/gpu_detect


### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact


### Notes to PKH###
1) the "pH-regulatesigma", the PNP.py and t1.geo are unnecessary.

2) BS added all latest 3D .geo files in main directory

3) BS added all latest fenics scripts in main directory
