# README #

Fenics-based poisson-nernst-planck solver


5/5/17
PKH copied over files from /home/AD/bsu233/labscripts/poissonnernstplanck

### What is this repository for? ###

* Quick summary
Code in support of nanoporous paper 
* Version

### How do I get set up? ###

** Summary of set up
Compile .geo files
gmsh -3 UnitCellA.geo 
dolfin-convert UnitCellA.msh UnitCellB.xml 

** Sample run 
mpirun -np 16 python KCl_CaCl_PNP.py 


** WARNING: 
* Mesh files are hard-coded into the scripts
* Note also that the lengths/etc attributes of the mesh are also hardcoded in the scripts

** Computer vision/matched filering techniques 
* are available on bitbucket: https://bitbucket.org/pkhlab/gpu_detect


### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact
