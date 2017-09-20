# Poisson-Nernst-Planck equation..

from __future__ import division
from fenics import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import math

#parameters 
nm=1e-9
length = 90*nm # length of nanopore
radius = 4*nm

RevH=20*nm # depth of reservoir

F = 96485.3365 # C/mol
cm0 = 100 # mol/m^3 = 1mM
cp0 = 100 # mol/m^3 = 1mM
zm = 1 # cation
zp = -1 # anion
eps0 = 8.854187817e-12 # C/V-m
epsr = 78.5
kb = 1.380648e-23
T = 298 #K
charge = 1.60217e-19 #C
Dm = 1.96e-9 # m^2/s
Dp = 2.03e-9 # m^2/s

#add H+ and OH+
pH = 7
ch0 = 10**(-pH+3) #mol/m^3
coh0 = 10**(-11+pH) #mol/m^3
zh=1
zoh=-1
Dh=9.31e-9 #
Doh=5.3e-9 #

pKa = 7 
pKb = 1.9
Gamma = 5e-6 #mol/m^2 

class poresurface(SubDomain):
    def inside(self,x,on_boundary):
        indomain_x = x[2] >= DOLFIN_EPS and x[2] <= length - DOLFIN_EPS
        result = indomain_x and on_boundary
        return result

class bottom_boundary(SubDomain):
    def inside(self,x,on_boundary):
        result = near(x[2],-RevH) and on_boundary
        return result 

class top_boundary(SubDomain):
    def inside(self,x,on_boundary):
        result = near(x[2],length+RevH) and on_boundary
        result = near(x[2],1.1e-7) and on_boundary
        #print x[2]
        #if result: 
        #  print True
        return result 


import numpy as np

## trying to define a boundary that is within a finite distance from the 
# pore surface 
class inside_pore(SubDomain):
  def inside(self, x, on_boundary):
    #  halfway point 
    indomain_x = x[2] >= (2.5e-8) and x[2] <= (3.5e-8)                
    #if indomain_x:
    #  print x
    # midpoint in x,y 
    mp = 0.8e-8 
    rsqd = (x[0]-mp)**2 + (x[1]-mp)**2
    r = np.sqrt(rsqd) 
 
    alongWall = r >= 4e-9
    alongWall = r >-1       
    result = alongWall and indomain_x
    
    return result      



    

#meshfile = "/home/AD/bsu233/labscripts/poissonnernstplanck/contantSigma/onepore/3Dpore.xml"
meshfile = "onepore.xml"
mesh = Mesh(meshfile)
coords = mesh.coordinates()[:]

print "WARNING: I think I should be using a discontinuous gallerkin here, but was getting nan's"
V = FunctionSpace(mesh, 'Lagrange',1) 

# Define test and trial functions
u = TrialFunction(V)
v = TestFunction(V)

# Define normal vector and mesh size
n = FacetNormal(mesh)

# Define the source term f, Dirichlet term u0 and Neumann term g
f = Constant(0)

# Mark facets of the mesh
boundaries = FacetFunction('size_t', mesh, 0)
boundaries.set_all(0)
#inboundary = inner_boundary()
pore_S = poresurface()
bottomboundary = bottom_boundary()
topboundary = top_boundary()

bottomboundary.mark(boundaries,1) #mark the boundary
topboundary.mark(boundaries,2)
pore_S.mark(boundaries,3)

bc1 = DirichletBC(V,1,boundaries,1)
bc2 = DirichletBC(V,0,boundaries,2)
bcc = [bc1,bc2]



# Initialize mesh function for interior domains
domains = CellFunction("size_t", mesh)
domains.set_all(0)
bufferRegion = inside_pore()
bufferRegion.mark(domains, 1)


# Define outer surface measure aware of Dirichlet and Neumann boundaries
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)
dx = Measure('dx', domain=mesh, subdomain_data=domains)     
#plot(domains)
#interactive()

# define buffering 


# Define variational problem
## Give two diffusion constants, one for the non-buffered region, the
## other for the buffered region 
Jm = -Dm*(grad(v))

### Here I am updating the Ca diffusion coefficient with 'bufferRegion'
B_tot = 1e5 # mol/m^3 = 1mM                 
Ks = 1.    # mol/m^3 = 1mM        
print "Not quite correct, had prob. with compiling, dbl check weak form"
## This correction is based on 
## Eqns 14 and 15 of Wagner, J., & Keizer, J. (1994). Effects of rapid buffers on Ca2+ diffusion and Ca2+ oscillations. Biophysical Journal, 67(1), 447. Retrieved from http://www.sciencedirect.com/science/article/pii/S0006349594805004
betaDenom = 1 + B_tot / (1 + cm0/Ks)**2
beta = 1/betaDenom
print "beta", beta
Db = Dm * beta
Jmb= -Db*(grad(v))



#a = dot(Jm, grad(u))*dx
a = dot(Jm, grad(u))*dx(0)
a+= dot(Jmb, grad(u))*dx(1)
#L = v*f*dx - u0*dot(grad(v), n)*ds(1) 
L = v*f*dx

# Compute solution
u = Function(V)
solve(a == L, u,bcs=bcc)


File("buff.pvd")<<u

print u.vector()[100]




