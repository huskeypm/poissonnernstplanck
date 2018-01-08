#################################
#
# Bin Sun's PNP solver-rama --- CaCl2 Conductance Validation Script...
# 
# Revisions
#       10.08.10 inception

##################################

#---------------------
# For unit cell B 
# Four ions are present: K Cl H OH
# Constant surface charge density is appiled on nanopore wall
# 

from __future__ import division
from fenics import *
import numpy as np
from scipy.optimize import fsolve

import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import math
#import cPickle as pickle
Dca = 1.1e-9 #m^2/s --Ca2+
Dcl = 2.03e-9 # m^2/s  --Cl-
zh=1
zoh=-1
Dh=9.31e-9 # ---H+
Doh=5.3e-9 # ---OH-

  #---------------------------
  # Basic constants
  #---------------------------
F = 96485.3365 # C/mol
eps0 = 8.854187817e-12 # C/V-m
epsr = 78.5
kb = 1.380648e-23
T = 298 #K
charge = 1.60217e-19 #C --- electron charge
cKCl = 0

nm=1e-9

cCaCl2 = CCC #mol/m^3 == 1mM ----Bulk [KCl]
phi0 = PHI

spacing = 4*nm # the distance between the line and slit mouth..
length = 5000*nm # length of nanoslit 
height = 18*nm # height of nanoslit
RevL= 400*nm #length of reservoir

meshfile = "/net/share/bsu233/Nanopores/CaCl2conductance/5000nmsig3mc/CaCl2conductance.xml"
pH = 5

  
  #Adding ability to change KCl concentration applied to the top for Figure 4
cca0 = cCaCl2 #initial Ca concentration
ccl0 = 2*cCaCl2 #initial Cl- concentration
   
 #--------------------------------------------------
sigmaS = -3e-3 #C/m^2  --- Surface charge density 
#ck0 = cKCl #initial K+ concentration                SRB -moved to inside function to manipulate cKCl input
#ccl0 = cKCl #initial Cl- concentration
zca = 2 # Ca2+
zcl = -1 # Cl
  
ch0 = 10**(-pH+3) #mol/m^3  initial [H]
coh0 = 10**(-11+pH) #mol/m^3 inital [OH]

#-- Boundary definition---
class Leftboundary(SubDomain):
	  def inside(self,x,on_boundary):
  		indomain_x = x[0] <= -RevL +  DOLFIN_EPS 
  		result = indomain_x and on_boundary
  		return result
  
class Rightboundary(SubDomain):
	  def inside(self,x,on_boundary):
  		indomain_x = x[0] >= RevL + length  -  DOLFIN_EPS 
  		result = indomain_x and on_boundary
  		return result

class slitwall(SubDomain):
 	 def inside(self,x,on_boundary):
  		indomain_x = x[0] >= DOLFIN_EPS and x[0] <= length  -  DOLFIN_EPS 
  		result = indomain_x and on_boundary
  		return result
#
class midline(SubDomain):
        def inside(self,x,on_boundary):
                return x[0] <= length/2 + DOLFIN_EPS and x[0] >= length/2 - DOLFIN_EPS and on_boundary
    
  
mesh = Mesh(meshfile)
subdomain = MeshFunction("size_t",mesh,1)
subdomain.set_all(0)
  
wall = slitwall()
Lboundary = Leftboundary()
Rboundary = Rightboundary()
  
Lboundary.mark(subdomain,1) #mark the boundary
Rboundary.mark(subdomain,2)
wall.mark(subdomain,3)
  
  # Mark Facet-- will be used for NeumannBC and normal flux calculations
facet_domains = FacetFunction("size_t",mesh)
facet_domains.set_all(0)
Midline=midline()
Lboundary.mark(facet_domains,5)
Rboundary.mark(facet_domains,6)
Midline.mark(facet_domains,7)
 
  
dS=Measure('dS',subdomain_data=facet_domains) # Exterior surface integration
ds=Measure('ds',subdomain_data=facet_domains) # Interior surface integration
  
  
  #efine MixedFunctionSpace--
  #--- for dolfin 1.6.0 (grimm)
P1 = FunctionSpace(mesh,"CG",1)
P2 = FunctionSpace(mesh,"CG",1)
P3 = FunctionSpace(mesh,"CG",1)
P4 = FunctionSpace(mesh,"CG",1)
P5 = FunctionSpace(mesh,"CG",1)
V = MixedFunctionSpace([P1,P2,P3,P4,P5])
  
u = Function(V)
u.interpolate(Constant((cca0,ccl0,ch0,coh0,0.0)))
  
  
cca, ccl, ch, coh, v = split(u)
ccam, cclp, chp, cohp, vv = TestFunctions(V)
  
Jm = -Dca*(grad(cca) + charge*zca/(kb*T)*cca*grad(v))
Jp = -Dcl*(grad(ccl) + charge*zcl/(kb*T)*ccl*grad(v))
Jh = -Dh*(grad(ch) + charge*zh/(kb*T)*ch*grad(v))
Joh = -Doh*(grad(coh) + charge*zoh/(kb*T)*coh*grad(v))
  
aJm = inner(Jm,grad(ccam))*dx
aJp = inner(Jp,grad(cclp))*dx
aJh = inner(Jh,grad(chp))*dx
aJoh = inner(Joh,grad(cohp))*dx
  
  #--LHS and RHS of poisson equation
aPoissonL = inner(eps0*epsr*grad(v),grad(vv))*dx
aPoissonR = F*(cca*zca + ccl*zcl+ch*zh +coh*zoh)*vv*dx #+ sigmaS*vv*ds(1) #+ avg(inner(sigmaS,vv))*dS(1)
  
  
FF = aJm + aJp + aJh + aJoh + aPoissonL - aPoissonR
J = derivative(FF, u)
  
  
  #--------Boundary Conditions--------------------------
  #-- Ground Potential at the two ends of reservoirs
bc1 = DirichletBC(V.sub(4),0.05,subdomain,1)
bc2 = DirichletBC(V.sub(4),0,subdomain,2) #currently testing conductance
  
  #---------------------------------------
  # assigin boundary condition for K+ and Cl-
bc3 = DirichletBC(V.sub(0),cca0,subdomain,1) #----> Assign a 0 [K+] at the botton reservor
bc4 = DirichletBC(V.sub(0),cca0,subdomain,2)
bc5 = DirichletBC(V.sub(1),ccl0,subdomain,1)
bc6 = DirichletBC(V.sub(1),ccl0,subdomain,2)
  
  # assign boundary condition for H+ and OH-
bc7 = DirichletBC(V.sub(2),ch0,subdomain,1)
bc8 = DirichletBC(V.sub(2),ch0,subdomain,2)
bc9 = DirichletBC(V.sub(3),coh0,subdomain,1)
bc10 = DirichletBC(V.sub(3),coh0,subdomain,2)
  
  # Now most important: Surface charge density at the nanopore surface
  #----------------------------------------------------------------
  # convert surface charge density to potential
  # use Grahame euqation
  #  sigma = sqrt(8*e0*er*kT)sinh(e*v0/2kT){[Na]inf + [Ca2+]_inf*(2+exp(-e*v0/kT)}^0.5
  # at 25^o (T = 298k) : sigma = 0.117*sqrt([NaCl])*sinh(v0/51.4)  ----for 1:1 electrolyte and [NaCl] is in M
  
def Grahamequation(x):
        FF = 0.117*math.sinh(x/51.4)*(cKCl/1000 + (cCaCl2/1000)*(2+math.exp(-x/25.7)))**0.5 - sigmaS
        return FF

#phi0 = fsolve(Grahamequation,0)[0]
#print phi0
bcx = DirichletBC(V.sub(4),Constant(phi0/1000),subdomain,3) # electric potential on the surface of sphere



  
  #-------------------------------------
  # TODO: Apply surface charge density as NeumanBC
  # Now works as DBC using Grahame euqaiton to convert sigma to psi
  #
  
bcc = [bc1,bc2,bc3,bc4,bc5,bc6,bc7,bc8,bc9,bc10,bcx]
  
  #-------------------
  # Solve the problem
  #--------------------
problem = NonlinearVariationalProblem(FF, u, bcs=bcc,J=J)
solver = NonlinearVariationalSolver(problem)
  #solver.parameters["newton_solver"]["linear_solver"] = "gmres"
  #solver.parameters["newton_solver"]["preconditioner"] = "ilu"
solver.solve()
  
cca_u,ccl_u,ch_u,coh_u,v_u = u.split(True) 
  
v1file = File("cca.pvd")
v1file << cca_u
  
v2file = File("v.pvd")
v2file << v_u


  # Measure [Ca+] flux at the reservoir boundary
n=FacetNormal(mesh)
TT = cca_u.function_space()
degree = TT.ufl_element().degree()
  # define new functionspace to store flux
W = VectorFunctionSpace(mesh,'P',degree)
flux_ca = project(grad(cca_u)*Constant(Dca),W)
flux_cl = project(grad(ccl_u)*Constant(Dcl),W)
flux_h = project(grad(ch_u)*Constant(Dh),W)
flux_oh = project(grad(coh_u)*Constant(Doh),W)

#  total flux at the mid slit
lfca = assemble(dot(flux_ca,n)*ds(5))
rfca = assemble(dot(flux_ca,n)*ds(6))
mfca = assemble(dot(flux_ca,n)*ds(7))

v3file = File("fca.pvd")
v3file << flux_ca

v4file = File("fcl.pvd")
v4file << flux_cl

v5file = File("fh.pvd")
v5file << flux_h

v6file = File("foh.pvd")
v6file << flux_oh
