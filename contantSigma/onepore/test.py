#---------------------
# For one single pore
# Four ions are present: Ca Cl H OH
# Constant surface charge density is appiled on nanopore wall
# 

from __future__ import division
from fenics import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.optimize import fsolve
import math

#parameters 
#----------------Domain Size--------------

nm=1e-9
spacing = 5.2*nm # the spacing between nanopores
length = 40*nm # length of nanopore
radius = 5.4*nm # radius of nanopore

RevL= spacing + 2*radius #length of reservoir
RevH=20*nm # depth of reservoir
#--------------------------------------------------
sigmaS = -40e-3 #C/m^2  --- Surface charge density 
cKCl = 0 #mol/m^3 == 1mM ----Bulk [KCl]
cCaCl2 = 100 # Bulk [CaCl2]

cca0 = cCaCl2 #initial Ca2+ concentration
ccl0 = 2*cCaCl2 + cKCl #initial Cl- concentration
zca = 2 # Ca2+ 
zcl = -1 # Cl

pH = 7
ch0 = 10**(-pH+3) #mol/m^3  initial [H]
coh0 = 10**(-11+pH) #mol/m^3 inital [OH]

Dk = 1.96e-9 # m^2/s --Diffusion constant for K+
Dca = 0.15e-9 #m^2/s --Ca2+
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
#-----------------------------
# TODO:apply a pH regulated surface charge density
#pKa = 7 
#pKb = 1.9
#Gamma = 5e-6 #mol/m^2 


#meshfile = "/home/AD/bsu233/labscripts/poissonnernstplancca/contantSigma/unitCell/UnitCellB/UnitCellB.xml"
#mesh = Mesh(meshfile)
class bottom_boundary(SubDomain):
    def inside(self,x,on_boundary):
        return near(x[1],0) and on_boundary

class left_boundary(SubDomain):
        def inside(self,x,on_boundary):
                return near(x[0],0) and on_boundary

class right_boundary(SubDomain):
        def inside(self,x,on_boundary):
                return near(x[0],length) and on_boundary
class  top_boundary(SubDomain):
        def inside(self,x,on_boundary):
                return near(x[1],length) and on_boundary



#meshfile = "/home/AD/bsu233/labscripts/poissonnernstplanck/t1.xml"
#mesh = Mesh(meshfile)
mesh = UnitSquareMesh(300,300)
mesh.coordinates()[:] = mesh.coordinates()[:]*length

subdomain = MeshFunction("size_t",mesh,1)

subdomain.set_all(0)
#inboundary = inner_boundary()
bottomboundary = bottom_boundary()
topboundary = top_boundary()
leftboundary = left_boundary()
rightboundary = right_boundary()



bottomboundary.mark(subdomain,1) #mark the boundary
topboundary.mark(subdomain,2)
leftboundary.mark(subdomain,3)
rightboundary.mark(subdomain,4)



#efine MixedFunctionSpace--
#--- for dolfin 1.6.0 (grimm)
P1 = FunctionSpace(mesh,"CG",1)
P2 = FunctionSpace(mesh,"CG",1)
P3 = FunctionSpace(mesh,"CG",1)
P4 = FunctionSpace(mesh,"CG",1)
P5 = FunctionSpace(mesh,"CG",1)
V = MixedFunctionSpace([P1,P2,P3,P4,P5])

#--- for dolfin 1.6.2 (kafka)
#P1 = FiniteElement('P',tetrahedron,1)
#P = MixedElement([P1,P1,P1,P1,P1])
#V = FunctionSpace(mesh,P)

# assign initial values for v,cca,ccl
# two ways to do that:
u = Function(V)
u.interpolate(Constant((cca0,ccl0,ch0,coh0,0)))

#another way to assigin initial values:
# 
# newSpace = V.sub(0).collapse()
# v0 = interpolate("expression",newSpace) --> interpolate is used for expression
# cca0 = Function(newSpace)
# cca0.vector()[:] = cca00*numpy.exp(-v0.vector()[:]/(kb*T))

# ccl0 = Function(newSpace)
# ccl0.vector()[:] = ccl00*numpy.exp(-v0.vector()[:]/(kb*T))

# assign(u,[v0,cca0,ccl0])

cca, ccl, ch, coh, v = split(u)
ccam, cclp, chp, cohp, vv = TestFunctions(V)

## TODO: apply a pH-regulated surface charge density
##
#define pH regulated surface charge i?
#chhh = Function(V.sub(2).collapse())
#sigmas = Expression("-F*Gamma*((pow(10,-pKa+3) - pow(10,-pKb-3)*ch*ch)/(pow(10,-pKa+3) +ch+ pow(10,-pKb-3)*ch*ch))",F=F,Gamma=Gamma,ch=chhh,pKa=pKa,pKb=pKb,degree=1)


#-----------------------------------
# Write the variational forms
#----------------------------------
# flux of ions
Jm = -Dca*(grad(cca) + charge*zca/(kb*T)*cca*grad(v))
Jp = -Dcl*(grad(ccl) + charge*zcl/(kb*T)*ccl*grad(v))
Jh = -Dh*(grad(ch) + charge*zh/(kb*T)*ch*grad(v))
Joh = -Doh*(grad(coh) + charge*zoh/(kb*T)*coh*grad(v))

# get the LHS 
#  -div(J) = 0 ==> J*grad(testf)*dx = 0
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
bc1 = DirichletBC(V.sub(4),0,subdomain,2)

#---------------------------------------
# assigin boundary condition for Ca2+ and Cl-
bc4 = DirichletBC(V.sub(0),cca0,subdomain,2)
bc6 = DirichletBC(V.sub(1),ccl0,subdomain,2)

# assign boundary condition for H+ and OH-
bc8 = DirichletBC(V.sub(2),ch0,subdomain,2)
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

phi0 = fsolve(Grahamequation,0)[0]
bcx = DirichletBC(V.sub(4),Constant(phi0/1000),subdomain,1) # electric potential on the surface of sphere

#-------------------------------------
# TODO: Apply surface charge density as NeumanBC
# Now works as DBC using Grahame euqaiton to convert sigma to psi
#

bcc = [bc1,bc4,bc6,bc8,bc10,bcx]

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

v1file = File("ccl.pvd")
v1file << ccl_u


