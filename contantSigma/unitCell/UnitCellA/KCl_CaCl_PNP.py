#---------------------
# Five ions are present: K Ca Cl H OH
# Constant surface charge density is appiled on nanopore wall
# 

from __future__ import division
from fenics import *
from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import math

#parameters 
#----------------Domain Size--------------

nm=1e-9
spacing = 5.2*nm # the spacing between nanopores
length = 90*nm # length of nanopore
radius = 5.4*nm # radius of nanopore

RevL= spacing + 2*radius #length of reservoir
RevH=20*nm # depth of reservoir
#--------------------------------------------------
sigmaS = -40e-3 #C/m^2  --- Surface charge density 
cKCl = 100 #mol/m^3 == 1mM ----Bulk [KCl]
cCaCl2 = 100 # Bulk [CaCl2] ---Bulk [CaCl2]

ck0 = cKCl # initial K+ concentration
cca0 = cCaCl2 #initial Ca2+ concentration
ccl0 = 2*cCaCl2 + cKCl #initial Cl- concentration
zca = 2 # Ca2+ 
zcl = -1 # Cl
zk = 1 # K

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

#-- Boundary definition---
class poresurface_for_one_pore(SubDomain):
    def inside(self,x,on_boundary):
        indomain_x = x[2] >= DOLFIN_EPS and x[2] <= length - DOLFIN_EPS
        result = indomain_x and on_boundary
        return result

class poresurface_for_cellA(SubDomain):
    def inside(self,x,on_boundary):
        indomain_x = ((x[0] - RevL/2)**2 + (x[1] - RevL/2)**2) <= (RevL/2)**2 + (RevL/2 - radius)**2
        indomain_y = x[2] >= DOLFIN_EPS and x[2] <= length - DOLFIN_EPS
        result = indomain_x and indomain_y and on_boundary
        return result
 
class poresurface_for_cellB(SubDomain):
    def inside(self,x,on_boundary):
        indomain_x = x[0] > DOLFIN_EPS and x[0] < RevL - DOLFIN_EPS
        indomain_y = x[1] > DOLFIN_EPS and x[1] < RevL - DOLFIN_EPS
        #indomain_x = ((x[0] - RevL/2)**2 + (x[1] - RevL/2)**2) <= (RevL/2)**2 + (radius)**2 
        indomain_z = x[2] >= DOLFIN_EPS and x[2] <= length - DOLFIN_EPS
        result = indomain_x and indomain_y and indomain_z and on_boundary
        return result

class bottom_boundary(SubDomain):
    def inside(self,x,on_boundary):
        return near(x[2],-RevH) and on_boundary

class top_boundary(SubDomain):
	def inside(self,x,on_boundary):
		return near(x[2],length+RevH) and on_boundary

meshfile = "/home/AD/bsu233/labscripts/poissonnernstplanck/contantSigma/unitCell/UnitCellA/UnitCellA.xml"
mesh = Mesh(meshfile)
subdomain = MeshFunction("size_t",mesh,2)
subdomain.set_all(0)

pore_S = poresurface_for_cellA()
bottomboundary = bottom_boundary()
topboundary = top_boundary()


bottomboundary.mark(subdomain,1) #mark the boundary
topboundary.mark(subdomain,2)
pore_S.mark(subdomain,3)


# Mark Facet-- will be used for NeumannBC and normal flux calculations
facet_domains = FacetFunction("size_t",mesh)
facet_domains.set_all(0)
bottomboundary.mark(facet_domains,1)
topboundary.mark(facet_domains,2)
dS=Measure('dS',subdomain_data=facet_domains) # Exterior surface integration
ds=Measure('ds',subdomain_data=facet_domains) # Interior surface integration


#efine MixedFunctionSpace--
#--- for dolfin 1.6.0 (grimm)
P1 = FunctionSpace(mesh,"CG",1)
P2 = FunctionSpace(mesh,"CG",1)
P3 = FunctionSpace(mesh,"CG",1)
P4 = FunctionSpace(mesh,"CG",1)
P5 = FunctionSpace(mesh,"CG",1)
P6 = FunctionSpace(mesh,"CG",1)
V = MixedFunctionSpace([P1,P2,P3,P4,P5,P6])

#--- for dolfin 1.6.2 (kafka)
#P1 = FiniteElement('P',tetrahedron,1)
#P = MixedElement([P1,P1,P1,P1,P1])
#V = FunctionSpace(mesh,P)

# assign initial values 
# two ways to do that:
u = Function(V)
u.interpolate(Constant((ck0,cca0,ccl0,ch0,coh0,0)))

#another way to assigin initial values:
# 
# newSpace = V.sub(0).collapse()
# v0 = interpolate("expression",newSpace) --> interpolate is used for expression
# cca0 = Function(newSpace)
# cca0.vector()[:] = cca00*numpy.exp(-v0.vector()[:]/(kb*T))

# ccl0 = Function(newSpace)
# ccl0.vector()[:] = ccl00*numpy.exp(-v0.vector()[:]/(kb*T))

# assign(u,[v0,cca0,ccl0])

ck, cca, ccl, ch, coh, v = split(u)
ckm, ccam, cclp, chp, cohp, vv = TestFunctions(V)

## TODO: apply a pH-regulated surface charge density
##
#define pH regulated surface charge i?
#chhh = Function(V.sub(2).collapse())
#sigmas = Expression("-F*Gamma*((pow(10,-pKa+3) - pow(10,-pKb-3)*ch*ch)/(pow(10,-pKa+3) +ch+ pow(10,-pKb-3)*ch*ch))",F=F,Gamma=Gamma,ch=chhh,pKa=pKa,pKb=pKb,degree=1)


#-----------------------------------
# Write the variational forms
#----------------------------------
# flux of ions
Jk = -Dk*(grad(ck) + charge*zk/(kb*T)*ck*grad(v))
Jm = -Dca*(grad(cca) + charge*zca/(kb*T)*cca*grad(v))
Jp = -Dcl*(grad(ccl) + charge*zcl/(kb*T)*ccl*grad(v))
Jh = -Dh*(grad(ch) + charge*zh/(kb*T)*ch*grad(v))
Joh = -Doh*(grad(coh) + charge*zoh/(kb*T)*coh*grad(v))

# get the LHS 
#  -div(J) = 0 ==> J*grad(testf)*dx = 0
aJk = inner(Jk,grad(ckm))*dx
aJm = inner(Jm,grad(ccam))*dx
aJp = inner(Jp,grad(cclp))*dx
aJh = inner(Jh,grad(chp))*dx
aJoh = inner(Joh,grad(cohp))*dx

#--LHS and RHS of poisson equation
aPoissonL = inner(eps0*epsr*grad(v),grad(vv))*dx
aPoissonR = F*(ck*zk + cca*zca + ccl*zcl+ch*zh +coh*zoh)*vv*dx #+ sigmaS*vv*ds(1) #+ avg(inner(sigmaS,vv))*dS(1)


FF = aJk + aJm + aJp + aJh + aJoh + aPoissonL - aPoissonR
J = derivative(FF, u)


#--------Boundary Conditions--------------------------
#-- Ground Potential at the two ends of reservoirs
bc1 = DirichletBC(V.sub(5),0,subdomain,1)
bc2 = DirichletBC(V.sub(5),0,subdomain,2)

#---------------------------------------
# assigin boundary condition for Ca2+ and Cl-
bc3 = DirichletBC(V.sub(0),0,subdomain,1) #----> Assign a 0 K+ at the botton reservor
bc4 = DirichletBC(V.sub(0),ck0,subdomain,2)

bc13 = DirichletBC(V.sub(1),0,subdomain,1) #----> Assign a 0 Ca2+ at the botton reservor
bc14 = DirichletBC(V.sub(1),cca0,subdomain,2)

bc5 = DirichletBC(V.sub(2),ccl0,subdomain,1)
bc6 = DirichletBC(V.sub(2),ccl0,subdomain,2)

# assign boundary condition for H+ and OH-
bc7 = DirichletBC(V.sub(3),ch0,subdomain,1)
bc8 = DirichletBC(V.sub(3),ch0,subdomain,2)
bc9 = DirichletBC(V.sub(4),coh0,subdomain,1)
bc10 = DirichletBC(V.sub(4),coh0,subdomain,2)

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
print "phi0",phi0
bcx = DirichletBC(V.sub(5),Constant(phi0/1000),subdomain,3) # electric potential on the surface of sphere


#-------------------------------------
# TODO: Apply surface charge density as NeumanBC
# Now works as DBC using Grahame euqaiton to convert sigma to psi
#

bcc = [bc1,bc2,bc3,bc4,bc5,bc6,bc7,bc8,bc9,bc10,bcx,bc13,bc14]

#-------------------
# Solve the problem
#--------------------
problem = NonlinearVariationalProblem(FF, u, bcs=bcc,J=J)
solver = NonlinearVariationalSolver(problem)
#solver.parameters["newton_solver"]["linear_solver"] = "gmres"
#solver.parameters["newton_solver"]["preconditioner"] = "ilu"
solver.solve()

ck_u,cca_u,ccl_u,ch_u,coh_u,v_u = u.split(True) 

v1file = File("cca.pvd")
v1file << cca_u

v2file = File("ck.pvd")
v2file << ck_u
# Measure [K+] flux at the reservoir boundary
TT = cca_u.function_space()
degree = TT.ufl_element().degree()
# define new functionspace to store flux
W = VectorFunctionSpace(mesh,'P',degree)
flux_ca = project(grad(cca_u)*Constant(Dca),W)
flux_ck = project(grad(ck_u)*Constant(Dk),W)

n=FacetNormal(mesh)
area = assemble(Constant(1.0)*ds(2,domain=mesh))
flux_top_ca = assemble(dot(flux_ca,n)*ds(2))
flux_top_ck = assemble(dot(flux_ck,n)*ds(2))
avgf_ca = flux_top_ca/area
avgf_ck = flux_top_ck/area
Deff_ca = avgf_ca*(length + 2*RevH)/cca0/Dca
Deff_ck = avgf_ck*(length + 2*RevH)/ck0/Dk

print "Average Flux of Ca2+ is: ",avgf_ca
print "Effective Diffusion constant is of Ca2+ is: ",Deff_ca

print "Average Flux of K+ is: ",avgf_ck
print "Effective Diffusion constant is of K+ is: ",Deff_ck
V2file = File("flux_ca.pvd")
V2file << flux_ca

V3file = File("flux_ck.pvd")
V3file << flux_ck

