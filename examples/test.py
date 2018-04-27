#---------------------
# 

from __future__ import division
from fenics import *
import ufl
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import math
import sys

arg1=sys.argv[0]
arg2=np.float(sys.argv[1]) # [KCl] concentration
arg3=np.float(sys.argv[2]) # surface charge density at the bottom boarder of the domain

print arg1, arg2, arg3


#parameters 
#----------------Domain Size--------------

nm=1e-9
length = 100*nm # length of domain
width = 50*nm
gap = 30*nm

#--------------------------------------------------
cKCl = arg2 #mol/m^3 == 1mM ----Bulk [KCl]
ck0 = cKCl #initial K+ concentration
ccl0 = cKCl #initial Cl- concentration
zk = 1 # K+ 
zcl = -1 # Cl
sigmaS = arg3 # C/m^2
Dk = 1.96e-9 # m^2/s --Diffusion constant for K+
Dcl = 2.03e-9 # m^2/s  --Cl-

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

#-- Boundary definition---

class bottomboundary(SubDomain):
    def inside(self,x,on_boundary):
        range = x[0] >= length*0.25 and x[0] <= length*0.75
        return (near(x[1],0.5*(width-gap)) or near(x[1],0.5*(width+gap))) and range and on_boundary

class topboundary(SubDomain):
	def inside(self,x,on_boundary):
                range = x[0] >= length*0.25 -DOLFIN_EPS and x[0] <= length*0.75 - DOLFIN_EPS
		return near(x[1],length) and range and on_boundary

class leftboundary(SubDomain):
	def inside(self,x,on_boundary):
		return near(x[0],0) and on_boundary

class rightboundary(SubDomain):
	def inside(self,x,on_boundary):
		return near(x[0],length) and on_boundary

meshfile = "/home/AD/sbh239/pnptest/mesh10.xml"
mesh = Mesh(meshfile)
subdomain = MeshFunction("size_t",mesh,1)
subdomain.set_all(0)

botboundary = bottomboundary()
topboundary = topboundary()
lboundary = leftboundary()
rboundary = rightboundary()


botboundary.mark(subdomain,1) #mark the boundary
topboundary.mark(subdomain,2)
lboundary.mark(subdomain,3)
rboundary.mark(subdomain,4)




#efine MixedFunctionSpace--
#--- for dolfin 1.6.0 (grimm)
P = FunctionSpace(mesh,"P",1)
P1 = FunctionSpace(mesh,"P",1)
P2 = FunctionSpace(mesh,"P",1)
P3 = FunctionSpace(mesh,"P",1)
V = MixedFunctionSpace([P1,P2,P3])

u = Function(V)
u.interpolate(Constant((ck0,ccl0,0)))


ck, ccl, v = split(u)
ckm, cclp, vv = TestFunctions(V)

Jm = -Dk*(grad(ck) + charge*zk/(kb*T)*ck*grad(v))
Jp = -Dcl*(grad(ccl) + charge*zcl/(kb*T)*ccl*grad(v))

aJm = inner(Jm,grad(ckm))*dx
aJp = inner(Jp,grad(cclp))*dx

aPoissonL = inner(eps0*epsr*grad(v),grad(vv))*dx
aPoissonR = F*(ck*zk + ccl*zcl)*vv*dx




FF = aJm + aJp + aPoissonL - aPoissonR
J = derivative(FF, u)


#--------Boundary Conditions--------------------------
#-- concentrations at the top boundary
bc1 = DirichletBC(V.sub(0),ck0,subdomain,3)
bc11 = DirichletBC(V.sub(0),0,subdomain,4)

bc2 = DirichletBC(V.sub(1),ccl0,subdomain,3)
bc22 = DirichletBC(V.sub(1),0,subdomain,4)

# Now most important: Surface charge density at bottom boundary
#----------------------------------------------------------------
# convert surface charge density to potential
# use Grahame euqation
#  sigma = sqrt(8*e0*er*kT)sinh(e*v0/2kT){[Na]inf + [Ca2+]_inf*(2+exp(-e*v0/kT)}^0.5
# at 25^o (T = 298k) : sigma = 0.117*sqrt([NaCl])*sinh(v0/51.4)  ----for 1:1 electrolyte and [NaCl] is in M
phi0 = math.asinh((sigmaS/(0.117*(ck0/1000)**0.5)))*51.4/1000
print phi0
bcx = DirichletBC(V.sub(2),Constant(phi0),subdomain,1) # electric potential at the bottom
bcx1 = DirichletBC(V.sub(2),0,subdomain,3) # electric potential at the bottom
bcx2 = DirichletBC(V.sub(2),0,subdomain,4) # electric potential at the bottom


bcc = [bc1,bc11,bc2,bc22,bcx,bcx1,bcx2]

#-------------------
# Solve the problem
#--------------------
problem = NonlinearVariationalProblem(FF, u, bcs=bcc,J=J)
#solver = NonlinearVariationalSolver(problem)
c_low = Constant(0.0)
c_up = Constant(3000.0)
lower = Function(V)
upper = Function(V)
ninfty = Function(P) ; ninfty.vector()[:] = -np.infty
pinfty = Function(P) ; pinfty.vector()[:] = np.infty
fa = FunctionAssigner(V,[P,P,P])
fa.assign(lower,[interpolate(c_low,P),interpolate(c_low,P),ninfty])
fa.assign(upper,[interpolate(c_up,P),interpolate(c_up,P),pinfty])
problem.set_bounds(lower,upper)
solver = NonlinearVariationalSolver(problem)
snes_solver_parameters = {"nonlinear_solver": "snes",
                          "snes_solver": {"linear_solver": "lu",
                                          "maximum_iterations": 20,
                                          "report": True,
                                          "error_on_nonconvergence": False}}

solver.parameters.update(snes_solver_parameters)
info(solver.parameters, True)

(iter, converged) = solver.solve()


#solver.parameters["newton_solver"]["preconditioner"] = "ilu"
#solver.solve()

ck_u,ccl_u,v_u = u.split(True) 

TT = ck_u.function_space()
degree = TT.ufl_element().degree()
W = VectorFunctionSpace(mesh,'P',degree)
fluxck = project(grad(ck_u)*Constant(-Dk),W)
fluxccl = project(grad(ccl_u)*Constant(-Dcl),W)



# save concentration of K+
v1file = File("ck.pvd")
v1file << ck_u

# save concentration of Cl-
v1file = File("ccl.pvd")
v1file << ccl_u

#save K+ and Cl- flux
v1file = File("ck_flux.pvd")
v1file << fluxck

v1file = File("cl_flux.pvd")
v1file << fluxccl

# save electric potential v
v1file = File("v.pvd")
v1file << v_u
