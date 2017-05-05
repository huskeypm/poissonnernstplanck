# Poisson-Nernst-Planck equation..

from __future__ import division
from fenics import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import math

#parameters 
nm=1e-9
length = 40*nm
radius = 5*nm
interdis = 20*nm
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

#add H+ and OH-
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
#--------------------------------


class inner_boundary(SubDomain):
    def inside(self,x,on_boundary):
        indomain_x = x[0] > DOLFIN_EPS and x[0] < length - DOLFIN_EPS
        indomain_y = x[1] > DOLFIN_EPS and x[1] < length - DOLFIN_EPS
        result = indomain_x and indomain_y and on_boundary
        return result

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
mesh = UnitSquareMesh(400,400)
mesh.coordinates()[:] = mesh.coordinates()[:]*length
#define element
P1 = FiniteElement('P',triangle,1)
P = MixedElement([P1,P1,P1,P1,P1])
V = FunctionSpace(mesh,P)


# assign initial values for v,cm,cp
# two ways to do that:
u = Function(V)
u.interpolate(Constant((cm0,cp0,ch0,coh0,0)))

#another way to assigin initial values:
#
# newSpace = V.sub(0).collapse()
# v0 = interpolate("expression",newSpace) --> interpolate is used for expression
# cm0 = Function(newSpace)
# cm0.vector()[:] = cm00*numpy.exp(-v0.vector()[:]/(kb*T))

# cp0 = Function(newSpace)
# cp0.vector()[:] = cp00*numpy.exp(-v0.vector()[:]/(kb*T))

# assign(u,[v0,cm0,cp0])



cm, cp, ch, coh, v = split(u)
cmm, cpp, chp, cohp, vv = TestFunctions(V)

# flux of ions

Jm = -Dm*(grad(cm) + charge*zm/(kb*T)*cm*grad(v)) # K+
Jp = -Dp*(grad(cp) + charge*zp/(kb*T)*cp*grad(v)) # Cl-
Jh = -Dh*(grad(ch) + charge*zh/(kb*T)*ch*grad(v)) # H+
Joh = -Doh*(grad(coh) + charge*zoh/(kb*T)*coh*grad(v)) #OH-


# get the form 
#  -d Jm = 0 ==> Jm.cmm.dx = 0
aJm = inner(Jm,grad(cmm))*dx
aJp = inner(Jp,grad(cpp))*dx
aJh = inner(Jh,grad(chp))*dx
aJoh = inner(Joh,grad(cohp))*dx

aPoissonL = inner(eps0*epsr*grad(v),grad(vv))*dx


subdomain = MeshFunction("size_t",mesh,1)

subdomain.set_all(0)
#inboundary = inner_boundary()
bottomboundary = bottom_boundary()
topboundary = top_boundary()
leftboundary = left_boundary()
rightboundary = right_boundary()

# For Neumann boundary condition
# facet_domains = FacetFunction("size_t",mesh,0)
# facet_domains.set_all(0)
# bottomboundary.mark(facet_domains,1)
# dS=Measure('dS',subdomain_data=facet_domains)
# ds=Measure('ds',subdomain_data=facet_domains)



bottomboundary.mark(subdomain,1) #mark the boundary
topboundary.mark(subdomain,2)
leftboundary.mark(subdomain,3)
rightboundary.mark(subdomain,4)
#inboundary.mark(subdomain,5)

#-----------------------------------
# DrichletBC
#----------------------------------
bc2 = DirichletBC(V.sub(4),0,subdomain,2)

# assigin boundary condition for K+ and Cl-
bc3 = DirichletBC(V.sub(0),cm0,subdomain,2)

bc6 = DirichletBC(V.sub(1),cp0,subdomain,2)
# assign boundary condition for H+ and OH-

bc9 = DirichletBC(V.sub(2),ch0,subdomain,2)

bc12 = DirichletBC(V.sub(3),coh0,subdomain,2)
#------------------------------------------

#------------------------
# Variational Form
#-----------------------
aPoissonR = F*(cm*zm + cp*zp+ch*zh +coh*zoh)*vv*dx #+ sigmaS*vv*ds(1) #+ avg(inner(sigmaS,vv))*dS(1)
L = aPoissonR  # lhs and rhs are only suitable for bilinear form ???
a = aJm + aJp + aPoissonL + aJh + aJoh
FF = a - L
J = derivative(FF, u)


def Myexpression(Experssion):
	def __init__(self,F):
		self.F = F
	def eval(self,val,x):
             val=math.asinh((-F*Gamma*((10**(-pKa+3) - 10**(-pKb -3)*self.F(x)*self.F(x))/(10**(-pKa+3) +self.F(x) + 10**(-pKb -3)*self.F(x)*self.F(x))) )/(0.117*(cm0/1000)**0.5))*51.4/1000
	

#---------------------------------------------
# Do iterations --> update [H] and sigmas
#---------------------------------------------
TT=1
VV=V.sub(0).collapse()
chh = interpolate(Constant(ch0),VV) # construct a functon to store [H] after each iteration

while TT<10:
  TT += 1

# convert surface charge density to potential
# use Grahame euqation
#  sigma = sqrt(8*e0*er*kT)sinh(e*v0/2kT){[Na]inf + [Ca2+]_inf*(2+exp(-e*v0/kT)}^0.5
# at 25^o (T = 298k) : sigma = 0.117*sqrt([NaCl])*sinh(v0/51.4)  ----for 1:1 electrolyte and [NaCl] is in M
# See detailed describtion of Grahame equation on p308 of Israelachvili's book 
#phi0 = math.asinh((sigmaS/(0.117*(cm0/1000)**0.5)))*51.4/1000

#bc1 = DirichletBC(V.sub(4),Constant(phi0),subdomain,1) # electric potential on the surface of sphere
#  phi0 = Expression("asinh((-F*Gamma*((pow(10,-pKa+3) - pow(10,-pKb -3)*chh*chh)/(pow(10,-pKa+3) +chh + pow(10,-pKb -3)*chh*chh)) \
 #)/(0.117*pow(cm0/1000,0.5)))*51.4/1000",F=F,Gamma=Gamma,chh=chh,pKa=pKa,pKb=pKb,cm0=cm0,degree=1)

  #chh.vector()[:] = math.asinh((-F*Gamma*(10**(-pKa+3) - 10**(-pKb-3)*chh.vector()[:]*chh.vector()[:])/(10**(-pKa+3) +chh.vector()[:] + 10**(-pKb -3)*chh.vector()[:]*chh.vector()[:]))/(0.117*(cm0/1000)**0.5))*51.4/1000

  #print chh.vector()[100]
  bc1 = DirichletBC(V.sub(4),Myexpression(chh),subdomain,1) # electric potential on the surface of sphere

  bcc = [bc2,bc3,bc6,bc9,bc12,bc1]
  problem = NonlinearVariationalProblem(FF, u, bcs=bcc,J=J)
  solver = NonlinearVariationalSolver(problem)
  #solver.parameters["newton_solver"]["linear_solver"] = "gmres"
  #solver.parameters["newton_solver"]["preconditioner"] = "ilu"
  solver.solve()

  c1,c2,c3,c4,v1 = u.split(True)
  chh.assign(c3)

V1file = File("c1.pvd")
V1file << c1
V2file = File("c2.pvd")
V2file << c2
V3file = File("c3.pvd")
V3file << c3
V4file = File("c4.pvd")
V4file << c4
