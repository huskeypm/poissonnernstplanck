# Poisson-Nernst-Planck equation..

from __future__ import division
from fenics import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import math

#parameters 
nm=1e-9
length = 30*nm # length of nanopore
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
        return near(x[2],-RevH) and on_boundary

class top_boundary(SubDomain):
	def inside(self,x,on_boundary):
		return near(x[2],length+RevH) and on_boundary

meshfile = "/home/AD/bsu233/labscripts/poissonnernstplanck/contantSigma/onepore/3Dpore.xml"
mesh = Mesh(meshfile)
#mesh = UnitSquareMesh(400,400)
#mesh.coordinates()[:] = mesh.coordinates()[:]*length
#define element
P1 = FiniteElement('P',tetrahedron,1) 
P = MixedElement([P1,P1,P1,P1,P1])
V = FunctionSpace(mesh,P)


# Define test functions
#v_1, v_2, v_3 = TestFunctions(V)

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


#define pH regulated surface charge i?
#chhh = Function(V.sub(2).collapse())
#sigmas = Expression("-F*Gamma*((pow(10,-pKa+3) - pow(10,-pKb-3)*ch*ch)/(pow(10,-pKa+3) +ch+ pow(10,-pKb-3)*ch*ch))",F=F,Gamma=Gamma,ch=chhh,pKa=pKa,pKb=pKb,degree=1)

cmm, cpp, chp, cohp, vv = TestFunctions(V)

# flux of two ions

Jm = -Dm*(grad(cm) + charge*zm/(kb*T)*cm*grad(v))
Jp = -Dp*(grad(cp) + charge*zp/(kb*T)*cp*grad(v))
Jh = -Dh*(grad(ch) + charge*zh/(kb*T)*ch*grad(v))
Joh = -Doh*(grad(coh) + charge*zoh/(kb*T)*coh*grad(v))


# get the form 
#  -d Jm = 0 ==> Jm.cmm.dx = 0
aJm = inner(Jm,grad(cmm))*dx
aJp = inner(Jp,grad(cpp))*dx
aJh = inner(Jh,grad(chp))*dx
aJoh = inner(Joh,grad(cohp))*dx

aPoissonL = inner(eps0*epsr*grad(v),grad(vv))*dx


subdomain = MeshFunction("size_t",mesh,2) #since it is 3D, so the third parapeter is 3 -1 = 2

subdomain.set_all(0)
#inboundary = inner_boundary()
pore_S = poresurface()
bottomboundary = bottom_boundary()
topboundary = top_boundary()
#leftboundary = left_boundary()
#rightboundary = right_boundary()

#-------------------------------------
# TODO: Apply surface charge density as NeumanBC
# Now works as DBC using Grahame euqaiton to convert sigma to psi
#
#facet_domains = FacetFunction("size_t",mesh,0)
#facet_domains.set_all(0)
#bottomboundary.mark(facet_domains,1)
#dS=Measure('dS',subdomain_data=facet_domains)
#ds=Measure('ds',subdomain_data=facet_domains)



bottomboundary.mark(subdomain,1) #mark the boundary
topboundary.mark(subdomain,2)
pore_S.mark(subdomain,3)
#rightboundary.mark(subdomain,4)
#inboundary.mark(subdomain,5)

bc1 = DirichletBC(V.sub(4),0,subdomain,1)
bc2 = DirichletBC(V.sub(4),0,subdomain,2)

#---------------------------------------
# assigin boundary condition for K+ and Cl-
bc3 = DirichletBC(V.sub(0),cm0,subdomain,1)
bc4 = DirichletBC(V.sub(0),cm0,subdomain,2)

bc5 = DirichletBC(V.sub(1),cp0,subdomain,1)
bc6 = DirichletBC(V.sub(1),cp0,subdomain,2)
# assign boundary condition for H+ and OH-

bc7 = DirichletBC(V.sub(2),ch0,subdomain,1)
bc8 = DirichletBC(V.sub(2),ch0,subdomain,2)

bc9 = DirichletBC(V.sub(3),coh0,subdomain,1)
bc10 = DirichletBC(V.sub(3),coh0,subdomain,2)



#TT=1
#VV=V.sub(0).collapse()
#VV=FunctionSpace(mesh,P1)
#chh = interpolate(Constant(ch0),VV)

#while TT<10:
#  TT += 1
sigmaS = -20e-3


# convert surface charge density to potential
# use Grahame euqation
#  sigma = sqrt(8*e0*er*kT)sinh(e*v0/2kT){[Na]inf + [Ca2+]_inf*(2+exp(-e*v0/kT)}^0.5
# at 25^o (T = 298k) : sigma = 0.117*sqrt([NaCl])*sinh(v0/51.4)  ----for 1:1 electrolyte and [NaCl] is in M
phi0 = math.asinh((sigmaS/(0.117*(cm0/1000)**0.5)))*51.4/1000

bcx = DirichletBC(V.sub(4),Constant(phi0),subdomain,3) # electric potential on the surface of sphere

bcc = [bc1,bc2,bc3,bc4,bc5,bc6,bc7,bc8,bc9,bc10,bcx]
#sigmaS = Expression("-F*Gamma*((pow(10,-pKa+3) - pow(10,-pKb -3)*chh*chh)/(pow(10,-pKa+3) +chh + pow(10,-pKb -3)*chh*chh))",F=F,Gamma=Gamma,chh=chh,pKa=pKa,pKb=pKb,degree=1)
aPoissonR = F*(cm*zm + cp*zp+ch*zh +coh*zoh)*vv*dx #+ sigmaS*vv*ds(1) #+ avg(inner(sigmaS,vv))*dS(1)

L = aPoissonR  # lhs and rhs are only suitable for bilinear form ???

a = aJm + aJp + aPoissonL + aJh + aJoh
FF = a - L
J = derivative(FF, u)
problem = NonlinearVariationalProblem(FF, u, bcs=bcc,J=J)
solver = NonlinearVariationalSolver(problem)
#solver.parameters["newton_solver"]["linear_solver"] = "gmres"
#solver.parameters["newton_solver"]["preconditioner"] = "ilu"
solver.solve()

c1,c2,c3,c4,v1 = u.split(True)
print c3.vector()[100]
print "--"
print c1.vector()[100]

#  chh.assign(c3)

V1file = File("c1.pvd")
V1file << c1
V2file = File("c2.pvd")
V2file << c2
V3file = File("c3.pvd")
V3file << c3
V4file = File("c4.pvd")
V4file << c4
