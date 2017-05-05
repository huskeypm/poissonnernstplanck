# Poisson-Nernst-Planck equation..

from fenics import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

#parameters 
nm=1e-9
length = 100*nm
radius = 5*nm
interdis = 20*nm
F = 96485.3365 # C/mol
cm0 = 100 # mol/m^3 = 1mM
cp0 = 100 # mol/m^3 = 1mM
zm = 1 # cation
zp = -1 # anion
eps0 = 8.854187817e-12 # C^2/Jm
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
		return near(x[0],100*nm) and on_boundary
class  top_boundary(SubDomain):
	def inside(self,x,on_boundary):
		return near(x[1],100*nm) and on_boundary



meshfile = "/home/AD/bsu233/labscripts/poissonnernstplanck/t1.xml"
mesh = Mesh(meshfile)
#mesh = UnitSquareMesh(800,400)
#mesh.coordinates()[:] = mesh.coordinates()[:]*100*nm
#define element
P1 = FiniteElement('P',triangle,1)
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

Jm = -Dm*(grad(cm) + charge*zm/(kb*T)*grad(v))
Jp = -Dp*(grad(cp) + charge*zp/(kb*T)*grad(v))
Jh = -Dh*(grad(ch) + charge*zh/(kb*T)*grad(v))
Joh = -Doh*(grad(coh) + charge*zoh/(kb*T)*grad(v))


# get the form 
#  -d Jm = 0 ==> Jm.cmm.dx = 0
aJm = inner(Jm,grad(cmm))*dx
aJp = inner(Jp,grad(cpp))*dx
aJh = inner(Jh,grad(chp))*dx
aJoh = inner(Joh,grad(cohp))*dx

aPoissonL = inner(eps0*epsr*grad(v),grad(vv))*dx


subdomain = MeshFunction("size_t",mesh,1)

subdomain.set_all(0)
inboundary = inner_boundary()
bottomboundary = bottom_boundary()
topboundary = top_boundary()
leftboundary = left_boundary()
rightboundary = right_boundary()

facet_domains = FacetFunction("size_t",mesh,0)
facet_domains.set_all(10)
inboundary.mark(facet_domains,11)
ds=Measure('ds')[facet_domains]


bottomboundary.mark(subdomain,1) #mark the boundary
topboundary.mark(subdomain,2)
leftboundary.mark(subdomain,3)
rightboundary.mark(subdomain,4)
inboundary.mark(subdomain,5)

#bc1 = DirichletBC(V.sub(4),-20e-3,subdomain,5) # electric potential on the surface of sphere
bc2 = DirichletBC(V.sub(4),0,subdomain,2)
bcx = DirichletBC(V.sub(4),0,subdomain,1)


#---------------------------------------
# assigin boundary condition for K+ and Cl-
bc3 = DirichletBC(V.sub(0),cm0,subdomain,3)
bc5 = DirichletBC(V.sub(0),cm0,subdomain,4)

bc6 = DirichletBC(V.sub(1),cp0,subdomain,3)
bc8 = DirichletBC(V.sub(1),cp0,subdomain,4)
# assign boundary condition for H+ and OH-

bc9 = DirichletBC(V.sub(2),ch0,subdomain,3)
bc11 = DirichletBC(V.sub(2),ch0,subdomain,4)

bc12 = DirichletBC(V.sub(3),coh0,subdomain,3)
bc14 = DirichletBC(V.sub(3),coh0,subdomain,4)

bcc = [bc2,bc3,bc5,bc6,bc8,bc9,bc11,bc12,bc14,bcx]

eos = 1.0
tol = 1e-5

TT=1
#VV=V.sub(0).collapse()
VV=FunctionSpace(mesh,P1)
chh = interpolate(Constant(ch0),VV)

#while TT<10:
#  TT += 1
sigmaS = -20e-3
#sigmaS = Expression("-F*Gamma*((pow(10,-pKa+3) - pow(10,-pKb -3)*chh*chh)/(pow(10,-pKa+3) +chh + pow(10,-pKb -3)*chh*chh))",F=F,Gamma=Gamma,chh=chh,pKa=pKa,pKb=pKb,degree=1)
aPoissonR = F*(cm*zm + cp*zp+ch*zh +coh*zoh)*vv*dx + (sigmaS/(eps0*epsr))*vv*ds(11)
L = aPoissonR  # lhs and rhs are only suitable for bilinear form ???

a = aJm + aJp + aPoissonL + aJh + aJoh
FF = a - L

solve(FF == 0,u,bcs=bcc)

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
#solver_parameters={'linear_solver':'cg',
#'preconditioner':'ilu'})
#problem = PNPequation(a,L)

#solver = NewtonSolver()
#solver.parameters["linear_solver"] = "lu"
#solver.parameters["convergence_criterion"] = "incremental"
#solver.parameters["relative_tolerance"] = 1e-6
#solver.solve(problem,u,bcs)
