from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

#parameters 
nm=1e-9
length = 80*nm
radius = 10*nm
interdis = 20*nm
F = 96485.3365 # C/mol
c0 = 1 # mol/m^3 = 1mM
z1 = 1 # cation
z2 = -1 # anion
eps0 = 8.854187817e-12 # C^2/Jm
epsr = 78.5
kb = 1.380648e-23
T = 298 #K
charge = 1.60217e-19 #C
Dk = 1.96e-9 # m^2/s
Dcl = 2.03e-9 # m^2/s


class inner_boundary(SubDomain):
    def inside(self,x,on_boundary):
        indomain_x = x[0] > DOLFIN_EPS and x[0] < length - DOLFIN_EPS
        indomain_y = x[1] > DOLFIN_EPS and x[1] < length - DOLFIN_EPS
        result = indomain_x and indomain_y and on_boundary
        return result
class outside_boundary(SubDomain):
    def inside(self,x,on_boundary):
        indomain_x = x[0] > DOLFIN_EPS and x[0] < length - DOLFIN_EPS
        indomain_y = x[1] > DOLFIN_EPS and x[1] < length - DOLFIN_EPS
        within = indomain_x and indomain_y
        if within:
		result = False
        else:
                result = True
        return result and on_boundary

# Now try to do iteration
def solvePBproblem(cation_conc,anion_conc,bcs):
        bcs = bcs
	c1 = Function(V)
	c1.vector()[:] = cation_conc.vector()[:]
	c2 = Function(V)
	c2.vector()[:] = anion_conc.vector()[:]
	u = TrialFunction(V)
        v = TestFunction(V)
	a = inner(grad(u),grad(v))*dx
	L = ((F*c1*z1+F*c2*z2)/(eps0*epsr))*v*dx
        phi = Function(V)
	pbproblem = LinearVariationalProblem(a,L,phi,bcs)
	solver = LinearVariationalSolver(pbproblem)
        solver.parameters["linear_solver"] = "gmres"
        solver.parameters["preconditioner"] = "amg"
	solver.solve()
    	return phi


        
def solveNPproblem(phi,ze,diffu_constant):
        pmf = Function(V) # store electric potential
        pmf.vector()[:] = phi.vector()[:]
        zz = ze
        beta = charge*zz/(kb*T)
        DD=diffu_constant
        prefact = Expression("exp(-beta*pmf)",beta=beta,pmf=pmf)
        f = Constant(0.0)
        u = TrialFunction(V)
        v = TestFunction(V)
        a = DD*prefact*inner(grad(u),grad(v))*dx
        L = f*v*dx
        bc = DirichletBC(V,1*exp(beta*pmf),subdomain,2)
        u = Function(V)
        npproblem = LinearVariationalProblem(a,L,u,bc)
        solver = LinearVariationalSolver(npproblem)
        solver.parameters["linear_solver"] = "gmres"
        solver.parameters["preconditioner"] = "amg"
        solver.solve()
        return project(prefact*u)



meshfile = "/u1/bsu233/FEniCS/PNP/t1.xml"
mesh = Mesh(meshfile)
V = FunctionSpace(mesh,"CG",1)
subdomain = MeshFunction("size_t",mesh,1)
subdomain.set_all(0)

innerboundary = inner_boundary()
innerboundary.mark(subdomain,1) #mark the boundary
outboundary = outside_boundary()
outboundary.mark(subdomain,2)
bcs = []
bc1 = DirichletBC(V,20e-3,subdomain,1) # electric potential on the surface of sphere
bcs.append(bc1)
bc2 = DirichletBC(V,0,subdomain,2)
bcs.append(bc2)

K = 1
cp1 = Function(V)
cp1.vector()[:] = c0
cp2 = Function(V)
cp2.vector()[:] = c0
while K < 100:
        print "Iteration --", K
	phi = solvePBproblem(cp1,cp2,bcs)
        print "Phi",phi.vector()[88] # print electric potentail near sphere surface
        cp1 = solveNPproblem(phi,z1,Dk)
        cp2 = solveNPproblem(phi,z2,Dcl)
        K += 1
        
   
        for i in range(4):
           if K == (i+1)*20 :
               name1 = "phi" + str(K) + ".pvd"
               name2 = "pcon" + str(K) + ".pvd"
               name3 = "ncon" + str(K) + ".pvd"
               pp = File(name1)
               cc = File(name2)
               cc2 = File(name3)
               phi.vector()[:] = phi.vector()[:]
               pp << phi
               cc << cp1
               cc2 << cp2
            
        # calculate flux 
        #==================================
        #J = -Dk*(grad(c1) + c1*beta*grad(phi))
        #flux = assemble(J)
        #print J.vector()[2]
        #print "interation %s " %T
