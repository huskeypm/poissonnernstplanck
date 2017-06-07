# 
# This code perLs a reaction-diffusion simulation within a cube 
# General reaction-diffusion example for buffered Ca2+ 
# Largely borrowed from Cahn-Hilliard example 
#
# Validation
# -Buffering seems to be correct (set cInit to 45 everywhere, checked that cb final agreed with ExcessBuffer output)
# -anistripic diffusion is at least qualitiatively correct



from __future__ import division
from fenics import *
import numpy as np
from scipy.optimize import fsolve

import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import math
import cPickle as pickle

## My reaction system 
# dc/dt = del D del c - R(c,cB)
# dcB/dt = del D del cB + R(c,cB)
# R(c,CB) = kon * (B - cB)*c - koff*cB


Dca = 0.15e-9 #m^2/s --Ca2+
Dcl = 2.03e-9 # m^2/s  --Cl-
#zh=1
#zoh=-1
#Dh=9.31e-9 # ---H+
#Doh=5.3e-9 # ---OH-

  #---------------------------
  # Basic constants
  #---------------------------
#F = 96485.3365 # C/mol
#eps0 = 8.854187817e-12 # C/V-m
#epsr = 78.5
koff = 1.380648e-23
T = 298 #K
#charge = 1.60217e-19 #C --- electron charge
cKCl = 0
cInit = 5.0
nm=1e-9

cCaCl2 = 100 #mol/m^3 == 1mM ----Bulk [KCl]


def runPNP(
  spacing = 5.2*nm, # the spacing between nanopores
  length = 34*nm, # length of nanopore SRB changed fromn 90 to 34 to match fig 4
  radius = 5.1*nm, # radius of nanopore
#meshfile = "/home/AD/bsu233/labscripts/poissonnernstplanck/contantSigma/unitCell/UnitCellA/UnitCellA.xml"
# PKH 
  meshfile = "/net/share/shared/papers/nanoporous/meshes/UAL34R5-1.xml",
  cCaCl2 = cCaCl2,
  Km = 6.5
  ):

  if radius>2e-7:
    spacing = spacing*nm
    length = length*nm
    radius = radius*nm



  #Adding ability to change KCl concentration applied to the top for Figure 4
  cca0 = cCaCl2 #initial K+ concentration
  ccl0 = 2*cCaCl2 #initial Cl- concentration




  ########### SRB-Changed length from 90 to 34 

  RevL= spacing + 2*radius #length of reservoir
  RevH=20*nm # depth of reservoir
  #--------------------------------------------------
  Gamma = 5e-6 #mol/m^2 
  tubeWallArea = length*(2*radius*pi)
  bT = Gamma*tubeWallArea
  
  cb0 = bT/10
  Kd = (10**(-Km))
  kon = 1
  koff = Kd/kon
  class poresurface_for_one_pore(SubDomain):
          def inside(self,x,on_boundary):
                indomain_x = x[2] >= DOLFIN_EPS and x[2] <= length - DOLFIN_EPS  #SRB - is this not indomain_z, checking to see if within "pore height"?
                result = indomain_x and on_boundary
                return result

  class poresurface_for_cellA(SubDomain):
         def inside(self,x,on_boundary):
                indomain_x = ((x[0] - RevL/2)**2 + (x[1] - RevL/2)**2) <= (RevL/2)**2 + (RevL/2 - radius)**2
                indomain_y = x[2] >= DOLFIN_EPS and x[2] <= length - DOLFIN_EPS  #SRB y or z?
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
  class tubeTop(SubDomain):
        def inside(self,x,on_boundary):
                return near(x[2], length) and on_boundary

  class tubeBottom(SubDomain):
        def inside(self,x,on_boundary):
                return near(x[2], 0) and on_boundary

  class bottom_boundary(SubDomain):
         def inside(self,x,on_boundary):
                return near(x[2],-RevH) and on_boundary

  class top_boundary(SubDomain):
        def inside(self,x,on_boundary):
                return near(x[2],length+RevH) and on_boundary
  class pore_wall(SubDomain):
        def inside(self,x,on_boundary):
          wall=near(((x[0] - RevL/2)**2 + (x[1] - RevL/2)**2),(RevL/2)**2 + (RevL/2 - radius)**2)
          z = x[2] >= DOLFIN_EPS and x[2] <= length - DOLFIN_EPS
          return wall and z 

  mesh = Mesh(meshfile)
  subdomain = MeshFunction("size_t",mesh,2)
  subdomain.set_all(0)

  pore_S = poresurface_for_cellA()
  bottomboundary = bottom_boundary()
  topboundary = top_boundary()
  pore_Top = tubeTop()
  pore_Bottom = tubeBottom()

  bottomboundary.mark(subdomain,1) #mark the boundary
  topboundary.mark(subdomain,2)
  pore_S.mark(subdomain,3)

  # Mark Facet-- will be used for NeumannBC and normal flux calculations
  facet_domains = FacetFunction("size_t",mesh)
  facet_domains.set_all(0)
  bottomboundary.mark(facet_domains,1)
  topboundary.mark(facet_domains,2)

  pore_Top.mark(facet_domains, 3)
  pore_Bottom.mark(facet_domains, 4)


  poreWall = pore_wall()
  poreWall.mark(facet_domains,5)



  dS=Measure('dS',subdomain_data=facet_domains) # Exterior surface integration
  ds=Measure('ds',subdomain_data=facet_domains) # Interior surface integration






  ## Params 
  field=False
  if(field==True):
    Dii = Constant((1.0,1.0,1.0))
    Dii = Constant((5.,0.1,0.1))
    Dij = diag(Dii)
  else:
    Dij = Constant(0.15e-9)
  
  #Kd = KD # [uM] 
    
  
  ## equiliubrium equation 
  def CalccB(c):
    cB = c * bT / (Kd + c) 
    return cB
  
  ## initiation conditions 
  # left hand side of box has Ca2+, Ca2+/Buff present 
  class InitialConditions(Expression):
  
    def eval(self, values, x):
      # edge  
      #if (x[0] < 0.5 and np.linalg.norm(x[1:2]) < 0.5):
      # corner 
      if (np.linalg.norm(x -np.zeros(3) ) < 0.5):
        values[0] = cInit
        values[1] = CalccB( values[0] )
      else:
        values[0] = 0         
        values[1] = 0
    def value_shape(self):
      return (2,)
  
  class InitialConditionsMyo(Expression):
  
    def eval(self, values, x):
      # edge  
      #if (x[0] < 0.5 and np.linalg.norm(x[1:2]) < 0.5):
      # corner 
      if (x[2]> length+RevH-DOLFIN_EPS):
	values[0] = cca0
       #    (np.linalg.norm(x[0:1] -np.array([0,0]) ) < 0.15) or 
        #   (np.linalg.norm(x[0:1] -np.array([xMax,0]) ) < 0.15) or 
         #  (np.linalg.norm(x[0:1] -np.array([xMax,yMax]) ) < 0.15) or 
          # (np.linalg.norm(x[0:1] -np.array([0,yMax]) ) < 0.15) 
        #):
	
        values[1] = 0
	values[2] = 0
      elif (x[2]<(length+DOLFIN_EPS) and near(((x[0] - RevL/2)**2 + (x[1] - RevL/2)**2),(RevL/2)**2 + (RevL/2 - radius)**2) and x[2]>(0-DOLFIN_EPS)):
        values[0] = 0         
        values[1] = 0
	values[2] = Gamma
      else:
	values[0] = 0
	values[1] = 0
	values[2] = 0

    
      #print x
      #print values[0]
    def value_shape(self):
      return (3,)
  
  
  # Class for interfacing with the Newton solver
  class MyEqn(NonlinearProblem):
      def __init__(self, a, L):
          NonlinearProblem.__init__(self)
          self.L = L
          self.a = a
          self.reset_sparsity = True
      def F(self, b, x):
          assemble(self.L, tensor=b)
      def J(self, A, x):
          assemble(self.a, tensor=A, reset_sparsity=self.reset_sparsity)
          self.reset_sparsity = False
  
  
  
  # Define mesh and function space 

  mesh = Mesh(meshfile)
  subdomain = MeshFunction("size_t",mesh,2)
  subdomain.set_all(0)

  pore_S = poresurface_for_cellA()
  bottomboundary = bottom_boundary()
  topboundary = top_boundary()
  pore_Top = tubeTop()
  pore_Bottom = tubeBottom()

  bottomboundary.mark(subdomain,1) #mark the boundary
  topboundary.mark(subdomain,2)
  pore_S.mark(subdomain,3)

  # Mark Facet-- will be used for NeumannBC and normal flux calculations
  facet_domains = FacetFunction("size_t",mesh)
  facet_domains.set_all(0)
  bottomboundary.mark(facet_domains,1)
  topboundary.mark(facet_domains,2)

  pore_Top.mark(facet_domains, 3)
  pore_Bottom.mark(facet_domains, 4)

  dS=Measure('dS',subdomain_data=facet_domains) # Exterior surface integration
  ds=Measure('ds',subdomain_data=facet_domains) # Interior surface integration


  #efine MixedFunctionSpace--
  #--- for dolfin 1.6.0 (grimm)
  P1 = FunctionSpace(mesh,"CG",1)
  P2 = FunctionSpace(mesh,"CG",1)
  P3 = FunctionSpace(mesh,"CG",1)
  #P4 = FunctionSpace(mesh,"CG",1)
  #P5 = FunctionSpace(mesh,"CG",1)
  V = MixedFunctionSpace([P1,P2,P3])#,P4,P5])
  

  du = TrialFunction(V)
  q,v, Help  = TestFunctions(V)  # I dont know what to do with test function for unoccupied sites on the wall



  u = Function(V) # current soln
  u0 = Function(V)
  # split mixed functions
  c,cb,gamma = split(u)
  c0,cb0,gamma0 = split(u0)
  ccam,ccbp,gammaP = TestFunctions(V)
  #c0,cb0 = split(u0)
  
  #u.interpolate(Constant((cca0,0)))
  # Init conts
  #init_cond = InitialConditions()
  init_cond = InitialConditionsMyo()
  print "init", init_cond
  u.interpolate(init_cond)
  u0.interpolate(init_cond)
  

####BCS##########3
   
  bc1 = DirichletBC(V.sub(0),0,subdomain,1) #----> Assign a 0 [K+] at the botton reservor
  bc2 = DirichletBC(V.sub(0),cca0,subdomain,2)
  bc3 = DirichletBC(V.sub(2),Gamma,subdomain,3)
  ## Weak Ls for RHS  
  # See notetaker 121213 notes for details 
  bcs = [bc1,bc2,bc3] 
  # Diffusion
  #if(field==False):
  #  RHS1 = -Dij*inner(grad(c),grad(q))*dx  
  #  RHS2 = -Dij*inner(grad(cb),grad(v))*dx 
  #else:
  RHS1 = inner(Dij*grad(c),grad(q))*dx  #only considering grad(c) not grad(cb) b/c grad(cb) is bound to wall
  #RHS2 = -inner(Dij*grad(cb),grad(v))*dx 
  RHS1 += -(kon*-c*q*(bT-cb)*(gamma-cb/tubeWallArea))*v*ds # add in rxn rates 
  RHS1 = koff*cb*q*ds 					  # rxn rates so KD does somwthing
  RHS2 = c*(gamma-cb/Constant(tubeWallArea))*v*ds - (-koff)*cb*v*ds
  L = RHS1 + RHS2
  # Reaction: b + c --kon--> cb,  
  #           b + c <-koff--- cb
  #R = np.array([
    #[-kon,koff],     # s
    #[kon,-koff]      # p
    #])
  
  
  # operator splitting 
  #opSplit=False
  #if(opSplit==False):
    #RHS1 += (R[0,0]*(bT-cb)*c*q + R[0,1]*cb*q)*dx
    #RHS2 = (R[1,0]*(bT-cb)*c*v + R[1,1]*cb*v)*dx
  
  
  # Add in time dependence 
  # (dc/dt) = RHS  --> c1-c0 - dt * RHS = 0
  #L1 = c*q*dx - c0*q*dx - RHS1 
  #L2 = cb*v*dx - cb0*v*dx - RHS2 
  #L = L1 + L2
  
######## my try ##########

  #RHS1 = (inner(grad(c),grad(ccam)))*dx  
  #RHS2 = L1 -(c*ccam)*(bT-cb)*ds(5, domain=mesh)
  #L3 = form2 + cb*ccam*ds(5, domain=mesh)
  #L4 = form3 + c*(bT-cb)*ccbp*ds(5, domain=mesh) 
  #L = form4 -cb*ccbp*ds(5, domain= mesh)

  #J = derivative(L,u)


  # Compute directional derivative about u in the direction of du (Jacobian)
  # (for Newton iterations) 
  a = derivative(L, u, du)
  
  
  # Create nonlinear problem and Newton solver
  problem = MyEqn(a, L)
  solver = NewtonSolver("lu")
  solver.parameters["convergence_criterion"] = "incremental"
  solver.parameters["relative_tolerance"] = 1e-6
  
  # Output file
  file = File("output.pvd", "compressed")
  
  # Step in time
  t = 0.0
  T = 25*dt
  while (t < T):
      t += dt
      u0.vector()[:] = u.vector()
      solver.solve(problem, u.vector())
  
      #if(opSplit==True):
        # Check Johans 121213/4 email regarding using dof maps and sympy functions here 
        
        
      #file << (u.split()[0], t)
      file << (u,t)
  
      # check values
      for i,ele in enumerate(split(u)):
        tot = assemble(ele*dx,mesh=mesh)
        vol = assemble(Constant(1.)*dx,mesh=mesh)
        print "Conc(%d) %f " % (i,tot/vol)
  
  #roblem = NonlinearVariationalProblem(L, u, bcs=bcs,J=J)
  #olver = NonlinearVariationalSolver(problem)
  #solver.parameters["newton_solver"]["linear_solver"] = "gmres"
  #solver.parameters["newton_solver"]["preconditioner"] = "ilu"
  #olver.solve()


  
  print "It did something!"



if __name__ == "__main__":
  import sys
#  msg = helpmsg()
  remap = "none"

  if len(sys.argv) < 2:
      raise RuntimeError(msg)

  #fileIn= sys.argv[1]
  #if(len(sys.argv)==3):
  #  1
  #  #print "arg"

  # Loops over each argument in the command line 
  for i,arg in enumerate(sys.argv):
    # calls 'doit' with the next argument following the argument '-validation'
    if(arg=="-validation"):
      runPNP()
      quit()

