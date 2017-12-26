##################################
#
# Bin Sun's PNP solver-rama
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
import ufl
from scipy.optimize import fsolve

import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import math
import cPickle as pickle


def wrapf(f):
# takes either Function/uflExpression or float/tuple and wraps with Constant in the latter case
# for easy specification of Dirichlet or Neumann data
    if isinstance(f, GenericFunction) or isinstance(f, ufl.core.expr.Expr):
        return f
    elif isinstance(f, float) or isinstance(f, int) or isinstance(f, tuple):
        return Constant(f)
    else:
        dolfin_error(__name__+".py", "use given function",
            "Function is of unexpected type '%s'." % (type(f),))

def _invert_dict(d):
    if isinstance(d,dict):
        return {i:s for s,t in d.items() for i in t}




Dk = 1.96e-9
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
cKCl = 1
phi0 = -76.7
pH=7.5
nm=1e-9
sigmaS = -0.02 #C/m^2
cCaCl2 = 0 #mol/m^3 == 1mM ----Bulk [KCl]
meshfile = "/home/bsu233/scratch/project/NanoPNP/KClconductance/onepore.xml"

def runPNP(
  meshfile = meshfile,
  cKCl = cKCl,
  phi0 = phi0
  ):

  length = 34*nm
  radius = 5.1*nm
  

  
  #Adding ability to change KCl concentration applied to the top for Figure 4
  #cca0 = cCaCl2 #initial K+ concentration
  #ccl0 = 2*cCaCl2 #initial Cl- concentration
   


  
  ########### SRB-Changed length from 90 to 34 
  
  #RevL= spacing + 2*radius #length of reservoir
  RevL= 40*nm #length of reservoir
  RevH= 100*nm # depth of reservoir
  #--------------------------------------------------
  #sigmaS = -20e-3 #C/m^2  --- Surface charge density 
  #ck0 = cKCl #initial K+ concentration                SRB -moved to inside function to manipulate cKCl input
  zca = 2 # Ca2+
  zk = 1
  zcl = -1 # Cl
  
  ch0 = 10**(-pH+3) #mol/m^3  initial [H]
  coh0 = 10**(-11+pH) #mol/m^3 inital [OH]
  ccl0 = cKCl 
  ck0 = cKCl + coh0 - ch0
  #---------------------------
  # TODO:apply a pH regulated surface charge density
  #pKa = 7 
  #pKb = 1.9
  #pKm = pKm
  #Gamma = 5e-6 #mol/m^2 
  
  #-- Boundary definition---
  class poresurface_for_one_pore(SubDomain):
	  def inside(self,x,on_boundary):
  		indomain_x = x[2] >= DOLFIN_EPS and x[2] <= length - DOLFIN_EPS  #SRB - is this not indomain_z, checking to see if within "pore height"?
  		result = indomain_x and on_boundary
  		return result
  
  class poresurface_for_cellA(SubDomain):
 	 def inside(self,x,on_boundary):
                indomain_x = x[0] >= DOLFIN_EPS and x[0] <= RevL - DOLFIN_EPS
                indomain_y = x[1] >= DOLFIN_EPS and x[1] <= RevL - DOLFIN_EPS
  		indomain_z = x[2] >= DOLFIN_EPS and x[2] <= length - DOLFIN_EPS  #SRB y or z?
  		result = indomain_x and indomain_y and indomain_z and on_boundary
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
  
  
  mesh = Mesh(meshfile)
  subdomain = MeshFunction("size_t",mesh,2)
  subdomain.set_all(0)
  
  pore_S = poresurface_for_one_pore()
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
  pore_S.mark(facet_domains,3)
 
  pore_Top.mark(facet_domains, 3)
  pore_Bottom.mark(facet_domains, 4)
  
  dS=Measure('dS',subdomain_data=facet_domains) # Exterior surface integration
  ds=Measure('ds',subdomain_data=facet_domains) # Interior surface integration
  dx = Measure("dx",subdomain_data=subdomain,domain=mesh) 
  
  #efine MixedFunctionSpace--
  #--- for dolfin 1.6.0 (grimm)
  P = FunctionSpace(mesh,"P",1)
  #V = MixedFunctionSpace([P,P,P,P,P])
  
  #--- for dolfin 1.6.2 (kafka)
  ele = FiniteElement('P',tetrahedron,1)
  mele = MixedElement([ele,ele,ele,ele,ele])
  V = FunctionSpace(mesh,mele)
  
  # assign initial values for v,ck,ccl
  # two ways to do that:
  u = Function(V)
  u.interpolate(Constant((ck0,ccl0,ch0,coh0,0.1)))
  
  #another way to assigin initial values:
  # 
  # newSpace = V.sub(0).collapse()
  # v0 = interpolate("expression",newSpace) --> interpolate is used for expression
  # ck0 = Function(newSpace)
  # ck0.vector()[:] = ck00*numpy.exp(-v0.vector()[:]/(kb*T))
  
  # ccl0 = Function(newSpace)
  # ccl0.vector()[:] = ccl00*numpy.exp(-v0.vector()[:]/(kb*T))
  
  # assign(u,[v0,ck0,ccl0])
  
  ck, ccl, ch, coh, v = split(u)
  ckm, cclp, chp, cohp, vv = TestFunctions(V)
  
  ## TODO: apply a pH-regulated surface charge density
  ##
  #define pH regulated surface charge i?
  #caaa = Function(V.sub(0).collapse())
  #chhh = Function(V.sub(2).collapse())
  #sigmas = Expression("-F*Gamma*((pow(10,-pKa+3) - pow(10,-pKb-3)*ch*ch)/(pow(10,-pKa+3) +ch+ pow(10,-pKb-3)*ch*ch))",F=F,Gamma=Gamma,ch=chhh,pKa=pKa,pKb=pKb,degree=1)
  #sigmaS = Expression("-F*Gamma*((pow(10,-pKa+3) - pow(10,-pKb-3)*ch*ch-pow(10,-pKm-3)*cca)/(pow(10,-pKa+3) +ch+ pow(10,-pKb-3)*ch*ch+pow(10,-pKm-3)*cca))",F=F,Gamma=Gamma,ch=chhh,cca=caaa,pKa=pKa,pKb=pKb,pKm=pKm,degree=1) 

  #-----------------------------------
  # Write the variational forms
  #----------------------------------
  # flux of ions
  Jm = -Dk*(grad(ck) + charge*zk/(kb*T)*ck*grad(v))
  Jp = -Dcl*(grad(ccl) + charge*zcl/(kb*T)*ccl*grad(v))
  Jh = -Dh*(grad(ch) + charge*zh/(kb*T)*ch*grad(v))
  Joh = -Doh*(grad(coh) + charge*zoh/(kb*T)*coh*grad(v))
  
  # get the LHS 
  #  -div(J) = 0 ==> J*grad(testf)*dx = 0
  aJm = inner(Jm,grad(ckm))*dx
  aJp = inner(Jp,grad(cclp))*dx
  aJh = inner(Jh,grad(chp))*dx
  aJoh = inner(Joh,grad(cohp))*dx
  
  #--LHS and RHS of poisson equation
  sigmaS = -0.16
  lscale=1e-9
  boundaryS=lscale*sigmaS/(eps0*epsr)
  aPoissonL = inner(grad(v),grad(vv))*dx
  aPoissonR = F/(eps0*epsr)*(ck*zk + ccl*zcl+ch*zh +coh*zoh)*vv*dx  #+ inner(wrapf(boundaryS),vv)*ds(3) 
  
  
  FF = aJm + aJp + aJh + aJoh + aPoissonL - aPoissonR
  J = derivative(FF, u)
  
  
  #--------Boundary Conditions--------------------------
  #-- Ground Potential at the two ends of reservoirs
  bc1 = DirichletBC(V.sub(4),0,subdomain,1)
  bc2 = DirichletBC(V.sub(4),0.1,subdomain,2)
  
  #---------------------------------------
  # assigin boundary condition for K+ and Cl-
  bc3 = DirichletBC(V.sub(0),ck0,subdomain,1) #----> Assign a 0 [K+] at the botton reservor
  bc4 = DirichletBC(V.sub(0),ck0,subdomain,2)
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
  c_low = Constant(0.0)
  c_up = Constant(3000.0)
  lower = Function(V)
  upper = Function(V)
  ninfty = Function(P) ; ninfty.vector()[:] = -np.infty
  pinfty = Function(P) ; pinfty.vector()[:] = np.infty
  fa = FunctionAssigner(V,[P,P,P,P,P])
  fa.assign(lower,[interpolate(c_low,P),interpolate(c_low,P),interpolate(c_low,P),interpolate(c_low,P),ninfty])
  fa.assign(upper,[interpolate(c_up,P),interpolate(c_up,P),interpolate(c_up,P),interpolate(c_up,P),pinfty])
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

  ck_u,ccl_u,ch_u,coh_u,v_u = u.split(True) 
  
# TT = ck_u.function_space()
# degree = TT.ufl_element().degree()
# define new functionspace to store flux
# W = VectorFunctionSpace(mesh,'P',degree)
# fluxck = project(grad(ck_u)*Constant(-Dk),W)
# fluxccl = project(grad(ccl_u)*Constant(-Dcl),W)
# fluxch = project(grad(ch_u)*Constant(-Dh),W)
# fluxoh = project(grad(coh_u)*Constant(-Doh),W)
 
  R = radius/nm
  v1file = File("ck_%s_%s.pvd"%(str(cKCl),str(sigmaS)))
  v1file << ck_u

  v1file = File("ccl_%s_%s.pvd"%(str(cKCl),str(sigmaS)))
  v1file << ccl_u
# v2file = File("ckflux_%s_%s.pvd"%(str(cKCl),str(sigmaS)))
# v2file << fluxck
# v3file = File("cclflux_%s_%s.pvd"%(str(cKCl),str(sigmaS)))
# v3file << fluxccl
# v4file = File("chflux_%s_%s.pvd"%(str(cKCl),str(sigmaS)))
# v4file << fluxch
# v5file = File("cohflux_%s_%s.pvd"%(str(cKCl),str(sigmaS)))
# v5file << fluxoh
# v6file = File("v_%s_%s.pvd"%(str(cKCl),str(sigmaS)))
# v6file << v_u
  """
  rs = np.linspace(0.1,R,R*1000)
  for i in range(int(R*1000)):
    if i==0:
	voltsTop= v_u(rs[i]*nm,rs[i]*nm,length)*(rs[i]*nm)**2*pi
        voltsBottom=v_u(rs[i]*nm,rs[i]*nm,0)*(rs[i]*nm)**2*pi
    else:
        voltsTop+=v_u(rs[i]*nm,rs[i]*nm,length)*((nm*rs[i])**2 - (nm*rs[i-1])**2)*pi
        voltsBottom+=v_u(rs[i]*nm,rs[i]*nm,0)*((nm*rs[i])**2 - (rs[i-1]*nm)**2)*pi

  voltsTop /= pi*radius**2
  voltsBottom /= pi*radius**2
  """
  
###############
# Commenting out hte above until I get sigmaS working
  #pickle.dump(Results, open("Ca_%s_%s_%s.p"%(str(R),str(length/nm),str(cCaCl2)),"wb"))

#!/usr/bin/env python
import sys
#
# Message printed when program run without arguments 
#
def helpmsg():
  scriptName= sys.argv[0]
  msg="""
Purpose: 
  Bin's solver
 
Usage:
"""
  msg+="mpirun -np #proc python  %s -validation" % (scriptName)
  msg+="""
  
 
Notes:
  Pay attention to the units in the input arguments!

"""
  return msg

#
# MAIN routine executed when launching this script from command line 
#
if __name__ == "__main__":
  import sys
  msg = helpmsg()
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
    if(arg=="-run"):          
      arg1=sys.argv[i+1]
      arg2=np.float(sys.argv[i+2])
      arg3=np.float(sys.argv[i+3])
      #arg7=(np.float(sys.argv[i+7]))
      runPNP(meshfile=arg1, cKCl=arg2, phi0=arg3) #same as at the pickle line, justing doing radii
      quit()  
    if(arg=="-runConcs"):
      arg1=np.float(sys.argv[i+1])
      runPNP(cCaCl2=arg1) #value provided is in mM = mol/m^3
      quit() 




  raise RuntimeError("Arguments not understood")


