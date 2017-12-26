##################################
#
# Bin Sun's PNP solver-rama
# 
# Revisions
#       10.08.10 inception

##################################

#---------------------
# CF Diffusion  CaCl2 + KCl
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
Dcf = 4.14e-10 # m^2/s
Dca= 1.1e-9
zcf = -1
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
IonicS = 100 #background ionic strength
cCF = 1 # concentration of CF
ratio=0.001
phi0 = -25.7
pH=7.4
nm=1e-9
rr = 2 # radius of pore
dd = 1 # thickness of buffering shell
Btot=100 # ratio of [buffer]/[CF]

meshfile = "/home/bsu233/scratch/project/NanoPNP/CFdiffusion/UnitCell/tt.xml"

def runPNP(
  meshfile = meshfile,
  rr = rr,
  dd = dd,
  Btot=Btot,
  ):
  cKCl=IonicS
  
  cCl = cKCl
  cK = cKCl

  length = 90*nm
  radius = rr*nm
  thickness = dd*nm
  #Adding ability to change KCl concentration applied to the top for Figure 4
  #cca0 = cCaCl2 #initial K+ concentration
  #ccl0 = 2*cCaCl2 #initial Cl- concentration
  
  hsp = 6*nm
  vsp = 21*nm

  RevW = hsp + 2*5e-9
  RevH = vsp + 2*5e-9 
  RevD = 40*nm
  #--------------------------------------------------
  #sigmaS = -20e-3 #C/m^2  --- Surface charge density 
  #ck0 = cKCl #initial K+ concentration                SRB -moved to inside function to manipulate cKCl input
  zca = 2 # Ca2+
  zcf = -1
  zk = 1
  zcl = -1 # Cl
  
  ch0 = 10**(-pH+3) #mol/m^3  initial [H]
  coh0 = 10**(-11+pH) #mol/m^3 inital [OH]
  ccl0 = cCl 
  ck0 = cK
  ccf0 = cCF
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
  
  class poresurface_for_cellB(SubDomain):
         def inside(self,x,on_boundary):
                indomain_x = x[0] >= DOLFIN_EPS and x[0] <= RevW - DOLFIN_EPS
                indomain_y = x[1] >= DOLFIN_EPS and x[1] <= RevH - DOLFIN_EPS
  #indomain_x = ((x[0] - RevL/2)**2 + (x[1] - RevL/2)**2) <= (RevL/2)**2 + (radius)**2 
                indomain_z = x[2] >= DOLFIN_EPS and x[2] <= length - DOLFIN_EPS
                result = indomain_x and indomain_y and indomain_z and on_boundary
                return result

  class poresurface_for_fusedpore(SubDomain):
 	 def inside(self,x,on_boundary):
  		indomain_x = near(x[0],0) or near(x[0], 15e-9)
  		indomain_z = x[2] >= DOLFIN_EPS and x[2] <= length - DOLFIN_EPS
  		result = indomain_x and indomain_z and on_boundary
 		return result
###
## --Define a shell with thickness d inside the pore, CF with in this shell has a smaller Dcf_bulk
## df_buffer = 1/(1 + Ks*B_tot/(Ks+cf0)**2))*Dcf_bulk 
## equation is from  Wagner, J., & Keizer, J. (1994). Effects of rapid buffers on Ca2+ diffusion and Ca2+ oscillations. Biophysical Journal, 67(1), 447.
## here suppose B_tot ( buffer concentration is one half of [CF]bulk = 0.5 mol/m^2 (mM)
## Ks is dissociation constant = 1
#
## 
  class buffer_shell(SubDomain):
       def inside(self, x, on_boundary):
    #  halfway point 
               indomain_x = x[2] >= DOLFIN_EPS and x[2] <= length - DOLFIN_EPS
               #indomain_x = x[2] >= 10e-9 and x[2] <= 50e-9
    # x,y coordinates for pore centers 
               pc1x = RevW/2
               pc1y = 0
               pc2x = RevW
               pc2y = RevH/2
               pc3x = RevW/2
               pc3y = RevH
               pc4x = 0
               pc4y = RevH/2

               r1sq = (x[0]-pc1x)**2 + (x[1]-pc1y)**2
               r1 = np.sqrt(r1sq) 
               p1_shell = r1 >= radius - thickness and r1 <= radius

               r2sq = (x[0]-pc2x)**2 + (x[1]-pc2y)**2
               r2 = np.sqrt(r2sq) 
               p2_shell = r2 >= radius - thickness and r2 <= radius

               r3sq = (x[0]-pc3x)**2 + (x[1]-pc3y)**2
               r3 = np.sqrt(r3sq) 
               p3_shell = r3 >= radius - thickness and r3 <= radius

               r4sq = (x[0]-pc4x)**2 + (x[1]-pc4y)**2
               r4 = np.sqrt(r4sq) 
               p4_shell = r4 >= radius - thickness and r4 <= radius
    
               inshell = p1_shell or p2_shell or p3_shell or p4_shell
               result = inshell and indomain_x
               return result  


  class bottom_boundary(SubDomain):
 	 def inside(self,x,on_boundary):
  		return near(x[2],-RevD) and on_boundary
    
  class top_boundary(SubDomain):
  	def inside(self,x,on_boundary):
  		return near(x[2],length+RevD) and on_boundary
  
  
  mesh = Mesh(meshfile)
  subdomain = FacetFunction("size_t",mesh,2)
  subdomain.set_all(0)
  
  bottomboundary = bottom_boundary()
  topboundary = top_boundary()
  
  bottomboundary.mark(subdomain,1) #mark the boundary
  topboundary.mark(subdomain,2)
  

# buffered region
  regions = CellFunction("size_t",mesh,2)
  regions.set_all(0)
  bufferRegion = buffer_shell()
  bufferRegion.mark(regions, 1)


  dx = Measure("dx",subdomain_data=regions,domain=mesh) 
  
  #efine MixedFunctionSpace--
  #--- for dolfin 1.6.0 (grimm)
  P = FunctionSpace(mesh,"P",1)
  V = MixedFunctionSpace([P,P,P,P,P,P])
  
  u = Function(V)
  u.interpolate(Constant((ck0,ccf0,ccl0,ch0,coh0,0.0)))
  
  #another way to assigin initial values:
  # 
  # newSpace = V.sub(0).collapse()
  # v0 = interpolate("expression",newSpace) --> interpolate is used for expression
  # ck0 = Function(newSpace)
  # ck0.vector()[:] = ck00*numpy.exp(-v0.vector()[:]/(kb*T))
  
  # ccl0 = Function(newSpace)
  # ccl0.vector()[:] = ccl00*numpy.exp(-v0.vector()[:]/(kb*T))
  
  # assign(u,[v0,ck0,ccl0])
  
  ck,ccf, ccl, ch, coh, v = split(u)
  ckm,ccfp, cclp, chp, cohp, vv = TestFunctions(V)
  
  #-----------------------------------
  # Write the variational forms
  #----------------------------------
  # flux of ions
  Jm = -Dk*(grad(ck) + charge*zk/(kb*T)*ck*grad(v))
  
  #---------------
  Ks = 1 # dissociation constant
  B_tot = Btot*ccf0 # concentration of buffering substrate, suppose = 1/2*[CF] = 0.5 mM
  Dcf_buffer = 1/(1 + Ks*B_tot/(Ks+ccf0)**2)*Dcf # diffusion constant in buffer region

  Jcf_bulk  = -Dcf*(grad(ccf) + charge*zcf/(kb*T)*ccf*grad(v))
  Jcf_buffer = -Dcf_buffer*(grad(ccf) + charge*zcf/(kb*T)*ccf*grad(v))
 

  Jp = -Dcl*(grad(ccl) + charge*zcl/(kb*T)*ccl*grad(v))
  Jh = -Dh*(grad(ch) + charge*zh/(kb*T)*ch*grad(v))
  Joh = -Doh*(grad(coh) + charge*zoh/(kb*T)*coh*grad(v))
  
  # get the LHS 
  #  -div(J) = 0 ==> J*grad(testf)*dx = 0
  # here, we consider dx(0) is unbuffered region and dx(1) is bufferd region
  aJm = inner(Jm,grad(ckm))*dx
  aJcf = inner(Jcf_bulk,grad(ccfp))*dx(0) + inner(Jcf_buffer,grad(ccfp))*dx(1) 
  aJp = inner(Jp,grad(cclp))*dx
  aJh = inner(Jh,grad(chp))*dx
  aJoh = inner(Joh,grad(cohp))*dx
  
  #--LHS and RHS of poisson equation
  aPoissonL = inner(grad(v),grad(vv))*dx
  aPoissonR = F/(eps0*epsr)*(ck*zk +ccf*zcf + ccl*zcl+ch*zh +coh*zoh)*vv*dx
  
  
  FF = aJm + aJp + aJh + aJoh + aJcf + aPoissonL - aPoissonR
  J = derivative(FF, u)
  
  
  #--------Boundary Conditions--------------------------
  #-- Ground Potential at the two ends of reservoirs
  bc1 = DirichletBC(V.sub(5),0,subdomain,1)
  bc2 = DirichletBC(V.sub(5),0,subdomain,2)
  
  #---------------------------------------
  # assigin boundary condition for K+ and Cl-
  bc3 = DirichletBC(V.sub(0),0,subdomain,1) #----> Assign a 0 [K+] at the botton reservor
  bc4 = DirichletBC(V.sub(0),ck0,subdomain,2)
  
  bc5x = DirichletBC(V.sub(1),0,subdomain,1)
  bc6x = DirichletBC(V.sub(1),ccf0,subdomain,2)
  bc5xx = DirichletBC(V.sub(2),0,subdomain,1)
  bc6xx = DirichletBC(V.sub(2),ccl0,subdomain,2)
  # assign boundary condition for H+ and OH-
  bc7 = DirichletBC(V.sub(3),0,subdomain,1)
  bc8 = DirichletBC(V.sub(3),ch0,subdomain,2)
  bc9 = DirichletBC(V.sub(4),0,subdomain,1)
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




  
  #-------------------------------------
  # TODO: Apply surface charge density as NeumanBC
  # Now works as DBC using Grahame euqaiton to convert sigma to psi
  #
  
  bcc = [bc1,bc2,bc3,bc4,bc7,bc8,bc9,bc10,bc5x,bc6x,bc5xx,bc6xx]
  
  #-------------------
  # Solve the problem
  #--------------------
  problem = NonlinearVariationalProblem(FF, u, bcs=bcc,J=J)
  solver = NonlinearVariationalSolver(problem)
  #solver.parameters["newton_solver"]["linear_solver"] = "gmres"
  #solver.parameters["newton_solver"]["preconditioner"] = "ilu"
# c_low = Constant(0.0)
# c_up = Constant(3000.0)
# lower = Function(V)
# upper = Function(V)
# ninfty = Function(P) ; ninfty.vector()[:] = -np.infty
# pinfty = Function(P) ; pinfty.vector()[:] = np.infty
# fa = FunctionAssigner(V,[P,P,P,P,P,P])
# fa.assign(lower,[interpolate(c_low,P),interpolate(c_low,P),interpolate(c_low,P),interpolate(c_low,P),interpolate(c_low,P),ninfty])
# fa.assign(upper,[interpolate(c_up,P),interpolate(c_up,P),interpolate(c_up,P),interpolate(c_up,P),interpolate(c_up,P),pinfty])
# problem.set_bounds(lower,upper)
# solver = NonlinearVariationalSolver(problem)
# snes_solver_parameters = {"nonlinear_solver": "snes",
#                         "snes_solver": {"linear_solver": "lu",
#                                         "maximum_iterations": 20,
#                                         "report": True,
#                                         "error_on_nonconvergence": False}}

# solver.parameters.update(snes_solver_parameters)
# info(solver.parameters, True)

  (iter, converged) = solver.solve()
  #solver.solve()

  ck_u,ccf_u,ccl_u,ch_u,coh_u,v_u = u.split(True) 
  
  TT = ck_u.function_space()
  degree = TT.ufl_element().degree()
  W = VectorFunctionSpace(mesh,'P',degree)
# fluxck = project(grad(ck_u)*Constant(-Dk),W)
# fluxccl = project(grad(ccl_u)*Constant(-Dcl),W)
# fluxch = project(grad(ch_u)*Constant(-Dh),W)
# fluxoh = project(grad(coh_u)*Constant(-Doh),W)
  fluxcf = project(grad(ccf_u)*Constant(-Dcf),W)
 
  #R = radius/nm
#  v1file = File("ccf_%s_%s.pvd"%(str(RRR),str(phi0)))
#  v1file << ccf_u

#  v3file = File("ck_%s_%s_%s.pvd"%(str(IonicS),str(ratio)))
#  v3file << ck_u
  v1file = File("cf_%s_%s_%s.pvd"%(str(rr),str(dd),str(Btot)))
  v1file << ccf_u
  v2file = File("ccfflux_%s_%s_%s.pvd"%(str(rr),str(dd),str(Btot)))
  v2file << fluxcf
# v3file = File("cclflux_%s_%s.pvd"%(str(cKCl),str(phi0)))
# v3file << fluxccl
# v4file = File("chflux_%s_%s.pvd"%(str(cKCl),str(phi0)))
# v4file << fluxch
# v5file = File("cohflux_%s_%s.pvd"%(str(cKCl),str(phi0)))
# v5file << fluxoh
# v6file = File("v_%s_%s.pvd"%(str(cKCl),str(phi0)))
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
      arg4=np.float(sys.argv[i+4])
      #arg7=(np.float(sys.argv[i+7]))
      runPNP(meshfile=arg1, rr=arg2,dd=arg3,Btot=arg4) #same as at the pickle line, justing doing radii
      quit()  
    if(arg=="-runConcs"):
      arg1=np.float(sys.argv[i+1])
      runPNP(cCaCl2=arg1) #value provided is in mM = mol/m^3
      quit() 




  raise RuntimeError("Arguments not understood")


