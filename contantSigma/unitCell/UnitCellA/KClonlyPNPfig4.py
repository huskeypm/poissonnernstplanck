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
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import math

#parameters 
#----------------Domain Size--------------

nm=1e-9
spacing = 5.2*nm # the spacing between nanopores
########### SRB-Changed length from 90 to 34 
length = 34*nm # length of nanopore   #SRB - changed for Fig 4
radius = 5.1*nm # radius of nanopore

RevL= spacing + 2*radius #length of reservoir
RevH=20*nm # depth of reservoir
#--------------------------------------------------
#sigmaS = -20e-3 #C/m^2  --- Surface charge density 

sigmaS = -20e-3 #C/m^2  --- Surface charge density 

cKCl = 100 #mol/m^3 == 1mM ----Bulk [KCl]
#ck0 = cKCl #initial K+ concentration       SRB -moved to inside function to manipulate cKCl input
#ccl0 = cKCl #initial Cl- concentration
zk = 1 # K+ 
zcl = -1 # Cl

pH = 7.5
ch0 = 10**(-pH+3) #mol/m^3  initial [H]  # reckless changes incoming
coh0 = 10**(-11+pH) #mol/m^3 inital [OH]


#ch0 = 10**(-pH) #mol/m^3  initial [H]  # reckless changes incoming
#coh0 = 10**(-14+pH) #mol/m^3 inital [OH]

Dk = 1.96e-9 # m^2/s --Diffusion constant for K+
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
#class top_pore(SubDomain):
     
class bottom_boundary(SubDomain):
    def inside(self,x,on_boundary):
        return near(x[2],-RevH) and on_boundary

class top_boundary(SubDomain):
	def inside(self,x,on_boundary):
		return near(x[2],length+RevH) and on_boundary



def runPNP(
  spacing = 5.2*nm, # the spacing between nanopores
  length = 34*nm, # length of nanopore SRB changed fromn 90 to 34 to match fig 4
  radius = 5.1*nm, # radius of nanopore
#meshfile = "/home/AD/bsu233/labscripts/poissonnernstplanck/contantSigma/unitCell/UnitCellA/UnitCellA.xml"
# PKH 
  meshfile = "/net/share/shared/papers/nanoporous/meshes/UAL34R5-1.xml",
  cKCl = cKCl  
  ):

  

  
  #Adding ability to change KCl concentration applied to the top for Figure 4
  ck0 = cKCl #initial K+ concentration
  ccl0 = cKCl #initial Cl- concentration
   
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
  V = MixedFunctionSpace([P1,P2,P3,P4,P5])
  
  #--- for dolfin 1.6.2 (kafka)
  #P1 = FiniteElement('P',tetrahedron,1)
  #P = MixedElement([P1,P1,P1,P1,P1])
  #V = FunctionSpace(mesh,P)
  
  # assign initial values for v,ck,ccl
  # two ways to do that:
  u = Function(V)
  u.interpolate(Constant((ck0,ccl0,ch0,coh0,0)))
  
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
  #chhh = Function(V.sub(2).collapse())
  #sigmas = Expression("-F*Gamma*((pow(10,-pKa+3) - pow(10,-pKb-3)*ch*ch)/(pow(10,-pKa+3) +ch+ pow(10,-pKb-3)*ch*ch))",F=F,Gamma=Gamma,ch=chhh,pKa=pKa,pKb=pKb,degree=1)
  

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
  aPoissonL = inner(eps0*epsr*grad(v),grad(vv))*dx
  aPoissonR = F*(ck*zk + ccl*zcl+ch*zh +coh*zoh)*vv*dx #+ sigmaS*vv*ds(1) #+ avg(inner(sigmaS,vv))*dS(1)
  
  
  FF = aJm + aJp + aJh + aJoh + aPoissonL - aPoissonR
  J = derivative(FF, u)
  
  
  #--------Boundary Conditions--------------------------
  #-- Ground Potential at the two ends of reservoirs
  bc1 = DirichletBC(V.sub(4),0,subdomain,1)
  bc2 = DirichletBC(V.sub(4),0.02,subdomain,2)
  
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
  phi0 = math.asinh((sigmaS/(0.117*(ck0/1000)**0.5)))*51.4/1000
  bcx = DirichletBC(V.sub(4),Constant(phi0),subdomain,3) # electric potential on the surface of sphere
  
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
  solver.solve()
  
  ck_u,ccl_u,ch_u,coh_u,v_u = u.split(True) 
  
  v1file = File("ck.pvd")
  v1file << ck_u
  voltsfile = File("v.pvd")
  voltsfile << v_u
  # Measure [K+] flux at the reservoir boundary
  TT = ck_u.function_space()
  degree = TT.ufl_element().degree()
  # define new functionspace to store flux
  W = VectorFunctionSpace(mesh,'P',degree)
  flux = project(grad(ck_u)*Constant(Dk),W)
  
  n=FacetNormal(mesh)
  #SRB adding print statement to check why area = 0 when L = 34nm to match fig 4A data
  # area is commented out to quiet down the mismatch for ds*n=0 where length=/= xml data
  area = assemble(Constant(1.0)*ds(2,domain=mesh))
  print area
  #vtopPore = 
  flux_top = assemble(dot(flux,n)*ds(2))
  avgf = flux_top/area
  Deff = avgf*(length + 2*RevH)/ck0/Dk
  
  print "Average Flux of K+ is", flux_top/area
  print "Effective Diffusion constant is", Deff
  I = (flux_top)*F
  
  print "In volts phi0 is ", phi0 
  print "I is K+ * F = ", I #this should be mmol/s of K+ * F to get C/s for I
  V2file = File("flux.pvd")
  V2file << flux

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
      arg1=np.float(sys.argv[i+1]) 
      runPNP(spacing=arg1) 
      quit()  
    if(arg=="-runConcs"):
      arg1=np.float(sys.argv[i+1])
      runPNP(cKCl=arg1) #value provided is in mM = mol/m^3
      quit() 




  raise RuntimeError("Arguments not understood")




