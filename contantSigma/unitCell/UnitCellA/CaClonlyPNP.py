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
from scipy.optimize import fsolve

import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import math
import cPickle as pickle
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
cKCl = 0

nm=1e-9

cCaCl2 = 100 #mol/m^3 == 1mM ----Bulk [KCl]

def runPNP(
  spacing = 5.2*nm, # the spacing between nanopores
  length = 34*nm, # length of nanopore SRB changed fromn 90 to 34 to match fig 4
  radius = 5.4*nm, # radius of nanopore
#meshfile = "/home/AD/bsu233/labscripts/poissonnernstplanck/contantSigma/unitCell/UnitCellA/UnitCellA.xml"
# PKH 
  meshfile = "/net/share/shared/papers/nanoporous/meshes/UnitCellA.xml",
  cCaCl2 = cCaCl2,
  pH = 7
  ):

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
  sigmaS = -20e-3 #C/m^2  --- Surface charge density 
  #ck0 = cKCl #initial K+ concentration                SRB -moved to inside function to manipulate cKCl input
  #ccl0 = cKCl #initial Cl- concentration
  zca = 2 # Ca2+
  zcl = -1 # Cl
  
  ch0 = 10**(-pH+3) #mol/m^3  initial [H]
  coh0 = 10**(-11+pH) #mol/m^3 inital [OH]
  #---------------------------
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
  u.interpolate(Constant((cca0,ccl0,ch0,coh0,0)))
  
  #another way to assigin initial values:
  # 
  # newSpace = V.sub(0).collapse()
  # v0 = interpolate("expression",newSpace) --> interpolate is used for expression
  # ck0 = Function(newSpace)
  # ck0.vector()[:] = ck00*numpy.exp(-v0.vector()[:]/(kb*T))
  
  # ccl0 = Function(newSpace)
  # ccl0.vector()[:] = ccl00*numpy.exp(-v0.vector()[:]/(kb*T))
  
  # assign(u,[v0,ck0,ccl0])
  
  cca, ccl, ch, coh, v = split(u)
  ccam, cclp, chp, cohp, vv = TestFunctions(V)
  
  ## TODO: apply a pH-regulated surface charge density
  ##
  #define pH regulated surface charge i?
  #chhh = Function(V.sub(2).collapse())
  #sigmas = Expression("-F*Gamma*((pow(10,-pKa+3) - pow(10,-pKb-3)*ch*ch)/(pow(10,-pKa+3) +ch+ pow(10,-pKb-3)*ch*ch))",F=F,Gamma=Gamma,ch=chhh,pKa=pKa,pKb=pKb,degree=1)
  

  #-----------------------------------
  # Write the variational forms
  #----------------------------------
  # flux of ions
  Jm = -Dca*(grad(cca) + charge*zca/(kb*T)*cca*grad(v))
  Jp = -Dcl*(grad(ccl) + charge*zcl/(kb*T)*ccl*grad(v))
  Jh = -Dh*(grad(ch) + charge*zh/(kb*T)*ch*grad(v))
  Joh = -Doh*(grad(coh) + charge*zoh/(kb*T)*coh*grad(v))
  
  # get the LHS 
  #  -div(J) = 0 ==> J*grad(testf)*dx = 0
  aJm = inner(Jm,grad(ccam))*dx
  aJp = inner(Jp,grad(cclp))*dx
  aJh = inner(Jh,grad(chp))*dx
  aJoh = inner(Joh,grad(cohp))*dx
  
  #--LHS and RHS of poisson equation
  aPoissonL = inner(eps0*epsr*grad(v),grad(vv))*dx
  aPoissonR = F*(cca*zca + ccl*zcl+ch*zh +coh*zoh)*vv*dx #+ sigmaS*vv*ds(1) #+ avg(inner(sigmaS,vv))*dS(1)
  
  
  FF = aJm + aJp + aJh + aJoh + aPoissonL - aPoissonR
  J = derivative(FF, u)
  
  
  #--------Boundary Conditions--------------------------
  #-- Ground Potential at the two ends of reservoirs
  bc1 = DirichletBC(V.sub(4),0,subdomain,1)
  bc2 = DirichletBC(V.sub(4),.02,subdomain,2)
  
  #---------------------------------------
  # assigin boundary condition for K+ and Cl-
  bc3 = DirichletBC(V.sub(0),0,subdomain,1) #----> Assign a 0 [K+] at the botton reservor
  bc4 = DirichletBC(V.sub(0),cca0,subdomain,2)
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

  phi0 = fsolve(Grahamequation,0)[0]
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
  solver.solve()
  
  cca_u,ccl_u,ch_u,coh_u,v_u = u.split(True) 
  
  v1file = File("cca.pvd")
  v1file << cca_u
  
  # Measure [Ca+] flux at the reservoir boundary
  TT = cca_u.function_space()
  degree = TT.ufl_element().degree()
  # define new functionspace to store flux
  W = VectorFunctionSpace(mesh,'P',degree)
  flux = project(grad(cca_u)*Constant(Dca),W)
  
  #print "Degree of cca", degree
  YY = v_u.function_space()
  newDegree = YY.ufl_element().degree()
 # print "Degree of v_u =", newDegree
  #VV = VectorFunctionSpace(mesh,'P',newDegree)
  #vProfile = project(grad(v_u),VV)
  #v_u.set_allow_extrapolation(True)
   
  R = radius/nm
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
  
  #voltsTop = u
  #mwe = v_u(1e-9,1e-9,length)
  #print " mwe = ", mwe
  n=FacetNormal(mesh)
  #SRB adding print statement to check why area = 0 when L = 34nm to match fig 4A data
  # area is commented out to quiet down the mismatch for ds*n=0 where length=/= xml data
  area = assemble(Constant(1.0)*ds(2,domain=mesh))
  #print area
  newVTop = assemble(v_u*ds(3,domain=mesh))
  newVBottom = assemble(v_u*ds(4,domain=mesh))
  flux_top = assemble(dot(flux,n)*ds(2))
 
  flux_topP = assemble(dot(flux,n)*ds(3))
 # print "Current top of pore:", flux_topP*F


  flux_botP = assemble(dot(flux,n)*ds(4))
  #print "Current bottom of pore:", flux_botP*F
  


  flux_midP = (flux_topP+flux_botP)/2
  midI = flux_midP*F
  #print "Current mid of pore:", midI

 
  #avgf = flux_top/area
  #Deff = avgf*(length + 2*RevH)/cca0/Dca
  tubeArea = assemble(Constant(1.0)*ds(3,domain=mesh)) 
  #vTop = assemble(dot(grad(v_u),n)*ds(3)) #assemble(dot(vProfile,n)*ds(3))
  #vBottom =assemble(dot((v_u),n)*ds(4)) # assemble(dot(vProfile,n)*ds(4))
  #print "vTop", newVTop
 # print "vBot", newVBottom 
 # print "Average Flux of K+ is",flux_top/area
  #print "Effective Diffusion constant is",Deff
  I = (flux_top)*F
  delVolts = (newVTop-newVBottom)/tubeArea
  avgf = flux_midP/(tubeArea*length + RevH*(spacing+radius)**2)
  Deff = avgf*(tubeArea)/cca0/Dca
  G = midI/delVolts
  #print "I is K+ * F = ", I #this should be mmol/s of K+ * F to get C/s for I
  #print "Voltage difference is", delVolts
  #print "G = ", midI/delVolts, "S"
  
  V2file = File("flux.pvd")
  V2file << flux
  V2file = File("v.pvd")
  V2file << v_u
  
  Results = { "Deff":"%s"%(str(Deff)),"G":"%s"%str(G)}
  pickle.dump(Results, open("Ca_%s_%s_%s.p"%(str(R),str(length/nm),str(cCaCl2)),"wb"))

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
      arg2=np.float(sys.argv[i+2])
      arg3=np.float(sys.argv[i+3])
      arg4=sys.argv[i+4]
      arg5=np.float(sys.argv[i+5])
      arg6=np.float(sys.argv[i+6])
      runPNP(spacing=arg1, length=arg2, radius=arg3, meshfile=arg4, cCaCl2=arg5, pH=arg6) 
      quit()  
    if(arg=="-runConcs"):
      arg1=np.float(sys.argv[i+1])
      runPNP(cCaCl2=arg1) #value provided is in mM = mol/m^3
      quit() 




  raise RuntimeError("Arguments not understood")


