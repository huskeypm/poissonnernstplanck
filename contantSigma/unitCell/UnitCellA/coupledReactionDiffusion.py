# 
# This code performs a reaction-diffusion simulation within a cube 
# General reaction-diffusion example for buffered Ca2+ 
# Largely borrowed from Cahn-Hilliard example 
#
# Validation
# -Buffering seems to be correct (set cInit to 45 everywhere, checked that cb final agreed with ExcessBuffer output)
# -anistripic diffusion is at least qualitiatively correct

from dolfin import *
import numpy as np

## My reaction system 
# dc/dt = del D del c - R(c,cB)
# dcB/dt = del D del cB + R(c,cB)
# R(c,CB) = kp * (B - cB)*c - km*cB

## Params 
field=True
if(field==True):
  Dii = Constant((1.0,1.0,1.0))
  Dii = Constant((5.,0.1,0.1))
  Dij = diag(Dii)
else:
  Dij = Constant(1.0)

kp = 1.0     
#km = 0.6 
cInit = 5.0 # [uM]
bT = 70  # [uM]  
Kd =10**(-6.5) *10**6 # km/kp; # [uM] 
dt     = 1.0e-03  
km = Kd*kp

xMax= 2
yMax = 2 * 0.59
zMax = 5 

# based on bers example for TnC on pg 46 in the Bers Bible
# NOTE: refer to Johan's email to use this more cleanly with code 
def ExcessBuffer(bMax=70,cFree=1.,Keq= 0.6,cBound=43.75):
  cb=sympy.Symbol('cb')
  c=sympy.Symbol('c')
  x=sympy.Symbol('x')
  Keqs = sympy.Symbol('Keqs')
  bt = sympy.Symbol('bt')
  # TESTeqn = sympy.Eq(c*bt / (Keqs+c),cb)
  # eqn.subs({bt:70,c:1.,Keqs:0.6})
  eqnAdded = sympy.Eq((c+x)*bt / (Keqs+(c+x)),(cb-x))

  sol  = sympy.solve(eqnAdded,x)
  #sol[0] - gives wacky numbers
  dx = sol[1].subs({bt:bMax,c:cFree,Keqs:Keq,cb:cBound})

  cFreeNew = cFree+dx
  cBoundNew = cBound-dx

  #print cFreeNew
  #print cBoundNew

  return(cFreeNew,cBoundNew)

#def ExcessBufferCheap():
#   # from soln1 of ExcessBuffer solution
#   dx = cb/2 - Keqs/2 - bt/2 - c/2 + (-2*bt*c - 2*bt*cb + 2*Keqs*bt + 2*Keqs*c + 2*Keqs*cb + 2*c*cb + Keqs**2 + bt**2 + c**2 + cb**2)**(1/2)/2

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
    if (
         (np.linalg.norm(x[0:1] -np.array([0,0]) ) < 0.15) or 
         (np.linalg.norm(x[0:1] -np.array([xMax,0]) ) < 0.15) or 
         (np.linalg.norm(x[0:1] -np.array([xMax,yMax]) ) < 0.15) or 
         (np.linalg.norm(x[0:1] -np.array([0,yMax]) ) < 0.15) 
      ):
      values[0] = cInit
      values[1] = CalccB( values[0] )
    else:
      values[0] = 0         
      values[1] = 0
  
   # print x
   # print values[0]
  def value_shape(self):
    return (2,)


# Class for interfacing with the Newton solver
class MyEqn(NonlinearProblem):
    def __init__(self, a, L):
        NonlinearProblem.__init__(self)
        self.L = L
        self.a = a
       # self.reset_sparsity = True
    def F(self, b, x):
        assemble(self.L, tensor=b)
    def J(self, A, x):
        assemble(self.a, tensor=A,)# reset_sparsity=self.reset_sparsity)
        #self.reset_sparsity = False



# Define mesh and function space 
#mesh = UnitSquare(16,16)
mesh = UnitCubeMesh(8,8,8)      
mesh.coordinates()[:]= np.array([xMax,yMax,zMax])* mesh.coordinates()
V = FunctionSpace(mesh,"CG",1)
ME = V *V 

# Trial and Test functions 
du = TrialFunction(ME) 
q,v  = TestFunctions(ME)

# Define function
u = Function(ME) # current soln
u0 = Function(ME) # prev soln

# split mixed functions
c,cb = split(u)
c0,cb0 = split(u0)

# Init conts
#init_cond = InitialConditions()
init_cond = InitialConditionsMyo()
u.interpolate(init_cond)
u0.interpolate(init_cond)


## Weak forms for RHS  
# See notetaker 121213 notes for details 

# Diffusion
#if(field==False):
#  RHS1 = -Dij*inner(grad(c),grad(q))*dx  
#  RHS2 = -Dij*inner(grad(cb),grad(v))*dx 
#else:
RHS1 = -inner(Dij*grad(c),grad(q))*dx  
#RHS2 = -inner(Dij*grad(cb),grad(v))*dx 

# Reaction: b + c --kp--> cb,  
#           b + c <-km--- cb
R = np.array([
  [-kp,km],     # s
  [kp,-km]      # p
  ])


# operator splitting 
opSplit=False
if(opSplit==False):
  RHS1 += (R[0,0]*(bT-cb)*c*q + R[0,1]*cb*q)*ds
  RHS2 = (R[1,0]*(bT-cb)*c*v + R[1,1]*cb*v)*ds


# Add in time dependence 
# (dc/dt) = RHS  --> c1-c0 - dt * RHS = 0
L1 = c*q*dx - c0*q*dx - dt * RHS1 
L2 = cb*v*dx - cb0*v*dx - dt * RHS2 
L = L1 + L2

# Compute directional derivative about u in the direction of du (Jacobian)
# (for Newton iterations) 
a = derivative(L, u, du)


# Create nonlinear problem and Newton solver
problem = MyEqn(a, L)
solver = NewtonSolver()
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

    if(opSplit==True):
      # Check Johans 121213/4 email regarding using dof maps and sympy functions here 
      1
      
    #file << (u.split()[0], t)
    file << (u,t)

    # check values
    for i,ele in enumerate(split(u)):
      tot = assemble(ele*dx(domain=mesh))
      vol = assemble(Constant(1.)*dx(domain=mesh))
      print "Conc(%d) %f " % (i,tot/vol)




#
