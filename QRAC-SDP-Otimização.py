#%% Import
import qutip as qt # Quantum Mechanics Lib
import numpy as np # Standard math lib
import picos as pc # Optimization lib
#import statistics

def Pi_random(): # Return a random projector for dimension d
    # Definir theta e phi para o exemplo
    d=2
    theta = np.random.rand(1)[0] * np.pi # Theta da esfera de Bloch
    phi = np.random.rand(1)[0] * 2 * np.pi # Phi da esfera de Bloch
    # criação da unitária U
    out = np.cos(theta)*np.sin(theta)*np.exp(1j*phi)
   
    
    return np.array([[np.cos(theta)**2,out],[out.conjugate(),np.sin(theta)**2]])
    
n = 2 # Number of letters
for d in range(2,3):
    Id = np.eye(d) # Defining the Identity operator
    w =np.exp(2j*np.pi/d) # Defining the phase
    Pc = 0.5*(1 + 1/d) # Classical probability of success
    # Lists fopr the first phase (encoding)
    Pq = [] # first quantum/classical probability of success
    S = [] # 1st quantum success
    RHO = []
    RHO_0 = []
    # Lists for the second phase (decoding)
    Pq1 = [] # 2nd quantum/classical probability of success
    S1 = [] # 2nd quantum success
    M1 = []
    DeltaPq = [] 
    
    # Lists for first phase
    A = [] # Cumputational basis
    B = [] # Fourier basis
    rho0 = [[0 for _ in range(d)] for _ in range(d)]# pure sate to be optimized
    rho = []   # mixed sate to be optimized
    for x0 in range(d):
        for x1 in range(d):
            ketx=qt.Qobj(np.zeros(d)) # Constructiong one of the MUBs
            for j in range(d):
                ketx += d**(-1/2)*w**(x1*j)*qt.basis(d,j)
        
            b = ketx*ketx.dag()  # 2nd measurement operator
            B.append(b[:])
        a = qt.basis(d,x0)*qt.basis(d,x0).dag() # 1st measurement operator
        A.append(a[:])
    B=B[:d]
    
# OPTIMIZATION FOR THE STATE

F = pc.Problem() #Initiate first-phase solution

    
Success = 0 # Optimization variable
for x0 in range(d):
    for x1 in range(d):
        # State rho0 that will be optimizied 
        rho0[x0][x1]=pc.HermitianVariable(f"rho0_{x0}_{x1}", (d, d))
        # Restriction to the state be a state.
        F.add_constraint(rho0[x0][x1]>>0)
        F.add_constraint(pc.trace(rho0[x0][x1]) == 1)
        # Success definition
        Success += np.real((1/(2*d**2))*(pc.trace(A[x0]*rho0[x0][x1]) + pc.trace(B[x1]*rho0[x0][x1])))
                
# Our goal: maximize Success   
F.set_objective("max", Success) 
# Solve the problem: choose the solver (e.g: cvxopt or mosek).
F.solve(solver = "cvxopt")#, verbosity=1)

#storing the optimized success 
S=F.value
Pq=F.value/Pc

# storing the optimized states in a list
RHO = [[0 for _ in range(d)] for _ in range(d)]
for x0 in range(d):
    for x1 in range(d):
        RHO[x0][x1]=(np.array(rho0[x0][x1].value))
print("Success_RHO_opt = ",S) 
print("A = \n",np.round(A,3))
print("B = \n",np.round(B,3))
# print("RHO = \n",np.round(RHO,3))
        
#OPTIMIZATION FOR THE MEASUREMENT
RANDOM_RHO = 0 # 1 if random, 0 otherwise
d=2
if RANDOM_RHO ==1:
    RHO = [[0 for _ in range(d)] for _ in range(d)]
    for x0 in range(d):
        for x1 in range(d):
            RHO[x0][x1]=Pi_random()
                

n=[]
A1 = [0 for _ in range(d)]
B1 = [0 for _ in range(d)]
F1 = pc.Problem() # Initiate second-phase solution
"""
from now on we optimize the measurement using the optimized states
"""
# Defining the measurements.
for i in range(d):
    # Fisrt basis to decode x0
    A1[i]=pc.HermitianVariable(f"A1_{i}", (d, d)) 
    # Constraint to be a projector
    F1.add_constraint(A1[i]>>0)
    #F1.add_constraint(pc.trace(A1[i]) == 1)
    # Second basis to decode x1
    B1[i]=pc.HermitianVariable(f"B1_{i}", (d, d))
    # Constraint to be a projector
    F1.add_constraint(B1[i]>>0)
    #F1.add_constraint(pc.trace(B1[i]) == 1)
# Constraint to that the projectors of A1 (B1) form a complete basis.   
F1.add_constraint(pc.sum(A1)==Id)
F1.add_constraint(pc.sum(B1)==Id)

Success1 = 0 # Optimization variable
for x0 in range(d):
    for x1 in range(d):
        Success1 += np.real((1/(2*d**2))*(pc.trace(A1[x0]*RHO[x0][x1]) + pc.trace(B1[x1]*RHO[x0][x1])))


# Our goal: maximize Success   
F1.set_objective("max", Success1) 
# Solve the problem: choose the solver (e.g.: cvxopt or mosek).
F1.solve(solver = "cvxopt")#,verbosity=1)
S1=F1.value
print("Success_Measure_Opt = ",S1)
Am=[]
Bm=[]
for i in range(d):
    Am.append(np.array(A1[i].value))
    Bm.append(np.array(B1[i].value))
print("Am = \n",np.round(Am,3))
print("Bm = \n",np.round(Bm,3))
#print("RHO = \n",np.round(RHO,3))
