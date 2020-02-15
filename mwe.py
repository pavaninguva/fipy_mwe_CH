from fipy import CellVariable, Grid2D, GaussianNoiseVariable, DiffusionTerm, TransientTerm, ImplicitSourceTerm, VTKCellViewer, LinearLUSolver, parallelComm
from fipy.tools import numerix
import time

TIME_STRIDE = 500
TIME_MAX = 2000
DT = 2.0
chi_AB = 0.0075
N_A = 400
N_B = 400
NOISE_MAGNITUDE = 0.01
A_RAW = 0.5

print ("Yay")

# # Define mesh
mesh = Grid2D(nx=64.0, ny=64.0, dx=1.0, dy=1.0)
print ("mesh loaded")

# We need to define the relevant variables: 
x_a = CellVariable(name=r"x_a", mesh = mesh, hasOld=1)
mu_AB = CellVariable(name=r"mu_AB", mesh = mesh, hasOld=1)

x_a.setValue(GaussianNoiseVariable(mesh=mesh,
                                   mean=A_RAW,
                                   variance=NOISE_MAGNITUDE)
)


dgdx_a = ((1.0/N_A) - (1.0/N_B)) + (1.0/N_A)*numerix.log(x_a) - (1.0/N_B)*numerix.log(1.0 - x_a) + chi_AB*(1.0 - 2*x_a)
d2gdx_a2 = (1.0/(N_A*x_a)) + (1.0/(N_B*(1.0 - x_a))) - 2*chi_AB

# Define the equations

# evaluating kappa
kappa = (1.0/6.0)*chi_AB

# eqn 1 is the 2nd order transport equation
eq1 = (TransientTerm(var=x_a)) == DiffusionTerm(coeff = x_a * (1 - x_a), var=mu_AB)

# eqn 2 is the chemical potential definition
eq2 = (ImplicitSourceTerm(coeff=1. , var=mu_AB)) == ImplicitSourceTerm(coeff=d2gdx_a2, var=x_a) - d2gdx_a2 * x_a + dgdx_a - DiffusionTerm(coeff=kappa, var=x_a)

# Adding the equations together
eq = eq1 & eq2

elapsed = 0.
dt = DT
if __name__ == "__main__":
    duration = TIME_MAX


time_stride = TIME_STRIDE
timestep = 0

# Defining the solver 
solver = LinearLUSolver(tolerance=1e-10, iterations=25)
# solver = PETSc.KSP().create()
start = time.time()

while elapsed < duration: 
    if (timestep == 0):
        vw = VTKCellViewer(vars=(x_a, mu_AB))
        vw.plot(filename="0_output.%d.vtk" %(parallelComm.procID))
    elapsed += dt
    timestep += 1
    x_a.updateOld()
    mu_AB.updateOld()
    res = 1e+10
    while res > 1e-10:
        res = eq.sweep(dt=dt, solver=solver)
        print ("sweep!")
    print (elapsed)
    end = time.time()
    print(end-start)
    if (timestep % time_stride ==0):
        vw = VTKCellViewer(vars=(x_a, mu_AB))
        vw.plot(filename="%s_output.%d.vtk" %(elapsed, parallelComm.procID))