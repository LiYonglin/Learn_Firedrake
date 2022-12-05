# For problem with Dirichlet boundary condition
#   -div(grad(u)) = f in Omega;  u = g on Bdry.
# Let u = D*un+ud, where un is the same as u except the value on boundary, 
# and ud denotes the vector on boundary and is zero inside Omega. Then for K*u=f
#   D'*K*(D*un+ud) = D'*f
# which yields
#   (D'*K*D)*un = D'*f-D'K*ud
# Noting D'*K*D is singular, we define A=D'*K*D+I_d, where I_d denotes the 
# matrix whose value is zero if and only if on boundary. Then
#   A*un = D'*K*D*un, A*ud = I_d*ud = ud
# Finally,
#   A*u = D'*f-D'*K*ud+ud

from firedrake.petsc import PETSc
rank, size = COMM_WORLD.rank, COMM_WORLD.size

def printf(*args, **kwargs):
    PETSc.Sys.Print(*args, **kwargs)

def syncprintf(*args, **kwargs):
    PETSc.Sys.syncPrint(*args, **kwargs)
    PETSc.Sys.syncFlush()

# create identity matrix
D = PETSc.Mat().createAIJ((None,10),(None,10))
D.setUp()
D.assemble()
D.shift(1)
syncprintf(D.view())
