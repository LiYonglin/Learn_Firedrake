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

from firedrake import *
from firedrake.petsc import PETSc, OptionsManager
import numpy as np

rank, size = COMM_WORLD.rank, COMM_WORLD.size

def printf(*args, **kwargs):
    PETSc.Sys.Print(*args, **kwargs)

def syncprintf(*args, **kwargs):
    PETSc.Sys.syncPrint(*args, **kwargs)
    PETSc.Sys.syncFlush()

def exactu(x):
    return (1-x**2)/4

def rhsf(x):
    u = exactu(x)
    return -div(grad(u))+u

def rhsg(x):
    u = exactu(x)
    return u

n = PETSc.Options().getInt('-n', default=10)
mesh = UnitSquareMesh(nx=n, ny=n)
x = SpatialCoordinate(mesh)
V = FunctionSpace(mesh, 'CG', degree=1)
bc = DirichletBC(V, rhsg(x), 'on_boundary')

u = TrialFunction(V)
v = TestFunction(V)

a = inner(grad(u),grad(v))*dx + inner(u,v)*dx
f = inner(rhsf(x),v)*dx

a_mat = assemble(a).petscmat
# a_mat.view()
f_ = assemble(f)

uh = Function(V)
solve(a==f, uh, bcs=bc)
printf('Ref. solution is')
with uh.dat.vec as vec:
    vec.view()


D = PETSc.Mat().createAIJ(a_mat.getSizes())
rlgmap, clgmap = a_mat.getLGMap()
D.setLGMap(rlgmap, clgmap)
D.setUp()
D.assemble()
D.shift(1)
# D.view()


# Find boundary nodes
_is_bc_nodes = bc.nodes < len(f_.dat.data)
bc_nodes = bc.nodes[_is_bc_nodes]
# syncprintf(f'[{rank}/{size}] nodes_local={bc_nodes}')
for i in range(len(bc_nodes)):
    D.setValueLocal(bc_nodes[i], bc_nodes[i], 0)
# D.zeroRowsColumnsLocal(bc_nodes, diag=0)
D.assemble()
# D.view()

A = D.matMatMult(a_mat, D)
A.setLGMap(rlgmap, clgmap)
A.zeroRowsColumnsLocal(bc_nodes, diag=1)
# A.view()

# solution on Dirichlet boundary
b_fc = Function(V)
bc.apply(b_fc)

# _mid = D'*K-I, where Ku=f
_mid = D.matMult(a_mat)
_mid.shift(-1) 
# _mid.view()

with f_.dat.vec as f_vec:
    with b_fc.dat.vec as ud:
        b = ud.copy()
        _mid.mult(ud, b)
        b.scale(-1)
    D.multAdd(f_vec, b, b) # b=D*f+b=D*f-D*K*ud+ud
# b.view()

# KSP solver
ksp = PETSc.KSP().create()
ksp.setOperators(A)
om = OptionsManager({'ksp_type': 'preonly', 'pc_type': 'lu'}, \
    options_prefix='test')
om.set_from_options(ksp)
wh = Function(V)
with om.inserted_options(), wh.dat.vec as w_vec:
    ksp.solve(b, w_vec)

syncprintf(np.allclose(wh.dat.data_ro_with_halos, uh.dat.data_ro_with_halos))

err = errornorm(exactu(x), wh, norm_type="H1")
unorm = norm(exactu(x), norm_type="H1", mesh=mesh)
printf(f'H1 norm of u is {unorm}')
printf(np.real(err/unorm))
