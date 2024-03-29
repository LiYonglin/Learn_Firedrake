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


def add_DirichletBC(V, bc, a_mat, f_vec, ToSolve=False):
    """
    Add Dirichlet data (bc) to variational formulation a==f.
    Return: The solution is given by u=D*(inv(A)*b)+ud
        result.A: stiffness matrix
        result.b: RHS vector
        result.D: permutation matrix
        result.ud: Dirichlet value
        result.others
    """

    # Ku=f
    a_mat_size = a_mat.getSizes()
    rlgmap, clgmap = a_mat.getLGMap()
    a_bsize = a_mat.getBlockSize()

    # Find boundary nodes 
    last_index = bc.function_space().node_set.size
    _is_bc_nodes = bc.nodes < last_index
    bc_nodes = bc.nodes[_is_bc_nodes] # remove the ghost nodes
    len_bc_nodes = len(bc_nodes)

    # Create permutation matrix D: u=D*un+ud
    ncol_Local = a_mat_size[1][0]-len_bc_nodes*a_bsize
    nrow_Local = a_mat_size[0][0]
    D = PETSc.Mat().createAIJ((a_mat_size[0], (ncol_Local,None)), bsize=a_bsize)
    D.setUp()
    # Create the mapping from local to global for D
    _clgmap = PETSc.LGMap()
    owned_sz = np.array(ncol_Local//a_bsize, dtype=IntType)
    offset = np.empty_like(owned_sz)
    COMM_WORLD.Scan(owned_sz, offset)
    offset -= owned_sz
    indices = np.arange(offset, offset + owned_sz, dtype=IntType)
    _clgmap.create(indices, bsize=a_bsize, comm=COMM_WORLD)
    D.setLGMap(rlgmap, _clgmap)
    # Set D the identity matrix whose columns on boundary are removed
    not_bc_nodes = np.setdiff1d(np.arange(nrow_Local//a_bsize), bc_nodes)
    for i1 in range(ncol_Local//a_bsize):
        for i2 in range(a_bsize):
            D.setValueLocal(not_bc_nodes[i1]*a_bsize+i2, i1*a_bsize+i2, 1)
    D.assemble()

    # A = D'*K*D
    DTA = D.transposeMatMult(a_mat)
    A = DTA.matMult(D)
    A.setLGMap(_clgmap, _clgmap)

    # The values of solution on Dirichlet boundary
    b_fc = Function(V)
    bc.apply(b_fc)
    # Set the RHS vector b = D'*(f-K*ud) 
    b, _mid = D.createVecs() # _mid = D*b
    with b_fc.dat.vec_ro as ud:
        a_mat.mult(ud, _mid)
        _mid.aypx(-1, f_vec) # _mid=f-a*ud
        D.multTranspose(_mid, b) 

    # Solve the solution: u = D*(A\b)+ud
    if ToSolve:
        # KSP solver
        ksp = PETSc.KSP().create()
        ksp.setOperators(A)
        om = OptionsManager({'ksp_type': 'preonly', 'pc_type': 'lu'}, \
            options_prefix='solve_Dirichlet')
        om.set_from_options(ksp)
        uh = Function(V)
        un_vec = D.createVecRight() 
        with om.inserted_options():
            ksp.solve(b, un_vec)
        with uh.dat.vec as u_vec, b_fc.dat.vec_ro as ud:
            D.multAdd(un_vec, ud, u_vec)
    else:
        uh = None

    linear_system = {'A':A, 'b':b, 'D':D, 'ud':ud, 'uh':uh}
    return linear_system

# The following is for testing
def solve_Dirichlet_v1(V, bc, a, f):
    uh = Function(V)
    solve(a==f, uh, bcs=bc, \
        solver_parameters={"ksp_type": "preonly", "pc_type": "lu"}, \
            options_prefix='solve_v1')
    # printf('Ref. solution is')
    # with uh.dat.vec as vec:
    #     vec.view()

    # err = errornorm(exactu(x), uh, norm_type="H1")
    # unorm = norm(exactu(x), norm_type="H1", mesh=mesh)
    # # printf(f'H1 norm of u is {unorm}')
    # printf(np.real(err/unorm))

    return uh


def solve_Dirichlet_v2(V, bc, a, f):
    a_mat = assemble(a).petscmat
    # a_mat.view()
    f_ = assemble(f)

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
    # or use >>
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

    b, _mid = D.createVecs()
    with f_.dat.vec as f_vec:
        with b_fc.dat.vec as ud:
            a_mat.mult(ud, _mid)
            _mid.scale(-1)
            D.mult(_mid, b)
            b.axpy(1, ud)
        D.multAdd(f_vec, b, b) # b=D*f+b=D*f-D*K*ud+ud
    # b.view()

    # KSP solver
    ksp = PETSc.KSP().create()
    ksp.setOperators(A)
    om = OptionsManager({'ksp_type': 'preonly', 'pc_type': 'lu'}, \
        options_prefix='solve_v2')
    om.set_from_options(ksp)
    wh = Function(V)
    with om.inserted_options(), wh.dat.vec as w_vec:
        ksp.solve(b, w_vec)

    uh = solve_Dirichlet_v1(V, bc, a, f)
    syncprintf(np.allclose(wh.dat.data_ro_with_halos, uh.dat.data_ro_with_halos))

    err = errornorm(exactu(x), wh, norm_type="H1")
    unorm = norm(exactu(x), norm_type="H1", mesh=mesh)
    # printf(f'H1 norm of u is {unorm}')
    printf(np.real(err/unorm))


def solve_Dirichlet_v3(V, bc, a, f):
    a_mat = assemble(a).petscmat
    # a_mat.view()
    f_ = assemble(f)

    # Find boundary nodes
    _is_bc_nodes = bc.nodes < len(f_.dat.data)
    bc_nodes = bc.nodes[_is_bc_nodes]

    rlgmap, clgmap = a_mat.getLGMap()
    B = a_mat.copy()
    # B.setLGMap(rlgmap, clgmap)
    B.zeroRowsColumnsLocal(bc_nodes, diag=1)
    # B.view()

    # solution on Dirichlet boundary
    b_fc = Function(V)
    bc.apply(b_fc)

    b = B.createVecRight()
    g_nodes = rlgmap.apply(bc_nodes)
    with f_.dat.vec as f_vec:
        with b_fc.dat.vec as ud:
            a_mat.mult(ud, b)
            b.scale(-1)
            b.setValues(g_nodes, np.zeros(len(g_nodes)))
            f_vec.setValues(g_nodes, np.zeros(len(g_nodes)))
            # f_vec.setLGMap(rlgmap)
            # b.setValuesLocal(bc_nodes, np.zeros(len(bc_nodes)))
            # f_vec.setValuesLocal(bc_nodes, np.zeros(len(bc_nodes)))
            b.assemble()
            f_vec.assemble()
            b.axpy(1, ud) 
            b.axpy(1, f_vec) # b=D*f-D*K*ud+ud
    # b.view()

    # KSP solver
    ksp = PETSc.KSP().create()
    ksp.setOperators(B)
    om = OptionsManager({'ksp_type': 'preonly', 'pc_type': 'lu'}, \
        options_prefix='solve_v3')
    om.set_from_options(ksp)
    zh = Function(V)
    with om.inserted_options(), zh.dat.vec as z_vec:
        ksp.solve(b, z_vec)

    uh = solve_Dirichlet_v1(V, bc, a, f)
    syncprintf(np.allclose(zh.dat.data_ro_with_halos, uh.dat.data_ro_with_halos))

    err = errornorm(exactu(x), zh, norm_type="H1")
    unorm = norm(exactu(x), norm_type="H1", mesh=mesh)
    # printf(f'H1 norm of u is {unorm}')
    printf(np.real(err/unorm))

if __name__ == '__main__':
    # A simple test
    n = PETSc.Options().getInt('-n', default=10)
    mesh = UnitSquareMesh(nx=n, ny=n)
    x = SpatialCoordinate(mesh)
    V = FunctionSpace(mesh, 'CG', degree=1)
    bc = DirichletBC(V, rhsg(x), 'on_boundary')

    u = TrialFunction(V)
    v = TestFunction(V)

    a = inner(grad(u),grad(v))*dx + inner(u,v)*dx
    f = inner(rhsf(x),v)*dx

    solve_Dirichlet_v2(V, bc, a, f)
    solve_Dirichlet_v3(V, bc, a, f)
