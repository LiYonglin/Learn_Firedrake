from firedrake.petsc import PETSc

# create identity matrix
D = PETSc.Mat().createAIJ((None,10),(None,10))
D.setUp()
D.assemble()
D.shift(1)
syncprintf(D.view())
