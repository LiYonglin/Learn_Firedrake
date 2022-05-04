### Create a ksp solver
```
# solve Au=f
A = assemble(a)
A_ = A.petscmat # get the PETSc matrix
ksp = PETSc.KSP().create() 
ksp.setOperators(A)
om = OptionsManager({'ksp_type': 'preonly', 'pc_type': 'lu'}, \
    options_prefix='grad_g')
om.set_from_options(ksp)  
with om.inserted_options():
    ksp.setFromOptions()
# solve Ax=b
u = Function(V)
with u.dat.vec as x:
    with f.dat.vec as b:
        ksp.solve(b, x)
```
