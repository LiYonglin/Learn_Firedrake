### Print
```
from firedrake import *
from firedrake.petsc import PETSc, OptionsManager

rank, size = COMM_WORLD.rank, COMM_WORLD.size
PETSc.Sys.syncPrint(f'[{rank}/{size}]')
PETSc.Sys.syncFlush()
```
