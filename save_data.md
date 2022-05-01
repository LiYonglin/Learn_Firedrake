### Print mesh to Matlab file
```
from scipy import io
V = uh.ufl_function_space()
p = mesh.coordinates.dat.data.real
t = V.cell_node_list
io.savemat('mesh_pt.mat',{'p': p, 't': t})
```

### Print vector and matrix to dat file
```
viewer = PETSc.Viewer().createBinary('data.dat', mode='w')
A.view(viewer) # A is a matrix
b.view(viewer) # b is a vector
viewer.flush()

"""
Get data in MATLAB:
    fd = PetscOpenFile('data.dat');
    A = PetscBinaryRead(fd, 'complex', 1, 'indices', 'int64');
    b = PetscBinaryRead(fd, 'complex', 1, 'indices', 'int64');
"""
```
