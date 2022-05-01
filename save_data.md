### Print mesh to Matlab file
```
V = uh.ufl_function_space()
p = mesh.coordinates.dat.data.real
t = V.cell_node_list
io.savemat('data2mat/mesh_pt.mat',{'p': p, 't': t})

