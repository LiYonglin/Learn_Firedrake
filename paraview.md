### open Paraview-file on server
First, execute
```
mpiexec -n 5 pvserver
```
on server. Then
```
ssh -N -vvv -L 22222:127.0.0.1:11111 user_name@ip_address
```
on local machine. Finally, open Paraview and use `connect` to `22222`.
