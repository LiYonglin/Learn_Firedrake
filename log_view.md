### Create log files or informations for debug
To add your own events
```
from firedrake.petsc import PETSc

with PETSc.Log.Event("foo"):
    do_something_expensive()
```
or
```
from firedrake.petsc import PETSc

@PETSc.Log.EventDecorator("foo")
def do_something_expensive():
    ...
```

To generate the flame graph, you should execute programs on Terminal by
```
python file.py -ksp_view -dm_view -log_view :log_file.txt:ascii_flamegraph
```

The log_file.txt can be graphed on https://www.speedscope.app/

See https://www.firedrakeproject.org/optimising.html for more details.
