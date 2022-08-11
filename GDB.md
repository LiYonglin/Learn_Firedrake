# find the pid
```
ps -aux | grep command1 | grep command2 
pgrep -fla command
```
# gdb
```
gdb -p PID
(gdb) py-list
(gdb) py-bt
```
