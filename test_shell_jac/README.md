combine snes/tutorial/ex1f.F90 and mat/tutorial/ex6f.F90
mock original formJacobian with a shell matrix operation

2 methods possible:
(1) set Vec as ctx in mat shell. directly use solver vec x as jac shell mat when MatCreateShell. an empty FormShellJac is ready to go
(2) use a MatCtx derived type, and put base vec as its member. in FormShellJac one needs effort to update base vec
(3) use a MatCtx derived type, and put base vec as member. direct init this base vec as solver vec x. no need to update in FormShellJac
