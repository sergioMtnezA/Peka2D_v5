## Mode PeKa2D-v5.0 for hydro-morphodynamical surface flows
Open docs/html/index.html in any web-browser (Firefox, Chrome, Edge, etc)

## Numerical method
Find a detailed explanation on the numerical method implemented [here](./water_Riemann_solver.pdf).

## CPU compilation
Set SET_SIMGPU=0 in define.h

`make`

`./peka pathFolder/ caseName`

**DEBUG mode**

`make DEBUG=yes`

`valgrind --tool=memcheck --leak-check=full --track-origins=yes --show-reachable=yes ./peka pathFolder/ caseName`

**Parallelzed OMP**

Set the number of CPU cores SET_OPENMP in define.h

`make OMP=yes`


## GPU compilation
Set SET_SIMGPU=1 in define.h

`make -f MakefileCUDA`

`./gpeka pathFolder/ caseName`

**DEBUG mode**

`make -f MakefileCUDA DEBUG=yes`

`compute-sanitizer --tool memcheck --leak-check full ./peka pathFolder/ caseName`


 

