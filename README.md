## PeKa2D-v5.0 model for hydro-morphodynamical surface flows 

## Documentation
Code documentation for developers can be found [online](https://sergiomtneza.github.io/Peka2D_v5/) or within the local [docs](./docs/index.html)

## Numerical method
Find a detailed explanation on the numerical method implemented [here](./docs/water_Riemann_solver.pdf).

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


 

