# PeKa2D-v5.0 model
PeKa2D-v5.0 is a numerical model for the simulation of hydro-morphodynamical surface flows.

---

## Documentation
Code documentation for developers can be found:
* [online](https://sergiomtneza.github.io/Peka2D_v5/)
* [repo](./docs/index.html)

## Numerical method
Find a detailed explanation on the implemented numerical method:
* [repo](./docs/water_Riemann_solver.pdf).

---

## GPU compilation
To enable GPU computation, set `SET_SIMGPU=1` in define.h

Compile with:
`make -f MakefileCUDA`

Run with:
`./gpeka pathFolder/ caseName`

**Options:**
- `RECONSTRUC_ACTIVE=1` 
  Enable active cell/wall array completion within the time loop.
- `UPDATE_ACTIVE_ARRAYS=1` 
  Enable full reconstruction of active cell/wall arrays every `nIterArrangeActElem` iterations.  
 
## Solute transport computation
To enable solute transport computation, set `SET_SOLUTE=1` in define.h

**Options:**
- `SET_SOLUTE_UNROLL=1` 
  Enable unroll multi-solute computation in CUDA kernels

## DEBUG mode
Compile in debug mode with:
`make -f MakefileCUDA DEBUG=yes`

Run using NVIDIA Compute Sanitizer:
`compute-sanitizer --tool memcheck --leak-check full ./peka pathFolder/ caseName`

