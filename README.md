# **CPU compilation**
[IMP] Set SET_SIMGPU=0 in define.h
make
[EXE] ./peka pathFolder/ caseName
*DEBUG mode*
make DEBUG=yes
[EXE] valgrind --tool=memcheck --leak-check=full --track-origins=yes --show-reachable=yes ./peka pathFolder/ caseName
*Parallelzed OMP*
make OMP=yes
The number of CPU cores is set with SET_OPENMP in define.h

# **GPU compilation**
[IMP] Set SET_SIMGPU=1 in define.h
make -f MakefileCUDA
[EXE] ./gpeka pathFolder/ caseName
*DEBUG mode*
make -f MakefileCUDA DEBUG=yes
[EXE] compute-sanitizer --tool memcheck --leak-check full ./peka pathFolder/ caseName
 

