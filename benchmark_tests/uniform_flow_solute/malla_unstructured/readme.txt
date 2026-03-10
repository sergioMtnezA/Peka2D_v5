Instruciones actuales
######################
//execute triangle
triangle -pne -q32 -a2.0 mesh_triangle.poly

//for refining mesh
#./refine mesh_triangle.1 polygons.bln NO
#triangle -pnera -q32 mesh_triangle.1.poly NO

//convert Triangle to Peka

g++ convertTRIANGLE2PEKA.c -o convertT2P
./convertT2P


## Unstructured mesh creation

1. Execute TRIANGLE mesh generator
triangle -pne -q32 -a2.0 mesh_triangle.poly

2. Compile and execute PEKA file generator
g++ convertTriangle2Peka.c -o run2Peka
./run2Peka

3. Compile and execute FED file generator
g++ convertPeka2FED.c -o run2FED
./run2FED case1

4. Check the mesh topology in file case1.vtk

5. Copy FED file to the simulation folder




