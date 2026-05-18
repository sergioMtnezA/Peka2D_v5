## Square mesh creation

1. Compile and execute mesh generator
g++ structuredMesh.c -o mesher 
./mesher

2. Compile and execute FED file generator
g++ convertPEKA2FED.c -o convert
./convert case1

3. Check the mesh topology in file case1.vtk

3. Copy FED file to the simulation folder
