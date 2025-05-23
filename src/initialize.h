#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <omp.h>

#include "define.h"
#include "structs.h"
#include "loadData.h"


////////////////////////////////////////////////////////////////
EXPORT_DLL int initilizeComputationControls(
    t_parameters spar, 
    t_mesh *mesh,
    t_arrays *carrays,
    t_message *msg);
/*----------------------------*/  

////////////////////////////////////////////////////////////////
EXPORT_DLL int allocateMeshArraysMem(
    t_parameters spar, 
    t_mesh *mesh,
    t_arrays *carrays,
    t_message *msg);
/*----------------------------*/

////////////////////////////////////////////////////////////////
EXPORT_DLL int initilizeMeshArrays(
    t_parameters spar, 
    t_mesh *mesh,
    t_arrays *carrays,
    t_message *msg);
/*----------------------------*/

////////////////////////////////////////////////////////////////
EXPORT_DLL int allocateBoundaryArraysMem(
    t_parameters spar, 
    t_mesh *mesh,
    t_arrays *carrays,
    t_message *msg);
/*----------------------------*/

////////////////////////////////////////////////////////////////
EXPORT_DLL int initilizeBoundaryControlArrays(
    t_parameters spar, 
    t_mesh *mesh,
    t_arrays *carrays,
    t_message *msg);
/*----------------------------*/

////////////////////////////////////////////////////////////////
EXPORT_DLL int initilizeBoundaryMeshArrays(
    t_parameters spar, 
    t_mesh *mesh,
    t_arrays *carrays,
    t_message *msg);
/*----------------------------*/
