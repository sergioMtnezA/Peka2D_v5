#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <omp.h>

#include "define.h"
#include "structs.h"
#include "loadData.h"

////////////////////////////////////////////////////////////////
EXPORT_DLL int allocateArraysMemory(
    t_parameters spar, 
    t_mesh *mesh,
    t_arrays *carrays,
    t_message *msg);
/*----------------------------*/

////////////////////////////////////////////////////////////////
EXPORT_DLL int initilizeControlArrays(
    t_parameters spar, 
    t_mesh *mesh,
    t_arrays *carrays,
    t_message *msg);
/*----------------------------*/  

// ////////////////////////////////////////////////////////////////
// EXPORT_DLL int initilizeCommunicationArrays(
//     t_parameters spar, 
//     t_mesh *mesh,
//     t_arrays *carrays,   
//     t_cuPtr *cuPtr, 
//     t_message *msg);
// /*----------------------------*/

////////////////////////////////////////////////////////////////
EXPORT_DLL int initilizeMeshArrays(
    t_parameters spar, 
    t_mesh *mesh,
    t_arrays *carrays,
    t_message *msg);
/*----------------------------*/