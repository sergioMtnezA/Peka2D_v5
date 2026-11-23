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
/**
 * @brief Assigns mesh cells to rainfall polygons using their centers.
 * @param pksetup Structure containing rainfall polygon data.
 * @param mesh Computational mesh containing cell coordinates.
 * @param carrays Structure containing the rainRegion array.
 * @param msg Message handler for errors and notifications.
 * @return 1 on success, 0 on failure.
 */
EXPORT_DLL int assignRainPolygonsToCells(
    Peka2D_Setup *pksetup,
    t_mesh *mesh,
    t_arrays *carrays,
    t_message *msg);
/*----------------------------*/


////////////////////////////////////////////////////////////////
/**
 * @brief Releases CPU memory allocated for mesh computation arrays.
 * @param carrays Structure containing the allocated computation arrays.
 * @return 1 on success, 0 if the input pointer is NULL.
 */
EXPORT_DLL int freeMeshArraysMem(
    t_arrays *carrays);
/*----------------------------*/


////////////////////////////////////////////////////////////////
/**
 * @brief Allocates CPU memory for boundary computation arrays.
 * @param spar Simulation parameters.
 * @param mesh Computational mesh.
 * @param carrays Structure receiving the allocated arrays.
 * @param msg Message handler for errors and notifications.
 * @return 1 on success, 0 on failure.
 */
EXPORT_DLL int allocateBoundaryArraysMem(
    t_parameters spar, 
    t_mesh *mesh,
    t_arrays *carrays,
    t_message *msg);
/*----------------------------*/
//------------------------------------------//


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


#if SET_SOLUTE
////////////////////////////////////////////////////////////////
EXPORT_DLL int allocateSoluteArraysMem(
    t_parameters spar, 
    t_mesh *mesh,
    t_arrays *carrays,
    t_message *msg);
/*----------------------------*/

////////////////////////////////////////////////////////////////
EXPORT_DLL int initilizeSoluteArrays(
    t_parameters spar, 
    t_mesh *mesh,
    t_arrays *carrays,   
    t_message *msg);
/*----------------------------*/
#endif