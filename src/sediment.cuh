#include "define.h"
#include "structs.h"
#include "cuTilities.cuh"

#define FIXED_FRACTIONM 1
#define ACTIVE_FRACTIONM 2

#define EQUCONCF_BAGNOLD 1
#define EQUCONCF_WU 2
#define EQUCONCF_ZHANGANDXIE 3

#define WSF_RUBEY 1
#define WSF_ZHANG 2
#define WSF_ZANKE 3
#define WSF_VANRIJN 4 
#define WSF_RAUDKIVI 5
#define WSF_JULIEN 6
#define WSF_CHENG 7 
#define WSF_JIMENEZ 8
#define WSF_WU 9

#define ANISOTROPIC_DIFF 2
#define FRICTION_DIFF 3
#define EDDY_VISCOSITY 0
#define EDDY_VISCOSITY_LINEAR 0
#define EDDY_VISCOSITY_PARABOLIC 0
#define ws 0.00075


/**
 * @brief This function initializes the variation of conserved variables in the arrays structure to 0.0 and the array solidWallByCell to 0.
 * @param nTasks This integer variable passes the number of times the computation needs to be done.
 * @param arrays This pointer variable passes the arrays structure.
 */
__global__ void g_initialize_sediment_erosion_delta(int nTasks, t_arrays *arrays);


/**
 * @brief This function initializes the variation of conserved variables in the arrays structure to 0.0 and the array solidWallByCell to 0.
 * @param nTasks This integer variable passes the number of times the computation needs to be done.
 * @param arrays This pointer variable passes the arrays structure.
 */
__global__ void g_cell_sediment_Erosion_calculus(int nTasks, t_arrays *arrays);


/**
 * @brief This function initializes the variation of conserved variables in the arrays structure to 0.0 and the array solidWallByCell to 0.
 * @param nTasks This integer variable passes the number of times the computation needs to be done.
 * @param arrays This pointer variable passes the arrays structure.
 */
__global__ void g_update_sediment_erosion_cells(int nTasks, t_arrays *arrays);


/**
 * @brief This function initializes the variation of conserved variables in the arrays structure to 0.0 and the array solidWallByCell to 0.
 * @param nTasks This integer variable passes the number of times the computation needs to be done.
 * @param arrays This pointer variable passes the arrays structure.
 */
__global__ void g_multilayer_implicit_update_sediment_cells(int nTasks, t_arrays *arrays);

