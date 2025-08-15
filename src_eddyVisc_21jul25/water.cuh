#include "define.h"
#include "structs.h"
#include "cuTilities.cuh"

/**
 * @brief This function recounts and stores the total mass of all active cells.
 * @param nTasks This integer variable passes the number of times the computation needs to be done.
 * @param arrays This pointer variable passes the arrays structure.
 * @param mass This pointer variable passes the mass variable to store it.
 */
__global__ void g_compute_cell_mass(int nTasks, t_arrays *arrays, double *mass);

/**
 * @brief This function initializes the variation of conserved variables in the arrays structure to 0.0 and the array solidWallByCell to 0.
 * @param nTasks This integer variable passes the number of times the computation needs to be done.
 * @param arrays This pointer variable passes the arrays structure.
 */
__global__ void g_initialize_delta(int nTasks, t_arrays *arrays);

/**
 * @brief This function initializes the variation of conserved variables in the arrays structure to 0.0 and the array solidWallByCell to 0.
 * @param nTasks This integer variable passes the number of times the computation needs to be done.
 * @param localDt This pointer variable gets the values of the local time steps.
 */
__global__ void g_wall_rotated_calculus(int nTasks, t_arrays *arrays, double *localDt);

//__global__ void g_bound_calculus

/**
 * @brief This function finds the minimum time step out of the vector localDt. If the whole domain is dry, it is set to 1.0 seconds.
 * @param arrays This pointer variable passes the arrays structure. arrays->dt is modifed.
 * @param localDt This pointer variable contains the local time steps.
 * @param idmin This pointer variable indicates the component of the minimum time step.
 */
__global__ void g_get_dtmin(t_arrays *arrays, double *localDt, int *idmin);

/**
 * @brief This function computes the right-hand side of the equations to update each physical variable.
 * @param nTasks This integer variable passes the number of times the computation needs to be done.
 * @param arrays This pointer variable passes the arrays structure. arrays->dt is modifed.
 */
__global__ void g_update_contributions(int nTasks, t_arrays *arrays);

/**
 * @brief This function checks the positivity of the water depth. If there is a cell where the depth is negative (-tol), it write a notification flag in variable check.
 * @param nTasks This integer variable passes the number of times the computation needs to be done.
 * @param arrays This pointer variable passes the arrays structure.
 * @param check This integer acts as a flag to check the positivity of the water depth and is used in calculate.cpp to reduced the time step.
 */
__global__ void g_checkpos_h(int nTasks, t_arrays *arrays, int *check);

/**
 * @brief This function multiplies the time step by a factor of 0.8 if the there is a change to get negative water depths (as checked in function c_checkpos_h).
 * @param arrays This pointer variable passes the arrays structure.
 */
__global__ void g_reduce_dt(t_arrays *arrays);

/**
 * @brief This function sets the new time step according to the dumping and final times. If the precomputed time step is too big, it is reduced to match the desired time.
 * @param arrays This pointer variable passes the arrays structure.
 */
__global__ void g_set_new_dt(t_arrays *arrays);

//__global__ void g_add_rain_contributions
//__global__ void g_add_infiltration_contributions
//__global__ void g_add_evaporation_contributions

/**
 * @brief This function updates the water depth and water discharges (in wet cells).
 * @param nTasks This integer variable passes the number of times the computation needs to be done.
 * @param arrays This pointer variable passes the arrays structure.
 */
__global__ void g_update_cells(int nTasks, t_arrays *arrays);

/**
 * @brief This function checks the presence of wet/dry fronts and writes flags to correct them later.
 * @param nTasks This integer variable passes the number of times the computation needs to be done.
 * @param arrays This pointer variable passes the arrays structure.
 */
__global__ void g_check_wetdry(int nTasks, t_arrays *arrays);

/**
 * @brief This function corrects the wet/dry fronts found in function c_check_wetdry.
 * @param nTasks This integer variable passes the number of times the computation needs to be done.
 * @param arrays This pointer variable passes the arrays structure.
 */
__global__ void g_update_wetdry_cells(int nTasks, t_arrays *arrays);

//__global__ void c_update_boundaries

/**
 * @brief This function computes the mass error of the current state with respect to the previous one (in %).
 * @param arrays This pointer variable passes the arrays structure.
 */
__global__ void g_compute_mass_error(t_arrays *arrays);

/**
 * @brief This function computes the mass error of the current state with respect to the previous one (in %).
 * @param arrays This pointer variable passes the arrays structure.
 */
__global__ void g_reconstruct_active_elements(int nTasks, t_arrays *arrays);

/**
 * @brief This function adds the newly wetted cells to the active cells set.
 * @param arrays This pointer variable passes the arrays structure.
 * @param id This integer variable passes the index of the wetted cell.
 */
__device__ void g_add_active_cells(t_arrays *arrays, int id);

/**
 * @brief This function adds wall identified by index id to the active walls set.
 * @param arrays This pointer variable passes the arrays structure.
 * @param id This integer variable passes the index of the wetted cell.
 */
__device__ void g_add_active_walls(t_arrays *arrays, int id);





