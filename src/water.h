#include "define.h"
#include "structs.h"

/**
 * @brief This function recounts and stores the total mass of all active cells.
 * @param nTasks This integer variable passes the number of times the computation needs to be done.
 * @param arrays This pointer variable passes the arrays structure.
 */
EXPORT_DLL void c_compute_cell_mass(int nTasks, t_arrays *arrays);

/**
 * @brief This function initializes the variation of conserved variables in the arrays structure to 0.0 and the array solidWallByCell to 0.
 * @param nTasks This integer variable passes the number of times the computation needs to be done.
 * @param arrays This pointer variable passes the arrays structure.
 */
EXPORT_DLL void c_initialize_delta(int nTasks, t_arrays *arrays);

/**
 * @brief This function computes the variation of conserved variables (dh, dhu, dhv) using the ARoe method at active walls. It includes source terms (bed slope and bed friction), positivity, entropy and wet/dry fixes. Additionally, it updates aux arrays and the list of active cells.
 * @param nTasks This integer variable passes the number of times the computation needs to be done.
 * @param arrays This pointer variable passes the arrays structure.
 */
EXPORT_DLL void c_wall_rotated_calculus(int nTasks, t_arrays *arrays);

//EXPORT_DLL void c_bound_calculus

/**
 * @brief This function finds the minimum time step out of the vector localDt. If the whole domain is dry, it is set to 1.0 seconds.
 * @param nTasks This integer variable passes the number of times the computation needs to be done.
 * @param arrays This pointer variable passes the arrays structure. arrays->dt is modifed.
 */
EXPORT_DLL void c_get_dtmin(int nTasks, t_arrays *arrays);

/**
 * @brief This function finds the minimum time step out of the vector localDt. If the whole domain is dry, it is set to 1.0 seconds.
 * @param nTasks This integer variable passes the number of times the computation needs to be done.
 * @param arrays This pointer variable passes the arrays structure. arrays->dt is modifed.
 */
EXPORT_DLL void c_update_contributions(int nTasks, t_arrays *arrays);

/**
 * @brief This function checks the positivity of the water depth. If there is a cell where the depth is negative (-tol), it write a notification flag in variable check.
 * @param nTasks This integer variable passes the number of times the computation needs to be done.
 * @param arrays This pointer variable passes the arrays structure.
 * @param check This integer acts as a flag to check the positivity of the water depth and is used in calculate.cpp to reduced the time step.
 */
EXPORT_DLL void c_checkpos_h(int nTasks, t_arrays *arrays, int *check);

/**
 * @brief This function multiplies the time step by a factor of 0.8 if the there is a change to get negative water depths (as checked in function c_checkpos_h).
 * @param arrays This pointer variable passes the arrays structure.
 */
EXPORT_DLL void c_reduce_dt(t_arrays *arrays);

/**
 * @brief This function sets the new time step according to the dumping and final times. If the precomputed time step is too big, it is reduced to match the desired time.
 * @param arrays This pointer variable passes the arrays structure.
 */
EXPORT_DLL void c_set_new_dt(t_arrays *arrays);

//EXPORT_DLL void c_add_rain_contributions
//EXPORT_DLL void c_add_infiltration_contributions
//EXPORT_DLL void c_add_evaporation_contributions

/**
 * @brief This function updates the water depth and water discharges (in wet cells).
 * @param nTasks This integer variable passes the number of times the computation needs to be done.
 * @param arrays This pointer variable passes the arrays structure.
 */
EXPORT_DLL void c_update_cells(int nTasks, t_arrays *arrays);

/**
 * @brief This function checks the presence of wet/dry fronts and writes flags to correct them later.
 * @param nTasks This integer variable passes the number of times the computation needs to be done.
 * @param arrays This pointer variable passes the arrays structure.
 */
EXPORT_DLL void c_check_wetdry(int nTasks, t_arrays *arrays);

/**
 * @brief This function corrects the wet/dry fronts found in function c_check_wetdry.
 * @param nTasks This integer variable passes the number of times the computation needs to be done.
 * @param arrays This pointer variable passes the arrays structure.
 */
EXPORT_DLL void c_update_wetdry_cells(int nTasks, t_arrays *arrays);

//EXPORT_DLL void c_update_boundaries

/**
 * @brief This function corrects the wet/dry fronts found in function c_check_wetdry.
 * @param nTasks This integer variable passes the number of times the computation needs to be done.
 * @param arrays This pointer variable passes the arrays structure.
 */
EXPORT_DLL void c_compute_mass_error(t_arrays *arrays);

EXPORT_DLL void c_reconstruct_active_elements(int nTasks, t_arrays *arrays);

EXPORT_DLL void c_add_active_cells(t_arrays *arrays, int id);

EXPORT_DLL void c_add_active_walls(t_arrays *arrays, int id);



