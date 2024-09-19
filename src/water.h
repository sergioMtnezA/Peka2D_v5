#include "define.h"
#include "structs.h"

EXPORT_DLL void c_compute_cell_mass(int nTasks, t_arrays *arrays);

EXPORT_DLL void c_initialize_delta(int nTasks, t_arrays *arrays);

EXPORT_DLL void c_wall_rotated_calculus(int nTasks, t_arrays *arrays);

//EXPORT_DLL void c_bound_calculus

EXPORT_DLL void c_get_dtmin(int nTasks, t_arrays *arrays);

EXPORT_DLL void c_update_contributions(int nTasks, t_arrays *arrays);

EXPORT_DLL void c_checkpos_h(int nTasks, t_arrays *arrays, int *check);

EXPORT_DLL void c_reduce_dt(t_arrays *arrays);

EXPORT_DLL void c_set_new_dt(t_arrays *arrays);

//EXPORT_DLL void c_add_rain_contributions
//EXPORT_DLL void c_add_infiltration_contributions
//EXPORT_DLL void c_add_evaporation_contributions

EXPORT_DLL void c_update_cells(int nTasks, t_arrays *arrays);

EXPORT_DLL void c_check_wetdry(int nTasks, t_arrays *arrays);

EXPORT_DLL void c_update_wetdry_cells(int nTasks, t_arrays *arrays);

//EXPORT_DLL void c_update_boundaries

EXPORT_DLL void c_compute_mass_error(t_arrays *arrays);

EXPORT_DLL void c_reconstruct_active_elements(int nTasks, t_arrays *arrays);

EXPORT_DLL void c_add_active_cells(t_arrays *arrays, int id);

EXPORT_DLL void c_add_active_walls(t_arrays *arrays, int id);



