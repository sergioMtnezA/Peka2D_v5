#include "define.h"
#include "structs.h"
#include "cuTilities.cuh"

__global__ void g_compute_cell_mass(int nTasks, t_arrays *arrays, double *mass);

__global__ void g_initialize_delta(int nTasks, t_arrays *arrays);

__global__ void g_wall_rotated_calculus(int nTasks, t_arrays *arrays, double *localDt);

//__global__ void g_bound_calculus

__global__ void g_get_dtmin(t_arrays *arrays, double *localDt, int *idmin);

__global__ void g_update_contributions(int nTasks, t_arrays *arrays);

__global__ void g_checkpos_h(int nTasks, t_arrays *arrays, int *check);

__global__ void g_reduce_dt(t_arrays *arrays);

__global__ void g_set_new_dt(t_arrays *arrays);

//__global__ void g_add_rain_contributions
//__global__ void g_add_infiltration_contributions
//__global__ void g_add_evaporation_contributions

__global__ void g_update_cells(int nTasks, t_arrays *arrays);

__global__ void g_check_wetdry(int nTasks, t_arrays *arrays);

__global__ void g_update_wetdry_cells(int nTasks, t_arrays *arrays);

//__global__ void c_update_boundaries

__global__ void g_compute_mass_error(t_arrays *arrays);

__global__ void g_reconstruct_active_elements(int nTasks, t_arrays *arrays);

__device__ void g_add_active_cells(t_arrays *arrays, int id);

__device__ void g_add_active_walls(t_arrays *arrays, int id);





