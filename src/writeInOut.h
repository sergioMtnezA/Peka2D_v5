#include "define.h"
#include "structs.h"
#include "utilities.h"

EXPORT_DLL void Notify(char *msg, int type, t_message *e);

EXPORT_DLL int write_vtk_mesh(char *filename, t_mesh *mesh, t_message *msg);

EXPORT_DLL int create_computation_files(char *path, t_message *msg);

EXPORT_DLL int write_massBalance(char *path, t_arrays *arrays, t_message *msg);

EXPORT_DLL int write_timers(char *path, t_timers timers, t_message *msg);

EXPORT_DLL void dump_screen_info(t_arrays *arrays, t_message *msg);

EXPORT_DLL int write_vtk_state(char *filename, t_mesh *mesh, t_arrays *arrays, t_message *msg);

EXPORT_DLL int write_mesh_structure_in_vtk_binary(FILE *fp, t_mesh *mesh);

EXPORT_DLL int write_mesh_structure_in_vtk_ascii(FILE *fp, t_mesh *mesh);

EXPORT_DLL void write_Dscalar_in_vtk(FILE *fp, double val);

EXPORT_DLL void write_Iscalar_in_vtk(FILE *fp, int val);

EXPORT_DLL void write_Dvector_in_vtk(FILE *fp, double val1, double val2);