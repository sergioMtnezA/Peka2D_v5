#include <search.h>
#include <math.h>
#include "omp.h"

#include "define.h"
#include "structs.h"

#if SET_SIMGPU
    #include "writeInOut.cuh"
#else
    #include "writeInOut.h"
#endif

//#include "writeInOut.cuh"

#ifndef i_mesh_
#define i_mesh_

typedef struct t_pair_ t_pair;

struct t_pair_{
    int included;
	int id1;
	int id2;
	int idW;
	int idW2;
	int idC;
	int n1,n2;
};

int new_lg_cells(l_g_cells *lcells,int n);

int new_lc_cells(l_c_cells *lcells,int n);

int new_g_cell(t_g_cell *cell, int id);

int new_c_cell(t_c_cell *cell, int id);

int add_g_cell(t_g_cell *cell, l_g_cells *lcells);

int add_c_cell(t_c_cell *cell, l_c_cells *lcells);

int new_lnodes(l_nodes *lnodes,int n);

int new_node(t_node *node, int id);

int add_node(t_node *node, l_nodes *lnodes);

int new_lwalls(l_wall *lwalls, int n);

int calc_norm(l_g_cells *lcells, l_c_cells *lccells, int nwall, t_message *msg);

int cons_lpocells_qsort(t_mesh *mesh, t_message *msg);

int cmpp(const void *id1, const void *id2);

int check_angle(l_g_cells *lcells,t_message *msg);

int build_wall_inlet(t_mesh *mesh, t_bound *iinn, t_wall *w, int i, t_message *msg);

int build_inner_inlet(t_mesh *mesh, t_bound *iinn, int i, t_message *msg);

int build_wall_outlet(t_mesh *mesh, t_bound *outt, t_wall *w, int i,t_message *msg);

int build_inner_outlet(t_mesh *mesh,t_bound *outt, int i, t_message *msg);

#endif
