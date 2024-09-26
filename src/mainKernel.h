#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#ifdef _WIN32
	#if APP_MODE == 1
	#using <mscorlib.dll>
	#using <System.dll>
	#using <System.Windows.Forms.dll>
	#using <System.Drawing.dll>
	#include "..\dev\p2dGUI.h"
	#endif
#endif

#include <time.h>
#include "define.h"
#include "structs.h"
#include "loadData.h"
#include "initialize.h"

#if SET_SIMGPU
	#include "writeInOut.cuh"
	#include "calculate.cuh"
#else
	#include "writeInOut.h"
    #include "calculate.h"
#endif


extern Peka2D_Setup *pksetup;
extern t_parameters spar;
extern t_mesh *mesh;
extern t_message *msg;

int runMainKernel (int argc, char * argv[]);

/**
 * @brief This function reads raster files and stores the information in the struct-type variable called msg.
 */
int loadSimulationDomain(
	Peka2D_Setup *pksetup, 
	t_parameters *spar, 
	t_mesh *mesh,
	t_timers *timers,
	t_message *msg);

/**
 * @brief This function loads mesh-type structures to array-type structures.
 */
int initializeComputationArrays(
	t_parameters spar, 
	t_mesh *mesh,
    t_arrays *carrays,   
	t_timers *timers,
	t_message *msg);



