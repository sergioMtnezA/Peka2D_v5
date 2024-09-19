#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#include <time.h>

#ifdef _WIN32
    #include <sys/types.h>
    #include <assert.h>
#endif

#include "define.h"
#include "structs.h"
#include "mesh.h"

#if SET_SIMGPU
    #include "writeInOut.cuh"
#else
    #include "writeInOut.h"
#endif

#ifndef hy_mesh_
#define hy_mesh_

#define STR_SIZE 1024

// BOUNDARY CONDITIONS
#define HYD_N_OBT 200
#define HYD_N_OBC_TYPES 12

//Boundary labels
#define HYD_CLOSED 0

#define HYD_INFLOW_Q 1
#define HYD_INFLOW_HZ 2
#define HYD_INFLOW_QHZ 3

#define HYD_OUTFLOW_GAUGE 11
#define HYD_OUTFLOW_HZ 12
#define HYD_OUTFLOW_FREE 13
#define HYD_OUTFLOW_FR 14
#define HYD_OUTFLOW_NORMAL 15


typedef struct Peka2D_Setup_ Peka2D_Setup;
typedef struct Peka2D_Run_ Peka2D_Run;
typedef struct Peka2D_NodeBoundary_ Peka2D_NodeBoundary;
typedef struct Peka2D_OBCWalls_ Peka2D_OBCWalls;


struct Peka2D_Run_{
    int release;      // Release number

    //components
    int rain;         // 1 if there is rain
    int sediment;     // 1 if there is sediment
    int weirs;        /**< 1 if there are weirs*/
    int sources;      /**< 1 if there are source/sink terms*/
    int dambreach;    /**< 1 if there are dambreach*/
    int wind;

    //controls
    int writeMass;
    int extremes;
    int crossSection; // 1 if there are cross sections
    int profile;      // 1 if there are profiles
    int obs;          // 1 to print out probes

    int nIterOut;
    double CFL;
    double dtDump;    // Should not be used by PEKA
    double dtOut;      // Output times
    double tLimit;    // Max Simulation time
    double tInitial;    // Initial Simulation time

    int iniType;      /**< Initial condition switch*/
    double initialWSE;
    int hotStart;     /**< Hot start switch*/
    int indexHS;

    double XnMan;     // Manning's n multiplier
    double hMin;      // Minimum depth for computations

    int ncores;
    int gdevice;

    char hsfile[STR_SIZE]; 

    int itrash;
    double dtrash;   

};

struct Peka2D_NodeBoundary_{
      int n;      // Nnode
      int countInlet;
      int countOutlet;
      int *type;
      int *obcId;
};

struct Peka2D_OBCWalls_{
      int *wall;
      int n;
};

struct Peka2D_Setup_{
    Peka2D_Run pkrun;
    Peka2D_NodeBoundary *pknode;
    Peka2D_OBCWalls *IOBC;
    Peka2D_OBCWalls *OOBC;    
};


////////////////////////////////////////////////////////////////
EXPORT_DLL int loadControlParameters(
    Peka2D_Setup *pksetup, 
    t_parameters *spar, 
    t_mesh *mesh, 
    t_message *msg);
/*----------------------------*/

////////////////////////////////////////////////////////////////
EXPORT_DLL int readControlDataFile(
    char *filename,
    Peka2D_Setup *pksetup,
    t_message *msg);
/*----------------------------*/

////////////////////////////////////////////////////////////////
EXPORT_DLL int setControlParameters(
    Peka2D_Setup *pksetup, 
    t_parameters *spar, 
    t_message *msg);
/*----------------------------*/

////////////////////////////////////////////////////////////////
EXPORT_DLL int loadMeshData(
    Peka2D_Setup *pksetup, 
    t_parameters *spar, 
    t_mesh *mesh, 
    t_message *msg);
/*----------------------------*/

////////////////////////////////////////////////////////////////
int IsMeshAllocated(
    t_mesh *mesh, 
    t_message *e);
/*----------------------------*/

////////////////////////////////////////////////////////////////
int readMeshFile(
    char *filename,
    Peka2D_Setup *pksetup,
    t_mesh *mesh, 
    t_message *e);
/*----------------------------*/

////////////////////////////////////////////////////////////////
int setInitialState(
    Peka2D_Setup *pksetup,
    t_parameters *spar,
    t_mesh *mesh, 
    t_message *e);
/*----------------------------*/

////////////////////////////////////////////////////////////////
int loadBoundaryConditions(
    Peka2D_Setup *pksetup, 
    t_parameters *spar,    
    t_mesh *mesh, 
    t_message *e);
/*----------------------------*/

////////////////////////////////////////////////////////////////
int readOpenBoundaryNodes(
    char *filename,
    int nnodes,
    Peka2D_NodeBoundary *nbc, 
    t_message *e);
/*----------------------------*/

////////////////////////////////////////////////////////////////
int createOpenBounds(
    Peka2D_NodeBoundary *nbc,
    Peka2D_OBCWalls *IOBC,
    Peka2D_OBCWalls *OOBC,
    t_parameters *spar, 
    t_mesh *mesh,    
    t_message *e);
/*----------------------------*/

////////////////////////////////////////////////////////////////
int readOpenBoundaryFile(
    char *filename,
    Peka2D_NodeBoundary *nbc,
    Peka2D_OBCWalls *IOBC,
    Peka2D_OBCWalls *OOBC,
    t_parameters *spar, 
    t_mesh *mesh,
    t_message *e);
/*----------------------------*/


#endif //end #ifndef hy_mesh_
