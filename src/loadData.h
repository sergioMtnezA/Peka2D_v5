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

// OUR LABELS. REALMENTE CREO QUE NO SE UTILIZAN
#define HYD_CLOSED 0

#define HYD_INFLOW_Q 1
#define HYD_INFLOW_HZ 2
#define HYD_INFLOW_QHZ 3

#define HYD_OUTFLOW_GAUGE 11
#define HYD_OUTFLOW_HZ 12
#define HYD_OUTFLOW_FREE 13
#define HYD_OUTFLOW_FR 14
#define HYD_OUTFLOW_NORMAL 15

// HYDRONIA LABELS FOR BOUNDARY CONDITIONS
#define HYD_OBC_Q 6
#define HYD_OBC_QHZT 5
#define HYD_OBC_QHZ 9
#define HYD_OBC_HZ 1
#define HYD_OBC_HZD 17 //mesh orientation
#define HYD_OBC_FREE 11
#define HYD_OBC_FR 13
#define HYD_OBC_NORMAL 12
#define HYD_OBC_RT 14
#define HYD_OBC_HZS 18 ///// NUEVA CONDICIÓN DE SEDIMENTO/SOLUTO, LES AÑADO UN 8
#define HYD_OBC_QQS 26
#define HYD_OBC_WEIR 30 //cc vertedero
#define HYD_OBC_2D0D 31//cc modelo 2d0d
#define HYD_OBC_WEIRPID 32 //cc weir with pid
#define HYD_OBC_2D0DPID 33 //cc modelo2d0d with pid
#define HYD_OBC_QWEIR 40


// INITIAL CONDITIONS
#define HYD_INIT_WSE 0
#define HYD_INIT_DRY 1
#define HYD_INIT_VAR 2
#define HYD_INIT_WET -9999


typedef struct Hydronia_OBCWalls_ Hydronia_OBCWalls;
typedef struct Hydronia_Run_ Hydronia_Run;
typedef struct Hydronia_NodeBoundaryConditions_ Hydronia_NodeBoundaryConditions;
typedef struct Hydronia_NodeInitialConditions_ Hydronia_NodeInitialConditions;
typedef struct Hydronia_Culvert_ Hydronia_Culvert;
typedef struct Hydronia_CulvertGroup_ Hydronia_CulvertGroup;
typedef struct Hydronia_Pier_ Hydronia_Pier;
typedef struct Hydronia_PiersGroup_ Hydronia_PiersGroup;
typedef struct Hydronia_RainEvap_  Hydronia_RainEvap;
typedef struct Hydronia_LocalRainEvap_  Hydronia_LocalRainEvap;
typedef struct Hydronia_LocalInf_  Hydronia_LocalInf;
typedef struct Hydronia_Wind_  Hydronia_Wind;
typedef struct Hydronia_ManZone_  Hydronia_ManZone;
typedef struct Hydronia_Sediment_ Hydronia_Sediment;
typedef struct Hydronia_SedimentGroup_ Hydronia_SedimentGroup;
typedef struct Hydronia_SourceSink_ Hydronia_SourceSink;
typedef struct Hydronia_SourceSinkGroup_ Hydronia_SourceSinkGroup;
typedef struct Hydronia_Weir_ Hydronia_Weir;
typedef struct Hydronia_WeirsGroup_ Hydronia_WeirsGroup;
typedef struct Hydronia_Rill_ Hydronia_Rill;
typedef struct Hydronia_RillsGroup_ Hydronia_RillsGroup;
typedef struct Hydronia_IRT_ Hydronia_IRT;
typedef struct Hydronia_IRTGroup_ Hydronia_IRTGroup;
typedef struct Hydronia_Gate_ Hydronia_Gate;
typedef struct Hydronia_GatesGroup_ Hydronia_GatesGroup;
typedef struct Hydronia_Dambreach_ Hydronia_Dambreach;
typedef struct Hydronia_DambreachGroup_ Hydronia_DambreachGroup;
typedef struct Hydronia_Probe_ Hydronia_Probe;
typedef struct Hydronia_ProbeGroup_ Hydronia_ProbeGroup;
typedef struct Hydronia_NodeLink_ Hydronia_NodeLink;   // NEW
typedef struct Hydronia_NodeLinkGroup_ Hydronia_NodeLinkGroup;   // NEW
typedef struct Hydronia_Plt_ Hydronia_Plt;
typedef struct Hydronia_Profile_ Hydronia_Profile;
typedef struct Hydronia_Solute_ Hydronia_Solute;
typedef struct Hydronia_ProfileGroup_ Hydronia_ProfileGroup;
typedef struct Hydronia_CrossSection_ Hydronia_CrossSection;
typedef struct Hydronia_CrossSectionGroup_ Hydronia_CrossSectionGroup;
typedef struct Hydronia_BoundaryConditionData_ Hydronia_BoundaryConditionData;
typedef struct Hydronia_StageDischarge_ Hydronia_StageDischarge;
typedef struct Hydronia_CulvertDepthDischarge_ Hydronia_CulvertDepthDischarge;
typedef struct Hydronia_MudSetup_ Hydronia_MudSetup;
typedef struct Hydronia_SoluteGroup_ Hydronia_SoluteGroup;
typedef struct Hydronia_WQMGroup_ Hydronia_WQMGroup;
typedef struct Hydronia_Bridge_ Hydronia_Bridge;
typedef struct Hydronia_BridgesGroup_ Hydronia_BridgesGroup;
typedef struct Hydronia_Setup_ Hydronia_Setup;

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


// Lee la malla
EXPORT_DLL int ReadRunControlDataFile(char *filename, Hydronia_Run *run,t_message *e, char *nameFileTemp);

EXPORT_DLL int createDateOfCompilation(char *datecomp);

EXPORT_DLL int verifyHydroniaInitial(char *directory, char *filename,t_mesh *mesh,Hydronia_Setup *setup, t_message *e);

EXPORT_DLL int LoadHydroniaMesh(char *directory, char *filename,Hydronia_Run *run,t_mesh *mesh, t_message *e);

EXPORT_DLL int LoadHydroniaSetup(char *directory, char *name,t_mesh *mesh, Hydronia_Setup *setup, t_message *e);

EXPORT_DLL int ConvertHydroniaToPEKASetup(Hydronia_Setup *hs, t_parameters *Pm,t_mesh *mesh, t_message *e);

EXPORT_DLL int write_mass_balance_hydronia(t_mesh *mesh, double t, int flagIni, int units);

int Read_Hydronia_Mesh_Size(char *filename, int *Nnode, int *Ncell, t_message *e);

int Read_Hydronia_Mesh(char *filename,Hydronia_Run *run, t_mesh *mesh, Hydronia_NodeBoundaryConditions *nbc, Hydronia_NodeInitialConditions *nic, t_message *e);

int IsMeshAllocated(t_mesh *mesh,t_message *e);

int BuildBoundariesFromNodeDataOBCP(char *directory, char *basename,Hydronia_Run *run, t_mesh *mesh, Hydronia_NodeBoundaryConditions *nbc, t_message *e);
int ReadManningDataFile(char *directory, char *filename, Hydronia_Run *run, Hydronia_ManZone *MZ, t_message *e);

EXPORT_DLL int removeOldOUTFILESAndWriteNew(char *directory, char *proj, int index);
EXPORT_DLL int write_mesh_vtk_hydronia(char *filename, t_mesh *mesh, t_message *msg, int units);
EXPORT_DLL int write_cell_time_textout_files(char *directory, t_mesh *mesh, double t, int units);

#endif //end #ifndef hy_mesh_
