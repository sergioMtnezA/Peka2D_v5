#include "define.h"

#ifndef s_mesh_
#define s_mesh_

#define STR_SIZE 1024

typedef struct t_message_ t_message;

typedef struct t_parameters_ t_parameters;
typedef struct t_mesh_ t_mesh;

typedef struct t_c_cell_ t_c_cell;
typedef struct l_c_cells_ l_c_cells;

typedef struct t_g_cell_ t_g_cell;
typedef struct l_g_cells_ l_g_cells;

typedef struct t_node_ t_node;
typedef struct l_nodes_ l_nodes;

typedef struct t_wall_ t_wall;
typedef struct l_wall_ l_wall;

typedef struct t_edge_ t_edge;
typedef struct l_edges_ l_edges;

typedef struct t_bound_ t_bound;

typedef struct t_solute_ t_solute;
typedef struct l_solutes_ l_solutes;

typedef struct t_arrays_ t_arrays;
typedef struct t_cuPtr_ t_cuPtr;

typedef struct t_timers_ t_timers;

/**
 * @brief Compilation and execution messages
*/
struct t_message_{
	char logFile[1024]; /**< @brief Log file path*/
	int error; /**< If Error > 0, Execution will be aborted */
	char errorProperty[1024]; /**< @brief Error string */	
	int nWarning;
	char **warning;
	int nMsgL0;
	char **MsgL0;
	int nMsgL1;
	char **MsgL1;
	int nMsgL2;
	char **MsgL2;
	int nMsgL3;
	char **MsgL3;
};


/**
 * @brief Simulation execution paramenters
*/
struct t_parameters_{
	char proj[1024];
	char dir[1024];
	int ncores, gpuid;
	double ti,tf,CFL;

	int writeMass;
	int writeExtremes;
	int indexWriteHotstart;

    int nIterOut;
	double dtOut;
	double dtDump;
	double minh;

};


/**
 * @brief Geometrical mesh and flow data
*/
struct t_mesh_{
	int id;
    int NCwall;
    int nnodes;
	int ncells;
	int ncellsBound;
	int nw_calc;
	int nw_bound;

    //mesh structure
	l_c_cells *c_cells;
	l_g_cells *g_cells;
    l_nodes *nodes;

	l_wall *w_bound;
	l_wall *w_calc;
	
    l_c_cells *b_cells; // Only the boundary Cells

	int nInlet;
	int nTotalCellsIn;
	t_bound *in;
	int nOutlet;
	int nTotalCellsOut;
	t_bound *out;

    double minZ, maxZ;

	//solutes
	int nSolutes;
	l_solutes *solutes;

};


/**
 * @brief Calculus cells: flow variables in cells
*/
struct t_c_cell_{
	int id;
    double z;
	double h;
	double hu;
	double hv;
	double u;
	double v;
	double modulou;

	double hini;
    double zini;

	double nMan;

	t_g_cell *geom;

	//Advection-Diffusion Transport 
	double *hphi;
	double *phi;
};


/**
 * @brief List of calculus cells @ref t_c_cell_
*/
struct l_c_cells_{
	int size;
	int n;
	t_c_cell* cells;
	t_c_cell** pcells;
};


/**
 * @brief Geometry cells: mesh characteristics in cells
*/
struct t_g_cell_{
	int id;
    int isBound; 	// type cell identifier 
					// 0 = inner cell
					// 99 = closed boundary cell
					// negative = inlet cell (- ID inlet)
					// positive = outlet cell (+ ID outlet)
	double area;
	double center[3];
    int nneig; //número de vecinas, descontando los contornos
	int nbound; //número de contornos
	t_node *nodes[4];
	t_c_cell *neigcell[4];
	t_wall *neigwall[4];
	t_edge *wall;
    double Lwallmax;
};


/**
 * @brief List of geometry cells @ref t_g_cell_
*/
struct l_g_cells_{
	int size;
	int n;
	t_g_cell *cells;
};


/**
 * @brief Geometry mesh node
*/
struct t_node_{
	int id;
	double x,y,z;
};


/**
 * @brief List of geometry nodes @ref t_node_
*/
struct l_nodes_{
	int size;
	int n;
	double xmin,ymin,xmax,ymax;
	t_node *nodes;
};


/**
 * @brief Wall characteristics for computation
*/
struct t_wall_{
    int idWall;

    int node1;
    int node2;
    // This value represents, for c1, the identifier {0,1,2,{3}} of the wall
    // in order to select where to write at dh,dhu,dhv of the cell
    int i;
    // The same for c2
    int j;

    int id; // Identificador de la celda i
    int id2;// Identificaador de la celda j
    double lwall;
    double normal[2];
    double deltaX;
    double distNormal;//length between the centers of the cells
    double length;//length of the edge of each cell  
    double distCenter;
    double ncx, ncy;  
   
    // Punteros a las celdas i,j (0,1) de la pared
    t_c_cell *ccells[2];
    t_g_cell *gcells[2];

    int idBound; 
    int typeOfBound;	//type wall identifier
						//-1 = inner wall
						// 0 = closed bound wall
						// 2 = inlet bound wall
						// 3 = outlet bound wall
};


/**
 * @brief List of walls for computation @ref t_wall_
*/
struct l_wall_{
	int size;
	int n;
	t_wall *wall;
};


/**
 * @brief Geometrical edges for mesh topology construction
*/
struct t_edge_{
	int id;
	double normal[3];
	//t_c_cell *cells[2];
	double length;
	int isBoundary;
};


/**
 * @brief List of geometrical edges for mesh construction
*/
struct l_edges_{
	int size;
	int n;
	t_edge *walls;
};


/**
 * @brief Open boundary geometry and flow conditions
*/
struct t_bound_{
    char idname[1024];
    int type;               /**< Hydraulic boundary type*/
    int n;                  /**< Number of data points*/
    int ncellsBound;        /**< Number of cells in boundary*/
    int ncellsInner;        /**< Number of cells that are neighbors to boundary cells*/
    t_c_cell *cellzMin;     /**< Pointer to boundary with minimum elevation*/
    t_wall **wallBound;     /**< Array of pointers to boundary walls (of size ncellsBound)*/
    t_c_cell **cellBound;   /**< Array of pointers to boundary cells associated to boundary walls (of size ncellsBound) */
    t_c_cell **cellInner;   /**< Array of pointers to boundary cell neighbors (of size ncellsInner) */
    double *t;
    double *q;
    double *hZ;
    double *Fr;
    double **phi;
    double maxFroude;
    double totalLength;     /**< totalLength is the length of the boundary*/
    double totalArea;       /**< totalArea is the area (in plant view) of the boundary*/
    double zmin;            /**< zmin is the minimum value of z along the boundary*/
    double levelControl;    /**< levelControl accounts for the water level to be imposed in the case of extreme boundaries (dry domain i.e.) */
    int flagInitialize;    /**< flag to determine if the inlet bound has been got inside the initialize inlet outlet function */
    double qMass,qDischarge; //qMass is the mass of water that enters/leaves in form of water depth and qDischarge in form of discharge. The sum of this two quantities (times dt) represents the water volume that enters/leaves the boundary
    t_node node[2],normal;
    // *(n1)-------------*
    // |		     |
    // |     	     |
    // *-------------(n2)*
    int iAct;	// Nos dice el indice, en la tabla de tiempos, correspondiente al tiempo actual
    double qliquid;
    double qsolid;
};


/**
 * @brief Passive solute
*/
struct t_solute_{
	char name[STR_SIZE]; /**< @brief Solute name*/
	int typeDiff; /**< @brief Diffusion model*/
	double k_xx,k_yy; /**< @brief Longitudinal and transversal diffusion coefficients*/
	double maxConc;
};


/**
 * @brief List of passive solutes
*/
struct l_solutes_{
	int n;
	int flagDiffussion; //flag to determine if it is necessary to go inside diffusion subroutines
	t_solute *solute;
};



/**
 * @brief Arrange all the domain parameters, run controls and data arrays required for computation
*/
struct t_arrays_{

    //parameters
	int ncores; /**< @brief Number of CPU threads selected */
	int gpuid; /**< @brief GPU device ID selected by user*/
	double ti; /**< @brief Initial time for simulation*/
	double tf; /**< @brief Final time for simulation*/
	double CFL; /**< @brief CFL restriction for time step limitation */

	int writeMass; /**< @brief Write in output file mass balance information */
	int writeExtremes; /**< @brief Compute and write maximum values for the variables */
	int indexWriteHotstart; /**< @brief Initial output index for the hortstart initialization */

    int nIterOut; /**< @brief Number of iterations for sreen output */
	double dtOut; /**< @brief Temporal frequency for writing state resutls */  
	double dtDump; /**< @brief Temporal frequency for dump components integrations */
	double minh; /**< @brief Minimun flow depth for neglecting cell velocity */

	int NCwall; /**< @brief Number of walls per cell :: NCWALL*/
	int ncells; /**< @brief Number of cells in the domain :: NCELLS */
	int nw_calc; /**< @brief Number of internal (calculus) walls :: NWCALC */
	int nw_bound; /**< @brief Number of boundary walls */
	int nInlet; /**< @brief Number of inlets in the domain */
	int nOutlet; /**< @brief Number of outlets in the domain */


    //computation controls
    double t; /**< @brief Run-control time mark */
	double dt; /**< @brief Run-control time step */

	int nIter; /**< @brief Run-control iteration*/
	int indexOut; /**< @brief Run-control state output index*/
	int indexDump; /**< @brief Run-control component output index*/
    int dumpComponent; /**< @brief Run-control flag to dump components in the current iteration*/
	int dumpState; /**< @brief Run-control flag to dump state in the current iteration*/

	double massOld; /**< @brief Run-control integrated volume at the beginning of the time step*/
	double massNew; /**< @brief Run-control current integrated volume */
    double massError; /**< @brief Run-control mass error computed at the end of the time stop */

	double massIn; /**< @brief Run-control total inflow volume accumulated during the simulation*/
	double massOut; /**< @brief Run-control total outflow volume accumulated during the simulation*/
	double qTotalIn; /**< @brief Run-control inflow discharge inegrated for all the inlet boundaries*/
	double qTotalOut; /**< @brief Run-control outflow discharge inegrated for all the outlet boundaries*/

	//solute controls
	int nSolutes; /**< @brief Number of solutes */

	
	//ARRAY DE CELDA
    int nActCells; /**< @brief Run-control number of active cells*/
    int *actCells; /**< @brief [NCELLS] Active cell indexes*/
    int *activeC; /**< @brief [NCELLS] Local flags for active cells :: [0] Unactive cell - [1] Active cell*/	

    int *cidx; /**< [NCELLS] @brief Cell indexes*/
    int *nneig; /**< [NCELLS] @brief Number of neighboring cells*/

	double *z; /**< @brief [NCELLS] Bed elevation in cells*/
	double *h; /**< @brief [NCELLS] Flow depth in cells*/
	double *hu; /**< @brief [NCELLS] X-discharge in cells*/
	double *hv; /**< @brief [NCELLS] Y-discharge in cells*/
    double *u; /**< @brief [NCELLS] X-velocity in cells*/ 
	double *v; /**< @brief [NCELLS] Y-velocity in cells*/
    double *modulou; /**< @brief [NCELLS] Velocity modulus in cells*/
	double *sqrh; /**< @brief [NCELLS] Square depth of flow depth in cells*/

	double *area; /**< @brief [NCELLS] Cell area*/
	double *nman; /**< @brief [NCELLS] Manning roughness coefficient nMan [s m^(-1/3)] in cells*/
    double *SOX; /**< @brief [NCELLS] X-direction bed slope*/
	double *SOY; /**< @brief [NCELLS] Y-direction bed slope*/

	double *mass; /**< @brief [NCELLS] Mass in cells*/
	

	//ARRAY DE CELDAS*NPAREDES
    int nWallCell;
	double *dh; /**< @brief [NCELLS*NCWALL] Wall-contributions to cell mass*/
	double *dhu;/**< @brief [NCELLS*NCWALL] Wall-contributions to cell X-momentum*/
	double *dhv;/**< @brief [NCELLS*NCWALL] Wall-contributions to cell X-momentum*/

    int *solidWallByCell; /**< @brief [NCELLS*NCWALL] Local flag for solid wall behaviour*/
    int *neighCell; /**< @brief [NCELLS*NCWALL] Neigboring cell indexes. Boundary neigbour marked with index -1*/
    int *neighWall; /**< @brief [NCELLS*NCWALL] Neigboring wall IDs :: Internal calculus wall - Bound walls*/
	
	double *normalXbyCell; /**< @brief [NCELLS*NCWALL] X-nomall neighboring wall*/
	double *normalYbyCell; /**< @brief [NCELLS*NCWALL] Y-nomall neighboring wall*/
	

    //ARRAY DE PAREDES INTERNAS DE CALCULO
    int nActWalls; /**< @brief Run-control number of active calculus walls*/
    int *actWalls; /**< @brief [NWCALC] Active calculus walls indexes*/
    int *activeW; /**< @brief [NWCALC] Local flags for active calculus walls :: [0] Unactive wall - [1] Active wall*/

    int *widx; /**< @brief [NWCALC] Calculus walls IDs*/

	int *idx1; /**< @brief [NWCALC] Neighboring left-cell ID for calculus wall*/
	int *idx2; /**< @brief [NWCALC] Neigboring right-cell ID for calculus wall*/
	int *idw1; /**< @brief [NWCALC] Wall-ordering flag in left-cell*/
	int *idw2; /**< @brief [NWCALC] Wall-ordering flag in right-cell*/
	    
	double *normalX; /**< @brief [NWCALC] X-normal for calculus wall*/
	double *normalY; /**< @brief [NWCALC] Y-normal for calculus wall*/
	double *deltaX; /**< @brief [NWCALC] Distance factor for time step computation*/
	double *length; /**< @brief [NWCALC] Wall length*/
	double *distNormal; /**< @brief [NWCALC] Distance between neighboring cell centers projected on the wall-normal direction*/
    double *distCentX;  /**< @brief [NWCALC] Distance between neighboring cell centers projected on the X-direction*/
	double *distCentY; /**< @brief [NWCALC] Distance between neighboring cell centers projected on the Y-direction*/
    
    double *nman2wall; /**< @brief [NWCALC] Square of the averaged Manning parameter at the calculus wall*/
    double *gp; /**< @brief [NWCALC] Projected gravity acceleration along the bed-normal direction :: gp=g*cos^2(B)*/

	int *typeOfBound; /**< @brief [NWCALC] Type of bound condition at the wall*/
    int *solidWall; /**< @brief [NWCALC] Local flag of solid wall behaviour for calculus wall*/
	
    double *qnormalL; /**< @brief [NWCALC] Upwind-left volume rate for the calculus wall*/
	double *localDt; /**< @brief [NWCALC] Local time step limitation for the calculus wall*/


	// BOUNDARY ARRAYS ///////////////////////////	
	double *normalXIn, *normalYlIn;
	double *normalXOut, *normalYOut;	

	int nTotalCellsIn;	
	int *cidxIn;
	int *idBoundIn;	
	double *lengthwallIn;

	int nTotalCellsOut;	
	int *cidxOut;
	int *idBoundOut;
	double *lengthwallOut;		


	// SOLUTE ARRAYS ///////////////////////////	
	#if SET_SOLUTE
		int flagDiffusion; /**< @brief Solute diffusion activation flag */

		//solute data
		int *typeDiff; //nsol
		double *k_xx, *k_yy; //nsol	

		//solute conservative
		double *hcsol; /**< @brief [NSOL x NCELLS] Solute mass in cells*/
		double *csol; /**< @brief [NSOL x NCELLS] Solute concentration in cells*/	

		//solute contributions
		double *dhcsol;/**< @brief [NSOL*NCELLS*NCWALL] Wall-contributions to cell solute mass*/
	#endif

};


/**
 * @brief Arrange the pointers to CUDA memory for all the run-control and data-arrays in @ref t_arrays
 * 
 * The definition of the variables and memory size is the same as that in the @ref t_arrays description.
 * This structure is created in the RAM memory but the CUDA memory for data-arrays is allocated using these pointers.
*/
struct t_cuPtr_{

	int *index, *check;

    double *t, *dt;

	int *nActCells, *nActWalls;

    int *nIter, *indexOut, *indexDump;
    int *dumpComponent, *dumpState;

    double *massOld, *massNew;
    double *massError;

	double *qTotalIn, *qTotalOut;	

	// GPU COMPUTATION ARRAYS ////////////////////////////////////////////
    //cells
	int *cidx, *nneig;
	int *activeC, *actCells;
	double *z, *h, *hu, *hv, *u, *v, *modulou, *sqrh;
	double *area, *nman, *SOX, *SOY;
	double *mass;

	//cells*NCwall
	double *dh, *dhu, *dhv;
	int *solidWallByCell, *neighCell, *neighWall;
	double *normalXbyCell, *normalYbyCell;

	//internal walls
	int *activeW, *actWalls;
	int *widx;
	int *idx1, *idx2, *idw1, *idw2;
	double *normalX, *normalY, *deltaX, *length, *distNormal, *distCentX, *distCentY;
	double *nman2wall, *gp;
	int *typeOfBound, *solidWall;
	double *qnormalL, *localDt;


	// BOUNDARY ARRAYS ///////////////////////////	
	int *cidxIn;
	int *idBoundIn;	
	double *normalXwallIn, *normalYwallIn;
	double *lengthwallIn;

	int *cidxOut;
	int *idBoundOut;
	double *normalXwallOut, *normalYwallOut;
	double *lengthwallOut;		


	// SOLUTE ARRAYS ///////////////////////////
	#if SET_SOLUTE
		//solutes
		int *typeDiff;
		double *k_xx, *k_yy; 

		//nsolutes*cells
		double *hcsol, *csol;

		//nsolutes*cells*NCwall
		double *dhcsol; 
	#endif

};


/**
 * @brief Timers for accounting computational time
*/
struct t_timers_{
	double total=0.0;
	double loading=0.0;
	double init=0.0;
	double initGPU=0.0;
	double computeSim=0.0;
	double wallCalculus=0.0;
	double cellUpdating=0.0;
	double boundConditon=0.0;
	double wetDryFix=0.0;
	double memoryTransfer=0.0;
	double writeOut=0.0;
};


#endif
