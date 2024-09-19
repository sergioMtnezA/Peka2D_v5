#include "define.h"

#ifndef s_mesh_
#define s_mesh_

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

typedef struct t_arrays_ t_arrays;
typedef struct t_cuPtr_ t_cuPtr;

typedef struct t_timers_ t_timers;


struct t_message_{
	char logFile[1024];     /**< Log file */
	int error;                    /**< If Error > 0, Execution will be aborted */
	char errorProperty[1024];     /**< Error string */	
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

struct t_mesh_{
	int id;
    int NCwall;
    int nnodes;
	int ncells;
	int ncellsBound;
	int nw_calc;
	int nw_bound;

    int nInlet;
    int nOutlet;

    //mesh structure
	l_c_cells *c_cells;
	l_g_cells *g_cells;
    l_nodes *nodes;

	l_wall *w_bound;
	l_wall *w_calc;
	
    l_c_cells *b_cells; // Only the boundary Cells

	t_bound *in;
	t_bound *out;

    double minZ, maxZ;
};

struct t_c_cell_{
	int id;
    double z;
	double h;
	double hu;
	double hv;
	double u;
	double v;
	double modulou;

	double *hphi;
	double *phi;

	double hini;
    double zini;

	double nMan;

	t_g_cell *geom;
};

struct l_c_cells_{
	int size;
	int n;
	t_c_cell* cells;
	t_c_cell** pcells;
};

struct t_g_cell_{
	int id;
    int isBound;//isBound is and identifier to know if a cell is bound or not. 0->normal cell, negative=inlet cell, positive->outletcell
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

struct l_g_cells_{
	int size;
	int n;
	t_g_cell *cells;
};

struct t_node_{
	int id;
	double x,y,z;
};

struct l_nodes_{
	int size;
	int n;
	double xmin,ymin,xmax,ymax;
	t_node *nodes;
};

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
    int typeOfBound;
    // 1 si es esquina (velocidad 0)
    // 0 si es pared (qx o qy =0) Incluye normal
    // 2 si es celda de entrada
    // 3 si es celda de salida
    // 4 si es celda de contorno interna
    // 5 si es puente

};

struct l_wall_{
	int size;
	int n;
	t_wall *wall;
};

struct t_edge_{
	int id;
	double normal[3];
	//t_c_cell *cells[2];
	double length;
	int isBoundary;
};

struct l_edges_{
	int size;
	int n;
	t_edge *walls;
};


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


struct t_arrays_{

    //parameters
	int ncores, gpuid;
	double ti,tf,CFL;

	int writeMass;
	int writeExtremes;
	int indexWriteHotstart;

    int nIterOut;
	double dtOut;
	double dtDump;
	double minh;  

	int NCwall;
	int ncells;		//-->NC
	int nw_calc;	//-->NW
	int nw_bound;
	int nInlet,nOutlet;


    //computation controls
    double t, dt;

	int nIter, indexOut, indexDump;
    int dumpComponent, dumpState;

	double massOld, massNew;
    double massError;

	double massIn, massOut;
	double qTotalIn, qTotalOut;


	//ARRAY DE CELDA
    int nActCells;
    int *actCells;   //NC
    int *activeC;    //NC	

    int *cidx;
    int *nneig;     //NC

	double *z;		//NC
	double *h;		//NC
	double *hu;		//NC
	double *hv;		//NC	
    double *u,*v; 		//NC equivalent to c1->u,c1->v respectively
    double *modulou;
	double *sqrh;   //NC equivalent to c1->raizh	

    
	double *area;	//NC
	double *nman;	//NC
    double *SOX,*SOY; //NC

	double *mass;	//NC
	

	//ARRAY DE CELDAS*NPAREDES
    int nWallCell;
	double *dh;		//NC * NCwall
	double *dhu;	//NC * NCwall
	double *dhv;	//NC * NCwall

    int *solidWallByCell;	//NC * NCwall
    int *neighCell;		//NC * NCwall It specifies the idenfifier of the neigcell of i at edge j by means of neigOfCell[j*ncell+i]
    int *neighWall;		//NC * NCwall It specifies the idenfifier of the neigcell of i at edge j by means of neigOfCell[j*ncell+i]
	
	double *normalXbyCell; //NC * NCwall
	double *normalYbyCell; //NC * NCwall
	

    //ARRAY DE PAREDES INTERNAS DE CALCULO
    int nActWalls;
    int *actWalls;      //NW
    int *activeW;       //NW

    int *widx;			//NW

	int *idx1;			//NW
	int *idx2;			//NW
	int *idw1;			//NW
	int *idw2;			//NW
	    
	double *normalX, *normalY;	//NW
	double *deltaX;		//NW
	double *length;		//NW
	double *distNormal; //NW
    double *distCentX, *distCentY;
    
    double *nman2wall; 		//NW
    double *gp;			//NW

	int *typeOfBound;	//NW
    int *solidWall;		//NW

    double *qnormalL;	//NW
	double *localDt;	//NW

};


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

    //boundaries

};


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
