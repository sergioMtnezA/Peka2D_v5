#include "loadData.h"

////////////////////////////////////////////////////////////////
EXPORT_DLL int loadControlParameters(
    Peka2D_Setup *pksetup, 
    t_parameters *spar, 
    t_mesh *mesh, 
    t_message *msg){
/*----------------------------*/

    char filename[1024], temp[1024];

	//Read run control file
    sprintf(filename,"%s%s.DAT",spar->dir,spar->proj);
	if(readControlDataFile(filename, pksetup, msg)){
		sprintf(temp,"Reading control file completed");
		Notify(temp,MSG_L1,msg);		
	}

	if(setControlParameters(pksetup, spar, mesh, msg)){
		sprintf(temp,"Set control parameters completed");
		Notify(temp,MSG_L1,msg);		
	}    

	return 1;
	
}

////////////////////////////////////////////////////////////////
EXPORT_DLL int readControlDataFile(
    char *filename,
	Peka2D_Setup *pksetup, 
	t_message *msg){
/*----------------------------*/

	FILE *fp;
	int nc;
	double aux;
	char temp[1024];
	Peka2D_Run *run;
	
	//Load .DAT file
	fp = fopen(filename,"r");
	if(!fp){
		sprintf(temp,"%s file not found",filename);
		Notify(temp,MSG_ERROR,msg);
		return 0;
	}

	//Assign local pointer to pksetup structure
	run = &(pksetup->pkrun);

	////////////////////////////////////////////////////////////////////////////////////////
	// Line 1: Internal program version number
	fscanf(fp,"%d",&run->release);

	////////////////////////////////////////////////////////////////////////////////////////
	// Line 2: Model selector switch
	fscanf(fp,"%d",&run->itrash);

	////////////////////////////////////////////////////////////////////////////////////////
	// Line 3: Physical processes of components switches. The number depends on the release version.
	switch(run->release){
		case 202407:
			fscanf(fp,"%d %d %d %d %d %d %d %d %d %d %d",
				&run->rain,
				&run->sediment,
				&run->itrash,
				&run->weirs,
				&run->itrash,
				&run->sources,
				&run->itrash,
				&run->itrash,
				&run->itrash,
				&run->dambreach,
				&run->itrash);
			break;
		default:
			sprintf(temp,"Release %d not valid",run->release);
			Notify(temp,MSG_ERROR,msg);
			return 0;
			break;
	}

	////////////////////////////////////////////////////////////////////////////////////////
	// Line 4: Wet-dry
	fscanf(fp,"%d",&run->itrash);

	////////////////////////////////////////////////////////////////////////////////////////
	// Line 5: Output control switches
	fscanf(fp,"%d %d %d %d %d",
		&run->extremes, 
		&run->crossSection, 
		&run->profile, 
		&run->itrash, 
		&run->obs);

	////////////////////////////////////////////////////////////////////////////////////////
	// Line 6: Time control data
	fscanf(fp,"%d %lf %lf %lf %lf",
		&run->nIterOut, 
		&run->CFL, 
		&run->dtDump, 
		&run->dtOut, 
		&run->tLimit);

	run->dtDump*=3600.0;
	run->dtOut*=3600.0;
	run->tLimit*=3600.0;

	////////////////////////////////////////////////////////////////////////////////////////
	// Line 7: Initial conditions and hot start control switches
	fscanf(fp,"%d %d",&run->iniType, &run->hotStart);

	////////////////////////////////////////////////////////////////////////////////////////
	// Line 8: Manning's n variable with depth switch
	fscanf(fp,"%lf",&run->dtrash);

	////////////////////////////////////////////////////////////////////////////////////////
	// Line 9: Manning's n values global multiplication factor
	fscanf(fp,"%lf",&run->XnMan);

	////////////////////////////////////////////////////////////////////////////////////////
	// Line 10: Mass balance reporting switch
	fscanf(fp,"%d",&run->writeMass);

	////////////////////////////////////////////////////////////////////////////////////////
	// Line 11: Unit system definition switches
	fscanf(fp,"%d",&run->itrash);

	////////////////////////////////////////////////////////////////////////////////////////
	// Line 12: MinDepth
	fscanf(fp,"%lf",&run->hMin);

	////////////////////////////////////////////////////////////////////////////////////////
	// Line 13: Initial water surface elevation
	fscanf(fp,"%lf",&run->initialWSE);

	////////////////////////////////////////////////////////////////////////////////////////
	// Line 14: Pollutant transport model switch
	fscanf(fp,"%d",&run->solutes);

	////////////////////////////////////////////////////////////////////////////////////////
	// Line 15: Wind stress switch
	fscanf(fp,"%d",&run->wind);

	////////////////////////////////////////////////////////////////////////////////////////
	// Line 16: Mudflow model switch
	fscanf(fp,"%d",&run->itrash);

	////////////////////////////////////////////////////////////////////////////////////////
	// Line 17: Number of cores or GPU ID
	fscanf(fp,"%d",&nc);
	#if SET_SIMGPU
		run->gdevice=nc;
		run->ncores=1;
		sprintf(temp,"GPU device ID set by user: %d",run->gdevice);
	#else
		run->gdevice=0;
		#if SET_OPENMP
			run->ncores=nc;
			if(run->ncores>omp_get_num_procs() || run->ncores<1){
				run->ncores=omp_get_num_procs();
			}
		#else
			run->ncores=1;
		#endif
		sprintf(temp,"Number of CPU cores set by user: %d",run->ncores);
	#endif
	Notify(temp,MSG_WARN,msg);

	////////////////////////////////////////////////////////////////////////////////////////
	// Line 18: time*.EXP files writing
	fscanf(fp,"%d",&run->itrash);

	////////////////////////////////////////////////////////////////////////////////////////
	// Line 19: Additional components
    //none

    fclose(fp);

	return 1;
	
}




////////////////////////////////////////////////////////////////
EXPORT_DLL int setControlParameters(
    Peka2D_Setup *pksetup, 
    t_parameters *spar, 
    t_mesh *mesh, 
    t_message *msg){
/*----------------------------*/

    Peka2D_Run *run;

    //Assign local pointer to pksetup structure
	run = &(pksetup->pkrun);

    //run parameters
	spar->ncores=run->ncores;
    spar->gpuid=run->gdevice;

	spar->ti=0.0;
    spar->tf=run->tLimit;
    spar->CFL=run->CFL;

	spar->writeMass=run->writeMass;
	spar->writeExtremes=run->extremes;
	spar->indexWriteHotstart=run->hotStart;
	
    spar->nIterOut=run->nIterOut;
    spar->dtOut=run->dtOut;
	spar->dtDump=run->dtDump;

	spar->minh=run->hMin;

    //initialize configuration parameters
    mesh->nSolutes=0;

    return 1;

}

////////////////////////////////////////////////////////////////
EXPORT_DLL int loadMeshData(
    Peka2D_Setup *pksetup, 
    t_parameters *spar, 
    t_mesh *mesh, 
    t_message *msg){
/*----------------------------*/

	char meshfile[1024];
	int Nnode = 0;
	int Ncell = 0;

	t_g_cell *myg_Cell;
	t_c_cell *myc_Cell;
	t_node *myNode;
	l_g_cells *lgcells;
	l_c_cells *lccells;
	l_nodes *lnodes;
	
    l_wall *w_calc;
	l_wall *w_bound;

	FILE *fp,*fp1;
	char temp[1024];
    char filename[1024];
    int i,index,nsolutes,nseds,nSolutesMax;
    int flag;

    // Read mesh size
    sprintf(meshfile,"%s%s.FED",spar->dir,spar->proj);
    fp = fopen(meshfile,"r");
	if(!fp){
		sprintf(temp,"%s file not found",meshfile);
		Notify(temp,MSG_ERROR,msg);
		return 0;
	}
    fscanf(fp,"%d %d %*d %*d",&Ncell,&Nnode);
    fclose(fp);


    /*****************************************************************
    // Create the cell structure
    *****************************************************************/
    mesh->ncells = Ncell;  

	lgcells = (l_g_cells*) malloc( sizeof(l_g_cells) );
	//lgcells->cells = (t_g_cell*) malloc( Ncell * sizeof(t_g_cell) );
	new_lg_cells(lgcells,Ncell);

	lccells = (l_c_cells*) malloc( sizeof(l_c_cells));
	//lccells->cells = (t_c_cell*) malloc( Ncell * sizeof(t_c_cell));
	new_lc_cells(lccells,Ncell);

	// Create/initialize new cells and its correspound!
	for (i=0; i<Ncell; i++) {
		myc_Cell=(t_c_cell*)malloc(sizeof(t_c_cell));
		myg_Cell=(t_g_cell*)malloc(sizeof(t_g_cell));
		new_g_cell(myg_Cell,i);
		new_c_cell(myc_Cell,i);
		myc_Cell->geom = myg_Cell;
		add_g_cell(myg_Cell,lgcells);
		add_c_cell(myc_Cell,lccells);
	}

    mesh->g_cells = lgcells;
	mesh->c_cells = lccells;
    

    /*****************************************************************
    // Create the node structure
    *****************************************************************/
    mesh->nnodes = Nnode;

	// Allocate nodes list
	lnodes=(l_nodes*)malloc(sizeof(l_nodes));
	//lnodes->nodes=(t_node*)malloc(Nnode*sizeof(t_node));
	new_lnodes(lnodes,Nnode);

	// Create/initialize new nodes
	for(i=0; i<Nnode; i++){
		myNode=(t_node*)malloc(sizeof(t_node));
		new_node(myNode,i);
		add_node(myNode,lnodes);
	}

    mesh->nodes = lnodes;
	
    /*****************************************************************
    // Generate wall structure
    *****************************************************************/
	w_calc=(l_wall*)malloc(sizeof(l_wall));
	//w_calc->wall=(t_wall*)malloc( 2 * mesh->ncells *sizeof(t_wall));
    new_lwalls(w_calc, 2*mesh->ncells);
	mesh->w_calc = w_calc;

	w_bound=(l_wall*)malloc(sizeof(l_wall));
	//w_bound->wall=(t_wall*)malloc(mesh->ncells*sizeof(t_wall));
    new_lwalls(w_bound, mesh->ncells);
	mesh->w_bound = w_bound;


    /*****************************************************************
    // Load nodes and cells
    *****************************************************************/
    if(!IsMeshAllocated(mesh,msg)){
        sprintf(temp,"Mesh memory allocation failed");
		Notify(temp,MSG_ERROR,msg);
        return 0;
    }

    if(!readMeshFile(meshfile,pksetup,mesh,msg)){
        sprintf(temp,"Mesh file reading failed");
		Notify(temp,MSG_ERROR,msg);
        return 0;
    }


    /*****************************************************************
    // Calculate normal vectors of the cells
    *****************************************************************/
    if(!calc_norm(lgcells,lccells,mesh->NCwall,msg)){
        sprintf(temp,"Mesh geometry calculation failed");
		Notify(temp,MSG_ERROR,msg);
        return 0;
    }    

    /*****************************************************************
    // Create the neighbouring cells structure
    *****************************************************************/
	// The pocells list!
	if(!cons_lpocells_qsort(mesh,msg)){
        sprintf(temp,"Mesh conectivity generation failed");
		Notify(temp,MSG_ERROR,msg);
		return 0;
	}

    /*****************************************************************
    // Check mesh quality
    *****************************************************************/
	if(!check_angle(lgcells,msg)){
        sprintf(temp,"Mesh quality checking failed");
		Notify(temp,MSG_ERROR,msg);
		return 0;
 	}

    /*****************************************************************
    // Set initial cell conditions
    *****************************************************************/
	if(!setInitialState(pksetup,spar,mesh,msg)){
        sprintf(temp,"Set initial state failed");
		Notify(temp,MSG_ERROR,msg);
		return 0;
 	}    

    return 1;
}


////////////////////////////////////////////////////////////////
int IsMeshAllocated(
    t_mesh *mesh, 
    t_message *e){
/*----------------------------*/
    char temp[1024];

    // Check if main array are allocated
    if(mesh == NULL){
        Notify("Mesh has not been allocated due to memory space",MSG_ERROR,e);
        return(0);
    }
    // Check if primary arrays are allocated
    if(mesh->nodes == NULL){
        Notify("Nodes list has not been allocated due to memory space",MSG_ERROR,e);
        return(0);
    }
    if(mesh->g_cells == NULL){
        Notify("Geometric cell list has not been allocated due to memory space",MSG_ERROR,e);
        return(0);
    }
    if(mesh->c_cells == NULL){
        Notify("Computational cell list has not been allocated due to memory space",MSG_ERROR,e);
        return(0);
    }
    if(mesh->w_calc== NULL){
        Notify("Wall list has not been allocated due to memory space",MSG_ERROR,e);
        return(0);
    }
    if(mesh->w_bound== NULL){
        Notify("Wall-Bound list has not been allocated due to memory space",MSG_ERROR,e);
        return(0);
    }

    // Check if secondary arrays are allocated
    if(mesh->nodes->nodes == NULL){
        Notify("Nodes list has not been allocated due to memory space",MSG_ERROR,e);
        return(0);
    }
    if(mesh->g_cells->cells == NULL){
        Notify("Geometric cell list has not been allocated due to memory space",MSG_ERROR,e);
        return(0);
    }
    if(mesh->c_cells->cells== NULL){
        Notify("Computational cell list has not been allocated due to memory space",MSG_ERROR,e);
        return(0);
    }
    if(mesh->w_calc-> wall== NULL){
        Notify("Wall list has not been allocated due to memory space",MSG_ERROR,e);
        return(0);
    }
    if(mesh->w_bound->wall== NULL){
        Notify("Wall-Bound list has not been allocated due to memory space",MSG_ERROR,e);
        return(0);
    }

    sprintf(temp,"Mesh memory allocation completed");
    Notify(temp,MSG_L1,e);
    return(1);
}


////////////////////////////////////////////////////////////////
int readMeshFile(
    char *filename,
    Peka2D_Setup *pksetup,
    t_mesh *mesh, 
    t_message *e){
/*----------------------------*/

    FILE *fp;
    int i;
    int Nnode,Ncell,NCwall;
    int n1,n2,n3,n4;

    double aux1, wse;
    char temp[1024];
	Peka2D_Run *run;
	
	//Assign local pointer to pksetup structure
	run = &(pksetup->pkrun);    

    //Open FED file
    fp = fopen(filename,"r");

    // Initialize counters
    mesh->nodes->n=0;
    mesh->g_cells->n=0;
    mesh->c_cells->n=0;

    // Read the header
    fscanf(fp,"%d %d %d %*d",&Ncell,&Nnode,&NCwall);
    if(NCwall!=3 && NCwall!=4){
        sprintf(temp,"Number of wall-per-cell not allowed: %d",NCwall);
        Notify(temp,MSG_ERROR,e);
        return(0);       
    }
    mesh->NCwall=NCwall;

    // Loop through the nodes
    for(i=0; i<Nnode; i++){
        // id,x,y, z,initial wse, max erotion depth, bctype, bcfilename
        fscanf(fp,"%d %lf %lf %lf %*f %*f %*d %*s",
            &mesh->nodes->nodes[i].id,
            &mesh->nodes->nodes[i].x, 
            &mesh->nodes->nodes[i].y,
            &mesh->nodes->nodes[i].z);

        mesh->nodes->nodes[i].id--;
        mesh->nodes->n++;
    }
    sprintf(temp,"Nodes read %d",mesh->nodes->n);
    Notify(temp,MSG_L0,e);


    // Loop through the cells
    for(i=0; i<Ncell; i++){
        if(NCwall==3){
            //id, n1, n2, n3, nMan, z, wse, aux1, aux2);
            fscanf(fp,"%d %d %d %d %lf %lf %lf %lf %*lf",
                &mesh->c_cells->cells[i].id,
                &n1,
                &n2,
                &n3,
                &mesh->c_cells->cells[i].nMan,
                &mesh->c_cells->cells[i].z,
                &mesh->c_cells->cells[i].h,
                &aux1);
        }else if(NCwall==4){
            //id, n1, n2, n3, n4, nMan, z, wse, aux1, aux2);
            fscanf(fp,"%d %d %d %d %d %lf %lf %lf %lf %*lf",
                &mesh->c_cells->cells[i].id,
                &n1,
                &n2,
                &n3,
                &n4,
                &mesh->c_cells->cells[i].nMan,
                &mesh->c_cells->cells[i].z,
                &mesh->c_cells->cells[i].h,
                &aux1);
        }            

        //Assing nodes
        mesh->c_cells->cells[i].id--;
        mesh->g_cells->cells[i].id = mesh->c_cells->cells[i].id;
        n1--;
        mesh->g_cells->cells[i].nodes[0]=&(mesh->nodes->nodes[n1]);
        n2--;
        mesh->g_cells->cells[i].nodes[1]=&(mesh->nodes->nodes[n2]);
        n3--;
        mesh->g_cells->cells[i].nodes[2]=&(mesh->nodes->nodes[n3]);
        if(NCwall==4){
            n4--;
            mesh->g_cells->cells[i].nodes[3]=&(mesh->nodes->nodes[n4]);            
        }

        mesh->g_cells->n++;
        mesh->c_cells->n++;
    }
    sprintf(temp,"Cells read %d",mesh->g_cells->n);
    Notify(temp,MSG_L0,e);

    fclose(fp);

    sprintf(temp,"Mesh file reading completed");
    Notify(temp,MSG_L1,e);
    return(1);
}


////////////////////////////////////////////////////////////////
int setInitialState(
    Peka2D_Setup *pksetup,
    t_parameters *spar,
    t_mesh *mesh, 
    t_message *e){
/*----------------------------*/

    FILE *fp;
    int i;
    // Here, c1 is defined as a local pointer
    t_c_cell *c1;

    double aux1, wse;
    char temp[1024];
	Peka2D_Run *run;
	
	//Assign local pointer to pksetup structure
	run = &(pksetup->pkrun);

    mesh->minZ=1e6;
    mesh->maxZ=-1e6;    

	for(i=0;i<mesh->ncells;i++){
        // Here, c1 is pointed to the memory position of the i-th cell contained in mesh.
		c1=&(mesh->c_cells->cells[i]);

        //Initial bed level
        c1->zini=c1->z;

        //Initialize water depth
        wse=c1->h;
        if(run->iniType==0){
            c1->h = run->initialWSE - c1->z;
        }else if(run->iniType==1){
            c1->h=0.0;
        }else if(run->iniType==2){
            c1->h = wse - c1->z;
        } 

        //negative flow depth fix
        if(c1->h<0.0){ 
            c1->h=0.0;
        }     

        //initial wetted area
        if(c1->h>0.0){
            c1->hini=c1->h;
        }else{
            c1->hini=-1.;
        }
        
        //Scale nMan
        c1->nMan*=run->XnMan;

        // Initialize momentum
        c1->hu=0.0;
        c1->hv=0.0;
        c1->u=0.;
        c1->v=0.;
        c1->modulou=0.0;

        //Altitude range
        mesh->minZ = MIN(mesh->minZ,c1->z);
        mesh->maxZ = MAX(mesh->maxZ,c1->z);          

	} // end cell loop


    sprintf(temp,"Set initial state completed");
    Notify(temp,MSG_L1,e);
    return(1);
}


////////////////////////////////////////////////////////////////
int loadBoundaryConditions(
    Peka2D_Setup *pksetup, 
    t_parameters *spar,    
    t_mesh *mesh, 
    t_message *e){
/*----------------------------*/

    int i;
    int j;
    int k;
    
    int ibc,type;
    int n1,n2,n3,n4;

    Peka2D_NodeBoundary *nbc;
    int *listClosedBC;

    t_wall *w;
    Peka2D_OBCWalls *IOBC;
    Peka2D_OBCWalls *OOBC;

    char filename[1024], temp[1024];   

      
    //BOUNDARY NODES
    nbc=(Peka2D_NodeBoundary*) malloc(sizeof(Peka2D_NodeBoundary));
    nbc->type = (int*) malloc(mesh->nnodes*sizeof(int));
    nbc->obcId = (int*) malloc(mesh->nnodes*sizeof(int));

    sprintf(filename,"%s%s.OBCP",spar->dir,spar->proj);
    if(!readOpenBoundaryNodes(filename,mesh->nnodes,nbc,e)){
        sprintf(temp,"Read open boundary nodes failed");
        Notify(temp,MSG_ERROR,e);
        return 0;
    }      
    pksetup->pknode=nbc;

 
    //CLOSED BOUNDARY WALLS
    listClosedBC = (int*) malloc(mesh->w_bound->n*sizeof(int));
	k=0;
    for(i=0; i<mesh->w_bound->n; i++){
        listClosedBC[k] = 0;

        // Store nodes of each wall
        n1 = mesh->w_bound->wall[i].node1;
        n2 = mesh->w_bound->wall[i].node2;

        // If any of the wall nodes are CLOSED BOUNDARY, the wall is SURELY a closed boundary
        if(nbc->type[n1] == HYD_CLOSED || nbc->type[n2] == HYD_CLOSED){
            listClosedBC[k] = i;
            //printf("w_bound index: %d, cell id: %d, nodes: %d %d, type %d\n", i, mesh->w_bound->wall[i].ccells[0]->id,n1,n2,0);
            k++;
        }

        // If one of the nodes is a CLOSED BOUNDARY (type 0), the wall is assumed to be type 0.
        // This is already considered by the earlier conditional, but we will warn the user
        if(nbc->type[n1] != HYD_CLOSED || nbc->type[n2] != HYD_CLOSED){
            if(nbc->type[n1] == HYD_CLOSED){
                sprintf(temp,"Boundary ID %d: Nodes %d-%d are set to closed bound wall",nbc->type[n2],n2,n1);
                Notify(temp,MSG_WARN,e);
            }
            if(nbc->type[n2] == HYD_CLOSED){
                sprintf(temp,"Boundary ID %d: Nodes %d-%d are set to closed bound wall",nbc->type[n1],n1,n2);
                Notify(temp,MSG_WARN,e);
            }
        }

        // If the nodes of the wall have different OPEN types (> 0) there is an inconsistency in the definitions
        if(nbc->type[n1] != HYD_CLOSED && nbc->type[n2] != HYD_CLOSED ){
            if( nbc->obcId[n1] != nbc->obcId[n2]){
                sprintf(temp,"Inconsistent boundary definition in nodes %d-%d: Two consecutive nodes must have the same ID open boundary",n1+1,n2+1);
                Notify(temp,MSG_ERROR,e);
                return(0);
            }
        }
    }
    sprintf(temp,"Found %d/%d closed boundary walls",k,mesh->nw_bound);
    Notify(temp,MSG_L0,e);    

    //OPEN BOUNDARIES WALLS
    IOBC = (Peka2D_OBCWalls*) malloc(sizeof(Peka2D_OBCWalls) * HYD_N_OBT);
    OOBC = (Peka2D_OBCWalls*) malloc(sizeof(Peka2D_OBCWalls) * HYD_N_OBT);

    for(i=0; i < HYD_N_OBT; i++){
        IOBC[i].wall = (int*) malloc(sizeof(int)*mesh->w_bound->n);
        IOBC[i].n  = 0;
        OOBC[i].wall = (int*) malloc(sizeof(int)*mesh->w_bound->n);
        OOBC[i].n  = 0;
        for(j=0; j < mesh->w_bound->n; j++){
            IOBC[i].wall[j] = 0;
            OOBC[i].wall[j] = 0;
        }
    }

    j=0;
    for(i=0; i < mesh->w_bound->n; i++){
        // Store nodes of each wall
        n1 = mesh->w_bound->wall[i].node1;
        n2 = mesh->w_bound->wall[i].node2;

        //If both nodes of the wall share the same boundary type, the wall is a single boundary
        if( (nbc->type[n1] != HYD_CLOSED) && 
            (nbc->type[n2] != HYD_CLOSED) && 
            (nbc->type[n1] == nbc->type[n2]) ){

            if(nbc->type[n1] < 0){
                k = abs(nbc->type[n1])-1;
                j = OOBC[k].n;
                OOBC[k].wall[j] = i;      // Store in position j the boundary wall index i
                OOBC[k].n++;

            }else if(nbc->type[n1] > 0){
                k = abs(nbc->type[n1])-1;
                j = IOBC[k].n;
                IOBC[k].wall[j] = i;      // Store in position j the boundary wall index i
                IOBC[k].n++;
            }

        }
    }
    pksetup->IOBC=IOBC;
    pksetup->OOBC=OOBC;

    //count inlets-outlets from wall list
    k=0;
    mesh->nInlet = 0;
    mesh->nOutlet = 0;
    for(i=0; i<HYD_N_OBT; i++){
        // INLETS
        if(IOBC[i].n > 0){
            sprintf(temp,"Found %d/%d open boundary walls in inlet %d",IOBC[i].n,mesh->nw_bound,i);
            Notify(temp,MSG_L0,e);
            mesh->nInlet++;
        }
        // OUTLETS
        if(OOBC[i].n > 0){
            sprintf(temp,"Found %d/%d open boundary walls in outlet %d",OOBC[i].n,mesh->nw_bound,i);
            Notify(temp,MSG_L0,e);
            mesh->nOutlet++;
        }
    }
    if(nbc->countInlet != mesh->nInlet){
        sprintf(temp,"Number of inlets file OBCP (%d) is inconsistent with number of inlets found in mesh (%d)",nbc->countInlet,mesh->nInlet);
        Notify(temp,MSG_ERROR,e);
        return(0);
    }
    if(nbc->countOutlet != mesh->nOutlet){
        sprintf(temp,"Number of inlets file OBCP (%d) is inconsistent with number of inlets found in mesh (%d)",nbc->countOutlet,mesh->nOutlet);
        Notify(temp,MSG_ERROR,e);
        return(0);
    }            

    if(mesh->nInlet == 0){
        sprintf(temp,"No inlet OBC defined");
        Notify(temp,MSG_WARN,e);
    }else{
        sprintf(temp,"Found %d inlet OBC",mesh->nInlet);
        Notify(temp,MSG_L0,e);
    }

    if(mesh->nOutlet == 0){
        sprintf(temp,"No outlet OBC defined");
        Notify(temp,MSG_WARN,e);
    }else{
        sprintf(temp,"Found %d outlet OBC",mesh->nOutlet);
        Notify(temp,MSG_L0,e);
    }

    //OPEN BOUNDARY CONDITIONS
    if(mesh->nInlet+mesh->nOutlet>0){

        //create open bounds structures
        if(!createOpenBounds(nbc, IOBC, OOBC, spar, mesh, e)){
            sprintf(temp,"Creating open boundary structures failed");
            Notify(temp,MSG_ERROR,e);
            return 0;
        }              
 
        //read open boundary conditions
        sprintf(filename,"%s%s.OBCP",spar->dir,spar->proj);
        if(!readOpenBoundaryFile(filename, nbc, IOBC, OOBC, spar, mesh, e)){
            sprintf(temp,"Read open boundary conditions data failed");
            Notify(temp,MSG_ERROR,e);
            return 0;
        }        
    }  

    free(listClosedBC);

    Notify("Boundaries ready",MSG_L1,e);

    return(1);
}


////////////////////////////////////////////////////////////////
int readOpenBoundaryNodes(
    char *filename,
    int nnodes,
    Peka2D_NodeBoundary *nbc, 
    t_message *e){
/*----------------------------*/

	FILE *f;
	int nb=0;;
    int i,j;
	int npb,idn,obcId,np;
	int type,inlet;
	char temp[1024];
	int countInlet;
	int countOutlet;
	int release=0;

    //initialize boundary nodes
    nbc->n=0;
    for(idn=0;idn<nnodes;idn++){
        nbc->type[idn]=HYD_CLOSED;
        nbc->obcId[idn]=-1;
    }

    //read OBCP file
	f=fopen(filename,"r");
	if(f==NULL){
		sprintf(temp,"File %s with Open Boundary Conditions not found",filename);
		Notify(temp,MSG_ERROR,e);
		return(0);
	}    

    fscanf(f,"%d",&release);
    if(release!=202407){
        sprintf(temp,"The Solute file version should be 202407. Current version: %d",release);
        Notify(temp,MSG_ERROR,e);
        return(0);
    }

    fscanf(f,"%d",&nb);
    sprintf(temp,"Number of open boundaries defined %d",nb);
    Notify(temp,MSG_L0,e);


    countInlet=0;
    countOutlet=0;
    if(nb>0){
        obcId=0;    
        for(i=0;i<nb;i++){ //Boundaries loop
            fscanf(f,"%*s"); 
            fscanf(f,"%d",&type); 
            inlet=0;
            if( type==HYD_INFLOW_Q ||
                type==HYD_INFLOW_HZ ||  
                type==HYD_INFLOW_QHZ){
                countInlet++;
                inlet=1;
            }else if( type==HYD_OUTFLOW_GAUGE ||
                type==HYD_OUTFLOW_HZ ||
                type==HYD_OUTFLOW_FREE ||
                type==HYD_OUTFLOW_FR||
                type==HYD_OUTFLOW_NORMAL){
                countOutlet--;
                inlet=0;
            }else{
                sprintf(temp,"Undefined open boundary condition type %d",type);
                Notify(temp,MSG_ERROR,e);
                fclose(f);
                return(0);
            }

            fscanf(f,"%*s"); // filename
            fscanf(f,"%d",&npb); // Number of nodes
            //printf("inlet %d npb %d\n", inlet, npb);
            nbc->n+=npb;
            for(j=0;j<npb;j++){
                fscanf(f,"%d",&idn);
                if(inlet){
                    nbc->type[idn-1]=countInlet;
                }else{
                    nbc->type[idn-1]=countOutlet;
                }
                nbc->obcId[idn-1]=obcId;
                //printf("node %d obcId %d in-out (%d)-(%d)\n",idn,obcId,countInlet,countOutlet);
            }
            sprintf(temp,"Number of nodes in ID %d open boundary: %d",obcId,npb);
            Notify(temp,MSG_L0,e);

            obcId++;
            
        }

        sprintf(temp,"Read open boundary nodes completed");
        Notify(temp,MSG_L1,e);

    }
    fclose(f);

    nbc->countInlet=countInlet;
    nbc->countOutlet=fabs(countOutlet);

	return 1;
}


////////////////////////////////////////////////////////////////
int createOpenBounds(
    Peka2D_NodeBoundary *nbc,
    Peka2D_OBCWalls *IOBC,
    Peka2D_OBCWalls *OOBC,
    t_parameters *spar, 
    t_mesh *mesh,    
    t_message *e){
/*----------------------------*/

    int i,j,k;
    char temp[1024];

    double nx,ny;
    double module;

    //create inlets structure
    mesh->in=(t_bound*) malloc(mesh->nInlet*sizeof(t_bound));
    for(i=0; i<mesh->nInlet; i++){
        mesh->in[i].wallBound=(t_wall**) malloc(sizeof(t_wall*));
        mesh->in[i].ncellsBound = 0;
        mesh->in[i].ncellsInner = 0;
        mesh->in[i].totalArea = 0.0;
        mesh->in[i].totalLength = 0.0;
        mesh->in[i].zmin = 1E6;
        mesh->in[i].maxFroude=0.9;

        for(j=0; j < IOBC[i].n; j++){
            k = IOBC[i].wall[j];

            if(!build_wall_inlet(mesh, &mesh->in[i], &mesh->w_bound->wall[k], i, e)){
                sprintf(temp,"Build inlet walls failed");
                Notify(temp,MSG_ERROR,e);
                return 0;
            }

        }
        sprintf(temp,"Boundary cells assigned to ID %d inlet: %d",i,mesh->in[i].ncellsBound);
        Notify(temp,MSG_L0,e);          

        if(!build_inner_inlet(mesh, &mesh->in[i], i, e)){
            sprintf(temp,"Build inlet boundary failed");
            Notify(temp,MSG_ERROR,e);
            return 0;
        }
        sprintf(temp,"Inner cells assigned to ID %d inlet: %d",i,mesh->in[i].ncellsInner);
        Notify(temp,MSG_L0,e);  

        //global inlet normal direction
        nx=0.0;
        ny=0.0;
        for(j=0;j<mesh->in[i].ncellsBound;j++){
            nx+=mesh->in[i].wallBound[j]->normal[_X_];
            ny+=mesh->in[i].wallBound[j]->normal[_Y_];
        }
        nx/=mesh->in[i].ncellsBound;
        ny/=mesh->in[i].ncellsBound;
        nx=-nx;
        ny=-ny;

        module=sqrt(nx*nx+ny*ny);
        if(module>TOL12){
            nx*=1./module;
            ny*=1./module;
        }

        mesh->in[i].normal.x=nx;
        mesh->in[i].normal.y=ny;         

    }
    sprintf(temp,"Inlet boundary cell structures completed");
    Notify(temp,MSG_L1,e);       


    // Outlets
    mesh->out=(t_bound*) malloc(mesh->nOutlet*sizeof(t_bound));
    for(i=0; i<mesh->nOutlet; i++){
        mesh->out[i].wallBound=(t_wall**) malloc(sizeof(t_wall*));
        mesh->out[i].ncellsBound = 0;
        mesh->out[i].ncellsInner = 0;
        mesh->out[i].totalArea = 0.0;
        mesh->out[i].totalLength = 0.0;
        mesh->out[i].zmin = 1E6;
        mesh->out[i].maxFroude=0.9;

        //bound cells
        for(j=0; j < OOBC[i].n; j++){
            k = OOBC[i].wall[j];

            if(!build_wall_outlet(mesh, &mesh->out[i], &mesh->w_bound->wall[k], i, e)){
                sprintf(temp,"Build outlet walls failed");
                Notify(temp,MSG_ERROR,e);
            }

        }
        sprintf(temp,"Boundary cells assigned to ID %d outlet: %d",i,mesh->out[i].ncellsBound);
        Notify(temp,MSG_L0,e);   

        //inner cells
        if(!build_inner_outlet(mesh, &mesh->out[i], i, e)){
            sprintf(temp,"Build outlet boundary failed");
            Notify(temp,MSG_ERROR,e);
            return 0;
        } 
        sprintf(temp,"Inner cells assigned to ID %d outlet: %d",i,mesh->out[i].ncellsInner);
        Notify(temp,MSG_L0,e); 

        //global outlet normal direction
        nx=0.0;
        ny=0.0;
        for(j=0;j<mesh->out[i].ncellsBound;j++){
            nx+=mesh->out[i].wallBound[j]->normal[_X_];
            ny+=mesh->out[i].wallBound[j]->normal[_Y_];
        }
        nx/=mesh->out[i].ncellsBound;
        ny/=mesh->out[i].ncellsBound;

        module=sqrt(nx*nx+ny*ny);
        if(module>TOL12){
            nx*=1./module;
            ny*=1./module;
        }

        mesh->out[i].normal.x=nx;
        mesh->out[i].normal.y=ny;

    }
    sprintf(temp,"Outlet boundary cell structures completed");
    Notify(temp,MSG_L1,e); 

    return 1;
}


////////////////////////////////////////////////////////////////
int readOpenBoundaryFile(
    char *filename,
    Peka2D_NodeBoundary *nbc,
    Peka2D_OBCWalls *IOBC,
    Peka2D_OBCWalls *OOBC,
    t_parameters *spar, 
    t_mesh *mesh,    
    t_message *e){
/*----------------------------*/

    int i,j,k;
    int release;
    int nobc;
    int l,ii,jj,nn;
    int type,inlet;
    int countInlet;
    int countOutlet;
    char temp[1024];
    FILE *fp,*fdata;
    char fdataname[1024], label[1024];
    char bcidname[1024];


    //read OBCP file
	fp=fopen(filename,"r");
	if(fp==NULL){
		sprintf(temp,"File %s with Open Boundary Conditions not found",filename);
		Notify(temp,MSG_ERROR,e);
		return(0);
	} 

    fscanf(fp,"%d",&release);
    fscanf(fp,"%d",&nobc);
    //printf("release %d np %d\n", release, np);

    countInlet=0;
    countOutlet=0;
    for(i=0;i<nobc;i++){
        fscanf(fp,"%s",&(bcidname));
        fscanf(fp,"%d",&(type));

        if( type==HYD_INFLOW_Q ||
            type==HYD_INFLOW_HZ || 
            type==HYD_INFLOW_QHZ){
            inlet=1;
            strcpy(mesh->in[countInlet].idname, bcidname);
            mesh->in[countInlet].type=type;

            sprintf(temp,"Inlet boundary idname %s of type %d",mesh->in[countInlet].idname,type);
            Notify(temp,MSG_L0,e); 

        }else if( type==HYD_OUTFLOW_GAUGE ||
            type==HYD_OUTFLOW_HZ ||
            type==HYD_OUTFLOW_FREE ||
            type==HYD_OUTFLOW_FR||
            type==HYD_OUTFLOW_NORMAL){
            inlet=0;
            strcpy(mesh->out[countOutlet].idname,bcidname);
            mesh->out[countOutlet].type=type;

            sprintf(temp,"Outlet boundary idname %s of type %d",mesh->out[countOutlet].idname,type);
            Notify(temp,MSG_L0,e); 

        }else{
            sprintf(temp,"Boundary condition %d not implemented",type);
            Notify(temp,MSG_ERROR,e);
            fclose(fp);
            return(0);
        }

        //boundary table file
        fscanf(fp,"%s",label);

        // saltamos los nodos
        fscanf(fp,"%d",&nn);
        for(j=0;j<nn;j++){
            fscanf(fp,"%*d");
        }

        // Read boundary table file
        sprintf(fdataname,"%s%s",spar->dir,label);
        fdata=fopen(fdataname,"r");
        if(fdata==NULL && type!=HYD_OUTFLOW_FREE){
            sprintf(temp,"File with boundary data %s not found",fdataname);
            Notify(temp,MSG_ERROR,e);
            fclose(fp);
            return(0);
        }

        if(inlet==1){ //Entrada
            switch (mesh->in[countInlet].type){
                case HYD_INFLOW_Q:
                    fscanf(fdata,"%d\n",&(mesh->in[countInlet].n));
                    mesh->in[countInlet].t=(double*) malloc(mesh->in[countInlet].n*sizeof(double));
                    mesh->in[countInlet].q=(double*) malloc(mesh->in[countInlet].n*sizeof(double));

                    for(j=0;j<mesh->in[countInlet].n;j++){
                        fscanf(fdata,"%lf %lf",&(mesh->in[countInlet].t[j]),&(mesh->in[countInlet].q[j]));
                        
                        mesh->in[countInlet].t[j] *= 3600.0;
                    }
                    fclose(fdata);
                    break;

                case HYD_INFLOW_HZ://h+z(t)
                    fscanf(fdata,"%d\n",&(mesh->in[countInlet].n));
                    mesh->in[countInlet].t=(double*) malloc(mesh->in[countInlet].n*sizeof(double));
                    mesh->in[countInlet].hZ=(double*) malloc(mesh->in[countInlet].n*sizeof(double));

                    for(j=0;j<mesh->in[countInlet].n;j++){
                        fscanf(fdata,"%lf %lf",&(mesh->in[countInlet].t[j]),&(mesh->in[countInlet].hZ[j]));
                        
                        mesh->in[countInlet].t[j] *= 3600.0;
                    }
                    fclose(fdata);
                    break;

                case HYD_INFLOW_QHZ:
                    fscanf(fdata,"%d\n",&(mesh->in[countInlet].n));
                    mesh->in[countInlet].t=(double*) malloc(mesh->in[countInlet].n*sizeof(double));
                    mesh->in[countInlet].q=(double*) malloc(mesh->in[countInlet].n*sizeof(double));
                    mesh->in[countInlet].hZ=(double*) malloc(mesh->in[countInlet].n*sizeof(double));

                    for(j=0;j<mesh->in[countInlet].n;j++){
                        fscanf(fdata,"%lf %lf %lf",&(mesh->in[countInlet].t[j]),&(mesh->in[countInlet].q[j]),&(mesh->in[countInlet].hZ[j]));
                        
                        mesh->in[countInlet].t[j] *= 3600.0;
                    }
                    fclose(fdata);
                    break;

                default:
                    sprintf(temp,"Invalid boundary condition (%d) for inlet %d",mesh->in[countInlet].type,countInlet+1);
                    Notify(temp,MSG_ERROR,e);
                    return(0);
                    break;

            }// End case
            mesh->in[countInlet].iAct=0;
            
            sprintf(temp,"Inlet %d boundary conditions ready",countInlet);
            Notify(temp,MSG_L1,e);

            countInlet++;

        }else{ // Salida
            switch (mesh->out[countOutlet].type){
                case HYD_OUTFLOW_GAUGE: //q(h+z)
                    fscanf(fdata,"%d",&(mesh->out[countOutlet].n));
                    mesh->out[countOutlet].hZ=(double*) malloc(mesh->out[countOutlet].n*sizeof(double));
                    mesh->out[countOutlet].q=(double*) malloc(mesh->out[countOutlet].n*sizeof(double));

                    for(j=0;j<mesh->out[countOutlet].n;j++){
                        fscanf(fdata,"%lf %lf",&(mesh->out[countOutlet].hZ[j]),&(mesh->out[countOutlet].q[j]));
                    }
                    fclose(fdata);
                    break;

            case HYD_OUTFLOW_HZ: //h+z(t)
                fscanf(fdata,"%d",&(mesh->out[countOutlet].n));
                mesh->out[countOutlet].hZ=(double*) malloc(mesh->out[countOutlet].n*sizeof(double));
                mesh->out[countOutlet].t=(double*) malloc(mesh->out[countOutlet].n*sizeof(double));

                for(j=0;j<mesh->out[countOutlet].n;j++){
                    fscanf(fdata,"%lf %lf",&(mesh->out[countOutlet].t[j]),&(mesh->out[countOutlet].hZ[j]));
                    
                    mesh->out[countOutlet].t[j] *= 3600.0;
                }
                fclose(fdata);
                break;

            case HYD_OUTFLOW_FREE: //salida libre
                break;

            case HYD_OUTFLOW_FR: //Froude cte
                mesh->out[countOutlet].n=1;
                mesh->out[countOutlet].Fr=(double*) malloc(sizeof(double));
                fscanf(fdata,"%lf",&(mesh->out[countOutlet].Fr[0]));
                fclose(fdata);
                break;

            case HYD_OUTFLOW_NORMAL:
                fscanf(fdata,"%d",&(mesh->out[countOutlet].n));
                mesh->out[countOutlet].hZ=(double*) malloc(mesh->out[countOutlet].n*sizeof(double));
                mesh->out[countOutlet].q=(double*) malloc(mesh->out[countOutlet].n*sizeof(double));

                for(j=0;j<mesh->out[countOutlet].n;j++){
                    fscanf(fdata,"%lf %lf",&(mesh->out[countOutlet].hZ[j]),&(mesh->out[countOutlet].q[j]));
                }
                fclose(fdata);
                break;

            default:
                sprintf(temp,"Invalid boundary condition (%d) for outlet %d",mesh->out[i].type,i+1);
                Notify(temp,MSG_ERROR,e);
                return(0);
                break;
            } // end case
            mesh->out[countOutlet].iAct=0;

            sprintf(temp,"Outlet %d boundary conditions ready",countOutlet);
            Notify(temp,MSG_L1,e);

            countOutlet++;
        } // End if inlet==1

    }// End for(i=0;i<nobc;i++){

    sprintf(temp,"Open boundary conditions data ready");
    Notify(temp,MSG_L1,e);

    return 1;
}



#if SET_SOLUTE
////////////////////////////////////////////////////////////////
EXPORT_DLL int loadSoluteData(
    Peka2D_Setup *pksetup, 
    t_parameters *spar, 
    t_mesh *mesh, 
    t_message *msg){
/*----------------------------*/

    char filename[1024], temp[1024];
    int solute_enabled_by_run=0;

    Peka2D_SoluteGroup *soluteGroup;
    soluteGroup = (Peka2D_SoluteGroup*) malloc(sizeof(Peka2D_SoluteGroup));

    //initialize default
    soluteGroup->nSolutes=0;
    sprintf(soluteGroup->initialFile,"%s%s.SOLINITIAL",spar->dir,spar->proj);

    //Load solutes if activated
    solute_enabled_by_run = pksetup->pkrun.solutes;
    if(solute_enabled_by_run){

        //Read solute file
        sprintf(filename,"%s%s.SOLUTES",spar->dir,spar->proj);
        if(readSoluteFile(filename, soluteGroup, msg)){
            sprintf(temp,"Reading solute file completed");
            Notify(temp,MSG_L1,msg);		
        }

        //associted soluteGroup to pksetup
        pksetup->soluteGroup = soluteGroup; 

        //Create solute structures
        if(createSoluteStructures(soluteGroup, spar, mesh, msg)){
            sprintf(temp,"Solute structures completed");
            Notify(temp,MSG_L1,msg);		
        }

        //Initialize solute concentration
        if(setInitialSoluteState(soluteGroup, spar, mesh, msg)){
            sprintf(temp,"Set initial solute state completed");
            Notify(temp,MSG_L1,msg);		
        }  

    }

	return 1;
	
}


////////////////////////////////////////////////////////////////
int readSoluteFile(
    char *filename,
    Peka2D_SoluteGroup *soluteGroup,   
    t_message *e){
/*----------------------------*/

    int i;
    FILE *fp;
    char temp[1024];
    int nSolutes=0;
    int release=0;


    fp = fopen(filename,"r");
    if(!fp){
        sprintf(temp,"%s file not found",filename);
        Notify(temp,MSG_ERROR,e);
        return 0;
    }

    fscanf(fp,"%d",&release); //line 1
    if(release!=202407){
        sprintf(temp,"The Solute file version should be 202407. Current version: %d",release);
        Notify(temp,MSG_ERROR,e);
        return(0);
    }    

    fscanf(fp,"%d",&nSolutes); //line 2

    if(nSolutes>0){

        soluteGroup->nSolutes=nSolutes;
        soluteGroup->solute = (Peka2D_Solute*) malloc(nSolutes*sizeof(Peka2D_Solute));

        fscanf(fp,"%d",&soluteGroup->flagDiffussion); //line 3 

        for(i=0;i<nSolutes;i++){ //line 4
            fscanf(fp,"%d",&soluteGroup->solute[i].typeDiff);
        }

        for(i=0;i<nSolutes;i++){ //line 5 to +nSolutes
            fscanf(fp,"%lf %lf",
                &soluteGroup->solute[i].k_xx,
                &soluteGroup->solute[i].k_yy);
        }

        for(i=0;i<nSolutes;i++){ //line 6 to +nSolutes
            fscanf(fp,"%s",&soluteGroup->solute[i].name);
            //printf("Cname %s\n",soluteGroup->solute[i].name);
        }

        for(i=0;i<nSolutes;i++){ //line 7 to +nSolutes - skip
            for(i=0;i<nSolutes;i++){
                fscanf(fp,"%*f");
            }
        }

        for(i=0;i<nSolutes;i++){ //line 8
            fscanf(fp,"%lf",&soluteGroup->solute[i].maxConc);
            //printf("Cmax %lf\n",soluteGroup->solute[i].maxConc);
        }

    }else{
        sprintf(temp,"Solutes are not defined in file %s - Runing simulation without solutes",filename);
        Notify(temp,MSG_WARN,e);
        return(0);
    }

	fclose(fp);
    
    return(1);

}





////////////////////////////////////////////////////////////////
int createSoluteStructures(
    Peka2D_SoluteGroup *soluteGroup, 
    t_parameters *spar, 
    t_mesh *mesh,    
    t_message *e){
/*----------------------------*/

    int i,j;
    int nSolutes = soluteGroup->nSolutes;
    t_c_cell *c1;

    //solute flag in mesh
    mesh->nSolutes = nSolutes;

    if(mesh->nSolutes){
        
        //solute data in mesh
        mesh->solutes = (l_solutes*) malloc(sizeof(l_solutes));

        mesh->solutes->n = nSolutes;
        mesh->solutes->flagDiffussion = soluteGroup->flagDiffussion;
        //printf("flagDiffusion %d\n",mesh->solutes->flagDiffussion);

        mesh->solutes->solute=(t_solute*) malloc(nSolutes*sizeof(t_solute));
        for(j=0;j<nSolutes;j++){
            sprintf(mesh->solutes->solute[j].name,"%s",soluteGroup->solute[j].name);
            //printf("Cname %s\n",mesh->solutes->solute[j].name);

            mesh->solutes->solute[j].typeDiff=soluteGroup->solute[j].typeDiff;

            mesh->solutes->solute[j].k_xx = soluteGroup->solute[j].k_xx;
            mesh->solutes->solute[j].k_yy = soluteGroup->solute[j].k_yy;

            mesh->solutes->solute[j].maxConc=soluteGroup->solute[j].maxConc;
            //printf("Cmax %lf\n",mesh->solutes->solute[j].maxConc);
            
        }

        //solute variables in c_cells
        for(i=0;i<mesh->ncells;i++){
            c1=&(mesh->c_cells->cells[i]);
            c1->hphi=(double*) malloc(nSolutes*sizeof(double));
            c1->phi=(double*) malloc(nSolutes*sizeof(double));
        }
    
    }

	return 1;
	
}




////////////////////////////////////////////////////////////////
int setInitialSoluteState(
    Peka2D_SoluteGroup *soluteGroup,
    t_parameters *spar,
    t_mesh *mesh, 
    t_message *e){
/*----------------------------*/

    FILE *fp;
    int i,j;

    // Here, c1 is defined as a local pointer
    t_c_cell *c1;

    double aux1, dataSolute;
    char temp[1024];

    int nSolutes = soluteGroup->nSolutes;

    fp=NULL;
    fp=fopen(soluteGroup->initialFile,"r");
    if(fp){

        for(i=0;i<mesh->ncells;i++){
            c1=&(mesh->c_cells->cells[i]);

            for(j=0;j<nSolutes;j++){  
                fscanf(fp,"%lf",&dataSolute);         

                if(c1->h > TOL12){
                    c1->phi[j] = MAX(0.0,dataSolute);
                    c1->hphi[j] = c1->h * c1->phi[j];
                }else{
                    c1->phi[j] = 0.0;
                    c1->hphi[j] = 0.0;
                }
            
            } 

        } // end cell loop
        fclose(fp);

        sprintf(temp,"Solute initial concentration set from file %s",soluteGroup->initialFile);
        Notify(temp,MSG_L1,e);
    }else{

        for(i=0;i<mesh->ncells;i++){
            c1=&(mesh->c_cells->cells[i]);

            for(j=0;j<nSolutes;j++){  
                c1->phi[j] = 0.0;
                c1->hphi[j] = 0.0;
            } 

        } // end cell loop        
        sprintf(temp,"File %s not found - Initial concentration set to 0.0",soluteGroup->initialFile);
        Notify(temp,MSG_WARN,e);
    }        

    return 1;

}
#endif