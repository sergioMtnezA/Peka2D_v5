#include "initialize.h"

////////////////////////////////////////////////////////////////
EXPORT_DLL int allocateArraysMemory(
    t_parameters spar, 
    t_mesh *mesh,
    t_arrays *carrays,
    t_message *msg){
/*----------------------------*/

	int i,j;
    int NCwall;
	int ncells,nwc,nwb;
    int nWallCell;

	size_t free_mem, total_mem;

	//Local variables just for allocation
    NCwall=mesh->NCwall;	//walls per cell
	ncells=mesh->ncells;
    nWallCell=NCwall*ncells;
	nwc=mesh->nw_calc;
	nwb=mesh->nw_bound;


    // COMPUTATION ARRAYS ////////////////////////////////////////////
    //cells
    carrays->activeC=(int*)malloc(ncells*sizeof(int));
    carrays->actCells=(int*)malloc(ncells*sizeof(int));

    carrays->cidx=(int*)malloc(ncells*sizeof(int));
    carrays->nneig=(int*)malloc(ncells*sizeof(int));

    carrays->z=(double*)malloc(ncells*sizeof(double));
	carrays->h=(double*)malloc(ncells*sizeof(double));
	carrays->hu=(double*)malloc(ncells*sizeof(double));
	carrays->hv=(double*)malloc(ncells*sizeof(double));
	carrays->u=(double*)malloc(ncells*sizeof(double));
	carrays->v=(double*)malloc(ncells*sizeof(double));
    carrays->modulou=(double*)malloc(ncells*sizeof(double));
    carrays->sqrh=(double*)malloc(ncells*sizeof(double));
	
	carrays->area=(double*)malloc(ncells*sizeof(double));
	carrays->nman=(double*)malloc(ncells*sizeof(double));
	carrays->SOX=(double*)malloc(ncells*sizeof(double));
	carrays->SOY=(double*)malloc(ncells*sizeof(double));  

    carrays->mass=(double*)malloc(ncells*sizeof(double)); 


    //cells*NCwall
	carrays->dh=(double*)malloc(nWallCell*sizeof(double));
	carrays->dhu=(double*)malloc(nWallCell*sizeof(double));
	carrays->dhv=(double*)malloc(nWallCell*sizeof(double));

	carrays->solidWallByCell=(int*)malloc(nWallCell*sizeof(int));
	carrays->neighCell=(int*)malloc(nWallCell*sizeof(int));
    carrays->neighWall=(int*)malloc(nWallCell*sizeof(int));

    carrays->normalXbyCell=(double*)malloc(nWallCell*sizeof(double));
    carrays->normalYbyCell=(double*)malloc(nWallCell*sizeof(double));
	
    //internal walls
    carrays->activeW=(int*)malloc(nwc*sizeof(int));
    carrays->actWalls=(int*)malloc(nwc*sizeof(int));

    carrays->widx=(int*)malloc(nwc*sizeof(int));

	carrays->idx1=(int*)malloc(nwc*sizeof(int));
	carrays->idx2=(int*)malloc(nwc*sizeof(int));
	carrays->idw1=(int*)malloc(nwc*sizeof(int));
	carrays->idw2=(int*)malloc(nwc*sizeof(int));  

	carrays->normalX=(double*)malloc(nwc*sizeof(double));
	carrays->normalY=(double*)malloc(nwc*sizeof(double));
	carrays->deltaX=(double*)malloc(nwc*sizeof(double));
	carrays->length=(double*)malloc(nwc*sizeof(double));
    carrays->distNormal=(double*)malloc(nwc*sizeof(double));
    carrays->distCentX=(double*)malloc(nwc*sizeof(double));
    carrays->distCentY=(double*)malloc(nwc*sizeof(double));

	carrays->nman2wall=(double*)malloc(nwc*sizeof(double));
	carrays->gp=(double*)malloc(nwc*sizeof(double));

    carrays->typeOfBound=(int*)malloc(nwc*sizeof(int)); 
	carrays->solidWall=(int*)malloc(nwc*sizeof(int));
 
    carrays->qnormalL=(double*)malloc(nwc*sizeof(double));
    carrays->localDt=(double*)malloc(nwc*sizeof(double));


    //boundaries




    return 1;

}




////////////////////////////////////////////////////////////////
EXPORT_DLL int allocateBoundaryArraysMemory(
    t_parameters spar, 
    t_mesh *mesh,
    t_arrays *carrays,
    t_message *msg){
/*----------------------------*/

	int i,j;
    int NCwall;
    int nInlet, nOutlet;
    int nTotalCellsIn, nTotalCellsOut;
	int nwb;


	size_t free_mem, total_mem;

	//Local variables just for allocation
    NCwall=mesh->NCwall;	//walls per cell
    nInlet=mesh->nInlet;
    nOutlet=mesh->nOutlet;
    nTotalCellsIn=mesh->nTotalCellsIn;
    nTotalCellsOut=mesh->nTotalCellsOut;
	nwb=mesh->nw_bound;


    //bounds
    carrays->normalXIn=(double*)malloc(nInlet*sizeof(double));
    carrays->normalYIn=(double*)malloc(nInlet*sizeof(double));

    carrays->normalXOut=(double*)malloc(nOutlet*sizeof(double));
    carrays->normalYOut=(double*)malloc(nOutlet*sizeof(double));


    //cells
    carrays->cidxIn=(int*)malloc(nTotalCellsIn*sizeof(int));
    carrays->idBoundIn=(int*)malloc(nTotalCellsIn*sizeof(int));

    carrays->cidxOut=(int*)malloc(nTotalCellsOut*sizeof(int));
    carrays->idBoundOut=(int*)malloc(nTotalCellsOut*sizeof(int));    
	
    //bound walls
    carrays->lengthwallIn=(double*)malloc(nTotalCellsIn*sizeof(double));

    carrays->lengthwallOut=(double*)malloc(nTotalCellsOut*sizeof(double));

    return 1;

}




////////////////////////////////////////////////////////////////
EXPORT_DLL int initilizeControlArrays(
    t_parameters spar, 
    t_mesh *mesh,
    t_arrays *carrays,
    t_message *msg){
/*----------------------------*/

    size_t free_mem, total_mem;

    //////////////////////////////////COMPUTATION ARRAYS
  	carrays->ncores=spar.ncores;
    carrays->gpuid=spar.gpuid;
	carrays->ti=spar.ti;
    carrays->tf=spar.tf;
    carrays->CFL=spar.CFL;

	carrays->writeMass=spar.writeMass;
	carrays->writeExtremes=spar.writeExtremes;
	carrays->indexWriteHotstart=spar.indexWriteHotstart;

    carrays->nIterOut=spar.nIterOut;
	carrays->dtOut=spar.dtOut;
	carrays->dtDump=spar.dtDump;
	carrays->minh=spar.minh;

    //mesh metadata
	carrays->NCwall=mesh->NCwall;
	carrays->ncells=mesh->ncells;
    carrays->nWallCell=mesh->NCwall*mesh->ncells;

	carrays->nw_calc=mesh->nw_calc;
	carrays->nw_bound=mesh->nw_bound;

	carrays->nInlet=mesh->nInlet;
	carrays->nOutlet=mesh->nOutlet;

    //computation controls
	carrays->t = spar.ti;
    carrays->dt = 0.0;

    carrays->nActCells=0;
    carrays->nActWalls=0;

	carrays->nIter=1;
    carrays->indexOut=1;
    carrays->indexDump=1;
    carrays->dumpComponent=0;
    carrays->dumpState=0;    

	carrays->massOld=0.0;
	carrays->massNew=0.0;
    carrays->massError=0.0; 

    //boundary control
    carrays->massIn=0.0;
    carrays->massOut=0.0;
    carrays->qTotalIn=0.0;
    carrays->qTotalOut=0.0; 


    //solute permanent flag to switch on/off memory & computation
	carrays->nSolutes = mesh->nSolutes;


    return 1;
 
}




// ////////////////////////////////////////////////////////////////
// EXPORT_DLL int initilizeCommunicationArrays(
//     t_parameters spar, 
//     t_mesh *mesh,
//     t_arrays *carrays,  
//     t_cuPtr *cuPtr, 
//     t_message *msg){
// /*----------------------------*/

// 	int i,j;
//     int cidx;
// 	int NCwall,ncells,nwc,nwb,nWallCell;
//     int nActCells, nActWalls;

//     double daux;

//     t_c_cell *c1,*c2;
//     t_g_cell *g1;
//     t_wall *w1;
	
// 	//Local variables just for the function
//     NCwall=mesh->NCwall;
// 	ncells=mesh->ncells;
//     nwc=mesh->nw_calc;
//     nwb=mesh->nw_bound;

//     //////////////////////////////////COMMUNICATION ARRAYS
//     cuPtr->t=carrays->ti;
//     cuPtr->dt=carrays->dt;

// 	cuPtr->nIter=carrays->nIter;
//     cuPtr->indexOut=carrays->indexOut;
//     cuPtr->indexDump=carrays->indexDump;
//     cuPtr->dumpComponent=carrays->dumpComponent;
//     cuPtr->dumpState=carrays->dumpState; 

//     cuPtr->nActCells=carrays->nActCells;
//     cuPtr->nActWalls=carrays->nActWalls; 

//     cuPtr->massOld=carrays->massOld;
//     cuPtr->massNew=carrays->massNew;
//     cuPtr->massError=carrays->massError; 

//     cuPtr->qTotalIn=carrays->qTotalIn;
//     cuPtr->qTotalOut=carrays->qTotalOut;            
	
//     return 1;

// }




////////////////////////////////////////////////////////////////
EXPORT_DLL int initilizeMeshArrays(
    t_parameters spar, 
    t_mesh *mesh,
    t_arrays *carrays,   
    t_message *msg){
/*----------------------------*/

	int i,j;
	int NCwall,ncells,nwc,nwb,nWallCell;
    int nActCells, nActWalls;

    double daux;

    t_c_cell *c1,*c2;
    t_g_cell *g1;
    t_wall *w1;

    size_t free_mem, total_mem;
	
	//Local variables just for the function
    NCwall=mesh->NCwall;
	ncells=mesh->ncells;
    nWallCell=carrays->nWallCell;
    nwc=mesh->nw_calc;
    nwb=mesh->nw_bound;

	//Copy array values by cells
    nActCells=0;
	for(i=0;i<ncells;i++){
        g1=&(mesh->g_cells->cells[i]);
        c1=&(mesh->c_cells->cells[i]);

        carrays->cidx[i]=g1->id;
        carrays->nneig[i]=g1->nneig;
        //printf("nneig %d \n",carrays->nneig[i]);        

        carrays->z[i]=c1->z;

		carrays->h[i]=c1->h;
        if(c1->h > TOL12){
            carrays->sqrh[i]=sqrt(c1->h);
            if(c1->h > carrays->minh){
                carrays->hu[i]=c1->hu;
                carrays->hv[i]=c1->hv;
                carrays->u[i]=c1->u;
                carrays->v[i]=c1->v;
                carrays->modulou[i]=c1->modulou;
            }else{
                carrays->hu[i]=0.0;
                carrays->hv[i]=0.0;
                carrays->u[i]=0.0;
                carrays->v[i]=0.0;
                carrays->modulou[i]=0.0;             
            }
        }else{
            carrays->h[i]=0.0;
            carrays->sqrh[i]=0.0;
            carrays->hu[i]=0.0;
            carrays->hv[i]=0.0;
            carrays->u[i]=0.0;
            carrays->v[i]=0.0;
            carrays->modulou[i]=0.0;              
        }

		carrays->area[i]=g1->area;
		carrays->nman[i]=c1->nMan;
		//carrays->SOX[i]=mesh->c_cells->cells[i].SOX;
		//carrays->SOY[i]=mesh->c_cells->cells[i].SOY;

        carrays->mass[i]=c1->h*g1->area;		

        // initial active cells -------------------------------------------
        carrays->activeC[i]=0;
        carrays->actCells[i]=-1;
        if(c1->h > TOL12){
            carrays->activeC[i]=1;
            carrays->actCells[nActCells]=g1->id;
            nActCells++;
        }
        // ---------------------------------------------------------------

    }
    carrays->nActCells=nActCells;
	

	for(i=0;i<nWallCell;i++){
		carrays->dh[i]=0.0;
		carrays->dhu[i]=0.0;
		carrays->dhv[i]=0.0;
	}
    for(i=0;i<ncells;i++){
        for(j=0;j<NCwall;j++){
            g1=&(mesh->g_cells->cells[i]);
            c1=&(mesh->c_cells->cells[i]);

            carrays->solidWallByCell[j*ncells+i]=0;  

            if(j<g1->nneig){
                carrays->neighCell[j*ncells+i]=g1->neigcell[j]->id;
                carrays->neighWall[j*ncells+i]=g1->neigwall[j]->idWall;
                carrays->normalXbyCell[j*ncells+i]=g1->neigwall[j]->normal[_X_];
                carrays->normalYbyCell[j*ncells+i]=g1->neigwall[j]->normal[_Y_];     
            }else if(j<(g1->nneig+g1->nbound)){
                carrays->neighCell[j*ncells+i]=-1; //no neighcell
                carrays->neighWall[j*ncells+i]=g1->neigwall[j]->idWall;; //boundary wall
                carrays->normalXbyCell[j*ncells+i]=g1->neigwall[j]->normal[_X_];
                carrays->normalYbyCell[j*ncells+i]=g1->neigwall[j]->normal[_Y_];                
            }
            //printf("wall %d cell %d neighCell %d neighWall %d\n",j,i,carrays->neighCell[j*ncells+i],carrays->neighWall[j*ncells+i]);
		}
	}

	//Copy array values by walls
    nActWalls=0;
	for(i=0;i<nwc;i++){
        w1=&(mesh->w_calc->wall[i]);

        carrays->widx[i]=w1->idWall;

		carrays->idx1[i]=w1->gcells[0]->id;	
		carrays->idx2[i]=w1->gcells[1]->id;	
		carrays->idw1[i]=w1->i;	
		carrays->idw2[i]=w1->j;

		carrays->normalX[i]=w1->normal[_X_];
		carrays->normalY[i]=w1->normal[_Y_];
		carrays->deltaX[i]=w1->deltaX;
		carrays->length[i]=w1->length;
		carrays->distNormal[i]=w1->distNormal;
        carrays->distCentX[i]=w1->distCenter*w1->ncx;
        carrays->distCentY[i]=w1->distCenter*w1->ncy;

        //square of the wall-averaged nMan
        c1=&(mesh->c_cells->cells[w1->gcells[0]->id]);
        c2=&(mesh->c_cells->cells[w1->gcells[1]->id]);
        daux=0.5*(c1->nMan+c2->nMan);
        carrays->nman2wall[i]=daux*daux;

		carrays->gp[i]=_g_;

        carrays->typeOfBound[i]=w1->typeOfBound;
		carrays->solidWall[i]=0;
        
        carrays->qnormalL[i]=0.0;
        carrays->localDt[i]=1e6;

        // initial active walls -------------------------------------------
        carrays->activeW[i]=0;
        carrays->actWalls[i]=-1;
        if((w1->ccells[0]->h > TOL12) || (w1->ccells[1]->h > TOL12)){
            carrays->activeW[i]=1;
            carrays->actWalls[nActWalls]=w1->idWall;
            nActWalls++;
        }  
        // ----------------------------------------------------------------   

	}
    carrays->nActWalls=nActWalls;
    //printf("nActWalls %d \n",nActWalls);

	
	/*for(i=0;i<nw_calc;i++){	
		carrays->normalXbc[carrays->idw1[i]*ncells+mesh->idx1[i]]=mesh->normalX[i];
		carrays->normalYbc[carrays->idw1[i]*ncells+mesh->idx1[i]]=mesh->normalY[i];
		
		carrays->normalXbc[carrays->idw2[i]*ncells+mesh->idx2[i]]=-mesh->normalX[i];
		carrays->normalYbc[carrays->idw2[i]*ncells+mesh->idx2[i]]=-mesh->normalY[i];
	}*/
	
    return 1;

}





#if SET_SOLUTE
////////////////////////////////////////////////////////////////
EXPORT_DLL int allocateArraysSoluteMemory(
    t_parameters spar, 
    t_mesh *mesh,
    t_arrays *carrays,
    t_message *msg){
/*----------------------------*/

	int i,j;
    int NCwall;
	int ncells,nwc,nwb;
    int nWallCell;
    int nSolutes;

	size_t free_mem, total_mem;

	//Local variables just for allocation
    NCwall=mesh->NCwall;	//walls per cell
	ncells=mesh->ncells;
    nWallCell=NCwall*ncells;
	nwc=mesh->nw_calc;
	nwb=mesh->nw_bound;
    nSolutes=mesh->nSolutes;

    //solutes
    if(mesh->nSolutes){
        //solute arrays
        carrays->typeDiff=(int*)malloc(nSolutes*sizeof(int));
        carrays->k_xx=(double*)malloc(nSolutes*sizeof(double)); 
        carrays->k_yy=(double*)malloc(nSolutes*sizeof(double));  

        //cell arrarys
        carrays->hcsol=(double*)malloc(nSolutes*ncells*sizeof(double));
        carrays->csol=(double*)malloc(nSolutes*ncells*sizeof(double)); 

        //nWallCell arrays
        carrays->dhcsol=(double*)malloc(nSolutes*nWallCell*sizeof(double));    
    }

    return 1;

}


////////////////////////////////////////////////////////////////
EXPORT_DLL int initilizeControlSoluteArrays(
    t_parameters spar, 
    t_mesh *mesh,
    t_arrays *carrays,
    t_message *msg){
/*----------------------------*/

    int j;
    int nSolutes;

    nSolutes=mesh->nSolutes;

    if(mesh->nSolutes){
        carrays->flagDiffusion=mesh->solutes->flagDiffussion;

        //solute arrays
        for(j=0;j<nSolutes;j++){  
            carrays->typeDiff[j] = mesh->solutes->solute[j].typeDiff;
            carrays->k_xx[j] = mesh->solutes->solute[j].k_xx; 
            carrays->k_yy[j] = mesh->solutes->solute[j].k_yy;          
        }

    } 


    return 1;
 
}




////////////////////////////////////////////////////////////////
EXPORT_DLL int initilizeSoluteArrays(
    t_parameters spar, 
    t_mesh *mesh,
    t_arrays *carrays,   
    t_message *msg){
/*----------------------------*/

	int i,j;
	int NCwall,ncells,nwc,nwb,nWallCell;
    int nActCells, nActWalls;
    int nSolutes;
    int idx;


    double daux;

    t_c_cell *c1,*c2;
    t_g_cell *g1;
    t_wall *w1;

    size_t free_mem, total_mem;
	
	//Local variables just for the function
    NCwall=mesh->NCwall;
	ncells=mesh->ncells;
    nWallCell=carrays->nWallCell;
    nwc=mesh->nw_calc;
    nwb=mesh->nw_bound;
    nSolutes=mesh->nSolutes;

    if(mesh->nSolutes){
        //soltute*cell arrays
        for(j=0;j<nSolutes;j++){        
            for(i=0;i<ncells;i++){
                idx = j*ncells+i;
                c1=&(mesh->c_cells->cells[i]);

                carrays->hcsol[idx] = c1->hphi[j];
                carrays->csol[idx] = c1->phi[j];  
                //printf("Phi %d - Cell %d : %lf\n",j,i,carrays->csol[idx])  ;        
            }
        }

        //solute*cell*NCwall arrays
        for(i=0;i<nSolutes*nWallCell;i++){
		    carrays->dhcsol[i]=0.0;
	    }
      
    } 

    return 1;

}


#endif



