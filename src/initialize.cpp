#include "initialize.h"


////////////////////////////////////////////////////////////////
EXPORT_DLL int initilizeComputationControls(
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
	carrays->nInlet = mesh->nInlet; //permanent flag to switch on/off memory & computation
	carrays->nOutlet = mesh->nOutlet; //permanent flag to switch on/off memory & computation
    //carrays->nOBC = mesh->nInlet + mesh->nOutlet; //permanent flag to switch on/off memory & computation

    //solute permanent flag to switch on/off memory & computation
	carrays->nSolutes = mesh->nSolutes;

    return 1;
 
}



////////////////////////////////////////////////////////////////
EXPORT_DLL int allocateMeshArraysMem(
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
    carrays->typeWallByCell=(int*)malloc(nWallCell*sizeof(int));

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

    return 1;

}



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
                carrays->typeWallByCell[j*ncells+i]=g1->neigwall[j]->typeOfBound;

                carrays->normalXbyCell[j*ncells+i]=g1->neigwall[j]->normal[_X_];
                carrays->normalYbyCell[j*ncells+i]=g1->neigwall[j]->normal[_Y_]; 

            }else if(j<(g1->nneig+g1->nbound)){
                carrays->neighCell[j*ncells+i]=-1; //no neighcell
                carrays->neighWall[j*ncells+i]=g1->neigwall[j]->idWall;; //boundary wall
                carrays->typeWallByCell[j*ncells+i]=g1->neigwall[j]->typeOfBound;

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
	
    return 1;

}




////////////////////////////////////////////////////////////////
EXPORT_DLL int allocateBoundaryArraysMem(
    t_parameters spar, 
    t_mesh *mesh,
    t_arrays *carrays,
    t_message *msg){
/*----------------------------*/

	int i,j;
    int NCwall;
    int nInlet, nOutlet;
    int nTotalCellsIn, nTotalCellsOut, nTotalBoundCells;
    int nTotalInnerCells;
    int nTotalPointSeries;

    int nOBC;

    int nSolutes;

	//Local variables just for allocation
    NCwall = mesh->NCwall;	//walls per cell
    
    nInlet = mesh->nInlet;
    nOutlet = mesh->nOutlet;

    nTotalCellsIn = mesh->nTotalCellsIn;
    nTotalCellsOut = mesh->nTotalCellsOut;
    nTotalBoundCells = nTotalCellsIn+nTotalCellsOut;

    //nTotalInnerCells = mesh->nTotalInnerIn+mesh->nTotalInnerOut;

    nTotalPointSeries = mesh->nTotalSeriesIn+mesh->nTotalSeriesOut;

    nSolutes = carrays->nSolutes;

    //bound blocks
    nOBC=0;
    if(nInlet){
        for(j=0;j<nInlet;j++){  
            nOBC += mesh->in[j].ncellsBound/threadsPerOBC + 1;
        }
    }
    if(nOutlet){
        for(j=0;j<nOutlet;j++){
            nOBC += mesh->out[j].ncellsBound/threadsPerOBC + 1;
        }  
    }      
    carrays->nOBC = nOBC;    


    if(nOBC){
        //open bounds
        carrays->nCellsOBC=(int*)malloc(nOBC*sizeof(int));
        carrays->iniIndexOBC=(int*)malloc(nOBC*sizeof(int));
        carrays->idBoundOBC=(int*)malloc(nOBC*sizeof(int));
        carrays->typeOBC=(int*)malloc(nOBC*sizeof(int)); 
        carrays->flagInitializeOBC=(int*)malloc(nOBC*sizeof(int));

        carrays->blockSectionOBC=(double*)malloc(nOBC*sizeof(double));
        carrays->normalXOBC=(double*)malloc(nOBC*sizeof(double));
        carrays->normalYOBC=(double*)malloc(nOBC*sizeof(double));
        carrays->totalLengthOBC=(double*)malloc(nOBC*sizeof(double));
        carrays->totalAreaOBC=(double*)malloc(nOBC*sizeof(double));        
        
        carrays->cellZminOBC=(int*)malloc(nOBC*sizeof(int));

        // carrays->nInnerCellsOBC=(int*)malloc(nOBC*sizeof(int));
        // carrays->iniInnerIndexOBC=(int*)malloc(nOBC*sizeof(int));        

        //open bound cells
        //cells
        carrays->cidxBound=(int*)malloc(nTotalBoundCells*sizeof(int));
        carrays->zCellBound=(double*)malloc(nTotalBoundCells*sizeof(double));
        carrays->areaCellBound=(double*)malloc(nTotalBoundCells*sizeof(double));
        //walls
        carrays->nxWallBound=(double*)malloc(nTotalBoundCells*sizeof(double));
        carrays->nyWallBound=(double*)malloc(nTotalBoundCells*sizeof(double));
        carrays->lWallBound=(double*)malloc(nTotalBoundCells*sizeof(double));

        //inner cells
        // if(nTotalInnerCells){
        //     carrays->cidxInner=(int*)malloc(nTotalInnerCells*sizeof(int));
        // }
    
        //time series
        carrays->nPointsSeriesOBC=(int*)malloc(nOBC*sizeof(int));
        carrays->iniIndexSeriesOBC=(int*)malloc(nOBC*sizeof(int));

        carrays->tSeriesOBC=(double*)malloc(nTotalPointSeries*sizeof(double)); 
        carrays->qSeriesOBC=(double*)malloc(nTotalPointSeries*sizeof(double)); 
        carrays->hzSeriesOBC=(double*)malloc(nTotalPointSeries*sizeof(double));
        carrays->frSeriesOBC=(double*)malloc(nTotalPointSeries*sizeof(double)); 
        #if SET_SOLUTE
        carrays->phiSeriesOBC=(double*)malloc(nSolutes*nTotalPointSeries*sizeof(double)); 
        #endif

        //mass balance arrays
        carrays->qBoundByCell=(double*)malloc(nTotalBoundCells*sizeof(double));
        carrays->mBoundByCell=(double*)malloc(nTotalBoundCells*sizeof(double));

        // if(nTotalInnerCells){
        //     carrays->mInnerByCell=(double*)malloc(nTotalInnerCells*sizeof(double));
        // }

        if(nInlet){
            carrays->qInByInlet=(double*)malloc(nInlet*sizeof(double));
            carrays->mInByInlet=(double*)malloc(nInlet*sizeof(double));
        }
        if(nOutlet){
            carrays->qOutByOutlet=(double*)malloc(nOutlet*sizeof(double));
            carrays->mOutByOutlet=(double*)malloc(nOutlet*sizeof(double));
        }

    }

    return 1;

}



////////////////////////////////////////////////////////////////
EXPORT_DLL int initilizeBoundaryControlArrays(
    t_parameters spar, 
    t_mesh *mesh,
    t_arrays *carrays,
    t_message *msg){
/*----------------------------*/

    int j, k, m, l;
    int idx;
    int nInlet, nOutlet;
    int countOBC, countIdx0, countInnerIdx0;
    int bcBlocks, ncellsBlock; 
    int nTotalPointSeries;

    char temp[1024];

    nInlet=carrays->nInlet;
    nOutlet=carrays->nOutlet;

    //mandatory control initialization
    carrays->nTotalCellsIn = 0;
    carrays->nTotalCellsOut = 0;
    carrays->nTotalBoundCells = 0;

    carrays->nMaxBoundCells = 0;

    carrays->nTotalInnerCells = 0;

    carrays->nTotalPointSeries = 0;

    //open boundaries
    if(carrays->nOBC){
        carrays->nTotalCellsIn = mesh->nTotalCellsIn;
        carrays->nTotalCellsOut = mesh->nTotalCellsOut;
        carrays->nTotalBoundCells = mesh->nTotalCellsIn+mesh->nTotalCellsOut;
        //carrays->nTotalInnerCells = mesh->nTotalInnerIn+mesh->nTotalInnerOut;

        //bound data
        countOBC=0;
        countIdx0=0;
        //countInnerIdx0=0;
        if(nInlet){
            for(j=0;j<nInlet;j++){ 
                bcBlocks = mesh->in[j].ncellsBound/threadsPerOBC + 1; 
                for(m=0;m<bcBlocks;m++){

                    if(m < (bcBlocks-1)) {
                        ncellsBlock=threadsPerOBC;
                    }else{
                        ncellsBlock = mesh->in[j].ncellsBound-m*threadsPerOBC;
                    }
                    carrays->nCellsOBC[countOBC] = ncellsBlock;
                    carrays->blockSectionOBC[countOBC] = (double)ncellsBlock/(double)mesh->in[j].ncellsBound;
                    carrays->nMaxBoundCells = fmax(ncellsBlock,carrays->nMaxBoundCells);
                    //carrays->nCellsOBC[countOBC] = mesh->in[j].ncellsBound;

                    //OBC geometry features
                    carrays->iniIndexOBC[countOBC] = countIdx0;
                    carrays->idBoundOBC[countOBC] = -(j+1);
                    carrays->typeOBC[countOBC] = mesh->in[j].type;
                    carrays->flagInitializeOBC[countOBC] = 0;

                    //OBC direction
                    carrays->normalXOBC[countOBC] = mesh->in[j].normal.x; 
                    carrays->normalYOBC[countOBC] = mesh->in[j].normal.y;
                    carrays->totalLengthOBC[countOBC] = mesh->in[j].totalLength;
                    carrays->totalAreaOBC[countOBC] = mesh->in[j].totalArea;                 

                    carrays->totalLengthOBC[countOBC] *= carrays->blockSectionOBC[countOBC]; 
                    carrays->totalAreaOBC[countOBC] *= carrays->blockSectionOBC[countOBC];                

                    //OBC zmin cell index
                    carrays->cellZminOBC[countOBC] = mesh->in[j].cellzMin->id;

                    //OBC inner cells
                    //carrays->nInnerCellsOBC[countOBC] = mesh->in[j].ncellsInner;
                    //carrays->iniInnerIndexOBC[countOBC] = countInnerIdx0;

                    printf("bid %d block %d nbc %d blockSectionOBC %lf id0 %d\n",
                        carrays->idBoundOBC[countOBC],
                        m,
                        ncellsBlock,
                        carrays->blockSectionOBC[countOBC],
                        carrays->iniIndexOBC[countOBC]);    


                    //update index count
                    countOBC++;
                    countIdx0 += ncellsBlock;
                    //countIdx0 += mesh->in[j].ncellsBound;
                    //countInnerIdx0 += mesh->in[j].ncellsInner;
                }
            }
        }
        if(nOutlet){
            for(j=0;j<nOutlet;j++){ 
                bcBlocks = mesh->out[j].ncellsBound/threadsPerOBC + 1; 
                for(m=0;m<bcBlocks;m++){
                    //OBC geometry features
                    if(m < (bcBlocks-1)) {
                        ncellsBlock=threadsPerOBC;
                    }else{
                        ncellsBlock = mesh->out[j].ncellsBound-m*threadsPerOBC;
                    }
                    carrays->nCellsOBC[countOBC] = ncellsBlock;  
                    carrays->blockSectionOBC[countOBC] = (double)ncellsBlock/(double)mesh->out[j].ncellsBound;
                    carrays->nMaxBoundCells = fmax(ncellsBlock,carrays->nMaxBoundCells); 
                    //carrays->nCellsOBC[countOBC] = mesh->out[j].ncellsBound;

                    //OBC geometry features
                    carrays->iniIndexOBC[countOBC] = countIdx0;
                    carrays->idBoundOBC[countOBC] = (j+1);
                    carrays->typeOBC[countOBC] = mesh->out[j].type;
                    carrays->flagInitializeOBC[countOBC] = 0;

                    //OBC direction
                    carrays->normalXOBC[countOBC] = mesh->out[j].normal.x; 
                    carrays->normalYOBC[countOBC] = mesh->out[j].normal.y;
                    carrays->totalLengthOBC[countOBC] = mesh->out[j].totalLength;
                    carrays->totalAreaOBC[countOBC] = mesh->out[j].totalArea;

                    carrays->totalLengthOBC[countOBC] *= carrays->blockSectionOBC[countOBC]; 
                    carrays->totalAreaOBC[countOBC] *= carrays->blockSectionOBC[countOBC];                

                    //OBC zmin cell index
                    carrays->cellZminOBC[countOBC] = mesh->out[j].cellzMin->id;    

                    //OBC inner cells
                    //carrays->nInnerCellsOBC[countOBC] = mesh->out[j].ncellsInner;
                    //carrays->iniInnerIndexOBC[countOBC] = countInnerIdx0;

                    printf("bid %d block %d nbc %d blockSectionOBC %lf id0 %d\n",
                        carrays->idBoundOBC[countOBC],
                        m,
                        ncellsBlock,
                        carrays->blockSectionOBC[countOBC],
                        carrays->iniIndexOBC[countOBC]);                                                    

                    //update index count
                    countOBC++;
                    countIdx0 += ncellsBlock;
                    //countIdx0 += mesh->out[j].ncellsBound;
                    //countInnerIdx0 += mesh->out[j].ncellsInner;                
                }
            }  
        } 
        if((countOBC == carrays->nOBC) && (countIdx0 == carrays->nTotalBoundCells)){
			sprintf(temp,"Number of boundary block calculus: %d", carrays->nOBC);
			Notify(temp,MSG_L0,msg);	
			sprintf(temp,"Maximun number of cell in boundary block: %d", carrays->nMaxBoundCells);
			Notify(temp,MSG_L0,msg);                        
        }else{
            sprintf(temp,"Boundary block construction fails");
			Notify(temp,MSG_ERROR,msg);
            return 0;         
        }           

        
        //inlet time series
        nTotalPointSeries = mesh->nTotalSeriesIn+mesh->nTotalSeriesOut;
        carrays->nTotalPointSeries = nTotalPointSeries;

        countOBC=0;
        countIdx0=0;        
        if(nInlet){        
            for(j=0;j<nInlet;j++){
                bcBlocks = mesh->in[j].ncellsBound/threadsPerOBC + 1; 
                for(m=0;m<bcBlocks;m++){ 
                    carrays->nPointsSeriesOBC[countOBC] = mesh->in[j].n;
                    carrays->iniIndexSeriesOBC[countOBC] = countIdx0;

                    //update block count
                    countOBC++;
                }

                //time series                          
                switch(mesh->in[j].type){
                    case HYD_INFLOW_Q:
                        for(k=0;k<mesh->in[j].n;k++){
                            idx = countIdx0 + k;
                            carrays->tSeriesOBC[idx] = mesh->in[j].t[k];
                            carrays->qSeriesOBC[idx] = mesh->in[j].q[k];
                            carrays->hzSeriesOBC[idx] = 0.0;
                            carrays->frSeriesOBC[idx] = 0.0;
                        }
                        #if SET_SOLUTE
                        for(l=0;l<mesh->nSolutes;l++){
                            for(k=0;k<mesh->in[j].n;k++){
                                idx = l*nTotalPointSeries + countIdx0 + k;
                                carrays->phiSeriesOBC[idx] = mesh->in[j].phi[l][k];
                                //printf("cpu id %d sol %d phi %lf \n",j,l,carrays->phiSeriesOBC[idx]);

                            }
                        }
                        #endif
                        break;

                    case HYD_INFLOW_HZ://h+z(t)
                        for(k=0;k<mesh->in[j].n;k++){
                            idx = countIdx0 + k;
                            carrays->tSeriesOBC[idx] = mesh->in[j].t[k];
                            carrays->qSeriesOBC[idx] = 0.0;   
                            carrays->hzSeriesOBC[idx] = mesh->in[j].hZ[k];
                            carrays->frSeriesOBC[idx] = 0.0;
                        }
                        #if SET_SOLUTE
                        for(l=0;l<mesh->nSolutes;l++){
                            for(k=0;k<mesh->in[j].n;k++){
                                idx = l*nTotalPointSeries + countIdx0 + k;
                                carrays->phiSeriesOBC[idx] = mesh->in[j].phi[l][k];
                                //printf("cpu id %d sol %d phi %lf \n",j,l,carrays->phiSeriesOBC[idx]);

                            }
                        }
                        #endif
                        break;

                    case HYD_INFLOW_QHZ:
                        for(k=0;k<mesh->in[j].n;k++){
                            idx = countIdx0 + k;
                            carrays->tSeriesOBC[idx] = mesh->in[j].t[k];
                            carrays->qSeriesOBC[idx] = mesh->in[j].q[k];
                            carrays->hzSeriesOBC[idx] = mesh->in[j].hZ[k];
                            carrays->frSeriesOBC[idx] = 0.0;
                        }
                        #if SET_SOLUTE
                        for(l=0;l<mesh->nSolutes;l++){
                            for(k=0;k<mesh->in[j].n;k++){
                                idx = l*nTotalPointSeries + countIdx0 + k;
                                carrays->phiSeriesOBC[idx] = mesh->in[j].phi[l][k];
                                //printf("cpu id %d sol %d phi %lf \n",j,l,carrays->phiSeriesOBC[idx]);

                            }
                        }
                        #endif
                        break;
                }// End case

                //update time series index
                countIdx0 += mesh->in[j].n; 
            }
        } 
        if(nOutlet){
            for(j=0;j<nOutlet;j++){ 
                bcBlocks = mesh->out[j].ncellsBound/threadsPerOBC + 1; 
                for(m=0;m<bcBlocks;m++){ 
                    carrays->nPointsSeriesOBC[countOBC] = mesh->out[j].n;
                    carrays->iniIndexSeriesOBC[countOBC] = countIdx0;

                    //update block count
                    countOBC++;
                }                          

                //time series
                switch(mesh->out[j].type){
                    case HYD_OUTFLOW_GAUGE:
                        for(k=0;k<mesh->out[j].n;k++){
                            idx = countIdx0 + k;
                            carrays->tSeriesOBC[idx] = 0.0;
                            carrays->qSeriesOBC[idx] = mesh->out[j].q[k];
                            carrays->hzSeriesOBC[idx] = mesh->out[j].hZ[k];
                            carrays->frSeriesOBC[idx] = 0.0;
                        }
                        break;

                    case HYD_OUTFLOW_HZ:
                        for(k=0;k<mesh->out[j].n;k++){
                            idx = countIdx0 + k;
                            carrays->tSeriesOBC[idx] = mesh->out[j].t[k];
                            carrays->qSeriesOBC[idx] = 0.0;
                            carrays->hzSeriesOBC[idx] = mesh->out[j].hZ[k];
                            carrays->frSeriesOBC[idx] = 0.0;
                        }
                        break;
                    case HYD_OUTFLOW_FREE:
                        break;

                    case HYD_OUTFLOW_FR:
                            idx = countIdx0;
                            carrays->tSeriesOBC[idx] = 0.0;
                            carrays->qSeriesOBC[idx] = 0.0;
                            carrays->hzSeriesOBC[idx] = 0.0;
                            carrays->frSeriesOBC[idx] = mesh->out[j].Fr[0];                    
                        break;

                    case HYD_OUTFLOW_NORMAL:
                        for(k=0;k<mesh->out[j].n;k++){
                            idx = countIdx0 + k;
                            carrays->tSeriesOBC[idx] = 0.0;
                            carrays->qSeriesOBC[idx] = mesh->out[j].q[k];
                            carrays->hzSeriesOBC[idx] = mesh->out[j].hZ[k];
                            carrays->frSeriesOBC[idx] = 0.0;
                        }
                        break;
                }// End case  

                //update time series index
                countIdx0 += mesh->out[j].n;           
            }

        }           

    } 


    //mass balance initialization
    if(nInlet){
        for(j=0;j<nInlet;j++){
            carrays->qInByInlet[j]=0.0;
            carrays->mInByInlet[j]=0.0;
        }
    }
    if(nOutlet){
        for(j=0;j<nOutlet;j++){
            carrays->qOutByOutlet[j]=0.0;
            carrays->mOutByOutlet[j]=0.0;
        }
    }
    
    carrays->qTotalIn=0.0;
    carrays->qTotalOut=0.0;  
    carrays->mTotalIn=0.0;
    carrays->mTotalOut=0.0; 

    carrays->massTotalIn=0.0;
    carrays->massTotalIn=0.0;       

    return 1;
 
}



////////////////////////////////////////////////////////////////
EXPORT_DLL int initilizeBoundaryMeshArrays(
    t_parameters spar, 
    t_mesh *mesh,
    t_arrays *carrays,
    t_message *msg){
/*----------------------------*/

    int i,j;
    int nInlet, nOutlet;
    int nTotalCellsOBC;
    int countIdx;

	//Local variables just for allocation
    nInlet=mesh->nInlet;
    nOutlet=mesh->nOutlet;

    nTotalCellsOBC=mesh->nTotalCellsIn+mesh->nTotalCellsOut;

    if(carrays->nOBC){

        //Bound cells
        countIdx=0;
        if(nInlet){
            for(i=0;i<nInlet;i++){
                for(j=0;j<mesh->in[i].ncellsBound;j++){
                    carrays->cidxBound[countIdx] = mesh->in[i].wallBound[j]->ccells[0]->id;
                    carrays->zCellBound[countIdx] = mesh->in[i].wallBound[j]->ccells[0]->z;
                    carrays->areaCellBound[countIdx] = mesh->in[i].wallBound[j]->gcells[0]->area;

                    carrays->nxWallBound[countIdx] = mesh->in[i].wallBound[j]->normal[_X_];
                    carrays->nyWallBound[countIdx] = mesh->in[i].wallBound[j]->normal[_Y_];
                    carrays->lWallBound[countIdx] = mesh->in[i].wallBound[j]->length;

                    countIdx++;
                }
            }            
        }
        if(mesh->nOutlet){
            for(i=0;i<nOutlet;i++){
                for(j=0;j<mesh->out[i].ncellsBound;j++){
                    carrays->cidxBound[countIdx] = mesh->out[i].wallBound[j]->ccells[0]->id;
                    carrays->zCellBound[countIdx] = mesh->out[i].wallBound[j]->ccells[0]->z;
                    carrays->areaCellBound[countIdx] = mesh->out[i].wallBound[j]->gcells[0]->area;

                    carrays->nxWallBound[countIdx] = mesh->out[i].wallBound[j]->normal[_X_];
                    carrays->nyWallBound[countIdx] = mesh->out[i].wallBound[j]->normal[_Y_];                    
                    carrays->lWallBound[countIdx] = mesh->out[i].wallBound[j]->length;

                    countIdx++;
                }
            }
        } 

        //Inner cells
        // countIdx=0;
        // if(nInlet){
        //     for(i=0;i<nInlet;i++){
        //         for(j=0;j<mesh->in[i].ncellsInner;j++){
        //             carrays->cidxInner[countIdx] = mesh->in[i].cellInner[j]->id;
        //             countIdx++;
        //         }
        //     }            
        // }
        // if(mesh->nOutlet){
        //     for(i=0;i<nOutlet;i++){
        //         for(j=0;j<mesh->out[j].ncellsInner;j++){
        //             carrays->cidxInner[countIdx] = mesh->out[i].cellInner[j]->id;
        //             countIdx++;
        //         }
        //     }
        // }          

    } 

    return 1;
 
}


#if SET_SOLUTE
////////////////////////////////////////////////////////////////
EXPORT_DLL int allocateSoluteArraysMem(
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
        carrays->hphi=(double*)malloc(nSolutes*ncells*sizeof(double));
        carrays->phi=(double*)malloc(nSolutes*ncells*sizeof(double)); 
        carrays->localDtd=(double*)malloc(nSolutes*ncells*sizeof(double));
        carrays->BTcell=(double*)malloc(nSolutes*ncells*sizeof(double));

        //nWallCell arrays
        carrays->dhphi=(double*)malloc(nSolutes*nWallCell*sizeof(double));  
        carrays->Bwall=(double*)malloc(nSolutes*nWallCell*sizeof(double));    
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

        carrays->flagDiffusion=mesh->solutes->flagDiffussion;
        carrays->Dtd = 0.0; 

        //solute arrays
        for(j=0;j<nSolutes;j++){  
            carrays->typeDiff[j] = mesh->solutes->solute[j].typeDiff;
            carrays->k_xx[j] = mesh->solutes->solute[j].k_xx; 
            carrays->k_yy[j] = mesh->solutes->solute[j].k_yy;          
        }        

        
        //soltute*cell arrays
        for(j=0;j<nSolutes;j++){        
            for(i=0;i<ncells;i++){
                idx = j*ncells+i;
                c1=&(mesh->c_cells->cells[i]);

                carrays->hphi[idx] = c1->hphi[j];
                carrays->phi[idx] = c1->phi[j];  
                carrays->localDtd[idx]= 1e6;
                carrays->BTcell[idx]=0.0;
                //printf("Phi %d - Cell %d : %lf\n",j,i,carrays->csol[idx])  ;        
            }
        }

        //solute*cell*NCwall arrays
        for(i=0;i<nSolutes*nWallCell;i++){
		    carrays->dhphi[i]=0.0;
            carrays->Bwall[i]=0.0;
	    }
      
    } 

    return 1;

}

#endif



