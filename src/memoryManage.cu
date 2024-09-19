#include "memoryManage.cuh"

////////////////////////////////////////////////////////////////
EXPORT_DLL int createArraysCudaMemory(
    t_arrays *carrays,
    t_arrays *garrays,
    t_cuPtr *cuPtr,
    t_message *msg){
/*----------------------------*/

	int i,j;
    int NCwall;
	int ncells,nwc,nwb;
    int nWallCell;

	size_t free_mem, total_mem;

    int nTasks, blocksPerGrid;

	//Local variables just for allocation
    NCwall=carrays->NCwall;	//walls per cell
	ncells=carrays->ncells;
    nWallCell=NCwall*ncells;
	nwc=carrays->nw_calc;
	nwb=carrays->nw_bound;

    //allocate arrays memory
    if(!allocateArraysCudaMem( 
        NCwall, ncells, nWallCell, nwc, nwb, 
        cuPtr) ){
        return(0);
    };
    cudaDeviceSynchronize();
    //getCudaMemoryState(&free_mem, &total_mem);   



    //allocate  control arrays memory
    if(!copyControlArraysCudaMem(carrays, garrays, cuPtr)){
        return(0);
    };
    cudaDeviceSynchronize();
    //getCudaMemoryState(&free_mem, &total_mem); 


    //copy mesh arrays memory
    if(!copyMeshArraysCudaMem(carrays, garrays, cuPtr)){
        return(0);
    };
    cudaDeviceSynchronize();
    //getCudaMemoryState(&free_mem, &total_mem);     


    //assign CUDA garrays to CUDA cuPtr 
    assignGArraysToCudaMem <<<1,1>>> (garrays,
        //------------------------cells
        cuPtr->activeC,
        cuPtr->actCells,
        cuPtr->cidx,
        cuPtr->nneig,
        cuPtr->z,
        cuPtr->h,
        cuPtr->hu,
        cuPtr->hv,
        cuPtr->u,
        cuPtr->v,
        cuPtr->modulou,
        cuPtr->sqrh,
        cuPtr->area,
        cuPtr->nman,
        cuPtr->SOX,
        cuPtr->SOY,
        cuPtr->mass,
        //--------------------- cells*NCwall
        cuPtr->dh,
        cuPtr->dhu,
        cuPtr->dhv,
        cuPtr->solidWallByCell,
        cuPtr->neighCell,
        cuPtr->neighWall,
        cuPtr->normalXbyCell,
        cuPtr->normalYbyCell,
        //---------------------- internal walls
        cuPtr->activeW,
        cuPtr->actWalls,
        cuPtr->widx,
        cuPtr->idx1,
        cuPtr->idx2,
        cuPtr->idw1,
        cuPtr->idw2,
        cuPtr->normalX,
        cuPtr->normalY,
        cuPtr->deltaX,
        cuPtr->length,
        cuPtr->distNormal,
        cuPtr->distCentX,
        cuPtr->distCentY,
        cuPtr->nman2wall,
        cuPtr->gp,
        cuPtr->typeOfBound,
        cuPtr->solidWall,
        cuPtr->qnormalL,
        cuPtr->localDt
    );
    cudaDeviceSynchronize();
    getCudaMemoryState(&free_mem, &total_mem); 

    return 1;
}




////////////////////////////////////////////////////////////////
EXPORT_DLL int allocateArraysCudaMem(
    int NCwall, int ncells, int nWallCell, int nwc, int nwb,
    t_cuPtr *cuPtr){
/*----------------------------*/

    //computation controls
    cudaMalloc((void**) &(cuPtr->index), sizeof(int));
    cudaMalloc((void**) &(cuPtr->check), sizeof(int));

    cudaMalloc((void**) &(cuPtr->t), sizeof(double));
    cudaMalloc((void**) &(cuPtr->dt), sizeof(double));

    cudaMalloc((void**) &(cuPtr->nActCells), sizeof(int));
    cudaMalloc((void**) &(cuPtr->nActWalls), sizeof(int));

    cudaMalloc((void**) &(cuPtr->nIter), sizeof(int));
    cudaMalloc((void**) &(cuPtr->indexOut), sizeof(int));
    cudaMalloc((void**) &(cuPtr->indexDump), sizeof(int));
    cudaMalloc((void**) &(cuPtr->dumpComponent), sizeof(int));
    cudaMalloc((void**) &(cuPtr->dumpState), sizeof(int));

    cudaMalloc((void**) &(cuPtr->massOld), sizeof(double));
    cudaMalloc((void**) &(cuPtr->massNew), sizeof(double));
    cudaMalloc((void**) &(cuPtr->massError), sizeof(double));

    cudaMalloc((void**) &(cuPtr->qTotalIn), sizeof(double));
    cudaMalloc((void**) &(cuPtr->qTotalOut), sizeof(double));


    //cells
    cudaMalloc((void**) &(cuPtr->activeC), ncells*sizeof(int));
    cudaMalloc((void**) &(cuPtr->actCells), ncells*sizeof(int));

    cudaMalloc((void**) &(cuPtr->cidx), ncells*sizeof(int));
    cudaMalloc((void**) &(cuPtr->nneig), ncells*sizeof(int));

    cudaMalloc((void**) &(cuPtr->z), ncells*sizeof(double));
    cudaMalloc((void**) &(cuPtr->h), ncells*sizeof(double));
    cudaMalloc((void**) &(cuPtr->hu), ncells*sizeof(double));
    cudaMalloc((void**) &(cuPtr->hv), ncells*sizeof(double));
    cudaMalloc((void**) &(cuPtr->u), ncells*sizeof(double));
    cudaMalloc((void**) &(cuPtr->v), ncells*sizeof(double));
    cudaMalloc((void**) &(cuPtr->modulou), ncells*sizeof(double));
    cudaMalloc((void**) &(cuPtr->sqrh), ncells*sizeof(double));

    cudaMalloc((void**) &(cuPtr->area), ncells*sizeof(double));
    cudaMalloc((void**) &(cuPtr->nman), ncells*sizeof(double));
    cudaMalloc((void**) &(cuPtr->SOX), ncells*sizeof(double));
    cudaMalloc((void**) &(cuPtr->SOY), ncells*sizeof(double));    

    cudaMalloc((void**) &(cuPtr->mass), ncells*sizeof(double));


    //cells*NCwall
    cudaMalloc((void**) &(cuPtr->dh), nWallCell*sizeof(double));
    cudaMalloc((void**) &(cuPtr->dhu), nWallCell*sizeof(double));
    cudaMalloc((void**) &(cuPtr->dhv), nWallCell*sizeof(double));

    cudaMalloc((void**) &(cuPtr->solidWallByCell), nWallCell*sizeof(int));
    cudaMalloc((void**) &(cuPtr->neighCell), nWallCell*sizeof(int));
    cudaMalloc((void**) &(cuPtr->neighWall), nWallCell*sizeof(int));

    cudaMalloc((void**) &(cuPtr->normalXbyCell), nWallCell*sizeof(double));
    cudaMalloc((void**) &(cuPtr->normalYbyCell), nWallCell*sizeof(double));


    //internal walls
    cudaMalloc((void**) &(cuPtr->activeW), nwc*sizeof(int));
    cudaMalloc((void**) &(cuPtr->actWalls), nwc*sizeof(int));

    cudaMalloc((void**) &(cuPtr->widx), nwc*sizeof(int));

    cudaMalloc((void**) &(cuPtr->idx1), nwc*sizeof(int));
    cudaMalloc((void**) &(cuPtr->idx2), nwc*sizeof(int));
    cudaMalloc((void**) &(cuPtr->idw1), nwc*sizeof(int));
    cudaMalloc((void**) &(cuPtr->idw2), nwc*sizeof(int)); 

    cudaMalloc((void**) &(cuPtr->normalX), nwc*sizeof(double));
    cudaMalloc((void**) &(cuPtr->normalY), nwc*sizeof(double));
    cudaMalloc((void**) &(cuPtr->deltaX), nwc*sizeof(double));
    cudaMalloc((void**) &(cuPtr->length), nwc*sizeof(double));
    cudaMalloc((void**) &(cuPtr->distNormal), nwc*sizeof(double));
    cudaMalloc((void**) &(cuPtr->distCentX), nwc*sizeof(double));
    cudaMalloc((void**) &(cuPtr->distCentY), nwc*sizeof(double)); 

    cudaMalloc((void**) &(cuPtr->nman2wall), nwc*sizeof(double));
    cudaMalloc((void**) &(cuPtr->gp), nwc*sizeof(double));

    cudaMalloc((void**) &(cuPtr->typeOfBound), nwc*sizeof(int));
    cudaMalloc((void**) &(cuPtr->solidWall), nwc*sizeof(int));

    cudaMalloc((void**) &(cuPtr->qnormalL), nwc*sizeof(double));
    cudaMalloc((void**) &(cuPtr->localDt), nwc*sizeof(double));


    //boundaries

    
      
    return(1);

}



////////////////////////////////////////////////////////////////
EXPORT_DLL int copyControlArraysCudaMem(
    t_arrays *carrays,
    t_arrays *garrays,
    t_cuPtr *cuPtr){
/*----------------------------*/

    int iaux=0;

    // COMPUTATION CONTROLS ////////////////////////////////
    cudaMemcpy(cuPtr->index, &(iaux), sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy(cuPtr->check, &(iaux), sizeof(int), cudaMemcpyHostToDevice );

    cudaMemcpy(cuPtr->t , &(carrays->t), sizeof(double), cudaMemcpyHostToDevice );
    cudaMemcpy(cuPtr->dt, &(carrays->dt), sizeof(double), cudaMemcpyHostToDevice );

    cudaMemcpy(cuPtr->nActCells, &(carrays->nActCells), sizeof(int), cudaMemcpyHostToDevice ); //initial active cells
    cudaMemcpy(cuPtr->nActWalls, &(carrays->nActWalls), sizeof(int), cudaMemcpyHostToDevice ); //initial active walls

    cudaMemcpy(cuPtr->nIter, &(carrays->nIter), sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy(cuPtr->indexOut, &(carrays->indexOut), sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy(cuPtr->indexDump, &(carrays->indexDump), sizeof(int), cudaMemcpyHostToDevice ); 
    cudaMemcpy(cuPtr->dumpComponent, &(carrays->dumpComponent), sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy(cuPtr->dumpState, &(carrays->dumpState), sizeof(int), cudaMemcpyHostToDevice );

    cudaMemcpy(cuPtr->massOld, &(carrays->massOld), sizeof(double), cudaMemcpyHostToDevice );
    cudaMemcpy(cuPtr->massNew, &(carrays->massNew), sizeof(double), cudaMemcpyHostToDevice );
    cudaMemcpy(cuPtr->massError, &(carrays->massError), sizeof(double), cudaMemcpyHostToDevice );

    cudaMemcpy(cuPtr->qTotalIn, &(carrays->qTotalIn), sizeof(double), cudaMemcpyHostToDevice );
    cudaMemcpy(cuPtr->qTotalOut, &(carrays->qTotalOut), sizeof(double), cudaMemcpyHostToDevice );




    // COMPUTATION ARRAYS ////////////////////////////////
    cudaMemcpy(&(garrays->ncores), &(carrays->ncores), sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy(&(garrays->gpuid), &(carrays->gpuid), sizeof(int), cudaMemcpyHostToDevice );

    cudaMemcpy(&(garrays->ti), &(carrays->ti), sizeof(double), cudaMemcpyHostToDevice );
    cudaMemcpy(&(garrays->tf), &(carrays->tf), sizeof(double), cudaMemcpyHostToDevice );
    cudaMemcpy(&(garrays->CFL), &(carrays->CFL), sizeof(double), cudaMemcpyHostToDevice );

    cudaMemcpy(&(garrays->writeMass), &(carrays->writeMass), sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy(&(garrays->writeExtremes), &(carrays->writeExtremes), sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy(&(garrays->indexWriteHotstart), &(carrays->indexWriteHotstart), sizeof(int), cudaMemcpyHostToDevice );

    cudaMemcpy(&(garrays->nIterOut), &(carrays->nIterOut), sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy(&(garrays->dtOut), &(carrays->dtOut), sizeof(double), cudaMemcpyHostToDevice );
    cudaMemcpy(&(garrays->dtDump), &(carrays->dtDump), sizeof(double), cudaMemcpyHostToDevice );
    cudaMemcpy(&(garrays->minh), &(carrays->minh), sizeof(double), cudaMemcpyHostToDevice );

    //mesh
    cudaMemcpy(&(garrays->NCwall), &(carrays->NCwall), sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy(&(garrays->ncells), &(carrays->ncells), sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy(&(garrays->nWallCell), &(carrays->nWallCell), sizeof(int), cudaMemcpyHostToDevice );

    cudaMemcpy(&(garrays->nw_calc), &(carrays->nw_calc), sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy(&(garrays->nw_bound), &(carrays->nw_bound), sizeof(int), cudaMemcpyHostToDevice );

    cudaMemcpy(&(garrays->nInlet), &(carrays->nInlet), sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy(&(garrays->nOutlet), &(carrays->nOutlet), sizeof(int), cudaMemcpyHostToDevice );

    //computation controls
    cudaMemcpy(&(garrays->t), &(carrays->t), sizeof(double), cudaMemcpyHostToDevice );
    cudaMemcpy(&(garrays->dt), &(carrays->dt), sizeof(double), cudaMemcpyHostToDevice );

    cudaMemcpy(&(garrays->nActCells), &(carrays->nActCells), sizeof(int), cudaMemcpyHostToDevice ); //initial active cells
    cudaMemcpy(&(garrays->nActWalls), &(carrays->nActWalls), sizeof(int), cudaMemcpyHostToDevice ); //initial active walls

    cudaMemcpy(&(garrays->nIter), &(carrays->nIter), sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy(&(garrays->indexOut), &(carrays->indexOut), sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy(&(garrays->indexDump), &(carrays->indexDump), sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy(&(garrays->dumpComponent), &(carrays->dumpComponent), sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy(&(garrays->dumpState), &(carrays->dumpState), sizeof(int), cudaMemcpyHostToDevice );

    cudaMemcpy(&(garrays->massOld), &(carrays->massOld), sizeof(double), cudaMemcpyHostToDevice );
    cudaMemcpy(&(garrays->massNew), &(carrays->massNew), sizeof(double), cudaMemcpyHostToDevice );
    cudaMemcpy(&(garrays->massError), &(carrays->massError), sizeof(double), cudaMemcpyHostToDevice );

    //boundary control
    cudaMemcpy(&(garrays->qTotalIn), &(carrays->qTotalIn), sizeof(double), cudaMemcpyHostToDevice );
    cudaMemcpy(&(garrays->qTotalOut), &(carrays->qTotalOut), sizeof(double), cudaMemcpyHostToDevice );


    return 1;

}



////////////////////////////////////////////////////////////////
int copyMeshArraysCudaMem(
    t_arrays *carrays,
    t_arrays *garrays,
    t_cuPtr *cuPtr){
/*----------------------------*/
    
	int i,j;
    int NCwall;
	int ncells,nwc,nwb;
    int nWallCell;

	//Local variables just for allocation
    NCwall=carrays->NCwall;	//walls per cell
	ncells=carrays->ncells;
    nWallCell=NCwall*ncells;
	nwc=carrays->nw_calc;
	nwb=carrays->nw_bound;

    // COMPUTATION ARRAYS ////////////////////////////////////////////
    //cells    
    cudaMemcpy((cuPtr->activeC), (carrays->activeC), ncells*sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy((cuPtr->actCells), (carrays->actCells), ncells*sizeof(int), cudaMemcpyHostToDevice );    

    cudaMemcpy((cuPtr->cidx), (carrays->cidx), ncells*sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy((cuPtr->nneig), (carrays->nneig), ncells*sizeof(int), cudaMemcpyHostToDevice );

    cudaMemcpy((cuPtr->z), (carrays->z), ncells*sizeof(double), cudaMemcpyHostToDevice );
    cudaMemcpy((cuPtr->h), (carrays->h), ncells*sizeof(double), cudaMemcpyHostToDevice );
    cudaMemcpy((cuPtr->hu), (carrays->hu), ncells*sizeof(double), cudaMemcpyHostToDevice );
    cudaMemcpy((cuPtr->hv), (carrays->hv), ncells*sizeof(double), cudaMemcpyHostToDevice );
    cudaMemcpy((cuPtr->u), (carrays->u), ncells*sizeof(double), cudaMemcpyHostToDevice );
    cudaMemcpy((cuPtr->v), (carrays->v), ncells*sizeof(double), cudaMemcpyHostToDevice );
    cudaMemcpy((cuPtr->modulou), (carrays->modulou), ncells*sizeof(double), cudaMemcpyHostToDevice );
    cudaMemcpy((cuPtr->sqrh), (carrays->sqrh), ncells*sizeof(double), cudaMemcpyHostToDevice );
	
    cudaMemcpy((cuPtr->area), (carrays->area), ncells*sizeof(double), cudaMemcpyHostToDevice );
    cudaMemcpy((cuPtr->nman), (carrays->nman), ncells*sizeof(double), cudaMemcpyHostToDevice );
    cudaMemcpy((cuPtr->SOX), (carrays->SOX), ncells*sizeof(double), cudaMemcpyHostToDevice );
    cudaMemcpy((cuPtr->SOY), (carrays->SOY), ncells*sizeof(double), cudaMemcpyHostToDevice );    

    cudaMemcpy((cuPtr->mass), (carrays->mass), ncells*sizeof(double), cudaMemcpyHostToDevice );


    //cells*NCwall
    cudaMemcpy((cuPtr->dh), (carrays->dh), nWallCell*sizeof(double), cudaMemcpyHostToDevice );
    cudaMemcpy((cuPtr->dhu), (carrays->dhu), nWallCell*sizeof(double), cudaMemcpyHostToDevice );
    cudaMemcpy((cuPtr->dhv), (carrays->dhv), nWallCell*sizeof(double), cudaMemcpyHostToDevice );

    cudaMemcpy((cuPtr->solidWallByCell), (carrays->solidWallByCell), nWallCell*sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy((cuPtr->neighCell), (carrays->neighCell), nWallCell*sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy((cuPtr->neighWall), (carrays->neighWall), nWallCell*sizeof(int), cudaMemcpyHostToDevice );

    cudaMemcpy((cuPtr->normalXbyCell), (carrays->normalXbyCell), nWallCell*sizeof(double), cudaMemcpyHostToDevice );
    cudaMemcpy((cuPtr->normalYbyCell), (carrays->normalYbyCell), nWallCell*sizeof(double), cudaMemcpyHostToDevice );    


    //internal walls
    cudaMemcpy((cuPtr->activeW), (carrays->activeW), nwc*sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy((cuPtr->actWalls), (carrays->actWalls), nwc*sizeof(int), cudaMemcpyHostToDevice );

    cudaMemcpy((cuPtr->widx), (carrays->widx), nwc*sizeof(int), cudaMemcpyHostToDevice );

    cudaMemcpy((cuPtr->idx1), (carrays->idx1), nwc*sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy((cuPtr->idx2), (carrays->idx2), nwc*sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy((cuPtr->idw1), (carrays->idw1), nwc*sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy((cuPtr->idw2), (carrays->idw2), nwc*sizeof(int), cudaMemcpyHostToDevice );

    cudaMemcpy((cuPtr->normalX), (carrays->normalX), nwc*sizeof(double), cudaMemcpyHostToDevice );
    cudaMemcpy((cuPtr->normalY), (carrays->normalY), nwc*sizeof(double), cudaMemcpyHostToDevice ); 
    cudaMemcpy((cuPtr->deltaX), (carrays->deltaX), nwc*sizeof(double), cudaMemcpyHostToDevice ); 
    cudaMemcpy((cuPtr->length), (carrays->length), nwc*sizeof(double), cudaMemcpyHostToDevice ); 
    cudaMemcpy((cuPtr->distNormal), (carrays->distNormal), nwc*sizeof(double), cudaMemcpyHostToDevice ); 
    cudaMemcpy((cuPtr->distCentX), (carrays->distCentX), nwc*sizeof(double), cudaMemcpyHostToDevice ); 
    cudaMemcpy((cuPtr->distCentY), (carrays->distCentY), nwc*sizeof(double), cudaMemcpyHostToDevice ); 

    cudaMemcpy((cuPtr->nman2wall), (carrays->nman2wall), nwc*sizeof(double), cudaMemcpyHostToDevice ); 
    cudaMemcpy((cuPtr->gp), (carrays->gp), nwc*sizeof(double), cudaMemcpyHostToDevice ); 

    cudaMemcpy((cuPtr->typeOfBound), (carrays->typeOfBound), nwc*sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy((cuPtr->solidWall), (carrays->solidWall), nwc*sizeof(int), cudaMemcpyHostToDevice );
 
    cudaMemcpy((cuPtr->qnormalL), (carrays->qnormalL), nwc*sizeof(double), cudaMemcpyHostToDevice ); 
    cudaMemcpy((cuPtr->localDt), (carrays->localDt), nwc*sizeof(double), cudaMemcpyHostToDevice ); 

    //boundaries
  

    return(1);
    
}



////////////////////////////////////////////////////////////////
__global__ void assignGArraysToCudaMem(t_arrays *garrays,
	//------------------------cells
	int *activeC,
	int *actCells,
	int *cidx,
	int *nneig,
	double *z,
	double *h,
	double *hu,
	double *hv,
	double *u,
	double *v,
	double *modulou,
	double *sqrh,
	double *area,
	double *nman,
	double *SOX,
	double *SOY,
	double *mass,
	//--------------------- cells*NCwall
	double *dh,
	double *dhu,
	double *dhv,
	int *solidWallByCell,
	int *neighCell,
	int *neighWall,
	double *normalXbyCell,
	double *normalYbyCell,
	//---------------------- internal walls
	int *activeW,
	int *actWalls,
	int *widx,
	int *idx1,
	int *idx2,
	int *idw1,
	int *idw2,
	double *normalX,
	double *normalY,
	double *deltaX,
	double *length,
	double *distNormal,
	double *distCentX,
	double *distCentY,
	double *nman2wall,
	double *gp,
	int *typeOfBound,
	int *solidWall,
	double *qnormalL,
	double *localDt){
/*----------------------------*/
	//mesh
	garrays->activeC=activeC;
	garrays->actCells=actCells;

	garrays->cidx=cidx;
	garrays->nneig=nneig;

	garrays->z=z;
	garrays->h=h;
	garrays->hu=hu;
	garrays->hv=hv;
	garrays->u=u;
	garrays->v=v;
	garrays->modulou=modulou;
	garrays->sqrh=sqrh;

	garrays->area=area;
	garrays->nman=nman;
	garrays->SOX=SOX;
	garrays->SOY=SOY;

	garrays->mass=mass;

	//cells*NCwall
	garrays->dh=dh;
	garrays->dhu=dhu;
	garrays->dhv=dhv;
	garrays->solidWallByCell=solidWallByCell;
	garrays->neighCell=neighCell;
	garrays->neighWall=neighWall;
	garrays->normalXbyCell=normalXbyCell;
	garrays->normalYbyCell=normalYbyCell;

	//internal walls
	garrays->activeW=activeW;
	garrays->actWalls=actWalls;

	garrays->widx=widx;

	garrays->idx1=idx1;
	garrays->idx2=idx2;
	garrays->idw1=idw1;
	garrays->idw2=idw2;

	garrays->normalX=normalX;
	garrays->normalY=normalY;
	garrays->deltaX=deltaX;
	garrays->length=length;
	garrays->distNormal=distNormal;
	garrays->distCentX=distCentX;
	garrays->distCentY=distCentY;

	garrays->nman2wall=nman2wall;
	garrays->gp=gp;

	garrays->typeOfBound=typeOfBound;
	garrays->solidWall=solidWall;

	garrays->qnormalL=qnormalL;
	garrays->localDt=localDt;

}




////////////////////////////////////////////////////////////////
EXPORT_DLL int freeCudaMemory(t_cuPtr *cuPtr){
/*----------------------------*/

    //computation controls
    cudaFree(cuPtr->index);
    cudaFree(cuPtr->check);

    cudaFree(cuPtr->t);
    cudaFree(cuPtr->dt);

    cudaFree(cuPtr->nActCells);
    cudaFree(cuPtr->nActWalls);

    cudaFree(cuPtr->nIter);
    cudaFree(cuPtr->indexOut);
    cudaFree(cuPtr->indexDump);
    cudaFree(cuPtr->dumpComponent);
    cudaFree(cuPtr->dumpState);

    cudaFree(cuPtr->massOld);
    cudaFree(cuPtr->massNew);
    cudaFree(cuPtr->massError);

    cudaFree(cuPtr->qTotalIn);
    cudaFree(cuPtr->qTotalOut);
    


    //cells
    cudaFree(cuPtr->activeC);
    cudaFree(cuPtr->actCells);

    cudaFree(cuPtr->cidx);
    cudaFree(cuPtr->nneig);

    cudaFree(cuPtr->z);
    cudaFree(cuPtr->h);
    cudaFree(cuPtr->hu);
    cudaFree(cuPtr->hv);
    cudaFree(cuPtr->u);
    cudaFree(cuPtr->v);
    cudaFree(cuPtr->modulou);
    cudaFree(cuPtr->sqrh);

    cudaFree(cuPtr->area);
    cudaFree(cuPtr->nman);
    cudaFree(cuPtr->SOX);
    cudaFree(cuPtr->SOY);  

    cudaFree(cuPtr->mass);


    //cells*NCwall
    cudaFree(cuPtr->dh);
    cudaFree(cuPtr->dhu);
    cudaFree(cuPtr->dhv);

    cudaFree(cuPtr->solidWallByCell);
    cudaFree(cuPtr->neighCell);
    cudaFree(cuPtr->neighWall);

    cudaFree(cuPtr->normalXbyCell);
    cudaFree(cuPtr->normalYbyCell);


    //internal walls
    cudaFree(cuPtr->activeW);
    cudaFree(cuPtr->actWalls);

    cudaFree(cuPtr->widx);

    cudaFree(cuPtr->idx1);
    cudaFree(cuPtr->idx2);
    cudaFree(cuPtr->idw1);
    cudaFree(cuPtr->idw2);

    cudaFree(cuPtr->normalX);
    cudaFree(cuPtr->normalY);
    cudaFree(cuPtr->deltaX);
    cudaFree(cuPtr->length);
    cudaFree(cuPtr->distNormal);
    cudaFree(cuPtr->distCentX);
    cudaFree(cuPtr->distCentY);

    cudaFree(cuPtr->nman2wall);
    cudaFree(cuPtr->gp);

    cudaFree(cuPtr->typeOfBound);
    cudaFree(cuPtr->solidWall);

    cudaFree(cuPtr->qnormalL);
    cudaFree(cuPtr->localDt);


    //boundaries

    
      
    return(1);

}

