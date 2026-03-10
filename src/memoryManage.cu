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
    int nSolutes;

    int nInlet, nOutlet, nOBC;
    int nTotalBoundCells, nTotalInnerCells;
    int nTotalPointSeries;

	size_t free_mem, total_mem;

    int nTasks, blocksPerGrid;

	//Mesh variables for allocation
    NCwall=carrays->NCwall;	//walls per cell
	ncells=carrays->ncells;
    nWallCell=NCwall*ncells;
	nwc=carrays->nw_calc;
	nwb=carrays->nw_bound;
   
    //Bound variables for allocation
    nInlet = carrays->nInlet;
    nOutlet = carrays->nOutlet;
    nOBC = carrays->nOBC;
    nTotalBoundCells = carrays->nTotalBoundCells;
    nTotalInnerCells = carrays->nTotalInnerCells;
    nTotalPointSeries = carrays->nTotalPointSeries;

    //Solute variables for allocation
    nSolutes=carrays->nSolutes;

    //Transfer computation arrays 
    if(!copyComputationControls(
        carrays, 
        garrays, 
        cuPtr)
    ) return 0;

    if(!allocateArraysCudaMem( 
        NCwall, ncells, nWallCell, nwc, nwb, 
        cuPtr) 
    ) return 0;

    if(!copyMeshArraysCudaMem(
        carrays, 
        garrays, 
        cuPtr)
    ) return 0;      

    assignMeshArraysToCudaMem <<<1,1>>> (garrays,
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
        cuPtr->typeWallByCell,
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
    //getCudaMemoryState(&free_mem, &total_mem);        



    //Transfer boundary arrays 
    if(!allocateBoundArraysCudaMem( 
        nOBC, nInlet, nOutlet,
        nTotalBoundCells, nTotalInnerCells, 
        nTotalPointSeries,
        nSolutes,
        cuPtr) 
    ) return 0; 

    if(!copyBoundSetupArraysCudaMem(
        carrays, 
        garrays, 
        cuPtr)
    ) return 0;  

    if(!copyBoundMeshArraysCudaMem(
        carrays, 
        garrays, 
        cuPtr)
    ) return 0; 

    assignBoundArraysToCudaMem <<<1,1>>> (garrays,
        //------------------------ bound geometry
        cuPtr->nCellsOBC,
        cuPtr->iniIndexOBC,
        cuPtr->idBoundOBC,
        cuPtr->typeOBC,
        cuPtr->flagInitializeOBC,
        cuPtr->blockSectionOBC,
        cuPtr->normalXOBC,
        cuPtr->normalYOBC,
        cuPtr->totalLengthOBC,
        cuPtr->totalAreaOBC,
        cuPtr->cellZminOBC,
        cuPtr->nInnerCellsOBC,
        cuPtr->iniInnerIndexOBC,
        //----------------------- bound cells
        cuPtr->cidxBound,
        cuPtr->zCellBound,
        cuPtr->areaCellBound,        
        cuPtr->nxWallBound,
        cuPtr->nyWallBound,
        cuPtr->lWallBound,
        //----------------------- inner cells
        cuPtr->cidxInner,
        //----------------------- time series
        cuPtr->nPointsSeriesOBC,
        cuPtr->iniIndexSeriesOBC,
        cuPtr->tSeriesOBC,
        cuPtr->qSeriesOBC,
        cuPtr->hzSeriesOBC,
        cuPtr->frSeriesOBC,
        cuPtr->phiSeriesOBC,
        //----------------------- mass balance pointers
        cuPtr->qBoundByCell,
        cuPtr->mBoundByCell,
        cuPtr->mInnerByCell,
        cuPtr->qInByInlet,
        cuPtr->qOutByOutlet,
        cuPtr->mInByInlet, 
        cuPtr->mOutByOutlet
    );     

    cudaDeviceSynchronize();
    //getCudaMemoryState(&free_mem, &total_mem);  


    #if SET_SOLUTE
    if(!allocateSoluteArraysCudaMem( 
        nSolutes,
        NCwall, ncells, nWallCell, nwc, nwb, 
        cuPtr)
    ) return 0;


    if(!copySoluteArraysCudaMem(
        nSolutes, 
        carrays, garrays, cuPtr)
    ) return 0;        

    assignSoluteArraysToCudaMem <<<1,1>>> (nSolutes, garrays,
        //------------------------solutes
        cuPtr->typeDiff,
        cuPtr->k_xx,
        cuPtr->k_yy,
        //------------------------solutes*cells
        cuPtr->hphi,
        cuPtr->phi,
        cuPtr->localDtd,
        cuPtr->BTcell,
        //------------------------solutes*cells*NCwall
        cuPtr->dhphi,
        cuPtr->Bwall
    );
    // cudaDeviceSynchronize();
    // //getCudaMemoryState(&free_mem, &total_mem); 
    #endif

    //Chech final CUDA memory
    cudaDeviceSynchronize();
    getCudaMemoryState(&free_mem, &total_mem); 

    return 1;
}






////////////////////////////////////////////////////////////////
EXPORT_DLL int copyComputationControls(
    t_arrays *carrays,
    t_arrays *garrays,
    t_cuPtr *cuPtr){
/*----------------------------*/

    int iaux=0;

    // CONTROL ARRAYS ////////////////////////////////
    //simulation
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

    //execution 
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

    //boundary permanent flag
    cudaMemcpy(&(garrays->nInlet), &(carrays->nInlet), sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy(&(garrays->nOutlet), &(carrays->nOutlet), sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy(&(garrays->nOBC), &(carrays->nOBC), sizeof(int), cudaMemcpyHostToDevice );

    //solute permanent flag
    cudaMemcpy(&(garrays->nSolutes), &(carrays->nSolutes), sizeof(int), cudaMemcpyHostToDevice );



    // CONTROL POINTERS ////////////////////////////////
    cudaMalloc((void**) &(cuPtr->index), sizeof(int));
    cudaMalloc((void**) &(cuPtr->check), sizeof(int));
    cudaMemcpy(cuPtr->index, &(iaux), sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy(cuPtr->check, &(iaux), sizeof(int), cudaMemcpyHostToDevice );    

    cudaMalloc((void**) &(cuPtr->t), sizeof(double));
    cudaMalloc((void**) &(cuPtr->dt), sizeof(double));
    cudaMemcpy(cuPtr->t , &(carrays->t), sizeof(double), cudaMemcpyHostToDevice );
    cudaMemcpy(cuPtr->dt, &(carrays->dt), sizeof(double), cudaMemcpyHostToDevice );    

    cudaMalloc((void**) &(cuPtr->nActCells), sizeof(int));
    cudaMalloc((void**) &(cuPtr->nActWalls), sizeof(int));
    cudaMemcpy(cuPtr->nActCells, &(carrays->nActCells), sizeof(int), cudaMemcpyHostToDevice ); //initial active cells
    cudaMemcpy(cuPtr->nActWalls, &(carrays->nActWalls), sizeof(int), cudaMemcpyHostToDevice ); //initial active walls

    cudaMalloc((void**) &(cuPtr->nIter), sizeof(int));
    cudaMalloc((void**) &(cuPtr->indexOut), sizeof(int));
    cudaMalloc((void**) &(cuPtr->indexDump), sizeof(int));
    cudaMalloc((void**) &(cuPtr->dumpComponent), sizeof(int));
    cudaMalloc((void**) &(cuPtr->dumpState), sizeof(int));
    cudaMemcpy(cuPtr->nIter, &(carrays->nIter), sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy(cuPtr->indexOut, &(carrays->indexOut), sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy(cuPtr->indexDump, &(carrays->indexDump), sizeof(int), cudaMemcpyHostToDevice ); 
    cudaMemcpy(cuPtr->dumpComponent, &(carrays->dumpComponent), sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy(cuPtr->dumpState, &(carrays->dumpState), sizeof(int), cudaMemcpyHostToDevice );    

    cudaMalloc((void**) &(cuPtr->massOld), sizeof(double));
    cudaMalloc((void**) &(cuPtr->massNew), sizeof(double));
    cudaMalloc((void**) &(cuPtr->massError), sizeof(double));
    cudaMemcpy(cuPtr->massOld, &(carrays->massOld), sizeof(double), cudaMemcpyHostToDevice );
    cudaMemcpy(cuPtr->massNew, &(carrays->massNew), sizeof(double), cudaMemcpyHostToDevice );
    cudaMemcpy(cuPtr->massError, &(carrays->massError), sizeof(double), cudaMemcpyHostToDevice );    

    return 1;

}





////////////////////////////////////////////////////////////////
EXPORT_DLL int allocateArraysCudaMem(
    int NCwall, int ncells, int nWallCell, int nwc, int nwb,
    t_cuPtr *cuPtr){
/*----------------------------*/

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
    cudaMalloc((void**) &(cuPtr->typeWallByCell), nWallCell*sizeof(int));

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
      
    return(1);

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
    cudaMemcpy((cuPtr->typeWallByCell), (carrays->typeWallByCell), nWallCell*sizeof(int), cudaMemcpyHostToDevice );

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

    return(1);
    
}


////////////////////////////////////////////////////////////////
__global__ void assignMeshArraysToCudaMem(t_arrays *garrays,
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
    int *typeWallByCell,
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
    garrays->typeWallByCell=typeWallByCell;

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
    cudaFree(cuPtr->typeWallByCell);

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
      
    return(1);

}





////////////////////////////////////////////////////////////////
EXPORT_DLL int allocateBoundArraysCudaMem(
    int nOBC, int nInlet, int nOutlet, 
    int nTotalBoundCells, int nTotalInnerCells,
    int nTotalPointSeries, 
    int nSolutes,
    t_cuPtr *cuPtr){
/*----------------------------*/


    if(nOBC){

        //bound geometry
        cudaMalloc((void**) &(cuPtr->nCellsOBC), nOBC*sizeof(int));
        cudaMalloc((void**) &(cuPtr->iniIndexOBC), nOBC*sizeof(int));
        cudaMalloc((void**) &(cuPtr->idBoundOBC), nOBC*sizeof(int));
        cudaMalloc((void**) &(cuPtr->typeOBC), nOBC*sizeof(int));
        cudaMalloc((void**) &(cuPtr->flagInitializeOBC), nOBC*sizeof(int));

        cudaMalloc((void**) &(cuPtr->blockSectionOBC), nOBC*sizeof(double));
        cudaMalloc((void**) &(cuPtr->normalXOBC), nOBC*sizeof(double));
        cudaMalloc((void**) &(cuPtr->normalYOBC), nOBC*sizeof(double));
        cudaMalloc((void**) &(cuPtr->totalLengthOBC), nOBC*sizeof(double));
        cudaMalloc((void**) &(cuPtr->totalAreaOBC), nOBC*sizeof(double));

        cudaMalloc((void**) &(cuPtr->cellZminOBC), nOBC*sizeof(int));
        cudaMalloc((void**) &(cuPtr->nInnerCellsOBC), nOBC*sizeof(int));
        cudaMalloc((void**) &(cuPtr->iniInnerIndexOBC), nOBC*sizeof(int));

        //bound cells
        cudaMalloc((void**) &(cuPtr->cidxBound), nTotalBoundCells*sizeof(int));
        cudaMalloc((void**) &(cuPtr->zCellBound), nTotalBoundCells*sizeof(double));
        cudaMalloc((void**) &(cuPtr->areaCellBound), nTotalBoundCells*sizeof(double));
        cudaMalloc((void**) &(cuPtr->nxWallBound), nTotalBoundCells*sizeof(double));
        cudaMalloc((void**) &(cuPtr->nyWallBound), nTotalBoundCells*sizeof(double));
        cudaMalloc((void**) &(cuPtr->lWallBound), nTotalBoundCells*sizeof(double));

        // cudaMalloc((void**) &(cuPtr->localh), nTotalBoundCells*sizeof(double));
        // cudaMalloc((void**) &(cuPtr->localhu), nTotalBoundCells*sizeof(double));
        // cudaMalloc((void**) &(cuPtr->localhv), nTotalBoundCells*sizeof(double));
        // #if SET_SOLUTE
        // cudaMalloc((void**) &(cuPtr->localhphi), nSolutes*nTotalBoundCells*sizeof(double));
        // #endif
        
        //inner cells
        // if(nTotalInnerCells){
        //     cudaMalloc((void**) &(cuPtr->cidxInner), nTotalInnerCells*sizeof(int));
        // }

        //time series
        cudaMalloc((void**) &(cuPtr->nPointsSeriesOBC), nOBC*sizeof(int)); 
        cudaMalloc((void**) &(cuPtr->iniIndexSeriesOBC), nOBC*sizeof(int));

        cudaMalloc((void**) &(cuPtr->tSeriesOBC), nTotalPointSeries*sizeof(double));
        cudaMalloc((void**) &(cuPtr->qSeriesOBC), nTotalPointSeries*sizeof(double));
        cudaMalloc((void**) &(cuPtr->hzSeriesOBC), nTotalPointSeries*sizeof(double));
        cudaMalloc((void**) &(cuPtr->frSeriesOBC), nTotalPointSeries*sizeof(double));

        #if SET_SOLUTE
        cudaMalloc((void**) &(cuPtr->phiSeriesOBC), nSolutes*nTotalPointSeries*sizeof(double));
        #endif

        //mass balance pointers
        cudaMalloc((void**) &(cuPtr->qBoundByCell), nTotalBoundCells*sizeof(double));
        cudaMalloc((void**) &(cuPtr->mBoundByCell), nTotalBoundCells*sizeof(double));

        // if(nTotalInnerCells){
        //     cudaMalloc((void**) &(cuPtr->mInnerByCell), nTotalInnerCells*sizeof(double));
        // }

        if(nInlet){
            cudaMalloc((void**) &(cuPtr->qInByInlet), nInlet*sizeof(double));
            cudaMalloc((void**) &(cuPtr->mInByInlet), nInlet*sizeof(double));
            cudaMalloc((void**) &(cuPtr->aux1sByInlet), nInlet*sizeof(double));
        }
        if(nOutlet){
            cudaMalloc((void**) &(cuPtr->qOutByOutlet), nOutlet*sizeof(double));
            cudaMalloc((void**) &(cuPtr->mOutByOutlet), nOutlet*sizeof(double));
            cudaMalloc((void**) &(cuPtr->aux1sByOutlet), nOutlet*sizeof(double));
        }

        cudaMalloc((void**) &(cuPtr->qTotalIn), sizeof(double));
        cudaMalloc((void**) &(cuPtr->qTotalOut), sizeof(double)); 
        cudaMalloc((void**) &(cuPtr->mTotalIn), sizeof(double));
        cudaMalloc((void**) &(cuPtr->mTotalOut), sizeof(double)); 

        cudaMalloc((void**) &(cuPtr->massTotalIn), sizeof(double));
        cudaMalloc((void**) &(cuPtr->massTotalOut), sizeof(double));

    }  

    return 1;    

}


////////////////////////////////////////////////////////////////
EXPORT_DLL int copyBoundSetupArraysCudaMem(
    t_arrays *carrays,
    t_arrays *garrays,
    t_cuPtr *cuPtr){
/*----------------------------*/

    int nOBC = carrays->nOBC;
    int nInlet = carrays->nInlet;
    int nOutlet = carrays->nOutlet;
    int nTotalPointSeries = carrays->nTotalPointSeries;
    int nSolutes = carrays->nSolutes;

    int i;
    double *aux1s;

    //ARRAYS ----------------------------------
    cudaMemcpy(&(garrays->nTotalCellsIn), &(carrays->nTotalCellsIn), sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy(&(garrays->nTotalCellsOut), &(carrays->nTotalCellsOut), sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy(&(garrays->nTotalBoundCells), &(carrays->nTotalBoundCells), sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy(&(garrays->nTotalBoundCells), &(carrays->nTotalBoundCells), sizeof(int), cudaMemcpyHostToDevice );

    cudaMemcpy(&(garrays->nMaxBoundCells), &(carrays->nMaxBoundCells), sizeof(int), cudaMemcpyHostToDevice );  
 
    cudaMemcpy(&(garrays->nTotalPointSeries), &(carrays->nTotalPointSeries), sizeof(int), cudaMemcpyHostToDevice );

    cudaMemcpy(&(garrays->qTotalIn), &(carrays->qTotalIn), sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy(&(garrays->qTotalOut), &(carrays->qTotalOut), sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy(&(garrays->massTotalIn), &(carrays->massTotalIn), sizeof(int), cudaMemcpyHostToDevice );  
    cudaMemcpy(&(garrays->massTotalIn), &(carrays->massTotalIn), sizeof(int), cudaMemcpyHostToDevice );
  
    // POINTERS ----------------------------------
    if(nOBC){

        //bound geometry
        cudaMemcpy((cuPtr->nCellsOBC), (carrays->nCellsOBC), nOBC*sizeof(int), cudaMemcpyHostToDevice );
        cudaMemcpy((cuPtr->iniIndexOBC), (carrays->iniIndexOBC), nOBC*sizeof(int), cudaMemcpyHostToDevice );
        cudaMemcpy((cuPtr->idBoundOBC), (carrays->idBoundOBC), nOBC*sizeof(int), cudaMemcpyHostToDevice );
        cudaMemcpy((cuPtr->typeOBC), (carrays->typeOBC), nOBC*sizeof(int), cudaMemcpyHostToDevice );
        cudaMemcpy((cuPtr->flagInitializeOBC), (carrays->flagInitializeOBC), nOBC*sizeof(int), cudaMemcpyHostToDevice );

        cudaMemcpy((cuPtr->blockSectionOBC), (carrays->blockSectionOBC), nOBC*sizeof(double), cudaMemcpyHostToDevice );
        cudaMemcpy((cuPtr->normalXOBC), (carrays->normalXOBC), nOBC*sizeof(double), cudaMemcpyHostToDevice );
        cudaMemcpy((cuPtr->normalYOBC), (carrays->normalYOBC), nOBC*sizeof(double), cudaMemcpyHostToDevice );
        cudaMemcpy((cuPtr->totalLengthOBC), (carrays->totalLengthOBC), nOBC*sizeof(double), cudaMemcpyHostToDevice );
        cudaMemcpy((cuPtr->totalAreaOBC), (carrays->totalAreaOBC), nOBC*sizeof(double), cudaMemcpyHostToDevice );

        cudaMemcpy((cuPtr->cellZminOBC), (carrays->cellZminOBC), nOBC*sizeof(int), cudaMemcpyHostToDevice );
        //cudaMemcpy((cuPtr->nInnerCellsOBC), (carrays->nInnerCellsOBC), nOBC*sizeof(int), cudaMemcpyHostToDevice );
        //cudaMemcpy((cuPtr->iniInnerIndexOBC), (carrays->iniInnerIndexOBC), nOBC*sizeof(int), cudaMemcpyHostToDevice );

        //time series
        cudaMemcpy((cuPtr->nPointsSeriesOBC), (carrays->nPointsSeriesOBC), nOBC*sizeof(int), cudaMemcpyHostToDevice );
        cudaMemcpy((cuPtr->iniIndexSeriesOBC), (carrays->iniIndexSeriesOBC), nOBC*sizeof(int), cudaMemcpyHostToDevice );
        
        cudaMemcpy((cuPtr->tSeriesOBC), (carrays->tSeriesOBC), nTotalPointSeries*sizeof(double), cudaMemcpyHostToDevice );
        cudaMemcpy((cuPtr->qSeriesOBC), (carrays->qSeriesOBC), nTotalPointSeries*sizeof(double), cudaMemcpyHostToDevice );
        cudaMemcpy((cuPtr->hzSeriesOBC), (carrays->hzSeriesOBC), nTotalPointSeries*sizeof(double), cudaMemcpyHostToDevice );
        cudaMemcpy((cuPtr->frSeriesOBC), (carrays->frSeriesOBC), nTotalPointSeries*sizeof(double), cudaMemcpyHostToDevice );

        #if SET_SOLUTE
        cudaMemcpy((cuPtr->phiSeriesOBC), (carrays->phiSeriesOBC), nSolutes*nTotalPointSeries*sizeof(double), cudaMemcpyHostToDevice );
        #endif

        //mass balance pointers  
        if(nInlet){
            cudaMemcpy(cuPtr->qInByInlet, &(carrays->qInByInlet), nInlet*sizeof(double), cudaMemcpyHostToDevice );
            cudaMemcpy(cuPtr->mInByInlet, &(carrays->mInByInlet), nInlet*sizeof(double), cudaMemcpyHostToDevice );

            aux1s = (double*)malloc(nInlet*sizeof(double));
            for(i=0;i<nInlet;i++){
                aux1s[i]=1.;
            }
            cudaMemcpy(cuPtr->aux1sByInlet, aux1s, nInlet*sizeof(double), cudaMemcpyHostToDevice );
            free(aux1s);
        }
        if(nOutlet){
            cudaMemcpy(cuPtr->qOutByOutlet, &(carrays->qOutByOutlet), nOutlet*sizeof(double), cudaMemcpyHostToDevice );      
            cudaMemcpy(cuPtr->mOutByOutlet, &(carrays->mOutByOutlet), nOutlet*sizeof(double), cudaMemcpyHostToDevice );  

            aux1s = (double*)malloc(nOutlet*sizeof(double));
            for(i=0;i<nOutlet;i++){
                aux1s[i]=1.;
            }
            cudaMemcpy(cuPtr->aux1sByOutlet, aux1s, nOutlet*sizeof(double), cudaMemcpyHostToDevice );
            free(aux1s);                
        }

        cudaMemcpy(cuPtr->qTotalIn, &(carrays->qTotalIn), sizeof(double), cudaMemcpyHostToDevice );
        cudaMemcpy(cuPtr->qTotalOut, &(carrays->qTotalOut), sizeof(double), cudaMemcpyHostToDevice );
        cudaMemcpy(cuPtr->mTotalIn, &(carrays->mTotalIn), sizeof(double), cudaMemcpyHostToDevice );
        cudaMemcpy(cuPtr->mTotalOut, &(carrays->mTotalOut), sizeof(double), cudaMemcpyHostToDevice );

        cudaMemcpy(cuPtr->massTotalIn, &(carrays->massTotalIn), sizeof(double), cudaMemcpyHostToDevice );
        cudaMemcpy(cuPtr->massTotalOut, &(carrays->massTotalOut), sizeof(double), cudaMemcpyHostToDevice );        
    }

    return 1;

}


////////////////////////////////////////////////////////////////
EXPORT_DLL int copyBoundMeshArraysCudaMem(
    t_arrays *carrays,
    t_arrays *garrays,
    t_cuPtr *cuPtr){
/*----------------------------*/

    int nOBC = carrays->nOBC;
    int nTotalBoundCells = carrays->nTotalBoundCells;
    int nTotalInnerCells = carrays->nTotalInnerCells;

    if(nOBC){
        //bound cells
        cudaMemcpy((cuPtr->cidxBound), (carrays->cidxBound), nTotalBoundCells*sizeof(int), cudaMemcpyHostToDevice );
        cudaMemcpy((cuPtr->zCellBound), (carrays->zCellBound), nTotalBoundCells*sizeof(double), cudaMemcpyHostToDevice );
        cudaMemcpy((cuPtr->areaCellBound), (carrays->areaCellBound), nTotalBoundCells*sizeof(double), cudaMemcpyHostToDevice );

        //bound walls
        cudaMemcpy((cuPtr->nxWallBound), (carrays->nxWallBound), nTotalBoundCells*sizeof(double), cudaMemcpyHostToDevice );
        cudaMemcpy((cuPtr->nyWallBound), (carrays->nyWallBound), nTotalBoundCells*sizeof(double), cudaMemcpyHostToDevice );
        cudaMemcpy((cuPtr->lWallBound), (carrays->lWallBound), nTotalBoundCells*sizeof(double), cudaMemcpyHostToDevice );

        //inner cells
        // if(nTotalInnerCells){
        //     cudaMemcpy((cuPtr->cidxInner), (carrays->cidxInner), nTotalInnerCells*sizeof(int), cudaMemcpyHostToDevice );
        // }

    }

    return 1;

}



////////////////////////////////////////////////////////////////
__global__ void assignBoundArraysToCudaMem(t_arrays *garrays,
	//------------------------ bound geometry
	int *nCellsOBC,
    int *iniIndexOBC,
    int *idBoundOBC,
    int *typeOBC,
    int *flagInitializeOBC,
    double *blockSectionOBC,
    double *normalXOBC,
    double *normalYOBC,
    double *totalLengthOBC,
    double *totalAreaOBC,
    int *cellZminOBC,
    int *nInnerCellsOBC,
    int *iniInnerIndexOBC,
    //----------------------- bound cells
    int *cidxBound,
    double *zCellBound,
    double *areaCellBound,
    double *nxWallBound,
    double *nyWallBound,
    double *lWallBound,
    //----------------------- inner cells
    int *cidxInner,
    //----------------------- time series
    int *nPointsSeriesOBC,
    int *iniIndexSeriesOBC,
    double *tSeriesOBC,
    double *qSeriesOBC,
    double *hzSeriesOBC,
    double *frSeriesOBC,
    double *phiSeriesOBC,
    //----------------------- mass balance pointers
    double *qBoundByCell,
    double *mBoundByCell,
    double *mInnerByCell,
    double *qInByInlet,
    double *qOutByOutlet,
    double *mInByInlet, 
    double *mOutByOutlet){
/*----------------------------*/

    if(garrays->nOBC){

        //bound geometry
        garrays->nCellsOBC=nCellsOBC;
        garrays->iniIndexOBC=iniIndexOBC;
        garrays->idBoundOBC=idBoundOBC;
        garrays->typeOBC=typeOBC;
        garrays->flagInitializeOBC=flagInitializeOBC;

        garrays->blockSectionOBC=blockSectionOBC;
        garrays->normalXOBC=normalXOBC;
        garrays->normalYOBC=normalYOBC;
        garrays->totalLengthOBC=totalLengthOBC;
        garrays->totalAreaOBC=totalAreaOBC;

        garrays->cellZminOBC=cellZminOBC;
        garrays->nInnerCellsOBC=nInnerCellsOBC;
        garrays->iniInnerIndexOBC=iniInnerIndexOBC;

        //bound cells
        garrays->cidxBound=cidxBound;
        garrays->zCellBound=zCellBound;
        garrays->areaCellBound=areaCellBound;
        garrays->nxWallBound=nxWallBound;
        garrays->nyWallBound=nyWallBound;
        garrays->lWallBound=lWallBound;  

        //inner cells 
        // if(garrays->nTotalInnerCells){
        //     garrays->cidxInner=cidxInner;
        // }     

        //time series
        garrays->nPointsSeriesOBC=nPointsSeriesOBC;
        garrays->iniIndexSeriesOBC=iniIndexSeriesOBC;

        garrays->tSeriesOBC=tSeriesOBC;
        garrays->qSeriesOBC=qSeriesOBC;
        garrays->hzSeriesOBC=hzSeriesOBC;
        garrays->frSeriesOBC=frSeriesOBC;

        #if SET_SOLUTE
        garrays->phiSeriesOBC=phiSeriesOBC;
        #endif

        //mass balance pointers
        garrays->qBoundByCell=qBoundByCell;
        garrays->mBoundByCell=mBoundByCell;
        // if(garrays->nTotalInnerCells){
        //     garrays->mInnerByCell=mInnerByCell;
        // }   

        if(garrays->nInlet){
            garrays->qInByInlet=qInByInlet;
            garrays->mInByInlet=mInByInlet;
        }
        if(garrays->nOutlet){
            garrays->qOutByOutlet=qOutByOutlet;
            garrays->mOutByOutlet=mOutByOutlet;
        }

    }

}


////////////////////////////////////////////////////////////////
EXPORT_DLL int freeBoundaCudaMemory(
    int nOBC, int nInlet, int nOutlet, 
    int nTotalBoundCells, int nTotalInnerCells,
    int nTotalPointSeries,
    t_cuPtr *cuPtr){
/*----------------------------*/

    if(nOBC){
        cudaFree(cuPtr->nCellsOBC);
        cudaFree(cuPtr->iniIndexOBC);
        cudaFree(cuPtr->idBoundOBC);
        cudaFree(cuPtr->typeOBC);
        cudaFree(cuPtr->flagInitializeOBC);


        cudaFree(cuPtr->blockSectionOBC);
        cudaFree(cuPtr->normalXOBC);
        cudaFree(cuPtr->normalYOBC);
        cudaFree(cuPtr->totalLengthOBC);
        cudaFree(cuPtr->totalAreaOBC);

        cudaFree(cuPtr->cellZminOBC);
        cudaFree(cuPtr->nInnerCellsOBC);
        cudaFree(cuPtr->iniInnerIndexOBC);

        cudaFree(cuPtr->cidxBound);
        cudaFree(cuPtr->zCellBound);
        cudaFree(cuPtr->areaCellBound);
        cudaFree(cuPtr->nxWallBound);
        cudaFree(cuPtr->nyWallBound);
        cudaFree(cuPtr->lWallBound);

        if(nTotalInnerCells){
            cudaFree(cuPtr->cidxInner);
        }

        cudaFree(cuPtr->nPointsSeriesOBC);
        cudaFree(cuPtr->iniIndexSeriesOBC);

        cudaFree(cuPtr->tSeriesOBC);
        cudaFree(cuPtr->qSeriesOBC);
        cudaFree(cuPtr->hzSeriesOBC);
        cudaFree(cuPtr->frSeriesOBC);


        #if SET_SOLUTE
        cudaFree(cuPtr->phiSeriesOBC);
        #endif        

        cudaFree(cuPtr->qBoundByCell);
        cudaFree(cuPtr->mBoundByCell);

        if(nTotalInnerCells){
            cudaFree(cuPtr->mInnerByCell);
        }

        if(nInlet){
            cudaFree(cuPtr->qInByInlet);
            cudaFree(cuPtr->mInByInlet);
            cudaFree(cuPtr->aux1sByInlet);
        }
        if(nOutlet){
            cudaFree(cuPtr->qOutByOutlet);
            cudaFree(cuPtr->mOutByOutlet);
            cudaFree(cuPtr->aux1sByOutlet);
        }

        cudaFree(cuPtr->qTotalIn);
        cudaFree(cuPtr->qTotalOut);
        cudaFree(cuPtr->mTotalIn);
        cudaFree(cuPtr->mTotalOut);

        cudaFree(cuPtr->massTotalIn);
        cudaFree(cuPtr->massTotalOut);         
    } 

    return(1);    

}





#if SET_SOLUTE
////////////////////////////////////////////////////////////////
EXPORT_DLL int allocateSoluteArraysCudaMem(
    int nSolutes, 
    int NCwall, int ncells, int nWallCell, int nwc, int nwb, 
    t_cuPtr *cuPtr){
/*----------------------------*/

    if(nSolutes){
        cudaMalloc((void**) &(cuPtr->Dtd), sizeof(double));
        cudaMalloc((void**) &(cuPtr->dtAux), sizeof(double));

        //solute control arrays
        cudaMalloc((void**) &(cuPtr->typeDiff), nSolutes*sizeof(int));
        cudaMalloc((void**) &(cuPtr->k_xx), nSolutes*sizeof(double));
        cudaMalloc((void**) &(cuPtr->k_yy), nSolutes*sizeof(double));

        //cells
        cudaMalloc((void**) &(cuPtr->hphi), nSolutes*ncells*sizeof(double));
        cudaMalloc((void**) &(cuPtr->phi), nSolutes*ncells*sizeof(double));
        cudaMalloc((void**) &(cuPtr->localDtd), nSolutes*ncells*sizeof(double));
        cudaMalloc((void**) &(cuPtr->BTcell), nSolutes*ncells*sizeof(double));

        //cells*NCwall
        cudaMalloc((void**) &(cuPtr->dhphi), nSolutes*nWallCell*sizeof(double));
        cudaMalloc((void**) &(cuPtr->Bwall), nSolutes*nWallCell*sizeof(double));
    }

    return(1);

}



////////////////////////////////////////////////////////////////
int copySoluteArraysCudaMem(
    int nSolutes,
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

    if(nSolutes){
        cudaMemcpy(&(garrays->flagDiffusion), &(carrays->flagDiffusion), sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(cuPtr->Dtd, &(carrays->Dtd), sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(&(garrays->Dtd), &(carrays->Dtd), sizeof(double), cudaMemcpyHostToDevice);

        //solutes
        cudaMemcpy((cuPtr->typeDiff), (carrays->typeDiff), nSolutes*sizeof(int), cudaMemcpyHostToDevice );
        cudaMemcpy((cuPtr->k_xx), (carrays->k_xx), nSolutes*sizeof(double), cudaMemcpyHostToDevice );
        cudaMemcpy((cuPtr->k_yy), (carrays->k_yy), nSolutes*sizeof(double), cudaMemcpyHostToDevice );        

        //solutes*cells    
        cudaMemcpy((cuPtr->hphi), (carrays->hphi), nSolutes*ncells*sizeof(double), cudaMemcpyHostToDevice );
        cudaMemcpy((cuPtr->phi), (carrays->phi), nSolutes*ncells*sizeof(double), cudaMemcpyHostToDevice );
        cudaMemcpy((cuPtr->localDtd), (carrays->localDtd), nSolutes*ncells*sizeof(double), cudaMemcpyHostToDevice );
        cudaMemcpy((cuPtr->BTcell), (carrays->BTcell), nSolutes*ncells*sizeof(double), cudaMemcpyHostToDevice );
    
        //solutes*cells*NCwall
        cudaMemcpy((cuPtr->dhphi), (carrays->dhphi), nSolutes*nWallCell*sizeof(double), cudaMemcpyHostToDevice );
        cudaMemcpy((cuPtr->Bwall), (carrays->Bwall), nSolutes*nWallCell*sizeof(double), cudaMemcpyHostToDevice );
    }

    return(1);

}  



////////////////////////////////////////////////////////////////
__global__ void assignSoluteArraysToCudaMem(int nSolutes, t_arrays *garrays,
	//------------------------cells
	int *typeDiff,
	double *k_xx,
	double *k_yy,
	double *hphi,
	double *phi,
    double *localDtd,
    double *BTcell,
	double *dhphi,
    double *Bwall){
/*----------------------------*/

    if(nSolutes){    
        //solute controls
        garrays->typeDiff=typeDiff;
        garrays->k_xx=k_xx;
        garrays->k_yy=k_yy;

        //mesh
        garrays->hphi=hphi;
        garrays->phi=phi;
        garrays->localDtd=localDtd;
        garrays->BTcell=BTcell;

        //cells*NCwall
        garrays->dhphi=dhphi;
        garrays->Bwall=Bwall;
    }

}


////////////////////////////////////////////////////////////////
EXPORT_DLL int freeSoluteCudaMemory(
    int nSolutes, 
    t_cuPtr *cuPtr){
/*----------------------------*/

    if(nSolutes){
        
        cudaFree(cuPtr->Dtd);
        cudaFree(cuPtr->dtAux);

        //solute control arrays
        cudaFree(cuPtr->typeDiff);
        cudaFree(cuPtr->k_xx);
        cudaFree(cuPtr->k_yy);

        //cells
        cudaFree(cuPtr->hphi);
        cudaFree(cuPtr->phi);
        cudaFree(cuPtr->localDtd);
        cudaFree(cuPtr->BTcell);

        //cells*NCwall
        cudaFree(cuPtr->dhphi);
        cudaFree(cuPtr->Bwall);
    }

    return(1);

}

#endif
