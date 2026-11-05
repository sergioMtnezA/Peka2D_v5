#include "calculate.cuh"

//--- Include preCICE lib
#include <precice/precice.hpp>
#include <cmath>

cublasHandle_t cuHandle;

////////////////////////////////////////////////////////////////
EXPORT_DLL int computeSimulation(
    t_parameters spar,
    t_mesh *mesh, 
    t_arrays *carrays,
    t_timers *timers,  
	t_message *msg){
/*----------------------------*/

	FILE *fp;
	char temp[1024], filename[1024];

    double t;
    int nIter, ncc;

    clock_t stime1, stime2;
    clock_t start0, end0;

    t_arrays *garrays;
    t_cuPtr cuPtr;
    int nTasks, blocksPerGrid;

    // --- for preCICE
    int k;
    std::vector<double> pekaOFbbox;
    ncc = carrays->nCoupPEKACells;
    
    t=carrays->ti;
    nIter=carrays->nIter;
 	//Start GPU initialization time .....................................
	stime1=clock();

    //Select GPU device
    int num_devices, device;
    cudaGetDeviceCount(&num_devices);
    sprintf(temp,"Num. GPU devices: %d",num_devices);
    Notify(temp,MSG_L0,msg);	

    cudaDeviceProp prop;
    if(num_devices > 1){
        int max_multiprocessors = 0, max_device = 0;
        for (device = 0; device < num_devices; device++){
            cudaGetDeviceProperties(&prop, device);
            if (max_multiprocessors < prop.multiProcessorCount) {
                max_multiprocessors = prop.multiProcessorCount;
                max_device = device;
            }
        }
        device=max_device;
    }else{
        device = 0;
    }
    cudaSetDevice(device);
    sprintf(temp,"Selected device ID: %d",device);
    Notify(temp,MSG_L3,msg);

    cudaGetDeviceProperties(&prop, device);
    sprintf(temp,"Compute capability %d.%d",prop.major, prop.minor);
    Notify(temp,MSG_L0,msg);     
    int coresPerSM = getCoresPerSM(prop.major, prop.minor);
    sprintf(temp,"CUDA sm %d || threads %d", prop.multiProcessorCount, prop.multiProcessorCount*coresPerSM);
    Notify(temp,MSG_L0,msg);  
    sprintf(temp,"Shared memory: Maximum %zu KB - Limit %zu KB",(prop.sharedMemPerBlock/1024),(prop.sharedMemPerBlockOptin/1024));
    Notify(temp,MSG_L0,msg);      


    //Initialize CUBLAS
	cublasCreate(&cuHandle);
	cublasSetPointerMode_v2(cuHandle, CUBLAS_POINTER_MODE_DEVICE);
    sprintf(temp,"CUBLAS methods initialized");
    Notify(temp,MSG_L3,msg);	


    //Set HEAP memory size for local malloc
    //size_t heapsize = 64*1024*1024; //64MB 
    //cudaDeviceSetLimit(cudaLimitMallocHeapSize, heapsize);     
    //sprintf(temp,"HEAP memory set");
    //Notify(temp,MSG_L3,msg);

    //create CUDA memory and initialize
    cudaMalloc( (void **) &garrays, sizeof(t_arrays));
    if( createArraysCudaMemory(carrays, garrays, &(cuPtr), msg) ){
        sprintf(temp,"CUDA memory allocated and initialized");
        Notify(temp,MSG_L3,msg);	
    }

	stime2=clock();
	timers->initGPU += double(stime2-stime1)/CLOCKS_PER_SEC;
	//sprintf(temp,"Initializtion time %.6lf",timers->init);
	//Notify(temp,MSG_L0,msg);
	//End GPU initialization time .....................................      




 	//Start IO time .....................................
	stime1=clock();

    //initialize computation files
    if(create_computation_files(spar.dir, msg)){
        sprintf(temp,"Computation files created");
        Notify(temp,MSG_L2,msg);	
    } 


    //initialize boundary conditions
    if( computeInitialBoundaryConditions(carrays, garrays, &(cuPtr), msg) ){
        sprintf(temp,"Initial boundary conditions set");
        Notify(temp,MSG_L2,msg);	
    }     


    //initialize mass balance
    if( computeInitialMassBalance(carrays, garrays, &(cuPtr), msg) ){
        write_massBalance(spar.dir, carrays, msg);
        sprintf(temp,"Initial mass in domain %lf m3",carrays->massNew);
        Notify(temp,MSG_L2,msg);	
    }       


    //write initial condition      
    sprintf(filename,"%sstate%d.vtk",spar.dir,0);
    if(write_vtk_state(filename,mesh,carrays,msg)){
        sprintf(temp,"State%d VTK file written",0);
        Notify(temp,MSG_L0,msg);	
    } 

	stime2=clock();
	timers->init += double(stime2-stime1)/CLOCKS_PER_SEC;
	//sprintf(temp,"Initializtion time %.6lf",timers->init);
	//Notify(temp,MSG_L0,msg);
	//End IO time .....................................   


	/////////////////////////////////////////////////////////////////////////
	// TIME LOOP
	/////////////////////////////////////////////////////////////////////////  
    #if SET_COUPLING
    // --- Creates preCICE participant
    precice::Participant sweParticipant(mesh->participantName, mesh->configFile, 0, 1);

    // --- Coupling interface informations
    std::vector<t_coupInterface> ifaces(mesh->nCoupInterfaces);

    // --- Define one bounding box for each interface
    for (k=0; k<mesh->nCoupInterfaces; k++){
        // Check for ncells == 0
        if (mesh->coupInfo[k].nCells == 0) continue;

        pekaOFbbox = CouplingBbox(mesh, k, 1.0);
        sweParticipant.setMeshAccessRegion(mesh->meshNameReceived[k], pekaOFbbox);
    }

    // --- Initialize preCICE participant
    sweParticipant.initialize();

    // --- OpenFOAM faces, coordinates, buffers and map for each interface
    for (k=0; k<mesh->nCoupInterfaces; k++){
        // Access the coupling interface
        t_coupInterface &ci = ifaces[k];

        // Check if nCells == 0
        if (mesh->coupInfo[k].nCells == 0) continue;

        // Get the number of OpenFOAM faces and mesh dimension
        ci.OFFaces = sweParticipant.getMeshVertexSize(mesh->meshNameReceived[k]);
        ci.meshDim = sweParticipant.getMeshDimensions(mesh->meshNameReceived[k]);

        // Vertex IDs and coordinates
        ci.vertexIDs.resize(ci.OFFaces);
        ci.OFCoord.resize(ci.OFFaces*ci.meshDim);
        sweParticipant.getMeshVertexIDsAndCoordinates(mesh->meshNameReceived[k], ci.vertexIDs, ci.OFCoord);

        // Buffers OpenFOAM -> Peka
        ci.alpha2D.resize(ci.OFFaces);
        ci.velocity2D.resize(ci.OFFaces*ci.meshDim);
        ci.pressure2D.resize(ci.OFFaces);

        // Buffers Peka -> OpenFOAM
        ci.alpha3D.resize(ci.OFFaces);
        ci.velocity3D.resize(ci.OFFaces*ci.meshDim);
        ci.pressure3D.resize(ci.OFFaces);

        // Manual neares-neighbor mapping
        mapOFwithPEKA(mesh, k, ci);

        sprintf(temp, "Interface %d (%s): # %d, %d Peka Cells, %d OpenFOAM Faces, max map dist %.4lf",
            k, mesh->meshNameReceived[k], mesh->coupInfo[k].interfaceIndex,
            mesh->coupInfo[k].nCells, ci.OFFaces, ci.mapMaxDist);
        Notify(temp, MSG_L1, msg);
    }

    sprintf(temp,"Simulation starts ...............");
    Notify(temp,MSG_L1,msg);

    //Start computation time .....................................
    start0=clock();

    while(t < carrays->tf){
        // --- preCICE max dt
        carrays->preciceDt = sweParticipant.getMaxTimeStepSize();
        
        // --- Send the dt to GPU
        cudaMemcpy(&(garrays->preciceDt), &(carrays->preciceDt), sizeof(double), cudaMemcpyHostToDevice);

        //update carrays to next time step
        generateTimeStep(&t, carrays, garrays, &(cuPtr), timers, msg);

        // --- After updating the fields
        nTasks = carrays->nCoupPEKACells;
        blocksPerGrid = nTasks/threadsPerBlock + 1;

        // --- Update hCoupInt, huCoupInt and hvCoupInt with the mesh values
        readWriteCouplingValues <<<blocksPerGrid,threadsPerBlock>>> (nTasks, 'r', garrays);
        updateCouplingArrays(carrays,  &(cuPtr), 'r');

        // --- Read the arrays from preCICE
        for (k=0; k<mesh->nCoupInterfaces; k++){
            t_coupInterface &ci = ifaces[k];
            if(mesh->coupInfo[k].nCells == 0) continue;

            // Alpha
            sweParticipant.readData(
                mesh->meshNameReceived[k],
                mesh->alphaWater2D,
                ci.vertexIDs,
                carrays->dt,
                ci.alpha2D);
            
            // Velocity
            sweParticipant.readData(
                mesh->meshNameReceived[k],
                mesh->velocity2D,
                ci.vertexIDs,
                carrays->dt,
                ci.velocity2D);
            
            // Pressure
            sweParticipant.readData(
                mesh->meshNameReceived[k],
                mesh->pressure2D,
                ci.vertexIDs,
                carrays->dt,
                ci.pressure2D);
        }




/*
        // --- Pass the write arrays to preCICE
        sweParticipant.writeData(mesh->meshName, mesh->alphaWater2D, preCICEcellsIDs,
            precice::span<const double>(carrays->hCoupInt, carrays->nCouplingCells)); // <---- UPDATE THIS
        sweParticipant.writeData(mesh->meshName, mesh->velocity2D, preCICEcellsIDs,
            precice::span<const double>(carrays->huCoupInt, carrays->nCouplingCells)); // <---- UPDATE THIS
        sweParticipant.writeData(mesh->meshName, mesh->pressure2D, preCICEcellsIDs,
            precice::span<const double>(carrays->hvCoupInt, carrays->nCouplingCells)); // <---- UPDATE THIS
*/

        // FUNCTION HERE 3D -> 2D

        //Function performing average values


        // --- Update hCoup, huCoup and hvCoup with the received values
        //readWriteCouplingValues <<<blocksPerGrid,threadsPerBlock>>> (nTasks, 'w', garrays);
        //updateCouplingArrays(carrays,  &(cuPtr), 'w');
        


        
        //Start IO time .....................................
        stime1=clock();

        //print screen info
        if(carrays->nIter%carrays->nIterOut==0){
            dump_screen_info(carrays, msg);
        }

        //write components results
        if(carrays->dumpComponent){ 
            //mass balance file
            write_massBalance(spar.dir, carrays, msg);  

            carrays->indexDump++;	
        }        

        //write output results
        if(carrays->dumpState){
            //VTK file
            sprintf(filename,"%sstate%d.vtk",spar.dir,carrays->indexOut);
            if(write_vtk_state(filename,mesh,carrays,msg)){
                sprintf(temp,"State%d VTK file written",carrays->indexOut);
                Notify(temp,MSG_L1,msg);	
            } 

            carrays->indexOut++;	
        }

        carrays->nIter++;

        stime2=clock();
        timers->writeOut += double(stime2-stime1)/CLOCKS_PER_SEC;
        //End IO time .....................................          


        // sprintf(filename,"%sprueba.vtk",spar.dir);
        // if(write_vtk_state(filename,mesh,carrays,msg)){
        //     sprintf(temp,"Prueba VTK file written");
        //     Notify(temp,MSG_L1,msg);	
        // }  
        // getchar();

        //--- advance preCICE dt
        sweParticipant.advance(carrays->dt);
    }
    // --- preCICE end communication
    sweParticipant.finalize();
    
    #endif
	/////////////////////////////////////////////////////////////////////////
	// END TIME LOOP
	/////////////////////////////////////////////////////////////////////////        
       
	end0=clock();
	timers->computeSim += double(end0-start0)/CLOCKS_PER_SEC;
	//sprintf(temp,"Computation time %.6lf",timers->computeSim);
	//Notify(temp,MSG_L0,msg);
	//End computation time .....................................

    sprintf(temp,"Computation loop completed");
    Notify(temp,MSG_L1,msg);



    //Start closing time .....................................
	start0=clock();

    //write hotstart initialization file      
    sprintf(filename,"%shotstart.out",spar.dir);
    if(write_hotstart_file(filename,carrays,msg)){
		sprintf(temp,"Hotstart initialization file written");
		Notify(temp,MSG_L1,msg);        
    }    

    // Free CUDA memory
    #if SET_SOLUTE
    freeSoluteCudaMemory(carrays->nSolutes,&(cuPtr));
    #endif
    freeBoundaCudaMemory(carrays->nOBC, carrays->nInlet, carrays->nOutlet,
        carrays->nTotalBoundCells, carrays->nTotalInnerCells, 
        carrays->nTotalPointSeries,
        &(cuPtr));
    freeCudaMemory(&(cuPtr));
    cudaFree(garrays);
    cublasDestroy(cuHandle);

    end0=clock();
	timers->closeSim += double(end0-start0)/CLOCKS_PER_SEC;
    //End closing time .....................................

	return 1;
}




////////////////////////////////////////////////////////////////
EXPORT_DLL int computeInitialBoundaryConditions(
    t_arrays *carrays,
    t_arrays *garrays,     
    t_cuPtr *cuPtr,
    t_message *msg){
/*----------------------------*/

    int i;
    int icount, ocount;

    int nTasks, blocksPerGrid;
    int obcPerGrid;
    size_t memPerOBC;    

    char temp[1024];

    //Initialize open boundaries
    carrays->qTotalIn=0.0;
    carrays->qTotalOut=0.0;
    cudaMemset(&(garrays->qTotalOut), 0, sizeof(double));
    cudaMemset(&(garrays->qTotalOut), 0, sizeof(double));

    carrays->mTotalIn=0.0;
    carrays->mTotalOut=0.0;
    cudaMemset(&(garrays->mTotalOut), 0, sizeof(double));
    cudaMemset(&(garrays->mTotalOut), 0, sizeof(double));   

    if(carrays->nOBC){

        nTasks=carrays->nTotalBoundCells;
        obcPerGrid = carrays->nOBC; 
        memPerOBC = 4*carrays->nMaxBoundCells*sizeof(double);
        #if SET_SOLUTE
        memPerOBC += carrays->nSolutes*sizeof(double);
        #endif
        sprintf(temp,"Shared memory: nMaxBoundCells %d - Reserved %zu KB",carrays->nMaxBoundCells,(memPerOBC/1024));
        Notify(temp,MSG_L0,msg);

        cudaFuncSetCacheConfig(g_update_open_boundary, cudaFuncCachePreferShared);
        g_update_open_boundary <<<obcPerGrid,threadsPerOBC,memPerOBC>>> (nTasks, garrays, 
            cuPtr->qBoundByCell, cuPtr->mBoundByCell, cuPtr->mInnerByCell,
            cuPtr->qInByInlet,cuPtr->mInByInlet,
            cuPtr->qOutByOutlet,cuPtr->mOutByOutlet);          


        //-------------------------------------------------------------
        //cublasDasum(cuHandle, carrays->nInlet, cuPtr->mInByInlet, 1, cuPtr->mTotalIn);
        cublasDdot(cuHandle, carrays->nInlet, cuPtr->aux1sByInlet, 1, cuPtr->mInByInlet, 1, cuPtr->mTotalIn);
        cudaMemcpy(&(garrays->mTotalIn), cuPtr->mTotalIn, sizeof(double), cudaMemcpyDeviceToDevice );
        cudaMemcpy(&(carrays->mTotalIn), cuPtr->mTotalIn, sizeof(double), cudaMemcpyDeviceToHost );

        //cublasDasum(cuHandle, carrays->nOutlet, cuPtr->mOutByOutlet, 1, cuPtr->mTotalOut);
        cublasDdot(cuHandle, carrays->nOutlet, cuPtr->aux1sByOutlet, 1, cuPtr->mOutByOutlet, 1, cuPtr->mTotalOut);
        cudaMemcpy(&(garrays->mTotalOut), cuPtr->mTotalOut, sizeof(double), cudaMemcpyDeviceToDevice );
        cudaMemcpy(&(carrays->mTotalOut), cuPtr->mTotalOut, sizeof(double), cudaMemcpyDeviceToHost );

        sprintf(temp,"Initial added mass: Inlet %lf m3 - Outlet %lf m3",carrays->mTotalIn, carrays->mTotalOut);
        Notify(temp,MSG_L2,msg);   

        //-------------------------------------------------------------------
        //cublasDasum(cuHandle, carrays->nInlet, cuPtr->qInByInlet, 1, cuPtr->qTotalIn);
        cublasDdot(cuHandle, carrays->nInlet, cuPtr->aux1sByInlet, 1, cuPtr->qInByInlet, 1, cuPtr->qTotalIn);
        cudaMemcpy(&(garrays->qTotalIn), cuPtr->qTotalIn, sizeof(double), cudaMemcpyDeviceToDevice );
        cudaMemcpy(&(carrays->qTotalIn), cuPtr->qTotalIn, sizeof(double), cudaMemcpyDeviceToHost );

        //cublasDasum(cuHandle, carrays->nOutlet, cuPtr->qOutByOutlet, 1, cuPtr->qTotalOut);
        cublasDdot(cuHandle, carrays->nOutlet, cuPtr->aux1sByOutlet, 1, cuPtr->qOutByOutlet, 1, cuPtr->qTotalOut);        
        cudaMemcpy(&(garrays->qTotalOut), cuPtr->qTotalOut, sizeof(double), cudaMemcpyDeviceToDevice );
        cudaMemcpy(&(carrays->qTotalOut), cuPtr->qTotalOut, sizeof(double), cudaMemcpyDeviceToHost );

        sprintf(temp,"Initial discharge: Inlet %lf m3/s - Outlet %lf m3/s",carrays->qTotalIn, carrays->qTotalOut);
        Notify(temp,MSG_L2,msg);          
    }    

    //Inititalize mass balance
    carrays->massTotalIn=0.0;
    carrays->massTotalOut=0.0;
    cudaMemset(&(garrays->massTotalIn), 0, sizeof(double));
    cudaMemset(&(garrays->massTotalOut), 0, sizeof(double));    

    return 1;
}




////////////////////////////////////////////////////////////////
EXPORT_DLL int computeInitialMassBalance(
    t_arrays *carrays,
    t_arrays *garrays,     
    t_cuPtr *cuPtr,
    t_message *msg){
/*----------------------------*/

    int nTasks, blocksPerGrid;

    nTasks=carrays->nActCells;
    blocksPerGrid = nTasks/threadsPerBlock + 1; 
    g_compute_cell_mass <<<blocksPerGrid,threadsPerBlock>>> (nTasks, garrays, cuPtr->mass);
    
    cublasDasum(cuHandle, carrays->ncells, cuPtr->mass, 1, cuPtr->massNew);
    cudaMemcpy(&(garrays->massNew), cuPtr->massNew, sizeof(double), cudaMemcpyDeviceToDevice );
    cudaMemcpy(&(carrays->massNew), cuPtr->massNew, sizeof(double), cudaMemcpyDeviceToHost );

    return 1;
}






////////////////////////////////////////////////////////////////
EXPORT_DLL void generateTimeStep(
    double *t,
    t_arrays *carrays,
    t_arrays *garrays,     
    t_cuPtr *cuPtr,
    t_timers *timers, 
    t_message *msg){
/*----------------------------*/

    int i;
	int checkpos;
    int ncells=carrays->ncells;
    int nwc=carrays->nw_calc;
    int nSolutes=carrays->nSolutes;
    int nSteps;
    double dtDifR, dtAux;

    int icount, ocount;

    int nTasks, blocksPerGrid;
    int obcPerGrid;
    size_t memPerOBC;

    clock_t stime1, stime2;

    carrays->massOld = carrays->massNew;
    cudaMemcpy(&(garrays->massOld), &(garrays->massNew), sizeof(double), cudaMemcpyDeviceToDevice );

    cudaMemcpy(&(garrays->nIter), &(carrays->nIter), sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy(&(garrays->indexOut), &(carrays->indexOut), sizeof(int), cudaMemcpyHostToDevice );
    cudaMemcpy(&(garrays->indexDump), &(carrays->indexDump), sizeof(int), cudaMemcpyHostToDevice );



    //Start wallCalculus time .....................................
	stime1=clock();

    nTasks=carrays->nWallCell;
    blocksPerGrid = nTasks/threadsPerBlock + 1; 
	g_initialize_delta <<<blocksPerGrid,threadsPerBlock>>> (nTasks, garrays);


    #if SET_SOLUTE  
    nTasks=carrays->nWallCell*carrays->nSolutes;  
    blocksPerGrid = nTasks/threadsPerBlock + 1; 
    g_initialize_solute_delta <<<blocksPerGrid,threadsPerBlock>>> (nTasks, garrays);
    #endif


    nTasks=carrays->nActWalls;
    blocksPerGrid = nTasks/threadsPerBlock + 1; 
    g_wall_rotated_calculus <<<blocksPerGrid,threadsPerBlock>>> (nTasks, garrays, cuPtr->localDt);
    cudaMemcpy(&(carrays->nActCells), &(garrays->nActCells), sizeof(int), cudaMemcpyDeviceToHost );


    #if SET_SOLUTE
    #if SET_SOLUTE_UNROLL==0  //compact 
    nTasks=carrays->nActWalls;
    #elif SET_SOLUTE_UNROLL==1  //unroll    
    nTasks=carrays->nActWalls*carrays->nSolutes;
    #endif
    blocksPerGrid = nTasks/threadsPerBlock + 1; 
    g_wall_solute_calculus <<<blocksPerGrid,threadsPerBlock>>> (nTasks, garrays, cuPtr->localDt);

    if(carrays->nOBC){
        #if SET_SOLUTE_UNROLL==0  //compact 
        nTasks=carrays->nTotalBoundCells;
        #elif SET_SOLUTE_UNROLL==1  //unroll  
        nTasks=carrays->nTotalBoundCells*carrays->nSolutes;
        #endif
        blocksPerGrid = nTasks/threadsPerBlock + 1; 
        g_bound_solute_calculus <<<blocksPerGrid,threadsPerBlock>>> (nTasks, garrays);
    }
    #endif
  
    
    nTasks=carrays->nw_calc;
    cublasIdamin(cuHandle, nTasks, cuPtr->localDt, 1, cuPtr->index);
    // --- preCICE dt is evaluated inside this function
	g_get_dtmin <<<1,1>>>  (garrays, cuPtr->localDt, cuPtr->index);


    nTasks=carrays->nActCells;
    blocksPerGrid = nTasks/threadsPerBlock + 1; 
    g_update_contributions <<<blocksPerGrid,threadsPerBlock>>> (nTasks, garrays);


    #if SET_SOLUTE
    #if SET_SOLUTE_UNROLL==0  //compact 
    nTasks=carrays->nActCells; 
    #elif SET_SOLUTE_UNROLL==1  //unroll  
    nTasks=carrays->nActCells*carrays->nSolutes;
    #endif 
    blocksPerGrid = nTasks/threadsPerBlock + 1; 
    g_update_solute_contributions <<<blocksPerGrid,threadsPerBlock>>> (nTasks, garrays);
    #endif

    // Sincronizar la CPU con la GPU
    cudaDeviceSynchronize();

	stime2=clock();
	timers->wallCalculus += double(stime2-stime1)/CLOCKS_PER_SEC;
    //End wallCalculus time .....................................  



    //Start cellUpdating time .....................................
	stime1=clock();

    checkpos=1;
    while(checkpos){
        cudaMemset(cuPtr->check, 0, sizeof(int)); //set the integer to 0
        //displayCudaIscalar <<<1,1>>> (cuPtr->check);

        nTasks=carrays->nActCells;
        blocksPerGrid = nTasks/threadsPerBlock + 1; 
        g_checkpos_h <<<blocksPerGrid,threadsPerBlock>>> (nTasks, garrays, cuPtr->check);
        cudaMemcpy(&(checkpos), cuPtr->check, sizeof(int), cudaMemcpyDeviceToHost );

        if(checkpos){
            // printf("checkpos %d\n",checkpos);
            g_reduce_dt <<<1,1>>> (garrays);
        }
    }


    g_set_new_dt <<<1,1>>> (garrays);
    cudaMemcpy(&(carrays->t), &(garrays->t), sizeof(double), cudaMemcpyDeviceToHost );
    cudaMemcpy(&(carrays->dt), &(garrays->dt), sizeof(double), cudaMemcpyDeviceToHost );
    cudaMemcpy(&(carrays->dumpComponent), &(garrays->dumpComponent), sizeof(int), cudaMemcpyDeviceToHost );
    cudaMemcpy(&(carrays->dumpState), &(garrays->dumpState), sizeof(int), cudaMemcpyDeviceToHost );
    (*t) = carrays->t;


    nTasks=carrays->nActCells;
    blocksPerGrid = nTasks/threadsPerBlock + 1; 
    // --- h, hu and hv is updated here g_update_cells
	g_update_cells <<<blocksPerGrid,threadsPerBlock>>> (nTasks, garrays);
    cudaMemcpy(&(carrays->nActWalls), &(garrays->nActWalls), sizeof(int), cudaMemcpyDeviceToHost );

    #if SET_SOLUTE
    #if SET_SOLUTE_UNROLL==0  //compact 
    nTasks=carrays->nActCells; 
    #elif SET_SOLUTE_UNROLL==1  //unroll  
    nTasks=carrays->nActCells*carrays->nSolutes;
    #endif     
    blocksPerGrid = nTasks/threadsPerBlock + 1; 
    g_update_solute_cells <<<blocksPerGrid,threadsPerBlock>>> (nTasks, garrays);
    #endif

    // Sincronizar la CPU con la GPU
    cudaDeviceSynchronize();

    stime2=clock();
	timers->cellUpdating += double(stime2-stime1)/CLOCKS_PER_SEC;
    //End cellUpdating time .....................................  


    //Start wetDryFix time .....................................
	stime1=clock();

    nTasks=carrays->nActWalls;
    blocksPerGrid = nTasks/threadsPerBlock + 1; 
	g_check_wetdry <<<blocksPerGrid,threadsPerBlock>>> (nTasks, garrays);
    
    
    nTasks=carrays->nActCells;
    blocksPerGrid = nTasks/threadsPerBlock + 1; 
	g_update_wetdry_cells <<<blocksPerGrid,threadsPerBlock>>> (nTasks, garrays);  

    // Sincronizar la CPU con la GPU
    cudaDeviceSynchronize();

    stime2=clock();
	timers->wetDryFix += double(stime2-stime1)/CLOCKS_PER_SEC;
    //End wetDryFix time .....................................    



    #if SET_SOLUTE
    //Start diffusion time .....................................
    stime1=clock();
    if(carrays->flagDiffusion){

        nTasks=carrays->ncells*carrays->nSolutes;
        blocksPerGrid = nTasks/threadsPerBlock + 1; 
        g_initialize_solute_diffusion_delta <<<blocksPerGrid,threadsPerBlock>>> (nTasks, garrays);

        #if SET_SOLUTE_UNROLL==0  //compact 
        nTasks=carrays->nActWalls;
        #elif SET_SOLUTE_UNROLL==1  //unroll  
        nTasks=carrays->nActWalls*carrays->nSolutes;
        #endif  
        blocksPerGrid = nTasks/threadsPerBlock + 1; 
        g_wall_solute_diffusion_calculus <<<blocksPerGrid,threadsPerBlock>>> (nTasks, garrays);

        #if SET_SOLUTE_UNROLL==0  //compact 
        nTasks=carrays->nActCells;
        #elif SET_SOLUTE_UNROLL==1  //unroll  
        nTasks=carrays->nActCells*carrays->nSolutes;
        #endif  
        blocksPerGrid = nTasks/threadsPerBlock + 1; 
        g_update_solute_diffusion_contributions <<<blocksPerGrid,threadsPerBlock>>> (nTasks, garrays,cuPtr->localDtd);

        nTasks=carrays->ncells*carrays->nSolutes;
        cublasIdamin(cuHandle, nTasks, cuPtr->localDtd, 1, cuPtr->index);
        g_get_solute_diffusion_dtmin <<<1,1>>> (garrays, cuPtr->localDtd, cuPtr->index); //update dtd in arrays 
        cudaMemcpy(&(carrays->Dtd), &(garrays->Dtd), sizeof(double), cudaMemcpyDeviceToHost);
        //printf("dt %lf Dtd %lf\n",carrays->dt, carrays->Dtd);
        
        nSteps=(int)(carrays->dt/carrays->Dtd);
        for(i=0;i<(nSteps+1);i++){
            dtDifR = carrays->dt - i*carrays->Dtd;
            dtAux = MIN(carrays->Dtd , dtDifR);
            cudaMemcpy((cuPtr->dtAux), &(dtAux), sizeof(double), cudaMemcpyHostToDevice );
            //printf("iter %d dtDifR %lf dtAux %lf\n",i, dtDifR, dtAux);

            #if SET_SOLUTE_UNROLL==0  //compact 
            nTasks=carrays->nActCells;
            #elif SET_SOLUTE_UNROLL==1  //unroll  
            nTasks=carrays->nActCells*carrays->nSolutes;
            #endif  
            blocksPerGrid = nTasks/threadsPerBlock + 1;                 
            g_update_solute_diffusion_cells <<<blocksPerGrid,threadsPerBlock>>> (nTasks, garrays, cuPtr->dtAux);            
        }   

    }

    // Sincronizar la CPU con la GPU
    cudaDeviceSynchronize();
    stime2=clock();
	timers->diffusion += double(stime2-stime1)/CLOCKS_PER_SEC;
    //End diffusion time .....................................    
    #endif




    //Start openBoundaries time .....................................   
    stime1=clock();
    if(carrays->nOBC){
        carrays->qTotalIn=0.0;
        carrays->qTotalOut=0.0;
        cudaMemset(&(garrays->qTotalIn), 0, sizeof(double));
        cudaMemset(&(garrays->qTotalOut), 0, sizeof(double));

        //cublasDasum(cuHandle, carrays->nInlet, cuPtr->qInByInlet, 1, cuPtr->qTotalIn);
        cublasDdot(cuHandle, carrays->nInlet, cuPtr->aux1sByInlet, 1, cuPtr->qInByInlet, 1, cuPtr->qTotalIn);
        cudaMemcpy(&(garrays->qTotalIn), cuPtr->qTotalIn, sizeof(double), cudaMemcpyDeviceToDevice );
        cudaMemcpy(&(carrays->qTotalIn), cuPtr->qTotalIn, sizeof(double), cudaMemcpyDeviceToHost );
         
        //cublasDasum(cuHandle, carrays->nOutlet, cuPtr->qOutByOutlet, 1, cuPtr->qTotalOut);
        cublasDdot(cuHandle, carrays->nOutlet, cuPtr->aux1sByOutlet, 1, cuPtr->qOutByOutlet, 1, cuPtr->qTotalOut);
        cudaMemcpy(&(garrays->qTotalOut), cuPtr->qTotalOut, sizeof(double), cudaMemcpyDeviceToDevice );
        cudaMemcpy(&(carrays->qTotalOut), cuPtr->qTotalOut, sizeof(double), cudaMemcpyDeviceToHost );


        //update boundaries
        cudaMemset(cuPtr->qInByInlet, 0, carrays->nInlet*sizeof(double));
        cudaMemset(cuPtr->mInByInlet, 0, carrays->nInlet*sizeof(double));
        cudaMemset(cuPtr->qOutByOutlet, 0, carrays->nOutlet*sizeof(double));
        cudaMemset(cuPtr->mOutByOutlet, 0, carrays->nOutlet*sizeof(double)); 
        
        nTasks=carrays->nTotalBoundCells;
        obcPerGrid = carrays->nOBC; 
        memPerOBC = 4*carrays->nMaxBoundCells*sizeof(double);
        #if SET_SOLUTE
        memPerOBC += carrays->nSolutes*sizeof(double);
        #endif
        cudaFuncSetCacheConfig(g_update_open_boundary, cudaFuncCachePreferShared);
        g_update_open_boundary <<<obcPerGrid,threadsPerOBC,memPerOBC>>> (nTasks, garrays, 
            cuPtr->qBoundByCell, cuPtr->mBoundByCell, cuPtr->mInnerByCell,
            cuPtr->qInByInlet,cuPtr->mInByInlet,
            cuPtr->qOutByOutlet,cuPtr->mOutByOutlet);


        carrays->mTotalIn=0.0;
        carrays->mTotalOut=0.0;
        cudaMemset(&(garrays->mTotalIn), 0, sizeof(double));
        cudaMemset(&(garrays->mTotalOut), 0, sizeof(double));        

        //cublasDasum(cuHandle, carrays->nInlet, cuPtr->mInByInlet, 1, cuPtr->mTotalIn);
        cublasDdot(cuHandle, carrays->nInlet, cuPtr->aux1sByInlet, 1, cuPtr->mInByInlet, 1, cuPtr->mTotalIn);
        cudaMemcpy(&(garrays->mTotalIn), cuPtr->mTotalIn, sizeof(double), cudaMemcpyDeviceToDevice );
        cudaMemcpy(&(carrays->mTotalIn), cuPtr->mTotalIn, sizeof(double), cudaMemcpyDeviceToHost ); 
            

        //cublasDasum(cuHandle, carrays->nOutlet, cuPtr->mOutByOutlet, 1, cuPtr->mTotalOut);
        cublasDdot(cuHandle, carrays->nOutlet, cuPtr->aux1sByOutlet, 1, cuPtr->mOutByOutlet, 1, cuPtr->mTotalOut);
        cudaMemcpy(&(garrays->mTotalOut), cuPtr->mTotalOut, sizeof(double), cudaMemcpyDeviceToDevice );     
        cudaMemcpy(&(carrays->mTotalOut), cuPtr->mTotalOut, sizeof(double), cudaMemcpyDeviceToHost ); 
  
    }

    // Sincronizar la CPU con la GPU
    cudaDeviceSynchronize();
    stime2=clock();
	timers->boundConditon += double(stime2-stime1)/CLOCKS_PER_SEC;   
    //End openBoundaries time .....................................   




    //Start memoryTransfer time .....................................
	stime1=clock();

    nTasks=carrays->nActCells;
    blocksPerGrid = nTasks/threadsPerBlock + 1; 
    g_compute_cell_mass <<<blocksPerGrid,threadsPerBlock>>> (nTasks, garrays, cuPtr->mass);    
    cublasDasum(cuHandle, carrays->ncells, cuPtr->mass, 1, cuPtr->massNew);
    cudaMemcpy(&(garrays->massNew), cuPtr->massNew, sizeof(double), cudaMemcpyDeviceToDevice );
    cudaMemcpy(&(carrays->massNew), cuPtr->massNew, sizeof(double), cudaMemcpyDeviceToHost );

    g_compute_mass_error <<<1,1>>> (garrays);
    cudaMemcpy(&(carrays->massError), &(garrays->massError), sizeof(double), cudaMemcpyDeviceToHost );

    if(carrays->dumpComponent){
        // Transfer massBalance from GPU to CPU 
        cudaMemcpy(&(carrays->qTotalIn), &(garrays->qTotalIn), sizeof(double), cudaMemcpyDeviceToHost );
        cudaMemcpy(&(carrays->qTotalOut), &(garrays->qTotalOut), sizeof(double), cudaMemcpyDeviceToHost );
        cudaMemcpy(&(carrays->mTotalIn), &(garrays->mTotalIn), sizeof(double), cudaMemcpyDeviceToHost );
        cudaMemcpy(&(carrays->mTotalOut), &(garrays->mTotalOut), sizeof(double), cudaMemcpyDeviceToHost );                
        cudaMemcpy(&(carrays->massTotalIn), &(garrays->massTotalIn), sizeof(double), cudaMemcpyDeviceToHost );
        cudaMemcpy(&(carrays->massTotalOut), &(garrays->massTotalOut), sizeof(double), cudaMemcpyDeviceToHost );
    }

    if(carrays->dumpState){
        // Transfer flow arrays from GPU to CPU 
        //cudaMemcpy((carrays->z), (cuPtr->z), ncells*sizeof(double), cudaMemcpyDeviceToHost );
        cudaMemcpy((carrays->h), (cuPtr->h), ncells*sizeof(double), cudaMemcpyDeviceToHost );
        //cudaMemcpy((carrays->hu), (cuPtr->hu), ncells*sizeof(double), cudaMemcpyDeviceToHost );
        //cudaMemcpy((carrays->hv), (cuPtr->hv), ncells*sizeof(double), cudaMemcpyDeviceToHost );
        cudaMemcpy((carrays->u), (cuPtr->u), ncells*sizeof(double), cudaMemcpyDeviceToHost );
        cudaMemcpy((carrays->v), (cuPtr->v), ncells*sizeof(double), cudaMemcpyDeviceToHost );
        cudaMemcpy((carrays->modulou), (cuPtr->modulou), ncells*sizeof(double), cudaMemcpyDeviceToHost );          

        
        #if SET_SOLUTE
        // Transfer solute arrays from GPU to CPU 
        if(carrays->nSolutes){
            cudaMemcpy((carrays->phi), (cuPtr->phi), nSolutes*ncells*sizeof(double), cudaMemcpyDeviceToHost );
        }
        #endif
    }

    // Sincronizar la CPU con la GPU
    cudaDeviceSynchronize();
    stime2=clock();
	timers->memoryTransfer += double(stime2-stime1)/CLOCKS_PER_SEC;
    //End memoryTransfer time ..................................... 




    // Reconstruct actCells and actWalls arrays 
    #if RECONSTRUC_ACTIVE && UPDATE_ACTIVE_ARRAYS
    //Start wetDryFix time .....................................
	stime1=clock();

    //if(carrays->dumpState){
    if(carrays->nIter%nIterArrangeActElem==0){

        cudaMemset(&(garrays->nActCells), 0, sizeof(int)); //set the integer to 0
        cudaMemset((cuPtr->activeC), 0, ncells*sizeof(int)); //set the array to 0
        cudaMemset((cuPtr->actCells), 0xFF, ncells*sizeof(int)); //set the array to -1

        cudaMemset(&(garrays->nActWalls), 0, sizeof(int)); //set the integer to 0
        cudaMemset((cuPtr->activeW), 0, nwc*sizeof(int)); //set the array to 0
        cudaMemset((cuPtr->actWalls), 0xFF, nwc*sizeof(int));  //set the array to -1

        nTasks=carrays->nw_calc;
        blocksPerGrid = nTasks/threadsPerBlock + 1; 
        g_reconstruct_active_elements <<<blocksPerGrid,threadsPerBlock>>> (nTasks, garrays);
        cudaMemcpy(&(carrays->nActWalls), &(garrays->nActWalls), sizeof(int), cudaMemcpyDeviceToHost );
        //cudaMemcpy((carrays->actWalls), (cuPtr->actWalls), nwc*sizeof(int), cudaMemcpyDeviceToHost );
        cudaMemcpy(&(carrays->nActCells), &(garrays->nActCells), sizeof(int), cudaMemcpyDeviceToHost );
        //cudaMemcpy((carrays->actCells), (cuPtr->actCells), ncells*sizeof(int), cudaMemcpyDeviceToHost );
           
    }

    // Sincronizar la CPU con la GPU
    cudaDeviceSynchronize();

    stime2=clock();
	timers->wetDryFix += double(stime2-stime1)/CLOCKS_PER_SEC;
    //End wetDryFix time .....................................  
    #endif


}

#if SET_COUPLING
////////////////////////////////////////////////////////////////
// --- Copy h, hu, hv of the coupling BC
__global__ void readWriteCouplingValues(
    int nTasks,
    char cond,
    t_arrays *arrays){
/*----------------------------*/

    int idx;
    int i = threadIdx.x+(blockIdx.x*blockDim.x);

    if(i<nTasks){
        if (cond == 'r'){
            idx = arrays->couplingPEKAIDCells[i];
            arrays->hCoupInt[i]  = arrays->h[idx];
            arrays->huCoupInt[i] = arrays->hu[idx];
            arrays->hvCoupInt[i] = arrays->hv[idx];
        }
        else if (cond == 'w'){
            idx = arrays->couplingPEKAIDCells[i];
            arrays->h[idx] = arrays->hCoupExt[i];
            arrays->hu[idx] = arrays->huCoupExt[i];
            arrays->hv[idx] = arrays->hvCoupExt[i];
        }
    }
}

////////////////////////////////////////////////////////////////
// --- Cells coordinates at the preCICE Coupling BC
EXPORT_DLL std::vector<double> getCellsCoordinatesAtBC(
    t_mesh *mesh,
    int k){
/*----------------------------*/
    
    int j;
    int cellID;
    int iOut;
    std::vector<double> coords;
    iOut = mesh->coupInfo[k].interfaceIndex;

    for(j=0; j<mesh->coupInfo[k].nCells; j++){
        cellID = mesh->cibc->cellsIDs[k].cellsIDIn2D[j];
        coords.push_back(mesh->g_cells->cells[cellID].center[0]); // x
        coords.push_back(mesh->g_cells->cells[cellID].center[1]); // y
    }

    return coords;
}

////////////////////////////////////////////////////////////////
// --- update the arrays in CPU <-> GPU
EXPORT_DLL void updateCouplingArrays(
    t_arrays *carrays,   
    t_cuPtr *cuPtr,
    char cond){
/*----------------------------*/
    
    int ncc;

    ncc = carrays->nCoupPEKACells;

    if (cond == 'r'){
        cudaMemcpy(carrays->hCoupInt, cuPtr->hCoupInt, ncc*sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(carrays->huCoupInt, cuPtr->huCoupInt, ncc*sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(carrays->hvCoupInt, cuPtr->hvCoupInt, ncc*sizeof(double), cudaMemcpyDeviceToHost);
    }
    else if (cond == 'w'){
        cudaMemcpy(cuPtr->hCoupExt, carrays->hCoupExt, ncc*sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(cuPtr->huCoupExt, carrays->huCoupExt, ncc*sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(cuPtr->hvCoupExt, carrays->hvCoupExt, ncc*sizeof(double), cudaMemcpyHostToDevice);
    }
}


// --- Create bounding boxes for preCICE direct mesh access
// The margin should be sufficent to include the OpenFOAM patches
// Returns an empty vector if no Cells
EXPORT_DLL std::vector<double> CouplingBbox(
    t_mesh *mesh,
    int k,
    double margin){

    int j;
    int nCells;

    nCells = mesh->coupInfo[k].nCells;
    
    // If nCells == 0 return empty vector
    if (nCells == 0) {
        return std::vector<double>();
    }

    // Get the coordinates
    std::vector<double> coords2D = getCellsCoordinatesAtBC(mesh, k);   // x,y
    
    // Create a bounding box
    std::vector<double> bbox = { 1e30, -1e30, 1e30, -1e30, -1e3, 1e3 };
   
    // Update the limits of the boundingbox
    for (j = 0; j < nCells; j++) {
        bbox[0] = fmin(bbox[0], coords2D[2*j]);   bbox[1] = fmax(bbox[1], coords2D[2*j]);
        bbox[2] = fmin(bbox[2], coords2D[2*j+1]); bbox[3] = fmax(bbox[3], coords2D[2*j+1]);
    }

    // Add and remove the margin value
    bbox[0] -= margin; bbox[1] += margin;
    bbox[2] -= margin; bbox[3] += margin;

    return bbox;
}

// --- Function to map the nearest PEKA coupling cell for each OF face in x-y plane
EXPORT_DLL void mapOFwithPEKA(
    t_mesh *mesh,
    int k,
    t_coupInterface &ci){
    
    int i, j;
    int offset;
    int nCells;
    double dx, dy, dist, distMin;
    double xOFCell, yOFCell, xPKCell, yPKCell;
    
    // Offset in array position from one to another interface.
    offset = mesh->coupInfo[k].offset;

    // Number of cells of the interface
    nCells = mesh->coupInfo[k].nCells;
    
    // Peka Coordinates
    std::vector<double> pcc = getCellsCoordinatesAtBC(mesh, k);
    
    // Resize the std::vector and assign -1 values
    ci.map.assign(ci.OFFaces, -1);

    // Assign the maxDist as 0
    ci.mapMaxDist = 0.0;

    // --- For every OF Cell
    for (i=0; i<ci.OFFaces; i++){
        distMin = 1e30;

        // --- OF cell coordinates
        xOFCell = ci.OFCoord[(i*3)];
        yOFCell = ci.OFCoord[(i*3)+1];
        
        // --- For every Peka Cell
        for (j=0; j<nCells; j++){

            // --- Peka cell coordinate
            xPKCell = pcc[(j*2)];
            yPKCell = pcc[(j*2)+1];
            dx = (xPKCell - xOFCell);
            dy = (yPKCell - yOFCell);
            dist = sqrt(dx*dx + dy*dy);

            // --- Map the index of Peka cell 
            if (dist < distMin){
                distMin = dist;
                ci.map[i] = offset + j;
            }
        }
        ci.mapMaxDist = fmax(ci.mapMaxDist, distMin);
    }
};

#endif
