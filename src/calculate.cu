#include "calculate.cuh"

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
    int nIter;

    clock_t stime1, stime2;
    clock_t start0, end0;

    t_arrays *garrays;
    t_cuPtr cuPtr; 
    int nTasks, blocksPerGrid;

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
    int coresPerSM = getCoresPerSM(prop.major, prop.minor);
    sprintf(temp,"CUDA sm %d || threads %d", prop.multiProcessorCount, prop.multiProcessorCount*coresPerSM);
    Notify(temp,MSG_L0,msg);    


    //Initialize CUBLAS
	cublasCreate(&cuHandle);
	cublasSetPointerMode_v2(cuHandle, CUBLAS_POINTER_MODE_DEVICE);
    sprintf(temp,"CUBLAS methods initialized");
    Notify(temp,MSG_L3,msg);	     


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


    //Initial mass balance
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



    getchar();



	/////////////////////////////////////////////////////////////////////////
	// TIME LOOP
	/////////////////////////////////////////////////////////////////////////
    sprintf(temp,"Simulation starts ...............");
    Notify(temp,MSG_L1,msg);

    //Start computation time .....................................
	start0=clock();

	while(t < carrays->tf){

        //update carrays to next time step
		generateTimeStep(&t, carrays, garrays, &(cuPtr), timers, msg);

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
            //getchar();
        }

        //sprintf(filename,"%sprueba.vtk",spar.dir);
        //if(write_vtk_state(filename,mesh,carrays,msg)){
        //    sprintf(temp,"Prueba VTK file written");
        //    Notify(temp,MSG_L1,msg);	
        //}         

        carrays->nIter++;

        stime2=clock();
        timers->writeOut += double(stime2-stime1)/CLOCKS_PER_SEC;
        //End IO time .....................................          



        //getchar();

	}
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

    // Free CUDA memory
    freeCudaMemory(&(cuPtr));
    cudaFree(garrays);
    cublasDestroy(cuHandle);

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

	int checkpos;
    int ncells=carrays->ncells;
    int nwc=carrays->nw_calc;

    int nTasks, blocksPerGrid;

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


    nTasks=carrays->nActWalls;
    blocksPerGrid = nTasks/threadsPerBlock + 1; 
    g_wall_rotated_calculus <<<blocksPerGrid,threadsPerBlock>>> (nTasks, garrays, cuPtr->localDt);
    cudaMemcpy(&(carrays->nActCells), &(garrays->nActCells), sizeof(int), cudaMemcpyDeviceToHost );


    //g_bound_calculus
    

    nTasks=carrays->nw_calc;
    cublasIdamin(cuHandle, nTasks, cuPtr->localDt, 1, cuPtr->index);
	g_get_dtmin <<<1,1>>>  (garrays, cuPtr->localDt, cuPtr->index);


    nTasks=carrays->nActCells;
    blocksPerGrid = nTasks/threadsPerBlock + 1; 
    g_update_contributions <<<blocksPerGrid,threadsPerBlock>>> (nTasks, garrays);

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


    //g_add_rain_contributions
    //g_add_infiltration_contributions
    //g_add_evaporation_contributions


    nTasks=carrays->nActCells;
    blocksPerGrid = nTasks/threadsPerBlock + 1; 
	g_update_cells <<<blocksPerGrid,threadsPerBlock>>> (nTasks, garrays);
    cudaMemcpy(&(carrays->nActWalls), &(garrays->nActWalls), sizeof(int), cudaMemcpyDeviceToHost );

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

    stime2=clock();
	timers->wetDryFix += double(stime2-stime1)/CLOCKS_PER_SEC;
    //End wetDryFix time .....................................    



    //g_update_boundaries



    //Start memoryTransfer time .....................................
	stime1=clock();

    nTasks=carrays->nActCells;
    blocksPerGrid = nTasks/threadsPerBlock + 1; 
    g_compute_cell_mass <<<blocksPerGrid,threadsPerBlock>>> (nTasks, garrays, cuPtr->mass);    
    cublasDasum(cuHandle, carrays->ncells, cuPtr->mass, 1, cuPtr->massNew);
    cudaMemcpy(&(garrays->massNew), cuPtr->massNew, sizeof(double), cudaMemcpyDeviceToDevice );
    cudaMemcpy(&(carrays->massNew), cuPtr->massNew, sizeof(double), cudaMemcpyDeviceToHost );


    if(carrays->nIter%carrays->nIterOut==0){
        g_compute_mass_error <<<1,1>>> (garrays);
        cudaMemcpy(&(carrays->massError), &(garrays->massError), sizeof(double), cudaMemcpyDeviceToHost );
    }


    if(carrays->dumpComponent || carrays->dumpState){
        // Transfer arrays from GPU to CPU
        //cudaMemcpy((carrays->z), (cuPtr->z), ncells*sizeof(double), cudaMemcpyDeviceToHost );
        cudaMemcpy((carrays->h), (cuPtr->h), ncells*sizeof(double), cudaMemcpyDeviceToHost );
        //cudaMemcpy((carrays->hu), (cuPtr->hu), ncells*sizeof(double), cudaMemcpyDeviceToHost );
        //cudaMemcpy((carrays->hv), (cuPtr->hv), ncells*sizeof(double), cudaMemcpyDeviceToHost );
        cudaMemcpy((carrays->u), (cuPtr->u), ncells*sizeof(double), cudaMemcpyDeviceToHost );
        cudaMemcpy((carrays->v), (cuPtr->v), ncells*sizeof(double), cudaMemcpyDeviceToHost );
        cudaMemcpy((carrays->modulou), (cuPtr->modulou), ncells*sizeof(double), cudaMemcpyDeviceToHost );          
    }

    stime2=clock();
	timers->memoryTransfer += double(stime2-stime1)/CLOCKS_PER_SEC;
    //End memoryTransfer time ..................................... 




    // Reconstruct actCells and actWalls arrays 
    #if RECONSTRUC_ACTIVE
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

    stime2=clock();
	timers->wetDryFix += double(stime2-stime1)/CLOCKS_PER_SEC;
    //End wetDryFix time .....................................  
    #endif


}
