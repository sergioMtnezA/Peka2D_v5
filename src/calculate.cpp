#include "calculate.h"

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

    t=carrays->ti;
    nIter=carrays->nIter;


 	//Start IO time .....................................
	stime1=clock();

    //initialize computation files
    if(create_computation_files(spar.dir, msg)){
        sprintf(temp,"Computation files created");
        Notify(temp,MSG_L2,msg);	
    } 


    //Initial mass balance
    if( computeInitialMassBalance(carrays, msg) ){
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
		generateTimeStep(&t, carrays, timers, msg);

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

	return 1;
}





////////////////////////////////////////////////////////////////
EXPORT_DLL int computeInitialMassBalance(
    t_arrays *carrays,
    t_message *msg){
/*----------------------------*/

    int nTasks;

    nTasks=carrays->nActCells;
    c_compute_cell_mass(nTasks, carrays);
    
    return 1;
}






////////////////////////////////////////////////////////////////
EXPORT_DLL void generateTimeStep(
    double *t,
    t_arrays *carrays,
    t_timers *timers, 
    t_message *msg){
/*----------------------------*/

	int checkpos;
    int ncells=carrays->ncells;
    int nwc=carrays->nw_calc;

    int nTasks;

    clock_t stime1, stime2;

    carrays->massOld = carrays->massNew;

    //Start wallCalculus time .....................................
	stime1=clock();

    nTasks=carrays->nWallCell;
	c_initialize_delta(nTasks, carrays);


    nTasks=carrays->nActWalls;
    c_wall_rotated_calculus(nTasks, carrays);


    //c_bound_calculus

    nTasks=carrays->nw_calc;
	c_get_dtmin(nTasks, carrays);


    nTasks=carrays->nActCells;
    c_update_contributions(nTasks, carrays);

	stime2=clock();
	timers->wallCalculus += double(stime2-stime1)/CLOCKS_PER_SEC;
    //End wallCalculus time .....................................  




    //Start cellUpdating time .....................................
	stime1=clock();

    checkpos=1;
    while (checkpos) {
        checkpos=0;

        nTasks=carrays->nActCells;
        c_checkpos_h(nTasks, carrays, &(checkpos));

        if(checkpos){
            // printf("checkpos %d\n",checkpos);
            c_reduce_dt(carrays);
        }
    }


    c_set_new_dt(carrays);
    (*t) = carrays->t;


    //c_add_rain_contributions
    //c_add_infiltration_contributions
    //c_add_evaporation_contributions


    nTasks=carrays->nActCells;
	c_update_cells(nTasks, carrays);

    stime2=clock();
	timers->cellUpdating += double(stime2-stime1)/CLOCKS_PER_SEC;
    //End cellUpdating time .....................................  




    //Start wetDryFix time .....................................
	stime1=clock();

    nTasks=carrays->nActWalls;
	c_check_wetdry(nTasks, carrays);
    
    
    nTasks=carrays->nActCells;
	c_update_wetdry_cells(nTasks, carrays);  

    stime2=clock();
	timers->wetDryFix += double(stime2-stime1)/CLOCKS_PER_SEC;
    //End wetDryFix time .....................................    



    //c_update_boundaries




    //Start memoryTransfer time .....................................
	stime1=clock();

    nTasks=carrays->nActCells;
    c_compute_cell_mass(nTasks, carrays);

    if(carrays->nIter%carrays->nIterOut==0){
        c_compute_mass_error(carrays);
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

        carrays->nActCells=0;
        memset((carrays->activeC), 0, ncells*sizeof(int)); //set the array to 0
        memset((carrays->actCells), 0xFF, ncells*sizeof(int)); //set the array to -1

        carrays->nActWalls=0;        
        memset((carrays->activeW), 0, nwc*sizeof(int)); //set the array to 0
        memset((carrays->actWalls), 0xFF, nwc*sizeof(int));  //set the array to -1

        nTasks=carrays->nw_calc;
        c_reconstruct_active_elements(nTasks, carrays);
           
    }

    stime2=clock();
	timers->wetDryFix += double(stime2-stime1)/CLOCKS_PER_SEC;
    //End wetDryFix time .....................................  
    #endif


}

