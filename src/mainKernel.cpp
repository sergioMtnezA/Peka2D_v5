#include "mainKernel.h"

#ifdef _WIN32
	#if APP_MODE == 1
	using namespace System;
	using namespace System::IO;
	using namespace System::Windows::Forms;
	using namespace System::Drawing;
	using namespace dev;
	#endif
#endif

int runMainKernel (int argc, char * argv[]) {

	strcpy(spar.dir,argv[1]);
	strcpy(spar.proj,argv[2]);

	char temp[1024];

    t_arrays *carrays;
    carrays = (t_arrays*) malloc(sizeof(t_arrays));

	t_timers timers;

	clock_t start0=clock();

	loadSimulationDomain(
        pksetup, 
        (&spar), 
        mesh,
		&(timers), 
        msg);

    getchar();    

    initializeComputationArrays(
        spar, 
        mesh, 
        carrays,
		&(timers), 
        msg);

	getchar(); 	

    computeSimulation(
        spar,
        mesh, 
        carrays,
		&(timers),
		msg);


	clock_t end0=clock();
	double dtime=double(end0-start0)/CLOCKS_PER_SEC;
	timers.total=dtime;
	sprintf(temp,"Total simulation time %.6lf",dtime);
	Notify(temp,MSG_L0,msg);

	//computaiton timers file
	if(write_timers(spar.dir, timers, msg)){
		sprintf(temp,"Timers file written");
		Notify(temp,MSG_L1,msg);	
	} 

	getchar();
	//freeMemory();

	return 1;
}



////////////////////////////////////////////////////////////////
int loadSimulationDomain(
	Peka2D_Setup *pksetup, 
	t_parameters *spar, 
	t_mesh *mesh,
	t_timers *timers,
	t_message *msg){
/*----------------------------*/	
	FILE *fp;
	char temp[1024], filename[1024];

	//Initializing exceptions
	msg->error=0;

	msg->nWarning=0;
	msg->warning=(char**)malloc(sizeof(char*));
	msg->warning[0]=(char*)malloc(sizeof(char)*1024);
	
	msg->nMsgL0=0;
	msg->MsgL0 = (char**) malloc(sizeof(char*));
	msg->MsgL0[0]=(char*)malloc(sizeof(char)*1024);

	msg->nMsgL1=0;
	msg->MsgL1 = (char**) malloc(sizeof(char*));
	msg->MsgL1[0]=(char*)malloc(sizeof(char)*1024);
	
	msg->nMsgL2=0;
	msg->MsgL2 = (char**) malloc(sizeof(char*));
	msg->MsgL2[0]=(char*)malloc(sizeof(char)*1024);
	
	msg->nMsgL3=0;
	msg->MsgL3 = (char**) malloc(sizeof(char*));
	msg->MsgL3[0]=(char*)malloc(sizeof(char)*1024);	


	//Create log file
	sprintf(msg->logFile,"%s%s.log",spar->dir,spar->proj);
	fp=fopen(msg->logFile,"w");
	fclose(fp);
	sprintf(temp,"Log file %s created",msg->logFile);
	Notify(temp,MSG_L0,msg);

	//Start loading time .....................................
	clock_t start0=clock();

	//Read run control data
	if(loadControlParameters(pksetup, spar, mesh, msg)){
		sprintf(temp,"Control data loading completed");
		Notify(temp,MSG_L1,msg);		
	}

	//Load mesh data
	if(msg->error){
		return 0;
	}else{
		if(loadMeshData(pksetup, spar, mesh, msg)){
			sprintf(temp,"Mesh data loading completed");
			Notify(temp,MSG_L1,msg);		
		}        
	}

    //Load boundary conditions
	if(msg->error){
		return 0;
	}else{
		 if(loadBoundaryConditions(pksetup, spar, mesh, msg)){
			sprintf(temp,"Mesh data loading completed");
			Notify(temp,MSG_L1,msg);		
		}        
	} 

    // Write mesh file
    if(msg->error){
		return 0;
	}else{
        sprintf(filename,"%s%s.vtk",spar->dir,spar->proj);
        if(write_vtk_mesh(filename,mesh,msg)){
			sprintf(temp,"Mesh VTK file written");
			Notify(temp,MSG_L0,msg);	
        }
    }

	clock_t end0=clock();
	timers->loading += double(end0-start0)/CLOCKS_PER_SEC;
	//sprintf(temp,"Loading time %.6lf",timers->loading);
	//Notify(temp,MSG_L0,msg);
	//End loading time .....................................

    sprintf(temp,"Loading simulation completed");
    Notify(temp,MSG_L1,msg);
	
	return 1;
}


////////////////////////////////////////////////////////////////
int initializeComputationArrays(
	t_parameters spar, 
	t_mesh *mesh,
    t_arrays *carrays,
	t_timers *timers,   
	t_message *msg){
/*----------------------------*/

	FILE *fp;
	char temp[1024], filename[1024];

	//Start initialization time .....................................
	clock_t start0=clock();

    //Create control arrays
	if(msg->error){
		return 0;
	}else{
		if(allocateArraysMemory(spar, mesh, carrays, msg)){
			sprintf(temp,"Arrays memory allocation completed");
			Notify(temp,MSG_L2,msg);		
		}        
	}    

    //Create control arrays
	if(msg->error){
		return 0;
	}else{
		if(initilizeControlArrays(spar, mesh, carrays, msg)){
			sprintf(temp,"Initialize control arrays completed");
			Notify(temp,MSG_L2,msg);		
		}        
	}

    //Create communication arrays
	// if(msg->error){
	// 	return 0;
	// }else{
	// 	if(initilizeCommunicationArrays(spar, mesh, carrays, cuPtr, msg)){
	// 		sprintf(temp,"Initialize communication arrays completed");
	// 		Notify(temp,MSG_L1,msg);		
	// 	}        
	// }     	

    //Create mesh arrays
	if(msg->error){
		return 0;
	}else{
		if(initilizeMeshArrays(spar, mesh, carrays, msg)){
			sprintf(temp,"Initialize mesh arrays completed");
			Notify(temp,MSG_L2,msg);		
		}        
	}     

    //Create boundary arrays
	/*if(msg->error){
		return 0;
	}else{
		 if(initilizeBoundaryArrays(spar, mesh, carrays, msg)){
			sprintf(temp,"Initialize mesh arrays completed");
			Notify(temp,MSG_L1,msg);		
		}        
	}*/  


	clock_t end0=clock();
	timers->init += double(end0-start0)/CLOCKS_PER_SEC;
	//sprintf(temp,"Initializtion time %.6lf",timers->init);
	//Notify(temp,MSG_L0,msg);
	//End initialization time .....................................       

    sprintf(temp,"Computation arrays completed");
    Notify(temp,MSG_L1,msg);

	return 1;
}




