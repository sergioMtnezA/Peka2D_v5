#include "cuTilities.cuh"


// ENDIAN ##########################################################################
static long _TestEndian=1;

EXPORT_DLL int IsLittleEndian(void) {
	return *(char*)&_TestEndian;
}

/******************************************************************************
  FUNCTION: SwapEndian
  PURPOSE: Swap the byte order of a structure
  EXAMPLE: float F=123.456;; SWAP_FLOAT(F);
******************************************************************************/

EXPORT_DLL void *SwapEndian(void* Addr, const int Nb) {
	static char Swapped[16];
	switch (Nb) {
		case 2:	Swapped[0]=*((char*)Addr+1);
				Swapped[1]=*((char*)Addr  );
				break;
		case 3:	// As far as I know, 3 is used only with RGB images
				Swapped[0]=*((char*)Addr+2);
				Swapped[1]=*((char*)Addr+1);
				Swapped[2]=*((char*)Addr  );
				break;
		case 4:	Swapped[0]=*((char*)Addr+3);
				Swapped[1]=*((char*)Addr+2);
				Swapped[2]=*((char*)Addr+1);
				Swapped[3]=*((char*)Addr  );
				break;
		case 8:	Swapped[0]=*((char*)Addr+7);
				Swapped[1]=*((char*)Addr+6);
				Swapped[2]=*((char*)Addr+5);
				Swapped[3]=*((char*)Addr+4);
				Swapped[4]=*((char*)Addr+3);
				Swapped[5]=*((char*)Addr+2);
				Swapped[6]=*((char*)Addr+1);
				Swapped[7]=*((char*)Addr  );
				break;
		case 16:Swapped[0]=*((char*)Addr+15);
				Swapped[1]=*((char*)Addr+14);
				Swapped[2]=*((char*)Addr+13);
				Swapped[3]=*((char*)Addr+12);
				Swapped[4]=*((char*)Addr+11);
				Swapped[5]=*((char*)Addr+10);
				Swapped[6]=*((char*)Addr+9);
				Swapped[7]=*((char*)Addr+8);
				Swapped[8]=*((char*)Addr+7);
				Swapped[9]=*((char*)Addr+6);
				Swapped[10]=*((char*)Addr+5);
				Swapped[11]=*((char*)Addr+4);
				Swapped[12]=*((char*)Addr+3);
				Swapped[13]=*((char*)Addr+2);
				Swapped[14]=*((char*)Addr+1);
				Swapped[15]=*((char*)Addr  );
				break;
	}
	return (void*)Swapped;
}
// ENDIAN ##########################################################################


////////////////////////////////////////////////////////////////
// Función que calcula el número de núcleos CUDA según la arquitectura GPU
EXPORT_DLL int getCoresPerSM(int major, int minor) {
/*----------------------------*/
	switch (major) {
		case 2: // Fermi
			if (minor == 1) return 48;  // GF10x (sm_2.1) = 48 CUDA cores por SM
			else return 32;              // GF100 (sm_2.0) = 32 CUDA cores por SM
		case 3: // Kepler
			return 192;                  // sm_3.x = 192 CUDA cores por SM
		case 5: // Maxwell
			return 128;                  // sm_5.x = 128 CUDA cores por SM
		case 6: // Pascal
			if (minor == 1) return 128;  // sm_6.1 = 128 CUDA cores por SM
			else if (minor == 0) return 64;  // sm_6.0 = 64 CUDA cores por SM
		case 7: // Volta y Turing
			return 64;                   // sm_7.x = 64 CUDA cores por SM
		case 8: // Ampere
			return 128;                  // sm_8.x = 128 CUDA cores por SM
		default:
			printf("Arquitectura GPU desconocida\n");
			return -1;
	}
}


////////////////////////////////////////////////////////////////
EXPORT_DLL void getCudaMemoryState(size_t *free_mem, size_t *total_mem ){
/*----------------------------*/	
    cudaError_t status = cudaMemGetInfo(free_mem, total_mem);
    printf("%s CUDA memory: Allocated %zu MB || Free %zu MB\n",MSGINFO, 
		((*total_mem)-(*free_mem)) / (1024 * 1024),
		(*free_mem) / (1024 * 1024), MSGINFO);
}


////////////////////////////////////////////////////////////////
__global__ void displayCudaDarray(int nTasks, double *array){
/*----------------------------*/
	int i = threadIdx.x+(blockIdx.x*blockDim.x);
	if(i<nTasks){
		printf("idx %d -- value %lf\n",i,array[i]);
	}	
}


////////////////////////////////////////////////////////////////
__global__ void displayCudaIarray(int nTasks, int *array){
/*----------------------------*/
	int i = threadIdx.x+(blockIdx.x*blockDim.x);
	if(i<nTasks){
		printf("idx %d -- value %d\n",i,array[i]);
	}	
}


////////////////////////////////////////////////////////////////
__global__ void displayCudaDscalar(double *scalar){
/*----------------------------*/
	printf("Double scalar %f\n",(*scalar));
}


////////////////////////////////////////////////////////////////
__global__ void displayCudaIscalar(int *scalar){
/*----------------------------*/
	printf("Integer scalar %d\n",(*scalar));
}


////////////////////////////////////////////////////////////////
__global__ void displayCudaArrayParameter(t_arrays *arrays){
/*----------------------------*/
	double d1 = arrays->dt;
	printf("Parameter 1 %f\n",d1);
}



