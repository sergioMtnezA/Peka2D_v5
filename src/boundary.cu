#include "boundary.cuh"


////////////////////////////////////////////////////
__device__ void d_initialize_inlet(t_arrays *arrays, 
	int i, int idb,
	double *qBoundByCell, double *mBoundByCell, double *mInnerByCell,
	double *qInByInlet, double *mInByInlet){
/*----------------------------*/

	int j,k,l;
	int cidx,cidxneig;
	int count;	

	int indexIN;
	int ncells, NCwall;
    int nbc,i0;
	int typebc; 
    int flagInitialize;	

	int nic,ii0;

	int npts, ip0;
	int tidx;
	double qt, hzt;
	double qsign;

	int cidxZmin;
	double zmin;
	double levelControl;
	double HLin, Uin;
	
	double totalDischarge, totalMass;;	

	double aux1,aux2,aux3,aux4,aux5;

    //typeOBC
    //HYD_INFLOW_Q 1
    //HYD_INFLOW_HZ 2
    //HYD_INFLOW_QHZ 3

	//Mesh data
	ncells=arrays->ncells;
    NCwall=arrays->NCwall;	

	//Bound index
	indexIN=(-idb)-1;
    nbc=arrays->nCellsOBC[i];
    i0=arrays->iniIndexOBC[i]; 
	typebc=arrays->typeOBC[i]; 
    flagInitialize=arrays->flagInitializeOBC[i];

	//Inner index
    nic=arrays->nInnerCellsOBC[i];
    ii0=arrays->iniInnerIndexOBC[i]; 

	//Bound geometry
	cidxZmin = arrays->cellZminOBC[i];
	zmin = arrays->z[cidxZmin];

	//Time series index
	npts=arrays->nPointsSeriesOBC[i];
	ip0=arrays->iniIndexSeriesOBC[i];	


	//Limitations
	double maxFroudeBound = 0.9;
	double inletFroude = 0.1;

	if(typebc==HYD_INFLOW_Q){ //Q(t)

		//interpolate time series value
		tidx = d_get_index(arrays->t, npts, ip0, arrays->tSeriesOBC);
		qt = d_interpolate_vector(arrays->t, npts, ip0, tidx, arrays->tSeriesOBC, arrays->qSeriesOBC);
		if(qt<0.0) qt=0.0;

		aux1 = qt/arrays->totalLengthOBC[i];
		aux2 = cbrt(aux1*aux1/inletFroude/inletFroude/_g_);
		if(aux2 > TOL9){

			// mass balance initialize ---------------------------------------------------
			for(j=0;j<nbc;j++){
				cidx=arrays->cidxBound[i0+j];
				qBoundByCell[i0+j] = 0.0;
				mBoundByCell[i0+j] = -(arrays->h[cidx]*arrays->area[cidx]);
			}

			for(j=0;j<nic;j++){
				cidx=arrays->cidxInner[ii0+j];
				mInnerByCell[ii0+j] = -(arrays->h[cidx]*arrays->area[cidx]);
			}						
		
			// set level at cells ---------------------------------------------------------
			levelControl = zmin+aux2;

			HLin=0.0;
			for(j=0;j<nbc;j++){
				cidx=arrays->cidxBound[i0+j];

				//update effective depth
				aux1=levelControl-arrays->z[cidx];
				if(aux1<TOL12) aux1=0.0;

				arrays->h[cidx] = aux1;
				arrays->sqrh[cidx]=sqrt(aux1);

				//effective discharge area
				if(aux1 >= arrays->minh){
					HLin += aux1*arrays->lWallBound[i0+j];
				}				
			}

			// Neighbouring inner cells
			for(j=0;j<nic;j++){
				cidx=arrays->cidxInner[ii0+j];

				//mass balance initialize
				mInnerByCell[ii0+j] = -(arrays->h[cidx]*arrays->area[cidx]);

				//update effective depth
				aux1=levelControl-arrays->z[cidx];
				if(aux1<TOL12) aux1=0.0;

				arrays->h[cidx] = aux1;
				arrays->sqrh[cidx]=sqrt(aux1);

			}

			//Set effective discharge
			if(HLin>TOL12) {
				Uin = qt/HLin;

				for(j=0;j<nbc;j++){
					cidx=arrays->cidxBound[i0+j];

					if(arrays->h[cidx] >= arrays->minh){
						aux1 = Uin*(-arrays->nxWallBound[i0+j]);
						aux2 = Uin*(-arrays->nyWallBound[i0+j]);

						if(fabs(aux1)<TOL14) aux1=0.0;
						if(fabs(aux2)<TOL14) aux2=0.0;

						arrays->hu[cidx] = arrays->h[cidx]*aux1;
						arrays->hv[cidx] = arrays->h[cidx]*aux2;

						arrays->u[cidx] = aux1;
						arrays->v[cidx] = aux2;
						arrays->modulou[cidx] = sqrt(aux1*aux1 + aux2*aux2);
					}else{
						arrays->hu[cidx] = 0.0;
						arrays->hv[cidx] = 0.0;
						arrays->u[cidx] = 0.0;
						arrays->v[cidx] = 0.0;
						arrays->modulou[cidx] = 0.0;						
					}
				}
			}


			// mass balance updating ---------------------------------------------------
			totalDischarge=0.0;
			totalMass=0.0;
			for(j=0;j<nbc;j++){
				cidx=arrays->cidxBound[i0+j];

				qBoundByCell[i0+j] = (arrays->hu[cidx]*arrays->nxWallBound[i] + arrays->hv[cidx]*arrays->nyWallBound[i]) * arrays->lWallBound[i];
				totalDischarge+=qBoundByCell[i0+j];

				mBoundByCell[i0+j] += arrays->h[cidx]*arrays->area[cidx];
				totalMass+=mBoundByCell[i0+j];	
			}

			for(j=0;j<nic;j++){
				cidx=arrays->cidxInner[ii0+j];

				mInnerByCell[ii0+j] += arrays->h[cidx]*arrays->area[cidx];
				totalMass+=mInnerByCell[ii0+j];
			}

			qInByInlet[indexIN]=totalDischarge;
			mInByInlet[indexIN]=totalMass;		

		}

	}

}











////////////////////////////////////////////////////
__global__ void g_update_open_boundary(int nTasks, t_arrays *arrays, 
	double *qBoundByCell, double *mBoundByCell, double *mInnerByCell,
	double *qInByInlet, double *mInByInlet,
	double *qOutByOutlet, double *mOutByOutlet){		
/*----------------------------*/

	int j,cidx;
	int k,l;

    // extern __shared__ double localh[256];
    // extern __shared__ double localhu[256];
    // extern __shared__ double localhv[256];
    // extern __shared__ double localz[256];

    //extern __shared__ double tempDa1[256];
    //extern __shared__ double tempDa2[256];

	extern __shared__ double sharedMem[];

    __shared__ double tempDa1[threadsPerOBC];
    __shared__ double tempDa2[threadsPerOBC];	

	__shared__ int ncells;
	__shared__ int nSolutes;

    __shared__ int idb;
	__shared__ int nbc,i0;
	__shared__ int npts, ip0;
    __shared__ int flagInitialized;
    __shared__ int typebc;
    __shared__ int indexIN, indexOUT;
	
    __shared__ double totalDepth;
    __shared__ double totalLength, totalArea;

    __shared__ double vol_over_hzMin1,area_over_hzMin1;
	__shared__ double hzMin1,hzMin2;

	__shared__ int tidx;
	__shared__ double qt, hzt;
    __shared__ double qcellin, hincrit; 
    __shared__ double inletFroude, outletFroude, maxFroudeBound;
    __shared__ int cidxZmin;
    __shared__ double zmin;
    __shared__ double HZave;

    __shared__ double areaT;
	__shared__ double qInflow, qOutflow;

    __shared__ double totalDischarge, totalMass;

	double levelControl;
	double areasC;
	double dir1,dir2;
    double aveVel;
    double hun;

    double u,v;
    
    double aux1,aux2,aux3,aux4,aux5;

	#if SET_SOLUTE
		//extern __shared__ double localhphi[1024];
		double phit;
	#endif


    int i = threadIdx.x+(blockIdx.x*blockDim.x); 
    int ithread = threadIdx.x; 
    int iblock = blockIdx.x;
    int idx;

    if(ithread==0){ //only master thread per block

		ncells = arrays->ncells;
		nSolutes = arrays->nSolutes;

        //Bound index
        idb=arrays->idBoundOBC[iblock]; //bound ID: (-id) inlet  (+id) outlet
        nbc=arrays->nCellsOBC[iblock];
        i0=arrays->iniIndexOBC[iblock];         
        //printf("iblock %d idb %d nbc %d i0 %d\n",iblock,idb,nbc,i0);

        //Time series index
        npts=arrays->nPointsSeriesOBC[iblock];
        ip0=arrays->iniIndexSeriesOBC[iblock];
        //printf("npts %d ip0 %d\n",npts,ip0);

        flagInitialized=arrays->flagInitializeOBC[iblock];

        cidxZmin=arrays->cellZminOBC[iblock];
        zmin = arrays->z[cidxZmin];  
        
        inletFroude=0.1;
        outletFroude = 0.1;
        maxFroudeBound=0.9;

        typebc=arrays->typeOBC[iblock]; 

        //initialize boun id
        if(idb<0){ 
            indexIN=(-idb)-1;
            indexOUT=-1;
        }
        if(idb>0){ 
            indexIN=-1;
            indexOUT=idb-1;        
        } 	
          
    }
    __syncthreads();


    //ibound index
    idx = i0+ithread;
	double *localh = sharedMem;
	double *localhu = &localh[nbc];
	double *localhv = &localhu[nbc];
	double *localz = &localhv[nbc];		
	#if SET_SOLUTE
	double *phiIn = &localz[nbc];
	#endif


    //Generate shared memory arrays per block 
    if(ithread<nbc){
        cidx=arrays->cidxBound[idx];

        localh[ithread]  = arrays->h[cidx];
        localhu[ithread] = arrays->hu[cidx];
        localhv[ithread] = arrays->hv[cidx];
        localz[ithread] = arrays->z[cidx];
        //printf("i %d iblock %d ithread %d idx %d h %lf\n",i,iblock,ithread,idx,localh[ithread]);

        //initialize mass bound account
        qBoundByCell[idx] = 0.0;
        mBoundByCell[idx] = 0.0; 
    }
    __syncthreads(); 

	#if SET_SOLUTE
	if(ithread<nSolutes){
		phiIn[ithread]=0.0;
	}
	#endif	   



    if(flagInitialized==0){
        if(ithread<nbc){
            tempDa1[ithread]=localh[ithread];
        }else{
            tempDa1[ithread]=0.0; 
        }
        __syncthreads();

        for(int stride=blockDim.x/2; stride>0; stride/=2){
            if(ithread < stride) {
                tempDa1[ithread] += tempDa1[ithread + stride];
            }
            __syncthreads(); // Sincronizar después de cada paso
        }            

        if(ithread==0){ //only master thread per block
            totalDepth=tempDa1[0];

            if(totalDepth > TOL12){
                flagInitialized = 1;
                arrays->flagInitializeOBC[iblock] = 1;
            }
			//printf("idb %d totalDepth %lf flagInitializeOBC %d \n",idb,totalDepth,arrays->flagInitializeOBC[iblock]);     
		}
        __syncthreads();
    }



    // Pre-condicioning step ///////////////////////////////////////////////////
    if(flagInitialized==1){

        // mass balance initialize ---------------------------------------------------
        if(ithread<nbc){
            qBoundByCell[idx] = 0.0;
            mBoundByCell[idx] -= (localh[ithread]*arrays->areaCellBound[idx]);
        }


        //minimum water level ---------------------------------------------------
        if(ithread<nbc){
            if(localh[ithread] >= arrays->minh){
                tempDa1[ithread]=localh[ithread]+localz[ithread];
            }else{
                tempDa1[ithread]=1e6; 
            }
        }else{
            tempDa1[ithread]=1e6; 
        }
        __syncthreads();        

        for(int stride=blockDim.x/2; stride>0; stride/=2){
            if(ithread < stride) {
                tempDa1[ithread] = MIN(tempDa1[ithread],tempDa1[ithread + stride]);
            }
            __syncthreads(); // Sincronizar después de cada paso
        } 

        if(ithread==0){ //only master thread per block
            hzMin1=tempDa1[0];
        }
        __syncthreads(); 
        

        //volume over minimum water level ---------------------------------------------------
        if(ithread<nbc){
            if(localh[ithread] >= arrays->minh){
                if(localz[ithread] <= hzMin1){
                    tempDa1[ithread]=arrays->areaCellBound[idx] * (localh[ithread]+localz[ithread] - hzMin1);
                    tempDa2[ithread]=arrays->areaCellBound[idx];
                }else{
                    tempDa1[ithread]=arrays->areaCellBound[idx] * localh[ithread];
                    tempDa2[ithread]=0.0;
                }
            }else{
                tempDa1[ithread]=0.0;
                tempDa2[ithread]=0.0; 
            }
        }else{
            tempDa1[ithread]=0.0; 
            tempDa2[ithread]=0.0; 
        }
        __syncthreads();        

        for(int stride=blockDim.x/2; stride>0; stride/=2){
            if(ithread < stride) {
                tempDa1[ithread] += tempDa1[ithread + stride];
                tempDa2[ithread] += tempDa2[ithread + stride];
            }
            __syncthreads(); // Sincronizar después de cada paso
        } 

        if(ithread==0){ //only master thread per block
            vol_over_hzMin1=tempDa1[0];
            area_over_hzMin1=tempDa2[0];

            if(area_over_hzMin1 > TOL12){
                hzMin2 = hzMin1 + vol_over_hzMin1/area_over_hzMin1;
            }else{
                hzMin2 = hzMin1;
            }          
        }
        __syncthreads(); 


        //distribibute volume over minimum water level ---------------------------------------------------
        if(vol_over_hzMin1 > TOL12){
            if(ithread<nbc){
                if(localh[ithread] >= arrays->minh){
                    if(localz[ithread] <= hzMin1){
                        localh[ithread]=MAX(hzMin2-localz[ithread],0.0);
                    }else{
                        localh[ithread]=0.0;
                    }
                }
            }    
        }
        __syncthreads(); 


        // flow regime --------------------------------------------------
        if(ithread<nbc){
            if(localh[ithread] >= arrays->minh){
                tempDa1[ithread]=arrays->lWallBound[idx];
                tempDa2[ithread]=arrays->lWallBound[idx]*localh[ithread];
            }else{
                tempDa1[ithread]=0.0;
                tempDa2[ithread]=0.0; 
            }
        }else{
            tempDa1[ithread]=0.0; 
            tempDa2[ithread]=0.0; 
        }
        __syncthreads();        

        for(int stride=blockDim.x/2; stride>0; stride/=2){
            if(ithread < stride) {
                tempDa1[ithread] += tempDa1[ithread + stride];
                tempDa2[ithread] += tempDa2[ithread + stride];
            }
            __syncthreads(); // Sincronizar después de cada paso
        } 

        if(ithread==0){ //only master thread per block
            totalLength=tempDa1[0];
            totalArea=tempDa2[0];       
        }
        __syncthreads();  

    }    
    

    // Application step /////////////////////////////////////////////////////////////////                
    //inlet
    if(idb<0){
        if(flagInitialized==1){

            //typeOBC
            //HYD_INFLOW_Q 1
            //HYD_INFLOW_HZ 2
            //HYD_INFLOW_QHZ 3

            if(typebc==HYD_INFLOW_Q){ //Q(t) 

                //Interpolate time series value
                if(ithread==0){ //Q(t)   
                    tidx = d_get_index(arrays->t, npts, ip0, arrays->tSeriesOBC);
                    qt = d_interpolate_vector(arrays->t, npts, ip0, tidx, arrays->tSeriesOBC, arrays->qSeriesOBC);
                    
					if(qt<0.0) qt=0.0;
					qt *= arrays->blockSectionOBC[iblock];

					//printf("qt %lf\n", qt);
                }
                __syncthreads();

                //distribute discharge by cells
                if(ithread<nbc){
                    // aux1 = sqrt(localhu[ithread]*localhu[ithread] + localhv[ithread]*localhv[ithread]);
					// if(aux1>TOL12){
                    //     localhu[ithread]/=aux1;
                    //     localhv[ithread]/=aux1;
                    // }else{
                    //     localhu[ithread]=0.0;
                    //     localhv[ithread]=0.0;
                    // }

                    //active inlet cells
                    dir1 = localhu[ithread]*arrays->normalXOBC[iblock] + localhv[ithread]*arrays->normalYOBC[iblock];
                    dir2 = localhu[ithread]*(-arrays->nxWallBound[idx]) + localhv[ithread]*(-arrays->nyWallBound[idx]);

                    if((localh[ithread] >= arrays->minh) && (dir1>0.0 && dir2>0.0) ){
                        areasC=arrays->lWallBound[idx]*localh[ithread];
                        tempDa1[ithread]=areasC;

                        // localhu[ithread]*=areasC;
                        // localhv[ithread]*=areasC;
                        localhu[ithread] = areasC*(-arrays->nxWallBound[idx]);
                        localhv[ithread] = areasC*(-arrays->nyWallBound[idx]);
                    }else{
                        tempDa1[ithread]=0.0;
                        localhu[ithread]=0.0;
                        localhv[ithread]=0.0;
                    } 
                }else{
                    tempDa1[ithread]=0.0;
                }
                __syncthreads();        

                for(int stride=blockDim.x/2; stride>0; stride/=2){
                    if(ithread < stride) {
                        tempDa1[ithread] += tempDa1[ithread + stride];
                    }
                    __syncthreads(); // Sincronizar después de cada paso
                } 

                if(ithread==0){ //only master thread per block
                    areaT=tempDa1[0];
                }
                __syncthreads(); 


                if(areaT>TOL12){
                    if(ithread<nbc){
                        aveVel = qt/areaT;
                        if(localh[ithread] >= arrays->minh){	
                            aux1=aveVel/arrays->lWallBound[idx];
                            localhu[ithread]*=aux1;
                            localhv[ithread]*=aux1;
                        }else{
                            localhu[ithread]=0.0;
                            localhv[ithread]=0.0;
                        }
                    }
                }else if(totalArea>TOL12){
                    if(ithread<nbc){
                        aveVel = qt/totalArea;
                        if(localh[ithread] >= arrays->minh){
                            aux1=aveVel*localh[ithread];
                            localhu[ithread] = aux1*arrays->normalXOBC[iblock];
                            localhv[ithread] = aux1*arrays->normalYOBC[iblock];
                        }else{
                            localhu[ithread]=0.0;
                            localhv[ithread]=0.0;
                        }						
                    }
                }else{	
                    if(ithread<nbc){
                        aux1 = qt/arrays->totalLengthOBC[iblock];
                        aux2 = aux1*aux1/(_g_*inletFroude*inletFroude);

                        localh[ithread] = cbrt(aux2);	
                        localhu[ithread] = aux1*arrays->normalXOBC[iblock];
                        localhv[ithread] = aux1*arrays->normalYOBC[iblock];										
                    }		
                }

                //qInflow
                if(ithread<nbc){
                    hun = localhu[ithread]*(-arrays->nxWallBound[idx]) + localhv[ithread]*(-arrays->nyWallBound[idx]);
                    tempDa1[ithread]=hun*arrays->lWallBound[idx]; 
                }else{
                    tempDa1[ithread]=0.0; 
                }
                __syncthreads();        
    
                for(int stride=blockDim.x/2; stride>0; stride/=2){
                    if(ithread < stride) {
                        tempDa1[ithread] += tempDa1[ithread + stride];
                    }
                    __syncthreads(); // Sincronizar después de cada paso
                } 
    
                if(ithread==0){ //only master thread per block
                    qInflow=tempDa1[0];        
                }
                __syncthreads();                 


                if(fabs(qInflow)>0.0){
                    if(ithread<nbc){
                        aux1 = qt/qInflow;
                        localhu[ithread]*=aux1;
                        localhv[ithread]*=aux1;	

						//printf("iblock %d ithread %d hu %lf hv %lf\n",iblock,ithread,localhu[ithread],localhv[ithread]);
                    }						
                }
                __syncthreads();

                //Froude number corrected flow discharge
                if(ithread<nbc){
                    if(localh[ithread] > TOL12){
                        aux1 = localhu[ithread]*localhu[ithread] + localhv[ithread]*localhv[ithread];
                        aux2 = _g_*maxFroudeBound*maxFroudeBound;
                        aux3 = cbrt(aux1/aux2) - localh[ithread]; 
                        if(aux3>0.0){
                            localh[ithread] += aux3;
                        }
                    }
                }
                __syncthreads();                              
        
            }	
		
		}
    }

    //outlet
    if(idb>0){ 
        if(flagInitialized==1){

            //typeOBC
            // HYD_OUTFLOW_GAUGE 11
            // HYD_OUTFLOW_HZ 12
            // HYD_OUTFLOW_FREE 13
            // HYD_OUTFLOW_FR 14
            // HYD_OUTFLOW_NORMAL 15            

            if(typebc==HYD_OUTFLOW_GAUGE){ //q(hz)

                // average surface level --------------------------------------------------
                if(ithread<nbc){
                    if(localh[ithread] >= arrays->minh){
                        tempDa1[ithread]=(localh[ithread]+localz[ithread])*arrays->lWallBound[idx];
                    }else{
                        tempDa1[ithread]=0.0;
                    }
                }else{
                    tempDa1[ithread]=0.0; 
                }
                __syncthreads();        

                for(int stride=blockDim.x/2; stride>0; stride/=2){
                    if(ithread < stride) {
                        tempDa1[ithread] += tempDa1[ithread + stride];
                    }
                    __syncthreads(); // Sincronizar después de cada paso
                } 

                if(ithread==0){ //only master thread per block
                    HZave=tempDa1[0]; 

                    if(totalLength>TOL12){
                        HZave/=totalLength;
                    }else{
                        HZave = arrays->h[cidxZmin]+arrays->z[cidxZmin];
                    }  					                            
                }
                __syncthreads();                  

                //interpolate time series value
                if(ithread==0){ //only master thread per block  
                    if(HZave >= arrays->hzSeriesOBC[ip0]){
                        tidx = d_get_index(HZave, npts, ip0, arrays->hzSeriesOBC);
                        qt = d_interpolate_vector(HZave, npts, ip0, tidx, arrays->hzSeriesOBC, arrays->qSeriesOBC);
                        
						if(qt<0.0) qt=0.0;
						qt *= arrays->blockSectionOBC[iblock];

                    }else{
                        qt=0.0;
                    } 
                }
                __syncthreads(); 


               //distribute discharge by cells
               if(ithread<nbc){
                    aux1 = sqrt(localhu[ithread]*localhu[ithread] + localhv[ithread]*localhv[ithread]);
                    if(aux1>TOL12){
                        localhu[ithread]/=aux1;
                        localhv[ithread]/=aux1;
                    }else{
                        localhu[ithread]=0.0;
                        localhv[ithread]=0.0;
                    }

                    //active outlet cells
                    dir1 = localhu[ithread]*arrays->normalXOBC[iblock] + localhv[ithread]*arrays->normalYOBC[iblock];
                    dir2 = localhu[ithread]*arrays->nxWallBound[idx] + localhv[ithread]*arrays->nyWallBound[idx];

                    if((localh[ithread] >= arrays->minh) && (dir1>0.0 && dir2>0.0) ){
                        areasC=arrays->lWallBound[idx]*localh[ithread];
                        tempDa1[ithread]=areasC;

                        localhu[ithread]*=areasC;
                        localhv[ithread]*=areasC;
                        //localhu[ithread] = areasC*arrays->nxWallBound[idx];
                        //localhv[ithread] = areasC*arrays->nyWallBound[idx];						
                    }else{
                        tempDa1[ithread]=0.0;
                        localhu[ithread]=0.0;
                        localhv[ithread]=0.0;
                    } 
                }else{
                    tempDa1[ithread]=0.0;
                }
                __syncthreads();        

                for(int stride=blockDim.x/2; stride>0; stride/=2){
                    if(ithread < stride) {
                        tempDa1[ithread] += tempDa1[ithread + stride];
                    }
                    __syncthreads(); // Sincronizar después de cada paso
                } 

                if(ithread==0){ //only master thread per block
                    areaT=tempDa1[0];
                }
                __syncthreads(); 
                

                if(ithread<nbc){
                    if(areaT>TOL12){
                        aveVel = qt/areaT;
                    }else{
                        aveVel = 0.0;
                    }
                    if(localh[ithread] >= arrays->minh){	
                        aux1=aveVel/arrays->lWallBound[idx];
                        localhu[ithread]*=aux1;
                        localhv[ithread]*=aux1;
                    }else{
                        localhu[ithread]=0.0;
                        localhv[ithread]=0.0;
                    }
                }

                //qOutflow
                if(ithread<nbc){
                    hun = localhu[ithread]*arrays->nxWallBound[idx] + localhv[ithread]*arrays->nyWallBound[idx];
                    tempDa1[ithread]=hun*arrays->lWallBound[idx]; 
                }else{
                    tempDa1[ithread]=0.0; 
                }
                __syncthreads();        
    
                for(int stride=blockDim.x/2; stride>0; stride/=2){
                    if(ithread < stride) {
                        tempDa1[ithread] += tempDa1[ithread + stride];
                    }
                    __syncthreads(); // Sincronizar después de cada paso
                } 
    
                if(ithread==0){ //only master thread per block
                    qOutflow=tempDa1[0];        
                }
                __syncthreads();                 


                if(qOutflow>0.0){
                    if(ithread<nbc){
                        aux1 = qt/qOutflow;
                        localhu[ithread]*=aux1;
                        localhv[ithread]*=aux1;	
                    }						
                }
                __syncthreads();                


            } //end if(typebc==HYD_INFLOW_Q){   
			else if(typebc==HYD_OUTFLOW_HZ){ //hz(t)

                //interpolate time series value
                if(ithread==0){ //only master thread per block  
                    tidx = d_get_index(arrays->t, npts, ip0, arrays->tSeriesOBC);
                    hzt = d_interpolate_vector(arrays->t, npts, ip0, tidx, arrays->tSeriesOBC, arrays->hzSeriesOBC);
                }
                __syncthreads(); 


               //set level at cells
               if(ithread<nbc){
					levelControl = localh[ithread]+localz[ithread]; 

					if(hzt <= localz[ithread]){
                        localh[ithread] = 0.0; 
					}
					else{
						aux1 = hzt - levelControl;
						localh[ithread] += aux1;
					}

                }
                __syncthreads();        



                //limit inflow Froude number
 				if(ithread<nbc){
					if(localh[ithread] >= arrays->minh){
						hun = localhu[ithread]*arrays->nxWallBound[idx] + localhv[ithread]*arrays->nyWallBound[idx];
						
						if(hun < 0.0){ //inflow
							aux1 = 0.15*localh[ithread]*sqrt(localh[ithread]*_g_); //discharge with Fr=0.15					
							if(hun < -aux1){
								hun = -aux1;
							}
						}

						localhu[ithread]=hun*arrays->nxWallBound[idx];
						localhv[ithread]=hun*arrays->nyWallBound[idx];
					}
				}
                __syncthreads();             


            } //end if(typebc==HYD_OUTFLOW_HZ){  			
            else if(typebc==HYD_OUTFLOW_FREE){ //free
				if(ithread<nbc){
					if(localh[ithread] >= arrays->minh){
						hun = localhu[ithread]*arrays->nxWallBound[idx] + localhv[ithread]*arrays->nyWallBound[idx];
						if(hun>TOL12){ //outflow
							aux1=localh[ithread]*sqrt(localh[ithread]*_g_); //critical discharge					
							if(hun<aux1){ //subcritical outflow
								hun=aux1;
							}
							
							localhu[ithread]=hun*arrays->nxWallBound[idx];
							localhv[ithread]=hun*arrays->nyWallBound[idx];

						}else{ // prevent back flow
							localhu[ithread]=0.0;
							localhv[ithread]=0.0;
						}
					}
				}
                __syncthreads();    

            } //end if(typebc==HYD_OUTFLOW_FREE){			       
				
        }
    }




    // Post-conditioning step ///////////////////////////////////////////////////////////
    if(flagInitialized==1){

        // mass balance updating ---------------------------------------------------
        if(ithread<nbc){
            if(fabs(localhu[ithread])<TOL14) localhu[ithread]=0.0;
            if(fabs(localhv[ithread])<TOL14) localhv[ithread]=0.0;

            aux1=arrays->lWallBound[idx]*(localhu[ithread]*arrays->nxWallBound[idx] + localhv[ithread]*arrays->nyWallBound[idx]);
            if(idb<0) qBoundByCell[idx] = -aux1; //inlet
            else if(idb>0) qBoundByCell[idx] = aux1; //outlet
            tempDa1[ithread] = qBoundByCell[idx];
            
            aux2 = localh[ithread]*arrays->areaCellBound[idx];
            mBoundByCell[idx] += aux2;
            tempDa2[ithread] = mBoundByCell[idx];     
        }else{
            tempDa1[ithread] = 0.0;
            tempDa2[ithread] = 0.0;
        }
        __syncthreads();        

        for(int stride=blockDim.x/2; stride>0; stride/=2){
            if(ithread < stride) {
                tempDa1[ithread] += tempDa1[ithread + stride];
                tempDa2[ithread] += tempDa2[ithread + stride];
            }
            __syncthreads(); // Sincronizar después de cada paso
        } 

        if(ithread==0){ //only master thread per block
            totalDischarge=tempDa1[0];
            totalMass=tempDa2[0];
            
            if(idb<0){ 
                //qInByInlet[indexIN] = totalDischarge;
                //mInByInlet[indexIN] = totalMass; 
				atomicAdd(&qInByInlet[indexIN], totalDischarge);
				atomicAdd(&mInByInlet[indexIN], totalMass);

				//printf("qin %lf\n", totalDischarge);
            }

            if(idb>0){ 
                //qOutByOutlet[indexOUT] = totalDischarge;
                //mOutByOutlet[indexOUT] = totalMass; 
				atomicAdd(&qOutByOutlet[indexOUT], totalDischarge);
				atomicAdd(&mOutByOutlet[indexOUT], totalMass);	

				//printf("qout %lf\n", totalDischarge);		           
            }

        }
        __syncthreads(); 
    
    

        //Restoring cell primitives ---------------------------------------------------
        if(ithread<nbc){
            cidx=arrays->cidxBound[idx];

            arrays->h[cidx] = localh[ithread];
            arrays->sqrh[cidx] = sqrt(localh[ithread]);

            if(localh[ithread] >= arrays->minh){
                arrays->hu[cidx] = localhu[ithread];
                arrays->hv[cidx] = localhv[ithread]; 

                u = localhu[ithread]/localh[ithread];
                v = localhv[ithread]/localh[ithread];
                arrays->u[cidx] = u;
                arrays->v[cidx] = v;
                arrays->modulou[cidx] = sqrt(u*u+v*v);
				
            }else{
                arrays->hu[cidx] = 0.0;
                arrays->hv[cidx] = 0.0;
                arrays->u[cidx] = 0.0;
                arrays->v[cidx] = 0.0;
                arrays->modulou[cidx] = 0.0;

            }        
        }
        __syncthreads(); 

	}


    // Solute inflow update ///////////////////////////////////////////////////////////
	#if SET_SOLUTE
	if(nSolutes){
	if(flagInitialized==1){
		if(idb<0){	
			if(typebc==HYD_INFLOW_Q){

				//interpolate time series value
				if(ithread<nSolutes){
					tidx = d_get_index(arrays->t, npts, ip0, arrays->tSeriesOBC);
					phit = d_interpolate_matrix(arrays->t, ithread, arrays->nTotalPointSeries, 
						npts, ip0, 
						tidx, 
						arrays->tSeriesOBC, arrays->phiSeriesOBC);
					
					if(phit<0.0) phit=0.0;
					phiIn[ithread]=phit;
                    //printf("sol %d tidx %d phit %lf phiIn %lf\n",ithread,tidx,phit,phiIn[ithread]);
				}
				__syncthreads();
				

				//impose solute concentration at cells
				if(ithread<nbc){	
					if(localh[ithread] > TOL12){
						for(j=0;j<nSolutes;j++){				
							arrays->phi[j*ncells+cidx] = phiIn[j];
						}
					}else{
						for(j=0;j<nSolutes;j++){
							arrays->phi[j*ncells+cidx] = 0.0;
						}
					}
				}			
				__syncthreads(); 

			}
		}

		//update solute at cells (always)
		if(ithread<nbc){	
			if(localh[ithread] > TOL12){
				for(j=0;j<nSolutes;j++){				
					arrays->hphi[j*ncells+cidx] = localh[ithread]*arrays->phi[j*ncells+cidx];
				}
			}else{
				for(j=0;j<nSolutes;j++){
					arrays->hphi[j*ncells+cidx] = 0.0;
				}
			}
		}			
		__syncthreads(); 

	}
	}		
	#endif	  



}




////////////////////////////////////////////////////
__device__ int d_get_index(double valor, int n, int idx0, double *x){
/*----------------------------*/

	int i;
	int pos;

	pos=0;
	if(valor >= x[idx0+n-1]){
		pos=idx0+n-1;
	}else{
		for (i=idx0;i<idx0+n-1;i++){
			if((fabs(valor-x[i])<TOL12) || (valor>x[i])&&(valor<x[i+1])){
				pos=i;
			}
		}
		if(fabs(valor-x[idx0+n-2])<TOL12){
			pos=idx0+n-2;
		}
	}
	return pos;

}


////////////////////////////////////////////////////
__device__ double d_interpolate_vector(double valor, 
	int n, int idx0,
	int tidx, 
	double *x, 
	double *y){
/*----------------------------*/
	
	double x0, y0;
	double x1, y1; 
	double data;

	if(tidx==idx0+n-1){
		data=y[tidx];
	}else{
		x0=x[tidx];
		x1=x[tidx+1];
		y0=y[tidx];
		y1=y[tidx+1];
		if(fabs(x1-x0)<=TOL6)
			data = y1;
		else 
			data = y0+(y1-y0)*(valor-x0)/(x1-x0);
	}
	return data;

}



////////////////////////////////////////////////////
__device__ double d_interpolate_matrix(double valor, 
	int line, int nPointsLine,
	int n, int idx0,
	int tidx, 
	double *x, 
	double *y){
/*----------------------------*/
	
	double x0, y0;
	double x1, y1; 
	double data;

	if(tidx==idx0+n-1){
		data=y[line*nPointsLine+tidx];
	}else{
		x0=x[tidx];
		x1=x[tidx+1];
		y0=y[line*nPointsLine+tidx];
		y1=y[line*nPointsLine+tidx+1];
		if(fabs(x1-x0)<=TOL6)
			data = y1;
		else 
			data = y0+(y1-y0)*(valor-x0)/(x1-x0);
	}	

    //printf("ithred %d phit %lf\n",line,data);

	return data;

}

