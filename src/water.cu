#include "water.cuh"


////////////////////////////////////////////////////
__global__ void g_compute_cell_mass(int nTasks, t_arrays *arrays, double *mass){
/*----------------------------*/
	int idx;
    int i = threadIdx.x+(blockIdx.x*blockDim.x);    
    if(i<nTasks){
        idx=arrays->actCells[i];
        mass[idx] = arrays->h[idx] * arrays->area[idx];
        //printf("cell %d mass %lf\n",idx,mass[idx]);
    }
}


////////////////////////////////////////////////////
__global__ void g_initialize_delta(int nTasks, t_arrays *arrays){
/*----------------------------*/
    int i = threadIdx.x+(blockIdx.x*blockDim.x);    
    if(i<nTasks){
        arrays->dh[i]=0.0;
        arrays->dhu[i]=0.0;
        arrays->dhv[i]=0.0;
        arrays->solidWallByCell[i]=0;
    }
}


////////////////////////////////////////////////////
__global__ void g_wall_rotated_calculus(int nTasks, t_arrays *arrays, double *localDt){
/*----------------------------*/

    int idx;
	int j,k;
    int id1,id2;
    int iw1,iw2; 

	double hL,hR;
	double zL,zR;
	double sqrhL,sqrhR;
    double huL,huR,hvL,hvR;
    double modUL, modUR;
	double unL, unR, vtL, vtR;
	double hbar, unbar,vtbar,cbar,inverseCbar;
	double nx,ny;
	double deltah,deltahu,deltahv;

    double aux1,aux2,aux3,aux4;

	double gp;
 	double deltaXl;	    
	double sqrghL,sqrghR;
	double landaL[3],landaR[3];
	double coc;
    double difqx,difqy;

    double beta_z;
	double dZ,dL,dR,dZs;
    double bedSlope_term,ps1,ps2,psmax;
	double hqi0,hqi1,hqi2;
	
	double beta_f;
	double maxh, minU, manning2;
	double modUbar, nux, nuy, dint, tau, stress, auxTau;
	double friction_term, friction_t;

	double hstar,hls,hrs;   		
	
	double eigel  [3][3]; // l de local
	double landl  [3]; 	
	double landEl [3]; 	
	double alfal  [3]; 
	double betal  [3]; 	 	
	double dUL    [3], dULrot  [3]; 	
	double dUR    [3], dURrot  [3];

	int    solidWall,wetWall;
    double qnormalL;

    double dtl;	

    int i = threadIdx.x+(blockIdx.x*blockDim.x);    
    if(i<nTasks){
        idx=arrays->actWalls[i];

        solidWall=0;
        wetWall = 1;	
        dtl = 1e6;
        qnormalL=0.0;

        nx = arrays->normalX[idx];
        ny = arrays->normalY[idx];		
        deltaXl = arrays->deltaX[idx];	
        gp = arrays->gp[idx];
        
        //cells index
        id1=arrays->idx1[idx];
        id2=arrays->idx2[idx];

        zL  = arrays->z[id1];
        zR  = arrays->z[id2];	        

        hL = arrays->h[id1];
        hR = arrays->h[id2];
		
	    //if ((hL>TOL12)||(hR>TOL12)){

	  
        if(hL<TOL12){ wetWall=0; }				
        if(hR<TOL12){ wetWall=0; }	

        sqrhL=arrays->sqrh[id1];
        sqrhR=arrays->sqrh[id2];

        huL = arrays->hu[id1];
        huR = arrays->hu[id2];
        
        hvL = arrays->hv[id1];
        hvR = arrays->hv[id2];	

        modUL  = arrays->modulou[id1];
        modUR  = arrays->modulou[id2];	                        	


        //---------------- Rotated reference system   	
        unL  = arrays->u[id1]*nx + arrays->v[id1]*ny;
        unR  = arrays->u[id2]*nx + arrays->v[id2]*ny;
				
        vtL  = -arrays->u[id1]*ny + arrays->v[id1]*nx;
        vtR  = -arrays->u[id2]*ny + arrays->v[id2]*nx;

        //ccccccccccccccccccccccccccccccccccccccccccccccccc Wall-averaged values		 
        hbar = 0.5*(hL+hR);
        cbar = sqrt(gp*hbar);
        inverseCbar=1./(cbar);

        aux1 = sqrhL + sqrhR;
        unbar = (unL*sqrhL + unR*sqrhR)/aux1;
        vtbar = (vtL*sqrhL + vtR*sqrhR)/aux1;

        if(fabs(unbar) < TOL12) unbar = 0.0;
        if(fabs(vtbar) < TOL12) vtbar = 0.0;


		//ccccccccccccccccccccccccccccccccccccccccccccccccc Valores propios
        landl[0]= unbar - cbar;
        landl[1]= unbar;
        landl[2]= unbar + cbar;



		//ccccccccccccccccccccccccccccccccccccccccccccccccc vectores propios		    
        eigel[0][0]=1.0;
        eigel[0][1]=landl[0];
        eigel[0][2]=vtbar;

        eigel[1][0]= 0.0;
        eigel[1][1]= 0.0;
        eigel[1][2]= cbar;

        eigel[2][0]=1.0;
        eigel[2][1]=landl[2];
        eigel[2][2]=vtbar;		


		//ccccccccccccccccccccccccccccccccccccccccccccccccc entropy correction
        sqrghL = sqrhL*sqrt(gp); 
        sqrghR = sqrhR*sqrt(gp);  
     
        landaL[0]=unL-sqrghL;
        landaR[0]=unR-sqrghR;

        landaL[1]=unL;
        landaR[1]=unR;
    
        landaL[2]=unL+sqrghL;
        landaR[2]=unR+sqrghR;

        landEl[0]=0.0;
        landEl[1]=0.0;
        landEl[2]=0.0;

        if( (wetWall==1) &&  //wet-wet wall
            (landl[0]*landl[2]<0.0) ){ //subcritical wall
            if( (landaL[0] < 0.0) && (landaR[0] > 0.0) ){
                coc=landaL[0]*(landaR[0]-landl[0])/(landaR[0]-landaL[0]);
                if(coc <= 0.0){
                    landEl[0]=landl[0] - coc;
                    landl [0]=coc;                    
                }

            }else if( (landaL[2] < 0.0) && (landaR[2] > 0.0) ){
                coc=landaR[2]*(landl[2]-landaL[2])/(landaR[2]-landaL[2]);
                if(coc >= 0.0){
                    landEl[2]=landl[2] - coc;
                    landl [2]=coc;
                }
            }               
        }
     	


		//cccccccccccccccccccccccccccccccccccccccccccccc alpha
        deltah = hR - hL ; 
        deltahu = hR*unR - hL*unL;
        deltahv = hR*vtR - hL*vtL;
		 
        aux1   = 0.5*deltah;
				
        difqx  = deltahu - unbar*deltah;
        difqy  = deltahv - vtbar*deltah;
        aux2   = 0.5*inverseCbar*difqx;
				
        alfal[0]= aux1 - aux2;
        alfal[1]= inverseCbar*difqy;	
        alfal[2]= aux1 + aux2;



        //ccccccccccccccccccccccccccccccccccccccccccccccc bed slope
        bedSlope_term = 0.0;
        dR    = zR + hR;
        dL    = zL + hL;       
        dZs   = dR - dL;

        dZ    = zR - zL;
        if(fabs(dZ)>TOL12){
            ps1 = -cbar*cbar * dZ;
            if(dZ<0.0){
                if(dR<zL){
                    aux1 = -hR;
                }else{
                    aux1 =  dZ;
                }
                ps2 = -gp*(hR-0.5*fabs(aux1))*aux1;
            }else{
                if(dL<zR){
                    aux1 =  hL;
                }else{
                    aux1 =  dZ;
                }
                ps2 = -gp*(hL-0.5*fabs(aux1))*aux1;
            }

            bedSlope_term = ps2;

            ////////////// Correccion /////////////////////
            if(fabs(ps1)>fabs(ps2)){
                psmax = ps1;
            }else{
                psmax = ps2;
            }

            if(dZs*dZ > 0.0){ // Free-surface y bed con las misma pendiente
                if(unbar*dZ > 0.0){ // Velocidad en contra de pendiente
                    bedSlope_term = psmax;	
                }
            }
            ////////////////////////////////////////////////
        }
        beta_z = - 0.5*inverseCbar*bedSlope_term; //beta1_z

        //inner state of the normal discharge -- homogeneous + bed-pressure
        hqi0=0.0;
        hqi1=0.0;
        if(landl[0]*landl[2] < 0.0){ // Subcritico + entropy fix
            if(landl[1] >= 0.0){
                hqi0 = hL*unL + alfal[0]*eigel[0][1];
                hqi1 = hqi0 - beta_z;
            }else{
                hqi0 = hR*unR - alfal[2]*eigel[2][1];
                hqi1 = hqi0 - beta_z; //+beta3_z=-beta1_z
            }
		}else if(landl[0] > 0.0){ // Supercritico derecha
            hqi0 = hR*unR - alfal[2]*eigel[2][1];
            hqi1 = hqi0 - beta_z; //+beta3_z=-beta1_z      
		}else if(landl[2] < 0.0){  // Supercritico izquierda
            hqi0 = hL*unL + alfal[0]*eigel[0][1];
            hqi1 = hqi0 - beta_z;
		}    
        if(fabs(hqi1)<TOL12){hqi1 = 0.0;}




        //ccccccccccccccccccccccccccccccccccccccccccccccccccccc friction_slope
        friction_term=0.0;

        if(wetWall==1){
            maxh=fmax(hL,hR);
            minU=fmin(modUR,modUL);
            manning2 = arrays->nman2wall[idx];

            modUbar=0.0;
            nux=0.0;
            nuy=0.0;
            if(modUL>TOL12 && modUR>TOL12){
                modUbar=0.5*(modUL+modUR);
                nux=0.5*(arrays->u[id1]/modUL + arrays->u[id2]/modUR);
                nuy=0.5*(arrays->v[id1]/modUL + arrays->v[id2]/modUR);
            }else if(modUL>TOL12){
                modUbar=modUL;
                nux=arrays->u[id1]/modUL;
                nuy=arrays->v[id1]/modUL;
            }else if(modUR>TOL12){
                modUbar=modUL;
                nux=arrays->u[id2]/modUR;
                nuy=arrays->v[id2]/modUR;
            }

            tau = 0.0;
            stress = 0.0;
            dint = 0.0;
            if(fabs(hqi1) > 0.0){

                //tau  = _rhow_ * gp * manning2 * modUbar * modUbar / cbrt(hbar);
                tau  = _rhow_ * gp * hbar * manning2 * modUbar * modUbar / (maxh*cbrt(maxh));

                aux1 = hqi1 / fabs(hqi1);
                stress = aux1 * tau;
                if(fabs(modUbar) > TOL3) {
                    dint = fabs( (arrays->distCentX[idx])*nux + (arrays->distCentY[idx])*nuy );
                } else {
                    dint = arrays->distNormal[idx];
                }

            }else{
                stress = 0.0;
                dint = arrays->distNormal[idx];
            }

            //ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
            // Integrated friction force
            friction_term = -1.* stress/_rhow_ * dint;
            //ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo	

            //*****************************************
            // turbulent friction energy limitation
            friction_t = fabs(friction_term)/(cbar*cbar);

            //Limite = fabs(un)*modUbar/(2g);
            aux1 = 10. * 0.5*fabs(unbar)*modUbar/gp;
            aux2 = fabs(friction_t);
            if (aux2>aux1){
                friction_term *= (aux1/aux2);
            }
            //**************************************

        }//end mojado==1
        beta_f = - 0.5*inverseCbar*friction_term; //beta1_f


        //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc friction fix
        hqi2=0.0;
        if(landl[0]*landl[2] < 0.0){ // Subcritico + entropy fix
            if(landl[1] >= 0.0){
                hqi2 = hqi1 - beta_f;
            }else{
                hqi2 = hqi1 - beta_f; //+beta3_f=-beta1_f 
            }
		}else if(landl[0] > 0.0){ // Supercritico derecha
            hqi2 = hqi1 - beta_f; //+beta3_f=-beta1_f         
		}else if(landl[2] < 0.0){  // Supercritico izquierda
            hqi2 = hqi1 - beta_f;
		} 
        if(fabs(hqi2)<TOL12){hqi2 = 0.0;}

        aux1 = hqi1*hqi2;
        if(aux1<=0.0){
            beta_f = hqi1;
        }

        betal[0] = beta_z + beta_f;
        betal[1] = 0.0;
        betal[2] = -betal[0];	 


        //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc positivity fix
        if( (wetWall == 1) && //wet-wet wall
            (landl[0]*landl[2] < 0.0) ){ // subcritical + entropy

            hls = hL + alfal[0] + (landEl[2]/landl[0])*alfal[2]; // homogeneous left h-inner state
            if(fabs(hls)<TOL12){hls=0.0;}

            hrs = hR - alfal[2] - (landEl[0]/landl[2])*alfal[0]; // homogenenous right h-inner state
            if(fabs(hrs)<TOL12){hrs=0.0;}

            // alternative 1
            // aux1 = hls - betal[0]/landl[0]; // augmented left h-inner state
            // if(hls>0.0 && aux1<0.0){
            //     betal[0]=  hls*landl[0];
            //     betal[2]= -betal[0];
            // }
            // aux1 = hrs + betal[2]/landl[2]; // augmented right h-inner state
            // if(hrs>0.0 && aux1<0.0){
            //     betal[0]=  hrs*landl[2];
            //     betal[2]= -betal[0];
            // }

            // alternative 2
            aux1 = landl[0]*hls; //minimum for beta1
            aux2 = landl[2]*hrs; //maximum for beta1
            if(aux1<aux2){
                if(betal[0]<aux1){
                    betal[0] = aux1;
                }else if(betal[0]>aux2){
                    betal[0] = aux2;
                }
                betal[2] = -betal[0];
            }
        }


        //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc time step
        //if ((hL>arrays->minh) && (hR>arrays->minh) ){
            dtl=fmin( dtl , deltaXl/fabs(landl[0]) );
            dtl=fmin( dtl , deltaXl/fabs(landl[2]) );
        //}


        //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc compute contributions	
        #pragma unroll
        for(j=0;j<3;j++){
            dULrot[j] = 0.0;
            dURrot[j] = 0.0;            
        }

        hls = hL + alfal[0] + (landEl[2]/landl[0])*alfal[2] - betal[0]/landl[0]; //left h-inner state
        if(fabs(hls)<TOL12){hls=0.0;}
        hrs = hR - alfal[2] - (landEl[0]/landl[2])*alfal[0] + betal[2]/landl[2]; //right h-inner state
        if(fabs(hrs)<TOL12){hrs=0.0;}        
				
        if(hR<TOL12 && hrs<0.0){ //right wet-dry
            solidWall=1;
            #pragma unroll
            for(k=0;k<3;k++){
                aux1 =  -(landl[k]*alfal[k] - betal[k]);
                dULrot[0] += aux1*eigel[k][0];
            }

        }else if(hL<TOL12 && hls<0.0){ //left wet-dry
            solidWall=1;
            #pragma unroll
            for(k=0;k<3;k++){
                aux1 =  -(landl[k]*alfal[k] - betal[k]);
                dURrot[0] += aux1*eigel[k][0];
            }

        }else{ //wet-wet
            #pragma unroll
            for(k=0;k<3;k++){
                aux1 =  -(landl[k]*alfal[k] - betal[k]);
                if(landl[k]<0.0){
                    #pragma unroll
                    for(j=0;j<3;j++){
                        dULrot[j] += aux1*eigel[k][j];
                    }
                }else{
                    #pragma unroll
                    for(j=0;j<3;j++){
                        dURrot[j] += aux1*eigel[k][j];
                    }
                }

                //enthropy fix contribution 
                aux1 =  -(landEl[k]*alfal[k]);
                if(landEl[k]<0.0){
                    #pragma unroll
                    for(j=0;j<3;j++){
                        dULrot[j] += aux1*eigel[k][j];
                    }
                }else{
                    #pragma unroll
                    for(j=0;j<3;j++){
                        dURrot[j] += aux1*eigel[k][j];
                    }
                }
            }
        }

        // Left cell X-Y contributions
        dUL[0] = dULrot[0];
        dUL[1] = dULrot[1]*nx - dULrot[2]*ny;
        dUL[2] = dULrot[1]*ny + dULrot[2]*nx;

        // Right cell Left cell X-Y contributions
        dUR[0] = dURrot[0];
        dUR[1] = dURrot[1]*nx - dURrot[2]*ny;
        dUR[2] = dURrot[1]*ny + dURrot[2]*nx;       


        //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc normal mass flux
        qnormalL = (huL*nx + hvL*ny) - dUL[0];

        // Neglect very small mass fluxes (wet-dry control)
        if(fabs(qnormalL)<TOL14){
            qnormalL=0.0;
            dUL[0] = (huL*nx + hvL*ny);
            dUR[0] = -(huR*nx + hvR*ny);
        }

        // Neglect very small momentum fluxes
        #pragma unroll
        for(k=1;k<3;k++){
            if(fabs(dUR[k])<TOL14){dUR[k]=0.0;}
            if(fabs(dUL[k])<TOL14){dUL[k]=0.0;}
        }

        //ccccccccccccccccccccccccccccccccccccccccccccccc Add contributions
        aux1 = arrays->length[idx]/arrays->area[id1];
        iw1=(arrays->idw1[idx]*arrays->ncells)+id1;
        arrays->dh[iw1]=dUL[0]*aux1;
        arrays->dhu[iw1]=dUL[1]*aux1;
        arrays->dhv[iw1]=dUL[2]*aux1;

        aux2 = arrays->length[idx]/arrays->area[id2];
        iw2=(arrays->idw2[idx]*arrays->ncells)+id2;
        arrays->dh[iw2]=dUR[0]*aux2;
        arrays->dhu[iw2]=dUR[1]*aux2;
        arrays->dhv[iw2]=dUR[2]*aux2;

        //} //end ((hL>TOL12)||(hR>TOL12)){ // Bucle de paredes mojadas


        //ccccccccccccccccccccccccccccccccccccccccccccccc update aux arrays
		arrays->qnormalL[idx] = qnormalL;
        localDt[idx]=dtl;

		arrays->solidWall[idx] = solidWall;
		arrays->solidWallByCell[iw1] = solidWall;
		arrays->solidWallByCell[iw2] = solidWall;


        //ccccccccccccccccccccccccccccccccccccccccccccccc update active cell list
 		if(wetWall==0){
            g_add_active_cells(arrays, id1);
            g_add_active_cells(arrays, id2);
		} 

	} // end iwall loop

}





////////////////////////////////////////////////////
//EXPORT_DLL void c_bound_calculus(t_mesh_cu *cmesh,t_parameters *sp);


////////////////////////////////////////////////////
__global__ void g_get_dtmin(t_arrays *arrays, double *localDt, int *idmin){
/*----------------------------*/	
	double minDt;

	minDt=localDt[(*idmin)-1];

    //dry domain
    //new rainDt implementation required
    if(minDt>1e5){
        minDt=1.0;
    }    
    minDt*=arrays->CFL;

    //storage time step
    arrays->dt=minDt;

    //printf("minDT %f idmin %d\n", minDt, (*idmin)-1);

}




////////////////////////////////////////////////////
__global__ void g_update_contributions(int nTasks, t_arrays *arrays){
/*----------------------------*/	
    int idx;
    int k;
    int NCwall,ncells;
	double total;

    NCwall=arrays->NCwall;
	ncells=arrays->ncells;

    int i = threadIdx.x+(blockIdx.x*blockDim.x);    
    if(i<nTasks){
        idx=arrays->actCells[i];

		total=0.0;
		total+=arrays->dh[idx];
		total+=arrays->dh[ncells+idx];
		total+=arrays->dh[2*ncells+idx];
		if(NCwall==4){
			total+=arrays->dh[3*ncells+idx];
		}
		if(fabs(total)<TOL14){
			total=0.0;
		}
		arrays->dh[idx]=total;

		total=0.0;
		total+=arrays->dhu[idx];
		total+=arrays->dhu[ncells+idx];
		total+=arrays->dhu[2*ncells+idx];
		if(NCwall==4){
			total+=arrays->dhu[3*ncells+idx];
		}
		if(fabs(total)<TOL14){
			total=0.0;
		}
		arrays->dhu[idx]=total;

		total=0.0;
		total+=arrays->dhv[idx];
		total+=arrays->dhv[ncells+idx];
		total+=arrays->dhv[2*ncells+idx];
		if(NCwall==4){
			total+=arrays->dhv[3*ncells+idx];
		}
		if(fabs(total)<TOL14){
			total=0.0;
		}
		arrays->dhv[idx]=total;

	}
}




////////////////////////////////////////////////////
__global__ void g_checkpos_h(int nTasks, t_arrays *arrays, int *check){
/*----------------------------*/
    int idx;
    int k;
	double updateh;

    int i = threadIdx.x+(blockIdx.x*blockDim.x);    
    if(i<nTasks){
        idx=arrays->actCells[i];

        updateh = arrays->h[idx] + arrays->dh[idx]*arrays->dt;

        if(updateh < (-TOL12)){   
            atomicAdd(check,1);
        }
    }

}


////////////////////////////////////////////////////
__global__ void g_reduce_dt(t_arrays *arrays){
/*----------------------------*/
    arrays->dt *= 0.8;
}




////////////////////////////////////////////////////
__global__ void g_set_new_dt(t_arrays *arrays){
/*----------------------------*/	
    double t0,tn;
	double minDt;
    
    t0=arrays->ti;
    tn=arrays->t;
	minDt=arrays->dt;
    
    arrays->dumpComponent=0;
    arrays->dumpState=0;

    //components
    if((tn+minDt) >= (t0 + arrays->indexDump*arrays->dtDump)){
        minDt = (t0 + arrays->indexDump*arrays->dtDump) - tn;
        arrays->dumpComponent=1;
    }

    //states
    //printf("indexOut %d\n",arrays->indexOut);
    if((tn+minDt) >= (t0 + arrays->indexOut*arrays->dtOut)){
        minDt = (t0 + arrays->indexOut*arrays->dtOut) - tn;
        arrays->dumpState=1;
        //printf("dumpState %d\n",arrays->dumpState);
    }    

    //last time
    if((tn+minDt) >= (arrays->tf)){
        minDt = (arrays->tf) - tn;
        arrays->dumpComponent=1;
        arrays->dumpState=1;
    }

    arrays->dt = minDt;
    arrays->t += minDt;

}



////////////////////////////////////////////////////
__global__ void g_update_cells(int nTasks, t_arrays *arrays){
/*----------------------------*/
	int idx;
    double dt;

    dt=arrays->dt;

    int i = threadIdx.x+(blockIdx.x*blockDim.x);    
    if(i<nTasks){
        idx=arrays->actCells[i];

		arrays->h[idx] += arrays->dh[idx]*dt;

        if(arrays->h[idx]>TOL12){ //wet cells
            arrays->sqrh[idx]=sqrt(arrays->h[idx]);

			if(arrays->h[idx] >= arrays->minh){
				arrays->hu[idx] += arrays->dhu[idx]*dt;
				arrays->hv[idx] += arrays->dhv[idx]*dt;

                if(fabs(arrays->hu[idx])<TOL14) arrays->hu[idx]=0.0;
                if(fabs(arrays->hv[idx])<TOL14) arrays->hv[idx]=0.0;

                arrays->u[idx] = arrays->hu[idx]/arrays->h[idx];
				arrays->v[idx] = arrays->hv[idx]/arrays->h[idx];
                arrays->modulou[idx] = sqrt(arrays->u[idx]*arrays->u[idx] + arrays->v[idx]*arrays->v[idx]);
			}else{
				arrays->hu[idx] = 0.0;
				arrays->hv[idx] = 0.0;
                arrays->u[idx] = 0.0;
				arrays->v[idx] = 0.0;
                arrays->modulou[idx] = 0.0;
			}


            //add walls to active-wall-array
            g_add_active_walls(arrays, idx);  


        }else{ //dry cell
            arrays->h[idx] = 0.0;
            arrays->sqrh[idx] = 0.0;
            arrays->hu[idx] = 0.0;
            arrays->hv[idx] = 0.0;
            arrays->u[idx] = 0.0;
            arrays->v[idx] = 0.0;
            arrays->modulou[idx] = 0.0;            
        }

    }

}


////////////////////////////////////////////////////
__global__ void g_check_wetdry(int nTasks, t_arrays *arrays){
/*----------------------------*/
	int idx;
	int id1, id2;
	int idw1,idw2;
	double h1,h2,z1,z2,hu,hv,un;
	double hu1,hu2,hv1,hv2;
	double nx,ny;

	int ncells;
	ncells=arrays->ncells;

    int i = threadIdx.x+(blockIdx.x*blockDim.x);    
    if(i<nTasks){
        idx=arrays->actWalls[i];

		id1=arrays->idx1[idx];
		id2=arrays->idx2[idx];
		idw1=arrays->idw1[idx];
		idw2=arrays->idw2[idx];

		h1=arrays->h[id1];
		h2=arrays->h[id2];

		z1=arrays->z[id1];
		z2=arrays->z[id2];

        //solidWalls from wall_calculus 
		if(arrays->solidWall[idx]==1){
			arrays->solidWallByCell[idw1*ncells+id1]=1;
			arrays->solidWallByCell[idw2*ncells+id2]=1;
		}

        //new solidWalls after update_cells
		if( (h1+z1<z2 && h2<TOL12 && h1>arrays->minh) || 
            (h2+z2<z1 && h1<TOL12 && h2>arrays->minh) ){
            arrays->solidWallByCell[idw1*ncells+id1]=1;
            arrays->solidWallByCell[idw2*ncells+id2]=1;
		}
	}

}


////////////////////////////////////////////////////
__global__ void g_update_wetdry_cells(int nTasks, t_arrays *arrays){
/*----------------------------*/
	int idx;
	int k, count, solidWall;
	double hun;
	double nx,ny,nx1,ny1;
    double daux;

	int ncells;
	ncells=arrays->ncells;

    int i = threadIdx.x+(blockIdx.x*blockDim.x);    
    if(i<nTasks){
        idx=arrays->actCells[i];

        if(arrays->h[idx] >= arrays->minh){

      		count=0;
      		for(k=0;k<arrays->NCwall;k++){ 

                solidWall=0;
      			if( (arrays->solidWallByCell[k*ncells+idx]==1) || //computational solidWalls
                    (arrays->neighCell[k*ncells+idx]==-1) ){ //closed boundary walls 
                    
                    solidWall=1;
      				count++;
      			}

      			// Modify velocity
                if(solidWall==1){

                    if(count==1){ //solidWall 1
                        nx=arrays->normalXbyCell[k*ncells+idx];
                        ny=arrays->normalYbyCell[k*ncells+idx];

                        hun=arrays->hu[idx]*nx + arrays->hv[idx]*ny;
                        arrays->hu[idx]-=hun*nx;
                        arrays->hv[idx]-=hun*ny;

                        //storage first solidWall normal
                        nx1=nx;
                        ny1=ny;

                    }else if(count==2){ //solidWall 2
                        if(arrays->NCwall==3){ //triangular cells
                            arrays->hu[idx] = 0.0;
                            arrays->hv[idx] = 0.0;
                        }else{ //square cells
                            nx=arrays->normalXbyCell[k*ncells+idx];
                            ny=arrays->normalYbyCell[k*ncells+idx];
                            if(fabs(nx1*nx+ny1*ny) < TOL9){ //n1*n2=0 orthogonal walls
                                arrays->hu[idx] = 0.0;
                                arrays->hv[idx] = 0.0;                           
                            }
                        }

                    }if(count>=3){
                        arrays->hu[idx] = 0.0;
                        arrays->hv[idx] = 0.0;                  
                    
                    }                      
                }

      		}

            if(fabs(arrays->hu[idx])<TOL14) arrays->hu[idx]=0.0;
            if(fabs(arrays->hv[idx])<TOL14) arrays->hv[idx]=0.0;

            arrays->u[idx] = arrays->hu[idx]/arrays->h[idx];
            arrays->v[idx] = arrays->hv[idx]/arrays->h[idx];
            arrays->modulou[idx] = sqrt(arrays->u[idx]*arrays->u[idx] + arrays->v[idx]*arrays->v[idx]);

        } //end if(arrays->h[idx] >= arrays->minh){
	}

}


////////////////////////////////////////////////////
__global__ void g_compute_mass_error(t_arrays *arrays){
/*----------------------------*/

	double massError;

    if(arrays->massOld>TOL3){
        massError=fabs(arrays->massNew-arrays->massOld)/(arrays->massOld)*100.0;
        if(massError<TOL14) massError=TOL14;
    }else{ 
        massError = TOL14;
    }

	arrays->massError=massError;
}


////////////////////////////////////////////////////
__global__ void g_reconstruct_active_elements(int nTasks, t_arrays *arrays){
/*----------------------------*/
	int idxw;
	int id1, id2;
	double h1,h2;

	int pos;

    int i = threadIdx.x+(blockIdx.x*blockDim.x);    
    if(i<nTasks){
        idxw=arrays->widx[i];

		id1=arrays->idx1[i];
		id2=arrays->idx2[i];
        
        h1=arrays->h[id1];
		h2=arrays->h[id2];

        if((h1 > TOL12) || (h2 > TOL12)){
            arrays->activeW[i]=1;
            pos = atomicAdd(&(arrays->nActWalls), 1); //obtenemos el índice para el nuevo elemento
            arrays->actWalls[pos] = idxw; //copiamos el nuevo elemento al array

            if(h1 > TOL12){
                g_add_active_cells(arrays, id1);
            }
            if(h2 > TOL12){
                g_add_active_cells(arrays, id2);
            }
        }

        arrays->solidWall[i]=0;

        arrays->qnormalL[i]=0.0;
        arrays->localDt[i]=1e6; 

    }

}




////////////////////////////////////////////////////
__device__ void g_add_active_cells(t_arrays *arrays, int id){
/*----------------------------*/
    int access, pos;

    access = atomicAdd(&(arrays->activeC[id]), 1); //dos cucores no acceden a la vez
    if(access==0){ //primer cucore que accede    

        pos = atomicAdd(&(arrays->nActCells), 1); //obtenemos el índice para el nuevo elemento
        arrays->actCells[pos] = id; //copiamos el nuevo elemento al array
       
        //printf("Cell %d added to active cell position %d\n",id,pos);         
    }
}

////////////////////////////////////////////////////
__device__ void g_add_active_walls(t_arrays *arrays, int id){
/*----------------------------*/
    int k,idxw;
    int access, pos;

    for(k=0;k<arrays->NCwall;k++){
        if(arrays->neighCell[k*arrays->ncells+id] >=0 ){ //internal wall
            idxw=arrays->neighWall[k*arrays->ncells+id];

            access = atomicAdd(&(arrays->activeW[idxw]), 1); //dos cucores no acceden a la vez
            if(access==0){ //primer cucore que accede

                pos = atomicAdd(&(arrays->nActWalls), 1); //obtenemos el índice para el nuevo elemento
                arrays->actWalls[pos] = idxw; //copiamos el nuevo elemento al array

                //printf("Cell %d :: Wall %d added to active wall position %d\n",id,idxw,pos);
            }
        }
    }

}



