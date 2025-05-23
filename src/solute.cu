#include "solute.cuh"

#if SET_SOLUTE
////////////////////////////////////////////////////
__global__ void g_initialize_solute_delta(int nTasks, t_arrays *arrays){
/*----------------------------*/
    int i = threadIdx.x+(blockIdx.x*blockDim.x);    
    if(i<nTasks){
        arrays->dhphi[i]=0.0;
        //arrays->Bwall[i]=0.0;
    }
}

////////////////////////////////////////////////////
__global__ void g_wall_solute_calculus(int nTasks, t_arrays *arrays, double *localDt){
/*----------------------------*/

    int idx;
    int ncells = arrays->ncells;
    int NCwall = arrays->NCwall;

    int id1,id2;
    int idw1,idw2; 

    //nSolutes*ncell indices
    int sid1, sid2;
    int siw1, siw2;
    
    //hydrodynamic parameters
    double qnormalL;
    double length;
    double areaL, areaR;

    // Solute variable
    double phiL, phiR;
    double dhphi;
    double dphi;

    double aux1,aux2;

    int nActWalls = arrays->nActWalls;
    int jphi;
    int iactWall;

    int i = threadIdx.x+(blockIdx.x*blockDim.x);    
    if(i<nTasks){
        
        //idx=arrays->actWalls[i]; //compact 

        //solute index
        jphi=(int)(i/nActWalls);
        
        //wall index
        iactWall=i-jphi*nActWalls;
        idx=arrays->actWalls[iactWall];
        
        //cells index
        id1=arrays->idx1[idx];
        id2=arrays->idx2[idx];

        //ccccccccccccccccccccccccccccccccccccccccccccccccc Wall-averaged values		 
        qnormalL = arrays->qnormalL[idx];

        length = arrays->length[idx];
        areaL = arrays->area[id1];
        areaR = arrays->area[id2];

        //cccccccccccccccccccccccccccccccccccccccccccccccccc Convective Solute Transport
        //for(jphi=0;jphi<arrays->nSolutes;jphi++){ //compact      

            //cccccccccccccccccccccccccccccccccccccccccccccc solute wall flux
            sid1 = jphi*ncells+id1;
            sid2 = jphi*ncells+id2;
            
            phiL = arrays->phi[sid1];
            phiR = arrays->phi[sid2];

            dphi=0.5*(phiL+phiR)-SIGN(qnormalL)*0.5*(phiR-phiL);
            dhphi=qnormalL*dphi;

            //ccccccccccccccccccccccccccccccccccccccccccccccc Add solute contributions
            aux1 = length/areaL;
            aux2 = length/areaR;

            idw1 = arrays->idw1[idx];
            idw2 = arrays->idw2[idx]; 
 
            siw1=jphi*(ncells*NCwall)+(idw1*ncells)+id1;
            siw2=jphi*(ncells*NCwall)+(idw2*ncells)+id2;
            
            arrays->dhphi[siw1] = - dhphi*aux1;
            arrays->dhphi[siw2] = dhphi*aux2;
            
        //} //compact 


	} // end iwall loop

}


////////////////////////////////////////////////////
__global__ void g_bound_solute_calculus(int nTasks, t_arrays *arrays){		
/*----------------------------*/

	int cidx;
    int ncells = arrays->ncells;
    int NCwall = arrays->NCwall;
	double nSolutes = arrays->nSolutes;

    int sid, siw0;

    double h, hu, hv;
    double hun;
    double length, area;

    int jphi;
    int iBoundCell;

    int nBoundCells=arrays->nTotalBoundCells;    

	int i = threadIdx.x+(blockIdx.x*blockDim.x);  
    if(i<nTasks){

        //cidx=arrays->cidxBound[i]; //compact
        //iBoundCell=i;

        //solute index
        jphi=(int)(i/nBoundCells);

        //cell index
        iBoundCell=i-jphi*nBoundCells;
        cidx=arrays->cidxBound[iBoundCell];

		hu=arrays->hu[cidx];
		hv=arrays->hv[cidx];
        area=arrays->area[cidx];

        hun = hu*arrays->nxWallBound[iBoundCell] + hv*arrays->nyWallBound[iBoundCell];
        length = arrays->lWallBound[iBoundCell];

        //for(jphi=0;jphi<nSolutes;jphi++){ //compact

            sid = jphi*ncells+cidx;
            siw0 = jphi*(ncells*NCwall)+cidx;

            arrays->dhphi[siw0] -= hun*arrays->phi[sid]*length/area;

        //} //compact

	}
	//__syncthreads(); // Sincronizar todos los hilos del bloque

}


////////////////////////////////////////////////////
__global__ void g_update_solute_contributions(int nTasks, t_arrays *arrays){
/*----------------------------*/	
    int idx;
	int ncells=arrays->ncells;
    int NCwall=arrays->NCwall;
	    
    int jphi;
    int iactCell;

    int nActCells=arrays->nActCells;
    
    double total;
    int siw0;

    int i = threadIdx.x+(blockIdx.x*blockDim.x);    
    if(i<nTasks){

        //idx=arrays->actCells[i]; //compact
        
        //solute index
        jphi=(int)(i/nActCells);

        //cell index
        iactCell=i-jphi*nActCells;
        idx=arrays->actCells[iactCell];

        //for(jphi=0;jphi<arrays->nSolutes;jphi++){ //compact

            siw0=jphi*(ncells*NCwall)+idx;

            total=0.0;
            total+=arrays->dhphi[siw0];
            total+=arrays->dhphi[siw0+ncells];
            total+=arrays->dhphi[siw0+2*ncells];
            if(NCwall==4){
                total+=arrays->dhphi[siw0+3*ncells];     
            }
            if(fabs(total)<TOL14){
                total=0.0;
            }
            arrays->dhphi[siw0]=total;                
                

        //} //compact

        
    }
}


////////////////////////////////////////////////////
__global__ void g_update_solute_cells(int nTasks, t_arrays *arrays){
/*----------------------------*/
	int idx;
    double dt;
    int siw0,sid;

    int NCwall=arrays->NCwall;
    int ncells=arrays->ncells;

    int jphi;
    int iactCell;

    int nActWalls=arrays->nActWalls;
    int nActCells=arrays->nActCells;    

    dt=arrays->dt;

    int i = threadIdx.x+(blockIdx.x*blockDim.x);    
    if(i<nTasks){

        //idx=arrays->actCells[i]; //compact

        //solute index
        jphi=(int)(i/nActCells);

        //cell index
        iactCell=i-jphi*nActCells;
        idx=arrays->actCells[iactCell];

        if(arrays->h[idx]>TOL12){ //wet cells

            //for(jphi=0;jphi<arrays->nSolutes;jphi++){ //compact

                siw0 = jphi*(ncells*NCwall)+idx;
                sid = jphi*ncells+idx;

                arrays->hphi[sid] += arrays->dhphi[siw0]*dt;
                if(fabs(arrays->hphi[sid])<TOL14){
                    arrays->hphi[sid]=0.0;
                }
                
                arrays->phi[sid] = arrays->hphi[sid]/arrays->h[idx];

            //} //compact

        }
    }

}


////////////////////////////////////////////////////
__global__ void g_initialize_solute_diffusion_delta(int nTasks, t_arrays *arrays){
/*----------------------------*/
    int k, jphi, idx, siw;
    int ncells = arrays->ncells;
    int NCwall = arrays->NCwall;

    int i = threadIdx.x+(blockIdx.x*blockDim.x);    
    if(i<nTasks){
        arrays->BTcell[i] = 0.0;
        arrays->localDtd[i] = 1e5;

        //solute index
        jphi=(int)(i/ncells);

        //cell index
        idx=i-jphi*ncells; 

        for(k=k;k<NCwall;k++){
            siw=jphi*(ncells*NCwall)+(k*ncells)+idx;
            arrays->dhphi[siw]=0.0;
            arrays->Bwall[siw]=0.0;
        }
    }
}


////////////////////////////////////////////////////
__global__ void g_wall_solute_diffusion_calculus(int nTasks, t_arrays *arrays){
/*----------------------------*/

    int idx;
	int k;
    int id1,id2;
    int siw1, siw2;
    int idw1, idw2;

    int ncells = arrays->ncells;
    int NCwall = arrays->NCwall;
    int typeDiff;

    int nActWalls = arrays->nActWalls;
    int jphi;
    int iactWall;

	double hL,hR;
	double sqrhL,sqrhR; 
	double uL, uR, vL, vR;
    double areaL;
    double areaR;
    double minh = arrays->minh;   

	double hbar,ubar,vbar,modU2; 
	double nx,ny;
	double gp; 
    double length;
    double nman2wall;
    double distNormal;  
 
    double aux1,aux2,aux3,aux4;

    // Solute variable
    double kL, kT;
    double ustar;
    double Cxx,Cxy,Cyy;
    double contrib=0.0;
    double Cnn;
    double unbar, Ca;

    double dt = arrays->dt;
	
    int i = threadIdx.x+(blockIdx.x*blockDim.x);    
    if(i<nTasks){

        //solute index
        jphi=(int)(i/nActWalls);
        
        //wall index
        iactWall=i-jphi*nActWalls;
        idx=arrays->actWalls[iactWall];        

        //cells index
        id1=arrays->idx1[idx];
        id2=arrays->idx2[idx];
    
        hL = arrays->h[id1];
        hR = arrays->h[id2];

        if(hL>minh && hR>minh){ //wet-wet walls
		
            sqrhL=arrays->sqrh[id1];
            sqrhR=arrays->sqrh[id2];

            uL  = arrays->u[id1];
            uR  = arrays->u[id2];
                    
            vL  = arrays->v[id1];
            vR  = arrays->v[id2];

            areaL = arrays->area[id1];
            areaR = arrays->area[id2];


            //ccccccccccccccccccccccccccccccccccccccccccccccccc Wall-averaged values		 
            hbar = 0.5*(hL+hR);

            aux1 = sqrhL + sqrhR;
            ubar = (uL*sqrhL + uR*sqrhR)/aux1;
            vbar = (vL*sqrhL + vR*sqrhR)/aux1;

            if(fabs(ubar) < TOL12) ubar = 0.0;
            if(fabs(vbar) < TOL12) vbar = 0.0;            

            modU2 = ubar*ubar+vbar*vbar;

            //ccccccccccccccccccccccccccccccccccccccccccccccccc Edge values
            nx = arrays->normalX[idx];
            ny = arrays->normalY[idx];		
            gp = arrays->gp[idx];
                    
            nman2wall = arrays->nman2wall[idx];
            distNormal = arrays->distNormal[idx];
            length = arrays->length[idx];


            //for(j=0;j<arrays->nSolutes;j++){

                //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc diffusion contribution calculation
                if(arrays->typeDiff[jphi]==NON_DIFF){
                    kL = 0.0;
                    kT = 0.0;

                }else if(arrays->typeDiff[jphi]==CONSTANT_DIFF){
                    kL = arrays->k_xx[jphi];
                    kT = arrays->k_yy[jphi];

                }else if(arrays->typeDiff[jphi]==ANISOTROPIC_DIFF){
                    ustar = sqrt(gp*nman2wall*modU2/cbrt(hbar));
                    kL = arrays->k_xx[jphi]*hbar*ustar;
                    kT = arrays->k_yy[jphi]*hbar*ustar;
                }
                //printf("kL %lf kT %lf\n", kL, kT);
    
                if(modU2>TOL9){
                    aux1 = kL*ubar*ubar/modU2 + kT*vbar*vbar/modU2;
                    aux2 = (kL-kT)*ubar*vbar/modU2;
                    aux3 = kT*ubar*ubar/modU2 + kL*vbar*vbar/modU2; 
                }else{
                    aux1 = 0.0;
                    aux2 = 0.0;
                    aux3 = 0.0;
                }
                
                Cxx = MAX(aux1,_Dm_);
                Cyy = MAX(aux3,_Dm_);
                Cxy = MAX(aux2,0.0);  
                
                Cnn = nx*(Cxx*nx + Cxy*ny) + ny*(Cxy*nx + Cyy*ny);
                Cnn = MAX(Cnn,0.0);

                // Numerical diff correction term
                unbar = fabs(ubar*nx + vbar*ny);
                Ca = 0.5*unbar*(distNormal - unbar*dt);
                Cnn = MAX(Cnn-Ca,0.0);

                //Wall contribution
                contrib = Cnn*hbar/distNormal;
                //printf("nx %lf ny %lf contrib %lf\n",nx, ny, contrib);
                

                //ccccccccccccccccccccccccccccccccccccccccccccccc Add solute contributions
                aux1 = length/areaL;
                aux2 = length/areaR;

                idw1 = arrays->idw1[idx];
                idw2 = arrays->idw2[idx];

                siw1=jphi*(ncells*NCwall)+(idw1*ncells)+id1;
                siw2=jphi*(ncells*NCwall)+(idw2*ncells)+id2;

                arrays->Bwall[siw1] = (contrib/hL)*aux1;
                arrays->Bwall[siw2] = (contrib/hR)*aux2;

            //}

        } //end if(hL>minh && hR>minh){ //wet-wet walls

	} // end iwall loop

}


////////////////////////////////////////////////////
__global__ void g_update_solute_diffusion_contributions(int nTasks, t_arrays *arrays, double *localDtd){
/*----------------------------*/	
    int idx;
    int siw0,sid;

    int ncells=arrays->ncells;
    int NCwall=arrays->NCwall;

    int nActCells=arrays->nActCells;
    int jphi;
    int iactCell;

    double Dtd;
	double total;

    Dtd = 10.*arrays->dt;
    
    int i = threadIdx.x+(blockIdx.x*blockDim.x);
    if(i<nTasks){
        
        //idx=arrays->actCells[i]; //compact

        //solute index
        jphi=(int)(i/nActCells);

        //cell index
        iactCell=i-jphi*nActCells;
        idx=arrays->actCells[iactCell]; 

        //for(jphi=0;jphi<arrays->nSolutes;jphi++){  //compact
            
            siw0=jphi*(ncells*NCwall)+idx;
            total = arrays->Bwall[siw0];
            total += arrays->Bwall[siw0+ncells];
            total += arrays->Bwall[siw0+2*ncells];
            if(NCwall==4){
                total += arrays->Bwall[siw0+3*ncells];
            }

            sid = jphi*ncells+idx;
            arrays->BTcell[sid] = total;

            //Diffusion time step
            if(fabs(total)>TOL12){
                arrays->localDtd[sid] = 1./fabs(total);
            }
            //printf("cell %d BTcell %lf Dtd %lf\n",idx, arrays->BTcell[sid], arrays->localDtd[sid]);

            //store solute concentration
            arrays->dhphi[sid] = arrays->phi[sid]; //phi is stored in [sid]-position of dhphi           

        //} //compact

    }

}


//////////////////////////////////////////////////////
__global__ void g_get_solute_diffusion_dtmin(t_arrays *arrays, double *localDtd, int *idmin){
/*----------------------------*/	
	double minDtd;

	minDtd=localDtd[(*idmin)-1];

    //storage time step
    arrays->Dtd=minDtd;
    

}


//////////////////////////////////////////////////////
__global__ void g_update_solute_diffusion_cells(int nTasks, t_arrays *arrays, double *Dtd){
/*----------------------------*/	

    int idx,neighid;
    int siw,sid;
    int k;
    int ncells=arrays->ncells;
    int NCwall=arrays->NCwall;
    int nw_calc=arrays->nw_calc;

    double dtd;
	double total;
    double contrib1, contrib2;

    int nActCells=arrays->nActCells;
    int jphi;
    int iactCell;
    
    int i = threadIdx.x+(blockIdx.x*blockDim.x);
    if(i<nTasks){

        //idx=arrays->actCells[i]; //compact

        //solute index
        jphi=(int)(i/nActCells);

        //cell index
        iactCell=i-jphi*nActCells;
        idx=arrays->actCells[iactCell]; 

        //for(jphi=0;jphi<arrays->nSolutes;jphi++){  //compact

            sid = jphi*ncells+idx;
            contrib1 = arrays->dhphi[sid]*(1.0-arrays->BTcell[sid]*(*Dtd)); //phi is stored in dhphi[sid]

            contrib2 = 0.0;
            for(k=0;k<NCwall;k++){
                siw=jphi*(ncells*NCwall)+(k*ncells)+idx;
                neighid=arrays->neighCell[k*ncells+idx];
                if(neighid>=0){
                    contrib2 += arrays->dhphi[jphi*ncells+neighid]*arrays->Bwall[siw]*(*Dtd); //phi is stored in dhphi[sid] 
                }
            }
            arrays->phi[sid]=contrib1+contrib2;

            //update conservative variable
            arrays->hphi[sid]=arrays->phi[sid]*arrays->h[idx]; 

        //}  //compact
    }

}



#endif 