#include "particle.cuh"

#if SET_SOLUTE || SET_SED

////////////////////////////////////////////////////
__global__ void g_initialize_particle_delta(int nTasks, t_arrays *arrays){
/*----------------------------*/
    int i = threadIdx.x+(blockIdx.x*blockDim.x);    
    if(i<nTasks){
        arrays->dhphi[i]=0.0;
        //arrays->Bwall[i]=0.0;
    }
    //printf("hphi%lf\n", arrays->dhphi[i]);
}

////////////////////////////////////////////////////
__global__ void g_wall_particle_calculus(int nTasks, t_arrays *arrays, double *localDt){
/*----------------------------*/

    int idx;
    int ncells = arrays->ncells;
    int NCwall = arrays->NCwall;

    int nSolutes=arrays->nSolutes;
    int nSediments=arrays->nSediments;
    int nParticles=nSolutes + nSediments;

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

    double u,v, moduloU;
    double hlayer;
    double nman2wall;
    double ustar;

    double aux1,aux2;

    int nActWalls = arrays->nActWalls;
    int jphi;
    int iactWall;

    int i = threadIdx.x+(blockIdx.x*blockDim.x);    

    
    //printf("nSolutes %d nSediments %d nParticles %d\n", nSolutes, arrays->nSediments, nParticles);
     

    if(i<nTasks){
        
        #if SET_SOLUTE_UNROLL==0 || SET_SED_UNROLL==0 //compact
        //wall index
        idx=arrays->actWalls[i]; 

        #elif SET_SOLUTE_UNROLL==1 || SET_SED_UNROLL==1 //unroll 
        //solute index
        jphi=(int)(i/nActWalls);
        
        //wall index
        iactWall=i-jphi*nActWalls;
        idx=arrays->actWalls[iactWall];
        #endif
        
        //cells index
        id1=arrays->idx1[idx];
        id2=arrays->idx2[idx];

        //ccccccccccccccccccccccccccccccccccccccccccccccccc Wall-averaged values		 
        qnormalL = arrays->qnormalL[idx];

        length = arrays->length[idx];
        areaL = arrays->area[id1];
        areaR = arrays->area[id2];

        //cccccccccccccccccccccccccccccccccccccccccccccccccc Convective Solute Transport
        #if SET_SOLUTE_UNROLL==0 || SET_SED_UNROLL==0 //compact
        for(jphi=0;jphi<(nParticles);jphi++){
        #endif    

            //cccccccccccccccccccccccccccccccccccccccccccccc solute wall flux
            sid1 = jphi*ncells+id1;
            sid2 = jphi*ncells+id2;
            
            phiL = arrays->phi[sid1];
            phiR = arrays->phi[sid2];

            // if(idx == 358){
            //     printf("phiL %lf phiR %lf\n", phiL, phiR);
            // }

            dphi=0.5*(phiL+phiR)-SIGN(qnormalL)*0.5*(phiR-phiL);


            #if SET_MULTILAYER
            #if SET_MULTILAYER_VELOCITY

            double k = 0.41;
            double z_normalized;
            double ks;
            double Re = 2800;
            double U_fonc, u_prom;
            double A[12];

            // MODEL LINEAR
            // A[0] = 0.0833333;
            // A[1] = 0.25;
            // A[2] = 0.416667;
            // A[3] = 0.583333;
            // A[4] = 0.75;
            // A[5] = 0.916667;
            // A[6] = 1.08333;
            // A[7] = 1.25;
            // A[8] = 1.41667;
            // A[9] = 1.58333;
            // A[10] = 1.75;
            // A[11] = 1.91667;

            // MODEL LOG
            // EXP BORJA LATORRE
            // A[0]=	0.6325546093710857;
            // A[1]=	0.8178943911224351;
            // A[2]=	0.9040724819724524;
            // A[3]=	0.9608365392691058;
            // A[4]=	1.0032341728737848;
            // A[5]=	1.0370880294931535;
            // A[6]=	1.0652706448812252;
            // A[7]=	1.089412263723802;
            // A[8]=	1.1105277289681712;
            // A[9]=	1.1292918873060287;
            // A[10]=	1.1461763210204554;
            // A[11]=	1.1615235820078802;

            //Guadalquivir
            A[0]=	0.8223322102639082;
            A[1]=	0.9119480013820693;
            A[2]=	0.953616971290203;
            A[3]=	0.9810636201081068;
            A[4]=	1.0015637925002303;
            A[5]=	1.0179328640221537;
            A[6]=	1.031559767808893;
            A[7]=	1.043232762408364;
            A[8]=	1.0534425463242685;
            A[9]=	1.0625154225207534;
            A[10]=	1.0706794112262679;
            A[11]=	1.0781001436879587;


            dhphi=(qnormalL*A[jphi]/nSolutes)*dphi;

            #else
            
            dhphi=(qnormalL/nSolutes)*dphi;

            #endif
            #else
            
            dhphi=qnormalL*dphi;

            #endif
            
            //dhphi=qnormalL*dphi;

            //ccccccccccccccccccccccccccccccccccccccccccccccc Add solute contributions
            aux1 = length/areaL;
            aux2 = length/areaR;

            idw1 = arrays->idw1[idx];
            idw2 = arrays->idw2[idx]; 
 
            siw1=jphi*(ncells*NCwall)+(idw1*ncells)+id1;
            siw2=jphi*(ncells*NCwall)+(idw2*ncells)+id2;
            
            arrays->dhphi[siw1] = - dhphi*aux1;
            arrays->dhphi[siw2] = dhphi*aux2;

        // if(idx == 80500){
        //     printf("dphi %lf, dhphi %lf\n", dphi, arrays->dhphi[siw2]);
        // }
            // if(idx == 358){
            //     printf("dhphi %lf jphi %d\n", arrays->dhphi[siw1], jphi);
            // }


        #if SET_SOLUTE_UNROLL==0 || SET_SED_UNROLL==0 //compact
        }
        #endif

	} // end iwall loop

}


////////////////////////////////////////////////////
__global__ void g_bound_particle_calculus(int nTasks, t_arrays *arrays){		
/*----------------------------*/

	int cidx;
    int ncells = arrays->ncells;
    int NCwall = arrays->NCwall;
	int nSolutes = arrays->nSolutes;
    int nSediments = arrays->nSediments;
    int nParticles = nSolutes + nSediments;

    int sid, siw0;

    double h, hu, hv;
    double hun;
    double length, area;

    int jphi;
    int iBoundCell;

    int nBoundCells=arrays->nTotalBoundCells;    

	int i = threadIdx.x+(blockIdx.x*blockDim.x);  
    if(i<nTasks){

        #if SET_SOLUTE_UNROLL==0 || SET_SED_UNROLL==0 //compact
        //cell index
        iBoundCell=i;
        cidx=arrays->cidxBound[i];

        #elif SET_SOLUTE_UNROLL==1  || SET_SED_UNROLL==1 //unroll 
        //solute index
        jphi=(int)(i/nBoundCells);

        //cell index
        iBoundCell=i-jphi*nBoundCells;
        cidx=arrays->cidxBound[iBoundCell];
        #endif

		hu=arrays->hu[cidx];
		hv=arrays->hv[cidx];
        area=arrays->area[cidx];

        #if SET_MULTILAYER
        hun = hu*arrays->nxWallBound[iBoundCell] + hv*arrays->nyWallBound[iBoundCell];  
        hun /=nSOlutes;
        #else
        hun = hu*arrays->nxWallBound[iBoundCell] + hv*arrays->nyWallBound[iBoundCell];  
        #endif

        length = arrays->lWallBound[iBoundCell];

        #if SET_SOLUTE_UNROLL==0 || SET_SED_UNROLL==0 //compact
        for(jphi=0;jphi<nParticles;jphi++){
        #endif

            sid = jphi*ncells+cidx;
            siw0 = jphi*(ncells*NCwall)+cidx;

            arrays->dhphi[siw0] -= hun*arrays->phi[sid]*length/area;

        #if SET_SOLUTE_UNROLL==0 || SET_SED_UNROLL==0 //compact
        }
        #endif

        
        
	}
	//__syncthreads(); // Sincronizar todos los hilos del bloque

}


////////////////////////////////////////////////////
__global__ void g_update_particle_contributions(int nTasks, t_arrays *arrays){
/*----------------------------*/	
    int idx;
	int ncells=arrays->ncells;
    int NCwall=arrays->NCwall;
	    
    int jphi;
    int iactCell;

    int nActCells=arrays->nActCells;

    int nSolutes = arrays->nSolutes;
    int nSediments = arrays->nSediments;
    int nParticles = nSolutes + nSediments;

    double total;
    int siw0;

    int i = threadIdx.x+(blockIdx.x*blockDim.x);    
    if(i<nTasks){

        #if SET_SOLUTE_UNROLL==0 || SET_SED_UNROLL==0  //compact
        //cell index
        idx=arrays->actCells[i];

        #elif SET_SOLUTE_UNROLL==1 || SET_SED_UNROLL==1  //unroll 
        //solute index
        jphi=(int)(i/nActCells);

        //cell index
        iactCell=i-jphi*nActCells;
        idx=arrays->actCells[iactCell];
        #endif        

        #if SET_SOLUTE_UNROLL==0 || SET_SED_UNROLL==0 //compact
        for(jphi=0;jphi<arrays->nParticles;jphi++){
        #endif

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
                
        #if SET_SOLUTE_UNROLL==0 || SET_SED_UNROLL==0 //compact
        }
        #endif

        
    }
}


////////////////////////////////////////////////////
__global__ void g_update_particle_cells(int nTasks, t_arrays *arrays){
/*----------------------------*/
	int idx;
    double dt;
    int siw0,sid;

    int NCwall=arrays->NCwall;
    int ncells=arrays->ncells;
    int nSolutes=arrays->nSolutes;
    int nSediments=arrays->nSediments;
    int nParticles=nSolutes+nSediments;

    int jphi;
    int iactCell;

    int nActWalls=arrays->nActWalls;
    int nActCells=arrays->nActCells;  

    double hlayer;

    dt=arrays->dt;

    int i = threadIdx.x+(blockIdx.x*blockDim.x);    
    if(i<nTasks){

        #if SET_SOLUTE_UNROLL==0 || SET_SED_UNROLL==0 //compact 
        //cell index
        idx=arrays->actCells[i];

        #elif SET_SOLUTE_UNROLL==1 || SET_SED_UNROLL==1 //unroll  
        //solute index
        jphi=(int)(i/nActCells);

        //cell index
        iactCell=i-jphi*nActCells;
        idx=arrays->actCells[iactCell];
        #endif 


        if(arrays->h[idx]>TOL12){ //wet cells

            hlayer=arrays->h[idx]/nSolutes;

            #if SET_SOLUTE_UNROLL==0 || SET_SED_UNROLL==0 //compact 
            for(jphi=0;jphi<nParticles;jphi++){
            #endif

                siw0 = jphi*(ncells*NCwall)+idx;
                sid = jphi*ncells+idx;

                arrays->hphi[sid] += arrays->dhphi[siw0]*dt;
                if(fabs(arrays->hphi[sid])<TOL14){
                    arrays->hphi[sid]=0.0;
                }

                #if SET_MULTILAYER
                arrays->phi[sid] = arrays->hphi[sid]/hlayer;
                #else
                arrays->phi[sid] = arrays->hphi[sid]/arrays->h[idx];
                #endif
                

            #if SET_SOLUTE_UNROLL==0  || SET_SED_UNROLL==0//compact 
            }
            #endif

        }
    }

}

////////////////////////////////////////////////////
__global__ void g_multilayer_implicit_update_solute_cells(int nTasks, t_arrays *arrays){
/*----------------------------*/

int idx1,idx2,idx;
int sid1,sid2,sid3;
int jphi;

int ncells=arrays->ncells;
int nInterfaces= arrays->nSolutes;

double phij1, phij2;
double hphi1, hphi2;
double aux1,aux2;
double hlayer;
double sqrhL, sqrhR;

double dt;
double epsis1 = 0.005;

//double A,C,B;
double Bi,Ci;
double Af,Bf;
double Bbis [21];
double u,v,moduloU;
double uL, uR, vL, vR, hL, hR;
double modU2;
double ustar;
double nman;
double gp;
double hbar, ubar,vbar;


double sigmaD, sigmap, sigmah, sigmaMax;
double As, Ai, Asigma;
double Ri;
double ra = 0.99;
double n = 1;
double k = 0.41;
double z_depth;

#if EDDY_VISCOSITY || EDDY_VISCOSITY_LINEAR || EDDY_VISCOSITY_PARABOLIC
double A [20];
double B [20];
double C [20];
double epsia;
double epsib;
double epsis2[20];
#else
double A;
double B;
double C;
#endif


dt=arrays->dt;
    
    int i = threadIdx.x+(blockIdx.x*blockDim.x);    
    if(i<nTasks){
        idx=arrays->actCells[i]; //compact
        
        double r [21];
        double rbis [21];


        if(arrays->h[idx] >= arrays->minh){

            u = arrays->u[idx];
            v = arrays->v[idx];
                
            moduloU = arrays->modulou[idx];

            hlayer = arrays->h[idx]/(nInterfaces);
            nman = arrays->nman[idx];
            ustar = nman*moduloU*sqrt(_g_/cbrt(arrays->h[idx]));

            // if(idx == 1035457){
            //     printf("depth %.12lf ustar %.12lf\n", arrays->h[idx], ustar);
            //     printf("manning %.12lf moduloU %.12lf\n", nman, moduloU);
            // } 

            //cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc NON CONSTANT VERTICAL DIFFUSION 
            //ccccccccccccccccccccccccccc ESTUARY VERSION
            #if EDDY_VISCOSITY
                
                sigmap = -0.3;
                sigmah = -0.8;

                //calculation of the parameters needed
                sigmaMax = ((1+sigmap)*(1+sigmap))/4.;
                
                As = ra*sigmaMax;
                Ai = (sigmah + 1)*(sigmap-sigmah);

                //first values of eddy viscosity
            
                //ustar = nman2wall*sqrt(_g_*moduloU*moduloU/cbrt(arrays->h[idx]));

                for(jphi=0;jphi<nInterfaces;jphi++){
                    sigmaD = (jphi*hlayer + hlayer/2.)/(arrays->h[idx]) - 1;

                    if(sigmaD <= sigmah){
                        Asigma = (sigmaD + 1)*(sigmap - sigmaD);
                    }else if (sigmaD >sigmah && n == 1.){
                        Asigma = (Ai-As)*fabs(sigmaD/sigmah) + As;
                    }else if (sigmaD >sigmah && n == 2.){
                        Asigma = (Ai-As)*fabs(sigmaD/sigmah)*fabs(sigmaD/sigmah) + As;
                    }

                    epsis2[jphi] = k*ustar*arrays->h[idx]*Asigma;


                    if(idx == 1035457){
                        printf("epsis %.12lf sigmaD %.12lf\n", epsis2[jphi], sigmaD);
                        printf("depth %.12lf ustar %.12lf\n", arrays->h[idx], ustar);
                    } 

                }

                // sigmaD = (hlayer/2.)/(arrays->h[idx]) - 1.;

                // if(sigmaD <= sigmah){
                //     Asigma = (sigmaD + 1)*(sigmap - sigmaD);

                // }else if (sigmaD >sigmah && n == 1.){
                //     Asigma = (Ai-As)*fabs(sigmaD/sigmah) + As;

                // }else if (sigmaD >sigmah && n == 2.){
                //     Asigma = (Ai-As)*fabs(sigmaD/sigmah)*fabs(sigmaD/sigmah) + As;

                // }else if (sigmaD >sigmah && n == 0.){
                //     Asigma = Ai;
                // }

                //ustar = nman2wall*sqrt(_g_*modU2/cbrt(hbar));

                //epsis2 = k*ustar*arrays->h[idx]*Asigma;

                //epsis2 = epsia;

                // Bi = 1+dt/(hlayer*hlayer)*epsis2;
                // Ci = -(dt*epsis2)/(hlayer*hlayer);

                Bi = 1.+dt/(hlayer*hlayer)*((epsis2[1]+epsis2[0])/2.);
                Ci = -((epsis2[1]+epsis2[0])/2.)*dt/(hlayer*hlayer);

                B[0] = Bi;
                C[0] = Ci;

                //sigmaD = (hlayer+hlayer/2.)/(arrays->h[idx]) - 1;

                // if(sigmaD <= sigmah){
                //     Asigma = (sigmaD + 1)*(sigmap - sigmaD);

                // }else if (sigmaD >sigmah && n == 1.){
                //     Asigma = (Ai-As)*fabs(sigmaD/sigmah) + As;

                // }else if (sigmaD >sigmah && n == 2.){
                //     Asigma = (Ai-As)*fabs(sigmaD/sigmah)*fabs(sigmaD/sigmah) + As;
                // }

               
                // //epsis2 = k*ustar*arrays->h[idx]*Asigma;
                // epsis2 = epsia;


                A[1] = -((epsis2[1]+epsis2[0])/2.)*dt/(hlayer*hlayer);
                B[1] = 1.+((epsis2[2]+2.*epsis2[1]+epsis2[0])/2.)*(dt/(hlayer*hlayer));
                C[1] = -((epsis2[2]+epsis2[1])/2.)*dt/(hlayer*hlayer);

                Bbis[0] = Bi;
                r[0] = arrays->phi[idx];
                r[1] = arrays->phi[ncells+idx];
                rbis[0] = r[0];
                rbis[1] = arrays->phi[ncells+idx] - (A[1]/Bi)*r[0];
                Bbis[1] = B[1]-(A[1]*Ci)/Bi;
            #else
            //cccccccccccccccccccccccccc LINEAR VERSION
            #if EDDY_VISCOSITY_LINEAR || EDDY_VISCOSITY_PARABOLIC

            for(jphi=0;jphi<nInterfaces;jphi++){
                z_depth = hlayer/2 + hlayer*jphi;
                #if EDDY_VISCOSITY_LINEAR
                epsis2[jphi] = 10*k*ustar*(z_depth+5e-6);
                #endif 
                #if EDDY_VISCOSITY_PARABOLIC
                epsis2[jphi] = k*ustar*z_depth*(1-z_depth/arrays->h[idx]);
                #endif

            }

            Bi = 1.+dt/(hlayer*hlayer)*((epsis2[1]+epsis2[0])/2.);
            Ci = -((epsis2[1]+epsis2[0])/2.)*dt/(hlayer*hlayer);

            B[0] = Bi;
            C[0] = Ci;

            A[1] = -((epsis2[1]+epsis2[0])/2.)*dt/(hlayer*hlayer);
            B[1] = 1.+((epsis2[2]+2.*epsis2[1]+epsis2[0])/2.)*(dt/(hlayer*hlayer));
            C[1] = -((epsis2[2]+epsis2[1])/2.)*dt/(hlayer*hlayer);

            Bbis[0] = Bi;
            r[0] = arrays->phi[idx];
            r[1] = arrays->phi[ncells+idx];
            rbis[0] = r[0];
            rbis[1] = arrays->phi[ncells+idx] - (A[1]/Bi)*r[0];
            Bbis[1] = B[1]-(A[1]*Ci)/Bi;
            
            #else
            
            A = -epsis1*dt/(hlayer*hlayer);
            B = 1+2*epsis1*(dt/(hlayer*hlayer));
            C = -epsis1*dt/(hlayer*hlayer);

            Bi = 1+dt/(hlayer*hlayer)*epsis1;
            Ci = -(dt*epsis1)/(hlayer*hlayer);

            Af = dt/(hlayer*hlayer)*(-epsis1);
            Bf = 1+dt/(hlayer*hlayer)*epsis1;

            Bbis[0] = Bi;
            r[0] = arrays->phi[idx];
            r[1] = arrays->phi[ncells+idx];
            rbis[0] = r[0];
            rbis[1] = arrays->phi[ncells+idx] - (A/Bi)*r[0];
            Bbis[1] = B-(A*Ci)/Bi;

            #endif

            #endif


            for(jphi=2;jphi<nInterfaces;jphi++){
                sid1 = jphi*ncells+idx; 
                r[jphi] = arrays->phi[sid1];

                #if EDDY_VISCOSITY || EDDY_VISCOSITY_LINEAR || EDDY_VISCOSITY_PARABOLIC

                    A[jphi] = -((epsis2[jphi]+epsis2[jphi-1])/2.)*dt/(hlayer*hlayer);
                    B[jphi] = 1.+((epsis2[jphi+1]+2.*epsis2[jphi]+epsis2[jphi-1])/2.)*(dt/(hlayer*hlayer));
                    C[jphi] = -((epsis2[jphi+1]+epsis2[jphi])/2.)*dt/(hlayer*hlayer);

                    Bbis[jphi] = B[jphi] - (A[jphi]*C[jphi-1])/Bbis[jphi-1] ;
                    rbis[jphi] = r[jphi] - (A[jphi]/Bbis[jphi-1])*rbis[jphi-1];  

                #else
                        
                    Bbis[jphi] = B - (A*C)/Bbis[jphi-1] ;
                    rbis[jphi] = r[jphi]- (A/Bbis[jphi-1])*rbis[jphi-1];  

                    
                #endif
                    
            } 

            //celdas internas

            //Version 3
            #if EDDY_VISCOSITY || EDDY_VISCOSITY_LINEAR || EDDY_VISCOSITY_PARABOLIC

            Af = -((epsis2[nInterfaces-1]+epsis2[nInterfaces-2])/2.)*dt/(hlayer*hlayer);
            Bf = 1+((epsis2[nInterfaces-1]+epsis2[nInterfaces-2])/2.)*dt/(hlayer*hlayer);

            A[nInterfaces-1] = Af;
            B[nInterfaces-1] = Bf;

            Bbis[nInterfaces-1] = Bf-(Af*C[nInterfaces-2])/Bbis[nInterfaces-2];
            #else
            Bbis[nInterfaces-1] = Bf-(Af*C)/Bbis[nInterfaces-2];
            #endif

            rbis[nInterfaces-1] = r[nInterfaces-1]-(Af/Bbis[nInterfaces-2])*rbis[nInterfaces-2];
            arrays->phi[(nInterfaces-1)*ncells+idx] = rbis[nInterfaces-1]/Bbis[nInterfaces-1];
            

            for(jphi=nInterfaces-2; jphi>0; jphi--){
                sid1 = jphi*ncells+idx;
                sid2 = (jphi+1)*ncells+idx;
                
                #if EDDY_VISCOSITY || EDDY_VISCOSITY_LINEAR || EDDY_VISCOSITY_PARABOLIC
                arrays->phi[sid1] = (rbis[jphi] - C[jphi]*arrays->phi[sid2])/Bbis[jphi];
                #else
                arrays->phi[sid1] = (rbis[jphi] - C*arrays->phi[sid2])/Bbis[jphi];
                #endif
                
            }

            jphi = 0;
            sid1 = jphi*ncells+idx;
            sid2 = (jphi+1)*ncells+idx;
            arrays->phi[sid1] = (rbis[jphi] - Ci*arrays->phi[sid2])/Bbis[jphi];

            
            for(jphi=0; jphi<nInterfaces; jphi++){
                sid1 = jphi*ncells+idx;
                arrays->hphi[sid1] = arrays->phi[sid1]*hlayer;

            }

            for(jphi=0;jphi<nInterfaces-1;jphi++){

                sid1 = jphi*ncells+idx;
                sid2 = (jphi+1)*ncells+idx; //j+1

                phij1 = arrays->phi[sid1];
                phij2 = arrays->phi[sid2]; //j-1 aux2

                hphi1 = arrays->hphi[sid1];
                hphi2 = arrays->hphi[sid2];

                aux1 = dt*(ws*phij2);
                
                if((fabs(aux1)>TOL14)){

                    
                    if(aux1>0 && aux1>(arrays->hphi[sid2])){
                        aux1 = arrays->hphi[sid2];
                    }else if(aux1<0 && fabs(aux1)>(arrays->hphi[sid1])){
                        aux1 = arrays->hphi[sid1];
                    }

                    arrays->hphi[sid1] += aux1;
                    arrays->hphi[sid2] -= aux1;
               
                }

                if(arrays->hphi[sid2]<0.0){
                    arrays->hphi[sid2]=0.0;
                    arrays->phi[sid2]=0.0;  
                } 
                
            }

            for(jphi=0;jphi<nInterfaces;jphi++){
                sid1 = jphi*ncells+idx;
                arrays->phi[sid1] = arrays->hphi[sid1]/hlayer;

                if(std::isnan(arrays->phi[sid1])){
                    printf("cell %d layer %f\n",idx, jphi);
                    //printf("h %lf\n",hlayer);
                    printf("hphi %lf phi %lf hlayer %lf\n ",arrays->hphi[sid1], arrays->phi[sid1], hlayer);
                }
            
            }

        
        }
        else{
            // for(jphi=0;jphi<nInterfaces;jphi++){
            //     sid1 = jphi*ncells+idx;
            //     arrays->phi[sid1] = 0.0;
            //     arrays->hphi[sid1] = 0.0;
            
            // }
        }

    }


}



////////////////////////////////////////////////////
__global__ void g_initialize_particle_diffusion_delta(int nTasks, t_arrays *arrays){
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
__global__ void g_wall_particle_diffusion_calculus(int nTasks, t_arrays *arrays){
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
    int nSolutes= arrays->nSolutes;
    int nInterface = arrays->nSolutes;
    int nSediments = arrays->nSediments;
    int nParticles = arrays->nSolutes + nSediments;

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

        #if SET_SOLUTE_UNROLL==0  //compact 
        //wall index
        idx=arrays->actWalls[i];

        #elif SET_SOLUTE_UNROLL==1  //unroll  
        //solute index
        jphi=(int)(i/nActWalls);
        
        //wall index
        iactWall=i-jphi*nActWalls;
        idx=arrays->actWalls[iactWall];
        #endif         

        //cells index
        id1=arrays->idx1[idx];
        id2=arrays->idx2[idx];
    
        hL = arrays->h[id1];
        hR = arrays->h[id2];

        #if SET_MULTILAYER
        hL = arrays->h[id1]/nInterfaces;
        hR = arrays->h[id2]/nInterfaces;
        #endif

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


            #if SET_SOLUTE_UNROLL==0 || SET_SED_UNROLL==0  //compact 
            for(jphi=0;jphi<nParticles;jphi++){
            #endif

            //printf("jbefore %d\n", jphi);

                
                #if SET_SOLUTE || SET_SED
                
                typeDiff = arrays->typeDiff[0];
                #else
                typeDiff = CONSTANT_DIFF;
                #endif
            
                //ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc diffusion contribution calculation
                if(typeDiff==NON_DIFF){
                    kL = 0.0;
                    kT = 0.0;

                }else if(typeDiff==CONSTANT_DIFF){
                    #if SET_SOLUTE && SET_SED

                    //printf("jafter%d\n", nParticles);
                    if(jphi<arrays->nSolutes){
                        kL = arrays->k_xx[jphi];
                        kT = arrays->k_yy[jphi];
                        //printf("kL %lf kT %lf particles%d\n",kL, kT,jphi);
                    }else{
                        kL = arrays->ks_xx[jphi-nSolutes];
                        kT = arrays->ks_yy[jphi-nSolutes];
                        //printf("kL %lf kT %lf particles%d\n",arrays->ks_xx[jphi-nSolutes], arrays->ks_yy[jphi-nSolutes],jphi);
                    }
                    //printf("kL %lf kT %lf particles%d\n",arrays->ks_xx[jphi], arrays->ks_yy[jphi],nSolutes);
                    #else
                    #if SET_SOLUTE
                    if(jphi<nParticles){
                        kL = arrays->k_xx[jphi];
                        kT = arrays->k_yy[jphi];
                    }
                    #endif

                    #if SET_SED
                    if(jphi<nParticles){
                        kL = arrays->ks_xx[jphi];
                        kT = arrays->ks_yy[jphi];
                    }
                    #endif
                    
                    #endif

                }else if(typeDiff==ANISOTROPIC_DIFF){
                    ustar = sqrt(gp*nman2wall*modU2/cbrt(hbar));
                    #if SET_SOLUTE && SET_SED
                    if(jphi<arrays->nSolutes){
                        kL = arrays->k_xx[jphi]*hbar*ustar;
                        kT = arrays->k_yy[jphi]*hbar*ustar;
                    }
                    if((jphi>=nSolutes) || (jphi<nParticles)){
                        kL = arrays->ks_xx[jphi]*hbar*ustar;
                        kT = arrays->ks_yy[jphi]*hbar*ustar;
                        //printf("kL %lf kT %lf j%d\n",arrays->ks_xx[jphi], arrays->ks_yy[jphi],jphi);
                    }
                    #else
                    #if SET_SOLUTE
                    if(jphi<nParticles){
                        kL = arrays->k_xx[jphi]*hbar*ustar;
                        kT = arrays->k_yy[jphi]*hbar*ustar;
                    }
                    #endif

                    #if SET_SED
                    if(jphi<nParticles){
                        kL = arrays->ks_xx[jphi]*hbar*ustar;
                        kT = arrays->ks_yy[jphi]*hbar*ustar;
                    }
                    #endif

                    #endif

                }
                //printf("kL %lf kT %lf j%d\n",arrays->ks_xx[jphi], arrays->ks_yy[jphi],jphi);
    
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

            #if SET_SOLUTE_UNROLL==0  || SET_SED_UNROLL==0//compact 
            }
            #endif

        } //end if(hL>minh && hR>minh){ //wet-wet walls

	} // end iwall loop

}


////////////////////////////////////////////////////
__global__ void g_update_particle_diffusion_contributions(int nTasks, t_arrays *arrays, double *localDtd){
/*----------------------------*/	
    int idx;
    int siw0,sid;

    int ncells=arrays->ncells;
    int NCwall=arrays->NCwall;

    int nActCells=arrays->nActCells;
    int jphi;
    int iactCell;
    
    int nParticles = arrays->nSolutes + arrays->nSediments;

    double Dtd;
	double total;

    Dtd = 10.*arrays->dt;
    
    int i = threadIdx.x+(blockIdx.x*blockDim.x);
    if(i<nTasks){

        #if SET_SOLUTE_UNROLL==0  //compact 
        //cell index
        idx=arrays->actCells[i];

        #elif SET_SOLUTE_UNROLL==1  //unroll  
        //solute index
        jphi=(int)(i/nActCells);

        //cell index
        iactCell=i-jphi*nActCells;
        idx=arrays->actCells[iactCell]; 
        #endif           
        
  
        #if SET_SOLUTE_UNROLL==0  //compact 
        for(jphi=0;jphi<nParticles;jphi++){
        #endif
            
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

        #if SET_SOLUTE_UNROLL==0  //compact 
        }
        #endif

    }

}


//////////////////////////////////////////////////////
__global__ void g_get_particle_diffusion_dtmin(t_arrays *arrays, double *localDtd, int *idmin){
/*----------------------------*/	
	double minDtd;

	minDtd=localDtd[(*idmin)-1];

    //storage time step
    arrays->Dtd=minDtd;
    

}


//////////////////////////////////////////////////////
__global__ void g_update_particle_diffusion_cells(int nTasks, t_arrays *arrays, double *Dtd){
/*----------------------------*/	

    int idx,neighid;
    int siw,sid;
    int k;
    int ncells=arrays->ncells;
    int NCwall=arrays->NCwall;
    int nw_calc=arrays->nw_calc;
    int nInterfaces= arrays->nSolutes;
    int nParticles = arrays->nSolutes + arrays->nSediments;

    double dtd;
	double total;
    double contrib1, contrib2;

    int nActCells=arrays->nActCells;
    int jphi;
    int iactCell;
    
    int i = threadIdx.x+(blockIdx.x*blockDim.x);
    if(i<nTasks){

        #if SET_SOLUTE_UNROLL==0  //compact 
        //cell index
        idx=arrays->actCells[i];

        #elif SET_SOLUTE_UNROLL==1  //unroll  
        //solute index
        jphi=(int)(i/nActCells);

        //cell index
        iactCell=i-jphi*nActCells;
        idx=arrays->actCells[iactCell]; 
        #endif           
        
  
        #if SET_SOLUTE_UNROLL==0  //compact 
        for(jphi=0;jphi<nParticles;jphi++){
        #endif

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

            #if SET_MULTILAYER
            arrays->hphi[sid]=arrays->phi[sid]*(arrays->h[idx]/nSolutes); 
            #endif

        #if SET_SOLUTE_UNROLL==0  //compact 
        }
        #endif
        

    }

}



#endif  