#include "sediment.cuh"

#if SET_SED

////////////////////////////////////////////////////
__global__ void g_initialize_sediment_erosion_delta(int nTasks, t_arrays *arrays){
/*----------------------------*/
    int k, jphi, idx, sid;
    int ncells = arrays->ncells;
    int NCwall = arrays->NCwall;

    int i = threadIdx.x+(blockIdx.x*blockDim.x);    
    if(i<nTasks){
        //solute index
        idx=arrays->actCells[i];
        arrays->Nb[idx] = 0.0;
        arrays->phiZero[idx] += arrays->phi[sid];

        //cell index
        for(jphi=0;jphi<arrays->nSediments;jphi++){
            sid = (arrays->nSolutes + jphi)*ncells+idx;
            arrays->Ns[sid] = 0.0;

        }
    }
}

////////////////////////////////////////////////////
__global__ void g_cell_sediment_Erosion_calculus(int nTasks, t_arrays *arrays){
/*----------------------------*/
    int idx;
    int sid;
    int jphi;
    int nActCells=arrays->nActCells;
    int iactCell;

    int ncells=arrays->ncells;
    int nSediments = arrays->nSediments;
    int nSolutes = arrays->nSolutes;
    int nParticles = nSediments + nSolutes;

    double aux1, aux2, aux3, aux4;

    double Ns;
    double Nb;
    double Ebj, Dbj;

    double pd;
    double WsFF, WsFs;
    double Wsmp;
    double EquConcFF, EquConcFs;
    double Csst;
    double dsp;
    double rhoS, rhoSW;
    double Fsp;
    double TdFsp;
    double Css;
    double fAngle;
    double ks_xx, ks_yy;
    double thob;
    double phiZero;
    double nman;
    double SsModulus;
    double Theta, Thetar;
    double h;
    double u,v;
    double moduloU;

    double dt=arrays->dt;

    int i = threadIdx.x+(blockIdx.x*blockDim.x);    
    if(i<nTasks){

        for(jphi=0;jphi<arrays->nSediments;jphi++){
            if(arrays->dsp[jphi] > TOL12){
                TdFsp += arrays->Fsp[jphi]/arrays->dsp[jphi];      
            }else{
                TdFsp += 0.0;
            }
        }
        
        #if SET_SOLUTE_UNROLL==0  //compact
        //wall index
        idx=arrays->actCells[i];

        #elif SET_SOLUTE_UNROLL==1  //unroll 
        //solute index
        jphi=(int)(i/nActCells);
        
        //wall index
        iactCell=i-(nSolutes+jphi)*nActCells;
        idx=arrays->actCells[iactCell];
        #endif

        h = arrays->h[idx];
        u = arrays->u[idx];
        v = arrays->v[idx];
                
        moduloU = sqrt(u*u + v*v);

        phiZero = arrays->phiZero[idx];
        nman = arrays->nman[idx];

        rhoS = arrays->rhoS;

        if(arrays->h[idx]>TOL12){ //wet cells

            #if SET_SOLUTE_UNROLL==0  //compact
            for(jphi=0;jphi<arrays->nSediments;jphi++){
            #endif

                sid = (nSolutes + jphi)*ncells+idx;
        
                Fsp = arrays->Fsp[jphi];
                Css = arrays->Css[jphi];
                dsp = arrays->dsp[jphi];
                EquConcFF = arrays->EquConcFF[jphi];


                rhoSW = rhoS/_rhow_;

                WsFs = arrays->WsFs[jphi];

                if(idx == 1035457){
                    printf("Fsp %.12lf Css %.12lf dsp %.12lf jphi %f\n", Fsp, Css, dsp, jphi);
                    printf("EquConcFF %.12lf WsFs %.12lf jphi %f\n", EquConcFF, WsFs, jphi);
                } 

                SsModulus = _rhow_*arrays->h[idx]*(nman*nman*(u*u + v*v))/(pow(arrays->h[idx],4./3.));
                Theta = abs(SsModulus)/((rhoS-_rhow_)*_g_*dsp);

                Thetar = Theta/Css; 
               
                if(WsFs >TOL12){
                    if(arrays->EquConcF==EQUCONCF_BAGNOLD){
                        aux1 = (1./(h*moduloU));
                        Csst = 0.01*(rhoS/(rhoS-_rhow_)*((thob*moduloU*moduloU)/WsFs));
                        EquConcFs = aux1*Csst;

                    }else if(arrays->EquConcF==EQUCONCF_WU){
                        aux1=1./21.1*pow(dsp,1./6.);
                        aux1 = aux1/nman;
                        aux2 = sqrt(aux1*aux1*aux1);
                        if(aux2<1.){
                            aux2 = 1.;
                        }
                        aux3 = aux2*Thetar-1.;
                        aux4 = (Thetar-1.)*moduloU/WsFs;

                        aux1 = 0.0053*pow(aux3,2.2)+0.0000262*pow(aux4,1.74);
                        aux2 = sqrt((rhoS/(rhoS-_rhow_))*_g_*dsp*dsp*dsp);

                        EquConcFs = aux2*aux1;
                    }
                }else{
                    EquConcFs = 0.0;
                }

                

                //Calculation of the settling velocity in the mixture
                //aux1 = EquConcFF*WsFs;

                // if(phiZero > 2.*sqrt(dsp)){
                //     if(dsp>TOL12){
                //         Wsmp = aux1 * (1-(phiZero/(2*sqrt(dsp))))*(1-(phiZero/(2*sqrt(dsp))))*(1-(phiZero/(2*sqrt(dsp))));
                //     }else{
                //         Wsmp = 0.0;
                //     }

                // }else{
                //     Wsmp = aux1*(1-phiZero)*(1-phiZero)*(1-phiZero)*(1-phiZero);
                // }

                //EROSION
                Ebj = EquConcFF*WsFs*EquConcFs;
                //DEPOSITION
                Dbj = WsFs*WsFF*arrays->phi[sid];
                
                arrays->Ns[sid] = (Ebj - Dbj);
                arrays->Nb[idx] += arrays->Ns[sid];

                if(idx == 1035457){
                    printf("Ns %.12lf Nb %.12lf sid %f\n", Ns, Nb, sid);
                } 


            #if SET_SOLUTE_UNROLL==0  //compact
            }
            #endif


        }    
    
    
    }

}

////////////////////////////////////////////////////
__global__ void g_update_sediment_erosion_cells(int nTasks, t_arrays *arrays){
/*----------------------------*/
int idx;
int sid;
int ncells=arrays->ncells;
int NCwall=arrays->NCwall;

int nActCells=arrays->nActCells;
int jphi;
int iactCell;

double dt;

double pd;
double EtaS;
double rhob, rhoBulk;
double rhoS;
double phiZero;

double bedExchangePos;
double bedExchangeNeg;
double aux1, aux2, aux3, aux4, aux5;
double bedExchange;
double mod_EtaS, EtaS_eff;
double cr;
double deltaz;


int i = threadIdx.x+(blockIdx.x*blockDim.x);
if(i<nTasks){

    pd = arrays->pd;
    rhoS = arrays->rhoS;
    dt=arrays->dt;

    if(pd<1.){
        EtaS = 1./(1-pd);
    }else{
        EtaS = 0;
    }
    

    
    //cell index
    idx=arrays->actCells[i];

    phiZero = arrays->phiZero[idx];

    rhob = _rhow_*pd + rhoS*(1-pd);
    rhoBulk = _rhow_*(1-phiZero) + rhoS*phiZero; 

    bedExchange = 0.0;

    if(arrays->h[idx]>TOL12){ //wet cells

        aux1 = 0.0;
        aux2 = 0.0;

        for(jphi=0;jphi<arrays->nSediments;jphi++){
            sid = (arrays->nSolutes + jphi)*ncells+idx;
            aux1 += arrays->hphi[sid]/arrays->h[idx];
            aux2 += (arrays->hphi[sid]/arrays->h[idx])*EtaS;
        }

        mod_EtaS = 0.0;

        if(aux2>1.0){
            mod_EtaS = 1./aux1;
        }

        aux2 = 0.0;
        for(jphi=0;jphi<arrays->nSediments;jphi++){
            sid = (arrays->nSolutes + jphi)*ncells+idx;
            if(mod_EtaS>0.0){
                EtaS_eff = mod_EtaS;
            }else{
                EtaS_eff = EtaS;
            }
            aux1 = arrays->Ns[sid]*dt;

            if(aux1<0.0){
                aux1 = fmax(-1.*arrays->hphi[sid], aux1);
            }

            aux2 += aux1*EtaS_eff;
        }

        cr = 1.0;

        if(abs(aux2)>=TOL12){
            aux3 = aux2;
            if(aux3<0.0){
                aux3 = fmax(-1.*arrays->h[idx], aux3);
                aux3 = fmin(aux3,0.0);
            }else{
                // aux3 = fmin((arrays->z[idx]-1e6), aux3);
                // aux3 = fmax(aux3,0.0); 
            }
            cr = aux3/aux2;
        }else{
            cr=0.0;
        }

        for(jphi=0;jphi<arrays->nSediments;jphi++){
            sid = (arrays->nSolutes + jphi)*ncells+idx;
            if(mod_EtaS>0.0){
                EtaS_eff = mod_EtaS;
            }else{
                EtaS_eff = EtaS;
            }

            aux1 = arrays->Ns[sid]*dt;

            if(aux1<0.0){
                aux1 = fmax(-1.*arrays->hphi[sid], aux1);
            }

            aux1 = cr*aux1;
            arrays->z[idx] += -aux1*(EtaS_eff);
            // if(arrays->z[idx] < 1e-6){
            //     arrays->z[idx] = 1e-6;
            // }

            arrays->hphi[sid] += aux1;
            if(arrays->hphi[sid]<0.0){
                arrays->hphi[sid]=0.0;
            }

            arrays->h[idx] += aux1*EtaS_eff;
            if(arrays->h[idx]<0.0){
                arrays->h[idx] = 0.0;
                arrays->hphi[sid] = 0.0;
            }

            arrays->phi[sid] = arrays->hphi[sid]/arrays->h[idx];

            if(idx == 1035457){
                printf("hphi %.12lf phi %.12lf\n", arrays->hphi[sid], arrays->phi[sid]);
                printf("porosityCoef %.12lf EtaS %.12lf rhoS %.12lf\n", pd, EtaS, rhoS);
            } 

            if(std::isnan(arrays->phi[sid])){
                printf("cell %d layer %f\n",idx, jphi);
                //printf("h %lf\n",hlayer);
                printf("hphi %lf phi %lf hlayer %lf\n ",arrays->hphi[sid], arrays->phi[sid], jphi);
            }

        }
        
        
    }
 
}
}

#endif 