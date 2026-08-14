#include "sediment.cuh"

#if SET_SED

////////////////////////////////////////////////////
__global__ void g_initialize_sediment_erosion_delta(int nTasks, t_arrays *arrays){
/*----------------------------*/
    int k, jphi, idx, sid, sid1;
    int ncells = arrays->ncells;
    int NCwall = arrays->NCwall;

    int i = threadIdx.x+(blockIdx.x*blockDim.x);    
    if(i<nTasks){
        //solute index
        idx=arrays->actCells[i];
        arrays->Nb[idx] = 0.0;
        arrays->phiZero[idx] = 0.0;

        //cell index
        
        #if SET_MULTILAYER_SED
        sid = (arrays->nSolutes)*ncells+idx;
        sid1 = 0;
        arrays->Ns[sid1] = 0.0;
        arrays->phiZero[idx] += arrays->phi[sid];
        #else  
        for(jphi=0;jphi<arrays->nSediments;jphi++){
            sid = (arrays->nSolutes + jphi)*ncells+idx;
            sid1 = jphi*ncells+idx;
            arrays->Ns[sid1] = 0.0;
            arrays->phiZero[idx] += arrays->phi[sid];
        } 
        #endif 

        
    }
}

////////////////////////////////////////////////////
__global__ void g_cell_sediment_Erosion_calculus(int nTasks, t_arrays *arrays){
/*----------------------------*/
    int idx;
    int sid, sid1;
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
    double EquConcFF;
    double EquConcFs;
    double Csst;
    double dsp;
    double rhoS, rhoSW;
    double rhob, rhoBulk;
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
    double u_comp, v_comp;
    double Ssu, Ssv;

    double dt=arrays->dt;

    int i = threadIdx.x+(blockIdx.x*blockDim.x);  

    if(i<nTasks){

        #if SET_MULTILAYER_SED
            TdFsp = arrays->Fsp[0]/arrays->dsp[0]; 
        #else  
        for(jphi=0;jphi<arrays->nSediments;jphi++){
            if(arrays->dsp[jphi] > TOL12){
                TdFsp += arrays->Fsp[jphi]/arrays->dsp[jphi];      
            }else{
                TdFsp += 0.0;
            }
        }
        #endif 

        
        //compact
        //wall index
        idx=arrays->actCells[i];

        #if SET_MULTILAYER_SED 
            h = arrays->h[idx]/nSediments;
        #else 
            h = arrays->h[idx];
        #endif 
        u = arrays->u[idx];
        v = arrays->v[idx];
                
        moduloU = sqrt(u*u + v*v);

        phiZero = arrays->phiZero[idx];
        nman = arrays->nman[idx];

        rhoS = arrays->rhoS;

        if(arrays->h[idx] >= arrays->minh){ //wet cells

            #if SET_MULTILAYER_SED
            sid = (nSolutes)*ncells+idx;
            sid1 = 0*ncells+idx;

            Fsp = arrays->Fsp[0];
            Css = arrays->Css[0];
            dsp = arrays->dsp[0];
            EquConcFF = arrays->EquConcFF[0];

            pd = arrays->pd;

            rhob = _rhow_*pd + rhoS*(1-pd);
            rhoBulk = _rhow_*(1-phiZero) + rhoS*phiZero; 


            rhoSW = rhoS/_rhow_;

            WsFs = arrays->WsFs[0];

            SsModulus = rhoBulk*_g_*arrays->h[idx] *(nman*nman*(u*u + v*v))/(pow(arrays->h[idx],4./3.));
            Theta = SsModulus/((rhoS-_rhow_)*_g_*dsp);

            if(Css > 0.0){
                    Thetar = Theta/Css;
                }else{
                    Thetar = 0.0; 
                }
            
            if(WsFs >TOL12){
                    if(arrays->EquConcF==EQUCONCF_BAGNOLD){
                       //if(moduloU > 0.0){
                        //    aux1 = (1./(h*moduloU));
                        //}else{
                        //    aux1 = 0.0;
                        //}
                        
                        //Csst = 0.01*((rhoS/(rhoS-_rhow_))*((SsModulus*moduloU*moduloU)/WsFs));
                        EquConcFs = 0.01*SsModulus*moduloU/((rhoS-_rhow_)*_g_*arrays->h[idx]*WsFs);
                        
                        aux1 = arrays->h[idx]*moduloU*(1-pd);

                        if (EquConcFs>aux1){
                            EquConcFs = aux1;
                        }

                        
                    }else if(arrays->EquConcF==EQUCONCF_WU){
                        aux1=(1./21.1)*pow(dsp,1./6.);
                        aux1 = aux1/nman;
                        aux2 = sqrt(aux1*aux1*aux1);
                        if(aux2<1.){
                            aux2 = 1.;
                        }
                        aux3 = aux2*Thetar-1.;
                        if(WsFs >0.0){
                            aux4 = (Thetar-1.)*moduloU/WsFs;
                        }else{
                            aux4 = 0.0;
                        }
                        if(aux4 >1 && aux3 > 1){
                            aux1 = 0.0053*pow(aux3,2.2)+0.0000262*pow(aux4,1.74);
                        }else if(aux4 < 1 && aux3 >1){
                            aux1 = 0.0053*pow(aux3,2.2);
                        }else if(aux4 > 1 && aux3 <1){
                            aux1 = 0.0000262*pow(aux4,1.74);
                        }
                        aux2 = sqrt((rhoS/(rhoS -_rhow_))*_g_*dsp*dsp*dsp);

                        EquConcFs = aux2*aux1;

                        aux1 = h*moduloU*(1-pd);

                        if (EquConcFs>aux1){
                            EquConcFs = aux1;
                        }
                    }
                }else{
                    EquConcFs = 0.0;
                }

             //EROSION
            Ebj = EquConcFF*Fsp*WsFs*EquConcFs;
            //DEPOSITION
            Dbj = arrays->WsFs[jphi]*arrays->WsFF[jphi]*arrays->phi[sid];

            arrays->Ns[sid1] = (Ebj - Dbj);
            arrays->Nb[idx] += arrays->Ns[sid1];

            //printf(" Ns %lf\n", arrays->Ns[sid1]);

            

            #else
            for(jphi=0;jphi<arrays->nSediments;jphi++){
             
                sid = (nSolutes + jphi)*ncells+idx;
                
                sid1 = jphi*ncells+idx;
        
                Fsp = arrays->Fsp[jphi];
                Css = arrays->Css[jphi];
                dsp = arrays->dsp[jphi];
                EquConcFF = arrays->EquConcFF[jphi];

                pd = arrays->pd;

                rhob = _rhow_*pd + rhoS*(1-pd);
                rhoBulk = _rhow_*(1-phiZero) + rhoS*phiZero; 

                rhoSW = rhoS/_rhow_;

                WsFs = arrays->WsFs[jphi];

                //printf("EquConcFF %.12lf WsFs %.12lf jphi %d\n", arrays->EquConcFF[jphi], WsFs, jphi);
            
                //SsModulus = _rhow_*aux2*moduloU;

                SsModulus = rhoBulk*_g_*h*(nman*nman*(u*u + v*v))/(pow(h,4./3.));

                //printf("SsModulus %.12lf rhoS %.12lf\n", SsModulus, rhoS);
                Theta = SsModulus/((rhoS-_rhow_)*_g_*dsp);

                if(Css > 0.0){
                    Thetar = Theta/Css;
                }else{
                    Thetar = 0.0; 
                }
               

                if(WsFs >TOL12){
                    if(arrays->EquConcF==EQUCONCF_BAGNOLD){
                       //if(moduloU > 0.0){
                        //    aux1 = (1./(h*moduloU));
                        //}else{
                        //    aux1 = 0.0;
                        //}
                        
                        //Csst = 0.01*((rhoS/(rhoS-_rhow_))*((SsModulus*moduloU*moduloU)/WsFs));
                        EquConcFs = 0.01*SsModulus*moduloU/((rhoS-_rhow_)*_g_*h*WsFs);
                        
                        aux1 = h*moduloU*(1-pd);

                        if (EquConcFs>aux1){
                            EquConcFs = aux1;
                        }

                        
                    }else if(arrays->EquConcF==EQUCONCF_WU){
                        aux1=(1./21.1)*pow(dsp,1./6.);
                        aux1 = aux1/nman;
                        aux2 = sqrt(aux1*aux1*aux1);
                        if(aux2<1.){
                            aux2 = 1.;
                        }
                        aux3 = aux2*Thetar-1.;
                        if(WsFs >0.0){
                            aux4 = (Thetar-1.)*moduloU/WsFs;
                        }else{
                            aux4 = 0.0;
                        }
                        if(aux4 >1 && aux3 > 1){
                            aux1 = 0.0053*pow(aux3,2.2)+0.0000262*pow(aux4,1.74);
                        }else if(aux4 < 1 && aux3 >1){
                            aux1 = 0.0053*pow(aux3,2.2);
                        }else if(aux4 > 1 && aux3 <1){
                            aux1 = 0.0000262*pow(aux4,1.74);
                        }
                        aux2 = sqrt((rhoS/(rhoS -_rhow_))*_g_*dsp*dsp*dsp);

                        EquConcFs = aux2*aux1;

                        aux1 = h*moduloU*(1-pd);

                        if (EquConcFs>aux1){
                            EquConcFs = aux1;
                        }
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
                Ebj = EquConcFF*Fsp*WsFs*EquConcFs;
                //DEPOSITION
                Dbj = arrays->WsFs[jphi]*arrays->WsFF[jphi]*arrays->phi[sid];

                arrays->Ns[sid1] = (Ebj - Dbj);
                arrays->Nb[idx] += arrays->Ns[sid1];

                //printf(" Ns %lf\n", arrays->Ns[sid1]);
            
            }
        #endif 

        }else{
            #if SET_MULTILAYER_SED
            sid1 = 0*ncells+idx;
            arrays->Ns[sid1] = 0.0;
        
            arrays->Nb[idx] = 0.0;
            #else  
            for(jphi=0;jphi<arrays->nSediments;jphi++){
                sid1 = jphi*ncells+idx;
                arrays->Ns[sid1] = 0.0;
            }

            arrays->Nb[idx] = 0.0;
            #endif 
        }    
    
    
    }

}

////////////////////////////////////////////////////
__global__ void g_update_sediment_erosion_cells(int nTasks, t_arrays *arrays){
/*----------------------------*/
int idx;
int sid, sid1;
int ncells=arrays->ncells;
int NCwall=arrays->NCwall;
int nSediments=arrays->nSediments;
int nSolutes = arrays->nSolutes;

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

double hlayer;


int i = threadIdx.x+(blockIdx.x*blockDim.x);

    if(i<nTasks){

        pd = arrays->pd;
        rhoS = arrays->rhoS;
        dt=arrays->dt;

        if(pd<1.){
            EtaS = 1./(1-pd);
            //printf("EtaS %lf\n", EtaS);
        }else{
            EtaS = 0;
        }
        
        
        //cell index
        idx=arrays->actCells[i];

        phiZero = arrays->phiZero[idx];

        rhob = _rhow_*pd + rhoS*(1-pd);
        rhoBulk = _rhow_*(1-phiZero) + rhoS*phiZero; 

        bedExchange = 0.0;

        if (idx == 109075){ 
        printf("phi %lf hphi %lf\n", arrays->phi[nSolutes*ncells+idx], arrays->hphi[nSolutes*ncells+idx]);
        } 
        if(arrays->h[idx] > arrays->minh){ //wet cells

            aux1 = 0.0;
            aux2 = 0.0;

            #if SET_MULTILAYER_SED 
            
            // sid = (arrays->nSolutes)*ncells+idx;
            // aux1 += arrays->hphi[sid]/arrays->h[idx]/nSediments;
            // aux2 += (arrays->hphi[sid]/arrays->h[idx]/nSediments)*EtaS;

            //for(jphi=0;jphi<arrays->nSediments;jphi++){
                sid = (arrays->nSolutes)*ncells+idx;
                aux1 += arrays->hphi[sid]/arrays->h[idx];
                aux2 += (arrays->hphi[sid]/arrays->h[idx]);
            //}
            
            #else
            for(jphi=0;jphi<arrays->nSediments;jphi++){
                    sid = (arrays->nSolutes + jphi)*ncells+idx;
                    aux1 += arrays->hphi[sid]/arrays->h[idx];
                    aux2 += (arrays->hphi[sid]/arrays->h[idx])*EtaS;
                }
            #endif 
            if (idx == 109075){ 
                printf("hphi %lf\n", arrays->hphi[sid]);
            }

            //printf("aux1 %lf aux2 %lf\n", aux1, aux2);

            mod_EtaS = 0.0;

            if(aux2>1.0){
                if(aux1>0.0){
                    mod_EtaS = 1./aux1;
                }
            }

            aux2 = 0.0;

            #if SET_MULTILAYER_SED
            
            // sid = (arrays->nSolutes)*ncells+idx;
            // sid1 = 0;


            // if(mod_EtaS>0.0){
            //     EtaS_eff = mod_EtaS;
            // }else{
            //     EtaS_eff = EtaS;
            // }
            // aux1 = arrays->Ns[sid1]*dt;

            // if(aux1<0.0){
            //     aux1 = fmax(-1.*arrays->hphi[sid]/nSediments, aux1);
            // }

            // aux2 += aux1*EtaS_eff;

            //for(jphi=0;jphi<arrays->nSediments;jphi++){
                sid = (arrays->nSolutes)*ncells+idx;
                sid1 = 0*ncells+idx;


                if(mod_EtaS>0.0){
                    EtaS_eff = mod_EtaS;
                }else{
                    EtaS_eff = EtaS;
                }
                aux1 = arrays->Ns[sid1]*dt;

                if(aux1<0.0){
                    aux1 = fmax(-1.*arrays->hphi[sid], aux1);
                }

                aux2 += aux1*EtaS_eff;
            //}
            
            #else  
            for(jphi=0;jphi<arrays->nSediments;jphi++){
                sid = (arrays->nSolutes + jphi)*ncells+idx;
                sid1 = jphi*ncells+idx;


                if(mod_EtaS>0.0){
                    EtaS_eff = mod_EtaS;
                }else{
                    EtaS_eff = EtaS;
                }
                aux1 = arrays->Ns[sid1]*dt;

                if(aux1<0.0){
                    aux1 = fmax(-1.*arrays->hphi[sid], aux1);
                }

                aux2 += aux1*EtaS_eff;
            }
            #endif 
            cr = 1.0;

            if(abs(aux2)>=TOL12){
                aux3 = aux2;
                if(aux3<0.0){ //deposition limited by h
                    #if SET_MULTILAYER_SED
                    aux3 = fmax(-1.*arrays->h[idx]/nSediments, aux3);
                    #else  
                    aux3 = fmax(-1.*arrays->h[idx], aux3);
                    #endif 
                    aux3 = fmin(aux3,0.0);
                }else{ // erosion limited by z
                    // aux3 = fmin((arrays->z[idx]-arrays->maxZ), aux3);
                    // //printf("z %.12lf zmax %.12lf", arrays->z[idx], arrays->maxZ);
                    // aux3 = fmax(aux3,0.0); 
                }
                cr = aux3/aux2; // 0<cr<1
            }else{
                cr=0.0;
            }

            #if SET_MULTILAYER_SED
            // sid = (arrays->nSolutes)*ncells+idx;
            // sid1 = 0;

            // if(mod_EtaS>0.0){
            //     EtaS_eff = mod_EtaS;
            // }else{
            //     EtaS_eff = EtaS;
            // }

            // aux1 = arrays->Ns[sid1]*dt;


            // if(aux1<0.0){
            //     aux1 = fmax(-1.*arrays->hphi[sid]/nSediments, aux1);
            // }

            // aux1 = cr*aux1;

            // arrays->z[idx] += -aux1*EtaS_eff;

            // arrays->hphi[sid] += aux1;

            //     if(arrays->hphi[sid]<TOL12){
            //         arrays->hphi[sid] = 0.0;
            //         arrays->phi[sid] = 0.0;
            //     }

            // arrays->h[idx] += aux1*EtaS_eff;

            //     if(arrays->h[idx] < 0.0){
            //         arrays->h[idx] = 0.0;
            //         arrays->hphi[sid] = 0.0;
            //         arrays->phi[sid] = 0.0;
            //         arrays->u[idx] = 0.0;
            //         arrays->v[idx] = 0.0;

            //     }
            
            // if(arrays->h[idx]>TOL12 && arrays->hphi[sid]>0.0){
            //     arrays->phi[sid] = arrays->hphi[sid]/(arrays->h[idx]/nSediments);
            // }else{
            //     arrays->phi[sid] = 0.0;
            // }

            //for(jphi=0;jphi<arrays->nSediments;jphi++){
                sid = (nSolutes)*ncells+idx;
                sid1 = 0*ncells+idx;

                if(mod_EtaS>0.0){
                    EtaS_eff = mod_EtaS;
                }else{
                    EtaS_eff = EtaS;
                }

                aux1 = arrays->Ns[sid1]*dt;


                if(aux1<0.0){
                    aux1 = fmax(-1.*arrays->hphi[sid], aux1);
                }
 
                aux1 = cr*aux1;

                arrays->z[idx] += -aux1*EtaS_eff;

                arrays->hphi[sid] += aux1;

                if(arrays->hphi[sid]<TOL12){
                    arrays->hphi[sid] = 0.0;
                    arrays->phi[sid] = 0.0;
                }
               
                arrays->h[idx] += aux1*nSediments;

                if (idx == 109075){ 
                printf("h %lf\n", arrays->h[idx]);
                } 

                if(arrays->h[idx] < 0.0){
                    arrays->h[idx] = 0.0;
                    arrays->hphi[sid] = 0.0;
                    arrays->phi[sid] = 0.0;
                    arrays->u[idx] = 0.0;
                    arrays->v[idx] = 0.0;
                }
                
                if(arrays->h[idx]>TOL12 && arrays->hphi[sid]>0.0){
                    arrays->phi[sid] = arrays->hphi[sid]/(arrays->h[idx]/nSediments);
                }else{
                    arrays->phi[sid] = 0.0;
                }
                
            //}

            #else   
            for(jphi=0;jphi<arrays->nSediments;jphi++){
                sid = (arrays->nSolutes + jphi)*ncells+idx;
                sid1 = jphi*ncells+idx;

                if(mod_EtaS>0.0){
                    EtaS_eff = mod_EtaS;
                }else{
                    EtaS_eff = EtaS;
                }

                aux1 = arrays->Ns[sid1]*dt;


                if(aux1<0.0){
                    aux1 = fmax(-1.*arrays->hphi[sid], aux1);
                }
 
                aux1 = cr*aux1;

                arrays->z[idx] += -aux1*EtaS_eff;

                 // if(abs(arrays->z[idx])> arrays->maxZ){
                //     //printf("z %.12lf zmax %.12lf", arrays->z[idx], arrays->maxZ);

                //     if(arrays->z[idx]< 0.0){
                //         arrays->z[idx] = -arrays->maxZ;
                //     }else{
                //         arrays->z[idx] = arrays->maxZ;
                //     }

                arrays->hphi[sid] += aux1;

                if(arrays->hphi[sid]<TOL12){
                    arrays->hphi[sid] = 0.0;
                    arrays->phi[sid] = 0.0;
                }

                // if(abs(arrays->z[idx])> arrays->maxZ){
                //     //printf("z %.12lf zmax %.12lf", arrays->z[idx], arrays->maxZ);

                //     if(arrays->z[idx]< 0.0){
                //         arrays->z[idx] = -arrays->maxZ;
                //     }else{
                //         arrays->z[idx] = arrays->maxZ;
                //     }

                // }
                arrays->h[idx] += aux1*EtaS_eff;

                if(arrays->h[idx] < 0.0){
                    arrays->h[idx] = 0.0;
                    arrays->hphi[sid] = 0.0;
                    arrays->phi[sid] = 0.0;
                    arrays->u[idx] = 0.0;
                    arrays->v[idx] = 0.0;

                }


                // if(arrays->h[idx] > 0.0){

                //     // if(rhoBulk > 0.0){
                //     //     aux1 = arrays->u[idx]*((rhob/rhoBulk)-1);
                //     //     aux2 = arrays->v[idx]*((rhob/rhoBulk)-1);
                        
                //     // }else{
                //     //     aux1 = -arrays->u[idx];
                //     //     aux2 = -arrays->v[idx];
                //     // }
                //     // arrays->u[idx] += -(1./arrays->h[idx])*aux1*aux3*EtaS_eff;
                //     // arrays->v[idx] += -(1./arrays->h[idx])*aux2*aux3*EtaS_eff;

                // }else{
                //     arrays->u[idx] = 0.0;
                //     arrays->v[idx] = 0.0;
                //     arrays->hphi[sid] = 0.0;
                //     arrays->phi[sid] = 0.0;
                // }

            
                // if(idx == 104477){
                //     printf("u %.12lf v %.12lf\n", arrays->u[idx], arrays->v[idx]);
                //     printf("porosityCoef %.12lf EtaS %.12lf rhoS %.12lf\n", pd, EtaS, rhoS);
                // } 
                
                if(arrays->h[idx]>TOL12 && arrays->hphi[sid]>0.0){
                    arrays->phi[sid] = arrays->hphi[sid]/arrays->h[idx];
                }else{
                    arrays->phi[sid] = 0.0;
                }
                
            }
            #endif
            

            // if(idx == 88511){
            //printf("hphi %.12lf phi %.12lf z %.12lf h %.12lf\n", arrays->hphi[sid], arrays->phi[sid], arrays->z[idx], arrays->h[idx]);
            //     //printf("porosityCoef %.12lf EtaS %.12lf rhoS %.12lf\n", pd, EtaS, rhoS);
            // } 

            
            
        }
    
    }
}

////////////////////////////////////////////////////
__global__ void g_multilayer_implicit_update_sediment_cells(int nTasks, t_arrays *arrays){
/*----------------------------*/

    int idx1,idx2,idx;
    int sid1,sid2,sid3;
    int jphi;

    int ncells=arrays->ncells;
    int nInterfaces= arrays->nSediments;

    int nSolutes=arrays->nSolutes;

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
                r[0] = arrays->phi[ncells*nSolutes+idx];
                r[1] = arrays->phi[ncells(nSolutes+1)+idx];
                rbis[0] = r[0];
                rbis[1] = arrays->phi[ncells(nSolutes+1)+idx] - (A[1]/Bi)*r[0];
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
            r[0] = arrays->phi[nSolutes*idx];
            r[1] = arrays->phi[ncells*(nSolutes+1)+idx];
            rbis[0] = r[0];
            rbis[1] = arrays->phi[ncells(nSolutes+1)+idx] - (A[1]/Bi)*r[0];
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
            r[0] = arrays->phi[nSolutes+idx];
            r[1] = arrays->phi[ncells*(nSolutes+1)+idx];
            rbis[0] = r[0];
            rbis[1] = arrays->phi[ncells*(nSolutes+1)+idx] - (A/Bi)*r[0];
            Bbis[1] = B-(A*Ci)/Bi;

            #endif

            #endif


            for(jphi=2;jphi<nInterfaces;jphi++){
                sid1 = (nSolutes+jphi)*ncells+idx; 
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
            arrays->phi[(nInterfaces+nSolutes-1)*ncells+idx] = rbis[nInterfaces-1]/Bbis[nInterfaces-1];
            

            for(jphi=nInterfaces-2; jphi>0; jphi--){
                sid1 = (nSolutes+jphi)*ncells+idx;
                sid2 = (jphi+nSolutes+1)*ncells+idx;
                
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

                sid1 = (nSolutes+jphi)*ncells+idx;
                sid2 = (jphi+nSolutes+1)*ncells+idx; //j+1

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
                sid1 = (nSolutes+jphi)*ncells+idx;
                arrays->phi[sid1] = arrays->hphi[sid1]/hlayer;

                if(std::isnan(arrays->phi[sid1])){
                    printf("cell %d layer %f\n",idx, jphi);
                    //printf("h %lf\n",hlayer);
                    printf("hphi %lf phi %lf hlayer %lf\n ",arrays->hphi[sid1], arrays->phi[sid1], hlayer);
                }
                //printf("hphi %lf phi %lf hlayer %lf\n ",arrays->hphi[sid1], arrays->phi[sid1], hlayer);
            
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



#endif
