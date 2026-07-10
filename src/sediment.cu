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
        for(jphi=0;jphi<arrays->nSediments;jphi++){
            sid = (arrays->nSolutes + jphi)*ncells+idx;
            sid1 = jphi*ncells+idx;
            arrays->Ns[sid1] = 0.0;
            arrays->phiZero[idx] += arrays->phi[sid];

        }
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

        for(jphi=0;jphi<arrays->nSediments;jphi++){
            if(arrays->dsp[jphi] > TOL12){
                TdFsp += arrays->Fsp[jphi]/arrays->dsp[jphi];      
            }else{
                TdFsp += 0.0;
            }
        }
        
        //compact
        //wall index
        idx=arrays->actCells[i];

        h = arrays->h[idx];
        u = arrays->u[idx];
        v = arrays->v[idx];
                
        moduloU = sqrt(u*u + v*v);

        // if(idx == 90526){
        //     printf("moduloU %.12lf u %.12lf v %.12lf\n", moduloU, u, v);
        // }

        phiZero = arrays->phiZero[idx];
        nman = arrays->nman[idx];

        rhoS = arrays->rhoS;

        if(arrays->h[idx] >= arrays->minh){ //wet cells

            
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

                //printf("Fsp %.12lf Css %.12lf dsp %.12lf jphi %d\n", Fsp, Css, dsp, jphi);

                //printf("EquConcFF %.12lf WsFs %.12lf jphi %d\n", arrays->EquConcFF[jphi], WsFs, jphi);
            
                aux1 = rhoBulk*_g_*h;
                aux2 = nman*nman/pow(h,4./3.);
                u_comp = moduloU*u;
                v_comp = moduloU*v;

                Ssu = aux1*aux2*u_comp;
                Ssv = aux1*aux2*v_comp;

                //SsModulus = sqrt(Ssu*Ssu + Ssv*Ssv);
                //SsModulus = _rhow_*aux2*moduloU;

                SsModulus = rhoBulk*_g_*h*(nman*nman*(u*u + v*v))/(pow(h,4./3.));

                //printf("SsModulus %.12lf\n", SsModulus);
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
                        //printf("aux1 %lf\n", aux1);
                        
                        //Csst = 0.01*((rhoS/(rhoS-_rhow_))*((SsModulus*moduloU*moduloU)/WsFs));
                        EquConcFs = 0.01*SsModulus*moduloU/((rhoS-_rhow_)*_g_*h*WsFs);
                        // if(idx == 90526){
                        //     printf("SsModulus %.12lf moduloU  %.12lf\n", SsModulus, moduloU);
                        // }

                        // if(idx == 104477){
                        //     printf("WsFs %.12lf moduloU %.12lf Csst %.12lf EquConcFs %.12lf\n", WsFs,moduloU, Csst, EquConcFs);
                        // }
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
                    }
                }else{
                    EquConcFs = 0.0;
                }

                if(idx == 88511){
                    printf("EquConcFs %.12lf\n", EquConcFs);
                }

                // if(idx == 260){
                //     printf("EquConcFs %.12lf WsFs %d\n", EquConcFs, WsFs);
                // }

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

                if(idx == 88511){
                    printf("Dbj %.12lf Ebj%.12lf\n", Dbj, Ebj);
                    printf("phi %.12lf\n",arrays->phi[sid]);
                }



                arrays->Ns[sid1] = (Ebj - Dbj);
                arrays->Nb[idx] += arrays->Ns[sid1];

                if(idx == 88511){
                  printf("Ns %.12lf Nb %.12lf sid1 %d\n", arrays->Ns[sid1], arrays->Nb[idx], sid1);
                } 


            
            }



        }else{

            for(jphi=0;jphi<arrays->nSediments;jphi++){
                sid1 = jphi*ncells+idx;
                arrays->Ns[sid1] = 0.0;
            }

            arrays->Nb[idx] = 0.0;
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

        //printf("rhoBulk %lf\n", rhoBulk);
        //printf("phiZero %lf\n", phiZero);

        if(arrays->h[idx] > arrays->minh){ //wet cells

            aux1 = 0.0;
            aux2 = 0.0;

            for(jphi=0;jphi<arrays->nSediments;jphi++){
                sid = (arrays->nSolutes + jphi)*ncells+idx;
                aux1 += arrays->hphi[sid]/arrays->h[idx];
                aux2 += (arrays->hphi[sid]/arrays->h[idx])*EtaS;
            }

            //printf("aux1 %lf aux2 %lf\n", aux1, aux2);

            mod_EtaS = 0.0;

            if(aux2>1.0){
                if(aux1>0.0){
                    mod_EtaS = 1./aux1;
                }
            }

            aux2 = 0.0;

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

            cr = 1.0;

            if(abs(aux2)>=TOL12){
                aux3 = aux2;
                if(aux3<0.0){ //deposition limited by h
                    aux3 = fmax(-1.*arrays->h[idx], aux3);
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

            for(jphi=0;jphi<arrays->nSediments;jphi++){
                sid = (arrays->nSolutes + jphi)*ncells+idx;
                sid1 = jphi*ncells+idx;

                if(mod_EtaS>0.0){
                    EtaS_eff = mod_EtaS;
                }else{
                    EtaS_eff = EtaS;
                }

                aux1 = arrays->Ns[sid1]*dt;

                if(idx ==104477){
                    printf("Ns %.12lf dt %.12lf\n", arrays->Ns[sid1], dt);
                }

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

                // if(idx == 104477){
                //     printf("aux1 %.12lf\n", aux1);
                // //printf("porosityCoef %.12lf EtaS %.12lf rhoS %.12lf\n", pd, EtaS, rhoS);
                // }

                //}

                //aux3 = arrays->Nb[idx]*dt;

                // if(idx ==104477){
                //     printf("Nb %.12lf dt %.12lf\n", arrays->Nb[idx], dt);
                // }

                // if(idx == 90526){
                //     printf("z %.12lf EtaS_eff %.12lf\n", arrays->z[idx], EtaS_eff);
                // }

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

                // if(idx == 88511){
                //     printf("u %.12lf v %.12lf\n", arrays->u[idx], arrays->v[idx]);
                //     //printf("porosityCoef %.12lf EtaS %.12lf rhoS %.12lf\n", pd, EtaS, rhoS);
                // }

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

                //     // if(idx == 88511){
                //     //     printf("u %.12lf v %.12lf\n", arrays->u[idx], arrays->v[idx]);
                //     // //printf("porosityCoef %.12lf EtaS %.12lf rhoS %.12lf\n", pd, EtaS, rhoS);
                //     // }

                // }else{
                //     arrays->u[idx] = 0.0;
                //     arrays->v[idx] = 0.0;
                //     arrays->hphi[sid] = 0.0;
                //     arrays->phi[sid] = 0.0;
                // }

            
                // if(idx == 104477){
                //     printf("u %.12lf v %.12lf\n", arrays->u[idx], arrays->v[idx]);
                //     //printf("porosityCoef %.12lf EtaS %.12lf rhoS %.12lf\n", pd, EtaS, rhoS);
                // } 


                // for(jphi=0;jphi<arrays->nSediments;jphi++){
                //     sid = (arrays->nSolutes + jphi)*ncells+idx;
                //     sid1 = jphi*ncells+idx;

                
                if(arrays->h[idx]>TOL12 && arrays->hphi[sid]>0.0){
                    arrays->phi[sid] = arrays->hphi[sid]/arrays->h[idx];
                }else{
                    arrays->phi[sid] = 0.0;
                }
                
            }
            

            if(idx == 88511){
                printf("hphi %.12lf phi %.12lf z %.12lf h %.12lf\n", arrays->hphi[sid], arrays->phi[sid], arrays->z[idx], arrays->h[idx]);
                //printf("porosityCoef %.12lf EtaS %.12lf rhoS %.12lf\n", pd, EtaS, rhoS);
            } 

            // if(std::isnan(arrays->phi[sid])){
            //     printf("cell %d layer %f\n",idx, jphi);
            //     //printf("h %lf\n",hlayer);
            //     printf("hphi %lf phi %lf hlayer %lf\n ",arrays->hphi[sid], arrays->phi[sid], jphi);
            // }

            
            
            
        }
    
    }
}



#endif 