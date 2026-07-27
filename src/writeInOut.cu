#include "writeInOut.cuh"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
EXPORT_DLL void Notify(char *msg, int type, t_message *e){
	FILE *fp;

	//write log file
	fp=fopen(e->logFile,"a+");
	if(type == MSG_ERROR){
		fprintf(fp," [FAIL] %s\n",msg);
	}
	if(type == MSG_WARN){
		fprintf(fp," [WARN] %s\n",msg);
	}
	if(type == MSG_L0){
		fprintf(fp," __ii__ %s\n",msg);
	}
	if(type == MSG_L1){
		fprintf(fp," [_OK_] %s\n",msg);
	}
	if(type == MSG_L2){
		fprintf(fp," [EXEC] %s\n",msg);
	}
	if(type == MSG_L3){
		fprintf(fp," [CUDA] %s\n",msg);
	}
	fclose(fp);



	#if APP_MODE
		if(type == MSG_ERROR){
			sprintf(e->errorProperty,"%s",msg);
			e->error = 1;
		}
		if(type == MSG_WARN){
			sprintf(e->warning[e->nWarning],"%s",msg);
			e->nWarning++;
			e->warning = (char**) realloc(e->warning,sizeof(char)*1024 * (e->nWarning+1));
			e->warning[e->nWarning]=(char*)malloc(sizeof(char)*1024);
		}
		if(type == MSG_L0){
			sprintf(e->MsgL0[e->nMsgL0],"%s",msg);
			e->nMsgL0++;
			e->MsgL0 = (char**) realloc(e->MsgL0,sizeof(char)*1024 * (e->nMsgL0+1));
			e->MsgL0[e->nMsgL0]=(char*)malloc(sizeof(char)*1024);
		}
		if(type == MSG_L1){
			sprintf(e->MsgL1[e->nMsgL1],"%s",msg);
			e->nMsgL1++;
			e->MsgL1 = (char**) realloc(e->MsgL1,sizeof(char)*1024 * (e->nMsgL1+1));
			e->MsgL1[e->nMsgL1]=(char*)malloc(sizeof(char)*1024);
		}
		if(type == MSG_L2){
			sprintf(e->MsgL2[e->nMsgL2],"%s",msg);
			e->nMsgL2++;
			e->MsgL2 = (char**) realloc(e->MsgL2,sizeof(char)*1024 * (e->nMsgL2+1));
			e->MsgL2[e->nMsgL2]=(char*)malloc(sizeof(char)*1024);
		}
		if(type == MSG_L3){
			sprintf(e->MsgL3[e->nMsgL3],"%s",msg);
			e->nMsgL3++;
			e->MsgL3 = (char**) realloc(e->MsgL3,sizeof(char)*1024 * (e->nMsgL3+1));
			e->MsgL3[e->nMsgL3]=(char*)malloc(sizeof(char)*1024);
		}
	#else
		if(type == MSG_ERROR){
			printf("%s %s\n",MSGERROR,msg);
			e->error = 1;
		}
		if(type == MSG_WARN){
			printf("%s %s\n",MSGWARN,msg);
		}
		if(type == MSG_L0){
			printf("%s %s\n",MSGINFO,msg);
		}
		if(type == MSG_L1){
			printf("%s %s\n",MSGOK,msg);
		}
		if(type == MSG_L2){
			printf("%s %s\n",MSGEXC,msg);
		}
		if(type == MSG_L3){
			printf("%s %s\n",MSGGPU,msg);
		}
	#endif

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
EXPORT_DLL int write_vtk_mesh(char *filename, t_mesh *mesh, t_message *msg){
	int i,j;
	t_g_cell *celda;
	t_c_cell *celdaC;
	t_node *nodo;
    t_wall *wall;
	
    int iaux;
    char temp[1024];

    FILE *fp;
	fp = fopen (filename,"w");

    //headers ----------------------------------------------------------------
	fprintf(fp,"# vtk DataFile Version 2.0\n");
	fprintf(fp,"Peka2D v5.0 mesh file\n"); 
	fprintf(fp,"ASCII\n");
	fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");

    //mesh structure ----------------------------------------------------------- 
    if(!write_mesh_structure_in_vtk_ascii(fp, mesh)){
        sprintf(temp,"VTK structure writing failed");
        Notify(temp,MSG_ERROR,msg);	
        return(0);
    }

    //cell data -----------------------------------------------------------------
    fprintf (fp, "CELL_DATA %d\n", mesh->g_cells->n);

	fprintf (fp, "SCALARS cid int\n");
    fprintf (fp, "LOOKUP_TABLE default\n");
    for(i=0;i<mesh->g_cells->n;i++){
        celda=&(mesh->g_cells->cells[i]);
        //celdaC=&(mesh->c_cells->cells[i]);
        fprintf (fp,"%d\n", celda->id);
    }
       
	fprintf (fp, "SCALARS isBoundG int\n");
    fprintf (fp, "LOOKUP_TABLE default\n");
    for(i=0;i<mesh->g_cells->n;i++){
        celda=&(mesh->g_cells->cells[i]);
        fprintf (fp,"%d\n", celda->isBound);
        //celdaC=&(mesh->c_cells->cells[i]);
        //fprintf (fp,"%d\n", celdaC->geom->isBound);
    } 

	fprintf (fp, "SCALARS isBoundC int\n");
    fprintf (fp, "LOOKUP_TABLE default\n");
    for(i=0;i<mesh->g_cells->n;i++){
        //celda=&(mesh->g_cells->cells[i]);
        //fprintf (fp,"%d\n", celda->isBound);
        celdaC=&(mesh->c_cells->cells[i]);
        fprintf (fp,"%d\n", celdaC->geom->isBound);
    }

	fprintf (fp, "SCALARS nneigG int\n");
    fprintf (fp, "LOOKUP_TABLE default\n");
    for(i=0;i<mesh->g_cells->n;i++){
        celda=&(mesh->g_cells->cells[i]);
        fprintf (fp,"%d\n", celda->nneig);
        //celdaC=&(mesh->c_cells->cells[i]);
        //fprintf (fp,"%d\n", celdaC->geom->nneig) 
    }         

	fprintf (fp, "SCALARS nneigC int\n");
    fprintf (fp, "LOOKUP_TABLE default\n");
    for(i=0;i<mesh->g_cells->n;i++){
        //celda=&(mesh->g_cells->cells[i]);
       // fprintf (fp,"%d\n", celda->nneig);
        celdaC=&(mesh->c_cells->cells[i]);
        fprintf (fp,"%d\n", celdaC->geom->nneig);
    }

	fprintf (fp, "SCALARS minNeighWall int\n");
    fprintf (fp, "LOOKUP_TABLE default\n");
    for(i=0;i<mesh->g_cells->n;i++){
        celda=&(mesh->g_cells->cells[i]);
        iaux=0;
        for(j=0;j<celda->nneig;j++){
            wall=celda->neigwall[j];
            iaux=MAX(iaux,wall->idWall);
        }
        fprintf (fp,"%d\n", iaux);
    }

	fprintf (fp, "SCALARS sumTypeWall int\n");
    fprintf (fp, "LOOKUP_TABLE default\n");
    for(i=0;i<mesh->g_cells->n;i++){
        celda=&(mesh->g_cells->cells[i]);
        iaux=0;
        for(j=0;j<celda->nneig;j++){
            wall=celda->neigwall[j];
            iaux+=wall->typeOfBound;
        }
        fprintf (fp,"%d\n", iaux);
    }          

    fclose(fp);        

	return(1);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
EXPORT_DLL int create_computation_files(char *path, t_message *msg){

	char filename[1024];
	FILE *fp;

	//mass balance file
	sprintf(filename,"%smassBalance.out",path);
	fp=fopen(filename,"w");
	fclose(fp);

	//solute mass balance file
	#if SET_SOLUTE
	sprintf(filename,"%ssoluteMassBalance.out",path);
	fp=fopen(filename,"w");
	fclose(fp);
	#endif

	return 1;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
EXPORT_DLL int write_massBalance(char *path, t_arrays *arrays, t_message *msg){

	char temp[1024], filename[1024];
	FILE *fp;

	//mass balance file
	sprintf(filename,"%smassBalance.out",path);
	fp=fopen(filename,"a+");

	if(fp){
		fprintf(fp,"%.3lf %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf\n",
			arrays->t,
			arrays->qTotalIn, arrays->qTotalOut, arrays->mTotalIn,arrays->mTotalOut,
			arrays->massTotalIn, arrays->massTotalOut, arrays->massNew);

		fclose(fp);		
	}else{
        sprintf(temp,"Mass balance file not reacheable");
        Notify(temp,MSG_ERROR,msg);		
	}

	return 1;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
EXPORT_DLL void dump_screen_info(t_arrays *arrays, t_message *msg){

	double aux1,aux2,aux3;

	printf("-----------------------------------------------\n");
	printf("%s Time: %.9f   dt: %.9f   nIter: %d\n",MSGOK,
		arrays->t,
		arrays->dt,
		arrays->nIter);
    printf("%s Active Cells: %d   Active Walls: %d\n",MSGEXC,
		arrays->nActCells,
		arrays->nActWalls);       
	printf("%s Inlet:  Discharge %.6e m3/s   Mass %.6e m3\n",MSGINFO,
		arrays->qTotalIn,
		arrays->mTotalIn);
	printf("%s Outlet: Discharge %.6e m3/s   Mass %.6e m3\n",MSGINFO,
		arrays->qTotalOut,
		arrays->mTotalOut);	
	// printf("%s Mass: In %.6e m3/s   Out %.6e m3\n",MSGINFO,
	// 	arrays->massTotalIn,
	// 	arrays->massTotalOut);			
    printf("%s Mass balance: Inner %.6e m3   Error %.6e\n",MSGINFO,arrays->massNew, arrays->massError);

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
EXPORT_DLL int write_vtk_state(char *filename, t_mesh *mesh, t_arrays *arrays, t_message *msg){
	int i,j;
	t_g_cell *celda;
	t_c_cell *celdaC;
	t_node *nodo;
	
	double auxd1,auxd2,auxd3;
	int auxi;
    char temp[1024];
	int idx;

    FILE *fp;
	fp = fopen (filename,"w");

    //headers ----------------------------------------------------------------
	fprintf(fp,"# vtk DataFile Version 2.0\n");
	fprintf(fp,"Peka2D v5.0 state output\n"); 
	#if BINARY_VTK
		fprintf(fp,"BINARY\n");
	#else
		fprintf(fp,"ASCII\n");
	#endif
	fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");

    //mesh structure ----------------------------------------------------------- 
    #if BINARY_VTK
	    if(!write_mesh_structure_in_vtk_binary(fp, mesh)){
	#else
	    if(!write_mesh_structure_in_vtk_ascii(fp, mesh)){
	#endif
            sprintf(temp,"VTK structure writing failed");
            Notify(temp,MSG_ERROR,msg);	
            return(0);
        }


    //cell data -----------------------------------------------------------------
    fprintf (fp, "CELL_DATA %d\n", mesh->g_cells->n);

	fprintf (fp, "SCALARS h double\n");
	fprintf (fp, "LOOKUP_TABLE default\n");
	for(i=0;i<arrays->ncells;i++){
		write_Dscalar_in_vtk(fp, arrays->h[i]);
	}

	// fprintf (fp, "SCALARS zs double\n");
	// fprintf (fp, "LOOKUP_TABLE default\n");
	// for(i=0;i<arrays->ncells;i++){
	// 	write_Dscalar_in_vtk(fp, arrays->z[i]+arrays->h[i]);
	// }

	fprintf (fp, "VECTORS vel double\n");
	for(i=0;i<arrays->ncells;i++){
		write_Dvector_in_vtk(fp, arrays->u[i], arrays->v[i]);
	}		

	fprintf (fp, "SCALARS modU double\n");
	fprintf (fp, "LOOKUP_TABLE default\n");
	for(i=0;i<arrays->ncells;i++){
		write_Dscalar_in_vtk(fp, arrays->modulou[i]);
	}

	fprintf (fp, "SCALARS z double\n");
	fprintf (fp, "LOOKUP_TABLE default\n");
	for(i=0;i<arrays->ncells;i++){
		write_Dscalar_in_vtk(fp, arrays->z[i]);
	}

	// fprintf (fp, "SCALARS actcell int\n");
	// fprintf (fp, "LOOKUP_TABLE default\n");
	// for(i=0;i<arrays->ncells;i++){
	// 	write_Iscalar_in_vtk(fp, arrays->activeC[i]);
	// } 

	#if SET_SOLUTE || SET_SED

	if(arrays->nSolutes + arrays->nSediments){
		for(j=0;j<arrays->nSolutes + arrays->nSediments;j++){
			if(j<arrays->nSolutes){
				fprintf (fp, "SCALARS solute_phi%d double\n",j);
			}else if(j>=arrays->nSolutes){
				fprintf (fp, "SCALARS sed_phi%d double\n",j-arrays->nSolutes);
			}
			fprintf (fp, "LOOKUP_TABLE default\n");
			for(i=0;i<arrays->ncells;i++){
				idx = j*arrays->ncells+i;
				write_Dscalar_in_vtk(fp, arrays->phi[idx]);
			}
		}	
	}
	#endif	          

    fclose(fp);

	return(1);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
EXPORT_DLL int write_mesh_structure_in_vtk_binary(FILE *fp, t_mesh *mesh){

    int i,j;
	t_g_cell *celda;
    t_node *nodo;

    double auxd;
	int auxi;

    fprintf(fp,"POINTS %d double\n",mesh->nodes->n);
	for (i=0;i<mesh->nodes->n;i++){
		nodo=&(mesh->nodes->nodes[i]);
        auxd=nodo->x;
        SWAP_DOUBLE(auxd);
        fwrite(&auxd,sizeof(double),1,fp);
        auxd=nodo->y;
        SWAP_DOUBLE(auxd);
        fwrite(&auxd,sizeof(double),1,fp);
        auxd=0.0;
        SWAP_DOUBLE(auxd);
        fwrite(&auxd,sizeof(double),1,fp);
	}

    fprintf(fp,"CELLS %d %d\n", mesh->g_cells->n, mesh->g_cells->n * (mesh->NCwall+1));
	for(i=0;i<mesh->g_cells->n;i++){
		celda=&(mesh->g_cells->cells[i]);
        auxi=mesh->NCwall;
        SWAP_INT(auxi);
        fwrite(&auxi,sizeof(int),1,fp);
        auxi=celda->nodes[0]->id;
        SWAP_INT(auxi);
        fwrite(&auxi,sizeof(int),1,fp);
        auxi=celda->nodes[1]->id;
        SWAP_INT(auxi);
        fwrite(&auxi,sizeof(int),1,fp);
        auxi=celda->nodes[2]->id;
        SWAP_INT(auxi);
        fwrite(&auxi,sizeof(int),1,fp);
        if(mesh->NCwall==4){
            auxi=celda->nodes[3]->id;
            SWAP_INT(auxi);
            fwrite(&auxi,sizeof(int),1,fp);
        }
	}

	fprintf(fp,"CELL_TYPES %d\n",mesh->g_cells->n);
	for(i=0;i<mesh->g_cells->n;i++){
        if(mesh->NCwall==5){
            auxi=5;
        }else{
            auxi=9;
        }
        SWAP_INT(auxi);
        fwrite(&auxi,sizeof(int),1,fp);
	}

    return(1);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
EXPORT_DLL int write_mesh_structure_in_vtk_ascii(FILE *fp, t_mesh *mesh){

    int i,j;
	t_g_cell *celda;
    t_node *nodo;

    double auxd;
	int auxi;

    fprintf(fp,"POINTS %d double\n",mesh->nodes->n);
	for (i=0;i<mesh->nodes->n;i++){
		nodo=&(mesh->nodes->nodes[i]);
		fprintf(fp,"%.9f %.9f 0.0\n", nodo->x,nodo->y);
	}

    fprintf(fp,"CELLS %d %d\n", mesh->g_cells->n, mesh->g_cells->n * (mesh->NCwall+1));
	for(i=0;i<mesh->g_cells->n;i++){
		celda=&(mesh->g_cells->cells[i]);
        fprintf(fp,"%d ",mesh->NCwall);
        fprintf(fp,"%d %d %d ",celda->nodes[0]->id,celda->nodes[1]->id,celda->nodes[2]->id);
        if(mesh->NCwall==4){
            fprintf(fp,"%d ",celda->nodes[3]->id);
        }
        fprintf(fp," \n");           
	}

	fprintf(fp,"CELL_TYPES %d\n",mesh->g_cells->n);
	for(i=0;i<mesh->g_cells->n;i++){
        if(mesh->NCwall==3){
            fprintf(fp,"5\n");
        }else{
            fprintf(fp,"9\n");
        }
	}

    return(1);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
EXPORT_DLL void write_Dscalar_in_vtk(FILE *fp, double val){
    
    double auxd;
   
    auxd=val;
    #if BINARY_VTK
        SWAP_DOUBLE(auxd);
        fwrite(&auxd,sizeof(double),1,fp);
    #else
        fprintf (fp,"%.9f\n", auxd);
    #endif
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
EXPORT_DLL void write_Iscalar_in_vtk(FILE *fp, int val){
    
    int auxi;
   
    auxi=val;
    #if BINARY_VTK
        SWAP_DOUBLE(auxi);
        fwrite(&auxi,sizeof(int),1,fp);
    #else
        fprintf (fp,"%d\n", auxi);
    #endif
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
EXPORT_DLL void write_Dvector_in_vtk(FILE *fp, double val1, double val2){
   
    double auxd1,auxd2,auxd3;
   
    auxd1=val1;
    auxd2=val2;
    auxd3=0.0;
    #if BINARY_VTK
        SWAP_DOUBLE(auxd1);
        fwrite(&auxd1,sizeof(double),1,fp);
        SWAP_DOUBLE(auxd2);
        fwrite(&auxd2,sizeof(double),1,fp);
        SWAP_DOUBLE(auxd3);
        fwrite(&auxd3,sizeof(double),1,fp);
    #else
        fprintf (fp,"%.9f %.9f 0.0\n",auxd1,auxd2);
    #endif
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
EXPORT_DLL int write_hotstart_file(char *filename, t_arrays *arrays, t_message *msg){

	int i,j;
	int ncells=arrays->ncells;
	int nSolutes=arrays->nSolutes;
	int nSediments= arrays->nSediments;

	char temp[1024];
	FILE *fp;

	fp=fopen(filename,"w");

	if(fp){
		//header
		fprintf(fp,"4 %d 0 0\n",(arrays->nSolutes+arrays->nSediments));

		//cells data
		for(i=0;i<ncells;i++){
			fprintf(fp,"%lf %lf %lf %lf",
				arrays->z[i],
				arrays->h[i],
				arrays->u[i],
				arrays->v[i]);
			#if SET_SOLUTE
			for(j=0;j<(nSolutes+nSediments);j++){
				fprintf(fp," %lf",arrays->phi[j*ncells+i]);
			}
			#endif
			fprintf(fp," \n");
		}
		fclose(fp);	
			
	}else{
		sprintf(temp,"Hotstart file not reacheable");
        Notify(temp,MSG_ERROR,msg);	
		return 0;	
	}

	return 1;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
EXPORT_DLL int dump_probes_gpu(t_mesh *mesh, t_arrays *arrays, char *path, t_message *msg, double t){
	FILE *fp;
	int i,j;
	int idc;
	t_c_cell *c;
	char temp[1024], filename[1024];
	double h,hz,u,v,moduloU;
	double phi;
	
	for(i=0;i<mesh->nProbes;i++){

		sprintf(filename,"%sprobe%d.out",path,i+1);

		fp=fopen(filename,"a+");


		//printf("idc %d/n", mesh->probe->probe[i].idc);
		if(fp){
			idc = mesh->probe->probe[i].idc;
			if(mesh->probe->probe[i].idc<0){// No habia sonda

				fprintf(fp,"%lf 0.0 0.0 0.0 0.0",t);
					if(mesh->nSolutes>0){
						for(j=0;j<mesh->nSolutes;j++){
							fprintf(fp," 0.0");
						}
					}

					if(mesh->nSediments>0){
						for(j=0;j<mesh->nSediments;j++){
							fprintf(fp,"0.0");
						}
					}
				fprintf(fp,"\n");

			}else{
				//c=mesh->probe->probe[i].cell;
				

				fprintf(fp,"%lf %.6lf %.6lf %.6lf %.6lf %.6lf",
					t,arrays->h[idc] + arrays->z[idc],
					arrays->h[idc],
					arrays->u[idc],
					arrays->v[idc],
					arrays->modulou[idc]);

				if(mesh->nSolutes>0 || mesh->nSediments>0){
					for(j=0;j<mesh->nSolutes + mesh->nSediments;j++){
						fprintf(fp," %.6lf",arrays->phi[j*arrays->ncells + idc]);
					}
				}
			
				fprintf(fp,"\n");
				fclose(fp);
			}
		
		}else{
			sprintf(temp,"obs file not reacheable");
        	Notify(temp,MSG_ERROR,msg);	
		}
		
	}
	return 1;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
EXPORT_DLL int dump_sections_gpu(t_mesh *mesh, t_arrays *arrays, double t,char* dir){
	FILE *f;
	int i,j,k;
	int jphi ;
	int idc, idc2;
	t_c_cell *c;
	t_section *sec;
	double sumQx,sumQy,sumH,sumZ,sumQ,sumVx,sumVy,sumphi, sumphisol, sumphised;
	double npoints,deltaX,npointsEfectivos,npointsMojados;
	double aux1;
	double phi;
	char nombre[1024];
	int nSolutes = arrays->nSolutes;
	int nSediments = arrays->nSediments;
	int ncells = arrays->ncells;

	sec=mesh->sec->sec;

	if(mesh->nSections>0){
		for(i=0;i<mesh->nSections;i++){
			sumQ=0.0;
			sec->zAverage=0.0;
			sec->hAverage=0.0;
			sec->phiAveragedSol=0.0;
			sec->phiAveragedSed=0.0;
		
			npointsEfectivos=0.0;
			npointsMojados=0.0;

			for(j=0;j<sec->npoints-1;j++){

				sumQx=0.;
				sumQy=0.;
				sumH=0.0;
				sumZ=0.0;
				deltaX=0.0;
				npoints=0.0;

				idc = sec->probeSec->probe[j].idc;

				if(sec->probeSec->probe[j].idc>0){
					//c=sec->probeSec->probe[j].cell;
					sumZ+=arrays->z[idc]*0.5;

					if(arrays->h[idc]>0.0){
						sumQx+=arrays->h[idc]*arrays->u[idc];
						sumQy+=arrays->h[idc]*arrays->v[idc];
						sumH+=arrays->h[idc];

						aux1 = 0.0;
						for(jphi=0;jphi<nSolutes;jphi++){
							aux1 += arrays->phi[jphi*ncells +idc]/nSolutes;
						}
						sumphisol+=aux1*0.5;

						aux1 = 0.0;
						for(jphi=nSolutes;jphi<nSolutes+nSediments;jphi++){
							aux1 += arrays->phi[jphi*ncells +idc]/nSediments;
						}
						sumphised+=aux1*0.5;

						deltaX+=(sec->dist/(sec->npoints-1))*0.5;
						npoints=npoints+0.5;
						npointsMojados+=0.5;
					}
					npointsEfectivos+=0.5;

				}
				idc2 = sec->probeSec->probe[j+1].idc;

				if(sec->probeSec->probe[j+1].idc>0){
					//c=sec->probeSec->probe[j+1].cell;
					sumZ+=arrays->z[idc2]*0.5;

					if(arrays->h[idc2]>0.0){
						sumQx+=arrays->h[idc]*arrays->u[idc2];
						sumQy+=arrays->h[idc]*arrays->v[idc2];
						sumH+=arrays->h[idc2];

						aux1 = 0.0;
						for(jphi=0;jphi<nSolutes;jphi++){
							aux1 += arrays->phi[jphi*ncells +idc2]/nSolutes;
						}
						sumphisol+=aux1*0.5;

						aux1 = 0.0;
						for(jphi=nSolutes;jphi<nSolutes+nSediments;jphi++){
							aux1 += arrays->phi[jphi*ncells +idc2]/nSediments;
						}
						sumphised+=aux1*0.5;

						deltaX+=(sec->dist/(sec->npoints-1))*0.5;
						npoints=npoints+0.5;
						npointsMojados+=0.5;
					}
					npointsEfectivos+=0.5;

				}

				//printf("sumQx %lf sumQy %lf\n", sumQx,sumQy);

				sec->zAverage+=sumZ;

				if(npoints>0.0){
					sec->hAverage+=0.5*sumH/npoints;
					sumQ+=0.5*(sumQx*sec->normal[_X_]+ sumQy*sec->normal[_Y_])*(deltaX/(npoints));
				}
				//printf("sumQ %lf\n", sumQ);

			}


			if(npointsEfectivos>0.0){
				sec->zAverage=sec->zAverage/npointsEfectivos;
			}else{
				sec->zAverage=0.0;
			}

			if(npointsMojados>0.0){
				sec->hAverage=sec->hAverage/npointsMojados;	
			}else{
				sec->hAverage=0.0;
			}

			//printf("sumQ %lf\n", sumQ);

			sec->discharge=sumQ;
			sec->phiAveragedSol=sumphisol/npoints;
			sec->phiAveragedSed=sumphised/npoints;


			//printf("discharge %lf\n", sec->discharge);
			
			sec++;

			//printf("discharge %lf\n", sec->discharge);

		}


		sprintf(nombre,"%sdischarge_section.out",dir);
		f=fopen(nombre,"a+");
		fprintf(f,"%lf",t);

		for(i=0;i<mesh->nSections;i++){
			fprintf(f," %lf",mesh->sec->sec[i].discharge);
		}

		fprintf(f,"\n");
		fclose(f);


		for(i=0;i<mesh->nSections;i++){
			sprintf(nombre,"%saverage_section%d.out",dir,i+1);
			f=fopen(nombre,"a+");
			fprintf(f,"%lf",t);

			fprintf(f," %lf %lf %lf %lf",
					mesh->sec->sec[i].zAverage,
					mesh->sec->sec[i].hAverage,
					mesh->sec->sec[i].zAverage + mesh->sec->sec[i].hAverage, 
					mesh->sec->sec[i].discharge);

			#if SET_SOLUTE 
			fprintf(f," %lf",
					mesh->sec->sec[i].phiAveragedSol);
			#endif

			#if SET_SED
			fprintf(f," %lf",
					mesh->sec->sec[i].phiAveragedSed);
			#endif

			fprintf(f,"\n");
			fclose(f);
		}


		sec=mesh->sec->sec;

		for(i=0;i<mesh->nSections;i++){

			sprintf(nombre,"%ssection%d.out",dir,i+1);
			f=fopen(nombre,"a+");
			
			fprintf(f,"TIME: %lf \t DISCHARGE: %lf\n",t,sec->discharge);


			for(j=0;j<sec->npoints;j++){
				idc = sec->probeSec->probe[j].idc;

				if(sec->probeSec->probe[j].idc>0){
					c=sec->probeSec->probe[j].cell;
					fprintf(f,"%d %.3f %lf %lf %lf %lf %lf",c->id,sec->deltaX*j,
						arrays->z[idc],
						arrays->h[idc],
						arrays->h[idc]+arrays->z[idc],
						arrays->u[idc],
						arrays->v[idc]);

					#if SET_SOLUTE 
					for(jphi=0;jphi<nSolutes;jphi++){
						fprintf(f," %lf",arrays->phi[jphi*ncells+idc]);
					}
					#endif

					#if SET_SED
					for(jphi=nSolutes;jphi<nSolutes+nSediments;jphi++){
						fprintf(f," %lf", arrays->phi[jphi*ncells+idc]);
					}
					#endif

				fprintf(f,"\n");
				}else{
					fprintf(f,"%d %.3f %lf %lf %lf %lf %lf",-1,sec->deltaX*j,-9999.0,-9999.0,-9999.0,-9999.0,-9999.0);

					#if SET_SOLUTE
					for(jphi=0;jphi<nSolutes;jphi++){
						fprintf(f,"%lf",-9999.0);
					}
					#endif 
					#if SET_SED
					for(jphi=nSolutes;jphi<nSolutes+nSediments;jphi++){
						fprintf(f,"%lf",-9999.0);
					}
					#endif

				fprintf(f,"\n");
				}

			}// end for(j=0;j<sec->npoints;j++){

			fclose(f);
			sec++;
		}

		#if SET_SOLUTE || SET_SED

		sec=mesh->sec->sec;

		for(i=0;i<mesh->nSections;i++){
			sprintf(nombre,"%saveragePHI_section%d.out",dir,i+1);
			f=fopen(nombre,"a+");
			fprintf(f,"%lf",t);
		
			if(mesh->nSections>0){

				for(jphi=0;jphi<nSolutes + nSediments;jphi++){
					
					sec->phiAveraged = 0.0;

					npointsMojados=0.0;

					for(j=0;j<sec->npoints-1;j++){

						sumphi = 0.0;
						npoints = 0.0;

						idc = sec->probeSec->probe[j].idc;

						if(sec->probeSec->probe[j].idc>0){
							if(arrays->h[idc]>0.0){
								sumphi+= arrays->phi[jphi *ncells + idc]*0.5;
								npoints=npoints+0.5;
								npointsMojados+=0.5;
							}
						}

						//printf("sum phi %lf\n", sumphi);


						idc2 = sec->probeSec->probe[j+1].idc;

						if(sec->probeSec->probe[j+1].idc>0){
							if(arrays->h[idc2]>0.0){
								sumphi+= arrays->phi[jphi *ncells + idc2]*0.5;
								npoints=npoints+0.5;
								npointsMojados+=0.5;
							}
						}

						//printf("sum phi %lf\n", sumphi);

					}
					//printf("sum phi %lf\n", sumphi);
					if(npoints>0.0){
						sec->phiAveraged += sumphi/npoints;
					}
					//printf("phiAve %lf\n", sec->phiAveraged);

					if(npointsMojados>0.0){
						sec->phiAveraged = sec->phiAveraged/npointsMojados;	
					}else{
						sec->phiAveraged = 0.0;
					}
				
					//fprintf(f," %lf", mesh->sec->sec[i].phiAveraged);
						
				}

				//fprintf(f," %lf", mesh->sec->sec[i].phiAveraged);
			fprintf(f,"\n");
			fclose(f);
			}

		sec++;
		}	
		#endif

	}
	return 1;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
EXPORT_DLL int write_timers(char *path, t_timers timers, t_message *msg){

	char temp[1024], filename[1024];
	FILE *fp;

	//mass balance file
	sprintf(filename,"%scomputationTime.out",path);
	fp=fopen(filename,"w");

	if(fp){
		fprintf(fp,"total              %12.6lf %12.3lf\n", timers.total, (timers.total/timers.total)*100.);
		fprintf(fp,"__loading          %12.6lf %12.3lf\n", timers.loading, (timers.loading/timers.total)*100.);
		fprintf(fp,"__init             %12.6lf %12.3lf\n", timers.init, (timers.init/timers.total)*100.);
		fprintf(fp,"__initGPU          %12.6lf %12.3lf\n", timers.initGPU, (timers.initGPU/timers.total)*100.);
		fprintf(fp,"__computeSim       %12.6lf %12.3lf\n", timers.computeSim, (timers.computeSim/timers.total)*100.);
		fprintf(fp,"____wallCalculus   %12.6lf %12.3lf\n", timers.wallCalculus, (timers.wallCalculus/timers.computeSim)*100.);
		fprintf(fp,"__multilayerCalc   %12.6lf %12.3lf\n", timers.multilayer_calculus, (timers.multilayer_calculus/timers.computeSim)*100.);
		fprintf(fp,"____cellUpdating   %12.6lf %12.3lf\n", timers.cellUpdating, (timers.cellUpdating/timers.computeSim)*100.);
		fprintf(fp,"____boundConditon  %12.6lf %12.3lf\n", timers.boundConditon, (timers.boundConditon/timers.computeSim)*100.);
		fprintf(fp,"____wetDryFix      %12.6lf %12.3lf\n", timers.wetDryFix, (timers.wetDryFix/timers.computeSim)*100.);
		fprintf(fp,"____diffusion      %12.6lf %12.3lf\n", timers.diffusion, (timers.diffusion/timers.computeSim)*100.);
		fprintf(fp,"____memoryTransfer %12.6lf %12.3lf\n", timers.memoryTransfer, (timers.memoryTransfer/timers.computeSim)*100.);
		fprintf(fp,"____writeOut       %12.6lf %12.3lf\n", timers.writeOut, (timers.writeOut/timers.computeSim)*100.);
		fprintf(fp,"__closeSim         %12.6lf %12.3lf\n", timers.closeSim, (timers.closeSim/timers.total)*100.);
		fclose(fp);		
	}else{
        sprintf(temp,"Computation time file not reacheable");
        Notify(temp,MSG_ERROR,msg);		
	}

	return 1;

}



