#include "mesh.h"

/*******************************************************/
// Creates a new t_c_cell element
// Make first
// lcells=(l_g_cells*)malloc(sizeof(l_g_cells));
// lcells->cells=(t_c_cell*)malloc(n*sizeof(t_c_cell));
int new_lg_cells(l_g_cells *lcells,int n){
	lcells->size=n;
	lcells->n=0;
    lcells->cells = (t_g_cell*) malloc( n * sizeof(t_g_cell) );
	return(1);
}

/*******************************************************/
// Creates a new t_c_cell element
// Make first
// lcells=(l_g_cells*)malloc(sizeof(l_g_cells));
// lcells->cells=(t_c_cell*)malloc(n*sizeof(t_c_cell));
int new_lc_cells(l_c_cells *lcells,int n){
	lcells->size=n;
	lcells->n=0;
    lcells->cells = (t_c_cell*) malloc( n * sizeof(t_c_cell));
	return(1);
}

/*******************************************************/
// Creates a new t_c_cell element
// Make first cell=(t_c_cell*)malloc(sizeof(t_c_cell));
int new_g_cell(t_g_cell *cell, int id){
	cell->id=id;
	cell->center[0]=0.0;
	cell->center[1]=0.0;
	cell->center[2]=0.0;
	cell->nneig=0;
	cell->nbound=0;
	cell->isBound=0;
    cell->wall = (t_edge*) malloc( 4 * sizeof(t_edge));
	//printf("%s Created new geom-cell [%d]\n",MSGOK,id);
	return(1);
}

/*******************************************************/
// Creates a new t_c_cell element
// Make first cell=(t_c_cell*)malloc(sizeof(t_c_cell));
int new_c_cell(t_c_cell *cell, int id){
	int i,j;
	cell->id=id;
	cell->h=0.0;
	cell->hu=0.0;
	cell->hv=0.0;
	cell->z=0.0;
	//printf("%s Created new compute-cell [%d]\n",MSGOK,id);
	return(1);
}

/*******************************************************/
int add_g_cell(t_g_cell *cell, l_g_cells *lcells){
	if(lcells->n<lcells->size){
		lcells->cells[(lcells->n)]=*cell;
		lcells->n++;
		return(1);
	}else{
		//printf("%s You've outnumbered geom-cells list\n",MSGERROR);
		return(0);
	}
}

/*******************************************************/
int add_c_cell(t_c_cell *cell, l_c_cells *lcells){
	if(lcells->n<lcells->size){
		lcells->cells[(lcells->n)]=*cell;
		lcells->n++;
		return(1);
	}else{
		//printf("%s You've outnumbered computed-cells list\n",MSGERROR);
		return(0);
	}
}

/*******************************************************/
// Creates a new t_c_cell element
// Make first
// lnodes=(l_nodes*)malloc(sizeof(l_nodes));
// lnodes->nodes=(t_node*)malloc(n*sizeof(t_node));
int new_lnodes(l_nodes *lnodes,int n){
	lnodes->size=n;
	lnodes->n=0;
	lnodes->nodes=(t_node*)malloc(n*sizeof(t_node));
	return(1);
}

/*******************************************************/
// Creates a new t_node element
// Make first node=(t_node*)malloc(sizeof(t_node));
int new_node(t_node *node, int id){
	node->id=id;
	node->x=0.0;
	node->y=0.0;
	node->z=0.0;
	//printf("%s Created new node [%d]\n",MSGOK,id);
	return(1);
}

/*******************************************************/
int add_node(t_node *node, l_nodes *lnodes){
	if(lnodes->n<lnodes->size){
		lnodes->nodes[(lnodes->n)]=*node;
		lnodes->n++;
		return(1);
	}else{
		//printf("%s You've outnumbered nodes list\n",MSGERROR);
		return(0);
	}
}

/*******************************************************/
// Creates a new l_wall element
// Make first
//lwall=(l_wall*)malloc(sizeof(l_wall));
//lwall->wall=(t_wall*)malloc( n *sizeof(t_wall));
int new_lwalls(l_wall *lwalls, int n){
	lwalls->size=n;
	lwalls->n=0;
	lwalls->wall=(t_wall*)malloc(n*sizeof(t_wall));
	return(1);
}

/*******************************************************/
// It Calculates the normal vectors of each wall on the cell
int calc_norm(l_g_cells *lcells, l_c_cells *lccells, int nwall, t_message *msg){
	int i,j;
	int k1,k2;
	t_g_cell *myCell;
	t_c_cell *mycCell;
	double kx1,kx2,ky1,ky2,kx3,ky3;
	double mX,mY;
	double vX,vY;
	double mwX,mwY;
	double cX,cY;
	double vcX,vcY;
	double vnX,vnY;
	double vax1,vax2,vay1,vay2,vnx1,vny1;
	double myArea;

    char temp[1024];

	for(i=0;i<(lcells->n);i++){
		myCell=&(lcells->cells[i]);
		mycCell=&(lccells->cells[i]);

		//cell center
		if(nwall==3){
			cX=myCell->nodes[0]->x+myCell->nodes[1]->x+myCell->nodes[2]->x;
			cX=cX/3.0;
			cY=myCell->nodes[0]->y+myCell->nodes[1]->y+myCell->nodes[2]->y;
			cY=cY/3.0;
		}else{
			if(myCell->nodes[0]->x==myCell->nodes[1]->x){
				cX=myCell->nodes[0]->x+myCell->nodes[2]->x;
				cX=cX*0.5;
			}else{
				cX=myCell->nodes[0]->x+myCell->nodes[1]->x;
				cX=cX/2.0;
			}
			if(myCell->nodes[0]->y==myCell->nodes[1]->y){
				cY=myCell->nodes[0]->y+myCell->nodes[2]->y;
				cY=cY*0.5;
			}else{
				cY=myCell->nodes[0]->y+myCell->nodes[1]->y;
				cY=cY/2.0;
			}
		}
		myCell->center[0]=cX;
		myCell->center[1]=cY;
        //printf("Center[%d]: x %.9f y %.9f\n",i,myCell->center[0],myCell->center[1]);

		//cell area
		if(nwall==3){
			kx1=myCell->nodes[0]->x;
			kx2=myCell->nodes[1]->x;
			kx3=myCell->nodes[2]->x;
			ky1=myCell->nodes[0]->y;
			ky2=myCell->nodes[1]->y;
			ky3=myCell->nodes[2]->y;
			vax1=(kx1-kx2);
			vax2=(kx3-kx2);
			vay1=(ky1-ky2);
			vay2=(ky3-ky2);
			vnx1=-vay1;
			vny1=vax1;
			myArea=(fabs((vnx1*vax2)+(vny1*vay2)))/2.0;
		}else{
			kx1=myCell->nodes[0]->x;
			kx2=myCell->nodes[1]->x;
			if(kx1==kx2){
				kx1=myCell->nodes[0]->x;
				kx2=myCell->nodes[2]->x;
			}
			myArea=fabs(kx2-kx1)*fabs(kx2-kx1);	            
		}
        myCell->area=myArea;
		//printf("Area[%d]:%.9f\n",i,myCell->area);
        
        //wall length and normal vector
		if(nwall==3){
            for(j=0;j<nwall;j++){
                if(j==0){
                    k1=1;
                    k2=2;
                }else{
                    if(j==1){
                        k1=2;
                        k2=0;
                    }else{
                        k1=0;
                        k2=1;
                    }
                }
                //printf("k1,k2:%d %d\n",k1,k2);
                kx1=myCell->nodes[k1]->x;
                kx2=myCell->nodes[k2]->x;
                ky1=myCell->nodes[k1]->y;
                ky2=myCell->nodes[k2]->y;
                mwX=(kx1+kx2)/2.0;
                mwY=(ky1+ky2)/2.0;

                // The vector between the center of the wall and the center of the cell
                vcX=cX-mwX;
                vcY=cY-mwY;
                mX=fabs(kx2-kx1);
                mY=fabs(ky2-ky1);
                myCell->wall[j].length=sqrt((mX*mX)+(mY*mY));
                
                // Vectores normales
                vX=(kx1-kx2);
                vY=(ky1-ky2);
                vnX=-vY;
                vnY=vX;
                // We must select the normal vector pointing outside the wall
                if((vnX*vcX)+(vnY*vcY)>0){//it points inside!
                    vnX=-vnX;
                    vnY=-vnY;
                }
                myCell->wall[j].normal[0]=vnX/(sqrt((vnX*vnX)+(vnY*vnY)));
                myCell->wall[j].normal[1]=vnY/(sqrt((vnX*vnX)+(vnY*vnY)));
                //printf("nx: %.9f ny: %.9f\n",myCell->wall[j].normal[0],myCell->wall[j].normal[1]);

            }
		}else{ //nwall==4
            for(j=0;j<nwall;j++){
                if(j==0){
                    k1=2;
                    k2=3;
                }else{
                    if(j==1){
                        k1=3;
                        k2=0;
                    }else{
                        if(j==2){
                        k1=0;
                        k2=1;
                        }else{
                            k1=1;
                            k2=2;
                        }
                    }
                }
                //printf("k1,k2:%d %d %.9f %.9f\n",k1,k2,cX,cY);
                kx1=myCell->nodes[k1]->x;
                kx2=myCell->nodes[k2]->x;
                ky1=myCell->nodes[k1]->y;
                ky2=myCell->nodes[k2]->y;
                //printf("%d\n",nwall);
                mwX=(kx1+kx2)/2.0;
                mwY=(ky1+ky2)/2.0;

                // The vector between the center of the wall and the center of the cell
                vcX=cX-mwX;
                vcY=cY-mwY;
                mX=fabs(kx2-kx1);
                mY=fabs(ky2-ky1);
                myCell->wall[j].length=MAX(fabs(kx1-kx2),fabs(ky1-ky2));

                // Vectores normales
                vX=(kx1-kx2);
                vY=(ky1-ky2);
                vnX=-vY;
                vnY=vX;

                // We must select the normal vector pointing outside the wall
                if((vnX*vcX)+(vnY*vcY)>0){//it points inside!
                    vnX=-vnX;
                    vnY=-vnY;
                }
                myCell->wall[j].normal[0]=vnX/(sqrt((vnX*vnX)+(vnY*vnY)));
                myCell->wall[j].normal[1]=vnY/(sqrt((vnX*vnX)+(vnY*vnY)));
                //printf("nx: %.9f ny: %.9f\n",myCell->wall[j].normal[0],myCell->wall[j].normal[1]);
            }   
		}
	}

    sprintf(temp,"Mesh geometry calculation completed");
    Notify(temp,MSG_L1,msg);

	return(1);
}



/*******************************************************/
// It constucts the lpocells structure from a mesh
int cons_lpocells_qsort(t_mesh *mesh, t_message *msg){
	int *celdasEnPunto;
	int found;
	FILE *fp;
	int *nCeldasEnPunto;
	t_pair *paredes,*parceldas;
	int *celdas;
	int i,j,k;
	int id;
	int p1,p2;
	int np;
	t_wall *tpc;
   	t_g_cell *c1,*c2;
	t_c_cell *cc1,*cc2;
	double x1,x2,y1,y2;
	double cx1,cx2,cy1,cy2;
	double maxlwall, aux1;
	int nparceldas;
    char temp[1024];

	paredes=(t_pair*)malloc(sizeof(t_pair)*mesh->ncells*mesh->NCwall);
	if(paredes==NULL){
        sprintf(temp,"Geometry allocation failes due to memory space");
		Notify(temp,MSG_ERROR,msg);        
        return(0);  
	}

	parceldas=(t_pair*)malloc(sizeof(t_pair)*mesh->ncells*mesh->NCwall);
	if(parceldas==NULL){
        sprintf(temp,"Geometry allocation failes due to memory space");
		Notify(temp,MSG_ERROR,msg);
        return(0);
	}

	celdas=(int*)malloc(sizeof(int)*mesh->ncells*mesh->NCwall);
	if(celdas==NULL){
        sprintf(temp,"Geometry allocation failes due to memory space");
		Notify(temp,MSG_ERROR,msg);        
        return(0);        
	}

	nCeldasEnPunto=(int*)malloc(sizeof(int)*mesh->nodes->n);
	if(nCeldasEnPunto==NULL){
        sprintf(temp,"Geometry allocation failes due to memory space");
		Notify(temp,MSG_ERROR,msg);        
        return(0);
	}

	// Creamos CeldasEnPunto y Celdas
	for(i=0;i<mesh->nodes->n;i++){
		nCeldasEnPunto[i]=0;
	}

    //cell loop
    k=0;
	for(i=0;i<mesh->ncells;i++){
		c1=&mesh->g_cells->cells[i];
		for(j=0;j<mesh->NCwall;j++){
			celdas[k]=c1->nodes[j]->id;

			p1=c1->nodes[j]->id;
			p2=c1->nodes[(j+1)%mesh->NCwall]->id;
			id=c1->nodes[j]->id;

			np=nCeldasEnPunto[id];
			nCeldasEnPunto[id]++;

			paredes[k].included=0;
			if(p1<p2){
				paredes[k].id1=p1;
				paredes[k].id2=p2;
				if(mesh->NCwall==3){
					paredes[k].idW=(j+mesh->NCwall-1)%mesh->NCwall;
				}else{ //mesh->NCwall==4
					paredes[k].idW=(j+mesh->NCwall-2)%mesh->NCwall;
				}
				paredes[k].idC=i;
			}else{
				paredes[k].id1=p2;
				paredes[k].id2=p1;
				if(mesh->NCwall==3){
					paredes[k].idW=(j+mesh->NCwall-1)%mesh->NCwall;
				}else{ //mesh->NCwall==4
					paredes[k].idW=(j+mesh->NCwall-2)%mesh->NCwall;
				}
				paredes[k].idC=i;
			}
			k++;
		}
	}

	mesh->w_calc->n=0;
	qsort(paredes,mesh->ncells*mesh->NCwall,sizeof(t_pair),cmpp);
	nparceldas=0;
	for(i=0;i<(mesh->ncells*mesh->NCwall-1);i++){
		if(	(paredes[i].id1==paredes[i+1].id1) &&
            (paredes[i].id2==paredes[i+1].id2) ){

			paredes[i].included=1;
			paredes[i+1].included=1;	

			parceldas[nparceldas].id1=paredes[i].idC;
			parceldas[nparceldas].id2=paredes[i+1].idC;
			parceldas[nparceldas].idW=paredes[i].idW;
			parceldas[nparceldas].idW2=paredes[i+1].idW;
			parceldas[nparceldas].n1=paredes[i].id1;
			parceldas[nparceldas].n2=paredes[i].id2;

			nparceldas++;
		}
	}
	qsort(parceldas,nparceldas,sizeof(t_pair),cmpp);  
	for(i=0;i<nparceldas;i++){
			tpc=&(mesh->w_calc->wall[mesh->w_calc->n]);

            tpc->idWall=mesh->w_calc->n;
			tpc->typeOfBound=-1;            

			tpc->id=parceldas[i].id1;
			tpc->id2=parceldas[i].id2;

			//nc1=cel
			c1=&mesh->g_cells->cells[parceldas[i].id1];
			c2=&mesh->g_cells->cells[parceldas[i].id2];
			cc1=&mesh->c_cells->cells[parceldas[i].id1];
			cc2=&mesh->c_cells->cells[parceldas[i].id2];
			tpc->length=c1->wall[parceldas[i].idW].length;
			tpc->lwall=c1->wall[parceldas[i].idW].length;
			tpc->normal[_X_]=c1->wall[parceldas[i].idW].normal[0];
			tpc->normal[_Y_]=c1->wall[parceldas[i].idW].normal[1];

			// Distancias normales
			cx1=c1->center[_X_];
			cx2=c2->center[_X_];
			cy1=c1->center[_Y_];
			cy2=c2->center[_Y_];
			tpc->distNormal=fabs((cx2-cx1)*tpc->normal[_X_]+(cy2-cy1)*tpc->normal[_Y_]);
			tpc->distCenter=sqrt((cx2-cx1)*(cx2-cx1)+(cy2-cy1)*(cy2-cy1));
			tpc->ncx = (cx2-cx1) / tpc->distCenter;
			tpc->ncy = (cy2-cy1) / tpc->distCenter;            

			tpc->deltaX=MIN(c1->area/c1->wall[0].length,c1->area/c1->wall[1].length);
			tpc->deltaX=MIN(tpc->deltaX,c1->area/c1->wall[2].length);
			tpc->deltaX=MIN(tpc->deltaX,c2->area/c2->wall[0].length);
			tpc->deltaX=MIN(tpc->deltaX,c2->area/c2->wall[1].length);
			tpc->deltaX=MIN(tpc->deltaX,c2->area/c2->wall[2].length);
			if(mesh->NCwall==4){
				tpc->deltaX=MIN(tpc->deltaX,c1->area/c1->wall[3].length);
				tpc->deltaX=MIN(tpc->deltaX,c2->area/c2->wall[3].length);	
			}						

			tpc->ccells[0]=cc1;
			tpc->ccells[1]=cc2;
			tpc->gcells[0]=c1;
			tpc->gcells[1]=c2;

			// Nuevo! Esto es para escribir en d_h[1] d_h[2] o d_h[3];
			tpc->i=c1->nneig;
			tpc->j=c2->nneig;
			//printf("c1: %d c2: %d i: %d j:%d a1: %lf a2: %lf\n",cc1->id,cc2->id,tpc->i,tpc->j,cc1->area,cc2->area);

			// Añadimos, a cada celda, el puntero a su vecina
			tpc->node1=parceldas[i].n1;
			tpc->node2=parceldas[i].n2;
			x1=mesh->nodes->nodes[parceldas[i].n1].x;
			x2=mesh->nodes->nodes[parceldas[i].n2].x;
			y1=mesh->nodes->nodes[parceldas[i].n1].y;
			y2=mesh->nodes->nodes[parceldas[i].n2].y;

			c1->neigcell[c1->nneig]=cc2;
			c1->neigwall[c1->nneig]=tpc;
			c1->nneig++;
			c2->neigcell[c2->nneig]=cc1;
			c2->neigwall[c2->nneig]=tpc;
			c2->nneig++;

			mesh->w_calc->n++;
			//printf("Pared %d/%d es de calculo (%d %d)\n",i,nparceldas,parceldas[i].idW1,parceldas[i].idW2);
	}
	mesh->nw_calc=mesh->w_calc->n;
    sprintf(temp,"Number of internal walls: %d",mesh->nw_calc);
    Notify(temp,MSG_L0,msg);  


	mesh->w_bound->n=0;
	for(i=0;i<(mesh->ncells*mesh->NCwall);i++){
		if(paredes[i].included==0){
			tpc=&(mesh->w_bound->wall[mesh->w_bound->n]);

            tpc->idWall = mesh->w_calc->n + mesh->w_bound->n;
            tpc->typeOfBound=0;

			c1=&mesh->g_cells->cells[paredes[i].idC];
			cc1=&mesh->c_cells->cells[paredes[i].idC];
			
			tpc->length=c1->wall[paredes[i].idW].length;
			tpc->lwall=c1->wall[paredes[i].idW].length;
			tpc->normal[_X_]=c1->wall[paredes[i].idW].normal[0];
			tpc->normal[_Y_]=c1->wall[paredes[i].idW].normal[1];

            tpc->ccells[0]=cc1;
			tpc->gcells[0]=c1;
            
            tpc->node1 = paredes[i].id1;  // Store nodes that define wall
            tpc->node2 = paredes[i].id2;

			
			c1->neigcell[c1->nneig+c1->nbound]=NULL;
			c1->neigwall[c1->nneig+c1->nbound]=tpc;
			c1->nbound++;

			mesh->w_bound->n++;
		}
	}
	mesh->nw_bound=mesh->w_bound->n;
    sprintf(temp,"Number of boundary walls: %d",mesh->nw_bound);
    Notify(temp,MSG_L0,msg);    


	//List of boundary cells
	mesh->ncellsBound=0;
	mesh->b_cells=(l_c_cells*)malloc(sizeof(l_c_cells));
	mesh->b_cells->pcells=(t_c_cell**)malloc(sizeof(t_c_cell*));
	for(i=0;i<mesh->nw_bound;i++){
        tpc=&(mesh->w_bound->wall[i]);

        //create typeOfBound flag for boundary cells
        if(mesh->NCwall==3){
            if(tpc->gcells[0]->nneig==1){ // corner!
                tpc->typeOfBound=1;
            }
        }else{//NCwall==4
            if(tpc->gcells[0]->nneig<=2){ // corner!
                tpc->typeOfBound=1;
            }
        }

        //add boundary cell to list
        found=0;
        for(j=0;j<mesh->ncellsBound;j++){
            if(mesh->b_cells->pcells[j]->id==tpc->gcells[0]->id){
                found=1;
            }
		}
        if(!found){ // We append the new found cell to the list of Boundary Cells
            tpc->ccells[0]->geom->isBound=99;
            mesh->b_cells->pcells[mesh->ncellsBound]=tpc->ccells[0];
            mesh->ncellsBound++;
            mesh->b_cells->pcells=(t_c_cell**)realloc(mesh->b_cells->pcells,(mesh->ncellsBound+1)*sizeof(t_c_cell*));
        }
	}
	mesh->b_cells->n=mesh->ncellsBound;     

	// Maximum wall length by cells
	for(i=0;i<mesh->ncells;i++){
		c1=&mesh->g_cells->cells[i];
		maxlwall=0.0;
		for(j=0;j<c1->nneig;j++){
			aux1=c1->neigwall[j]->length;
			maxlwall = MAX(aux1,maxlwall);
		}
		c1->Lwallmax = maxlwall;
		//printf("\nLwallmax %lf\n",maxlwall);
	}

	free(parceldas);
	free(nCeldasEnPunto);
	free(celdas);
	free(paredes);

    sprintf(temp,"Created new neighbouring-cells list");
    Notify(temp,MSG_L1,msg);       

	return(1);
}

/*******************************************************/
int cmpp(const void *id1, const void *id2){

    const t_pair *idx1,*idx2;

    idx1=(const t_pair*) id1;
    idx2=(const t_pair*) id2;

	if(idx1->id1<idx2->id1){
		return(-1);
	}else{
		if(idx1->id1>idx2->id1){
			return(1);
		}else{
			if(idx1->id2<=idx2->id2){
				return(-1);
			}else{
				return(1);
			}
		}
	}
}

/*******************************************************/
int check_angle(l_g_cells *lcells,t_message *msg){

	int i;
	t_g_cell *myCell;
	t_node node1,node2;
	char error[1024];
	double angle1,angle2,angle3, dot,mod1,mod2;

	for(i=0;i<(lcells->n);i++){
		myCell=&(lcells->cells[i]);
		node1.x=(myCell->nodes[1]->x-myCell->nodes[0]->x);
		node1.y=(myCell->nodes[1]->y-myCell->nodes[0]->y);
		node2.x=(myCell->nodes[1]->x-myCell->nodes[2]->x);
		node2.y=(myCell->nodes[1]->y-myCell->nodes[2]->y);

		dot=node1.x * node2.x + node1.y * node2.y;
		mod1=sqrt(node1.x*node1.x + node1.y*node1.y);
		mod2=sqrt(node2.x*node2.x + node2.y*node2.y);
		//printf("dot:%lf mod:%lf mod:%lf\n",dot,mod1,mod2);

		angle1=acos(dot/(mod1*mod2))*180/_PI_; //to centigrades
		//printf("angle1111111:%.9f\n",angle1);

		node1.x=(myCell->nodes[2]->x-myCell->nodes[0]->x);
		node1.y=(myCell->nodes[2]->y-myCell->nodes[0]->y);
		node2.x=(myCell->nodes[2]->x-myCell->nodes[1]->x);
		node2.y=(myCell->nodes[2]->y-myCell->nodes[1]->y);

		dot=node1.x * node2.x + node1.y * node2.y;
		mod1=sqrt(node1.x*node1.x + node1.y*node1.y);
		mod2=sqrt(node2.x*node2.x + node2.y*node2.y);

		angle2=acos(dot/(mod1*mod2))*180/_PI_; //to centigrades

		node1.x=(myCell->nodes[0]->x-myCell->nodes[1]->x);
		node1.y=(myCell->nodes[0]->y-myCell->nodes[1]->y);
		node2.x=(myCell->nodes[0]->x-myCell->nodes[2]->x);
		node2.y=(myCell->nodes[0]->y-myCell->nodes[2]->y);

		dot=node1.x * node2.x + node1.y * node2.y;
		mod1=sqrt(node1.x*node1.x + node1.y*node1.y);
		mod2=sqrt(node2.x*node2.x + node2.y*node2.y);

		angle3=acos(dot/(mod1*mod2))*180/_PI_; //to centigrades
		//printf("angle1:%lf angle2:%lf angle3:%lf\n",angle1,angle2,angle3);

		if(angle1<MAX_ANGLE_HYDRONIA_MESH || angle2<MAX_ANGLE_HYDRONIA_MESH || angle3 <MAX_ANGLE_HYDRONIA_MESH){
			sprintf(error,"Cell %d does not satisfy the angle constraint",i+1);
            Notify(error,MSG_ERROR,msg);
			return 0;
		}
	}

	return 1;
}


/*******************************************************/
int build_wall_inlet(t_mesh *mesh, t_bound *iinn, t_wall *w, int i, t_message *msg){
    int k,ncellsInlet;
    double vx1,vx2,vy1,vy2,vx3,vy3;
    double modulo;
    t_wall *w1,*w2;
    FILE *fp;
    int found;
    char error[1024];
    char temp[1024];

    ncellsInlet = iinn->ncellsBound;
    if(w->ccells[0]->z<iinn->zmin){
        iinn->zmin = w->ccells[0]->z;
        iinn->cellzMin = w->ccells[0];
    }
	// Si la celda que pertenece a la entrada es angulosa (Tiene solo
	// un vecino o dos en el caso de cuadradas), solo tenemos que añadir
	// a una de las paredes, y la información geometrica "equivalente" de
	// la condensación de las dos paredes en una sola, se modifica aquí.
	// OJO: Aun así, una de las paredes (la que no incluimos) debemos sacarla tambien
	// del calculo. No la incluimos pero modificamos el flag typeOfBound
	if(w->typeOfBound==0){
		ncellsInlet++;
		iinn->totalLength+=w->length;
		iinn->totalArea+=w->gcells[0]->area;
		w->typeOfBound=2;
		w->idBound=i+1;

		iinn->wallBound=(t_wall**)realloc(iinn->wallBound,ncellsInlet*sizeof(t_wall**));
		iinn->wallBound[ncellsInlet-1]=w;
		w->ccells[0]->geom->isBound = -1*(i+1);

        #if DEBUG_READ_BOUNDS
        sprintf(temp,"Inlet boundary cell %d assigned to open boundary ID %d",w->ccells[0]->id, w->idBound);
        Notify(temp,MSG_L0,msg);
        #endif
	
	}else{
		// Cuidado. Hay que incluir una. Solo una.
		// Verificamos primero si la celda no estaba incluida.
		k=0;
		found=0;
		while(k<ncellsInlet){
			if(iinn->wallBound[k]->ccells[0]->id == w->ccells[0]->id){// Ya la habiamos incluido antes
				found++;
				k=ncellsInlet;
			}
			k++;
		}
		// Ya hemos mirado si estaba o no estaba (found)
		if(found){ // Estaba ya incluida la celda
			w->typeOfBound=-2; // se pone esto para saber que es de entrada y evitar
								// que se hagan las correciones como si fuera de contorno
								// estándar. Se cambia a -2 por solutos en GPU

		}else{
			ncellsInlet++;
			// Modificamos la longitud y normal de la pared para adaptarlo a NORMALXI y
			// NORMALYI del codigo de javier
			// La longitud sera la de la pared de la celda que tiene una vecina
			if(mesh->NCwall==3){
				w->length=w->gcells[0]->neigwall[0]->length;
				// Necesitamos seleccionar la saliente. Para ello cambiamos el signo si
				// es necesario.
				if(w->ccells[0]->id==w->gcells[0]->neigwall[0]->gcells[0]->id){
					w->normal[_X_]=-w->gcells[0]->neigwall[0]->normal[_X_];
					w->normal[_Y_]=-w->gcells[0]->neigwall[0]->normal[_Y_];
				}else{
					w->normal[_X_]=w->gcells[0]->neigwall[0]->normal[_X_];
					w->normal[_Y_]=w->gcells[0]->neigwall[0]->normal[_Y_];
				}
			}else{
				if(w->gcells[0]->nneig==1){// Es una celda saliente
					// Necesitamos seleccionar la saliente. Para ello cambiamos el signo si
					// es necesario.
					if(w->ccells[0]->id==w->gcells[0]->neigwall[0]->gcells[0]->id){
						w->normal[_X_]=-w->gcells[0]->neigwall[0]->normal[_X_];
						w->normal[_Y_]=-w->gcells[0]->neigwall[0]->normal[_Y_];
					}else{
						w->normal[_X_]=w->gcells[0]->neigwall[0]->normal[_X_];
						w->normal[_Y_]=w->gcells[0]->neigwall[0]->normal[_Y_];
					}
					w->length=sqrt(3.0)*w->gcells[0]->neigwall[0]->length;
				}else{ // Es una celda "picuda". Tiene dos vecinas
                        // Vamos a sacar dos vectores (v1 y v2) de tal forma que
                        // el vector que buscamos, es la suma de ambos normalizado
                        /*
                        v2
                        ^
                        |
                        |
                    ---------
                    ||		|
                    ||		|-->v1
                    ||		|
                    =========
                    /
                    /
                    \/
                    v3
                    Las paredes con dos rayas, son contorno. v1 y v2 son los vectores normales salientes a las paredes
                    que tienen vecinas. El resultante de la suma de ambos, sera el normal que queremos
                    */
					w1=w->gcells[0]->neigwall[0];
					w2=w->gcells[0]->neigwall[1];
					// Sacamos el vector 1
					if(w1->ccells[0]->id==w->ccells[0]->id){
						// El vector uno, es el que queremos
						vx1=w1->normal[_X_];
						vy1=w1->normal[_Y_];
					}else{
						// El vector uno cambiado de signo, es el que queremos
						vx1=-w1->normal[_X_];
						vy1=-w1->normal[_Y_];
					}
					// Sacamos el vector 2
					if(w2->ccells[0]->id==w->ccells[0]->id){
						// El vector uno, es el que queremos
						vx2=w2->normal[_X_];
						vy2=w2->normal[_Y_];
					}else{
						// El vector uno cambiado de signo, es el que queremos
						vx2=-w2->normal[_X_];
						vy2=-w2->normal[_Y_];
					}
					vx3=-(vx1+vx2);
					vy3=-(vy1+vy2);
					modulo=sqrt(vx3*vx3+vy3*vy3);
					w->normal[_X_]=vx3/modulo;
					w->normal[_Y_]=vy3/modulo;
					w->length=sqrt(w1->length*w1->length + w2->length*w2->length);
				}
			}
            #if DEBUG_READ_BOUNDS
            sprintf(temp,"Modified boundary wall %d: New length %lf - New normal vector x(%lf) y(%lf)",w->ccells[0]->id, w->length,w->normal[_X_],w->normal[_Y_]);
            Notify(temp,MSG_WARN,msg);
            #endif

			iinn->totalLength+=w->length;
			iinn->totalArea+=w->gcells[0]->area;
			w->typeOfBound=2;
			w->idBound=i+1;

			iinn->wallBound=(t_wall**)realloc(iinn->wallBound,ncellsInlet*sizeof(t_wall**));
			iinn->wallBound[ncellsInlet-1]=w;
			w->ccells[0]->geom->isBound = -1*(i+1);

            #if DEBUG_READ_BOUNDS
            sprintf(temp,"Inlet boundary cell %d assigned to open boundary ID %d",w->ccells[0]->id, w->idBound);
            Notify(temp,MSG_L0,msg);
            #endif

		}

        //we have to check that the other walls in the boundary cell having typeofBound=1
        //are set to typeOfBound=-2
        found=0;
        w1=mesh->w_bound->wall;
        while(k<mesh->nw_bound && found<mesh->NCwall-2){
            if(found==0&&w1->ccells[0]->id==w->ccells[0]->id){//pared que es contorno;
                w1->typeOfBound=-2;
                found+=1;
            }
            k++;
            w1++;
        }

	}
	iinn->ncellsBound=ncellsInlet;

	//printf("in->totallength:%lf\n",iinn->totalLength);
	if(ncellsInlet==0){
		sprintf(error,"Cells in inlet %d not found",i);
     	Notify(error,MSG_ERROR,msg);
		return(0);
	}
 
    return(1);
}


/*******************************************************/
int build_inner_inlet(t_mesh *mesh, t_bound *iinn, int i, t_message *msg){

    int j;
    t_wall *w;
    t_c_cell *c1;
    t_g_cell *g1;
    t_c_cell *cAux;
    int k,l;
    int found;
    int nfound;
    int *founded;
    int ncellsInlet;

	founded=(int*)malloc(sizeof(int)* iinn->ncellsBound * 2);
	nfound=0;
	found=0;
	iinn->cellInner=0;
	iinn->cellInner=(t_c_cell**)malloc(sizeof(t_c_cell*));
	iinn->cellBound=(t_c_cell**)malloc(sizeof(t_c_cell*)* iinn->ncellsBound);

    ncellsInlet = iinn->ncellsBound;
    iinn->ncellsBound = 0;

	for(j=0; j<ncellsInlet; j++){

		w=iinn->wallBound[j];
		c1=w->ccells[0];
		g1=w->gcells[0];
		iinn->cellBound[iinn->ncellsBound]=c1;
		iinn->ncellsBound++;

		if(mesh->NCwall==3){
			for(k=0;k<g1->nneig;k++){ // Celda vecina de la celda de entrada
				cAux=g1->neigcell[k];

				if(g1->neigwall[k]->id==g1->id || g1->neigwall[k]->id2==g1->id){
					for(l=0;l<nfound;l++){
						if(cAux->id==founded[l]){
							found=1;
						}
					}
					if(	!found &&
						(cAux->geom->nneig==mesh->NCwall || 
						(cAux->geom->nneig!=mesh->NCwall && cAux->geom->isBound==0) )){

						cAux->geom->isBound = -1*(i);

						iinn->cellInner=(t_c_cell**)realloc(iinn->cellInner,sizeof(t_c_cell**)*(nfound+1));
						iinn->cellInner[nfound]=cAux;
						founded[nfound]=cAux->id;
						nfound++;
						iinn->ncellsInner=nfound;
					}
					found=0;
				}
			}
		}

	}
	free(founded);

    return(1);

}


/*******************************************************/
int build_wall_outlet(t_mesh *mesh, t_bound *outt, t_wall *w, int i,t_message *msg){
	int k,ncellsOutlet;
	t_wall *w1,*w2;
	FILE *fp;
	double vx1,vx2,vy1,vy2,vx3,vy3;
	double modulo;
	int found;
	char error[1024];
    char temp[1024];

    ncellsOutlet = outt->ncellsBound;
    if(w->ccells[0]->z<outt->zmin){
        outt->zmin=w->ccells[0]->z;
        outt->cellzMin=w->ccells[0];
    }
    // Si la celda que pertenece a la entrada es angulosa (Tiene solo
    // un vecino o dos en el caso de cuadradas), solo tenemos que añadir
    // a una de las paredes, y la información geometrica "equivalente" de
    // la condensación de las dos paredes en una sola, se modifica aquí.
    // OJO: Aun así, una de las paredes (la que no incluimos) debemos sacarla tambien
    // del calculo. No la incluimos pero modificamos el flag typeOfBound
    if(w->typeOfBound==0){
        ncellsOutlet++;
        outt->totalLength+=w->length;
        outt->totalArea+=w->gcells[0]->area;
        w->typeOfBound=3;
        w->idBound=mesh->nInlet+i+1;

        outt->wallBound=(t_wall**)realloc(outt->wallBound,ncellsOutlet*sizeof(t_wall**));
        outt->wallBound[ncellsOutlet-1]=w;
        w->ccells[0]->geom->isBound= i+1;

        #if DEBUG_READ_BOUNDS
        sprintf(temp,"Outlet boundary cell %d assigned to open boundary ID %d",w->ccells[0]->id, w->idBound);
        Notify(temp,MSG_L0,msg);
        #endif

    }else{
        // Cuidado. Hay que incluir una. Solo una.
        // Verificamos primero si la celda no estaba incluida.
        k=0;
        found=0;
        while(k<ncellsOutlet){
            if(outt->wallBound[k]->ccells[0]->id==w->ccells[0]->id){// Ya la habiamos incluido antes
                found++;
                k=ncellsOutlet;
            }
            k++;
        }
        // Ya hemos mirado si estaba o no estaba (found)
        if(found){ // Estaba ya incluida la celda
            w->typeOfBound=-3; // se pone esto para saber que es de salida y evitar
                                // que se hagan las correciones como si fuera de contorno
                                // estándar. Se cambia a -3 por solutos en GPU

        }else{
            ncellsOutlet++;
            // Modificamos la longitud y normal de la pared para adaptarlo a NORMALXI y
            // NORMALYI del codigo de javier
            // La longitud sera la de la pared de la celda que tiene una vecina
            if(mesh->NCwall==3){
                w->length=w->gcells[0]->neigwall[0]->length;
                // Necesitamos seleccionar la saliente. Para ello cambiamos el signo si
                // es necesario.
                if(w->ccells[0]->id==w->gcells[0]->neigwall[0]->gcells[0]->id){
                    w->normal[_X_]=-w->gcells[0]->neigwall[0]->normal[_X_];
                    w->normal[_Y_]=-w->gcells[0]->neigwall[0]->normal[_Y_];
                }else{
                    w->normal[_X_]=w->gcells[0]->neigwall[0]->normal[_X_];
                    w->normal[_Y_]=w->gcells[0]->neigwall[0]->normal[_Y_];
                        }
            }else{
                if(w->gcells[0]->nneig==1){// Es una celda saliente
                    // Necesitamos seleccionar la saliente. Para ello cambiamos el signo si
                    // es necesario.
                    if(w->ccells[0]->id==w->gcells[0]->neigwall[0]->gcells[0]->id){
                        w->normal[_X_]=-w->gcells[0]->neigwall[0]->normal[_X_];
                        w->normal[_Y_]=-w->gcells[0]->neigwall[0]->normal[_Y_];
                    }else{
                        w->normal[_X_]=w->gcells[0]->neigwall[0]->normal[_X_];
                        w->normal[_Y_]=w->gcells[0]->neigwall[0]->normal[_Y_];
                    }
                    w->length=sqrt(3.0)*w->gcells[0]->neigwall[0]->length;
                }else{ // Es una celda "picuda". Tiene dos vecinas
                        // Vamos a sacar dos vectores (v1 y v2) de tal forma que
                        // el vector que buscamos, es la suma de ambos normalizado
                        /*
                        v2
                        ^
                        |
                        |
                    ---------
                    ||		|
                    ||		|-->v1
                    ||		|
                    =========
                    /
                    /
                    \/
                    v3
                    Las paredes con dos rayas, son contorno. v1 y v2 son los vectores normales salientes a las paredes
                    que tienen vecinas. El resultante de la suma de ambos, sera el normal que queremos
                    */
                    w1=w->gcells[0]->neigwall[0];
                    w2=w->gcells[0]->neigwall[1];
                    // Sacamos el vector 1
                    if(w1->ccells[0]->id==w->ccells[0]->id){
                        // El vector uno, es el que queremos
                        vx1=w1->normal[_X_];
                        vy1=w1->normal[_Y_];
                    }else{
                        // El vector uno cambiado de signo, es el que queremos
                        vx1=-w1->normal[_X_];
                        vy1=-w1->normal[_Y_];
                    }
                    // Sacamos el vector 2
                    if(w2->ccells[0]->id==w->ccells[0]->id){
                        // El vector uno, es el que queremos
                        vx2=w2->normal[_X_];
                        vy2=w2->normal[_Y_];
                    }else{
                        // El vector uno cambiado de signo, es el que queremos
                        vx2=-w2->normal[_X_];
                        vy2=-w2->normal[_Y_];
                    }
                    vx3=-(vx1+vx2);
                    vy3=-(vy1+vy2);
                    modulo=sqrt(vx3*vx3+vy3*vy3);
                    w->normal[_X_]=vx3/modulo;
                    w->normal[_Y_]=vy3/modulo;
                    w->length=sqrt(w1->length*w1->length + w2->length*w2->length);

                }
            }//endif mesh nWall
            #if DEBUG_READ_BOUNDS
            sprintf(temp,"Modified boundary wall %d: New length %lf - New normal vector x(%lf) y(%lf)",w->ccells[0]->id,w->length,w->normal[_X_],w->normal[_Y_]);
            Notify(temp,MSG_WARN,msg);
            #endif

            outt->totalLength+=w->length;
            outt->totalArea+=w->gcells[0]->area;
            w->typeOfBound=3;
            w->idBound=mesh->nInlet+i+1;

            outt->wallBound=(t_wall**)realloc(outt->wallBound,ncellsOutlet*sizeof(t_wall**));
            outt->wallBound[ncellsOutlet-1]=w;
            w->ccells[0]->geom->isBound=i+1;
            
            #if DEBUG_READ_BOUNDS
            sprintf(temp,"Outlet boundary cell %d assigned to open boundary ID %d",w->ccells[0]->id, w->idBound);
            Notify(temp,MSG_L0,msg);
            #endif
        }

        //we have to check that the other walls in tghe boundary cell having typeofBound=1
        //are set to typeOfBound=3
        found=0;
        w1=mesh->w_bound->wall;
        while(k<mesh->nw_bound && found<mesh->NCwall-2){
            if(found==0&&w1->ccells[0]->id==w->ccells[0]->id){//pared que es contorno;
                w1->typeOfBound=-3;
                found+=1;
            }
            k++;
            w1++;
        }

    }
	outt->ncellsBound=ncellsOutlet;

	if(ncellsOutlet==0){
        sprintf(error,"Cells in outlet %d not found",i);
        Notify(error,MSG_ERROR,msg);
        return(0);
	}

   return(1);
}


/*******************************************************/
int build_inner_outlet(t_mesh *mesh,t_bound *outt, int i, t_message *msg){

    int j;
    t_wall *w;
    t_c_cell *c1;
    t_g_cell *g1;
    t_c_cell *cAux;
    int k,l;
    int found;
    int nfound;
    int *founded;
    int ncellsOutlet;


    founded=(int*)malloc(sizeof(int)* outt->ncellsBound*2);
    nfound=0;
    found=0;
    outt->cellInner=0;
    outt->cellInner=(t_c_cell**)malloc(sizeof(t_c_cell*));
    outt->cellBound=(t_c_cell**)malloc(sizeof(t_c_cell*)*outt->ncellsBound);

    ncellsOutlet = outt->ncellsBound;
    outt->ncellsBound = 0;

	for(j=0; j<ncellsOutlet ;j++){

		w=outt->wallBound[j];
		c1=w->ccells[0];
		g1=w->gcells[0];
		outt->cellBound[outt->ncellsBound]=c1;
		outt->ncellsBound++;
		
		if(mesh->NCwall==3){
			for(k=0;k<g1->nneig;k++){
				cAux=g1->neigcell[k];

				if(g1->neigwall[k]->id==g1->id||g1->neigwall[k]->id2==g1->id){
					for(l=0;l<nfound;l++){
						if(cAux->id==founded[l]){
							found=1;
						}
					}
					if(	!found &&
						(cAux->geom->nneig==mesh->NCwall ||
						(cAux->geom->nneig!=mesh->NCwall&&cAux->geom->isBound==0) )){

						cAux->geom->isBound=i+1;
						
						outt->cellInner=(t_c_cell**)realloc(outt->cellInner,sizeof(t_c_cell**)*(nfound+1));
						outt->cellInner[nfound]=cAux;
						founded[nfound]=cAux->id;
						nfound++;
						outt->ncellsInner=nfound;
					}
					found=0;
				}
			}
		}
	}
	free(founded);

    return(1);

}






