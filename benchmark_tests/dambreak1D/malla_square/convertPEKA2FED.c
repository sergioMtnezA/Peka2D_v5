#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#define OK " \033[1;35m[OK]\033[0m "
#define PI 3.1415926535
#define MIN(x,y) (x < y ? x : y)
#define MAX(x,y) (x > y ? x : y)




////////////////////////////////////////////////////
typedef struct{
	double x,y,z;
	int id;
}t_node;


typedef struct{
	int id;
	int n1,n2,n3,n4;
	double cx,cy; //center
	double area; //area
	double z,zR;
	double h, wsl;
	double nMan;
	double d50;
}t_cell;


typedef struct{
	int wallXcell;
	int ncells;
	int nnodes;
	t_cell *cell;
	t_node *node;
}t_mesh;


typedef struct{
	int nPolys;
	int *npoints;
	t_node **poly;
	double *manPoly;
} t_polys;




//////////////////////////////////////////////////////////////////////////
double triangle_size(t_node *n1, t_node *n2, t_node *n3){

	double x1, x2, x3, x4;

	x1 = n1->x - n2->x;
	x2 = n1->y - n2->y;
	x3 = n1->x - n3->x;
	x4 = n1->y - n3->y;

	return fabs(0.5*(x1*x4-x2*x3));

}




int dentro_poligono(t_node *poly,int npoints, t_node point) {

	int contador;
	int i, j;
	double x_inters;
	t_node *p1, *p2;
 	double coord_y, coord_x;

	contador=0;

   coord_x=point.x;
   coord_y=point.y;
	j=npoints;

	for (i=0;i<j;i++) {

		p1 = &(poly[i]);
		p2 = &(poly[(i+1)%j]);
		if ((coord_y > MIN(p1->y,p2->y)) && (coord_y <= MAX(p1->y,p2->y)) && (coord_x <= MAX(p1->x,p2->x))) {
        if (p1->y != p2->y) {
		    x_inters = p1->x + (coord_y-p1->y)*(p2->x-p1->x)/(p2->y-p1->y);
			 if (p1->x == p2->x || coord_x <= x_inters){
					contador++;
          }
		  }
		}
	}
	if (contador % 2 == 0)
		return(0);
	else
		return(1);

}









int read_peka_mesh(t_mesh *mesh, char *filename){

	int i;
	FILE *fp;
	t_node *n1,*n2,*n3,*n4;

	fp=fopen(filename,"r");

	if(!fp){
		printf("File %s not found\n",filename);
		return(0);
	}

	fscanf(fp,"%*s %d",&mesh->wallXcell);
	fscanf(fp,"%*s %d\n",&mesh->ncells);
	fscanf(fp,"%*s %d\n",&mesh->nnodes);

	mesh->cell=(t_cell*) malloc(mesh->ncells*sizeof(t_cell));
	mesh->node=(t_node*) malloc(mesh->nnodes*sizeof(t_node));

	for(i=0;i<mesh->nnodes;i++){
		mesh->node[i].id=i+1;
		fscanf(fp,"%lf %lf",&mesh->node[i].x,&mesh->node[i].y);
	}

	for(i=0;i<mesh->ncells;i++){
        if(mesh->wallXcell==3){
		    fscanf(fp,"%d %d %d",
                &mesh->cell[i].n1,
                &mesh->cell[i].n2,
                &mesh->cell[i].n3); // Los id de los nodos empiezan en 1
        }else if(mesh->wallXcell==4){
            fscanf(fp,"%d %d %d %d",
                &mesh->cell[i].n1,
                &mesh->cell[i].n2,
                &mesh->cell[i].n3,
                &mesh->cell[i].n4); // Los id de los nodos empiezan en 1
        }

	}

	for(i=0;i<mesh->ncells;i++){
		fscanf(fp,"%lf",&mesh->cell[i].zR);
	}

	fclose(fp);




	// Nodes ID
	for(i=0;i<mesh->nnodes;i++){
		mesh->node[i].id=i+1;
	}

	// Cell ID and nMan
	for(i=0;i<mesh->ncells;i++){
		mesh->cell[i].id=i+1;	
		mesh->cell[i].nMan=0.018;
	}




	return 1;

}




int read_peka_inicial(t_mesh *mesh, char *filename){

	int i;
	FILE *fp;
	t_node *n1,*n2,*n3;
	double aux;

	fp=fopen(filename,"r");

	if(!fp){
		printf("File %s not found\n",filename);
		return(0);
	}

	fscanf(fp,"%*s %*lf");
	fscanf(fp,"%*s %*lf");
	fscanf(fp,"%*s %*lf");
	fscanf(fp,"%*s %*d");

	//.inicial
	/*for(i=0;i<mesh->ncells;i++){
		fscanf(fp,"%lf %*lf %*lf %lf %*lf\n",&mesh->cell[i].h,&mesh->cell[i].z);
	}*/
	//.ini
	for(i=0;i<mesh->ncells;i++){
		fscanf(fp,"%lf %*lf %*lf %lf\n",&mesh->cell[i].h,&mesh->cell[i].z);
	}

	fclose(fp);

	// Cell WSL
	for(i=0;i<mesh->ncells;i++){
		mesh->cell[i].wsl=mesh->cell[i].h+mesh->cell[i].z;
	}

	return 1;

}


int read_EHC_ini(t_mesh *mesh, char *filename){

	int i;
	FILE *fp;
	t_node *n1,*n2,*n3;
	double aux;

	fp=fopen(filename,"r");

	if(!fp){
		printf("File %s not found\n",filename);
		return(0);
	}

	fscanf(fp,"%*s %*lf");
	fscanf(fp,"%*s %*lf");
	fscanf(fp,"%*s %*lf");
	fscanf(fp,"%*s %*d");


	//.ini
	for(i=0;i<mesh->ncells;i++){
		fscanf(fp,"%lf %*lf %*lf %lf\n",&mesh->cell[i].h,&mesh->cell[i].z);
	}

	fclose(fp);

	// Cell WSL
	for(i=0;i<mesh->ncells;i++){
		mesh->cell[i].wsl=mesh->cell[i].h+mesh->cell[i].z;
	}

	return 1;

}







int write_fed_mesh(t_mesh *mesh, char *filename){

	int i;
	FILE *fp;
	t_node *n1,*n2,*n3;

	fp=fopen(filename,"w");


	fprintf(fp,"%8d %8d %d 0\n",mesh->ncells,mesh->nnodes,mesh->wallXcell);


	for(i=0; i<mesh->nnodes; i++){
		fprintf(fp,"%8d %12.6lf %12.6lf    0.000    0.000    -9999.000    0    0\n",
            mesh->node[i].id,
            mesh->node[i].x,
            mesh->node[i].y);
	}

	printf("%s Nodes written %d\n",OK,mesh->nnodes);

	for(i=0; i<mesh->ncells; i++){
        if(mesh->wallXcell==3){
            fprintf(fp,"%8d %8d %8d %8d %12.3lf %12.6lf %12.6lf %12.6lf    0.000\n",
                mesh->cell[i].id,
                mesh->cell[i].n1,
                mesh->cell[i].n2,
                mesh->cell[i].n3,
                mesh->cell[i].nMan,
                mesh->cell[i].z,
                mesh->cell[i].wsl,
                mesh->cell[i].zR);
        }else if(mesh->wallXcell==4){
            fprintf(fp,"%8d %8d %8d %8d %8d %12.3lf %12.6lf %12.6lf %12.6lf    0.000\n",
                mesh->cell[i].id,
                mesh->cell[i].n1,
                mesh->cell[i].n2,
                mesh->cell[i].n3,
                mesh->cell[i].n4,
                mesh->cell[i].nMan,
                mesh->cell[i].z,
                mesh->cell[i].wsl,
                mesh->cell[i].zR);
        }
    }

	printf("%s Cells written %d\n",OK,mesh->ncells);

	fclose(fp);


	return 1;

}






int read_polys(t_polys *polys, char *filename){

	int i,j,nPolys;
	FILE *fp;

	fp=fopen(filename,"r");

	if(!fp){
		printf("Poly file %s not found",filename);
		return(0);
	}


	fscanf(fp,"%d",&polys->nPolys);

	polys->poly=(t_node**) malloc(polys->nPolys*sizeof(t_node*));
	polys->npoints=(int*) malloc(polys->nPolys*sizeof(int));
	polys->manPoly=(double*) malloc(polys->nPolys*sizeof(double));


	for(i=0;i<polys->nPolys;i++){

		fscanf(fp,"%lf %d",&polys->manPoly[i],&polys->npoints[i]);

		polys->poly[i]=(t_node*) malloc(polys->npoints[i]*sizeof(t_node));

		for(j=0;j<polys->npoints[i];j++){
			fscanf(fp,"%lf %lf",&polys->poly[i][j].x,&polys->poly[i][j].y);
		}


		printf("%s Poly %d read %d\n",OK,i+1,polys->npoints[i]);
	}

	fclose(fp);


	return 1;
}




void assign_Manning(t_mesh *mesh, t_polys *polys){


	int i,j,nPolys;
	t_cell *cell;
	t_node nodeAux;


	cell=mesh->cell;
	for(i=0;i<mesh->ncells;i++){
		nodeAux.x=cell->cx;
		nodeAux.y=cell->cy;

		for(j=0;j<polys->nPolys;j++){
			if(dentro_poligono(polys->poly[j],polys->npoints[j],nodeAux)==1){ // Nota: los poligonos deben estar ordenados en el archivo

				if(polys->manPoly[j] == 0.110) {
					cell->nMan = polys->manPoly[j];
					cell->d50 = 0.0074;
					cell->z = cell->zR + 5.0;
				} else if(polys->manPoly[j] == 0.040) {
					cell->nMan = 0.040;
					cell->d50 = 0.0;
					cell->z = cell->zR;
				} else if(polys->manPoly[j] == 1.0) {
					cell->nMan = 0.0283;
					cell->d50 = 0.0438;
					cell->z = cell->zR + 2.5;
				} else if(polys->manPoly[j] == 2.0) {
					cell->nMan = 0.0288;
					cell->d50 = 0.0486;
					cell->z = cell->zR + 2.5;
				} else if(polys->manPoly[j] == 3.0) {
					cell->nMan = 0.0297;
					cell->d50 = 0.0591;
					cell->z = cell->zR + 2.5;
				} else if(polys->manPoly[j] == 4.0) {
					cell->nMan = 0.0286;
					cell->d50 = 0.0469;
					cell->z = cell->zR + 2.5;
				} else if(polys->manPoly[j] == 5.0) {
					cell->nMan = 0.0295;
					cell->d50 = 0.0568;
					cell->z = cell->zR + 2.5;
				} else if(polys->manPoly[j] == 6.0) {
					cell->nMan = 0.0295;
					cell->d50 = 0.0567;
					cell->z = cell->zR + 2.5;
				}

			}
		}

		if(i%1000 == 0) printf("cell: %d\n",i);

		cell++;

	}



	// free polys
	for(j=0;j<polys->nPolys;j++){
		free(polys->poly[j]);
	}

	free(polys->poly);
	free(polys->manPoly);
	free(polys->npoints);



}





int write_manning_list(t_mesh *mesh, char *filename){

	int i;
	FILE *fp;

	fp=fopen(filename,"w");

	fprintf(fp,"%d\n",mesh->ncells);

	for(i=0;i<mesh->ncells;i++){
			fprintf(fp,"%lf\n",mesh->cell[i].nMan);
	}

	fclose(fp);

	return 1;

}



int write_d50_list(t_mesh *mesh, char *filename){

	int i;
	FILE *fp;

	fp=fopen(filename,"w");

	fprintf(fp,"%d\n",mesh->ncells);

	for(i=0;i<mesh->ncells;i++){
			fprintf(fp,"%lf\n",mesh->cell[i].d50);
	}

	fclose(fp);

	return 1;

}





void volcado_vtk(t_mesh *s,char *filename){


	int i,j;
	t_cell *celda;
	t_node *nodo;
	FILE *fp;

	fp = fopen (filename,"w");
	fprintf(fp,"# vtk DataFile Version 2.0\n");
	fprintf(fp,"mesh geometry\n");
	fprintf(fp,"ASCII\n");
	fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");

	// Imprimimos la informacion correspondiente a los nodos


   	fprintf(fp,"POINTS %d float \n",s->nnodes);
    nodo=s->node;
	for (i=0;i<s->nnodes;i++){
		fprintf(fp,"%lf %lf 0.0\n", nodo->x,nodo->y);
	    nodo++;
   	}

	// Imprimimos en el mismo fichero la informacion correspondiente a las celdas
    fprintf(fp,"CELLS %d %d \n", s->ncells, s->ncells*(s->wallXcell+1));
    celda=s->cell;
    for(i=0;i<s->ncells;i++){
        fprintf(fp,"%d ",s->wallXcell);
        if(s->wallXcell==3){
            fprintf(fp,"%d %d %d \n",celda->n1-1,celda->n2-1,celda->n3-1);
        }else if(s->wallXcell==4){
            fprintf(fp,"%d %d %d %d\n",celda->n1-1,celda->n2-1,celda->n3-1,celda->n4-1);
        }
        celda++;
    }

	fprintf(fp,"CELL_TYPES %d \n",s->ncells);
	for(i=0;i<s->ncells;i++){
        if(s->wallXcell==3){
            fprintf(fp,"5\n");
        }else{
            fprintf(fp,"9\n");
        }
	}


    //Imprimimos informacion sobre las variables fisicas del fluido
    fprintf (fp, "CELL_DATA %d \n", s->ncells);

    fprintf (fp, "SCALARS z float \n");
    fprintf (fp, "LOOKUP_TABLE default \n");
    celda=s->cell;
    for (i=0; i<s->ncells; i++){
    fprintf (fp,"%.6lf \n", celda->z);
    celda++;
    }

    fprintf (fp, "SCALARS h float \n");
    fprintf (fp, "LOOKUP_TABLE default \n");
    celda=s->cell;
    for (i=0; i<s->ncells; i++){
        fprintf (fp,"%.6lf \n", celda->h);
        celda++;
    }


	fclose(fp);

}







int main (int argc, char * argv[]){

	t_mesh *mesh;
	t_polys *polys;
	char filename[1024];

	if(argc!=2){
		printf("Usage %s projectName \n",argv[0]);
		printf("Example: ./convertPEKA2FED case1\n");
		return 0;
	}

	//inicializar la malla
	mesh=(t_mesh*) malloc(sizeof(t_mesh));

	//leer malla
	//sprintf(filename,"%s.txt",argv[1]);
	sprintf(filename,"%s.msh",argv[1]);
	read_peka_mesh(mesh,filename);

	printf("%s Mesh read\n",OK);


	//ller cond. iniciales
	sprintf(filename,"%s.ini",argv[1]);
	read_peka_inicial(mesh,filename);
	//sprintf(filename,"%s.ini",argv[1]);
	//read_EHC_ini(mesh,filename);
	

	printf("%s Initial conditions read\n",OK);

	// asociar manning a las celdas
	//polys=(t_polys*) malloc(sizeof(t_polys));
	//sprintf(filename,"new_manning_poligons/manning_pol_geo.txt");
	//read_polys(polys, filename);
	//assign_Manning(mesh, polys);

	sprintf(filename,"%s.FED",argv[1]);
	write_fed_mesh(mesh, filename);

	printf("%s FED mesh written\n",OK);

	//sprintf(filename,"%s.d50",argv[1]);
	//write_d50_list(mesh, filename);

	//printf("%s Data written\n",OK);


	//sprintf(filename,"%s.inicial",argv[1]);
	//write_inicial(mesh,filename);

	//printf("%s Inicial written\n",OK);


	sprintf(filename,"%s.vtk",argv[1]);
	volcado_vtk(mesh, filename);


	//printf("%s Dump results in vtk\n",OK);


	return 1;

}
