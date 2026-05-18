/*Given a mesh in triangle format, it generates the *.area file to refine with triangle)*/
/* Author: Mario Morales*/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
 
#define MIN(x,y) (x < y ? x : y)            
#define MAX(x,y) (x > y ? x : y)
#define PI 3.1415926

typedef struct{
	double x, y, z;
} point;
 
 
typedef struct{
	int n1, n2, n3;  //n1,n2 y n3 -> node's id
	double cx, cy, z; //cx = cell's center->x    cy = cell's center->y
   double size;
} cell;
 
typedef struct{
	int npoints; //number of points
	point *pPol; //points
	double ref;//required refinement for this polygon
} polygon;

typedef struct{
	int ncells;
	cell *cell;
	int nnodes;
	point *node;
	int nPolygons;
	polygon *pol;
	char meshName[1024];
	char polName[1024];
} mesh;


void read_mesh(mesh *m);
int read_polygons(mesh *m);
int dentro_poligono (mesh *m, double coord_x, double coord_y, polygon *polyG);
void refine(mesh *m);



int main (int argc, char *argv[]){
 
	mesh *meshMain;

	meshMain=(mesh*) malloc(sizeof(mesh));

	if(argc<3){
		printf("Error. Usage: ./refine meshFile polygonsFile\n");
		return 0;
	}

	sprintf(meshMain->meshName,"%s",argv[1]);
	sprintf(meshMain->polName,"%s",argv[2]);

	read_mesh(meshMain);
	read_polygons(meshMain);
	refine(meshMain);

	printf("%s.area generated\n",meshMain->meshName);


	return 0;
 
 
}





////////////////////////////////////////////////
//This function reads the triangle's output files and builds de mesh
void read_mesh(mesh *m){

	int i, j;
	double x1;
	cell *celda;
	point *nodo;
	point *n1, *n2, *n3;
	char fname2[1024];
	FILE *fp;
 
	sprintf(fname2, "%s.node", m->meshName);
	fp = fopen (fname2, "r");
 
	fscanf(fp, "%d %*d %*d %*d", &m->nnodes);
  
	m->node = (point*) malloc(m->nnodes*sizeof(point));
	
	nodo=m->node;
	for(i=0;i<m->nnodes;i++){
		fscanf(fp, "%*d %lf %lf %*lf", &nodo->x,&nodo->y); 
		nodo++;
	}
 
	fclose(fp);

	sprintf(fname2, "%s.ele", m->meshName);
	fp = fopen (fname2, "r");
 
	fscanf(fp, "%d %*d %*d", &m->ncells);
  
	m->cell = (cell*) malloc(m->ncells*sizeof(cell));
 
	celda=m->cell;
	for(i=0;i<m->ncells;i++){
 
		fscanf(fp, "%*d %d %d %d", &celda->n1,&celda->n2,&celda->n3); 
		celda->n1--;
		celda->n2--;
		celda->n3--;
 
		j=celda->n1;
		n1=&(m->node[j]);
		j=celda->n2;
		n2=&(m->node[j]);
		j=celda->n3;
		n3=&(m->node[j]);
 
		celda->cx=(n1->x+n2->x+n3->x)/3.;
		celda->cy=(n1->y+n2->y+n3->y)/3.;
		celda++;
 
	}
 
	fclose(fp);

}


int read_polygons(mesh *m){

	FILE *fp;
	int i,j;

	fp=fopen(m->polName,"r");

	fscanf(fp,"%d",&m->nPolygons);
	if(m->nPolygons<=0){
		printf("Error: Number of polygons is set to %d\n",m->nPolygons);
		return 0;
	}
	
	m->pol=(polygon*) malloc(m->nPolygons*sizeof(polygon));

	for(i=0;i<m->nPolygons;i++){
		fscanf(fp,"%d %lf",&m->pol[i].npoints,&m->pol[i].ref);
		//printf("%d npoints:%d ref:%d\n",i+1,m->pol[i].npoints,m->pol[i].ref);
		if(m->pol[i].npoints<3){
			printf("Error: Number of points defined for polygon %d is set to %d\n",i+1,m->pol[i].npoints);
			return 0;
		}
		if(m->pol[i].ref<1e-12){
			printf("Error: Refinement for polygon %d is set to %.12f\n",i+1,m->pol[i].ref);
			return 0;
		}

		m->pol[i].pPol=(point*) malloc(m->pol[i].npoints*sizeof(point));
		for(j=0;j<m->pol[i].npoints;j++){
			fscanf(fp,"%lf %lf",&m->pol[i].pPol[j].x,&m->pol[i].pPol[j].y);
		}
	}


	fclose(fp);

	return 1;


}


// Devuelve 0 si esta fuera y 1 si esta dentro
 
int dentro_poligono (mesh *m, double coord_x, double coord_y, polygon *polyG){

   int contador;
   int i,j,k;
   double x_inters,trash;
	FILE *fp;
	point *p1,*p2;
   
	j=polyG->npoints;
	
	contador=0;
	for (i=0;i<j;i++){
		p1 = polyG->pPol + i;	
		p2 = polyG->pPol + (i+1)%j;
		if ((coord_y > MIN(p1->y,p2->y)) && (coord_y <= MAX(p1->y,p2->y)) && (coord_x <= MAX(p1->x,p2->x))) {
			if (p1->y != p2->y) {
				x_inters = p1->x + (coord_y-p1->y)*(p2->x-p1->x)/(p2->y-p1->y);
				if (p1->x == p2->x || coord_x <= x_inters){
					contador++;
				}
			}		
		}
	}

	if (contador % 2 == 0){
		return(0);
	}else{
		return(1);
	}

}


void refine(mesh *m){

	int i,k;
	cell *celda;
	double cellRef; //cellRefinement
	FILE *fp1;
	char fname[1024];
	polygon *polyG;

	sprintf(fname, "%s.area", m->meshName);
	fp1 = fopen (fname, "w");

	celda=m->cell;
	fprintf(fp1,"%d\n",m->ncells);
	for(i=0;i<m->ncells;i++){
		cellRef=-1.0;
		polyG=m->pol;
		for(k=0;k<m->nPolygons;k++){
			if(dentro_poligono(m,celda->cx,celda->cy,polyG)==1){
				if(fabs(cellRef+1.0)<1e-12){
					cellRef=polyG->ref;
				}else{
					cellRef=MIN(polyG->ref,cellRef);
				}

			}
			polyG++;
		}
		fprintf(fp1,"%d %lf\n",i+1,cellRef);
		celda++;
	}

	fclose(fp1);


}
