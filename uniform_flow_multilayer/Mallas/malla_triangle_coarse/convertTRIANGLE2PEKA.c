#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define MIN(x,y) (x < y ? x : y)
#define MAX(x,y) (x > y ? x : y)
#define PI 3.1415926
#define tand(x) (tan(fmod((x),360) * M_PI / 180))
#define Nsed 1

typedef struct{
	double x, y, z;
} tipo_nodo;

typedef struct{
	long n1, n2, n3;
	double cx, cy, z, epsilon;
	double h, qx, qy, ux, uy, q, u;
   double size;
   double phi;
} tipo_celda;

typedef struct{
	long nceldas;
	tipo_celda *celda;
	long nnodos;
	tipo_nodo *nodo;
	long ncontorno;
	tipo_nodo *contorno;
} tipo_sim;

typedef struct{
	int n;
	tipo_nodo *nodes;
} l_nodes;


//////////////////////////////////


double triangle_size(tipo_nodo *n1, tipo_nodo *n2, tipo_nodo *n3);


///////////////////////////////////////

void lector_triangle(tipo_sim *s, char *fname){

	long i, j;
	double x1,R;
	long nceldas, nnodos;
	tipo_celda *celda;
	tipo_nodo *nodo;
	tipo_nodo *n1, *n2, *n3;
	char fname2[1024];
	double slope,lambdaX,lambdaY,Amplitude;
	double x0,y0,r,x0_2,y0_2;
	double aux1, aux2, aux3;
	FILE *fp;

	sprintf(fname2, "%s.node", fname);
	fp = fopen (fname2, "r");

	fscanf(fp, "%ld", &nnodos);
	s->nnodos=nnodos;

	fscanf(fp, "%ld", &j);
	fscanf(fp, "%ld", &j);
	fscanf(fp, "%ld", &j);

	s->nodo = (tipo_nodo*) malloc(nnodos*sizeof(tipo_nodo));
	nodo=s->nodo;
	for(i=0;i<nnodos;i++){

		fscanf(fp, "%ld", &j);
		fscanf(fp, "%lf", &(nodo->x));
		fscanf(fp, "%lf", &(nodo->y));
		fscanf(fp, "%ld", &j);
		//fscanf(fp, "%ld", &j);

		nodo++;

	}

	fclose(fp);


	sprintf(fname2, "%s.ele", fname);
	fp = fopen (fname2, "r");

	fscanf(fp, "%ld", &nceldas);

	s->nceldas=nceldas;

	fscanf(fp, "%ld", &j);
	fscanf(fp, "%ld", &j);

	s->celda = (tipo_celda*) malloc(nceldas*sizeof(tipo_celda));

	celda=s->celda;
	for(i=0;i<nceldas;i++){

		fscanf(fp, "%ld", &j);
		fscanf(fp, "%ld", &(celda->n1));
		fscanf(fp, "%ld", &(celda->n2));
		fscanf(fp, "%ld", &(celda->n3));

		celda->n1--;
		celda->n2--;
		celda->n3--;

		j=celda->n1;
		n1=&(s->nodo[j]);
		j=celda->n2;
		n2=&(s->nodo[j]);
		j=celda->n3;
		n3=&(s->nodo[j]);

		celda->cx=(n1->x+n2->x+n3->x)/3.;
		celda->cy=(n1->y+n2->y+n3->y)/3.;

		celda->size=triangle_size(n1,n2,n3);
		celda++;

	}


	fclose(fp);

}

///////////////////////////////////////

//Volcar malla formato JM

void generate_peka_mesh(tipo_sim *s,char *fname){

   long i;
   tipo_celda *celda;
   tipo_nodo *nodo;
   FILE *fp;
   celda=s->celda;

	fp=fopen(fname, "w");
	fprintf(fp,"NVDRT 3 \n");
	fprintf(fp,"NCDLL %ld \n",s->nceldas);
	fprintf(fp,"NNODD %ld \n",s->nnodos);

	nodo=s->nodo;
	for (i=0 ; i<s->nnodos ; i++){
		fprintf(fp,"%lf %lf \n", nodo->x, nodo->y);
		nodo++;
	}

	celda=s->celda;
	for(i=0 ; i<s->nceldas ; i++){
      fprintf(fp,"%ld %ld %ld \n", celda->n1+1, celda->n2+1, celda->n3+1);

	   celda++;
   }

	celda=s->celda;
	for(i=0 ; i<s->nceldas ; i++){

      fprintf(fp,"0.0 \n");

      celda++;
	}

   fclose(fp);


}
//////////////////////////////////////////////////////////////////////



void generate_peka_ini (tipo_sim *s, char *fname){

	long i;
	double r;
	tipo_celda *celda;
	FILE *fp;

	double aux1, aux2, aux3;

	fp=fopen(fname,"w");
	fprintf(fp,"NIVEL_SUPERFICIAL_INICIAL   0.0 \n");
	fprintf(fp,"VELOCIDADX_INICIAL          0.0 \n");
	fprintf(fp,"VELOCIDADY_INICIAL          0.0 \n");
	fprintf(fp,"LEER_CELDAS                 4\n");

	celda=s->celda;
	for (i=0 ; i<s->nceldas ; i++){

			//calado
			celda->h = 1.516;

			//pendiente
			celda->z = 0.005*1000.0-0.005*celda->cx;


			fprintf(fp,"%lf 0.0 0.0 %lf\n",celda->h, celda->z);


		//fprintf(fp,"%lf 1.0 1.0\n",celda->h);
		celda++;
	}

	fclose(fp);

}
/////////////////////////////////////////////////////////////////


void generate_peka_solini (tipo_sim *s, char *fname){

	long i;
	double r;
	tipo_celda *celda;
	FILE *fp;

	fp=fopen(fname,"w");

	celda=s->celda;
	for (i=0 ; i<s->nceldas ; i++){

			celda->phi=5.0;
			if(celda->cx>=100. && celda->cx<= 150.){
				if(celda->cy>=-25. && celda->cy<=25.){
					celda->phi=10.0;
				}
				}
							

			fprintf(fp,"%lf\n",celda->phi);


		celda++;
	}

	fclose(fp);

}
/////////////////////////////////////////////////////////////////////



void volcado_vtk(tipo_sim *s,char *fname){


	long i;
	tipo_celda *celda;
	tipo_nodo *nodo;
	FILE *fp;

	int j;

	celda=s->celda;
	nodo=s->nodo;
	fp = fopen (fname,"w");	fprintf(fp,"# vtk DataFile Version 2.0 \n");
	fprintf(fp,"Titulo \n");
	fprintf(fp,"ASCII \n");
	fprintf(fp,"DATASET UNSTRUCTURED_GRID \n");

	// Imprimimos la informacion correspondiente a los nodos


   	fprintf(fp,"POINTS %ld float \n",s->nnodos);
	for (i=0;i<s->nnodos;i++){
		fprintf(fp,"%lf %lf %lf \n", nodo->x,nodo->y,nodo->z);
	        nodo++;
   	}

	// Imprimimos en el mismo fichero la informacion correspondiente a las celdas

    	fprintf(fp,"CELLS %ld %ld \n", s->nceldas,s->nceldas * 4);

    	for(i=0;i<s->nceldas;i++){
        fprintf(fp,"3 %ld %ld %ld \n",celda->n1,celda->n2,celda->n3);
	     celda++;
    	}

	fprintf(fp,"CELL_TYPES %ld \n",s->nceldas);

	for(i=0;i<s->nceldas;i++){
	   fprintf(fp,"5 \n");
	}


//Imprimimos informacion sobre las variables fisicas del fluido

      fprintf (fp, "CELL_DATA %ld \n", s->nceldas);

      fprintf (fp, "SCALARS z-nonEro float \n");
      fprintf (fp, "LOOKUP_TABLE default \n");
      celda=s->celda;
      for (i=0; i<s->nceldas; i++){
	      fprintf (fp,"%e \n", celda->z);
	      celda++;
      }

        fprintf (fp, "SCALARS h float \n");
        fprintf (fp, "LOOKUP_TABLE default \n");
        celda=s->celda;
        for (i=0; i<s->nceldas; i++){
	        fprintf (fp,"%e \n", celda->h);
	        celda++;
        }

        fprintf (fp, "SCALARS eSED float \n");
        fprintf (fp, "LOOKUP_TABLE default \n");
        celda=s->celda;
        for (i=0; i<s->nceldas; i++){
	        fprintf (fp,"%lf \n", celda->epsilon);
	        celda++;
        }


        fprintf (fp, "SCALARS bed-surf float \n");
        fprintf (fp, "LOOKUP_TABLE default \n");
        celda=s->celda;
        for (i=0; i<s->nceldas; i++){
	        fprintf (fp,"%lf \n", celda->z+celda->epsilon);
	        celda++;
        }



        fclose(fp);

}
///////////////////////////////////////

double triangle_size(tipo_nodo *n1, tipo_nodo *n2, tipo_nodo *n3){

	double x1, x2, x3, x4;

	x1 = n1->x - n2->x;
	x2 = n1->y - n2->y;
	x3 = n1->x - n3->x;
	x4 = n1->y - n3->y;

	return fabs(0.5*(x1*x4-x2*x3));

}

//////////////////////////////////////



///////////////////////////////////////
int main(void){

	int i;
	char fname[1024];
	tipo_sim sim;

	lector_triangle(&sim, "mesh_triangle.1");

	generate_peka_mesh(&sim, "case1.msh");
	generate_peka_ini(&sim, "case1.ini");
	generate_peka_solini(&sim, "case1.solini");
	
	//volcado_vtk(&sim,"case1.vtk");

	return 0;

}
