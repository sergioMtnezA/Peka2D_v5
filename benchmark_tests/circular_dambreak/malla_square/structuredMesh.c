// STRUCTURED MESH GENERATOR (rectangular and triangular) for SFS2D
// Javier Fernández Pato
// Last modification: 13/11/12
// ----------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NVERT 4 		// Número de vértices de la celda

#define X_CELLS 500 	// Número de celdas en la dirección X
#define Y_CELLS 500		// Número de celdas en la dirección Y

#define dx 2
#define dy 2

#define xmin -500.0
#define ymin -500.0

#define Nsed 0

#define caso 0

	     /*  Caso 0: dry dambreak sobre flat bed*/

// ----------------------------------------------------------------

FILE *data; 
FILE *data2; 
FILE *data3; 
FILE *data4; 
FILE *data5;

int main()
{

	long a,b,c,d,i,j,NCELL,NNODE,X_NODES,Y_NODES;

	double x,y;
	double z,zr;
	double H,U,V;

	double *xc, *yc;
	double x1,x2,x3,y1,y2,y3;
	int count;
	double aux1;
	
	double slope, span, heigh;

	// Fracciones de sedimento en lecho
	double bed1, bed2, bed3;
	// Concentraciones de sedimento en suspension
	double susp1, susp2, susp3;


	data=fopen("case1.msh","w");
	data2=fopen("case1.ini","w");
	data3=fopen("case1.sedini","w");
	//data4=fopen("case1peka.txt","w");
	//data5=fopen("case1peka.inicial","w");
	

	printf("\n\nSTRUCTURED MESH GENERATOR\n");
	printf("-------------------------\n");
	if(NVERT==3)
	{
	   printf("\nMesh geometry: TRIANGULAR\n");
	   printf("\nNumber of cells: %d\n",2*X_CELLS*Y_CELLS);
	}
	if(NVERT==4)
	{
	   printf("\nMesh geometry: RECTANGULAR\n");
	   printf("\nNumber of cells: %d\n",X_CELLS*Y_CELLS);
	}

	/* Calculamos el número de celdas y nodos */
	X_NODES=X_CELLS+1; // Número de nodos en la dirección X
	Y_NODES=Y_CELLS+1; // Número de nodos en la dirección Y
	if(NVERT==3) NCELL=2*X_CELLS*Y_CELLS;
	if(NVERT==4) NCELL=X_CELLS*Y_CELLS;
	NNODE=X_NODES*Y_NODES;

	xc = (double*) malloc((NCELL)*sizeof(double));
	yc = (double*) malloc((NCELL)*sizeof(double));

	
	// Archivo malla.msh	
	fprintf(data,"NVERT          %d\n",NVERT);
	fprintf(data,"NCELL          %ld\n",NCELL);
	fprintf(data,"NNODE          %ld\n",NNODE);

	//fprintf(data4,"NVERT          %d\n",NVERT);
	//fprintf(data4,"NCELL          %ld\n",NCELL);
	//fprintf(data4,"NNODE          %ld\n",NNODE);

	/* Posición de los nodos */
	for(j=0;j<=Y_CELLS;j++) {
	   for(i=0;i<=X_CELLS;i++) {
		x = xmin+i*dx;
		y = ymin+j*dy;
	      fprintf(data,"%lf\t%lf\n",x,y);
	      //fprintf(data4,"%lf\t%lf\n",x,y);
	   }
	}



	//------------------ MALLA TRIANGULAR ------------------
	if(NVERT==3) {
	count=0;
	for(j=1;j<=Y_CELLS;j++) {
	   for(i=1;i<=X_CELLS;i++)
	   {
		a=i+(j-1)*X_NODES;
		b=a+1;
		c=b+X_NODES;
		d=a+X_NODES;

		if((i%2==0&&j%2==0) || (i%2!=0&&j%2!=0)) // X,Y misma paridad
		{
			/*	d---c
				| \ |
				a---b 
			*/
		   fprintf(data,"%ld\t%ld\t%ld\n",a,b,d);
		   //fprintf(data4,"%ld\t%ld\t%ld\n",a,b,d);

			x1 = xmin+(i-1)*dx;	y1 = ymin+(j-1)*dy;
			x2 = xmin+i*dx;		y2 = ymin+(j-1)*dy;
			x3 = xmin+(i-1)*dx;	y3 = ymin+j*dy;

			xc[count] = 1./3.*(x1+x2+x3);
			yc[count] = 1./3.*(y1+y2+y3);

			count++;
		
		   fprintf(data,"%ld\t%ld\t%ld\n",b,c,d);
		   //fprintf(data4,"%ld\t%ld\t%ld\n",b,c,d);

			x1 = xmin+i*dx;		y1 = ymin+(j-1)*dy;
			x2 = xmin+i*dx;		y2 = ymin+j*dy;
			x3 = xmin+(i-1)*dx;	y3 = ymin+j*dy;

			xc[count] = 1./3.*(x1+x2+x3);
			yc[count] = 1./3.*(y1+y2+y3);

			count++;
		}

		else // X,Y distinta paridad
		{
			/*	d---c
				| / |
				a---b 
			*/
		   fprintf(data,"%ld\t%ld\t%ld\n",a,b,c);
		   //fprintf(data4,"%ld\t%ld\t%ld\n",a,b,c);

			x1 = xmin+(i-1)*dx;	y1 = ymin+(j-1)*dy;
			x2 = xmin+i*dx;		y2 = ymin+(j-1)*dy;
			x3 = xmin+i*dx;		y3 = ymin+j*dy;

			xc[count] = 1./3.*(x1+x2+x3);
			yc[count] = 1./3.*(y1+y2+y3);

			count++;

		   fprintf(data,"%ld\t%ld\t%ld\n",c,d,a);
		   //fprintf(data4,"%ld\t%ld\t%ld\n",c,d,a);

			x1 = xmin+i*dx;		y1 = ymin+j*dy;
			x2 = xmin+(i-1)*dx;	y2 = ymin+j*dy;
			x3 = xmin+(i-1)*dx;	y3 = ymin+(j-1)*dy;

			xc[count] = 1./3.*(x1+x2+x3);
			yc[count] = 1./3.*(y1+y2+y3);

			count++;
		}


	   }
	
	}
	}
	//-------------------------------------------------------



	//------------------ MALLA RECTANGULAR ------------------
	if(NVERT==4) {
	for(j=1;j<=Y_CELLS;j++) {
	   for(i=1;i<=X_CELLS;i++)
	   {
		a=i+(j-1)*X_NODES;
		b=a+1;
		c=b+X_NODES;
		d=a+X_NODES;

		/* Conectividades celdas-nodos: d---c
							  |   |
							  a---b */

		fprintf(data,"%ld\t%ld\t%ld\t%ld\n",a,b,c,d);
		//fprintf(data4,"%ld\t%ld\t%ld\t%ld\n",a,b,c,d);

	   }
	}
	}
	//-------------------------------------------------------
	





	/* Altura de las celdas */


	/*  Caso 0: dry dambreak sobre flat bed*/
	#if caso == 0

	   H=0.0; // calado inicial
	   U=0.0; // Vx inicial
	   V=0.0; // Vy inicial

	   // Archivo malla.ini
	   fprintf(data2,"NIVEL_SUPERFICIAL_INICIAL   0\n");
	   fprintf(data2,"VELOCIDADX_INICIAL   0\n");
	   fprintf(data2,"VELOCIDADY_INICIAL   0\n");
	   fprintf(data2,"LEER_CELDAS          4\n");
	   //fprintf(data2,"LEER_CELDAS          3\n");

	   //fprintf(data5,"NIVEL_SUPERFICIAL_INICIAL	0\n");
	   //fprintf(data5,"VELOCIDADX_INICIAL	0\n");
	   //fprintf(data5,"VELOCIDADY_INICIAL	0\n");
	   //fprintf(data5,"SUSP_SED_CONCENTRATION	0\n");
	   //fprintf(data5,"LEER_CELDAS			4\n");


	   // Archivo malla.sedini
	   fprintf(data3,"NUMERO_SEDIMENTOS   %d\n",Nsed);
	   fprintf(data3,"LEER_CELDAS         %d\n",2*Nsed);

	   // Cota no erosionable
	   zr = 0.0;

		//Cota erosionable
		slope=0.05;
		span=sqrt(2.*500.*500.)-50.;
		heigh=slope*span;


	   if(NVERT == 4) {

		   for(j=0;j<=Y_CELLS-1;j++) {
			for(i=0;i<=X_CELLS-1;i++)
			{

				// Centro de la celda (malla cuadrada)
				x=xmin+i*dx+0.5*dx;
				y=ymin+j*dy+0.5*dy;

				// Superficie de
			  	aux1 = sqrt(x*x+y*y);

				//if( (x>-10 && x<=10) &&
              //    (y>-10 && y<=10)) {
              if(aux1<=50.){        
					z = heigh;
					H = 10.0;
					U = 0.0;
					V = 0.0;
					susp1 = 0.00;
					//susp2 = 0.0;  

				} else {
					z = heigh-slope*(aux1-50.);
					//H = H = 0.1 - z;
					H = 0.0;
					U = 0.0;
					V = 0.0;
					susp1 = 0.0;
					//susp2 = 0.0; 
				}
				// Fracciones de sedimento separadas en la diagonal
				bed1 = 1.0;
				//bed2 = 0.5;


				fprintf(data,"%lf\n",zr); // malla.msh
				fprintf(data2,"%lf\t%lf\t%lf\t%lf\n",H,U,V,z); // malla.ini
				//printf(data2," %lf   %lf   %lf\n",H,U,V); // malla.ini
				//fprintf(data3,"%lf\t%lf\t%lf\t%lf\n",bed1,bed2,susp1,susp2); // malla.sedini
				fprintf(data3,"%lf\t%lf\n",bed1,susp1); // malla.sedini


				//fprintf(data4,"%lf\n",zr); // malla.txt
				//fprintf(data5,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",H,U,V,z,bed1,bed2); // malla.inicial
				//fprintf(data5,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",H,U,V,z,bed1,susp1); // malla.inicial+solutos

			   }

			}

		}




	   if(NVERT == 3) {

		   count = 0;
		   for(i=0;i<=NCELL-1;i++) {

				// Centro de la celda
				x=xc[count];
				y=yc[count];

				// Superficie de
			  	aux1 = sqrt(x*x+y*y);

				//if( (x>-10 && x<=10) &&
                //    (y>-10 && y<=10)) {
                if(aux1<=50.){
					z = heigh;
					H = 10.0;
					U = 0.0;
					V = 0.0;
					susp1 = 0.00;
					//susp2 = 0.0;  

				} else {
					z = heigh-slope*(aux1-50.0);
					//H = 0.1 - z;
					H = 0.0;
					U = 0.0;
					V = 0.0;
					susp1 = 0.0;
					//susp2 = 0.0; 
				}
				// Fracciones de sedimento separadas en la diagonal
				bed1 = 1.0;
				//bed2 = 0.5;



				fprintf(data,"%lf\n",zr); // malla.msh
				fprintf(data2,"%lf\t%lf\t%lf\t%lf\n",H,U,V,z); // malla.ini
				//printf(data2," %lf   %lf   %lf\n",H,U,V); // malla.ini
				//fprintf(data3,"%lf\t%lf\t%lf\t%lf\n",bed1,bed2,susp1,susp2); // malla.sedini
				fprintf(data3,"%lf\t%lf%lf\t%lf\n",bed1,bed2,susp1,susp2); // malla.sedini


				//fprintf(data4,"%lf\n",zr); // malla.txt
				//fprintf(data5,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",H,U,V,z,bed1,bed2); // malla.inicial
				//fprintf(data5,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",H,U,V,z,bed1,susp1); // malla.inicial+solutos


				count++;

			   }

		}

	#endif
	// ------------------------------------------------------------------







	
	printf("\nMesh generated succesfully! :)\n\n");

	free(xc); free(yc);

	fclose(data); 
	fclose(data2); 
	fclose(data3);
	//fclose(data4); 
	//fclose(data5);


}





