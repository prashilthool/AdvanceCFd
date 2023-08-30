#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define N 128

int main (void)
{
	double u[N][N+1], un[N][N+1], uc[N][N],RHSU[N][N+1];
	double v[N+1][N], vn[N+1][N], vc[N][N],RHSV[N+1][N] ;
	double p[N+1][N+1], pn[N+1][N+1], pc[N][N];
	double m[N+1][N+1];
	int i, j, step;
	double dx, dy, dt, delta, error, Re,beta;
	step =1;
	dx = 1.0/(N-1);
	dy = 1.0/(N-1);
	dt = 0.001;
	error = 1.0;
	Re = 400.0;
    beta = dx/dy;
	
	// Initializing u
		for (i=0; i<=(N-1); i++)
		{
			for (j=0; j<=(N); j++)
			{
				u[i][j] = 0.0;
				u[i][N] = 1.0;
				u[i][N-1] = 1.0;
                //printf("%f\t",v[i][j]);
			}
		}
		
	// Initializing v
		for (i=0; i<=(N); i++)
		{
			for (j=0; j<=(N-1); j++)
			{
				v[i][j] = 0.0;
                //printf("%f\t",v[i][j]);
			}
		}
		
	// Initializing p
		for (i=0; i<=(N); i++)
		{
			for (j=0; j<=(N); j++)
			{
				p[i][j] = 1.0;
			}
		}
	
	while (error > 1e-6)
	{
		// Solve RHSU
		for (i=1; i<=(N-2); i++)
		{
			for (j=1; j<=(N-1); j++)
			{
			  RHSU[i][j] = u[i][j] - (dt/dx)*0.25*( pow((u[i+1][j]+u[i][j]),2)-pow( (u[i-1][j]+u[i][j]),2))-
               (dt/dy)*0.25*((u[i][j+1]+u[i][j])*(v[i+1][j]+v[i][j])-(u[i][j-1]+u[i][j])*(v[i][j-1]+v[i+1][j-1]))
               + (dt/(Re*dx*dx))*(u[i-1][j]-2*u[i][j]+u[i+1][j])+ (dt/(Re*dy*dy))*(u[i][j-1]-2*u[i][j]+u[i][j+1]);
									 // printf("%f",RHSU[i][j]);
			}
		}
		

		
		
		// Solves RHSV
		for (i=1; i<=(N-1); i++)
		{
			for (j=1; j<=(N-2); j++)
			{
				RHSV[i][j] = v[i][j] - (dt/dx)*0.25*((u[i][j+1]+u[i][j])*(v[i+1][j]+v[i][j])-(u[i-1][j]+u[i-1][j+1])*(v[i-1][j]+v[i][j]))
                                      -(dt/dy)*0.25*( pow((v[i][j]+v[i][j+1]),2)-pow((v[i][j]+v[i][j-1]),2))
                                      +(dt/(Re*dx*dx))*(v[i-1][j]-2*v[i][j]+v[i+1][j])+(dt/(Re*dy*dy))*(v[i][j-1]-2*v[i][j]+v[i][j+1]);
									
			}
		}
		

		
	
		// Solves pressure 
		for (i=1; i<=(N-1); i++)
		{
			for (j=1; j<=(N-1); j++)
			{
				pn[i][j] = (1.0/(2.0*(1+beta)))*(pn[i+1][j]+pn[i-1][j]+pow(beta,2)*(pn[i][j+1]+pn[i][j-1])
                            -(pow(dx,2)/dt)*((RHSU[i][j]-RHSU[i-1][j])/dx +(RHSV[i][j]-RHSV[i][j-1])/dy));
				
				
				
			}
		}
		
		
		// Boundary conditions
		for (i=1; i<=(N-1); i++)
		{
			pn[i][0] = pn[i][1];
			pn[i][N] = pn[i][N-1];
		}
		
		for (j=0; j<=(N); j++)
		{
			pn[0][j] = pn[1][j];
			pn[N][j] = pn[N-1][j];
		}

        // Solves u-momentum
        for (i=1; i<=(N-2); i++)
		{
			for (j=1; j<=(N-1); j++)
			{
			 un[i][j] =  - (dt/dx)*(pn[i+1][j]-pn[i][j])+RHSU[i][j];
			}
		}
        // Boundary conditions
		for (j=1; j<=(N-1); j++)
		{
			un[0][j] = 0.0;
			un[N-1][j] = 0.0;
		}
		
		for (i=0; i<=(N-1); i++)
		{
			un[i][0] = -un[i][1];
			un[i][N] = 2 - un[i][N-1];
		}	
        // Solves v-momentum
		for (i=1; i<=(N-1); i++)
		{
			for (j=1; j<=(N-2); j++)
			{
				vn[i][j] = - (dt/dy)*(pn[i][j+1]-pn[i][j])+RHSV[i][j];
									
			}
		}

        // Boundary conditions
		for (j=1; j<=(N-2); j++)
		{
			vn[0][j] = -vn[1][j];
			vn[N][j] = -vn[N-1][j];
		}		

		for (i=0; i<=(N); i++)
		{
			vn[i][0] = 0.0;
			vn[i][N-1] = 0.0;
		}	
		
		// Displaying error
		error = 0.0;
		
		for (i=1; i<=(N-1); i++)
		{
			for (j=1; j<=(N-1); j++)
			{
				m[i][j] = (  ( un[i][j]-un[i-1][j] )/dx + ( vn[i][j]-vn[i][j-1] )/dy  );
				error = error + fabs(m[i][j]);
			}
		}

        //printf("Error is %9.8lf for the step %d\n", error, step);
		
		if (step%1000 ==1)
		{
	    printf("Error is %.9lf for the step %d\n", error, step);
		}
		
		
		// Iterating u
		for (i=0; i<=(N-1); i++)
		{
			for (j=0; j<=(N); j++)
			{
				u[i][j] = un[i][j];
			}
		}
		
		// Iterating v
		for (i=0; i<=(N); i++)
		{
			for (j=0; j<=(N-1); j++)
			{
				v[i][j] = vn[i][j];
			}
		}
		
		// Iterating p
		for (i=0; i<=(N); i++)
		{
			for (j=0; j<=(N); j++)
			{
				p[i][j] = pn[i][j];
			}
		}

		step = step + 1;
	
	}
	
	for (i=0; i<=(N-1); i++)
	{
		for (j=0; j<=(N-1); j++)
		{	
			uc[i][j] = 0.5*(u[i][j]+u[i][j+1]);
			vc[i][j] = 0.5*(v[i][j]+v[i+1][j]);
			pc[i][j] = 0.25*(p[i][j]+p[i+1][j]+p[i][j+1]+p[i+1][j+1]);
		}
	}
	
	
	
	// OUTPUT DATA
	FILE *file2, *file3, *file1;
	file2 = fopen("UVP.plt","w+t");
	file3 = fopen("Central_U.plt","w+t");
	file1 = fopen("Central_V.plt","w+t");

	if ( file2 == NULL )
	{
    printf("\nERROR when opening file\n");
    fclose( file2 );
	}

  else
	{
	fprintf( file2, "VARIABLES=\"X\",\"Y\",\"U\",\"V\",\"P\"\n");
	fprintf( file2, "ZONE  F=POINT\n");
	fprintf( file2, "I=%d, J=%d\n", N, N );

	for ( j = 0 ; j < (N) ; j++ )
	{
    for ( i = 0 ; i < (N) ; i++ )
    {
		double xpos, ypos;
		xpos = i*dx;
		ypos = j*dy;

		fprintf( file2, "%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\n", xpos, ypos, uc[i][j], vc[i][j], pc[i][j] );
    }
	}
	}

	fclose( file2 );
	
	// CENTRAL --U
  fprintf(file3, "VARIABLES=\"U\",\"Y\"\n");
  fprintf(file3, "ZONE F=POINT\n");
  fprintf(file3, "I=%d\n", N );

  for ( j = 0 ; j < N ; j++ )
  {
	double ypos;
    ypos = (double) j*dy;

    fprintf( file3, "%5.8lf\t%5.8lf\n", (uc[N/2][j] + uc[(N/2)+1][j])/(2.), ypos );
  }
  	// CENTRAL --V
  fprintf(file1, "VARIABLES=\"X\",\"V\"\n");
  fprintf(file1, "ZONE F=POINT\n");
  fprintf(file1, "J=%d\n", N);
  for ( i = 0 ; i < N ; i++ )
  {
	double xpos;
    xpos = (double) i*dx;
    fprintf( file1, "%5.8lf\t%5.8lf\n", xpos,(vc[(N/2)+1][i] + vc[(N/2)][i])/(2.) );
}
}