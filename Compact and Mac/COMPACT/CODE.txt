#include<stdio.h>
#include<stdlib.h>
#include<math.h>
# define n 91
int main(void)
{
    double phi[n][n],df_anl[n][n],df[n][n],tem[n][n];
    double dx , dy,x,y,alpha, beta,a,b,c,cons,error;
    int i, j,order,iteration=1; 
    dx = 10.0/(n-1);
    dy = 10.0/(n-1);
    cons = (1.0/(12.0*dx));
   
    

    FILE *file2; 
    file2 = fopen("phi.plt","w");
    fprintf(file2,"VARIABLES= \"X\",\"Y\",\"PHI\"\n");
    fprintf(file2,"ZONE F=POINT\n ");
    fprintf(file2,"I=%d, J=%d\n",n,n);
     // value of phi at every point 
    for ( i = 0; i <n; i++)
    {
        for ( j = 0; j < n; j++)
        {   
            x= dx*i;
            y=dx*j;
            phi[i][j]= sin(x)*cos(y);
            fprintf(file2,"%f\t%f\t%f\n",j*dx,i*dx,phi[i][j]);
            //printf("%f\t",phi[i][j]);
        }
        
    }
    FILE*file1; 
    file1 = fopen("analytical.plt","w");
     fprintf(file1,"VARIABLES= \"X\",\"Y\",\"DF_ANA\"\n");
    fprintf(file1,"ZONE F=POINT\n ");
    fprintf(file1,"I=%d, J=%d\n",n,n);
    // analytical solution of problem
    for ( i = 0; i < n; i++)
    {
        for (  j= 0; j <n; j++)
        {
            x = dx*i;
            y = dx*j;
            df_anl [i][j] = cos(x)*cos(y);
            //printf("%f\t",df_anl[i][j]);
            fprintf(file1,"%f\t%f\t%f\n",j*dx,i*dx,df_anl[i][j]);
            df[i][j]=0.0;
           

        }
        
    } 
    printf("choose the order of acurracy\n");
    scanf("%d",& order);
    // if and else  conditions for the order of the accuracy
     if (order == 4)
     {
        alpha = 0.0;
        beta=0.0;
        a = (2.0/3.0)*(alpha+2);
        b= (1.0/3.0)*(4*alpha-1);
        c = 0.0;
        printf(" 4th order \n");
        
     }
     else if (order==6)
     {
      alpha = 1.0/3.0;
      beta = 0;
      a = (1.0/6.0)*(alpha+9);
      b= (1.0/15.0)*(32*alpha-9);
      c = (1.0/10.0)*(-3*alpha+1);
      printf(" 6th order \n");
     }
     else if (order == 8)
     {
      alpha = 3.0/8.0;
      beta = (1.0/20.0)*(-3+8*alpha);
      a = (1.0/6.0)*(12-7*alpha);
      b= (1.0/150.0)*(568*alpha-183);
      c = (1.0/50.0)*(9*alpha-4);
      printf(" 8th order \n");
     }
     else if (order == 10)
     {
      alpha = 1.0/2.0;
      beta = (1.0/20.0)*(-3+8*alpha);
      a = (1.0/6.0)*(12-7*alpha);
      b= (1.0/150.0)*(568*alpha-183);
      c = (1.0/50.0)*(9*alpha-4);
      printf(" 10th order \n");

      
     }


    for ( i = 0; i <n ; i++)
     {
        // bottom boundary
        df[i][0] = cons*(-25*phi[i][0]+48*phi[i][1]-36*phi[i][2]+16*phi[i][3]-3*phi[i][4]);
        //adjecent points 
        df[i][1] = cons*(-25*phi[i][1]+48*phi[i][2]-36*phi[i][3]+16*phi[i][4]-3*phi[i][5]);
        
     }

     
    
   //  printf("\n");
     for ( i = 0; i <n ; i++)
     {
        // top boundary
        df[i][n-1] = cons*(25*phi[i][n-1]-48*phi[i][n-2]+36*phi[i][n-3]-16*phi[i][n-4]+3*phi[i][n-5]);
        // adjecent points
        df[i][n-2] = cons*(25*phi[i][n-2]-48*phi[i][n-3]+36*phi[i][n-4]-16*phi[i][n-5]+3*phi[i][n-6]);
        

     }

     
     for ( j = 0; j < n; j++)
     {
        //left boundary
        df[0][j] = cons*(-25*phi[0][j]+48*phi[1][j]-36*phi[2][j]+16*phi[3][j]-3*phi[4][j]);
        //adjecent points 
        df[1][j] = cons*(-25*phi[1][j]+48*phi[2][j]-36*phi[3][j]+16*phi[4][j]-3*phi[5][j]);
        
     }
     

     for ( j =0; j <n; j++)
     {
            // right boundary
            df[n-1][j] = cons*(25*phi[n-1][j]-48*phi[n-2][j]+36*phi[n-3][j]-16*phi[n-4][j]+3*phi[n-5][j]);
            // adjecent points
            df[n-2][j] = cons*(25*phi[n-2][j]-48*phi[n-3][j]+36*phi[n-4][j]-16*phi[n-5][j]+3*phi[n-6][j]);
            
     } 
    
    // printing .m file for visualization

    FILE*file6;
    file6 = fopen("df_bound.m","w");
    fprintf(file6,"df_bound=[");
     // printing in plt file
    for ( i = 0; i <n; i++)
    {
       for ( j = 0; j <n; j++)
       {

        //printf("%f\t",df[i][j]);
        fprintf(file6,"%f\t",df[i][j]);
        
        
       }
       //printf("\n");
       fprintf(file6,"\n");
    }
   

  
 
   
    FILE*file3;
    file3 = fopen("error.txt","w");
    error =1.0;
    while (error>1e-6 )
    {
        // assining the df values to tem
        for ( i = 0; i <n; i++)
        {
            for ( j = 0; j < n; j++)
            {
                tem[i][j]= df[i][j];
            }
            
        }
        // point gausss seidel
        for ( i = 2; i < n-2; i++)
        {
            for ( j = 2; j < n-2; j++)
            {
                 df[i][j] =  a*0.5*(phi[i+1][j]-phi[i-1][j])/dx + b*0.25*(phi[i+2][j]-phi[i-2][j])/dx + c*(phi[i+1][j]-phi[i-1][j])/(6.0*dx)
                            -beta*df[i-2][j]-alpha*df[i-1][j]-alpha*df[i+1][j]-beta*df[i+2][j];
                             
            }
            
        }
        // error
        error=0.0;
        for ( i = 2; i < n-2; i++)
        {
            for ( j = 2; j < n-2; j++)
            {
                error = error+pow((df[i][j]-tem[i][j]),2);
            }
            
        }
        error = sqrt(error/(n*n));
        printf("iteration %d\t",iteration);
        printf("error %.9f\n",error);
        fprintf(file3,"%d\t%.9f\n",iteration,error);
        iteration=iteration+1;
        
        
    }
    FILE*file4;
    file4 = fopen("df.plt","w");
    
    fprintf(file4,"VARIABLES = \"X\",\"Y\",\"DF\"\n");
    fprintf(file4,"ZONE F=POINT\n");
    fprintf(file4,"I=%d, J=%d\n",n,n);

     // printing in plt file
    for ( i = 0; i <n; i++)
    {
       for ( j = 0; j <n; j++)
       {

        printf("%f\t",df[i][j]);
        fprintf(file4,"%f\t%f\t%f\n",j*dx,i*dx,df[i][j]);
        
        
       }
       
    }
    // line graph variation y at x= 5
     FILE*file7;
     file7 = fopen("df_variation_y.plt","w");
     fprintf(file7,"VARIABLES = \"DF\",\"Y\"\n");
     fprintf(file7, "ZONE F =POINT\n");
     fprintf(file7,"I=%d\n",n);
     for ( j = 0; j < n; j++)
     {
        double ypos;
        ypos = j*dx;
        fprintf(file7,"%f\t%f\n",df[(n-1)/2][j],ypos);
     }
    // line graph variation x at y= 5
    FILE*file8;
     file8 = fopen("df_variation_x.plt","w");
     fprintf(file8,"VARIABLES = \"X\",\"DF\"\n");
     fprintf(file8, "ZONE F =POINT\n");
     fprintf(file8,"J=%d\n",n);
     for ( i = 0; i < n; i++)
     {
         double xpos;
        xpos = i*dx;
        fprintf(file8,"%f\t%f\n",xpos,df[i][(n-1)/2]);
     }
     // line graph analytical variation y at x= 5
      FILE*file5;
      file5 = fopen("df_analytical_variation_y.plt","w");
      fprintf(file5,"VARIABLES = \"DF_anl\",\"Y\"\n");
      fprintf(file5, "ZONE F =POINT\n");
      fprintf(file5,"I=%d\n",n);
    
    for ( j = 0; j < n; j++)
    {
        double ypos;
        ypos = j*dx;
        fprintf(file5,"%f\t%f\n",df_anl[(n-1)/2][j],ypos);
    }
    // line graph  anlytical variation x at y= 5
     FILE*file9;
     file9 = fopen("df_analytical_variation_x.plt","w");
     fprintf(file9,"VARIABLES = \"X\",\"DF_anl\"\n");
     fprintf(file9, "ZONE F =POINT\n");
     fprintf(file9,"J=%d\n",n);
     for ( i = 0; i < n; i++)
     {
        double xpos;
        xpos = i*dx;
        fprintf(file9,"%f\t%f\n",xpos,df_anl[i][(n-1)/2]);
     }

     fclose(file9);
     // point A and B for numerical and anlytical and errror calculations
     file9 = fopen("POINT_A_AND_POINT_B.txt","w");
     fprintf(file9,"numA\tnumB\tanlA\tanlB\tdiffA\tdiffB\n");
     fprintf(file9,"%.9f\t%.9f\t%.9f\t%.9f\t%.9f\t%.9f",df[27][45],df[54][27],df_anl[27][45],df_anl[54][27],(df[27][45]-df_anl[27][45])*100,(df[54][27]-df_anl[54][27])*100);
        
   
     

   
}
