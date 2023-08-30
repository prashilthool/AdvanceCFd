# include<stdio.h>
# include<math.h>

# define nj 41
# define ni 41
int L(int i,int j)
{
    int l;
    l = (i-1) *nj + (j-1); 
    return l;
}
int main()
{
    double T[ni*nj],z[ni*nj],b[ni*nj],Tem[ni*nj],TA[ni][nj],Tana[ni][nj];
    double ap,as,aw,an,ae,dx,dy,beta,error,sum,sum1,pi=(22.0/7);
    double lw[ni*nj],ls[ni*nj],ue[ni*nj],un[ni*nj],lp[ni*nj],mnw[ni*nj],ztem[ni*nj],
            mse[ni*nj],p1[ni*nj],p2[ni*nj],B[ni*nj],r[ni][nj];
    double Lenght =1.0,width =1.0;

    int i,j,l,k,iter=0;
    dx = Lenght/(ni-1);
    dy = width/(nj-1);
    beta = dx/dy;
    ae = 1.0;
    aw = 1.0;
    as = beta * beta;
    an = beta * beta;
    ap = -2*(1+(beta*beta));

    //bounary conditions and U and L matrix
    for ( i = 1; i <= ni; i++)
    {
        for ( j = 1; j <=nj; j++)
        {
            l = L(i,j);
            lw[l] = aw;
            ls[l] = as;
            lp[l] = ap  - lw[l] * ue[l-nj]- ls[l] * un[l-1];
            un[l] = an/lp[l];
            ue[l] = ae/lp[l];
            z[l] =0.0;
            T[l] =0.0;
            b[l]=0.0;
            //printf("%f\t",lp[l]);

            if (i==1)
            {
                lp[l] = 1.0;
                b[l] = 0.0;
                lw[l] = 0.0;
                ls[l] = 0.0;
                ue[l] = 0.0;
                un[l] = 0.0;
                b[l] =0.0;

            }
            if (j== 1)
            {
                lp[l] = 1.0;
                b[l] = 0.0;
                lw[l] = 0.0;
                ls[l] = 0.0;
                ue[l] = 0.0;
                un[l] = 0.0;
                b[l] =0.0;
            }
            if (i==ni)
            {
                lp[l] = 1.0;
                b[l] = 0.0;
                lw[l] = 0.0;
                ls[l] = 0.0;
                ue[l] = 0.0;
                un[l] = 0.0; 
                b[l] =0.0;
            }
            if (j==nj)
            {
                lp[l] = 1.0;
                b[l] = 0.0;
                lw[l] = 0.0;
                ls[l] = 0.0;
                ue[l] = 0.0;
                un[l] = 0.0;
                b[l] = 1.0;
            } 

            //printf("%f\t",un[l]);  
        }
        
    }


    // creating the N matrix
    for ( i =1; i <ni; i++)
    {
        for ( j = 1; j < nj; j++)
        {
            l = L(i,j);
            mnw[l] = lw[l]*un[l-nj];
            mse[l] = ls[l]*ue[l-1];
            //printf("%f\t",mse[l]);
        }
        
    }

    
    error =1.0;
    
    while (error>1e-6)
    {
        //creating b prime matrix 
        for ( i = 1; i<= ni; i++)
        {
            for ( j = 1; j<= nj; j++)
            {
                l = L(i,j);
                Tem[l] = T[l];
            }
            
        }
        for ( i = 1; i<= ni; i++)
        {
            for ( j = 1; j<= nj; j++)
            {
                l = L(i,j);
                B[l] = b[l] + mnw[l] * T[l-nj+1] + mse[l] * T[l+nj-1];

            }
            
        }
        //z[0] = B[0]/lp[0];
        // forward substiturion
        for ( i = 1; i<= ni; i++)
        { 
            for ( j = 1; j<= nj; j++)
            {
                sum =0.0;
                l = L(i,j);
                z[l]  = (B[l]-(ls[l] * z[l-1] + lw[l] * z[l-nj]))/lp[l];
               

            }
            
        }
        // backward substitution
        //T[nj*ni]=z[nj*ni]; 
        for ( i = ni; i>= 1; i--)
        {
            
            for ( j =nj; j>=1; j--)
            {

                l = L(i,j);
                T[l] = z[l] -(un[l] * T[l+1] + ue[l] * T[l+nj]);
                

            }
            
            
        }
        error = 0.0;
        for ( i = 1; i<= ni; i++)
        {
            for ( j = 1; j<= nj; j++)
            {
                l = L(i,j);
                error = error + pow((T[l]-Tem[l]),2);

            }
            
        }
        error = sqrt(error);
        iter = iter + 1;
        printf("iter = %d\t error = %f\n",iter,error);
    }
    // numerical solution
    FILE*file;
    file = fopen("T.plt","w");
    fprintf(file," VARIABLES = \"X\",\"Y\",\"TA\"\n");
    fprintf(file," ZONE F = POINT\n");
    fprintf(file,"I=%d, J=%d\n",ni,nj);
        for ( i = 1; i<= ni; i++)
        {
            for ( j = 1; j<= nj; j++)
            {
                l = L(i,j);
                TA[i-1][j-1] = T[l];
                //TA[ni-1][j-1] =0.0;
                //TA[i-1][nj-1]=1.0;
                //TA[i-1][0] =0.0;
                //TA[0][j-1]=0.0;
                

                //printf("%f\t",TA[i-1][j-1]);
                fprintf(file,"%f\t%f\t%f\n",(i-1)*dx,(j-1)*dy,TA[i-1][j-1]);


            }
            
        }
    // analytical solution
    FILE*file1;
    file1 = fopen("Tanlytical.plt","w");
    fprintf(file1," VARIABLES = \"X\",\"Y\",\"Tana\"\n");
    fprintf(file1," ZONE F = POINT\n");
    fprintf(file1,"I=%d, J=%d\n",ni,nj);

    for ( i = 0; i < ni; i++)
    {
        for ( j = 0; j <nj; j++)
        {
          Tana[i][j] =0.0;
        }
        
    }
    for ( i = 0; i < ni; i++)
    {
        for ( j = 0; j <nj ; j++)
        {
            for ( k = 1; k <220; k++)
            {
                Tana[i][j] = Tana[i][j] + (2.0/pi) * ((pow(-1,k+1)+1)/k ) * sin(k*i*dx*pi)* (sinh(k*pi*j*dy)/sinh(k*pi));
            }
            
        }
        
    }
    for ( i = 0; i < ni; i++)
    {
        for ( j = 0; j <nj; j++)
        {
            printf("%f\t",Tana[i][j]);
            fprintf(file1,"%f\t%f\t%f\n",i*dx,j*dy,Tana[i][j]);
        }
        
    }
    

}
