#include "svd.h"

double* SVD::pseudo_inverse_SVD(double* V,double *S,double *U,int rows,int cols)
{
    double* temp=new double[rows*rows];

    for(int i=0;i<rows;i++)
        for(int j=0;j<rows;j++)
           temp[i*rows+j]=U[j*rows+i]/S[j];

    double* pinv=new double[rows*cols];

    float sum;
    for(int i=0;i<rows;i++)
        for(int j=0;j<cols;j++)
        {
           sum=0;
           for(int k=0;k<rows;k++)
               sum=sum+temp[i*rows+k]*V[j*rows+k];

           pinv[i*cols+j]=sum;
        }

    delete temp;

    return pinv;
}

int SVD::multiply (Matrix *m1, Matrix *m2, Matrix *y)
{
    if(m1->cols!=m2->rows)
        return 0;

    double sum;
    for(int i=0;i<m1->rows;i=i+1)
    {
        for(int j=0;j<m2->cols;j=j+1)
        {
            sum=0;
            for(int k=0;k<m1->cols;k=k+1)
                sum=sum+m1->mat[i*m1->cols+k]*m2->mat[k*m2->cols+j];

            y->mat[i*y->cols+j]=sum;
        }
    }

    return 1;
}

Matrix* SVD::pinv_compute(Matrix *A,int rows,int cols)
{
    A_pinv=new Matrix(cols,rows,"A_pinv",Matrix::as_matrix);
    U=new Matrix(rows,cols,"U",Matrix::as_matrix);
    VT=new Matrix(cols,cols,"VT",Matrix::as_matrix);
    S=new Matrix(cols,cols,"S",Matrix::as_matrix);

    double* M=new double[rows*cols];
    for(int i=0;i<rows*cols;i++)
        M[i]=A->mat[i];

    int m = cols;
    int n = rows;
    double* tU =new double[m*m];
    double* tS =new double[m];
    double* tVT=new double[m*n];

    svd_compute(&m, &n, &M[0], &m, &tS[0], &tU[0], &m, &tVT[0], &m);

    for(int i=0;i<rows;i++)
        for(int j=0;j<cols;j++)
            U->mat[i*cols+j]=tVT[i*cols+j];

    for(int i=0;i<cols;i++)
        for(int j=0;j<cols;j++)
            VT->mat[i*cols+j]=tU[i*cols+j];

    for(int i=0;i<cols;i++)
        for(int j=0;j<cols;j++)
        {
            if(i==j) S->mat[i*cols+j]=tS[i];
            else S->mat[i*cols+j]=0;
        }

    double* pinv;
    pinv=pseudo_inverse_SVD(tVT,tS,tU,m,n);
    A_pinv->mat=pinv;

    delete[] tU;
    delete[] tS;
    delete[] tVT;
    delete[] M;

    return A_pinv;
}

#define U(i,j) U_[(i)*dim[0]+(j)]
#define S(i,j) S_[(i)*dim[1]+(j)]
#define V(i,j) V_[(i)*dim[1]+(j)]

void SVD::GivensL(double* S_, const size_t dim[2], size_t m,double a, double b)
{
  double r=sqrt(a*a+b*b);
  double c=a/r;
  double s=-b/r;

  for(size_t i=0;i<dim[1];i++)
  {
    double S0=S(m+0,i);
    double S1=S(m+1,i);
    S(m,i)+=S0*(c-1);
    S(m,i)+=S1*(-s);
    S(m+1,i)+=S0*(s);
    S(m+1,i)+=S1*(c-1);
  }
}

void SVD::GivensR(double* S_, const size_t dim[2], size_t m, double a, double b)
{
  double r=sqrt(a*a+b*b);
  double c=a/r;
  double s=-b/r;

  for(size_t i=0;i<dim[0];i++)
  {
    double S0=S(i,m+0);
    double S1=S(i,m+1);
    S(i,m)+=S0*(c-1);
    S(i,m)+=S1*(-s);
    S(i,m+1)+=S0*(s);
    S(i,m+1)+=S1*(c-1);
  }
}

void SVD::svd_decomp(const size_t dim[2], double* U_, double* S_, double* V_, double eps)
{
  assert(dim[0]>=dim[1]);
  size_t n=std::min(dim[0],dim[1]);
  std::vector<double> house_vec(std::max(dim[0],dim[1]));

  for(size_t i=0;i<n;i++)
  {
      double x1=S(i,i);
      if(x1<0) x1=-x1;

      double x_inv_norm=0;
      for(size_t j=i;j<dim[0];j++)
      {
        x_inv_norm+=S(j,i)*S(j,i);
      }
      if(x_inv_norm>0) x_inv_norm=1/sqrt(x_inv_norm);

      double alpha=sqrt(1+x1*x_inv_norm);
      double beta=x_inv_norm/alpha;

      house_vec[i]=-alpha;
      for(size_t j=i+1;j<dim[0];j++)
      {
        house_vec[j]=-beta*S(j,i);
      }
      if(S(i,i)<0) for(size_t j=i+1;j<dim[0];j++)
      {
        house_vec[j]=-house_vec[j];
      }

      for(size_t k=i;k<dim[1];k++)
      {
        double dot_prod=0;
        for(size_t j=i;j<dim[0];j++)
        {
          dot_prod+=S(j,k)*house_vec[j];
        }
        for(size_t j=i;j<dim[0];j++)
        {
          S(j,k)-=dot_prod*house_vec[j];
        }
      }

      for(size_t k=0;k<dim[0];k++)
      {
        double dot_prod=0;
        for(size_t j=i;j<dim[0];j++)
        {
          dot_prod+=U(k,j)*house_vec[j];
        }
        for(size_t j=i;j<dim[0];j++)
        {
          U(k,j)-=dot_prod*house_vec[j];
        }
      }

      if(i>=n-1) continue;
      {
        double x1=S(i,i+1);
        if(x1<0) x1=-x1;

        double x_inv_norm=0;
        for(size_t j=i+1;j<dim[1];j++)
        {
          x_inv_norm+=S(i,j)*S(i,j);
        }
        if(x_inv_norm>0) x_inv_norm=1/sqrt(x_inv_norm);

        double alpha=sqrt(1+x1*x_inv_norm);
        double beta=x_inv_norm/alpha;

        house_vec[i+1]=-alpha;
        for(size_t j=i+2;j<dim[1];j++)
        {
          house_vec[j]=-beta*S(i,j);
        }
        if(S(i,i+1)<0) for(size_t j=i+2;j<dim[1];j++)
        {
          house_vec[j]=-house_vec[j];
        }
      }

      for(size_t k=i;k<dim[0];k++)
      {
        double dot_prod=0;
        for(size_t j=i+1;j<dim[1];j++)
        {
          dot_prod+=S(k,j)*house_vec[j];
        }
        for(size_t j=i+1;j<dim[1];j++)
        {
          S(k,j)-=dot_prod*house_vec[j];
        }
      }

      for(size_t k=0;k<dim[1];k++)
      {
        double dot_prod=0;
        for(size_t j=i+1;j<dim[1];j++)
        {
          dot_prod+=V(j,k)*house_vec[j];
        }
        for(size_t j=i+1;j<dim[1];j++)
        {
          V(j,k)-=dot_prod*house_vec[j];
        }
      }
  }

  size_t k0=0;
  if(eps<0)
  {
    eps=1.0;
    while(eps+(double)1.0>1.0) eps*=0.5;
    eps*=64.0;
  }

  while(k0<dim[1]-1)
  {
    double S_max=0.0;
    for(size_t i=0;i<dim[1];i++) S_max=(S_max>S(i,i)?S_max:S(i,i));

    while(k0<dim[1]-1 && fabs(S(k0,k0+1))<=eps*S_max) k0++;
    if(k0==dim[1]-1) continue;

    size_t n=k0+2;
    while(n<dim[1] && fabs(S(n-1,n))>eps*S_max) n++;

    double alpha=0;
    double beta=0;
    {
      double C[2][2];
      C[0][0]=S(n-2,n-2)*S(n-2,n-2);
      if(n-k0>2) C[0][0]+=S(n-3,n-2)*S(n-3,n-2);
      C[0][1]=S(n-2,n-2)*S(n-2,n-1);
      C[1][0]=S(n-2,n-2)*S(n-2,n-1);
      C[1][1]=S(n-1,n-1)*S(n-1,n-1)+S(n-2,n-1)*S(n-2,n-1);

      double b=-(C[0][0]+C[1][1])/2;
      double c=  C[0][0]*C[1][1] - C[0][1]*C[1][0];
      double d=0;
      if(b*b-c>0) d=sqrt(b*b-c);
      else
      {
        double b=(C[0][0]-C[1][1])/2;
        double c=-C[0][1]*C[1][0];
        if(b*b-c>0) d=sqrt(b*b-c);
      }

      double lambda1=-b+d;
      double lambda2=-b-d;

      double d1=lambda1-C[1][1]; d1=(d1<0?-d1:d1);
      double d2=lambda2-C[1][1]; d2=(d2<0?-d2:d2);
      double mu=(d1<d2?lambda1:lambda2);

      alpha=S(k0,k0)*S(k0,k0)-mu;
      beta=S(k0,k0)*S(k0,k0+1);
    }

    for(size_t k=k0;k<n-1;k++)
    {
      size_t dimU[2]={dim[0],dim[0]};
      size_t dimV[2]={dim[1],dim[1]};
      GivensR(S_,dim ,k,alpha,beta);
      GivensL(V_,dimV,k,alpha,beta);

      alpha=S(k,k);
      beta=S(k+1,k);
      GivensL(S_,dim ,k,alpha,beta);
      GivensR(U_,dimU,k,alpha,beta);

      alpha=S(k,k+1);
      beta=S(k,k+2);
    }

    {
      for(size_t i0=k0;i0<n-1;i0++)
      {
        for(size_t i1=0;i1<dim[1];i1++)
        {
          if(i0>i1 || i0+1<i1) S(i0,i1)=0;
        }
      }
      for(size_t i0=0;i0<dim[0];i0++)
      {
        for(size_t i1=k0;i1<n-1;i1++)
        {
          if(i0>i1 || i0+1<i1) S(i0,i1)=0;
        }
      }
      for(size_t i=0;i<dim[1]-1;i++)
      {
        if(fabs(S(i,i+1))<=eps*S_max)
        {
          S(i,i+1)=0;
        }
      }
    }
  }
}

void SVD::svd_compute(int *M, int *N, double *A, int *LDA,double *S, double *U, int *LDU, double *VT, int *LDVT){

  const size_t dim[2]={std::max(*N,*M), std::min(*N,*M)};
  double* U_=new double[dim[0]*dim[0]]; memset(U_, 0, dim[0]*dim[0]*sizeof(double));
  double* V_=new double[dim[1]*dim[1]]; memset(V_, 0, dim[1]*dim[1]*sizeof(double));
  double* S_=new double[dim[0]*dim[1]];

  const size_t lda=*LDA;
  const size_t ldu=*LDU;
  const size_t ldv=*LDVT;

  if(dim[1]==*M)
  {
    for(size_t i=0;i<dim[0];i++)
    for(size_t j=0;j<dim[1];j++)
    {
      S_[i*dim[1]+j]=A[i*lda+j];
    }
  }
  else
  {
    for(size_t i=0;i<dim[0];i++)
    for(size_t j=0;j<dim[1];j++)
    {
      S_[i*dim[1]+j]=A[j*lda+i];
    }
  }
  for(size_t i=0;i<dim[0];i++)
  {
    U_[i*dim[0]+i]=1;
  }
  for(size_t i=0;i<dim[1];i++)
  {
    V_[i*dim[1]+i]=1;
  }

  svd_decomp(dim, U_, S_, V_, -1.0);

  for(size_t i=0;i<dim[1];i++)
  {
    S[i]=S_[i*dim[1]+i];
  }
  if(dim[1]==*M)
  {
    for(size_t i=0;i<dim[1];i++)
    for(size_t j=0;j<*M;j++)
    {
      U[j+ldu*i]=V_[j+i*dim[1]]*(S[i]<0.0?-1.0:1.0);
    }
  }
  else
  {
    for(size_t i=0;i<dim[1];i++)
    for(size_t j=0;j<*M;j++)
    {
      U[j+ldu*i]=U_[i+j*dim[0]]*(S[i]<0.0?-1.0:1.0);
    }
  }
  if(dim[0]==*N)
  {
    for(size_t i=0;i<*N;i++)
    for(size_t j=0;j<dim[1];j++)
    {
      VT[j+ldv*i]=U_[j+i*dim[0]];
    }
  }
  else
  {
    for(size_t i=0;i<*N;i++)
    for(size_t j=0;j<dim[1];j++)
    {
      VT[j+ldv*i]=V_[i+j*dim[1]];
    }
  }
  for(size_t i=0;i<dim[1];i++)
  {
    S[i]=S[i]*(S[i]<0.0?-1.0:1.0);
  }

  delete[] U_;
  delete[] S_;
  delete[] V_;
}

SVD::SVD(){ }

#undef U
#undef S
#undef V
