#include "matrix.h"

Matrix::Matrix(int M,int N,string c_name,int c_type)
{
    rows=M;cols=N;
    name=c_name;
    type=c_type;

    mat=new double[rows*cols];
}

Matrix::~Matrix()
{
    delete []mat;
}

void Matrix::draw()
{
    if(type==as_matrix)
    {
        cout<<"Matrix "<<name<<": "<<endl;

        for(int i=0;i<rows;i=i+1)
        {
            for(int j=0;j<cols;j=j+1)
            {
                if(mat[i*cols+j]<10E-13&&mat[i*cols+j]>-10E-13)
                    mat[i*cols+j]=0;

                cout<<setprecision(3)<<setw(11)<<left<<mat[i*cols+j];
            }

            cout<<setw(0)<<endl;
        }
        cout<<endl;
    }

    else if(type==as_vector)
    {
        cout<<"Vector "<<name<<": "<<endl;

        for(int i=0;i<rows;i=i+1)
        {
            if(mat[i]<10E-13&&mat[i]>-10E-13)
                mat[i]=0;

            cout<<setprecision(3)<<setw(11)<<left<<mat[i];
        }
        cout<<setw(0)<<endl<<endl;
    }
}
