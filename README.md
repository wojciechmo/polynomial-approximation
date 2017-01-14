# Polynomial approximation
Discrete function polynomial approximation is determined based on solving overdetermined system Ax=b of linear equations with SVD matrix decomposition.
Pseudoinverse for overdetermined system makes it possible to find soution which minimizes mean squared error.
SVD matrix decomposition is done in order to determine pseudoinverse without calculating badly conditioned matrix.
Having SVD decomposition A=USVT (U,V - orthogonal matricies, S - diagonal matrix with matrix A singular values) 
it is extremely easy to calculate pseudoinverse: pinv(A)=VS^-1UT. Then searched soultion is defined as x=pinv(A)b.
Coefficients of system solution vector are polynomial coefficients. Matrix A is created based on discrete points x-coordinates, 
vector b, on discrete points y-coordinates.
