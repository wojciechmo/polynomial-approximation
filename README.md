# Polynomial approximation
Discrete function polynomial approximation is determined based on solving overdetermined system Ax=b of linear equations with SVD matrix decomposition.
Pseudoinverse for overdetermined system makes it possible to find soution which minimizes mean squared error.
SVD matrix decomposition is done in order to determine pseudoinverse without calculating badly conditioned matrix.
Having SVD decomposition A=USVT (U,V - orthogonal matricies, S - diagonal matrix with matrix A singular values) 
it is extremely easy to calculate pseudoinverse: pinv(A)=VS^-1UT. Then searched soultion is defined as x=pinv(A)b.
Coefficients of system solution vector are polynomial coefficients. Matrix A is created based on discrete points x-coordinates, 
vector b, on discrete points y-coordinates. See examples below (accurate explanation is below images):

<img src="https://s23.postimg.org/7s1c7bzsb/pol1.png" width="500">
<img src="https://s23.postimg.org/onl60sqp7/pol2.png" width="500">
<img src="https://s28.postimg.org/qro2sm925/pol3.png" width="500">

For discrete function { (-2.72,3.85), (-0.95,0.63), (3.13,-2.40), (6.22,0.68), (7.22,4.08) } detetrmined quadratic function (second degree polynomial): y=-1.14-1.22x+0.26x^2. Matricies for quadratic approximation:

<img src="https://s29.postimg.org/nrdnodtp3/pol5.png" width="390">

Matricies for linear approximation of another discrete function: 

<img src="https://s29.postimg.org/jwoytqf47/pol4.png" width="300">



