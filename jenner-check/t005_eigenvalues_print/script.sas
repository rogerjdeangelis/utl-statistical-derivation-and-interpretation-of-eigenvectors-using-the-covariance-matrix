/* EIGENVALUES of the 2x2 covariance matrix, sent to the listing.
   Same covariance-matrix literals (a,b,c,d) and the same quadratic-formula
   expressions as the repo's derivation block, computed in a DATA step and
   displayed with PROC PRINT so the two eigenvalues appear in the output
   listing (e1=1.4644714772, e2=0.1619589228 in the write-up). */

data eig;
  a=0.9351275;   b=0.6397438;
  c=0.6397438;   d=0.6913029;

  Eigenvalue1 = a/2 + d/2 + sqrt(a**2 - 2*a*d + 4*b**2 + d**2)/2;
  Eigenvalue2 = a/2 + d/2 - sqrt(a**2 - 2*a*d + 4*b**2 + d**2)/2;
  keep Eigenvalue1 Eigenvalue2;
run;

proc print data=eig noobs;
  var Eigenvalue1 Eigenvalue2;
run;
