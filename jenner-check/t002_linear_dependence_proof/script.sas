/* COMPUTATIONAL PROOF that the two eigenvector equations are linearly
   dependent for eigenvalue e1. Verbatim data _null_ block from the repo
   write-up (section 2: linear dependence). F is the proportionality
   constant; the two coefficient pairs are shown to be equal. */

data _null_;;

  a=0.9351275;   b=0.6397438;
  c=0.6397438;   d=0.6913029;

  e1=a/2 + d/2 + sqrt(a**2 - 2*a*d + 4*b**2 + d**2)/2;
  e2=a/2 + d/2 - sqrt(a**2 - 2*a*d + 4*b**2 + d**2)/2;

  F = b/(a-e1);

  /*
  F = b/(a-e1);

  equ1 = F*(a-e1)*v11  + b*v21)
  expanding

  equ1 = F*(a-e1)*v11 + F*b*v21
  equ2 =         bv11 + (d-e1)v21

  Need to show

  Note for equ1 to equal equ2 we need
      coef of v11 in equ1 = coef of v11 in equ2
      and
      coef of v21 in equ1 = coef of v12 in equ2

   F*(a-e1) = F*b
   F*b      = d-e

  */

   coef_v11_equ1 = F*(a-e1);
   coef_v11_equ2 = b;

   coef_v21_equ1 = F*b;
   coef_v21_equ2 = d-e1;


   put
       coef_v11_equ1= ' = '
       coef_v11_equ2=  /
       coef_v21_equ1= ' = '
       coef_v21_equ2=
   ;

   put
       coef_v11_equ1
       coef_v21_equ1 /
       coef_v11_equ2
       coef_v21_equ2
   ;
run;quit;
