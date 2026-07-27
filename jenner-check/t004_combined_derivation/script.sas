/* COMBINED DERIVATION: eigenvalues, raw eigenvectors and normalized
   eigenvectors in one pass. Verbatim data _null_ block from the repo
   write-up (section 2: "sas compute eigenvectors"). a,b,c,d are the
   covariance-matrix entries; the block prints all three result groups. */

 data _null_;;

    a=0.9351275;   b=0.6397438;
    c=0.6397438;   d=0.6913029;

    * eigenvalues;
    e1 = a/2 + d/2 - sqrt(a**2 - 2*a*d + 4*b**2 + d**2)/2;
    e2 = a/2 + d/2 + sqrt(a**2 - 2*a*d + 4*b**2 + d**2)/2;

    * non normalized eigenvectors;
    v11 =1;
    v21 =          -d/b + (a/2 + d/2 - sqrt(a**2 - 2*a*d + 4*b**2 + d**2)/2)/b;

    v12=1;
    v22 =          -d/b + (a/2 + d/2 + sqrt(a**2 - 2*a*d + 4*b**2 + d**2)/2)/b;

    * normalized eigervectors;
    norm11=v11   / sqrt(sum(1+v21**2));
    norm21=v21   / sqrt(sum(1+v21**2));
    norm12=v12   / sqrt(sum(1+v22**2));
    norm22=v22   / sqrt(sum(1+v22**2));

 put "Eigenvalues" / e1= / e2= //
     "Raw Eigenvectors" / v11= v21= / v12= v22= //
     "Normalized Eigenvectors" / norm11=  norm21= / norm12=  norm22= ;
 run;quit;
