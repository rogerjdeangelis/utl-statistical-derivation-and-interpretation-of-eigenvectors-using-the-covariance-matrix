%let pgm=utl-statistical-derivation-and-interpretation-of-eigenvectors-using-the-covariance-matrix;

Statistical-derivation-and-interpretation-of-eigenvectors-using-the-covariance-matrix

Eigenvector/Priciple Components Plot
https://tinyurl.com/bdfv4cu6
https://github.com/rogerjdeangelis/utl-statistical-derivation-and-interpretation-of-eigenvectors-using-the-covariance-matrix/blob/main/elp.pdf

github
https://tinyurl.com/379vv3yr
https://github.com/rogerjdeangelis/utl-statistical-derivation-and-interpretation-of-eigenvectors-using-the-covariance-matrix

I hope this is all correct?

     Sections

           1. theory
           2. sas conpute eigenvectors for any 2x2 covariance matrix
           3  r compute plot eigenvectors

/*     _     _
| |__ (_)___| |_ ___  _ __ _   _
| `_ \| / __| __/ _ \| `__| | | |
| | | | \__ \ || (_) | |  | |_| |
|_| |_|_|___/\__\___/|_|   \__, |
                           |___/
*/

DIFFERENTIAL EQUATIONS

The history of eigenvalues and eigenvectors dates back to the 18th century,
with their concepts gradually emerging through work on differential
equations and oscillatory phenomena in nature.

The concepts underlying eigenvalues and eigenvectors first appeared in the
18th century, long before the terms "matrix" and "vector" were in use.

They concept first emerged from work on solving differential
equations of the form y' = Ay, which describe various oscillatory
phenomena in nature like mechanical vibrations, light, and sound

STATISTICS

In the early 20th century, eigenvalues and eigenvectors began to be
applied to statistical problems.

Karl Pearson introduced the method of principal components in 1901,
which implicitly used eigenvectors, though he did not use that terminology.

The term "eigenvector" was introduced by David Hilbert in 1904,
 helping to standardize the terminology

Harold Hotelling formalized the use of eigenvectors in statistics
in his 1933 paper on principal component analysis.

/*               _     _
 _ __  _ __ ___ | |__ | | ___ _ __ ___
| `_ \| `__/ _ \| `_ \| |/ _ \ `_ ` _ \
| |_) | | | (_) | |_) | |  __/ | | | | |
| .__/|_|  \___/|_.__/|_|\___|_| |_| |_|
|_|
*/

/*************************************************************************************************************************/
/*                                                                                                                       */
/*                                                                                                                       */
/*          -----------------------------------------------------------------                                            */
/*          |                                                               |                                            */
/*          |    INTERPRETATION OF EIGENVALUES AND EIEGENVECTORS            |                                            */
/*          |                                                               |                                            */
/*          |    The major axis or first eigenvector corresponds to         |                                            */
/*          |    the direction of maximum variance in the data.             |                                            */
/*          |    This direction also defines the first principal            |                                            */
/*          |    component. The amount of explained variance                |                                            */
/*          |    accounted for by the first eigenvector is equal to         |                                            */
/*          |    the first eigenvalue.                                      |                                            */
/*          |                                                               |                                            */
/*          |    The minor axis or second                                   |                                            */
/*          |    eigenvector is perpendiculat to the first and              |                                            */
/*          |    accounts the maximum amount of remaining variance.         |                                            */
/*          |                                                               |                                            */
/*          |    The eigenvectors only provide direction. Positioning       |                                            */
/*          |    the head of the vectors at then means of x and y           |                                            */
/*          |    provides the major and minor axes.                         |                                            */
/*          |                                                               |                                            */
/*          |    INTERESTING FACTS                                          |                                            */
/*          |    =================                                          |                                            */
/*          |                                                               |                                            */
/*          |    1. The linear equation defining the major axis             |                                            */
/*          |       often leads to insights into the data.                  |                                            */
/*          |                                                               |                                            */
/*          |    2. Dimensionality based on variace reduction               |                                            */
/*          |                                                               |                                            */
/*          |    The total variance of the scatter plot is the trace of     |                                            */
/*          |    --                 --                                      |                                            */
/*          |    | var(x)    cov(x,y) |                                2    |                                            */
/*          |    |                    | = trace =var(x)+var(y)+cov(x,y)     |                                            */
/*          |    | cov(x,y)  var(y)   |                                     |                                            */
/*          |    --                 --                                      |                                            */
/*          |    --           --                                            |                                            */
/*          |    | 0.935 0.640 |                   2                        |                                            */
/*          |    |             |=.935 + .691 + .640  =2.0356 Total variance |                                            */
/*          |    | 0.640 0.691 |                                            |                                            */
/*          |    --           --                                            |                                            */
/*          |                                                               |                                            */
/*          |    Proportion of variance accounted for by eigenvalues        |                                            */
/*          |                                                               |                                            */
/*          |    eigenvalue1 = 1.464 => proportion = 1.464/2.0356 =72%      |                                            */
/*          |    eigenvalue2 = 0.162 => proportion = 0.162/2.0356 = 8%      |                                            */
/*          |                                                               |                                            */
/*          |    Total variance for e1 & e2 = (1.464+.162)/2.0356 =80%      |                                            */
/*          |                                                               |                                            */
/*          |    Non-Normalized Eigenvectors                                |                                            */
/*          |                                                               |                                            */
/*          |    From the graph you van see if v is a eigenvector           |                                            */
/*          |    -v is an eigenvector for the same eigenvalue               |                                            */
/*          |                                                               |                                            */
/*          |    Non Normalised Eigenvectors                                |                                            */
/*          |                                                               |                                            */
/*          |     --    --        --       --                               |                                            */
/*          |     | a  b |  then  | a  b    |                               |                                            */
/*          |     | b  d |        | b  -1/b |  perpendicular                |                                            */
/*          |     --    --        --       --                               |                                            */
/*          |                                                               |                                            */
/*          |     For given eigenvalues both are eigenvectors               |                                            */
/*          |                                                               |                                            */
/*          |     --    --       --     --                                  |                                            */
/*          |     | a  b |   or  | -a  -b |                                 |                                            */
/*          |     | b  d |       |- b  -d |                                 |                                            */
/*          |     --    --       --     --                                  |                                            */
/*          |                                                               |                                            */
/*          ---+-------------+-------------+--------------------------------|                                            */
/*          | -5             0             5                                |                                            */
/*        4 +                X                                              + 4                                          */
/*          |                                       INPUT                   |                                            */
/*          |                   Major Axis          -----                   |                                            */
/*          |                        . /     --                 --          |                                            */
/*          |                       ../      | var(x)    cov(x,y) |         |                                            */
/*          |                      ../       | cov(x,y)  var(y)   |         |                                            */
/*          |                     ../..      --                 --          |                                            */
/*        2 +                 ,..../..                                      + 2                                          */
/*          |                 ..../.         --            --               |                                            */
/*          |              ....../..         | 0.935  0.630 |               |                                            */
/*          |             ....../..          | 0.640  0.691 |               |                                            */
/*          |   Minor Axiss..../... .        --            --               |                                            */
/*          |           * ..../....                                         |                                            */
/*          |           ..*../...  .         mean(x) = -.002                |                                            */
/*      Y 0 +           ....C*...            mean(y) = -.009                + 0 Y                                        */
/*          |          ..../...* .                                          |                                            */
/*          |          .../......*               OUTPUT                     |                                            */
/*          |         .../,,,,                   ------                     |                                            */
/*          |         ../. . .                                              |                                            */
/*          |       .../....                 Eigenvalues                    |                                            */
/*          |        ./ ...                                                 |                                            */
/*       -2 +      ../.. .                    e1= 1.464                     + 2                                          */
/*          |      ./.                        e2= 0.162                     |                                            */
/*          |     ./                                                        |                                            */
/*          |     /                          Non Normal Eigenvectors        |                                            */
/*          |                                                               |                                            */
/*          |                                 [1  -1.209] 1/.827=1.209      |                                            */
/*          |                                 [1   0.827]                   |                                            */
/*       -4 +                                                               + 4                                          */
/*          |                                                               |                                            */
/*          |                           Normalized Eigenvectors             |                                            */
/*          |                                                               |                                            */
/*          |                   --            --    --            --        |                                            */
/*          |                   |  .770  -.637 | or | -.770  +.637 |        |                                            */
/*          |                   |  .637   .770 |    | -.637  -.770 |        |                                            */
/*       -5 |                   --            --    --            --        | 5                                          */
/*          |                                                               |                                            */
/*          |                                                               |                                            */
/*          ---+-------------+-------------+---------------------------------                                            */
/*            -5             0             5                                                                             */
/*                           X                                                                                           */
/*                                                                                                                       */
/*************************************************************************************************************************/

/*                   _
(_)_ __  _ __  _   _| |_
| | `_ \| `_ \| | | | __|
| | | | | |_) | |_| | |_
|_|_| |_| .__/ \__,_|\__|
        |_|
*/

libname sd1 "d:/sd1";

options validvarname=any;
libname sd1 "d:/sd1";
data sd1.have;
  call streaminit(54321);
  do ob=-4 to 4 by .01;
    x = rand('normal');
    y = 0.7 * x + rand('normal',0,.5);
    drop ob;
    output;
  end;
run;quit;

/**************************************************************************************************************************/
/*                                                                                                                        */
/*  SD1.HAVE total obs=801                                                                                                */
/*                                                                                                                        */
/*  Obs        x           y                                                                                              */
/*                                                                                                                        */
/*    1    -1.36487     0.02943                                                                                           */
/*    2     0. 51233     0.05964                                                                                          */
/*    3    -0.89565     0.17128                                                                                           */
/*    4    -0.41032    -0.41653                                                                                           */
/*    ....                                                                                                                */
/*  797     0.14316     0.57919                                                                                           */
/*  798    -1.48765    -1.59593                                                                                           */
/*  799     1.44784     1.00016                                                                                           */
/*  800    -0.18555    -0.99234                                                                                           */
/*  801     0.04708    -0.59225                                                                                           */
/*                                                                                                                        */
/**************************************************************************************************************************/

/*   _   _
/ | | |_| |__   ___  ___  _ __ _   _
| | | __| `_ \ / _ \/ _ \| `__| | | |
| | | |_| | | |  __/ (_) | |  | |_| |
|_|  \__|_| |_|\___|\___/|_|   \__, |
                               |___/
*/

DERIVATION OF EIGENVALUES AND EIGENVECTORS
--------------------------------------------

/*    _                            _
  ___(_) __ _  ___ _ ____   ____ _| |_   _  ___  ___
 / _ \ |/ _` |/ _ \ `_ \ \ / / _` | | | | |/ _ \/ __|
|  __/ | (_| |  __/ | | \ V / (_| | | |_| |  __/\__ \
 \___|_|\__, |\___|_| |_|\_/ \__,_|_|\__,_|\___||___/
        |___/
*/

Let

 --                 --     --    --
 | var(x)    cov(x,y) |    | a  b |
 |                    |  = | b  d |
 | cov(x,y)  var(y)   |    --    --
 --                 --
 --    --   --            --
 | a  b |   | 0.935  0.630 |
 | b  d | = | 0.640  0.691 |
 --    --   --            --

EIGENVALUES

Definition
  Determinant = 0

--     --
| a-e  b |                2
|        | = (a-e)(d-e)- b   = 0
| b  d-e |
--     --

Quadratic formula for eigenvalues
               __________________
              / 2          2   2
    a   d   \/ a - 2ad + 4b + d
e1  - + - + ---------------------
    2   2           2
               __________________
              / 2          2   2
    a   d   \/ a - 2ad + 4b + d
e2  - + - - ---------------------
    2   2           2

data _null_;;

  a=0.9351275;   b=0.6397438;
  c=0.6397438;   d=0.6913029;

  Eigenvalue1 = a/2 + d/2 + sqrt(a**2 - 2*a*d + 4*b**2 + d**2)/2;
  Eigenvalue2 = a/2 + d/2 - sqrt(a**2 - 2*a*d + 4*b**2 + d**2)/2;

  put / Eigenvalue1 = / Eigenvalue2 = /;

run;quit;

Eigenvalue1=1.4644714772
Eigenvalue2=0.1619589228

/*    _                                _              _ _                            _                           _
  ___(_) __ _  ___ _ ____   _____  ___| |_ ___  _ __ | (_)_ __   ___  __ _ _ __   __| | ___ _ __   ___ _ __   __| | ___ _ __   ___ ___
 / _ \ |/ _` |/ _ \ `_ \ \ / / _ \/ __| __/ _ \| `__|| | | `_ \ / _ \/ _` | `__| / _` |/ _ \ `_ \ / _ \ `_ \ / _` |/ _ \ `_ \ / __/ _ \
|  __/ | (_| |  __/ | | \ V /  __/ (__| || (_) | |   | | | | | |  __/ (_| | |   | (_| |  __/ |_) |  __/ | | | (_| |  __/ | | | (_|  __/
 \___|_|\__, |\___|_| |_|\_/ \___|\___|\__\___/|_|   |_|_|_| |_|\___|\__,_|_|    \__,_|\___| .__/ \___|_| |_|\__,_|\___|_| |_|\___\___|
        |___/                                                                              |_|

*/

DEFINITION OF EIGENVECTORS

--     --  --   --   --                 --    -- --
| a-e1 b | | v11 |   | (a-e1)v11 + bv21  |    | 0 |
|        | |     |   |                   |  = |   |
| b  d-e1| | v21 |   |  bv11 + (d-e1)v21 |    | 0 |
--     --  --   --   --                 --    -- --

 (a-e1)v11 + bv21  = 0

 bv11 + (d-e1)v21  = 0

The equations above have an infinite number of solutions.
The equations are linearly dependent.
[0,0] is one of them, but [0,0] is not a valid eigenvector.

 For a given eigenvalue if one row is a constant times the
 other row then the eigenvectors are linearly dependent.

 We will show that there is a constant F where F*equ1=equ2

 Let F be the constant

     F ((a-e1)*v11 + b*v21) = b*v11 + (d-e1)*v21

    expanding

    F*(a-e1)*v11 +      F*b*v21 =
           b*v11  +  (d-e1)*v21

    this requires

    F*(a-e1) = b          (v11 terms)
    F*b     =  d-e1       (v21 terms)

    solving First equation

    F = b/(a-e1)

    subsituting in second equation  (Fb=d-e1)

    b**2/(a-e1) = d-e1

    Rearranging

    (a-e1)*(d-e1) - b**2 = 0

 This was condition we used to derive the
 eigenvalues so it must be true. and the
 rows are linearly depedent, See the dterminant
 requirement above.

* COMPUTATIONAL PROOF;
* ===================;

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

The coeficients are equal so F*equ1 = equ2 and
we have linead dependence

coef_v11_equ1=0.6397438     =  coef_v11_equ2=0.6397438
coef_v21_equ1=-0.773168577  =  coef_v21_equ2=-0.773168577

/*   _           _                 _                                _
  __| | ___ _ __(_)_   _____   ___(_) __ _  ___ _ ____   _____  ___| |_ ___  _ __ ___
 / _` |/ _ \ `__| \ \ / / _ \ / _ \ |/ _` |/ _ \ `_ \ \ / / _ \/ __| __/ _ \| `__/ __|
| (_| |  __/ |  | |\ V /  __/|  __/ | (_| |  __/ | | \ V /  __/ (__| || (_) | |  \__ \
 \__,_|\___|_|  |_| \_/ \___| \___|_|\__, |\___|_| |_|\_/ \___|\___|\__\___/|_|  |___/
                                     |___/
*/

 Since the rows are lineraly dependent we only need one
 equation to derive for both eigenvectors.

 Knowing just the eigenvector is like knowing just the slope
 and trying to place the line on a graph.

 We only need to assume a value for V11 to find the value of v21

 --     --  --   --   --                --   -- --
 | a-e2 b | | v12 |   | (a-e2)v12 + bv22 |   | 0 |
 |        | |     | = |                  | = |   |
 | b  d-e2| | v22 |   |  bv12 (d-e2)v22  |   | 0 |
 --     --  --   --   --                --   -- --

 Let v12=1 and solve for v22

     (a-e1)v12 + bv22 = 0

     v22 = -(a-e1)/b

     Remember

     e1 = 0.1619589228
     e2 = 1.4644714772

      --                 --   --    --   --            --
      | var(x)    cov(x,y) |= | a  b | = | 0.935  0.630 |
      | cov(x,y)  var(y)   |  | b  d |   | 0.640  0.691 |
      --                 --   --    --   --            --

      substuting

      v22 = -(0.9351275-0.1619589228 )/0.6397438 = -1.208559703

      v2 = [v12 v22] = [1,-1.208559703]


  --     --  --   --   --                --    -- --
  | a-e2 b | | v11 |   | (a-e2)v11 + bv21 |    | 0 |
  |        | |     |   |                  |  = |   |
  | b  d-e2| | v21 |   |  bv11 (d-e2)v21  |    | 0 |
  --     --  --   --   --                --    -- --

  Let v11 = 1 and solve for v21


  Let v11=1

     (a-e2)v11 + bv21 = 0

     v21 = -(a-e2)/b

    v21 = -(0.9351275 - 1.4644714772)/0.6397438 = 0.8274311954

    v1 = [v11 v21] = [1, 0.8274311954]

/*                                _ _               _                                _
 _ __   ___  _ __ _ __ ___   __ _| (_)_______   ___(_) __ _  ___ _ ____   _____  ___| |_ ___  _ __ ___
| `_ \ / _ \| `__| `_ ` _ \ / _` | | |_  / _ \ / _ \ |/ _` |/ _ \ `_ \ \ / / _ \/ __| __/ _ \| `__/ __|
| | | | (_) | |  | | | | | | (_| | | |/ /  __/|  __/ | (_| |  __/ | | \ V /  __/ (__| || (_) | |  \__ \
|_| |_|\___/|_|  |_| |_| |_|\__,_|_|_/___\___| \___|_|\__, |\___|_| |_|\_/ \___|\___|\__\___/|_|  |___/
                                                      |___/
*/

We have both eigenvectors

    v1 = [v11 v21] = [1, 0.8274311954]    column 1
    v2 = [v12 v22] = [1,-1.208559703]      volumn 2

Lets normalize the eigenvectors

   len_v1  = sqrt(v11**2 + v12**2) = sqrt(1 + 0.8274311954**2)
           = 1.2979377424

   len_v2  = sqrt(v11**2 + v12**2) = sqrt(1 + 1.208559703**2)
           = 1.5686352526

   nrm_v11  = 1/1.2979377424             = 0.7704529789
   nrm_v21  = 0.8274311954/1.2979377424  = 0.637496829

   nrm_v12  = 1/1.5686352526             = 0.6374968294
   nrm_v22  = -1.208559703/1.5686352526  = -0.770452979

data _null_;

    v11 = 1;  v21=0.8274311954;
    v12 = 1;  v22=-1.208559703;

    len_v1  = sqrt(v11**2 + v21**2);
    len_v2  = sqrt(v12**2 + v22**2);

    put len_v1= / len_v2=;

     nrm_v11 = v11/len_v1;
     nrm_v21 = v21/len_v1;

     nrm_v12 = v12/len_v2;
     nrm_v22 = v22/len_v2;

    put
       nrm_v11
       nrm_v21
               /
       nrm_v12
       nrm_v22
    ;
run;quit;

/*___                                                _
|___ \   ___  __ _ ___    ___ ___  _ __  _ __  _   _| |_ ___
  __) | / __|/ _` / __|  / __/ _ \| `_ \| `_ \| | | | __/ _ \
 / __/  \__ \ (_| \__ \ | (_| (_) | | | | |_) | |_| | ||  __/
|_____| |___/\__,_|___/  \___\___/|_| |_| .__/ \__,_|\__\___|
                                        |_|
      _                                _
  ___(_) __ _  ___ _ ____   _____  ___| |_ ___  _ __ ___
 / _ \ |/ _` |/ _ \ `_ \ \ / / _ \/ __| __/ _ \| `__/ __|
|  __/ | (_| |  __/ | | \ V /  __/ (__| || (_) | |  \__ \
 \___|_|\__, |\___|_| |_|\_/ \___|\___|\__\___/|_|  |___/
        |___/
*/

 data _null_;;

    a=0.9351275;   b=0.6397438;
    c=0.6397438;   d=0.6913029;
                                                                                    9
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

Eigenvalues
e1=0.1619589228
e2=0.1619589228

Raw Eigenvectors
v11=1 v21=-0.827431195
v12=1 v22=1.2085597034

Normalized Eigenvectors
norm11=0.7704529789 norm21=-0.637496829
norm12=0.6374968293 norm22=0.7704529789

/*____                                          _               _       _
|___ /   _ __    ___ ___  _ __ ___  _ __  _   _| |_ ___   _ __ | | ___ | |_
  |_ \  | `__|  / __/ _ \| `_ ` _ \| `_ \| | | | __/ _ \ | `_ \| |/ _ \| __|
 ___) | | |    | (_| (_) | | | | | | |_) | |_| | ||  __/ | |_) | | (_) | |_
|____/  |_|     \___\___/|_| |_| |_| .__/ \__,_|\__\___| | .__/|_|\___/ \__|
                                   |_|                   |_|
      _                                _
  ___(_) __ _  ___ _ ____   _____  ___| |_ ___  _ __ ___
 / _ \ |/ _` |/ _ \ `_ \ \ / / _ \/ __| __/ _ \| `__/ __|
|  __/ | (_| |  __/ | | \ V /  __/ (__| || (_) | |  \__ \
 \___|_|\__, |\___|_| |_|\_/ \___|\___|\__\___/|_|  |___/
        |___/
*/

* SAME AS ABOVE;

libname sd1 "d:/sd1";

options validvarname=any;
libname sd1 "d:/sd1";
data sd1.have;
  call streaminit(54321);
  do ob=-4 to 4 by .01;
    x = rand('normal');
    y = 0.7 * x + rand('normal',0,.5);
    drop ob;
    output;
  end;
run;quit;


%utl_rbeginx;
parmcards4;
library(haven);
# Set a random seed for reproducibility
set.seed(123)

# Generate correlated data for an elliptical scatter plot
n <- 1000
x <- rnorm(n)
y <- 0.7 * x + rnorm(n, sd = 0.5)
data <- cbind(x, y)
data<-read_sas("d:/sd1/have.sas7bdat")[,c(1,2)]
data
# Create the scatter plot
pdf("d:/pdf/elp.pdf");
plot(data, asp = 1, main = "Elliptical Scatter Plot with Eigenvectors")

# Compute the covariance matrix
cov_matrix <- cov(data)
cov_matrix;
# Compute eigenvectors and eigenvalues
eigen_result <- eigen(cov_matrix)
eigenvectors <- eigen_result$vectors
eigenvalues <- eigen_result$values

# Calculate the center of the data
center <- colMeans(data)
# Function to plot eigenvectors
plot_eigenvector <- function(eigenvector, eigenvalue, color) {
  scaling_factor <- sqrt(eigenvalue) * 2  # Adjust scaling for visibility
  arrows(center[1], center[2],
         center[1] + scaling_factor * eigenvector[1],
         center[2] + scaling_factor * eigenvector[2],
         col = color, lwd = 2, length = 0.1)
}

# Plot eigenvectors
plot_eigenvector(eigenvectors[,1], eigenvalues[1], "red")
plot_eigenvector(eigenvectors[,2], eigenvalues[2], "blue")

# Add legend
legend("topright", legend = c("First Eigenvector", "Second Eigenvector"),
       col = c("red", "blue"), lwd = 2)

# Print eigenvectors and eigenvalues
print(center)
print("Eigenvectors:")
print(eigenvectors)
print("Eigenvalues:")
print(eigenvalues)
pdf()
;;;;
%utl_rendx;


> print(center)
-0.001784439  0.008892836

> print("Eigenvectors:")
[1] "Eigenvectors:"
> print(eigenvectors)
           [,1]       [,2]
[1,] -0.7704530  0.6374968
[2,] -0.6374968 -0.7704530
> print("Eigenvalues:")
[1] "Eigenvalues:"
> print(eigenvalues)
[1] 1.464471 0.161959
> pdf()
>

/*              _
  ___ _ __   __| |
 / _ \ `_ \ / _` |
|  __/ | | | (_| |
 \___|_| |_|\__,_|

*/
