SAS CODES

filename cardio '/folders/myfolders/Project3/cardio.csv';

DATA cardio; 
INFILE cardio DSD FIRSTOBS = 2 ;
INPUT uric dia hdl choles trig alco;
RUN;

PROC PRINT DATA=cardio (OBS=10);
RUN;

*.................................................................................................................................................................................................;
/*Question 1*/

PROC REG Data=cardio; 
MODEL uric = dia hdl choles trig alco; /*Fit a full model for 
predicting uric acid level using other explanatory variables*/
TEST hdl=choles=0 /*Test if the variables hdl and choles can be 
(jointly)dropped together from the full model */;
RUN;

*.................................................................................................................................................................................................;
/*Question 2*/

/* Selection based on Adj R^2 */
PROC REG DATA = cardio; 
MODEL uric = dia hdl choles trig alco / selection = adjrsq rsquare;
RUN;

/* Stepwise selection */
PROC REG DATA = cardio; 
MODEL uric = dia hdl choles trig alco / selection = stepwise; 
RUN;

*.................................................................................................................................................................................................;
/*Question 3*/

/* Regression model for selected variables */
PROC REG Data=cardio; 
MODEL uric = trig alco dia choles;
OUTPUT OUT=D RSTUDENT=R PREDICTED=P; 
RUN;

/* Lack of fit test for linearity*/
PROC REG DATA=cardio; 
MODEL uric = trig alco dia choles / lackfit ; 
RUN;

/* Brown Forsythe Test for homogeneity*/
PROC MEANS DATA=cardio median; /* Get Median for variables*/
VAR trig alco dia choles;
RUN;

/* BF for variable trig*/
DATA D; SET D; 
Group = (trig> 1.03); 
RUN;

PROC GLM Data=D; 
                class Group;
                model R=Group;
                means Group / hovtest=BF; 
run;

/* BF for variable alco*/
DATA D; SET D; 
Group = (alco> 0);
RUN;

PROC GLM Data=D; 
                class Group;
                model R=Group;
                means Group / hovtest=BF; 
run;

/* BF for variable dia*/
DATA D; SET D; 
Group = (dia> 87);
RUN;

PROC GLM Data=D; 
                class Group;
                model R=Group;
                means Group / hovtest=BF; 
run;

/* BF for variable choles*/
DATA D; SET D; 
Group = (choles> 5.5);
RUN;

PROC GLM Data=D; 
                class Group;
                model R=Group;
                means Group / hovtest=BF; 
run;

/* Breusch Pagan Test for homogeneity */
PROC MODEL DATA=D;
PARMS b0 b1 b2 b3 b4;
PSA=b0+b1*trig+b2*alco+b3*dia+b4*choles; 
fit PSA /WHITE BREUSCH=(trig alco dia choles);
fit PSA /BREUSCH=(trig);
fit PSA /BREUSCH=(alco);
fit PSA /BREUSCH=(dia);
fit PSA /BREUSCH=(choles);
RUN;

 /* absolute residuals vs fitted values to check homogeneity 
 assumption */

DATA D; SET D;
AbsR=ABS(R); 

PROC PLOT DATA = D;
PLOT AbsR*P;
RUN;

/* Check normality of Studentized residuals */
PROC UNIVARIATE DATA=D NORMAL PLOT; 
VAR R;
RUN;

/* Tests for outliers*/

/*calcutate hat matrix,Coock's distance, DFFITS and DFBETAS*/
PROC REG Data=cardio; 
MODEL uric = trig alco dia choles / INFLUENCE R;
ods output outputstatistics=results;
RUN;

PROC PRINT Data=results (obs=5);
RUN;

/* Test for outliers using Bonferroni method */
DATA results; set results; 
n=998; p= 5; alpha=0.05;
tvalue = tinv(1-alpha/(2*n), n-p-1);
if (abs(RStudent)) > tvalue then outlier=1;
else outlier=0;
RUN;

PROC PRINT data=results;
where outlier=1;
var RStudent;
RUN;

/* Test for outliers using R^2,hat matrix,Coock's distance,DFFITS and DFBETAS */
DATA results; SET results;

/* Checking if hii > 2*p/n */
if HatDiagonal > 2*(p/n) then hilev=1; 
else hilev=0;

/* Checking if DFFITS > 2*sqrt(p/n) */
if (abs(DFFITS) > 2*sqrt(p/n)) then dfflag=1;
else dfflag=0;

/* Calculating percentile for each Cook's D value using F(p,n-p) */
Fpercent = 100*probf(CooksD, p, n-p); 

/* Checking if each DFBETAS value > 1 */
if (abs(dfb_trig) > 1) then b1flag=1; 
else b1flag=0;
if (abs(dfb_alco) > 1) then b2flag=1;
else b2flag=0;
if (abs(dfb_dia) > 1) then b3flag=1;
else b3flag=0;
if (abs(dfb_choles) > 1) then b4flag=1;
else b4flag=0;
RUN;

PROC PRINT DATA = results (obs=5);
where hilev=1 or dfflag=1 or Fpercent>20 or b1flag=1 or b2flag=1 or b3flag=1 or b4flag=1;
var HatDiagonal hilev DFFITS dfflag CooksD Fpercent dfb_trig b1flag dfb_alco b2flag dfb_dia b3flag dfb_choles b4flag;
RUN;

PROC PRINT DATA = results;
where Fpercent>20;
var  Fpercent ;
RUN;

PROC PRINT DATA = results;
where dfflag=1;
var  DFFITS dfflag ;
RUN;

/* Collinearity diagnostics */
proc SGscatter data=cardio; 
matrix uric trig alco dia choles; 
run;

proc corr data=cardio plots=matrix;
var uric trig alco dia choles;
run;

PROC REG Data=cardio; 
MODEL uric = trig alco dia choles / collin tol vif;
RUN;

/*Transformations*/

/* Find Box-Cox transformation power */
PROC TRANSREG DATA=cardio; 
MODEL BoxCox(uric)=identity(trig alco dia choles);
RUN;

/* Tring log transformation */
DATA D; SET D;
loguric = log(uric);
RUN;

PROC REG Data=D;
MODEL loguric = trig alco dia choles/LACKFIT DWPROB;
OUTPUT OUT=E RSTUDENT=Rlog PREDICTED=Plog R=Res1;
RUN;

PROC MODEL DATA=D; /* Test for homogeneity assumption */
PARMS b0 b1 b2 b3 b4;
loguric = b0 + b1*trig + b2*alco + b3*dia + b4*choles;
fit loguric /WHITE BREUSCH=( trig alco dia choles);
RUN;

PROC UNIVARIATE DATA=E NORMAL PLOT; /* Checking normality of Studentized residuals */
VAR Rlog;
RUN;

/* Try INVSQRT transformation */
DATA D; SET D;
invsqrturic = 1/sqrt(uric);
RUN;
PROC REG Data=D;
MODEL invsqrturic = trig alco dia choles/LACKFIT DWPROB;
OUTPUT OUT=F RSTUDENT=Rsqinv PREDICTED=Psqinv;
RUN;

PROC MODEL DATA=D; /* Test for homogeneity assumption */
PARMS b0 b1 b2 b3 b4;
invsqrturic = b0 + b1*trig + b2*alco + b3*dia + b4*choles;
fit invsqrturic /WHITE BREUSCH=( trig alco dia choles);
RUN;

PROC UNIVARIATE DATA=F NORMAL PLOT; /* Checking normality of Studentized residuals */
VAR Rsqinv;
RUN;

*.................................................................................................................................................................................................;
/*Qestion 4*/

PROC REG Data = cardio;
MODEL uric = trig alco dia choles /R clb;
output out=results r =residual;
RUN;

DATA Step2;
SET results;
absresid = abs(residual);
RUN;

PROC PRINT DATA=Step2(obs=5);
RUN;

PROC REG Data = Step2;
MODEL absresid = trig alco dia choles/p; /* option p requests fitted values */
output out = Step3 p =ehat;
RUN;

DATA STEP3;
SET Step3;
wt = 1/(ehat**2);
RUN;

PROC PRINT DATA=STEP3 (obs=5);
RUN;

PROC REG Data=Step3; /* weighted least squares regression */
MODEL uric = trig alco dia choles/R clb;
WEIGHT wt;
output out=iteration2 r =residual2;
RUN;

PROC PRINT DATA=iteration2 (obs=5);
RUN;

/* Reiterate the process - Iteratively reweighted least squares */

DATA iteration2;
SET iteration2;
absresid2 = abs(residual2);
RUN;

PROC PRINT DATA=iteration2 (obs=5);
RUN;

PROC REG Data = iteration2;
MODEL absresid2 = trig alco dia choles/p; /* option p requests fitted values */
output out = results2 p =ehat2;
RUN;
PROC PRINT DATA=results2 (obs=5);
RUN;

DATA results2;
SET results2;
wt2 = 1/(ehat2**2);
RUN;
PROC PRINT DATA=results2 (obs=5);
RUN;

PROC REG Data=results2; /* weighted least squares regression */
MODEL uric = trig alco dia choles/R clb;
WEIGHT wt2;
output out=iteration3 r =residual3;
RUN;
PROC PRINT DATA=iteration3 (obs=5);
RUN;

/* Reiterate the process - Iteratively reweighted least squares */

DATA iteration3;
SET iteration3;
absresid3 = abs(residual3);
RUN;

PROC PRINT DATA=iteration3 (obs=5);
RUN;

PROC REG Data = iteration3;
MODEL absresid3 = trig alco dia choles/p; /* option p requests fitted values */
output out = results3 p =ehat3;
RUN;
PROC PRINT DATA=results3 (obs=5);
RUN;

DATA results3;
SET results3;
wt3 = 1/(ehat3**2);
RUN;
PROC PRINT DATA=results3 (obs=5);
RUN;

PROC REG Data=results3; /* weighted least squares regression */
MODEL uric = trig alco dia choles/R clb;
WEIGHT wt3;
output out=iteration4 r =residual4;
RUN;
PROC PRINT DATA=iteration4(obs=5);
RUN;

/* Reiterate the process - Iteratively reweighted least squares */

DATA iteration4;
SET iteration4;
absresid4 = abs(residual4);
RUN;

PROC PRINT DATA=iteration4 (obs=5);
RUN;

PROC REG Data = iteration4;
MODEL absresid4 = trig alco dia choles/p; /* option p requests fitted values */
output out = results4 p =ehat4;
RUN;
PROC PRINT DATA=results4 (obs=5);
RUN;

DATA results4;
SET results4;
wt4 = 1/(ehat4**2);
RUN;
PROC PRINT DATA=results4 (obs=5);
RUN;

PROC REG Data=results4; /* weighted least squares regression */
MODEL uric = trig alco dia choles/R clb;
WEIGHT wt4;
output out=iteration5 r =residual5;
RUN;
PROC PRINT DATA=iteration5(obs=5);
RUN;

*.................................................................................................................................................................................................;
/*Qestion 5*/

filename cardio '/folders/myfolders/Project3/cardio.csv';

DATA cardio; 
INFILE cardio DSD FIRSTOBS = 2 ;
INPUT uric dia hdl choles trig alco;
RUN;

PROC PRINT DATA=cardio (OBS=10);
RUN;

PROC REG Data = cardio;
MODEL uric = trig alco dia choles /R clb;
output out=E r =Res1;
RUN;


PROC MEANS Data=E Median;
Var Res1;
RUN;

DATA E;
Set E;
AD1=abs(Res1+7.6743106);
RUN;

PROC MEANS Data=E Median;
Var AD1;
RUN;

DATA E;
Set E;
MAD1=53.5500285/0.6745;
u1=Res1/MAD1; /* Calculating scaled residuals */
If abs(u1) le 4.685 then wt1=(1-(u1/4.685)**2)**2; Else wt1=0; /* Using Bisquare weight function */
RUN;

TITLE "Parameter Estimates from 1st Iteration";

PROC REG Data=E;
Model uric = trig alco dia choles / p;
Weight wt1;
Output Out=E2 R=Res2 P=P2;
RUN;

/* Iteratively Reweighted Least Squares - Second Iteration */
TITLE;
PROC MEANS Data=E2 Median;
Var Res2;
RUN;

DATA E2;
Set E2;
AD2=abs(Res2+5.1140754);
RUN;

PROC MEANS Data=E2 Median;
Var AD2;
RUN;

DATA E2;
Set E2;
MAD2=53.7849334/0.6745;
u2=Res2/MAD2; /* Calculating scaled residuals */
If abs(u2) le 4.685 then wt2=(1-(u2/4.685)**2)**2; Else wt2=0; /* Using Bisquare weight function */
RUN;

TITLE "Parameter Estimates from 2nd Iteration";
PROC REG Data=E2;
Model uric = trig alco dia choles / P;
Weight wt2;
Output Out=E3 R=Res3 P=P3;
RUN;

/* Iteratively Reweighted Least Squares - Third Iteration */

TITLE;
PROC MEANS Data=E3 Median;
Var Res3;
RUN;

DATA E3;
Set E3;
AD3=abs(Res3+5.0397940);
RUN;

PROC MEANS Data=E3 Median;
Var AD3;
RUN;

DATA E3;
Set E3;
MAD3=53.6012022/0.6745;
u3=Res3/MAD3; /* Calculating scaled residuals */
If abs(u3) le 4.685 then wt3=(1-(u3/4.685)**2)**2; Else wt3=0; /* Using Bisquare weight function */
RUN;

TITLE "Parameter Estimates from 3rd Iteration";
PROC REG Data=E3;
Model uric = trig alco dia choles / P;
Weight wt3;
Output Out=E4 R=Res4 P=P3;
RUN;


/* Iteratively Reweighted Least Squares - forth Iteration */

TITLE;
PROC MEANS Data=E4 Median;
Var Res4;
RUN;

DATA E4;
Set E4;
AD4=abs(Res4+5.1657584);
RUN;

PROC MEANS Data=E4 Median;
Var AD4;
RUN;

DATA E4;
Set E4;
MAD4=53.6156608/0.6745;
u4=Res4/MAD4; /* Calculating scaled residuals */
If abs(u4) le 4.685 then wt4=(1-(u4/4.685)**2)**2; Else wt4=0; /* Using Bisquare weight function */
RUN;

TITLE "Parameter Estimates from 4rd Iteration";
PROC REG Data=E4;
Model uric = trig alco dia choles / P;
Weight wt4;
Output Out=E5 R=Res5 P=P4;
RUN;

proc print data=E5 (obs=5);
var res1 wt1 res2 wt2 res3 wt3 res4 wt4 res5; 
run;
*.................................................................................................................................................................................................;
BONUS QUESTION

filename sample '/folders/myfolders/Project3/cardio.csv';

DATA sample; 
INFILE sample DSD FIRSTOBS = 2 ;
INPUT uric dia hdl choles trig alco;
RUN;

Proc means data=sample q1;
var uric; /*1st quartile 239.0000000*/
run;

%let NumSamples = 1000;       /* number of bootstrap resamples */
/* 2. Generate many bootstrap samples */
proc surveyselect data=sample NOPRINT seed=98638
     out=BootSSFreq(rename=(Replicate=SampleID))
     method=urs              /* resample with replacement */
     samprate=1              /* each bootstrap sample has N observations */
     /* OUTHITS                 option to suppress the frequency var */
     reps=&NumSamples;       /* generate NumSamples bootstrap resamples */
run;

/* 3. Compute the statistic for each bootstrap sample */
proc means data=BootSSFreq noprint;
   by SampleID;
   freq NumberHits;
   var uric;
   output out=OutStats Q1=Q1stat;  /* approx sampling distribution */
run;

proc print data=OutStats (obs=5);
run;

/*Visualize the bootstrap distribution*/
title "Bootstrap Distribution";
proc sgplot data=OutStats;
   histogram Q1stat;
run;

proc univariate data=OutStats normal plot;
var Q1stat;
run;


Proc means data=OutStats mean std q1 q3;
var Q1stat;
run;

/*mean=239.9860000, std = 3.9710922 
*/

/*CI*/
proc univariate data=OutStats noprint;
   var Q1stat;
   output out=Pctl pctlpre =CI95_
          pctlpts =2.5  97.5       /* compute 95% bootstrap confidence interval */
          pctlname=Lower Upper;
run;
 
proc print data=Pctl noobs; run;

data OutStats;
set outstats;
dif=Q1stat-239;
run;

proc print data=outstats (obs=5);
run;

proc univariate data=OutStats noprint;
   var dif;
   output out=Pct2 pctlpre =CI95_
          pctlpts =2.5  97.5       /* compute 95% bootstrap confidence interval */
          pctlname=Lower Upper;
run;

proc print data=Pct2 noobs; run;
