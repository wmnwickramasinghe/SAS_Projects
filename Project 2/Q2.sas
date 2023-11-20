filename prostate "/folders/myfolders/Project2/Prostate.dat";

DATA prostate;
INFILE prostate;
INPUT ID PSA CVOL WEIGHT AGE BPH SVI CP GS;
RUN;

PROC PRINT DATA=prostate (OBS=50);
RUN;

PROC REG Data=prostate; /*Model PSA vs CVOL and CP*/
MODEL PSA = CVOL CP;
OUTPUT OUT=D RSTUDENT=R PREDICTED=P; 
RUN;

PROC REG DATA=prostate; 
MODEL PSA=CVOL CP / lackfit ; /* Lack of fit test */
RUN;

DATA D; SET D;
AbsR=ABS(R); 

PROC PLOT DATA = D;
PLOT AbsR*P; /* absolute residuals vs fitted values to check homogeneity assumption */
RUN;

PROC UNIVARIATE DATA=prostate; /* Get Median of CVOL variable*/
VAR CVOL;
VAR CP;
RUN;

DATA D; SET D; 
Group = (CVOL> 4.2631); 
RUN;

/* Brown Forsythe Test for homogeneity for variable CVOL*/
PROC GLM Data=D; 
                class Group;
                model R=Group;
                means Group / hovtest=BF; 
run;

DATA D; SET D; 
Group = (CP> 0); /* Median = 0 from above */
RUN;

/* Brown Forsythe Test for homogeneity for variable CP*/
PROC GLM Data=D; 
                class Group;
                model R=Group;
                means Group / hovtest=BF; 
run;

PROC MODEL DATA=D;/* Breusch Pagan Test for homogeneity */
PARMS b0 b1 b2;
PSA=b0+b1*CVOL+b2*CP; /* specify model to be fitted */
fit PSA /WHITE BREUSCH=(CVOL CP);
fit PSA /BREUSCH=(CVOL);
fit PSA /BREUSCH=(CP);
RUN;

PROC UNIVARIATE DATA=D NORMAL PLOT; /* Check normality of Studentized residuals */
VAR R;
RUN;

********************************************************************************************;
/*Applying log transformation */

DATA prostate; SET prostate;
logPSA = log(PSA);
RUN;

PROC REG DATA=prostate; 
MODEL logPSA=CVOL CP / lackfit ; /* Lack of fit test */
RUN;

proc Print data=prostate (obs=25);
run;

PROC REG Data=prostate;
MODEL logPSA = CVOL CP; 
OUTPUT OUT=E RSTUDENT=Rlog PREDICTED=Plog UCLM=logUCI LCLM=logLCI; 
RUN;

DATA E; SET E;
AbslogR=ABS(Rlog); /* save absolute value of residuals */
RUN;

PROC PLOT DATA = E;
PLOT AbslogR*Plog; /* absolute residuals vs fitted values to check homogeneity assumption */
RUN;

PROC UNIVARIATE DATA=prostate; /* Get Median of CVOL variable*/
VAR CVOL;
VAR CP;
RUN;

DATA E; SET E; 
Group1 = (CVOL> 4.2631); /* Median = 4.2631 from above */
RUN;

/* Brown Forsythe Test for homogeneity*/
PROC GLM Data=E; 
                class Group1;
                model Rlog=Group1;
                means Group1 / hovtest=BF; 
run;

DATA E; SET E; 
Group2 = (CP> 0); /* Median = 0 from above */
RUN;

/* Brown Forsythe Test for homogeneity*/
PROC GLM Data=E; 
                class Group2;
                model Rlog=Group2;
                means Group2 / hovtest=BF; 
run;

PROC MODEL DATA=E;/* Breusch Pagan Test for homogeneity */
PARMS b0 b1 b2;
PSA=b0+b1*CVOL+b2*CP; /* specify model to be fitted */
fit logPSA /WHITE BREUSCH=(CVOL CP);
fit logPSA /BREUSCH=(CVOL);
fit logPSA /BREUSCH=(CP);
RUN;

PROC UNIVARIATE DATA=E NORMAL PLOT; /* Check normality of Studentized residuals */
VAR Rlog;
RUN;

******************************************************************;
/*Applying square root transformation */

DATA prostate; SET prostate;
sqrtPSA = sqrt(PSA);
RUN;

PROC REG DATA=prostate; 
MODEL sqrtPSA=CVOL CP / lackfit ; /* Lack of fit test */
RUN;

PROC REG Data=prostate;
MODEL sqrtPSA = CVOL CP; 
OUTPUT OUT=E RSTUDENT=R PREDICTED=P; 
RUN;

DATA E; SET E; 
Group1 = (CVOL> 4.2631); /* Median = 4.2631 from above */
RUN;

/* Brown Forsythe Test for homogeneity*/
PROC GLM Data=E; 
                class Group1;
                model R=Group1;
                means Group1 / hovtest=BF; 
run;

DATA E; SET E; 
Group2 = (CP> 0); /* Median = 0 from above */
RUN;

/* Brown Forsythe Test for homogeneity*/
PROC GLM Data=E; 
                class Group2;
                model R=Group2;
                means Group2 / hovtest=BF; 
run;


PROC UNIVARIATE DATA=E NORMAL PLOT; /* Check normality of Studentized residuals */
VAR R;
RUN;


****************************************************************************************;
/*Applying inverse transformation */

DATA prostate; SET prostate;
invPSA = (1/PSA);
RUN;

PROC REG DATA=prostate; 
MODEL invPSA=CVOL CP / lackfit ; /* Lack of fit test */
RUN;

PROC REG Data=prostate;
MODEL invPSA = CVOL CP; 
OUTPUT OUT=E RSTUDENT=R PREDICTED=P; 
RUN;

DATA E; SET E; 
Group1 = (CVOL> 4.2631); 
RUN;

/* Brown Forsythe Test for homogeneity of var CVOL*/
PROC GLM Data=E; 
                class Group1;
                model R=Group1;
                means Group1 / hovtest=BF; 
run;

DATA E; SET E; 
Group2 = (CP> 0); 
RUN;

/* Brown Forsythe Test for homogeneity of var CP*/
PROC GLM Data=E; 
                class Group2;
                model R=Group2;
                means Group2 / hovtest=BF; 
run;


PROC UNIVARIATE DATA=E NORMAL PLOT; /* Check normality of Studentized residuals */
VAR R;
RUN;

******************************************************************************;
/* Calculating Mean CI for Quantiles*/

PROC MEANS DATA=Prostate Q1 MEDIAN Q3 N; /*get summary statistics for Variable c Vol*/
VAR CVOL; 
OUTPUT OUT=stat Q1=Q1 MEDIAN=Q2 Q3=Q3 N=N;
RUN;

DATA stat (KEEP=Q1 Q2 Q3); SET stat; 
DO i=1 TO N; OUTPUT; END;
RUN;

DATA CI2;  
MERGE stat E; 
LCI=EXP(LogLCI); UCI=EXP(LogUCI); 
RUN;

PROC PRINT DATA=CI2; /*Dispaly the CI for mean response of 1st, 2nd, 3rd quartiles of CV*/
TITLE "CI for Quantile values of Variable CVOL";
VAR CVOL CP LCI UCI;
WHERE CVOL=Q1 OR CVOL=Q2 OR CVOL=Q3;
RUN;

