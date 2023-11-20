filename prostate "/folders/myfolders/Project2/Prostate.dat";

DATA prostate;
INFILE prostate;
INPUT ID PSA CVOL WEIGHT AGE BPH SVI CP GS;
RUN;

PROC PRINT DATA=prostate (OBS=50);
RUN;

PROC SGSCATTER DATA=prostate; /* draw scatter plots for PSA level vs other variables */
ods graphics on / width=9.5in height=2.0in;
  COMPARE y=PSA
          x=(CVOL WEIGHT AGE BPH SVI CP GS);
RUN;
ods graphics off;

PROC CORR DATA=prostate;  /* gives simple correlation coefficients */
VAR PSA;
WITH PSA CVOL WEIGHT AGE BPH SVI CP GS;
run;

/* draws the whole scatter plot matrix
proc SGscatter data=prostate; 
matrix ID PSA CVOL WEIGHT AGE BPH SVI CP GS; 
run;

proc corr data=prostate plots=matrix;
var ID PSA CVOL WEIGHT AGE BPH SVI CP GS;
run; */


PROC REG LINEPRINTER; /*Defining the model*/
MODEL PSA=CVOL / R ; 
PAINT RSTUDENT. >2 OR RSTUDENT.<-2 / SYMBOL='*'; 
PLOT PSA*CVOL RSTUDENT.*PREDICTED. ; 
OUTPUT OUT=D RSTUDENT=R PREDICTED=P; 
RUN;

PROC PRINT DATA=D;
RUN;


PROC REG DATA=prostate; 
MODEL PSA=CVOL / lackfit ; /* Lack of fit test */
RUN;

DATA D; SET D;
AbsR=ABS(R); 
RUN;

PROC PLOT DATA = D;
PLOT AbsR*P; /* absolute residuals vs fitted values to check homogeneity assumption */
RUN;

PROC UNIVARIATE DATA=prostate; /* Get Median of CVOL variable*/
VAR CVOL;
RUN;

DATA D; SET D; 
Group = (CVOL> 4.2631); 
RUN;

/* Brown Forsythe Test for homogeneity*/
PROC GLM Data=D; 
                class Group;
                model R=Group;
                means Group / hovtest=BF; 
run;


PROC MODEL DATA=D; /* Breusch Pagan Test for homogeneity */
PARMS b0 b1; 
PSA=b0+b1*CVOL; 
fit PSA / BREUSCH=(CVOL); 
RUN;

PROC UNIVARIATE DATA=D NORMAL PLOT; /* Check normality of Studentized residuals */
VAR R;
RUN;

**********************************************************************;
/*Applying transformation*/

PROC TRANSREG DATA=prostate; /* Find Box-Cox transformation power */
MODEL BoxCox(PSA)=identity(CVOL);
RUN;

/* Performing log transformation*/
DATA prostate; SET prostate;   
logPSA = log(PSA); 
RUN;

PROC REG Data = prostate;
MODEL logPSA = CVOL/ R;
OUTPUT OUT=E RSTUDENT=R1 Predicted=P1 LCLM=LogLCI UCLM=LogUCI; 
RUN;


PROC UNIVARIATE DATA=E NORMAL PLOT; /* Check normality of Studentized residuals */
VAR R1;
RUN;

PROC REG DATA=prostate; 
MODEL logPSA=CVOL  / lackfit ; /* Lack of fit test */
RUN;

DATA E; SET E;
AbsR1=ABS(R1); /* save absolute value of residuals */
RUN;

PROC PLOT DATA = E;
PLOT AbsR1*P1; /* absolute residuals vs fitted values to check homogeneity assumption */
RUN;

PROC UNIVARIATE DATA=prostate; /* Get Median of CVOL variable*/
VAR CVOL;
RUN;

DATA E; SET E; 
Group = (CVOL> 4.2631); 
RUN;

/* Brown Forsythe Test for homogeneity*/
PROC GLM Data=E; 
                class Group;
                model R1=Group;
                means Group / hovtest=BF; 
run;


PROC MODEL DATA=D; /* Breusch Pagan Test for homogeneity */
PARMS b0 b1;
logPSA=b0+b1*CVOL;
fit logPSA / BREUSCH=(CVOL); 
RUN;

***************************************************************************************;

/* Finding CI for 1st,2nd and 3rd Quartiles*/

PROC MEANS DATA=Prostate Q1 MEDIAN Q3 N; /*get summary statistics*/
VAR CVOL; 
OUTPUT OUT=stat Q1=Q1 MEDIAN=Q2 Q3=Q3 N=N;
RUN;

DATA stat (KEEP=Q1 Q2 Q3); SET stat; 
DO i=1 TO N; OUTPUT; END;
RUN;

DATA CI;  
MERGE stat E; 
LCI=EXP(LogLCI); UCI=EXP(LogUCI); 
RUN;

PROC PRINT DATA=CI; /*Dispaly the CI for mean response of 1st, 2nd, 3rd quartiles of CV*/
TITLE "CI for Quantile values of Variable CVOL";
VAR CVOL LCI UCI;
WHERE CVOL=Q1 OR CVOL=Q2 OR CVOL=Q3;
RUN;


