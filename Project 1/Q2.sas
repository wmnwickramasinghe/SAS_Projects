DATA values;
 p=0.1; n=100; alpha=0.05; z=QUANTILE('normal',1-alpha/2);
 nruns=4000; /* Because SE of the estimated coverage probability should not exceed 0.005 */  seed=98638;
RUN;

DATA montecarlo;
 SET values;
 CALL streaminit(seed);
 DO sample=1 TO nruns;
	x=RAND('binomial',p,n); /* sample */
	phat=x/n;  /* Estimate for p */
	lb=phat - z*sqrt((phat*(1-phat)/n));
	ub=phat + z*sqrt((phat*(1-phat)/n));
	indicator=(lb<=p<=ub); /* Indicator is 1 if p lies with the confidence interval, otherwise 0 */
	OUTPUT;
 END;
RUN;

PROC PRINT DATA=montecarlo (OBS=5);
 VAR sample n p x phat lb ub indicator;
 TITLE1 'Output for Question 2';
 TITLE2 'Part of the dataset generated for Monte Carlo';
RUN;

PROC MEANS DATA=montecarlo NOPRINT;
 VAR indicator;
 OUTPUT OUT=results MEAN=coverage; /* Estimating coverage probability using the proportion of Indicators */
RUN;

PROC PRINT DATA=results;
 VAR coverage; 
 TITLE 'Coverage probability for n = 100 and p = 0.1';
RUN;