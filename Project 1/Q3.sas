/* Question 3 */

data values;
 p=0.1; n=25; alpha=0.05; z=quantile('normal',1-alpha/2);
 nruns=4000; /* Because SE of the estimated coverage probability should not exceed 0.005 */  seed=69754;
run;

data montecarlo;
 set values;
 call streaminit(seed);
 do sample=1 to nruns;
	x=rand('binomial',p,n); /* sample */
	phat=x/n;  /* Estimate for p */
	lb=phat - z*sqrt((phat*(1-phat)/n));
	ub=phat + z*sqrt((phat*(1-phat)/n));
	indicator=(lb<=p<=ub); /* Indicator is 1 if p lies with the confidence interval, otherwise 0 */
	output;
 end;
run;

proc print data=montecarlo (obs=5);
 var sample n p x phat lb ub indicator;
 title1 'Output for Question 2';
 title2 'Part of the dataset generated for Monte Carlo';
run;

proc means data=montecarlo noprint;
 var indicator;
 output out=results mean=coverage; /* Estimating coverage probability using the proportion of Indicators */
run;

proc print data=results;
 var coverage; 
 title 'Coverage probability for n = &n and p = 0.1';
run;


