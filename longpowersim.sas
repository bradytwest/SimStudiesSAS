/* Simulation program that computes empirical power for GEE logistic regression models fitted to
longitudinal data sets, with two binary covariates (time and group) and their interaction */

/* Author: Brady T. West */

/* Date: June 5, 2014 */

/* Simulation macro */

/* INPUT PARAMETERS: 

nreps = number of simulations
nint = number of subjects in each group (e.g., number of intervention cases, number of control cases)
bsubvar = between-subject variance (defined as (rho x pi^2/3) / (1 - rho), rho is exchangeable correlation) 
int = desired intercept coefficient
betat = desired coefficient for time (post = 1, pre = 0) when group = 0 (time effect for control group)
betag = desired coefficient for group (int = 1, cont = 0) when time = 0 (group effect at time 1)
betatg = desired interaction coefficient for time x group
alpha = desired significance level

*/

%macro geesim(nreps, nint, bsubvar, int, betat, betag, betatg, alpha);

%let jplusone = %eval(&nint + 1);
%let total = %eval(2 * &nint);

/* Repeatedly simulate samples based on the given model */

%do reps = 1 %to &nreps;

data simsample (drop = j k u2);
   do j = 1 to &nint;
      intid = j;
	  u2 = sqrt(&bsubvar)*rannor(0);
	  do k = 1 to 2; /* NOTE: two time points per subject */
		 rand = ranuni(0);
		 group = 0;
		 if k = 1 then post = 0;
		 else if k = 2 then post = 1;
		 p2 = exp(&int + &betat * post + &betag * group + &betatg * group*post + u2) / (1 + exp(&int + &betat *post + &betag * group + &betatg *group*post + u2));
	     if rand <= p2 then trbz_bin = 1; /* binary outcome */
		 else trbz_bin = 0;
         output;
	  end;
   end;
   do j = &jplusone to &total;
      intid = j;
	  u2 = sqrt(&bsubvar)*rannor(0);
	  do k = 1 to 2;
		 rand = ranuni(0);
		 group = 1;
		 if k = 1 then post = 0;
		 else if k = 2 then post = 1;
		 p2 = exp(&int + &betat * post + &betag * group + &betatg * group*post + u2) / (1 + exp(&int + &betat *post + &betag * group + &betatg *group*post + u2));
	     if rand <= p2 then trbz_bin = 1;
		 else trbz_bin = 0;
		 output;
	  end;
   end;
run;

/* Fit the model of interest to a simulated sample */

proc genmod data = simsample desc;
   class intid;
   model trbz_bin = group post group*post / dist = bin link = logit;
   repeated subject = intid / corr = exch corrw;
   ods output GEEEmpPEst = parmest;
run;

/* Determine whether the interaction of interest was significant at the alpha level */

%if &reps = 1 %then %do;
data parmests (keep = probz power_bin);
   set parmest;
   if probz < &alpha then power_bin = 1;
   else power_bin = 0;
   if Parm = "group*post";
run;
%end;
%else %do;
data parmest (keep = probz power_bin);
   set parmest;
   if probz < &alpha then power_bin = 1;
   else power_bin = 0;
   if Parm = "group*post";
run;
data parmests;
   set parmests parmest;
run;
%end;

%end;

/* Compute empirical power based on the simulated samples */

proc means data = parmests mean;
   var power_bin;
run;

%mend geesim;

/* Example of macro call */

%geesim(nreps=100, nint=1000, bsubvar=0.822, int=-2.666, betat=-0.050, betag=0.626, betatg=-0.608, alpha=0.05);



/* Simulation program that computes empirical power for GEE linear regression models fitted to
longitudinal data sets, with general sets of normally distributed covariates */

/* Author: Brady T. West */

/* Date: June 5, 2014 */

/* Simulation macro */

/* INPUT PARAMETERS: 

nreps = number of simulations
nint = number of subjects on which the dependent variable will be measured
bsubvar = variance of random intercepts
residvar = variance of residuals
int = desired intercept coefficient
betat = desired coefficient for time 
betaPA = desired coefficient for a time-varying covariate (PA)
alpha = desired significance level

*/

%macro geenormsim(nreps, nint, bsubvar, residvar, int, betat, betaPA, alpha);

/* Repeatedly simulate samples based on the given model */

%do reps = 1 %to &nreps;

data simsample (drop = j u2);
   do j = 1 to &nint;
      intid = j;
	  u2 = sqrt(&bsubvar)*rannor(0); * Random subject intercept;
	  do k = 1 to 12; * 12 time points;
	     PA = 3.00 + 1.1475*rannor(0); * Note PA independent of time, with mean 3 and SD 1.1475;
		 engage = &int + &betat * k + &betaPA * PA + u2 + sqrt(&residvar)*rannor(0);
         output;
	  end;
   end;
run;

/* Fit the model of interest to a simulated sample */

proc genmod data = simsample;
   class intid; *intid = person ID;
   model engage = k PA;
   repeated subject = intid / corr = exch corrw;
   ods output GEEEmpPEst = parmest;
run;

/* Determine whether the coefficient of interest was significant at the alpha level */

%if &reps = 1 %then %do;
data parmests (keep = probz power_bin);
   set parmest;
   if probz < &alpha then power_bin = 1;
   else power_bin = 0;
   if Parm = "PA"; * predictor that you want to assess power for;
run;
%end;
%else %do;
data parmest (keep = probz power_bin);
   set parmest;
   if probz < &alpha then power_bin = 1;
   else power_bin = 0;
   if Parm = "PA"; * predictor that you want to assess power for;
run;
data parmests;
   set parmests parmest;
run;
%end;

%end;

/* Compute empirical power based on the simulated samples */

proc means data = parmests mean;
   var power_bin;
run;

%mend geenormsim;

/* Example of macro call */

%geenormsim(nreps=100, nint=50, bsubvar=0.015, residvar=0.112, int=1.38, betat=2, betaPA=0.10, alpha=0.05);
