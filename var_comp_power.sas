/************************************************************************************************************/
/* A SAS program to simulate power for a study comparing variance components between two groups of clusters */
/*																											*/
/* AUTHOR: Brady T. West																					*/
/* DATE: November 1, 2017																					*/
/************************************************************************************************************/

* nreps = number of simulations;
* nint = number of interviewers (clusters) in each group;
* numobs = number of respondents per interviewer;
* intvar1 = interviewer variance component in group 1 (continuous variable)
* pintvar1 = “” (binary variable)
* intvar2 = “” in group 2 (continuous variable)
* pintvar2 = “” (binary variable)
* sigeps1 = within-interviewer variance in group 1 (continuous variable)
* sigeps2 = within-interviewer variance in group 2 (continuous variable);

%macro varcompsim(nreps, nint, numobs, intvar1, pintvar1, intvar2, pintvar2, sigeps1, sigeps2);

%let jplusone = %eval(&nint + 1);
%let total = %eval(2 * &nint);

%do reps = 1 %to &nreps;

data simsample (drop = j k u u2);
   do j = 1 to &nint;
      intid = j;
	  u = sqrt(&intvar1)*rannor(0); /* u ~ N(0,intvar1) */
	  u2 = sqrt(&pintvar1)*rannor(0);
	  do k = 1 to &numobs;
         y1 = 1.25 + u + sqrt(&sigeps1)*rannor(0); /* e ~ N(0,sigeps1) */
		 p2 = exp(0.5 + u2) / (1 + exp(0.5 + u2));
		 rand = ranuni(0);
		 if rand <= p2 then y2 = 1;
		 else y2 = 0;
		 mode = 0;
		 output;
	  end;
   end;
   do j = &jplusone to &total;
      intid = j;
	  u = sqrt(&intvar2)*rannor(0); /* u ~ N(0,intvar2) */
	  u2 = sqrt(&pintvar2)*rannor(0);
	  do k = 1 to &numobs;
         y1 = 1.25 + u + sqrt(&sigeps2)*rannor(0); /* e ~ N(0,sigeps2) */
		 p2 = exp(0.5 + u2) / (1 + exp(0.5 + u2));
		 rand = ranuni(0);
		 if rand <= p2 then y2 = 1;
		 else y2 = 0;
		 mode = 1;
		 output;
	  end;
   end;
run;

proc glimmix data = simsample;
   class intid mode;
   model y1 = / solution;
   random int / subject = intid group = mode;
   covtest homogeneity;
   ods output covtests = lrtresult;
run;

proc glimmix data = simsample;
   class intid mode;
   model y2 = / solution link = logit dist = binomial;
   random int / subject = intid group = mode;
   covtest homogeneity;
   nloptions tech=nrridg;
   ods output covtests = lrtresult2;
run;

%if &reps = 1 %then %do;
data lrtresults (keep = ProbChiSq power_cont);
   set lrtresult;
   if ProbChiSq < 0.05 then power_cont = 1;
   else power_cont = 0;
run;
data lrtresults2 (keep = ProbChiSq power_bin);
   set lrtresult2;
   if ProbChiSq < 0.05 then power_bin = 1;
   else power_bin = 0;
run;
%end;
%else %do;
data lrtresult (keep = ProbChiSq power_cont);
   set lrtresult;
   if ProbChiSq < 0.05 then power_cont = 1;
   else power_cont = 0;
run;
data lrtresults;
   set lrtresults lrtresult;
run;
data lrtresult2 (keep = ProbChiSq power_bin);
   set lrtresult2;
   if ProbChiSq < 0.05 then power_bin = 1;
   else power_bin = 0;
run;
data lrtresults2;
   set lrtresults2 lrtresult2;
run;
%end;

%end;

proc means data = lrtresults mean;
   var power_cont;
run;

proc means data = lrtresults2 mean;
   var power_bin;
run;

%mend varcompsim;

%varcompsim(nreps=100, nint=30, numobs=30, intvar1=1, intvar2=7, pintvar1=0.03, pintvar2=0.33, sigeps1=64, sigeps2=64);
