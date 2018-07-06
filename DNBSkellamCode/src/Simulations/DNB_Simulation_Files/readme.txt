*******************************************************

The "readme" file for the replication codes used for the simulation study from the paper:
"Bayesian Dynamic Modeling of High-Frequency Integer Price Changes"
by István Barra, Agnieszka Borowska and Siem Jan Koopman,
to appear in JFEC 2018.

Uncommenting of three lines at the bottom of the main if/else statement in CodeDynamicDNB_Sim.cpp 
(for the Skellam or the DNB model) allows also for running the out-of-sample study.
Note, however, that the in-sample and the out-of-sample data must be provided. The former can be taken from the simulation study.

This file was created by Agnieszka Borowska on Mar 27, 2018.

*******************************************************

This folder contains the following files:
	CodeDynamicDNB_Sim.cpp		The main code executing a single simulation experiment 
								for the Skellam or the DNB model.
	MCMC.h						Definition of the MCMC class (objects and methods)
	MCMC.cpp					The methods of the MCMC class
	OutOfSample.h  				Definition of the OutOfSample class (objects and methods)
	OutOfSample.cpp				The methods of the OutOfSample class
	struct.h

The code requires the following (unprovided) library:
  		gsl (GNU Scientific Library)
(easily available for Linux).
  
  
To compile:
	0. Go to your desired folder to run the computations.
	1. Set and export the path to your gsl lib folder inside the gsl installation folder, typically given by ([username] is your username):
		LD_LIBRARY_PATH=/home/[username]/gsl/lib 
		export LD_LIBRARY_PATH
	2. Use e.g. the g++ compiler from the GNU Compiler Collection to compile the program DNB_sim_pr from the main simulation code CodeDynamicDNB_Sim.cpp simultaneously linking the required libraries:
		g++ -o DNB_sim MCMC.cpp CodeDynamicDNB_Sim.cpp -I/home/[username]/gsl/include -lm -L/home/[username]/gsl/lib -lgsl -lgslcblas

To run:		
	1. Call:
		./DNB_sim_pr T M BurnIn mu phi sigma2 gamma (nu) iSeedNo
		where the input arguments and their types are as follows: 
			int T - the required length of the generated time series 
			int M - the number of MCMC iterations
			int Burn - the length of the burn-in period
			double mu - the value of mu parameter (unconditional mean of h_t)
			double phi - the value of phi parameter (persistence of the AR(1) process x_t)
			double sigma2 - the value of sigma2_eta parameter (variance of the noise of the AR(1) process x_t)
			double gamma - the value of gamma parameter (zero-inflation)
			double nu - the value of nu parameter (optional: degrees of freedom for DNB)
			int iSeedNo - seed number for data generation (index of the simulation)
 
		e.g.
		./DNB_sim_pr 20000 100000 20000 -1.0 0.97 0.25 0.1 15 2
	2. Optionally you can run the program and save the screen output to a logfile:
		./DNB_sim_pr 20000 100000 20000 -1.7 0.97 0.02 0.001 15 2 |& tee DNB_logfile.txt
		./DNB_sim_pr 20000 100000 20000 -1.7 0.97 0.02 0.001 2 |& tee Sk_logfile.txt

Output:
	1. Saved in an automatically created folder OutPut.
	2. Output files indexed by the model type (with a prefix __Sk or __DNB for the Skellam and the DNB model, respectively) and by the iSeedNo used to generate the data (number added to the mode prefix e.g. __DNB32).
	3. Simulated data files:
		[prefix]_SimulatedData.csv
		[prefix]_SimulatedLogInt.csv
		[prefix]_Season.csv
	4. Markov chain of parameter draws:
		[prefix]_EstimationResults.csv
	5. Estimates for the latent states e.g. x_t and s_t (posterior mean, variance, 5% and 95% quantiles, respectively)
		[prefix]_vXEst.csv
		[prefix]_vXVol.csv
		[prefix]_vXLow.csv
		[prefix]_vXHigh.csv
		[prefix]_vSeasonEst.csv
		[prefix]_vSeaonEst.csv
		[prefix]_vSVol.csv		
		[prefix]_vSLow.csv
		[prefix]_vSHigh.csv 
	6. Simulation time:
		[prefix]_vTime.csv
		
		
		