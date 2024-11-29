# Instructions to replicate experiments
This file defines the Python codes used in the paper to generate the results for the experiments illustrated in this paper. In order to run the experiments, please create a virtual Python environment and install the dependencies from the '/scripts/requirements.txt' file using 'pip install -r requirements.txt'. 

## S-WST.py
This reads the arrival rates in the data folder ("ArrivalRates.xlsx") and uses the average service proportions specified in the lists 'reg_props' (for first-time service) and 'crisis_props' (for crisis service), on lines 345-346, to generate the optimal **S-WST** scheduling policies. The policy obtained using an average service proportion of 8% for first-time service and 4% for crisis service is shown in Figure 6 of the manuscript. The policies used for the results in Table 1 and Figure 9 of the manuscript are also obtained using this script. The key parameters are:\
line 345: crisis_props=[i for i in range(3,6)] # range of crisis average service proportions\
line 346: reg_props=[i for i in range(7,10)] # range of first-time average service proportions

## DPS.py
This reads the arrival rates in the data folder ("ArrivalRates.xlsx") and generates the DPS policy shown in Figure 10 of manuscript and used to generate the results in Figure 11 of the manuscript. Note that an active Gurobi license is required to run this script. The key parameters are:\
line 30: U_1 = 58 # upper bound of slots to be allocated to first-time service\
line 31: U_2 = 29 # upper bound of slots to be allocated to crisis service


## MCS_PSA.py 
This reads the arrival rates in the data folder ("ArrivalRates.xlsx") and generates the simulation experiment in Appendix C of the Online Supplement, in order to benchmark the performance of PSA. The key parameters are:\
line 373: n = 500 # Specify the total number of schedules desired (simulation runs)\
line 389: mcs_sims = 1000 # number of monte-carlo simulation runs to estimate the average system time for each schedule

## Solver_comparison.py
This reads the arrival rates in the data folder ("ArrivalRates.xlsx") and uses Gurobi to solve **S-SiST**, and compares the solution quality to that obtained using the solution scheme proposed in this paper. It also computes the total time to solve the problem using Gurobi, and also the propsed scheme. Note that an active Gurobi license is required to run this script. The key parameters are:\
line 294: planning_horizon = 200 # total length of the planning horizon (T)\
line 295: Ub = 50  # Upper bound (U)
