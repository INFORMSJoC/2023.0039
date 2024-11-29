#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 09:29:29 2023

@author: hebaish
"""

#%% IMPORT LIBRARIES
import nhpp 
import os
import numpy as np
import pandas as pd
import random
from tqdm.notebook import trange, tqdm
import math
from datetime import datetime
import matplotlib.pyplot as plt
from tqdm import tqdm
#%% INSTRUCTIONS TO GENERATE SCHEDULES
"""
To Generate random schedules do the following: 
    
    1. Generate the minimum schedule by calling MinViableSched(lam). Ouptut is
       an array of 0s and 1s indicating the allocated slots.
       
    2. Get service rates and times of that schedule by calling 
       get_ServTR_vecs(sched). Output is twofold: st & sr
       
    3. Test that minimum schedule for steady state by calling 
       checksteadystate(lam, sr).
       
    4. To generate n random schedules, run generate_scheds(MinSched, U, n), 
       where MinSched is the generated sched in step 1, U is an upper bound on
       the number of slots in the new schedules, n is the number of desired
       schedules. Note: U must be larger than sum(MinSched). 
"""

def get_ServTR_vecs(sched):
    """Takes a schedule as input and returns a vector of service times and 
       service rates, with all values after the last available slot being 0 for
       both vectors.

    Args:
        sched (list): binary vector representing the schedule. 

    Returns:
        list,list: a list of service times and service rates with all values 
        after the last available slot being 0 for both vectors.
    """   
    st = np.zeros(len(sched))
    index_last = next(i for i in range(len(sched) - 1, -1, -1) if sched[i] == 1)
    st[index_last]=1
    for i in range(index_last-1,-1,-1):
        if (sched[i]==1):
            st[i]=1
        if (sched[i]==0):
            st[i]= st[i+1]+1
    sr = np.zeros(len(sched))
    for i in range(index_last+1):
        sr[i] = 1/st[i]
    return st,sr


def checksteadystate(lam,sr):
    """Checks whether the service rate and lambda vectors conform to the steady
       state conditions. Returns 0 if it is VIOLATED and 1 if it CONFORMS to 
       steady state.

    Args:
        lam (list): vector of arrival rates.
        sr (list): vector of service rates.

    Returns:
        int: Returns 0 if it is VIOLATED and 1 if it CONFORMS to steady state.
    """    
    for i in range(len(lam)):
        if (sr[i]<lam[i]):
            return 0, i
    return 1


def MinViableSched(lam):
    """Creates a minimum feasible schedule that satisfy the steady state 
       condition.

    Args:
        lam (list): vector of arrival rates.
        
    Returns:
        Sched: An array of zeros and ones corresponding to the schedule.
    """
    week = 40
    sr = [x/0.9 for x in lam]
    st = [math.ceil(1/x) for x in sr]
    sched = [0 for i in range(len(lam))]

    for i in range(0, len(lam), week):
        if st[i]<=week:
            for j in range(0,week):
                if (j+1)%st[i]==0:
                    sched[i+j]=1
                sched[i+week-1]=1
        else:
            sched[i+week-1]=1
    return sched


def generate_scheds(MinSched, U, n):
    """Creates n random feasible schedules.

    Args:
        MinSched (list): A vector of the minimum feasible schedule.
        U (int): An upper bound of the number of slots to allocate.
        n (int): The required number of schedules
        
    Returns:
        results (list of list): An array of n lists, each containing a feasible
        schedule.
    """
    
    if U <= sum(MinSched):
        return f'Error: U should be greater than {sum(MinSched)}.'
    else:
        results = []
        for i in range(n):
            new_arr = MinSched.copy()
            change = random.randint(1, U-sum(MinSched))
            for j in range(change):
                index = random.choice([i for i, x in enumerate(new_arr) if x == 0])
                new_arr[index] = 1
            results.append(new_arr)
        return results

#%% INSTRUCTIONS TO RUN MONTE CARLO SIMULATION
"""
To run a Monte Carlo Simulation of the system, you will need to use a schedule
and an arrival rate array. Follow the below steps:
    
    1. Generate a schedule as previously discussed. 
    
    2. Run MonteCarloSim(lam, Sched, nreps), where lam is the arrival rate
       array, Sched is the an array of 0s and 1s of the schedule, nreps is the
       desired number of replications.    
"""

def MonteCarloSim(lam, Sched, nreps):
    """A Monte Carlo Simulation of the queueing system.

    Args:
        lam (list): An array of arrival rates.
        Sched (list): A list of zeros and ones indicating the schedule.
        nreps (int): The required number of replications.
        
    Returns:
        mean waiting time (float)
        mean system time (float)
    """
    
    Mean_Wait_Times = []
    Mean_Sys_Times = []
    
    for rep in range(nreps):
        
        # Creating the schedule from the array of zeros and ones
        Schedd = [i for i, x in enumerate(Sched) if x == 1]
        rejects = 0
        PatID = []
        ArrTim = []
        SchedTim = []
        
        # Creating a list of random arrivals according to a nonstationary 
        # Poisson process.
        ID = 0
        for t in range(len(lam)):
            arrivals = np.random.poisson(lam[t])
            if arrivals > 0:
                for a in range(arrivals):
                    PatID.append(ID)
                    ID+=1
                    ArrTim.append(t)
        
        # Scheduling patients
        Schedd = [i for i, x in enumerate(Sched) if x == 1]
        t_hat = Schedd[-1]
        for i in range(len(PatID)):
            if ArrTim[i]>t_hat:
                rejects +=1
                break
            
            if len(Schedd)<1:
                rejects +=1
                break
            
            time_to_sched = next((x for x in Schedd if x > ArrTim[i]), 0)
            SchedTim.append(time_to_sched)
            if time_to_sched in Schedd:
                Schedd.remove(time_to_sched)
        
        # Creating a list of patient IDs, arrival time, and scheduling time
        Full = [list(t) for t in zip(PatID, ArrTim, SchedTim)]
        
        wait_times = []
        sys_times = []
        for i in range(len(Full)):
            wait_times.append(Full[i][2]-Full[i][1])
            sys_times.append(Full[i][2]-Full[i][1]+1)
        
        mean_wait_time = np.mean(wait_times)
        mean_sys_time = np.mean(sys_times)
        Mean_Wait_Times.append(mean_wait_time)
        Mean_Sys_Times.append(mean_sys_time)
    
    return np.mean(Mean_Wait_Times), np.mean(Mean_Sys_Times)

def SimSystem(lam,x):
    """Takes a non-stationary arrival rate vector and a schedule and returns the mean system time of accepted items, total number of rejected items and the total number of arrivals (for debugging). 

    Args:
        lam (list): list of floats where each time is the corresponding arrival rate for that time slot.
        x (list): schedule, i.e., binary array of time slots with 1s in available slots.

    Returns:
        float, int, int: mean system time of accepted patients, total number of rejections and total number of arrivals.
    """    
    arr_vec =  getArrV(lam)
    tot_sysTime = getSysTime(arr_vec,x)
    rejs = getRejections(arr_vec, x)
    tot_arrs = sum(arr_vec)
    if (tot_arrs-rejs)==0:
        mean_sysTime = 0 
    else: 
        mean_sysTime = tot_sysTime/(tot_arrs-rejs)
    return mean_sysTime

def getRejections(a,x):
    """Takes an arrival vector and the schedule and returns the number of rejections. 

    Args:
        a (list): list of number of arrivals at the beginning of each time slot
        x (list): schedule, i.e., binary array of time slots with 1s in available slots

    Returns:
        int: number of rejections.
    """    
    pending_vec = [None]*len(a)
    for i in range(len(a)):
        if i==0:
            pending_vec[i]= max(a[i]-x[i],0)
        else: 
            pending_vec[i]= max( pending_vec[i-1]+ a[i]-x[i],0)
    return pending_vec[-1]
def getSysTime(a,x):
    """Takes an arrival vector and a schedule and returns the total system time. 

    Args:
        a (list): list of number of arrivals at the beginning of each time slot
        x (list): schedule, i.e., binary array of time slots with 1s in available slots

    Returns:
        float: total system time (sum accross all slots)
    """    
    a_lol=[]
    for timeslot in range(len(a)):
        if timeslot==0:
            appendage=[timeslot]*a[timeslot]
            appendage = appendage + [timeslot+1]
            a_lol.append(appendage)
        else:
            if x[timeslot-1]==0:
                appendage= appendage[:-1]+[timeslot]*a[timeslot]
                appendage = appendage + [timeslot+1]
                a_lol.append(appendage)
            elif x[timeslot-1]==1:
                appendage= appendage[1:-1]+[timeslot]*a[timeslot]
                appendage = appendage + [timeslot+1]
                a_lol.append(appendage)
    a_min=[]
    for ele in a_lol:
        a_min=a_min + [ele[0]]
    sysTime=0
    for i in range(len(a)):
        sysTime = sysTime + x[i]*(i-a_min[i]+1)
    return sysTime
def getArrV(lam):
    knots = dict()
    for i in range(len(lam)):
        knots[i] =  lam[i]
        knots[i+0.999999] = lam[i]
    knots[i+1]=0
    arrs = nhpp.get_arrivals(knots)
    a=[None]*(len(lam))
    # a[0]=0
    for i in range(len(lam)):
        a[i]= len(list(filter(lambda x: i <= x <= i+1, arrs)))
    return a

#%% INSTRUCTIONS TO RUN PSA
"""
To get PSA, run get_psaWT(lam_vec,serRt_vec), where lam_vec is the arrival rate
vector and serRT_vec is the service rate vector.    
   
"""


def PSA_MD1(lam_vec,serRt_vec):
    """Takes arrival rate vector and service rate vector and returns the PSA- based expected waiting time (not system time).

    Args:
        lam_vec (list): arrival rate vector (each element is stationary arrival rate at that time slot)
        serRt_vec (list): service rate vector (each element is service rate at that time slot calculated based on the schedule)

    Returns:
        float: Th PSA- based expected waiting time (NOT SYSTEM TIME) for the given arrival rate and service rate vectors.
    """    
    
    #empty utilization vector
    rho_vec = [None]* len(lam_vec)
    #update utilization vector
    rho_vec = [lam_vec[i]/serRt_vec[i] for i in range(len(lam_vec))]

    #empty PSA waiting time vector 
    waitTime_vec = [None]* len(lam_vec)
    #update PSA waiting time vector
    waitTime_vec = [lam_vec[i]*rho_vec[i]/(2*serRt_vec[i]*(1-rho_vec[i])) for i in range(len(lam_vec))]
    sysTtime_vec = [lam_vec[i]*((rho_vec[i]/(2*serRt_vec[i]*(1-rho_vec[i])))+(1/serRt_vec[i])) for i in range(len(lam_vec))]
    #total waiting time
    sum_WT = np.sum(waitTime_vec)
    sum_ST = np.sum(sysTtime_vec)
    # total expected number of arrivals
    sum_lam = np.sum(lam_vec)

    return sum_WT/sum_lam+1, sum_ST/sum_lam

#%% 
# READ ARRIVAL DATA
# Get the current working directory
current_directory = os.path.dirname(os.path.abspath(__file__))

# Go back one folder (.. means go up one directory)
parent_directory = os.path.join(current_directory, '..')

# Now, go to the 'data' folder
data_directory = os.path.join(parent_directory, 'data')

arrival_filename = "ArrivalRates.xlsx"

arrival_filepath = os.path.join(data_directory, arrival_filename)

arrival_data = pd.read_excel(arrival_filepath, index_col=0)

weekly_lam_r = arrival_data[arrival_data.columns[0]].tolist()

#convert weekly arrival rate vector to hourly arrival rate vector
lambda_r=[]
for i in weekly_lam_r:
    lambda_r=lambda_r+[i]*40
lambda_vec = lambda_r

#%% 
# TO GENERATE SCHEDULES

# Generate minimum schedule
MinSched = MinViableSched(lambda_vec)

# Get service rate array for minimum schedyle 
st, sr = get_ServTR_vecs(MinSched)

# Test the minimum schedule if it's viable
test = checksteadystate(lambda_vec, sr)

U = 720 #Specify the upper bound of the allocated number of slots
n = 500 # Specify the total number of schedules desired
Schedules = generate_scheds(MinSched, U, n)



#%% 
# Get PSA and MCS for the generated schedules

PSA_MD1_WT = [None]*len(Schedules)
PSA_MD1_ST = [None]*len(Schedules)
MCS_ST = [None]*len(Schedules)

for i in range(len(Schedules)):
    st, sr = get_ServTR_vecs(Schedules[i])
    PSA_MD1_WT[i], PSA_MD1_ST[i] = PSA_MD1(lambda_vec, sr)
    
    mcs_sims = 1000
    # total number of simulations (schedules)
    sims = 50
    #current sim initialize
    sim = 0 
    # empty sys times and rejection vectors      
    systimes=[None]*mcs_sims
    rejs=[None]*mcs_sims
    for mcs_sim in range(mcs_sims):
        try: 
            systimes[mcs_sim] = SimSystem(lambda_vec,Schedules[i])
        except Exception as e: 
            print("sim->",sim)
            print(str(e))
    MCS_ST[i] = np.mean(systimes)

        
    print(f'Schedule {i} is done.')
    



#%% 
# Generate an excel file with all results

sims = list(range(1,n+1))
full_data = {'Simulation': sims, 'MCS': MCS_ST, 'PSA': PSA_MD1_WT}
df_full = pd.DataFrame(full_data)
# Generate current timestamp for the filename
timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
df_full.to_excel(f'Full_Results_MCS_PSA_{timestamp}.xlsx', index=False)

# %%
