# -*- coding: utf-8 -*-
"""
Created on April '20'

@author: Sohom Chatterjee
"""
#%% 
# 1. Import Libraries
import matplotlib.pyplot as plt
import random
import os
from tqdm import tqdm
import numpy as np 
import pandas as pd
import math
#%% 
# 2. Function Definitions
def weight(lam,Node_1,Node_2):
    """Takes the arrival vector, the starting and 
        ending node and returns the edge weight for Problem S-SiST, using Eq. (7) of the manuscript

    Parameters
    ----------
    lam : list
        arrival rate vector
    Node_1 : int
        starting node
    Node_2 : int
        ending node

    Returns
    -------
    float
        the edge weight of the edge starting at Node 1 and ending at Node 2
    """
    sum=0
    for i in range(Node_1,Node_2):
        sum+=(lam[i]**2)*((Node_2-i)**2)/(1-(lam[i]*(Node_2-i)))  
    return sum


def feas_limit(lv):
    """This function takes as input the arrival vector and returns final_vec. 
        Each item in final_vec represents the node till which we construct edges for that index in the graph.
        eg., final_vec[i]=j represents that in the graph, each node from i to j-1 will have a node to j.

    Parameters
    ----------
    lambda_vec : list

    Returns
    -------
    final_vec : list
    """
    slots = len(lv)
    slot_vec=[i+math.floor(1/lv[i]) if lv[i]!=0 else 1000 for i in range(slots)]
    low_lim=0
    final_vec=[]
    while len(final_vec)<slots:
        mn=min(slot_vec[low_lim:])
        argmin=(slot_vec[low_lim:]).index(mn)+len(slot_vec[:low_lim])
        hm=argmin-low_lim+1
        final_vec=final_vec+[mn]*hm  
        low_lim=low_lim+hm
    return final_vec[0:slots]


def Bellman(cost,B):   
    """takes in an edge weight matrix and a value for B (maximum number of edges) 
    and returns the shortest path having at most B edges, and the optimal objective function value

    Parameters
    ----------
    cost : list of lists
        cost[i][j] is the edge weight from node i to node j
    B : int
        maximum number of edges for the shortest path

    Returns
    -------
    list, float
        node_list: a list representing the nodes in the shortest path, of: the optimal objective function value
    """
    d=[[None for i in range(pop+1)] for i in range(B+1)]
    d[0]=[float("inf") for i in range(pop+1)]
    d[0][0]=0.0
    pathNode=[[None for i in range(pop+1)] for i in range(B+1)]
    for k in range(B):
        for i in range(pop+1):
            tempMin=float("inf")
            tempIndex=-1.0
            for j in range(i):
                if d[k][j]+cost[j][i]<tempMin:
                    tempMin=d[k][j]+cost[j][i]
                    tempIndex=j
            if tempMin<d[k][i]:
                d[k+1][i]=tempMin
                pathNode[k+1][i]=tempIndex
            else:
                d[k+1][i]=d[k][i]
                pathNode[k+1][i]=pathNode[k][i]
    
    path=[]
    node=pop
    k=B
    while node>0:
        node=pathNode[k][node]
        path.append(node)
        k=k-1
        
    node_list=path
    of=0
    of+=cost[node_list[0]][slots]
    for ele in range(0,len(node_list)-1):
        of+=cost[node_list[ele+1]][node_list[ele]]
    return node_list,of

def calc_meas(cost,node_list):
    """takes in a edge weight matrix and a path and returns the total length of the path (sum of edge weights in the path)

    Parameters
    ----------
    cost : list of lists
        cost[i][j] is the edge weight from node i to node j
    node_list : list
        a list representing the nodes in the shortest path

    Returns
    -------
    float
        total path length (sum of edge weights in the path)
    """
    of=0
    of+=cost[node_list[0]][slots]
    for ele in range(0,len(node_list)-1):
        of+=cost[node_list[ele+1]][node_list[ele]]
    return of


def get_wts(lam,feas):
    """Constructs the edge weight matrix using the arrival rate vector and the feas vec 

    Parameters
    ----------
    lam : list
        arrival rate vector
    feas : list
        Each item in final_vec represents the node till which we construct edges for that index in the graph.
        eg., final_vec[i]=j represents that in the graph, each node from i to j-1 will have a node to j. This ensures that 
        the queues are stable.

    Returns
    -------
    list of lists
        the edge weight matrix
    """
    wt=[[100000000 for i in range(slots+1)]for i in range(slots)]
    for i in range(0,slots):
       for j in range(i+1,slots+1):
           if j>=feas[i]:
               break
           else:
               wt[i][j]=weight(lam,i,j)
    return wt

#%% 
# 3. BnB Algorithm

def chk_conf(a,b):
    """returns list of conflicting time slot

    Parameters
    ----------
    a : list
        list of nodes in first path
    b : list
        list of nodes in second path

    Returns
    -------
    list
        list of nodes common in both paths (conflicts in corresponding schedules)
    """
    # Remove last element (0) 
    a_set = set(a[:-1])
    b_set = set(b[:-1])
    if (a_set & b_set):
        return list(a_set & b_set)
    else:
        return []
    
def modifygraph(wt,conf):
    """modifies the edge weight matrix by setting a high value for all edges coming out of conflict node

    Parameters
    ----------
    wt : list of lists
        edge weight matrix
    conf : int
        conflict node

    Returns
    -------
    list of lits
        modified edge weight matrix such that every edge from conf node has very high weight (preventing them from being included in shortest path)
    """
    for i in range(len(wt)):
        wt[i][conf]=100000000
    return wt    
    


def Branch_Bound(wr,fr,wc,fc):
    """Recursive function to generate every node of the Branch and bound tree

    Parameters
    ----------
    wr : list of lists
        edge weight matrix for first-time
    fr : list
        feasibility vector for first-time
    wc : list of lists
        edge weight matrix for crisis
    fc : _type_
        feasibility vector for crisis
    """
    # Solve Shortest Path
    reg_sol=Bellman(wr,regular_slots)
    cris_sol=Bellman(wc,crisis_slots)
    # Check Conflicts
    conflict_list=chk_conf(reg_sol[0],cris_sol[0])

    # if no conflicts, add solution and OF to a global list and exit
    if conflict_list==[]:
        global solution_list
        solution_list=solution_list+[reg_sol,cris_sol]
        return
    # if conflicts, resolve conflicts, 
    else: 
        # solve P1 Disallow 1st regular conflict
        conf=conflict_list[0]
        
        feas_r_up=[i if i!=conf else i-1 for i in fr ]
        wt_reg_up=modifygraph(wr,conf)
        Branch_Bound(wt_reg_up,feas_r_up,wt_crisis,feas_c)

        # solve P2 Disallow 1st crisis conflict
        feas_c_up=[i if i!=conf else i-1 for i in fc ]
        wt_crisis_up=modifygraph(wc,conf)
        Branch_Bound(wt_reg,feas_r,wt_crisis_up,feas_c_up)


#%% 6. Set initial parameters

# Get the current working directory
current_directory = os.path.dirname(os.path.abspath(__file__))

# Go back one folder (.. means go up one directory)
parent_directory = os.path.join(current_directory, '..')

# Now, go to the 'data' folder
data_directory = os.path.join(parent_directory, 'data')

arrival_filename = "ArrivalRates.xlsx"

arrival_filepath = os.path.join(data_directory, arrival_filename)
print(f'{arrival_filepath}')

arrival_data = pd.read_excel(arrival_filepath, index_col=0)
weekly_lam_r = arrival_data[arrival_data.columns[0]].tolist()
weekly_lam_c = arrival_data[arrival_data.columns[1]].tolist()


#Crisis Arrival Rate
# weekly_lam_c=[0.0111111111111111,0.0118055555555556,0.0201388888888889,0.0125,0.0208333333333333,0.0229166666666667,0.0201388888888889,0.0194444444444444,0.0229166666666667,0.0131944444444444,0.0215277777777778,0.0166666666666667,0.0125,0.00138888888888889,0.000694444444444444,0.000694444444444444,0,0]
#convert weekly arrival rate vector to hourly arrival rate vector
lambda_c=[]
for i in weekly_lam_c:
    lambda_c=lambda_c+[i]*40
    
#First-time Arrival Rate
# weekly_lam_r=[0.0823295454463125,0.103177897995375,0.081846910096875,0.0801925233928875,0.084468831868575,0.078303409110675,0.0824852064396937,0.0711768749619375,0.0667458333100125,0.06103125001575,0.0570360668810625,0.0245395161167625,0.0070811831538,0.0170455645206,0.00590625,0.00039375,0.00039375,0.00118125]
#convert weekly arrival rate vector to hourly arrival rate vector
lambda_r=[]
for i in weekly_lam_r:
    lambda_r=lambda_r+[i]*40


slots=len(lambda_r)
pop=slots

#%%
        
# a= Branch_Bound(wt_reg,feas_r,wt_crisis,feas_c)
def getPolicy(reg, cris, change_freq):
    """Takes the shortest path solution for first-time (reg), crisis (cris), and the frequency of change and converts the schedule into a policy

    Parameters
    ----------
    reg : list
        first-time solution 
    cris : list
        crisis solution
    change_freq : int
        the parameter used to calculate the change frequency in policy, weekly or monthly

    Returns
    -------
    pandas Series
        input for dataframe having the policies
    """
    
    
    policy = pd.DataFrame(columns = ["StartDay","Reg_Prop","Cris_Prop"])
    start_days = [i for i in range(0,18*change_freq,change_freq)]
    policy["StartDay"]=pd.Series(start_days)

    
    reg = [i for i in reg if i != 0]
    cris = [i for i in cris if i != 0]

    
    end_day = 18*5*8
    start_days=start_days+[end_day]
    
    reg_policy=[]
    cris_policy=[]

    for i in range(len(start_days)-1):
        denom = (start_days[i+1]-start_days[i])
        reg_policy.append((len(list(x for x in reg if((x >= start_days[i])&(x < start_days[i+1])))))/denom)
        cris_policy.append((len(list(x for x in cris if((x >= start_days[i])&(x < start_days[i+1])))))/denom)

    
    policy['Reg_Prop']=pd.Series(reg_policy)
    policy['Cris_Prop']=pd.Series(cris_policy)

    
    
    return policy


#%% Final Simulation
# range of average proportion values for each service
crisis_props=[i for i in range(3,6)]
reg_props=[i for i in range(7,10)]

for rp in tqdm(reg_props):
    for cp in tqdm(crisis_props):
        #calculate upper bound for each service
        regular_slots=math.ceil(rp/100*slots)
        crisis_slots=math.ceil(cp/100*slots)
        

        solution_list=[]
        
        #calculate feas limit as a steady state constraint
        feas_r=feas_limit(lambda_r)
        feas_c=feas_limit(lambda_c)

        #calculate edge weight matrices
        wt_reg=get_wts(lambda_r,feas_r)
        wt_crisis=get_wts(lambda_c,feas_c)

        #recursive branch and bound
        Branch_Bound(wt_reg,feas_r,wt_crisis,feas_c)
        
        df=pd.DataFrame(solution_list)
        
        df.to_excel('Pf'+str(rp)+"_Cf"+str(cp)+".xlsx")
        reg_sol = df[0][0]+[719]
        cris_sol = df[0][1]+[719]

        change_freq = 5*8 #weekly policy
        policy = getPolicy(reg_sol, cris_sol, change_freq)    
        policy.to_excel('FinalPolicy_Pf'+str(rp)+"_Cf"+str(cp)+".xlsx")

  
        
        
