#%% 
# Import Libraries
import matplotlib.pyplot as plt
import random
from tqdm import tqdm
import numpy as np 
import pandas as pd
import gurobipy as gp
from gurobipy import *
from gurobipy import tupledict
import time
import cmath
import os
import math
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


lambda_vec=[]
for i in weekly_lam_r:
    lambda_vec=lambda_vec+[i]*40
    
M=100000


            
#%% 
# Helper functions for Gurobi solution
# Linearized Formulation
def minWT_new(lambda_vec,min_slots):
      M=9999999
      # Create the model within the Gurobi environment
      model = gp.Model("Test")
      model.params.NonConvex = 2
      model.params.MIPFocus = 3
      model.setParam('TimeLimit', 120*60)
      
      
      s = model.addMVar(len(lambda_vec), vtype= GRB.BINARY,name="s")
      h = model.addMVar(len(lambda_vec),lb=0, vtype= GRB.CONTINUOUS,name="h")
      z = model.addMVar(len(lambda_vec),lb=0, vtype= GRB.CONTINUOUS,name="z")
      zz = model.addMVar(len(lambda_vec),lb=0, vtype= GRB.CONTINUOUS,name="zz")

      model.setObjective(sum(zz[i] for i in range(len(lambda_vec))), GRB.MINIMIZE)




      #BOUNDARY CONDTIONS ON X AND H
      model.addConstr(h[-1]==1)
      model.addConstr(s[-1]==1)

      #NON LINEAR CONSTRAINTS
      # model.addConstrs(h[i]==s[i]+(1-s[i])*(1+h[i+1]) for i in range(len(lambda_vec)-1))
      
      
      # #LINEARIZING CONSTRAINTS
      model.addConstrs(h[i]==h[i+1]+1-z[i] for i in range(len(lambda_vec)-1))
      model.addConstrs(z[i]<=s[i]*M for i in range(len(lambda_vec)))
      model.addConstrs(z[i]<=h[i+1] for i in range(len(lambda_vec)-1))
      model.addConstrs(z[i]>=h[i+1]-M*(1-s[i]) for i in range(len(lambda_vec)-1))


      #UNDERUTILIZATION CONSTRAINT
      model.addConstrs(h[i]<=1/lambda_vec[i] for i in range(len(lambda_vec)))

      #MINPROP CONSTRAINT
      model.addConstr(sum(s)<=min_slots)
      
      #REMOVE DIVIDE TERM
      model.addConstrs(zz[i]*2*(1-(lambda_vec[i]*h[i]))>=(lambda_vec[i]**2)*(h[i]**2) for i in range(len(lambda_vec))) 

      #Solve
      model.update()
      #Disable console output
      model.setParam('OutputFlag', 0)
      model.optimize()

      return s.x, model.objVal    

def convert_list(sol):
    """converts the binary schedule into the exact slots in the shortest path

    Parameters
    ----------
    sol : list
        A binary list with 1 representing nodes that fall on the shortest path

    Returns
    -------
    list
        A list representing the nodes in the shortest path
    """
    list_sol = []
    for ele in range(len(sol)):
        if sol[ele]==1:
            list_sol.append(ele+1)
    list_sol=list_sol[:-1]
    list_sol=list_sol+[0]
    sorted_sol=list(np.sort(list_sol))
    sorted_sol = sorted_sol[::-1]
    sorted_sol = [int(ele) for ele in sorted_sol]
    return sorted_sol   
#%% 
# Helper functions for ssist
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


def compute_SSiSt(lambda_r,regular_slots):
    """solves ssist and returns solution and objective function value

    Parameters
    ----------
    lambda_r : list
        arrival rate vector
    regular_slots : int
        upper bound

    Returns
    -------
    list, float
        the shortest path node list, the optimal objective function value
    """
    feas_r=feas_limit(lambda_r)
    wt_reg=get_wts(lambda_r,feas_r)
    reg_sol, reg_of =Bellman(wt_reg,regular_slots)
    return reg_sol, reg_of


#%% 
# MAIN FUNCTION

if __name__ == '__main__':
    # set the planning horizon length and upper bound
    planning_horizon = 200
    Ub = 50
    slots = planning_horizon
    pop = slots
    #start time
    start_time = time.time()
    
    #slice the arrival vector
    lambda_vec_=lambda_vec[:planning_horizon]

    #solve using solver
    sol, sol_of = minWT_new(lambda_vec_,Ub) 

    #solver solution time 
    print("Gurobi time:--- %s mins ---" % ((time.time() - start_time)/60))
    
    #convert solution from binary to shortest path node list
    list_sol = convert_list(sol) 
    
    #reset start time
    start_time = time.time()

    feas_r=feas_limit(lambda_vec_)
    wt_reg=get_wts(lambda_vec_,feas_r)
    #ssist solution time 
    print("Ssist time: --- %s mins ---" % ((time.time() - start_time)/60))
    

    solver_of = calc_meas(wt_reg,list_sol)
    ssist_sol, ssist_of = compute_SSiSt(lambda_vec_,Ub)

    #Final solutions and OF values
    print(f'Final solver Solution: {list_sol}')
    print(f'Final ssist Solution: {ssist_sol}')
    print(f'Final solver OF: {solver_of}')
    print(f'Final ssist OF: {ssist_of}')

