#%% 
# 1. IMPORT LIBRARIES
from gurobipy import Model, GRB, quicksum
import pandas as pd
from datetime import datetime
import os
#%% 
# 2. READ ARRIVAL DATA
# Get the current working directory
current_directory = os.path.dirname(os.path.abspath(__file__))

# Go back one folder (.. means go up one directory)
parent_directory = os.path.join(current_directory, '..')

# Now, go to the 'data' folder
data_directory = os.path.join(parent_directory, 'data')

arrival_filename = "ArrivalRates.xlsx"

arrival_filepath = os.path.join(data_directory, arrival_filename)

arrival_data = pd.read_excel(arrival_filepath, index_col=0)
d_1 = arrival_data[arrival_data.columns[0]].tolist()
d_2 = arrival_data[arrival_data.columns[1]].tolist()
#%% 
# 3. CREATE OPTIMIZATION MODEL

# Parameters 
W = 40 
U_1 = 58
U_2 = 29

# Create a new Gurobi model
m = Model("DPS Model")

# Create decision variables
y_1 = m.addVars(len(d_1), vtype=GRB.INTEGER, lb = 0, name="y1")
y_2 = m.addVars(len(d_2), vtype=GRB.INTEGER, lb = 0, name="y2")

# Set the objective to minimize the sum of squared differences from the targets
objective = quicksum(((d_1[i]) - (y_1[i]/W))*((d_1[i]) - (y_1[i]/W)) for i in range(len(d_1))) + quicksum(((d_2[i]) - (y_2[i]/W))*((d_2[i]) - (y_2[i]/W)) for i in range(len(d_2)))
m.setObjective(objective, GRB.MINIMIZE)

#Add constraints
m.addConstr(y_1.sum() == U_1, "UB_constraint1")
m.addConstr(y_2.sum() == U_2, "UB_constraint2")
for i in range(len(d_2)):
    m.addConstr(y_1[i]+y_2[i] <= W, "Max slots Constraint")
#Disable console output
m.setParam('OutputFlag', 0)

# Optimize model
m.optimize()
#%% 
# 4. PRINT RESULTS
first_time=[]
crisis=[]
# Print results for x variables
print("Results for  variables:") 
print('First-time')
i=1
y_1sum =0
for v in y_1.values():
    print(f'Week {i}: {v.x}')
    first_time=first_time + [v.x]
    i+=1
    y_1sum+= v.x
print(f'R_sum: {y_1sum}')
print('---------------------------')

print('Crisis')
i=1
y_2sum =0
for v in y_2.values():
    print(f'Week {i}: {v.x}')
    crisis=crisis + [v.x]
    i+=1
    y_2sum+= v.x
print(f'C_sum: {y_2sum}')
print('---------------------------')


print(f'Obj: {m.objVal}')


# %%
# Create a DataFrame
df = pd.DataFrame({
    'first_time prop': [x / 40 for x in first_time],
    'crisis prop': [x / 40 for x in crisis]
})

# Set the index as Week, starting from 1
df.index = range(1, len(df) + 1)

# Generate current timestamp for the filename
timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')

# Save to Excel file
df.to_excel(f'DPS-Policy_{timestamp}.xlsx', index=True)

# %%
