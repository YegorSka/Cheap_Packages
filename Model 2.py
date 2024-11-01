from pulp import *
import csv
import networkx as nx
import matplotlib.pyplot as plt
import math

#excelbestand lezen
def read(flowname, costname, multname, buildcostname):
    """"leest excelbestand"""
    with open(flowname, newline='') as file:
        reader = csv.reader(file, delimiter=';')
        res = [list(row) for row in reader]
        res.pop(0)
        for i in range(len(res)):
            res[i].pop(0)
    flow = [[float(i) for i in lst] for lst in res]
    size = len(res)

    with open(costname, newline='') as file:
        reader = csv.reader(file, delimiter=';')
        res = [list(row) for row in reader]
        res.pop(0)
        for i in range(len(res)):
            res[i].pop(0)
    cost = [[float(i) for i in lst] for lst in res]

    with open(buildcostname, newline='') as file:
        reader = csv.reader(file, delimiter=';')
        tempres = [row for row in reader]
        for i in range(len(res)):
            tempres[i].pop(0)
        res = [cost[0] for cost in tempres]
    bc = [float(i) for i in res]

    with open(multname, newline='') as file:
        reader = csv.reader(file, delimiter=';')
        tempres = [row for row in reader]
        res = []
        for i in range(3):
            res.append(float(tempres[i][1]))
    mult = res

    return flow, cost, bc, mult, size

#input data
flow, cost, bc, mult, size = read("w_small.csv","c_small.csv","g_small.csv","f_small.csv")
PackageFlow = flow
EdgeCost = cost
HubFixedCost = bc
CollectionCoefficient = mult[0]
TransferCoefficient = mult[1]
DistributionCoefficient = mult[2]
n = size

#create LP problem
PostNL = LpProblem("PostNL", LpMinimize)

#introduce variables
x = LpVariable.dicts("HubConection", (range(n), range(n), range(n), range(n)), 0, 1, LpInteger)
z = LpVariable.dicts("IsHub", range(n), 0, 1, LpInteger)
y = LpVariable.dicts("Connection", (range(n), range(n)), 0, 1, LpInteger)

#add objective function
PostNL += lpSum([HubFixedCost[k]*z[k] for k in range(n)]) + lpSum([PackageFlow[i][j]*x[i][j][k][l]*(EdgeCost[i][k]*CollectionCoefficient + EdgeCost[k][l]*TransferCoefficient + EdgeCost[l][j]*DistributionCoefficient) for i in range(n) for j in range(n) for l in range(n) for k in range(n)])

#constraints
#PostNL += lpSum([z[k] for k in range(n)])>=1
for i in range(n):
    for j in range(n):
        PostNL += lpSum([x[i][j][k][l] for k in range(n) for l in range(n)])==1

for i in range(n):
    for j in range(n):
        for k in range(n):
            for l in range(n):
                PostNL += x[i][j][k][l]<=z[l]
                PostNL += x[i][j][k][l]<=z[k]

for i in range(n):
    for j in range(n):
        for k in range(n):
            for l in range(n):
                PostNL += x[i][j][k][l]<=y[i][k] #if i and k are not connected, then x[i][j][k][l]=0
                PostNL += x[i][j][k][l]<=y[j][l] #if j and l are not connected, then x[i][j][k][l]=0

for i in range(n): #every node has one hub
    PostNL += lpSum(y[i][k] for k in range(n))==1

#solve problem
PostNL.solve()
#print status
print("Status:", LpStatus[PostNL.status])
#print objective function value
print("Min cost:", value(PostNL.objective),"euros")
#print solution

if LpStatus[PostNL.status] == 'Optimal':
    established_hubs = []
    for k in range(n):
        if value(z[k]) == 1:
            established_hubs.append(k)
    print("List of established hubs:", established_hubs)

