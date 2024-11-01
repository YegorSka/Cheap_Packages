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

P = 200 #Max capacity
q = 300 #Cost to connect new hub to node

#create LP problem
PostNL = LpProblem("PostNL", LpMinimize)

#introduce variables
x = LpVariable.dicts("HubConection", (range(n), range(n), range(n), range(n)), 0, 1)
y = LpVariable.dicts("Connection", (range(n), range(n)), 0, 1, LpBinary)
z = LpVariable.dicts("IsHub", range(n), 0, 1, LpBinary)

#add objective function
PostNL += lpSum([HubFixedCost[k]*z[k] for k in range(n)])\
          \
          +lpSum([(lpSum([y[i][k] for k in range(n)]))*q for i in range(n)])\
          \
          +lpSum([PackageFlow[i][j]*x[i][j][k][l]*(EdgeCost[i][k]*CollectionCoefficient +\
          EdgeCost[k][l]*TransferCoefficient + EdgeCost[l][j]*DistributionCoefficient) \
          for i in range(n) for j in range(n) for l in range(n) for k in range(n)])

#constraints
PostNL += lpSum([z[k] for k in range(n)]) >= 1 #At least one hub
# with: 62216.99999971
# without: 62217.00000014
for i in range(n): #100% of the packages should arrive at their destination
    for j in range(n):
        PostNL += lpSum([x[i][j][k][l] for k in range(n) for l in range(n)]) == 1

for i in range(n): #non hubs nodes cannot act as hubs
    for j in range(n):
        for k in range(n):
            for l in range(n):
                PostNL += x[i][j][k][l] <= z[l]
                PostNL += x[i][j][k][l] <= z[k]

for k in range(n):
    PostNL+= lpSum([x[i][j][k][l]*PackageFlow[i][j]+x[i][j][l][k]*PackageFlow[i][j] for i in range(n) for l in range(n) for j in range(n)])\
    <= P*z[k]

for i in range(n):
    for k in range(n):
        PostNL+=y[i][k] >= (lpSum([x[i][j][k][l]+x[j][i][l][k] for l in range(n) for j in range(n)])/(2*n^2))-z[i]
        PostNL+=y[i][k] <= (lpSum([x[i][j][k][l] + x[j][i][l][k] for l in range(n) for j in range(n)]) / (2 * n ^ 2)) - z[i] + 1 - 10**-10

# PostNL += lpSum(z>=1, "there should be at least one hub"
# PostNL += (
#     lpSum(x for k in nodes)
# )
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

    G = nx.Graph()
    nodes_to_color = ['skyblue' for i in range(n)]
    node_sizes = [400 for i in range(n)]
    for i in range(n):
        for j in range(i + 1, n):  # Only iterate over upper triangular part
            if value(y[i][j]) != 0 and i != j:
                G.add_edge(i+1, j+1, weight=cost[i][j], edge_color="black", edge_width = 0.5)
    for i in established_hubs:
        nodes_to_color[i] = 'darkorange'
        node_sizes[i] = 800
        for j in established_hubs:
            if i != j:
                G.add_edge(i+1, j+1, weight=cost[i][j], edge_color="red", edge_width=2.5)
    plt.figure()
    # Determine edge widths
    # edge_widths = [2.5 if (edge[0] in hub_nodes and edge[1] in hub_nodes) else 1 for edge in G.edges()]
    # edge_colors = ['red' if (edge[0] in hub_nodes and edge[1] in hub_nodes) else 'black' for edge in G.edges()]
    nx.draw(G, with_labels=True, node_color=nodes_to_color, node_size=node_sizes, font_size=12, font_weight='bold')
    legend_text = "System Cost Information:\n\n" + str(value(PostNL.objective))
    plt.text(0.05, 0.05, legend_text, transform=plt.gca().transAxes, fontsize=10, verticalalignment='bottom')
    plt.show()