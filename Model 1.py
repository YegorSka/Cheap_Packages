import csv
import networkx as nx
import matplotlib.pyplot as plt
import copy

class transportation_system:

    def __init__(self,size = 0,flow = [[]],cost = [[]], hubs = set(), mult = [], bc = []): #initialize graphs its their data
        self.size = size
        self.flow = [[0 for j in range(size)] for i in range(size)]
        self.cost = [[0 for j in range(size)] for i in range(size)]
        self.connections = [[1 for j in range(size)] for i in range(size)]
        self.hubs = hubs
        self.mult = mult
        self.bc = bc

    def read(self,flowname,costname,multname,buildcostname):
        """"Builds a matrix graph from the entered flow,cost,multiplier and buildcost data"""
        with open(flowname, newline='') as file:
            reader = csv.reader(file, delimiter=';')
            res = [list(row) for row in reader]
            res.pop(0)
            for i in range(len(res)):
                res[i].pop(0)
        self.flow = [[float(i) for i in lst] for lst in res]
        self.size = len(res)

        with open(costname, newline='') as file:
            reader = csv.reader(file, delimiter=';')
            res = [list(row) for row in reader]
            res.pop(0)
            for i in range(len(res)):
                res[i].pop(0)
        self.cost = [[float(i) for i in lst] for lst in res]

        with open(buildcostname, newline='') as file:
            reader = csv.reader(file, delimiter=';')
            tempres = [row for row in reader]
            for i in range(len(res)):
                tempres[i].pop(0)
            res = [cost[0] for cost in tempres]
        self.bc =  [float(i) for i in res]

        with open(multname, newline='') as file:
            reader = csv.reader(file, delimiter=';')
            tempres = [row for row in reader]
            res = []
            for i in range(3):
                res.append(float(tempres[i][1]))
        self.mult = res

    def betweenness(self):
        """"Calculates the betweenness centrality of each node in the graph"""
        betweenness = {node:0 for node in range(self.size)} #initialize centrality of all nodes: 0
        for i in range(self.size): #compare all paths to find the minimal one between every pair (i,j) of nodes in graph
            for j in range(self.size):
                smallest_dist = (float("inf"))
                small = []
                for a in range(self.size):
                    if a != i and a != j: #Prevents paths where a would go to itself
                        if (self.mult[0] * self.cost[i][a] + self.mult[2] * self.cost[a][j]) == smallest_dist: #scenario with two shortest paths
                            small.append(a)
                        if (self.mult[0]*self.cost[i][a] + self.mult[2]*self.cost[a][j]) < smallest_dist:
                            small = [a] #node with smallest distance gets put in a list called small
                            smallest_dist = self.mult[0]*self.cost[i][a] + self.mult[2]*self.cost[a][j] #updates smallest distance
                for p in small:
                    betweenness[p] += self.flow[i][j] * (1/len(small))
        return betweenness
        
    def determine_alpha(self):
        build_cost = 0
        route_cost = 0
        flow = 0
        for i in range(self.size):
            for j in range(self.size):
                flow += self.flow[i][j]
                route_cost += self.cost[i][j]
        for k in range(self.size):
            build_cost += self.bc[k]
        route_cost *= 0.5
        alpha = build_cost / (flow * route_cost)
        return alpha

    def absolute_centrality(self,alpha):
        """"Every node now has a betweenness centrality.
            In order to get a better centrality (absolute centrality)
            the betweenness centrality of every node needs to be subtracted
            with the cost of making the node a hub (bc[node]).
            This new centrality applied to every node gets returned.
        """
        A = self.betweenness()
        res = []
        for node in range(self.size):
            res.append((node,A[node]-alpha*self.bc[node]))
        sorted_list = sorted(res, key=lambda x: x[1])
        return sorted_list

    def system_cost(self):
        cost = 0
        for hub in self.hubs: #add building costs to total cost
            cost+=self.bc[hub]
        for A in range(self.size): #calculate transport cost from A to B
            for B in range(self.size):
                routecost = 0
                if A not in self.hubs: #is A a hub? If yes then step 1 is free so we skip this step, else find it's hub.
                    h1 = 0
                    for i in range(self.size):
                        if self.connections[A][i] == 1: #found hub connected to A
                            routecost += self.mult[0]*self.flow[A][B]*self.cost[A][i]
                            h1=i #save hub 1
                            break
                if A in self.hubs:
                    h1 = A
                if B not in self.hubs: #is B a hub? If yes then step 2 is free so we skip this step, else find it's hub.
                    h2 = 0
                    for j in range(self.size):
                        if self.connections[B][j] == 1: #found hub connected to B
                            routecost += self.mult[2]*self.flow[A][B]*self.cost[B][j]
                            h2 = j #save hub 2
                            break
                if B in self.hubs:
                    h2 = B
                routecost+= self.mult[1]*self.flow[A][B]*self.cost[h1][h2]
                cost+=routecost
        return cost

    def determine_d_min(self): #Determines median distance between two nodes
        """Determines the minimum distance between hubs."""
        lst = []
        for i in range(self.size):  # makes a sorted list with all route costs
            for j in range(self.size):
                if i < j:
                    lst.append(self.cost[i][j])
        lst.sort()
        a = int(len(lst) * 0.5)
        return lst[a]

    def determine_2_hubs(self, d_min, alpha):
        """Determines which combination of 2 hubs is the best,
        based on absolute centrality. There must be a minimum distance
        between the two nodes."""
        A = self.absolute_centrality(alpha)
        A.reverse()
        for i in range(1, self.size):
            for j in range(i):
                if self.cost[A[i][0]][A[j][0]] >= d_min:
                    return {A[i][0], A[j][0]}

    def determine_3_hubs(self, d_min, alpha):
        """Determines which combination of 3 hubs is the best,
        based on absolute centrality. There must be a minimum distance
        between every pair of two hubs."""
        A = self.absolute_centrality(alpha)
        A.reverse()
        for i in range(2, self.size):
            for j in range(1, i):
                if self.cost[A[i][0]][A[j][0]] >= d_min:
                    for k in range(j):
                        if self.cost[A[k][0]][A[j][0]] >= d_min and self.cost[A[i][0]][A[k][0]] >= d_min:
                            return {A[i][0], A[j][0], A[k][0]}

    def determine_4_hubs(self, d_min, alpha):
        """Determines which combination of 4 hubs is the best,
        based on absolute centrality. There must be a minimum distance
        between every pair of two hubs."""
        A = self.absolute_centrality(alpha)
        A.reverse()
        for i in range(3, self.size):
            for j in range(2, i):
                if self.cost[A[i][0]][A[j][0]] >= d_min:
                    for k in range(1, j):
                        if self.cost[A[k][0]][A[j][0]] >= d_min and self.cost[A[i][0]][A[k][0]] >= d_min:
                            for l in range(k):
                                if self.cost[A[l][0]][A[j][0]] >= d_min and self.cost[A[i][0]][A[l][0]] >= d_min and \
                                        self.cost[A[k][0]][A[l][0]] >= d_min:
                                    return {A[i][0], A[j][0], A[k][0], A[l][0]}

    def determine_5_hubs(self, d_min, alpha):
        """Determines which combination of 5 hubs is the best,
        based on absolute centrality. There must be a minimum distance
        between every pair of two hubs."""
        A = self.absolute_centrality(alpha)
        A.reverse()
        for i in range(4, self.size):
            for j in range(3, i):
                if self.cost[A[i][0]][A[j][0]] >= d_min:
                    for k in range(2, j):
                        if self.cost[A[k][0]][A[j][0]] >= d_min and self.cost[A[i][0]][A[k][0]] >= d_min:
                            for l in range(1, k):
                                if self.cost[A[l][0]][A[j][0]] >= d_min and self.cost[A[i][0]][A[l][0]] >= d_min and \
                                        self.cost[A[k][0]][A[l][0]] >= d_min:
                                    for m in range(l):
                                        if self.cost[A[m][0]][A[j][0]] >= d_min and self.cost[A[i][0]][
                                            A[m][0]] >= d_min and self.cost[A[k][0]][A[m][0]] >= d_min and \
                                                self.cost[A[l][0]][A[m][0]] >= d_min:
                                            return {A[i][0], A[j][0], A[k][0], A[l][0], A[m][0]}

    def solution(self):
        alpha = self.determine_alpha()
        centrality = self.absolute_centrality(alpha)
        self.connections = self.build({centrality[-1][0]})
        self.hubs = {centrality[-1][0]}
        H = {centrality[-1][0]}
        min_cost = self.system_cost()
        d_min = self.determine_d_min()
        if self.size >= 2:
            h2 = self.determine_2_hubs(d_min, alpha)
            if h2 != None:
                self.connections = self.build(h2)
                self.hubs = h2
                if self.system_cost() < min_cost:
                    H = h2
                    min_cost = copy.copy(self.system_cost())
            if self.size >= 3:
                h3 = self.determine_3_hubs(d_min, alpha)
                if h3 != None:
                    self.connections = self.build(h3)
                    self.hubs = h3
                    if self.system_cost() < min_cost:
                        H = h3
                        min_cost = copy.copy(self.system_cost())
                if self.size >= 4:
                    h4 = self.determine_4_hubs(d_min, alpha)
                    if h4 != None:
                        self.connections = self.build(h4)
                        self.hubs = h4
                        if self.system_cost() < min_cost:
                            H = h4
                            min_cost = copy.copy(self.system_cost())
                    if self.size >= 5:
                        h5 = self.determine_5_hubs(d_min, alpha)
                        if h5 != None:
                            self.connections = self.build(h5)
                            self.hubs = h5
                            if self.system_cost() < min_cost:
                                H = h5
                                min_cost = copy.copy(self.system_cost())
        self.connections = self.build(H)
        self.hubs = H
        return min_cost

    def build(self,hubs):
        """"makes matrix with 1 if A is directly connected with B and one of them is a hub and 0 if A and B are not directly connected"""
        M = [[0 for j in range(self.size)] for i in range(self.size)] #makes matrix only zero's
        for i in range(self.size): #checks for every node is it a hub?
            if i in hubs:
                for hub in hubs: #add an 1 to the matrix for every contection between the hub and a other hub
                    M[i][hub] = 1
                    M[hub][i] = 1
            else:
                c = float("inf")
                hub = 0
                for k in range(self.size):
                    if k in hubs and self.cost[i][k] < c:
                        c = self.cost[i][k]
                        hub = k
                M[hub][i] = 1
                M[i][hub] = 1
        return M

    def show(self):
        """Creates a picture of the current graph"""
        num_nodes = self.size
        G = nx.Graph()
        for i in range(num_nodes):
            for j in range(i + 1, num_nodes):  # Only iterate over upper triangular part
                if self.connections[i][j] != 0:
                    G.add_edge(i+1, j+1, weight=self.cost[i][j])
        # Determine nodes with more than one edge
        nodes_to_color = ['darkorange' if G.degree(node) > 1 else 'skyblue' for node in G.nodes()]
        # Determine node sizes
        node_sizes = [800 if G.degree(node) > 1 else 400 for node in G.nodes()]
        # Draw the graph with a title and variable node sizes
        plt.figure()
        hub_nodes = [node for node in G.nodes() if G.degree(node) > 1]
        # Determine edge widths
        edge_widths = [2.5 if (edge[0] in hub_nodes and edge[1] in hub_nodes) else 1 for edge in G.edges()]
        edge_colors = ['red' if (edge[0] in hub_nodes and edge[1] in hub_nodes) else 'black' for edge in G.edges()]
        pos = nx.spring_layout(G, weight='weight')
        nx.draw(G, with_labels=True, node_color=nodes_to_color, node_size=node_sizes, font_size=12, font_weight='bold', width=edge_widths, edge_color=edge_colors)
        legend_text = "System Cost Information:\n\n" + "142707.0"
        plt.text(0.05, 0.05, legend_text, transform=plt.gca().transAxes, fontsize=10, verticalalignment='bottom')
        plt.show()

A = transportation_system()
A.read("w_small.csv","c_small.csv","g_small.csv","f_small.csv")
A.hubs = [1,5]
A.connections = A.build([1,5])
A.show()
print(A.system_cost())
# print(A.solution(0.5))
# A.show()
#
B = transportation_system()
B.read("w_large.csv","c_large.csv","g_large.csv","f_large.csv")
B.hubs = [1,9,13]
B.connections = B.build([1,9,13])
B.show()


