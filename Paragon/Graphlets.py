# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 21:26:21 2021

@author: mkaanarici
"""

#import graphlet_descriptor as descript
import networkx as nx, pandas as pd,pickle,numpy as np
import timeit

      




class GraphletSearch():
    def __init__(self,network):
        self.network=network
        self.neighbors=dict()
        for node in self.network.nodes:
            self.neighbors[node]=set([i for i in self.network.neighbors(node) if str(i)!="nan"])
        self.GD_vector=dict(zip(list(self.network.nodes),[np.zeros(15) for i in range(self.network.number_of_nodes())]))
        self.Graphlets=dict(zip([f'Graphlets{g}' for g in range(9)],[{} for i in range(9)]))
        
        self.Graphlets_Scores=dict(zip([f'Graphlets{g}' for g in range(9)],[{} for i in range(9)]))
        
        self.self_Graphlets=dict(zip([f'Graphlets{g}' for g in range(9)],[{} for i in range(9)]))
        self.Orbits=dict(zip([f'Orbits{o}' for o in range(15)],[set() for i in range(15)]))
        self.propagated_node_list=set()
        self.selected_edge_list=list()
        self.nodes=list()
        
        
        self.created_network=nx.Graph()
    def add_nodes_to_orbits(self, orbit, nodes):
        self.Orbits[f'Orbits{orbit}'].update(nodes)
        for node in nodes:
            self.GD_vector[node][orbit]+=1
    
    def get_graphlet_score(self, edges):
        Scores=0
        for edge in edges:
            try:
                Scores+=self.network[edge[0]][edge[1]]["Score"]
                return round(Scores/len(edges),3)
            except:
                return 1
        
    def get_graphlet_degree_vectors(self):
        return self.GD_vector
        
    def get_created_network(self):
        return self.created_network
        
    def get_selected_edge_list(self):
        return list(self.created_network.edges)
        
    def get_propagated_node_list(self):
        return list(self.propagated_node_list)
    
    def get_graphlets(self):
        return self.Graphlets
    
    def get_graphlets_scores(self):
        return self.Graphlets_Scores
    
    def get_orbits(self):
        return self.Orbits
     
    def write_created_network(self,file_name):
        DF=nx.to_pandas_edgelist(self.created_network)
        DF.to_csv(f'{file_name}.sif',sep="\t",index=False)
        return True
    
    def save_pickle_graphlets(self,file_name):
        with open(f'{file_name}.pickle',"wb") as f:
            pickle.dump(self.Graphlets,f)

    
    def save_pickle_graphlets_scores(self,file_name):
        with open(f'{file_name}.pickle',"wb") as f:
            pickle.dump(self.Graphlets_Scores,f)
        

    def save_pickle_orbits(self,file_name):
        with open(f'{file_name}.pickle',"wb") as f:
            pickle.dump(self.Orbits,f)
         
        
    def find_graphlets(self, Node_list,derivated_graphlets=False):
#        start=timeit.default_timer()
        self.nodes=set([i for i in Node_list if i in self.network.nodes])
#        print (f'{len(self.nodes)} of {len(Node_list)} input nodes have been found in the given network' )
#        start=timeit.default_timer()
        
        self.find_key_graphlets()
#        stop=timeit.default_timer()
#        print ("linear",stop-start)        
        
        self.G1_derivated=[]
        for key in self.Graphlets["Graphlets1"].keys():
            key_s=set(key)
            if len(key_s-self.nodes)==1:
                self.G1_derivated.append(self.Graphlets["Graphlets1"][key])     
 
        self.G2_derivated=[]   
        for key in self.Graphlets["Graphlets2"].keys():
            key_s=set(key)
            if len(key_s-self.nodes)==1:
                self.G2_derivated.append(self.Graphlets["Graphlets2"][key])
                
#        print("the number of G0 :",len(self.Graphlets["Graphlets0"].values()))
#        print("the number of G1 :",len(self.Graphlets["Graphlets1"].values()))
#        print("the number of G2 :",len(self.Graphlets["Graphlets2"].values()))

        if derivated_graphlets:
            self._derivate_graphlets()
            
        self.propagated_node_list=set()
        for i in range(15):
            self.propagated_node_list.update(self.Orbits[f'Orbits{i}'])
            
        for i in range (9):
            for key in self.Graphlets[f'Graphlets{i}']:
                self.created_network.add_edges_from( self.Graphlets[f'Graphlets{i}'][key])

        return self.Graphlets,self.Orbits

    def _derivate_graphlets(self):
#        start=timeit.default_timer()  
        Nodes=list(self.nodes)
        for i in range(len(self.nodes)):
            node1=Nodes[i]
            #N1=self.neighbors[node1] # set(i for i in self.network.neighbors(node1) if str(i)!="nan") #     
            for k in range(i+1,len(self.nodes)):
                node2=Nodes[k] 
                if node1==node2: 
                    continue 
                #N2=self.neighbors[node2] # set(i for i in self.network.neighbors(node2) if str(i)!="nan") # 
                for l in range(k+1,len(self.nodes)):
                    node3=Nodes[l]
                    if node3==node2 or node3==node1:
                        continue
                    #N3=self.neighbors[node3] # set(i for i in self.network.neighbors(node3) if str(i)!="nan") #
                    if ((node1,node2) in self.network.edges) and ((node1,node3) in self.network.edges) and ((node2,node3) in self.network.edges):
                        key=tuple(set([node1,node2,node3]))
                        self.self_Graphlets["Graphlets2"][key]=[(node1,node2),(node2,node3),(node3,node1)]
#                        self.add_nodes_to_orbits(3,[node1,node2,node3])
                        
                    elif ((node1,node2) in self.network.edges) and ((node1,node3) in self.network.edges) and ((node2,node3) not in self.network.edges):
                        key=tuple(set([node1,node2,node3]))
                        self.self_Graphlets["Graphlets1"][key]=[(node2,node1),(node1,node3)]
#                        self.add_nodes_to_orbits(1,[node2,node3])
#                        self.add_nodes_to_orbits(2,[node1])
                        
                    elif ((node1,node2) in self.network.edges) and ((node1,node3) not in self.network.edges) and ((node2,node3) in self.network.edges):
                        key=tuple(set([node1,node2,node3]))
                        self.self_Graphlets["Graphlets1"][key]=[(node1,node2),(node2,node3)]
#                        self.add_nodes_to_orbits(1,[node1,node3])
#                        self.add_nodes_to_orbits(2,[node2])                       
                        
                    elif ((node1,node2) not in self.network.edges) and ((node1,node3) in self.network.edges) and ((node2,node3) in self.network.edges):
                        key=tuple(set([node1,node2,node3]))
                        self.self_Graphlets["Graphlets1"][key]=[(node1,node3),(node3,node2)]
#                        self.add_nodes_to_orbits(1,[node1,node2])
#                        self.add_nodes_to_orbits(2,[node3])
                        
                        
                        
                        
                        
        self.find_graphlet_from_G1 ()            
        self.find_graphlet_from_G2 ()
#        print("the number of self_Graphlets3 :", len( self.self_Graphlets["Graphlets3"]))
#        print("the number of self_Graphlets4 :", len( self.self_Graphlets["Graphlets4"]))
#        print("the number of self_Graphlets5 :", len( self.self_Graphlets["Graphlets5"]))
#        print("the number of self_Graphlets6 :", len( self.self_Graphlets["Graphlets6"]))
#        print("the number of self_Graphlets7 :", len( self.self_Graphlets["Graphlets7"]))  
#        print("the number of self_Graphlets8 :", len( self.self_Graphlets["Graphlets8"]))       
#        
#        print("the number of G3 :", len( self.Graphlets["Graphlets3"]))
#        print("the number of G4 :", len( self.Graphlets["Graphlets4"]))
#        print("the number of G5 :", len( self.Graphlets["Graphlets5"]))
#        print("the number of G6 :", len( self.Graphlets["Graphlets6"]))
#        print("the number of G7 :", len( self.Graphlets["Graphlets7"])) 
#        print("the number of G8 :", len( self.Graphlets["Graphlets8"])) 
        
        
#        stop=timeit.default_timer()
#        print ("linear derivation",stop-start) 
        return (self.Graphlets, self.Orbits) 
    
    
 


    def find_graphlet_from_G1(self):

        G1_derivated=list(self.Graphlets["Graphlets1"].values())
        for i in range(len(G1_derivated)):
            g1=G1_derivated[i]
            A1,B1,i1=(g1[0][0],g1[1][1],g1[0][1])
            nodes=set([A1,B1,i1])
            
            NeighborsA1=self.neighbors[A1] - nodes
            NeighborsB1=self.neighbors[B1] - nodes
            Neighborsi1=self.neighbors[i1] - nodes
                
            Founded=(NeighborsA1 & self.nodes) - (NeighborsB1 | Neighborsi1)
            for i in Founded:
                key=tuple(set([A1,B1,i1,i]))
                self.Graphlets["Graphlets3"][key]=[(i,A1),(A1,i1),(i1,B1)]
                self.Graphlets_Scores["Graphlets3"][key]=self.get_graphlet_score([(i,A1),(A1,i1),(i1,B1)])
                self.add_nodes_to_orbits(4,[B1,i])
                self.add_nodes_to_orbits(5,[A1,i1])
                
            Founded=(NeighborsB1 & self.nodes) - (NeighborsA1 | Neighborsi1)
            for i in Founded:
                key=tuple(set([A1,B1,i1,i]))
                self.Graphlets["Graphlets3"][key]=[(i,B1),(B1,i1),(i1,A1)]
                self.Graphlets_Scores["Graphlets3"][key]=self.get_graphlet_score([(i,B1),(B1,i1),(i1,A1)])
                self.add_nodes_to_orbits(4,[A1,i])
                self.add_nodes_to_orbits(5,[B1,i1])
                     
            Founded=(Neighborsi1 & self.nodes) - (NeighborsA1 | NeighborsB1)
            for i in Founded:
                key=tuple(set([A1,B1,i1,i]))
                self.Graphlets["Graphlets4"][key]=[(A1,i1),(B1,i1),(i,i1)]
                self.Graphlets_Scores["Graphlets4"][key]=self.get_graphlet_score([(A1,i1),(B1,i1),(i,i1)])

                self.add_nodes_to_orbits(6,[A1,B1,i])
                self.add_nodes_to_orbits(7,[i1])         

                
            Founded=(NeighborsA1 & Neighborsi1 & self.nodes) - (NeighborsB1)
            for i in Founded:
                key=tuple(set([A1,B1,i1,i]))
                self.Graphlets["Graphlets6"][key]=[(A1,i1),(i1,B1),(i,A1),(i,i1)]
                self.Graphlets_Scores["Graphlets6"][key]=self.get_graphlet_score([(A1,i1),(i1,B1),(i,A1),(i,i1)])

                self.add_nodes_to_orbits(9,[B1])
                self.add_nodes_to_orbits(10,[A1,i])
                self.add_nodes_to_orbits(11,[i1])
                
            Founded=(NeighborsB1 & Neighborsi1 & self.nodes) - (NeighborsA1)
            for i in Founded:
                key=tuple(set([A1,B1,i1,i]))
                self.Graphlets["Graphlets6"][key]=[(B1,i1),(i1,A1),(i,B1),(i,i1)]
                self.Graphlets_Scores["Graphlets6"][key]=self.get_graphlet_score([(B1,i1),(i1,A1),(i,B1),(i,i1)])

                self.add_nodes_to_orbits(9,[A1])
                self.add_nodes_to_orbits(10,[B1,i])
                self.add_nodes_to_orbits(11,[i1])
                
                # self.Graphlets
        for graphlet in list(self.self_Graphlets["Graphlets1"].values()):
            A,B,C=graphlet[0][0],graphlet[0][1],graphlet[1][1]
            NeighborsA=self.neighbors[A] 
            NeighborsB=self.neighbors[B]
            NeighborsC=self.neighbors[C]
            
            Founded=(NeighborsA & NeighborsC) - (self.nodes |  NeighborsB)
            for each in Founded:
                key=tuple(set([A,B,C,each]))
                self.Graphlets["Graphlets5"][key]=[(A,B),(B,C),(C,each),(each,A)]
                self.Graphlets_Scores["Graphlets5"][key]=self.get_graphlet_score([(A,B),(B,C),(C,each),(each,A)])

                self.add_nodes_to_orbits(8,[A,B,C,each])
                
            Founded=(NeighborsA & NeighborsC & NeighborsB) - self.nodes
            for each in Founded:
                key=tuple(set([A,B,C,each]))
                self.Graphlets["Graphlets7"][key]=[(A,B),(B,C),(A,each),(B,each),(C,each)]
                self.Graphlets_Scores["Graphlets7"][key]=self.get_graphlet_score([(A,B),(B,C),(A,each),(B,each),(C,each)])
                self.add_nodes_to_orbits(12,[A,C])
                self.add_nodes_to_orbits(13,[B,each])
                
        self_Graphlet_1=list(self.self_Graphlets["Graphlets1"].values())
        for i in range(len(self_Graphlet_1)):
            g1=self_Graphlet_1[i]
            A1,B1,C1=(g1[0][0],g1[1][1],g1[0][1]) # A-->C-->B
            nodes=set([A1,B1,C1])
            
            NeighborsA1=(self.neighbors[A1] & self.nodes) - nodes
            NeighborsB1=(self.neighbors[B1] & self.nodes) - nodes
            NeighborsC1=(self.neighbors[C1] & self.nodes) - nodes
                
            Founded=(NeighborsA1) - (NeighborsB1 | NeighborsC1)
            for i in Founded:
                key=tuple(set([A1,B1,C1,i]))
#                self.self_Graphlets["Graphlets3"][key]=[(i,A1),(A1,C1),(C1,B1)]
                        
#                self.add_nodes_to_orbits(4,[B1,i])
#                self.add_nodes_to_orbits(5,[A1,C1])
                
            Founded=(NeighborsB1) - (NeighborsA1 | NeighborsC1)
            for i in Founded:
                key=tuple(set([A1,B1,C1,i]))
#                self.self_Graphlets["Graphlets3"][key]=[(i,B1),(B1,C1),(C1,A1)]
#                self.add_nodes_to_orbits(4,[A1,i])
#                self.add_nodes_to_orbits(5,[B1,C1])
                     
            Founded=(NeighborsC1) - (NeighborsA1 | NeighborsB1)
            for i in Founded:
                key=tuple(set([A1,B1,C1,i]))
#                self.self_Graphlets["Graphlets4"][key]=[(A1,C1),(C1,B1),(C1,i)]
#                self.add_nodes_to_orbits(6,[A1,B1,i])
#                self.add_nodes_to_orbits(7,[C1])
           
            Founded=(NeighborsA1 & NeighborsC1) - ( NeighborsB1)
            for i in Founded:
                key=tuple(set([A1,B1,C1,i]))
#                self.self_Graphlets["Graphlets5"][key]=[(A1,B1),(B1,C1),(C1,i),(i,A1)]
#                self.add_nodes_to_orbits(8,[A1,B1,C1,i])

                
            Founded=(NeighborsA1 & NeighborsC1) - (NeighborsB1)
            for i in Founded:
                key=tuple(set([A1,B1,C1,i]))
                self.self_Graphlets["Graphlets6"][key]=[(A1,C1),(C1,i),(i,A1),(C1,B1)]
#                self.Graphlets_Scores["Graphlets6"][key]=self.get_graphlet_score([(A1,C1),(C1,i),(i,A1),(C1,B1)])
#                self.add_nodes_to_orbits(9,[B1])
#                self.add_nodes_to_orbits(10,[A1,i])
#                self.add_nodes_to_orbits(11,[C1])
                
            Founded=(NeighborsB1 & NeighborsC1) - (NeighborsA1)
            for i in Founded:
                key=tuple(set([A1,B1,C1,i]))
#                self.self_Graphlets["Graphlets6"][key]=[(B1,C1),(C1,i),(i,C1),(C1,A1)]
#                self.add_nodes_to_orbits(9,[A1])
#                self.add_nodes_to_orbits(10,[B1,i])
#                self.add_nodes_to_orbits(11,[C1])

                
            Founded=NeighborsA1 & NeighborsC1 & NeighborsB1
            for i in Founded:
                key=tuple(set([A1,B1,C1,i]))
                self.self_Graphlets["Graphlets7"][key]=[(A1,B1),(B1,C1),(A1,i),(B1,i),(C1,i)]
#                self.Graphlets_Scores["Graphlets7"][key]=self.get_graphlet_score([(A1,B1),(B1,C1),(A1,i),(B1,i),(C1,i)])               
#                self.add_nodes_to_orbits(12,[A1,C1])
#                self.add_nodes_to_orbits(13,[B1,i])
                
        return (self.Graphlets, self.Orbits)             
    
                
    def find_graphlet_from_G2(self):  
        for graphlet in list(self.self_Graphlets["Graphlets2"].values()):
            A,B,C=graphlet[0][0],graphlet[0][1],graphlet[1][1]
            NeighborsA=self.neighbors[A] 
            NeighborsB=self.neighbors[B]
            NeighborsC=self.neighbors[C]
            
            results= (NeighborsA & NeighborsB & NeighborsC) - self.nodes
            for each in results:
                key=tuple(set([A,B,C,each]))
                self.Graphlets["Graphlets8"][key]=[(A,B),(B,C),(C,A),(A,each),(B,each),(C,each)]
                self.Graphlets_Scores["Graphlets8"][key]=self.get_graphlet_score([(A,B),(B,C),(C,A),(A,each),(B,each),(C,each)])               
#                self.Orbits["Orbits14"].update([A,B,C,each])
                self.add_nodes_to_orbits(14,[A,B,C,each])
                
        #self Graphlets        
        for graphlet in list(self.self_Graphlets["Graphlets2"].values()):
            A,B,C=graphlet[0][0],graphlet[0][1],graphlet[1][1]
            nodes={A,B,C}
            
            NeighborsA=(self.neighbors[A] & self.nodes) - nodes
            NeighborsB=(self.neighbors[B] & self.nodes) - nodes
            NeighborsC=(self.neighbors[C] & self.nodes) - nodes
            
            Founded= (NeighborsA & NeighborsB & NeighborsC)          
            for i in Founded:
                key=tuple(set([A,B,C,i]))
                self.self_Graphlets["Graphlets8"][key]=[(A,B),(B,C),(C,A),(A,i),(B,i),(C,i)]
#                self.add_nodes_to_orbits(14,[A,B,C,i])
            
            Founded= NeighborsA - (NeighborsB | NeighborsC)          
            for i in Founded:
                key=tuple(set([A,B,C,i]))
                self.self_Graphlets["Graphlets6"][key]=[(B,A),(A,C),(C,B),(A,i)]
#                self.add_nodes_to_orbits(9,[i])
#                self.add_nodes_to_orbits(10,[B,C])
#                self.add_nodes_to_orbits(11,[A])
            
            Founded= NeighborsB - (NeighborsA | NeighborsC)          
            for i in Founded:
                key=tuple(set([A,B,C,i]))
                self.self_Graphlets["Graphlets6"][key]=[(A,B),(B,C),(C,A),(B,i)]
#                self.add_nodes_to_orbits(9,[i])
#                self.add_nodes_to_orbits(10,[A,C])
#                self.add_nodes_to_orbits(11,[B])
            
            Founded= NeighborsC - (NeighborsB | NeighborsA)          
            for i in Founded:
                key=tuple(set([A,B,C,i]))
                self.self_Graphlets["Graphlets6"][key]=[(B,C),(C,A),(A,B),(C,i)]
#                self.add_nodes_to_orbits(9,[i])
#                self.add_nodes_to_orbits(10,[B,A])
#                self.add_nodes_to_orbits(11,[C])
            
            Founded= NeighborsA & NeighborsB - NeighborsC         
            for i in Founded:
                key=tuple(set([A,B,C,i]))
                self.self_Graphlets["Graphlets7"][key]=[(A,B),(B,C),(C,A),(A,i),(i,B)]
#                self.add_nodes_to_orbits(12,[A,B])
#                self.add_nodes_to_orbits(11,[C,i])
            
            Founded= NeighborsA & NeighborsC - NeighborsB         
            for i in Founded:
                key=tuple(set([A,B,C,i]))
                self.self_Graphlets["Graphlets7"][key]=[(A,C),(C,B),(B,A),(A,i),(i,C)]
#                self.add_nodes_to_orbits(12,[A,C])
#                self.add_nodes_to_orbits(11,[B,i])
            
            Founded= NeighborsC & NeighborsB - NeighborsA         
            for i in Founded:
                key=tuple(set([A,B,C,i]))
                self.self_Graphlets["Graphlets7"][key]=[(C,B),(B,A),(A,C),(C,i),(i,B)]
#                self.add_nodes_to_orbits(12,[C,B])
#                self.add_nodes_to_orbits(1                
               
##        print("the number of G8 :", len( self.Graphlets["Graphlets8"]))
        return (self.Graphlets, self.Orbits)        
 

    def find_key_graphlets (self):#,NodeA,NodeB): 
        Nodes=list(self.nodes)
        for i in range(len(self.nodes)):
            NodeA=Nodes[i]
            Neighbors_A=self.neighbors[NodeA]
            for k in range(i+1,len(self.nodes)):
                NodeB=Nodes[k]
                Neighbors_B=self.neighbors[NodeB]
        
                Propagated_nodes=(Neighbors_A & Neighbors_B) - self.nodes
                
                if NodeA in Neighbors_B:
                    key=tuple(set([NodeA,NodeB]))
                    self.Graphlets["Graphlets0"][key]=[(NodeA,NodeB)]
                    self.add_nodes_to_orbits(0,[NodeA,NodeB])
                    
                    self.Graphlets_Scores["Graphlets0"][key]=self.get_graphlet_score([(NodeA,NodeB)])
                        
#                        
#                    self.Orbits["Orbits0"].update(key)
                    for node in Propagated_nodes:
                        key=tuple(set([NodeA,node,NodeB]))
                        self.Graphlets["Graphlets2"][key]=[(NodeA,node),(node,NodeB),(NodeB,NodeA)]
                        self.Graphlets_Scores["Graphlets2"][key]=self.get_graphlet_score([(NodeA,node),(node,NodeB),(NodeB,NodeA)])
                        self.add_nodes_to_orbits(3,[NodeA,NodeB,node])
#                        self.Orbits["Orbits3"].update([NodeA,NodeB,node])
                    
                        
                        
                else:
                    for node in Propagated_nodes:
                        key=tuple(set([NodeA,node,NodeB]))
                        self.Graphlets["Graphlets1"][key]=[(NodeA,node),(node,NodeB)]
                        self.Graphlets_Scores["Graphlets1"][key]=self.get_graphlet_score([(NodeA,node),(node,NodeB)])
                        self.add_nodes_to_orbits(1,[NodeA,NodeB])
                        self.add_nodes_to_orbits(2,[node])                         
#                        self.Orbits["Orbits1"].update([NodeA,NodeB])
#                        self.Orbits["Orbits2"].add(node)
        
        return (self.Graphlets, self.Orbits)
         


        
if __name__=="__main__":
    
    
    
# =============================================================================
#     import argparse
#     
#     Epi=""
#     
#     description= """asdsds """#
# 
#     argument = argparse.ArgumentParser(description=description,
#                                        epilog=Epi)
#     
#     argument.add_argument("-r", "--reference_interactome_file", required=False)
#     
#     argument.add_argument("-i", "--input_node_file", required=False)
#     
#     argument.add_argument("-oi", "--output_interactome_file", required=False)
#     
#     
#     argument.add_argument("-og", "--output_graphlets_pickle", required=False)
#     
#     argument.add_argument("-oo", "--output_orbits_pickle", required=False)
#     
#     argument.add_argument("-d", "--derivate", action='store_true',required=False)
#     
#     
# 
# 
#     
#     args = vars(argument.parse_args())
#     
# 
#     if args["reference_interactome_file"] is None or args["input_node_file"] is None or args["output_interactome_file"] is None:
#         raise NameError("""
#               obligatory files
#               -r -----> reference_interactome_file  
#               -i -----> input_node_file
#              -oi -----> output_interactome_file
#               
#                optional graphlets
#               -d -----> find out derivated graphlets
#                          "G4 (Graphlet4),G5,G6,G7,G8
#                           G10,G11,G12,G13,G14,G15,G16
#                           G17,G18,G19,G20,G21,G22,G23
#                           G24,G25,G26,G27,G28,G29"
#              
#              optional output files
#              -og -----> output_graphlets_pickle
#              -oo -----> output_orbits_pickle
#         
#         """)
#     
#     Ref_df=pd.read_csv(args["reference_interactome_file"],sep="\t")
#     print("reference_interactome\n",args["reference_interactome_file"])
#     Ref_nx=nx.from_pandas_edgelist(Ref_df,"Gene_1","Gene_2")
#     
#     grf=Graphlet(Ref_nx)
#     
#     
#     Target_df=pd.read_csv(args["input_node_file"],sep="\t")[["Gene_name"]]
#     Target=Target_df.Gene_name
#     if args["derivate"]:
#         print("derivated graphlets from ",args["input_node_file"])
#         Res=grf.find_graphlets(Target,derivated_graphlets=True)
#     else:
#         print("underivated graphlets from ",args["input_node_file"])
#         Res=grf.find_graphlets(Target,derivated_graphlets=False)
#     
#     grf.write_created_network(args["output_interactome_file"])
#     print("the related subgraph was saved ",args["output_interactome_file"])   
# 
# 
#     if  args["output_graphlets_pickle"] :
#         print("the graphlets were stored at ",args["output_graphlets_pickle"])
#         grf.save_pickle_graphlets(args["output_graphlets_pickle"])
#     
#     if  args["output_orbits_pickle"] :
#         print(" the orbits were stored at ", args["output_orbits_pickle"])
#         grf.save_pickle_graphlets(args["output_orbits_pickle"])             
#     print ("\n"*3)     
# =============================================================================


                
#    Hippie_DF=pd.read_csv("../Source/Interactomes/HIPPIE.tab",sep="\t")
#    Hippie_DF=pd.read_csv("../Source/Filtered_HIPPIE/Breast_related_HIPPIE.sif",sep="\t")
#    Hippie_nx=nx.from_pandas_edgelist(Hippie_DF,"Gene_1","Gene_2","Score")
    #print (Hippie_DF["Gene_2"])
#    import timeit
#    
#    start=timeit.default_timer()
#    grf=GraphletSearch(Hippie_nx)
#    Target_df=pd.read_csv("../Source/NetPath__/Sampling_0_5/TCR_05A/TCR_05A_var_0.nodes",sep="\t")[["name"]]
#    Target=list(Target_df.name)
#    print (Target_df)
    
#    targets=["MP2K1","CXA1","ACTA"]
#    Res=grf.find_graphlets(Target,derivated_graphlets=True) 
#    grf.write_created_network("trial_TCR_05_node_0")
#    grf.save_pickle_graphlets("trial_TCR_05_graphlet_0")
#    grf.save_pickle_orbits("trial_TCR_05_orbits")
#    stop=timeit.default_timer()
#    print (stop-start)
    

    edges=[[0,1],[1,2],[2,3],[1,3],[2,4],[4,1],[4,3]]
    network=nx.Graph(edges)
    grf=GraphletSearch(network)
    nx.draw(network,with_labels=True, node_color="darkorange")
    grf.find_graphlets([2,4,3],derivated_graphlets=True)
    
    print(grf.get_selected_edge_list())
#    
#    print (grf.find_graphlets([2,4,3],derivated_graphlets=True))
    
    
#    grf.write_created_network("trial_egfr_100node")
    
# =============================================================================
#     
# =====================================================================