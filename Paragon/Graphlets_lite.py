# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 10:08:08 2022

@author: mkaanarici
"""
import networkx as nx

class GraphletSearch():
    def __init__(self,network,node_list=None):
        self.network=network
        self.neighbors=dict()
        for node in self.network.nodes:
            self.neighbors[node]=set([i for i in self.network.neighbors(node) if str(i)!="nan"])
        
        if node_list!=None:
            self.nodes=set([i for i in node_list if i in self.network.nodes])
            print (f'{len(self.nodes)} of {len(node_list)} input nodes have been found in the given network\n' )
            if len(set(node_list)-self.nodes)!=0:
                print(f'{set(node_list)-self.nodes} nonexisting nodes in the reference network\n' )
      
        else:
            self.nodes=None
        
        self.intermediate_nodes=set()#self.propagated_node_list=set()
        self.GGN=nx.Graph()
        self.nodes=set()
        
        self.targetted_graphlets=[f'Graphlets{g}' for g in range(9)]
        
        self.Key_Graphlets=dict(zip([f'Graphlets{g}' for g in range(3)],[{} for i in range(3)]))
        self.self_Graphlets=dict(zip([f'Graphlets{g}' for g in range(3)],[{} for i in range(3)]))
        
        
    def construct_GGN(self, node_list=None, Graphlets=None, extention=False):
        
        if Graphlets:
            self.targetted_graphlets=Graphlets
        
        
        #self.nodes=set([i for i in node_list if i in self.network.nodes])
        
        if node_list==None:
            print(' define a initial node set')
            return False
        
        else:
            self.nodes=set([i for i in node_list if i in self.network.nodes])
                        
        if len(self.nodes)>0:            
            self.nodes=set([i for i in node_list if i in self.network.nodes])
            print (f'{len(self.nodes)} of {len(node_list)} input nodes have been found in the given network\n' )
            if len(set(node_list)-self.nodes)!=0:
                print(f'{set(node_list)-self.nodes} nonexisting nodes in the reference network\n' )
        
        self.find_key_graphlets()
        self.find_self_graphlets()
        self.find_four_nodes_graphlets()
        
        
        
        self.intermediate_nodes= set(self.GGN.nodes) - self.nodes
        if extention:
            self.__extend_GGN()
        
        return self.GGN
   
    def write_guided_graphlet_network(self,file_name):
        DF=nx.to_pandas_edgelist(self.GGN)
        DF.to_csv(f'{file_name}.sif',sep="\t",index=False)
        return True
    
    
    def get_GGN_edges(self):
        return self.GGN.edges

    def __extend_GGN(self):
        intermediate_nodes=list(self.intermediate_nodes)
        print("run, intermediates",len(intermediate_nodes))
        for i in range(len(intermediate_nodes)):
            for k in range(i+1,len(intermediate_nodes)):
                if(intermediate_nodes[i],intermediate_nodes[k]) in self.network.edges:
                    self.GGN.add_edge(intermediate_nodes[i],intermediate_nodes[k])


    def find_self_graphlets(self):    
        Nodes=list(self.nodes)
        for i in range(len(self.nodes)):
            node1=Nodes[i]
            
            for k in range(i+1,len(self.nodes)):
                node2=Nodes[k] 
                if node1==node2: 
                    continue 
                
                for l in range(k+1,len(self.nodes)):
                    node3=Nodes[l]
                    if node3==node2 or node3==node1:
                        continue
 
                    if ((node1,node2) in self.network.edges) and ((node1,node3) in self.network.edges) and ((node2,node3) in self.network.edges):
                        key=tuple(set([node1,node2,node3]))
                        self.self_Graphlets["Graphlets2"][key]=[(node1,node2),(node2,node3),(node3,node1)]
                        
                    elif ((node1,node2) in self.network.edges) and ((node1,node3) in self.network.edges) and ((node2,node3) not in self.network.edges):
                        key=tuple(set([node1,node2,node3]))
                        self.self_Graphlets["Graphlets1"][key]=[(node2,node1),(node1,node3)]

                        
                    elif ((node1,node2) in self.network.edges) and ((node1,node3) not in self.network.edges) and ((node2,node3) in self.network.edges):
                        key=tuple(set([node1,node2,node3]))
                        self.self_Graphlets["Graphlets1"][key]=[(node1,node2),(node2,node3)]
                      
                        
                    elif ((node1,node2) not in self.network.edges) and ((node1,node3) in self.network.edges) and ((node2,node3) in self.network.edges):
                        key=tuple(set([node1,node2,node3]))
                        self.self_Graphlets["Graphlets1"][key]=[(node1,node3),(node3,node2)]
 
    
        
    def find_key_graphlets (self):#,NodeA,NodeB): 
        Nodes=list(self.nodes)
        for i in range(len(self.nodes)):
            NodeA=Nodes[i]
            Neighbors_A=self.neighbors[NodeA]
            for k in range(i+1,len(self.nodes)):
                NodeB=Nodes[k]
                Neighbors_B=self.neighbors[NodeB]
        
                intermediate_nodes=(Neighbors_A & Neighbors_B) - self.nodes
                
                if NodeA in Neighbors_B:
                    key=tuple(set([NodeA,NodeB]))
                    self.Key_Graphlets["Graphlets0"][key]=[(NodeA,NodeB)]
                    if "Graphlets0" in self.targetted_graphlets:
                        self.GGN.add_edge(NodeA,NodeB)
                    
                        
                    for node in intermediate_nodes :
                        key=tuple(set([NodeA,node,NodeB]))
                        self.Key_Graphlets["Graphlets2"][key]=[(NodeA,node),(node,NodeB),(NodeB,NodeA)]
                        if "Graphlets2" in self.targetted_graphlets:
                            self.GGN.add_edges_from([(NodeA,node),(node,NodeB),(NodeB,NodeA)])
                        
                else:
                    for node in intermediate_nodes:
                        key=tuple(set([NodeA,node,NodeB]))
                        self.Key_Graphlets["Graphlets1"][key]=[(NodeA,node),(node,NodeB)]
                        if "Graphlets1" in self.targetted_graphlets:
                            self.GGN.add_edges_from([(NodeA,node),(node,NodeB)])
        
    def find_four_nodes_graphlets(self):
        G1_derivated=list(self.Key_Graphlets["Graphlets1"].values())
        for i in range(len(G1_derivated)):
            g1=G1_derivated[i]
            A1,B1,i1=(g1[0][0],g1[1][1],g1[0][1]) # A1->i1->B1 i1 : intermediate
            nodes=set([A1,B1,i1])
            
            NeighborsA1=self.neighbors[A1] - nodes
            NeighborsB1=self.neighbors[B1] - nodes
            Neighborsi1=self.neighbors[i1] - nodes
            
            
            ## Graphlet3 
            Founded=(NeighborsA1 & self.nodes) - (NeighborsB1 | Neighborsi1)
            for C1 in Founded:
                if "Graphlets3" in self.targetted_graphlets:
                    self.GGN.add_edges_from([(C1,A1),(A1,i1),(i1,B1)]) # C1->A1->i1->B1 , i1 intermediate

            ### Graphlet3                 
            Founded=(NeighborsB1 & self.nodes) - (NeighborsA1 | Neighborsi1)
            for C1 in Founded:
                if "Graphlets3" in self.targetted_graphlets:
                    self.GGN.add_edges_from([(C1,B1),(B1,i1),(i1,A1)])  ## C1->B1->i1->A1 , i1 intermediate

            ### Graphlet4                          
            Founded=(Neighborsi1 & self.nodes) - (NeighborsA1 | NeighborsB1)
            for C1 in Founded:
                if "Graphlets4" in self.targetted_graphlets:
                    self.GGN.add_edges_from([(A1,i1),(B1,i1),(C1,i1)]) # A1->i1, B1->i1, C1->i1, i1 intermediate                
   

            ### Graphlet6                   
            Founded=(NeighborsA1 & Neighborsi1 & self.nodes) - (NeighborsB1)
            for C1 in Founded:
                if "Graphlets6" in self.targetted_graphlets:
                    self.GGN.add_edges_from([(B1,i1),(i1,A1),(A1,C1),(C1,i1)]) # B1->i1, i1->A1->C1->i1, i1 intermediate   

            ### Graphlet6                   
            Founded=(NeighborsB1 & Neighborsi1 & self.nodes) - (NeighborsA1)
            for C1 in Founded:
                if "Graphlets6" in self.targetted_graphlets:
                    self.GGN.add_edges_from([(A1,i1),(i1,B1),(B1,C1),(C1,i1)]) # A1->i1, i1->B1->C1->i1, i1 intermediate  
                
                #####################
                
                #### self.Graphlets
                
                ####################

                
        for graphlet in list(self.self_Graphlets["Graphlets1"].values()):
            A,B,C=graphlet[0][0],graphlet[0][1],graphlet[1][1]
            NeighborsA=self.neighbors[A] 
            NeighborsB=self.neighbors[B]
            NeighborsC=self.neighbors[C]
            
            
            ### Graphlet5
            Founded=(NeighborsA & NeighborsC) - (self.nodes |  NeighborsB)
            for i1 in Founded:
                if "Graphlets5" in self.targetted_graphlets:
                    self.GGN.add_edges_from([(A,B),(B,C),(C,i1),(i1,A)]) # A->B->Ci->1->A 

            ### Graphlet7                
            Founded=(NeighborsA & NeighborsC & NeighborsB) - self.nodes
            for i1 in Founded:
                if "Graphlets7" in self.targetted_graphlets:
                    self.GGN.add_edges_from([(A,B),(B,C),(A,i1),(B,i1),(C,i1)])  # A->B->C, A-i1, B-i1, C-i1, i1 intermediate
              
            
        for graphlet in list(self.self_Graphlets["Graphlets2"].values()):
            A,B,C=graphlet[0][0],graphlet[0][1],graphlet[1][1]
            NeighborsA=self.neighbors[A] 
            NeighborsB=self.neighbors[B]
            NeighborsC=self.neighbors[C]
            
            
            ### Graphlets8
            results= (NeighborsA & NeighborsB & NeighborsC) - self.nodes
            for i1 in results:
                if "Graphlets8" in self.targetted_graphlets:
                    self.GGN.add_edges_from([(A,B),(B,C),(C,A),(A,i1),(B,i1),(C,i1)]) ## A->B->C->A, A->i1,B->i1,C->i1, i1 intermediate Node 
          



if __name__=="__main__":
#    edges=[[0,1],[1,2],[2,3],[1,3],[2,4],[4,1],[4,3],[0,4]]
#    network=nx.Graph(edges)
#    grf=GraphletSearch(network)
#    nx.draw(network,with_labels=True, node_color="darkorange")
#    grf.construct_GGN([2,0,1])
    
    import pandas as pd
    HIPPIE_df=pd.read_csv("../Source/Interactomes/HIPPIE.tab", sep="\t")
    HIPPIE_nx=nx.from_pandas_edgelist(HIPPIE_df,"Gene_1","Gene_2","Score")
    
    tumor='GBM'
    Initial_df=pd.read_csv(f'../Source/TCGA/Input_from_2/{tumor}.freq',sep="\t")
    Initial_nodes=Initial_df.name
    len(Initial_nodes)
    
    output_name=f'{tumor}'
    GRF=GraphletSearch(HIPPIE_nx)
    returned_nx=GRF.construct_GGN(Initial_nodes,extention=True)
    
             