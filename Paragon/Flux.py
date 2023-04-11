# -*- coding: utf-8 -*-
"""
Created on Tue Dec 28 23:48:05 2021

@author: mkaanarici
"""
import networkx as nx, sys, numpy as np, pickle, math

class PageRankFlux():
    def __init__(self,network,motif_network,edge_attribute=None):


        self.network=network
        nan_nodes = []
        for node in self.network.nodes():
            if isinstance(node, float):
                if math.isnan(node):
                    nan_nodes.append(node)            
        for node in nan_nodes:
            self.network.remove_node(node)
        

        if edge_attribute==None:
            nx.set_edge_attributes(self.network,1,"Score")
            self.edge_attribute="Score"
            
        elif edge_attribute:
            self.edge_attribute=edge_attribute

        self.motif_nx=motif_network
        nan_nodes = []
        for node in self.motif_nx.nodes():
            if isinstance(node, float):
                if math.isnan(node):
                    nan_nodes.append(node)            
        for node in nan_nodes:
            self.motif_nx.remove_node(node)
    
        
    def load_initial_nodes(self,initial_nodes,weight=False):
        if weight==False:
            init=set(initial_nodes).intersection(set(self.network.nodes))
            self.initial_nodes=dict(zip(init,[1 for _ in range(len(initial_nodes))]))
        
        elif weight:#.any():
            self.initial_nodes=dict()
            for ind in range (len(initial_nodes)):
                if initial_nodes[ind] in self.network.nodes:
                    self.initial_nodes[initial_nodes[ind]]=weight[ind]
            #init=set(initial_nodes).intersection(set(self.network.nodes))
            #self.initial_nodes=dict(zip(init,[1 for _ in range(len(initial_nodes))]))
        else:
            print("Something wrong between weights and nodes")
        
        
    def reconstruct_subnetwork(self,alpha=0.8, threshold=0.5,max_edge_count=1000,intermediate_only=True): 
#    def pagerank(self,alpha=0.8, max_edge_count=500):
        self.pagerank_score= nx.pagerank(self.network,alpha=alpha,personalization=self.initial_nodes,weight=self.edge_attribute) 
        self.result_nx=nx.Graph()
        tot_flux=0
        for u in self.pagerank_score:
            denom = self.network.degree(u)
            for v in self.network.neighbors(u): #  neighbors
                
                flux1 = self.pagerank_score[u]*self.network[u][v][self.edge_attribute]/denom
                flux2 = self.pagerank_score[v]*self.network[u][v][self.edge_attribute]/denom
                flux_Score=min([flux1,flux2])
                self.network[u][v]['flux']=flux_Score
                if flux_Score==0:
                    self.network[u][v]['neglog_flux']=sys.float_info.max ## maximum float value
                else:
                    self.network[u][v]['neglog_flux'] = - (np.log10(flux_Score)) 
                    tot_flux+=flux_Score
        
        Motif_tot_flux=0
        for u,v in self.motif_nx.edges:
            self.motif_nx[u][v]["flux"]=self.network[u][v]['flux']
            self.motif_nx[u][v]['neglog_flux']=self.network[u][v]['neglog_flux']
            Motif_tot_flux+=self.motif_nx[u][v]["flux"]
        Returned_edges=dict()
        flux_sum=0
        
        init=list(self.initial_nodes.keys())
        for u,v,d in sorted(self.motif_nx.edges(data=True), key=lambda t: t[2]['neglog_flux']):
            Returned_edges[(u,v)] = self.motif_nx[u][v]['neglog_flux']
            flux_sum+=self.motif_nx[u][v]['flux']
            if u in init:
                init.remove(u)
            if v in init:
                init.remove(v)
            if  flux_sum/Motif_tot_flux>threshold or len(Returned_edges.keys())>max_edge_count-1:# or len(init)==0 
            
#            if flux_sum/Motif_tot_flux > threshold:
                print('theshold of %f limits predictions to %d edges'  %(round(flux_sum/Motif_tot_flux,3),len(Returned_edges)))
                break
        self.result_nx.add_edges_from(Returned_edges.keys())
        
        
        if intermediate_only==True:
            propagated_nodes=set(self.result_nx) -  set (self.initial_nodes.keys())
            delete_nodes=[]
            for node in propagated_nodes:
                if self.result_nx.degree(node)==1:
                    delete_nodes.append(node) 
            self.result_nx.remove_nodes_from(delete_nodes)
        self.result_nx=nx.from_edgelist(self.result_nx.edges())    
        return self.result_nx     
    
    

    def save_created_network(self, f_name):
        with open (f'{f_name}.pickle',"wb") as f:
            pickle.dump(self.result_nx,f)
 
    def get_created_network(self):
        return self.result_nx

    def write_created_network(self, f_name):
        DF=nx.to_pandas_edgelist(self.result_nx)
        DF.to_csv(f'{f_name}.sif', sep="\t",index=False)
        return True
         

if __name__=="__main__":
    import pandas as pd
    tumor='GBM'
    Initial_df=pd.read_csv(f'../Source/TCGA/Input_from_2/{tumor}.freq',sep="\t")
    Initial_nodes=Initial_df.name
    
    
    HIPPIE_df=pd.read_csv("../Source/Interactomes/HIPPIE.tab", sep="\t")
    HIPPIE_nx=nx.from_pandas_edgelist(HIPPIE_df,"Gene_1","Gene_2","Score")
    
    GGN_df=pd.read_csv("GBM_extended_GGN.sif",sep="\t")
    GGN_nx=nx.from_pandas_edgelist(GGN_df)
    
    pgrf=PageRankFlux(HIPPIE_nx,GGN_nx,edge_attribute="Score")
    
    pgrf.load_initial_nodes(Initial_nodes)
    returned_nx=pgrf.pagerank( max_edge_count=1000,alpha=0.8,threshold=0.8)
           
        # -*- coding: utf-8 -*-
"""
Created on Tue Dec 28 23:48:05 2021

@author: mkaanarici
"""
import networkx as nx, sys, numpy as np, pickle, math

class PageRankFlux():
    def __init__(self,network,motif_network,edge_attribute=None):


        self.network=network
        nan_nodes = []
        for node in self.network.nodes():
            if isinstance(node, float):
                if math.isnan(node):
                    nan_nodes.append(node)            
        for node in nan_nodes:
            self.network.remove_node(node)
        

        if edge_attribute==None:
            nx.set_edge_attributes(self.network,1,"Score")
            self.edge_attribute="Score"
            
        elif edge_attribute:
            self.edge_attribute=edge_attribute

        self.motif_nx=motif_network
        nan_nodes = []
        for node in self.motif_nx.nodes():
            if isinstance(node, float):
                if math.isnan(node):
                    nan_nodes.append(node)            
        for node in nan_nodes:
            self.motif_nx.remove_node(node)
    
        
    def load_initial_nodes(self,initial_nodes,weight=False):
        if weight==False:
            init=set(initial_nodes).intersection(set(self.network.nodes))
            self.initial_nodes=dict(zip(init,[1 for _ in range(len(initial_nodes))]))
        
        elif weight:#.any():
            self.initial_nodes=dict()
            for ind in range (len(initial_nodes)):
                if initial_nodes[ind] in self.network.nodes:
                    self.initial_nodes[initial_nodes[ind]]=weight[ind]
            #init=set(initial_nodes).intersection(set(self.network.nodes))
            #self.initial_nodes=dict(zip(init,[1 for _ in range(len(initial_nodes))]))
        else:
            print("Something wrong between weights and nodes")
        
        
    def reconstruct_subnetwork(self,alpha=0.8, threshold=0.5,max_edge_count=1000,intermediate_only=True): 
#    def pagerank(self,alpha=0.8, max_edge_count=500):
        self.pagerank_score= nx.pagerank(self.network,alpha=alpha,personalization=self.initial_nodes,weight=self.edge_attribute) 
        self.result_nx=nx.Graph()
        tot_flux=0
        for u in self.pagerank_score:
            denom = self.network.degree(u)
            for v in self.network.neighbors(u): #  neighbors
                
                flux1 = self.pagerank_score[u]*self.network[u][v][self.edge_attribute]/denom
                flux2 = self.pagerank_score[v]*self.network[u][v][self.edge_attribute]/denom
                flux_Score=min([flux1,flux2])
                self.network[u][v]['flux']=flux_Score
                if flux_Score==0:
                    self.network[u][v]['neglog_flux']=sys.float_info.max ## maximum float value
                else:
                    self.network[u][v]['neglog_flux'] = - (np.log10(flux_Score)) 
                    tot_flux+=flux_Score
        
        Motif_tot_flux=0
        for u,v in self.motif_nx.edges:
            self.motif_nx[u][v]["flux"]=self.network[u][v]['flux']
            self.motif_nx[u][v]['neglog_flux']=self.network[u][v]['neglog_flux']
            Motif_tot_flux+=self.motif_nx[u][v]["flux"]
        Returned_edges=dict()
        flux_sum=0
        
        init=list(self.initial_nodes.keys())
        for u,v,d in sorted(self.motif_nx.edges(data=True), key=lambda t: t[2]['neglog_flux']):
            Returned_edges[(u,v)] = self.motif_nx[u][v]['neglog_flux']
            flux_sum+=self.motif_nx[u][v]['flux']
            if u in init:
                init.remove(u)
            if v in init:
                init.remove(v)
            if  flux_sum/Motif_tot_flux>threshold or len(Returned_edges.keys())>max_edge_count-1:# or len(init)==0 
            
#            if flux_sum/Motif_tot_flux > threshold:
                print('theshold of %f limits predictions to %d edges'  %(round(flux_sum/Motif_tot_flux,3),len(Returned_edges)))
                break
        self.result_nx.add_edges_from(Returned_edges.keys())
        
        
        if intermediate_only==True:
            propagated_nodes=set(self.result_nx) -  set (self.initial_nodes.keys())
            delete_nodes=[]
            for node in propagated_nodes:
                if self.result_nx.degree(node)==1:
                    delete_nodes.append(node) 
            self.result_nx.remove_nodes_from(delete_nodes)
        self.result_nx=nx.from_edgelist(self.result_nx.edges())    
        return self.result_nx     
    
    

    def save_created_network(self, f_name):
        with open (f'{f_name}.pickle',"wb") as f:
            pickle.dump(self.result_nx,f)
 
    def get_created_network(self):
        return self.result_nx

    def write_created_network(self, f_name):
        DF=nx.to_pandas_edgelist(self.result_nx)
        DF.to_csv(f'{f_name}.sif', sep="\t",index=False)
        return True
         

if __name__=="__main__":
    import pandas as pd
    tumor='GBM'
    Initial_df=pd.read_csv(f'../Source/TCGA/Input_from_2/{tumor}.freq',sep="\t")
    Initial_nodes=Initial_df.name
    
    
    HIPPIE_df=pd.read_csv("../Source/Interactomes/HIPPIE.tab", sep="\t")
    HIPPIE_nx=nx.from_pandas_edgelist(HIPPIE_df,"Gene_1","Gene_2","Score")
    
    GGN_df=pd.read_csv("GBM_extended_GGN.sif",sep="\t")
    GGN_nx=nx.from_pandas_edgelist(GGN_df)
    
    pgrf=PageRankFlux(HIPPIE_nx,GGN_nx,edge_attribute="Score")
    
    pgrf.load_initial_nodes(Initial_nodes)
    returned_nx=pgrf.pagerank( max_edge_count=1000,alpha=0.8,threshold=0.8)
           
        