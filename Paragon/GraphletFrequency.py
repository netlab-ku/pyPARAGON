# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 16:18:16 2021

@author: mkaanarici
"""

import networkx as nx
from Paragon.Graphlets import GraphletSearch
import numpy as np,scipy.stats
import pandas as pd,pickle, timeit,os


import multiprocessing as mp 

class GraphletSelect():
    def __init__(self, network, node_list):
        self.network=network
        self.node_list=node_list
        self.GraphletSearch=GraphletSearch(network)
        self.GraphletSearch.find_graphlets(self.node_list,derivated_graphlets=True) # graphlet
        
  
        
        
        self.Graphlets_Scores=self.GraphletSearch.get_graphlets_scores()
        
        
        self.Graphlets=self.GraphletSearch.get_graphlets()
        self.Orbits=self.GraphletSearch.get_orbits()
        self.frequency=self.find_Graphlet_Frequency(self.Graphlets)
        
        self.selected_Graphlets=dict()
        self.selected_Graphlets_Scores=dict()
        self.selected_Orbits=dict()
#        with open ("frequencies.pickle","rb") as f:
#            self.frequecies_pool=pickle.load(f)
        self.frequecies_pool=[self.frequency]        
        
        self.created_network=nx.Graph()
        
#        print("the number of G0 :",len(self.Graphlets["Graphlets0"].values()))
#        print("the number of G1 :",len(self.Graphlets["Graphlets1"].values()))
#        print("the number of G2 :",len(self.Graphlets["Graphlets2"].values()))
#        print("the number of G3 :", len( self.Graphlets["Graphlets3"]))
#        print("the number of G4 :", len( self.Graphlets["Graphlets4"]))
#        print("the number of G5 :", len( self.Graphlets["Graphlets5"]))
#        print("the number of G6 :", len( self.Graphlets["Graphlets6"]))
#        print("the number of G7 :", len( self.Graphlets["Graphlets7"])) 
#        print("the number of G8 :", len( self.Graphlets["Graphlets8"]))
        
    def get_created_network(self):
        return self.created_network
        
    def get_selected_edge_list(self):
        return list(self.created_network.edges)
    
    def get_get_graphlets_scores(self):
        return self.selected_Graphlets_Scores

    def write_created_network(self,file_name):
        DF=nx.to_pandas_edgelist(self.created_network)
        DF.to_csv(f'{file_name}.sif',sep="\t",index=False)
        return True


    def find_Graphlet_Frequency(self,Graphlets):
        frequency=dict()
#        print(Graphlets.keys())
        total_graphlets_count=0
        for i in range(9):
            total_graphlets_count+=len(Graphlets[f'Graphlets{i}'])
#        print("graphlets count:",total_graphlets_count)
        for i in range(9):
            if total_graphlets_count==0:
                frequency[f'Graphlets{i}']=0
            else:               
                frequency[f'Graphlets{i}']=len(Graphlets[f'Graphlets{i}'])/total_graphlets_count
#        print(frequency)
        return frequency
    
    def get_Graphlet_Frequency(self):
        return self.frequency
    
    
    def get_Z_score(self):
        self.z_score=dict()
        rand_freq=[]
        for i in range(1,len(self.frequecies_pool)):
            temp=[]
            for k in range(9):
                temp.append(self.frequecies_pool[i][f'Graphlets{k}'])
            rand_freq.append(temp)
        rand_freq=np.array(rand_freq)
        for k in range(9):
            own_frq=self.frequency[f'Graphlets{k}']
            sample=rand_freq[:,k]
            Mean=np.mean(sample)
            std=np.std(sample)
            z_scores=(own_frq-Mean)/std
            p_values = scipy.stats.norm.sf(abs(z_scores))
            self.z_score[f'Graphlets{k}']=z_scores,p_values
        return self.z_score

    
    
    def permutate_network(self):
        Connected_components=[self.network.subgraph(c).copy() for c in nx.connected_components(self.network)]
        Permutated_nx=nx.Graph()
        for subgraph in Connected_components:
            if subgraph.number_of_nodes()<4:
                Permutated_nx=nx.union(Permutated_nx,subgraph)
                continue
            nx.connected_double_edge_swap(subgraph,subgraph.number_of_edges())
            Permutated_nx=nx.union(Permutated_nx,subgraph)
        return Permutated_nx
    
    def save_permuted_network(self,f_name):
        Permutated_nx=self.permutate_network()
        Permutated_df=nx.to_pandas_edgelist(Permutated_nx)
        Permutated_df=Permutated_df.rename(columns={"source":"Gene_1","target":"Gene_2"})
        Permutated_df.to_csv(f'{f_name}.tab',index=False, sep="\t")
        return f_name
        
        
        

    def permuted_frequencies(self):
        permutated_nx=self.permutate_network(self.network)
        Per_Grph_Src=GraphletSearch(permutated_nx)
        Per_Grph_Src.find_graphlets(self.node_list,derivated_graphlets=True)
        return self.find_Graphlet_Frequency(Per_Grph_Src.get_graphlets())
    
    def save_permuted_networks(self, f_name,size):
        argums=[[f'{f_name}_{i+1}'] for i in range(size)]
#        strt=timeit.default_timer()
        if __name__=="__main__":
            with mp.Pool() as pool:
                results=pool.starmap(self.save_permuted_network, argums)
            for result in results:
                print("permuted network saved in\t", result)
#            print("permuted network takes mp", timeit.default_timer()-strt)
            return 1
        else:
            for each in argums:
                self.save_permuted_network(each[0])
                print("permuted network saved in\t", each[0])
#            print("permuted network takes", timeit.default_timer()-strt)
            return 1
                
                
        
        with open (f'{f_name}.pickle',"wb") as f:
            pickle.dump(self.frequecies_pool,f)
            
    def select_significant_Graphlets(self):
        orbits={0:[0],1:[1,2],2:[3],3:[4,5],4:[6,7],5:[8],6:[9,10,11],7:[12,13],8:[14]}
        for i in range(9):
            if self.z_score[f'Graphlets{i}'][0]>3.27:
                self.selected_Graphlets[f'Graphlets{i}']=self.Graphlets[f'Graphlets{i}']
                self.selected_Graphlets_Scores[f'Graphlets{i}']=self.Graphlets_Scores[f'Graphlets{i}']
                for graphlet in self.Graphlets[f'Graphlets{i}'].values():
                    self.created_network.add_edges_from(graphlet)
                    for o in orbits[i]:
                        self.Orbits[f'Orbits{o}']=self.Orbits[f'Orbits{o}']
                
#            else:
#                self.selected_Graphlets[f'Graphlets{i}']={}
        return self.selected_Graphlets
    
    def calculate_frequency(self, network_file,node_list):
        temp_df=pd.read_csv(f'{network_file}',sep="\t")
        temp_nx=nx.from_pandas_edgelist(temp_df,"Gene_1","Gene_2")
        temp_GraphletSearch=GraphletSearch(temp_nx)
        temp_GraphletSearch.find_graphlets(node_list,derivated_graphlets=True) 
        temp_Graphlets=temp_GraphletSearch.get_graphlets()
        return self.find_Graphlet_Frequency(temp_Graphlets)
        
    
    def get_frequencies_from_pool(self, permuted_files_path):
        
#        if __name__=="__main__":
#            items=[(f'{permuted_files_path}/{i}',self.node_list) for i in os.listdir(f'{permuted_files_path}')]
#            with mp.Pool() as pool:
#                results=pool.starmap(self.calculate_frequency, items)
#            for result in results:
#                self.frequecies_pool.append(result)
#            return len(self.frequecies_pool)
#            
#        else:
        permuted_tabs=[i for i in os.listdir(f'{permuted_files_path}')]
        for file in permuted_tabs:
            temp_df=pd.read_csv(f'{permuted_files_path}/{file}',sep="\t")
            temp_nx=nx.from_pandas_edgelist(temp_df,"Gene_1","Gene_2")
            temp_GraphletSearch=GraphletSearch(temp_nx)
            temp_GraphletSearch.find_graphlets(self.node_list,derivated_graphlets=True) 
            temp_Graphlets=temp_GraphletSearch.get_graphlets()
            temp_frequency=self.find_Graphlet_Frequency(temp_Graphlets)
            self.frequecies_pool.append(temp_frequency)
        return self.frequecies_pool
                
    def save_pickle_zscore(self,file_name):
        with open(f'{file_name}.pickle',"wb") as f:
            pickle.dump(self.z_score,f)
            
    def save_pickle_graphlets_scores(self,file_name):
        with open(f'{file_name}.pickle',"wb") as f:
            pickle.dump(self.selected_Graphlets_Scores,f)            
    
    def save_pickle_graphlets(self,file_name):
        with open(f'{file_name}.pickle',"wb") as f:
            pickle.dump(self.selected_Graphlets,f)
    
    def save_pickle_orbits(self,file_name):
        with open(f'{file_name}.pickle',"wb") as f:
            pickle.dump(self.selected_Orbits,f)   

    def save_pickle_frequencies(self,file_name):
        with open(f'{file_name}.pickle',"wb") as f:
            pickle.dump(self.frequecies_pool,f)  
            
    def get_selected_graphlets(self):
        return self.selected_Graphlets
    
    def get_frequencies_from_pickle(self, frequency_file):   
        with open(f'{frequency_file}',"rb") as f:
            self.frequecies_pool=pickle.load(f)   
    
    def get_selected_orbits(self):
        return self.selected_Orbits
            
    def save_frequencies_into_pickle(self, f_name):
        with open (f'{f_name}.pickle',"wb") as f:
            pickle.dump(self.frequecies_pool,f)
            
            
    
    




if __name__=="__main__":
    
    
#    
#    # =============================================================================
#    import argparse
#     
#    Epi=""
#     
#    description= """asdsds """#
# 
#    argument = argparse.ArgumentParser(description=description,
#                                        epilog=Epi)
#     
#    argument.add_argument("-r", "--reference_interactome_file", required=False)
#     
#    argument.add_argument("-i", "--input_node_file", required=False)
#     
#    argument.add_argument("-oi", "--output_interactome_file", required=False)
#     
#    argument.add_argument("-of", "--output_frequency_pickle", required=False)
#     
#     
#    argument.add_argument("-og", "--output_graphlets_pickle", required=False)
#     
#    argument.add_argument("-oo", "--output_orbits_pickle", required=False)
#     
#    argument.add_argument("-oz", "--output_zscores", required=False)     
#     
#    argument.add_argument("-pn", "--permuted_network", required=False)
#     
#     
# 
# 
#     
#    args = vars(argument.parse_args())
#     
#
#
#    
#    start=timeit.default_timer()
#    
#    Hippie_DF=pd.read_csv(args["reference_interactome_file"],sep="\t")
#    Hippie_nx=nx.from_pandas_edgelist(Hippie_DF,"Gene_1","Gene_2")
#    
#    
#    Target_df=pd.read_csv(args["input_node_file"],sep="\t")[["name"]]
#    Target=list(Target_df.name)
#    
#    Select=GraphletSelect(Hippie_nx,Target)
#    Select.get_frequencies_from_pool(args["permuted_network"])
#    
#    print(Select.get_Z_score())
#    Select.save_pickle_zscore(args["output_zscores"])
#    print(Select.select_significant_Graphlets().keys())
#    
#    print(Select.get_selected_graphlets().keys())
#    print(Select.get_selected_orbits().keys())
#    
#    Select.save_pickle_graphlets(args["output_graphlets_pickle"])
#    Select.save_pickle_orbits(args["output_orbits_pickle"])
#    Select.save_pickle_frequencies(args["output_frequency_pickle"])
#    Select.write_created_network(args["output_interactome_file"])
#    stop=timeit.default_timer()
#    print (stop-start)

        
    
    

 

# =============================================================================
# =============================================================================
#      reselection of Graphlets (linux generated different graphlets )
# =============================================================================

     
     Hippie_DF=pd.read_csv("../Source/Interactomes/HIPPIE.tab",sep="\t")
     Hippie_nx=nx.from_pandas_edgelist(Hippie_DF,"Gene_1","Gene_2",  edge_attr='Score')
     
     
     pathway="TCR"
     i=0
     
     Pathways=["TCR","TGFbetaReceptor","TNFalpha","Wnt"]
 
     for pathway in Pathways:
         for i in range(5):
             print(pathway, i)
             
             os.makedirs(f'../Results/HIPPIE/NetPath__/Sampling_0_8/Selected/{pathway}_08A/Graphlet_Scores/',exist_ok=True)
        
             
             start=timeit.default_timer()
             Target_df=pd.read_csv(f'../Source/NetPath__/Sampling_0_8/{pathway}_08A/{pathway}_08A_var_{i}.nodes',sep="\t")[["name"]]
             Target=list(Target_df.name)
             
             Select=GraphletSelect(Hippie_nx,Target)     
             output=f'{pathway}_08A_var_{i}'
             
             Select.get_frequencies_from_pickle(f'../Results/HIPPIE/NetPath__/Sampling_0_8/Selected/{pathway}_08A/Frequencies/{output}.pickle')
        
        
             Select.get_Z_score()
             Select.save_pickle_zscore(f'../Results/HIPPIE/NetPath__/Sampling_0_8/Selected/{pathway}_08A/Zscores/{output}')
             Select.select_significant_Graphlets().keys()
             Select.save_pickle_graphlets_scores(f'../Results/HIPPIE/NetPath__/Sampling_0_8/Selected/{pathway}_08A/Graphlet_Scores/{output}')
             print(Select.get_selected_graphlets().keys())
#             print(Select.get_selected_orbits().keys())
             
             Select.save_pickle_graphlets(f'../Results/HIPPIE/NetPath__/Sampling_0_8/Selected/{pathway}_08A/Graphlets/{output}')
             Select.save_pickle_orbits(f'../Results/HIPPIE/NetPath__/Sampling_0_8/Selected/{pathway}_08A/Orbits/{output}')
             Select.write_created_network(f'../Results/HIPPIE/NetPath__/Sampling_0_8/Selected/{pathway}_08A/SIFs/{output}')
             stop=timeit.default_timer()
             print (stop-start)
 
         
     
     
 
 
# =============================================================================
