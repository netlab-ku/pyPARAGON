import pandas as pd, numpy as np

from community import community_louvain

from scipy.stats import hypergeom

class OverRepresentationAnalysis():
    def __init__(self):
        self.reference_network=None
        
    def define_reference_network(self,reference_network):
        self.reference_network=reference_network
        
        
        
    def define_prior_knowledge(self, prior_knowledge_df, prior_knowledge_on, name_on): 
        self.prior_knowledge_on=prior_knowledge_on
        self.name_on=name_on       
        self.prior_knowledge=prior_knowledge_df[prior_knowledge_df[f'{self.name_on}'].isin(self.reference_network.nodes)]
        self.prior_knowledge=prior_knowledge_df.drop_duplicates().reset_index(drop=True)
        
        prior_knowledge_in_group_A=self.prior_knowledge.groupby(f'{self.prior_knowledge_on}')[f'{self.name_on}'].apply(list).reset_index()
        prior_knowledge_in_group_A=prior_knowledge_in_group_A.rename(columns={f'{self.name_on}':"components of prior_knowledge"})
        prior_knowledge_in_group_B=self.prior_knowledge.groupby(f'{self.prior_knowledge_on}').count().reset_index()
        
        prior_knowledge_in_group_B=prior_knowledge_in_group_B.rename(columns={f'{self.name_on}':"The number of components of prior_knowledge"})
        self.prior_knowledge_in_group=prior_knowledge_in_group_A.merge(prior_knowledge_in_group_B, on=f'{self.prior_knowledge_on}',how="inner")
        #print (self.prior_knowledge_in_group)
        
        
        return self.prior_knowledge
    
    def __hypergeometric_test_formula (self, X,M,n,N):
        """ x --> Successfully identified target genes in target process
            M --> Population size, the number of genes in Reference Network
            n --> the number of genes in module
            N --> the number of genes in target parocess """

        p_value=hypergeom.sf(X-1, M, n, N)
        enrichment_score=np.log2(hypergeom.pmf(X-1, M, n, N) / hypergeom.sf(X-1, M, n, N))
        return [p_value,enrichment_score]
    
    def hypergeometric_test(self,community):
        Temp_GOA_bp_df=self.prior_knowledge[self.prior_knowledge[f'{self.name_on}'].isin(community)].reset_index(drop=True)
        
        
        Temp_GOA_bp_df_set=Temp_GOA_bp_df.groupby(f'{self.prior_knowledge_on}')[f'{self.name_on}'].apply(list).reset_index()
        Temp_GOA_bp_df_set=Temp_GOA_bp_df_set.rename(columns={f'{self.name_on}':"Intersecting Genes"})
        #print(Temp_GOA_bp_df_set)
    
    
        Temp_GOA_bp_df_num=Temp_GOA_bp_df.groupby(f'{self.prior_knowledge_on}').count().reset_index()
        Temp_GOA_bp_df_num=Temp_GOA_bp_df_num.rename(columns={f'{self.name_on}':"The number of intersecting genes"}).reset_index(drop=True)
        Temp_GOA_bp_df_num=Temp_GOA_bp_df_num[Temp_GOA_bp_df_num["The number of intersecting genes"]>2]
        #print (Temp_GOA_bp_df_num)
        
        Temp_GOA_bp_df=Temp_GOA_bp_df_set.merge(Temp_GOA_bp_df_num, how="inner")[[f'{self.prior_knowledge_on}',"Intersecting Genes","The number of intersecting genes"]]
        #print(Temp_GOA_bp_df)
        
        Temp_GOA_bp_df=Temp_GOA_bp_df.merge(self.prior_knowledge_in_group, on=f'{self.prior_knowledge_on}',how="inner")
        Temp_list=Temp_GOA_bp_df.T.to_dict()
        
        Result_temp=pd.DataFrame(columns=[f'{self.prior_knowledge_on}',"p-value","Erichment_Score","Genes in Module","Intersecting Genes","The number of intersecting genes","Process_Gene","The number of components of prior_knowledge"])

        
        for index in range(len(Temp_list)):
            go_id=Temp_list[index][f'{self.prior_knowledge_on}']
            #print (go_id)
            x=Temp_list[index]["The number of intersecting genes"]

            #M --> Length_ref_nx,    Population size, the number of genes in Reference Network
            n=len(community)#The_num_of_genes_in_community
            N=Temp_list[index]["The number of components of prior_knowledge"]

            #print (x,Length_ref_nx,n,N)
            Length_ref_nx=self.reference_network.number_of_nodes()
            P_value,Enrichment_Score=self.__hypergeometric_test_formula (X=x,M=Length_ref_nx,n=n,N=N)
            dictionary={f'{self.prior_knowledge_on}':go_id,"p-value":P_value,"Erichment_Score":Enrichment_Score,
                        "Genes in Module":community,"Intersecting Genes":Temp_list[index]["Intersecting Genes"],
                        "The number of intersecting genes":x,"Process_Gene":Temp_list[index]["components of prior_knowledge"],
                        "The number of components of prior_knowledge":Temp_list[index]["The number of components of prior_knowledge"]}


            Result_temp=Result_temp.append(dictionary, ignore_index=True)
           
        return Result_temp



class CommunityAnalysis(OverRepresentationAnalysis):
    def __init__(self, subnetwork):
        self.network=subnetwork
        self.find_louvain_communities(self.network)
        
    def find_louvain_communities(self,graph):
        output=[]
        comms = community_louvain.best_partition(graph,resolution=1,random_state=1)

        Commns_num=list(range(max(comms.values())+1))
        Commns_dict=dict(zip(Commns_num,[[]for i in Commns_num]))

        if max(comms.values())==0:
            if len(Commns_dict[0])>0:
                output.append(Commns_dict[0])
            return output


        for node,cluster in comms.items():
            #print(node,cluster)
            Commns_dict[cluster].append(node)

        for k in Commns_num:
            if len(Commns_dict[k])<5:
                #print (f' module, {Commns_dict[k]} was deleted due to the inadequate number of nodes')
                del Commns_dict[k]
                continue

            elif 20 >= len(Commns_dict[k]) >= 5:
                output.append(Commns_dict[k])

            else:
                #print (len(Commns_dict[k]))
                subgraph=graph.subgraph(Commns_dict[k]).copy()
                sub_ouput=self.find_louvain_communities(subgraph)
                output.extend(sub_ouput)
        self.returned_modules=output
        return self.returned_modules
    
    def get_communities_in_list(self):   
        return self.returned_modules
    
    def get_communities_in_DataFrames(self):
        returned_modules_list={"Community_name":[],"Community":[]}

        result_node_list={"Genes":[],"Community":[]}
        i=0

        for module in self.returned_modules:

            i+=1
            gene_list_str=(";").join(module)
            Community_name=f'Module_{i}'

            returned_modules_list["Community_name"].append(Community_name)
            returned_modules_list["Community"].append(gene_list_str)

            for gene in module:
                result_node_list["Genes"].append(gene)
                result_node_list["Community"].append(Community_name)


        module_name=[f'Module_{i}' for i in range(1,len(self.returned_modules)+1)]
        self.module_df=pd.DataFrame({"Community_name":returned_modules_list["Community_name"],
                                "Community":returned_modules_list["Community"]})


        self.Node_df=pd.DataFrame({"Genes":result_node_list["Genes"],
                                "Community":result_node_list["Community"]})
        
        return self.module_df,self.Node_df
    
    def hypergeometric_test_for_all_communities(self,reference_network=None,prior_knowledge_df=None, prior_knowledge_on=None, name_on=None):
        
        if reference_network!=None:
            self.define_reference_network(reference_network)
        if isinstance(prior_knowledge_df,pd.DataFrame) and prior_knowledge_on!=None  and name_on!=None:
            self.define_prior_knowledge(prior_knowledge_df, prior_knowledge_on, name_on)
        
        
        self.all_hypergeometric_results=pd.DataFrame(columns=["Community_name"])
        
        for index in self.module_df.index:
            Community=self.module_df.iloc[index]["Community"].split(";")
            temp_test_result=self.hypergeometric_test(Community)
            temp_test_result["Community_name"]=self.module_df.iloc[index]["Community_name"]
            self.all_hypergeometric_results=pd.concat([self.all_hypergeometric_results,temp_test_result],ignore_index=True)
                    
        return self.all_hypergeometric_results

    def hypergeometric_test_for_community(self,community_name,reference_network=None,prior_knowledge_df=None, prior_knowledge_on=None, name_on=None):
        
        if reference_network!=None:
            self.define_reference_network(reference_network)
        if isinstance(prior_knowledge_df,pd.DataFrame) and prior_knowledge_on!=None  and name_on!=None:
            self.define_prior_knowledge(prior_knowledge_df, prior_knowledge_on, name_on)
        
        Community=self.module_df[self.module_df["Community_name"]==community_name]["Community"].str.split(";").to_list()[0]
        temp_test_result=self.hypergeometric_test(Community)
        temp_test_result["Community_name"]=community_name
        return temp_test_result


