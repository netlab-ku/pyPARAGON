# PARAGON 

## Overview  

<p align="justify" font-weight="normal"><font size=3> PARAGON (PAgeRAnk-flux on Graphlet-guided network for multi-Omic data integratioN) is a novel method that infers a context specific subnetwork from a given reference network utilizing omic-hits as a initial node set and then identifies significant modules in specific subnetworks. During network inference, PARAGON run in three steps:  i. Graphlet-guided network (GGN) construction, ii. Propagation via  Personalized PageRank (PPR) algorithm, iii. Edge scoring and selection via flux calculation (Figure 1A). <font> </p> 


<p align="center">    
<img " src="Concept_Figs/Network_Inference.png" width="400">
<p align="center">
<p align="justify" font-weight="normal"><font size=2>Figure 1 Conceptual view of network inference in PARAGON
                                       
<p align="justify" font-weight="normal"><font size=3>Graphlet module of PARAGON identifies, an associated region of reference network through substantially frequent graphlets composed of 2-, 3-, and 4-nodes. Each graphlet includes an intermediate node that is the most interacting nodes in the graphlet. PARAGON collects the frequent graphlet motifs into GGN. In this way, PARAGON shrinks the size of the network into GGN and eliminates highly connected nodes and their unrelated interactions. <font> <p>

<p align="justify" font-weight="normal"><font size=3>Independent to the graphlet module, PARAGON propagates signals from omic-hits through the personalizaed PageRank algorithm in re-scoring all proteins in the reference. Then, flux calculation weights each edge of the GGN, considering PageRank scores of nodes, confidence scores of edges and the number of node interactions. Finally, PARAGON infer the context specific network by selecting the highly scored edges. <font> <p>

<p align="center"> 
<img src="Concept_Figs/Community__Analysis.png" width="400">
<p align="center"> 

<p align="justify" font-weight="normal"><font size=3>PARAGON goes beyond network inference by dividing the network into functional units. Community analysis module recruites the louvain community detection method based on network topology. Then, PARAGON can identify significant modules by appying hypergeometric test for a given prior knowledge such as biological processes, pathways. <font><p>     

## Installation 
<p align="justify" font-weight="normal"><font size=2>
    
## Citation <h2>      
    
<p align="justify" font-weight="normal"><font size=2>
Arici, M.K., Tuncbag, N. Discovering Hidden Connections in Omics Data: an Integrative Modeling Approach for Unveiling Cancer Networks, 2023, in submission <p>




