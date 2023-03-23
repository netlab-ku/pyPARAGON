# PARAGON 

## Overview  

<p style="text-align: justify;font-weight: normal;"><font size=3> PARAGON (PAgeRAnk-flux on Graphlet-guided network for multi-Omic data integratioN) is a novel method that infers a context specific subnetwork from a given reference network utilizing omic-hits as a initial node set and then identifies significant modules in specific subnetworks. During network inference, Paragon run in three steps:  i. Graphlet-guided network (GGN) construction, ii. Propagation via  Personalized PageRank (PPR) algorithm, iii. Edge scoring and selection via flux calculation. <font> <p> 


    
<img src="Concepts/Network_Inference.png" width="400">

<p style="text-align: justify;font-weight: normal;"><font size=3>Graphlet module of PARAGON identifies, an associated region of reference network through substantially frequent graphlets composed of 2-, 3-, and 4-nodes. Each graphlet includes an intermediate node that is the most interacting nodes in the graphlet. Paragon collects the frequent graphlet motifs into GGN. In this way, PARAGON shrinks the size of the network into GGN and eliminates highly connected nodes and their unrelated interactions. <font> <p>

<p style="text-align: justify;font-weight: normal;"><font size=3>Independent to the graphlet module, Paragon propagates signals from omic-hits through the personalizaed PageRank algorithm in re-scoring all proteins in the reference. Then, flux calculation weights each edge of the GGN, considering PageRank scores of nodes, confidence scores of edges and the number of node interactions. Finally, Paragon infer the context specific network by selecting the highly scored edges. <font> <p>

<img src="Concepts/Community__Analysis.png" width="400">
    
<p style="text-align: justify;font-weight: normal;"><font size=3>Paragon goes beyond network inference by dividing the network into functional units. Community analysis module recruites the louvain community detection method based on network topology. Then, Paragon can identify significant modules by appying hypergeometric test for a given prior knowledge such as biological processes, pathways. <font><p>     

## Installation 
<p style="text-align: justify; font-weight: normal;"> <font size=2>
    
## Citation <h2>      
    
<p style="text-align: justify; font-weight: normal;"> <font size=2>
Arici, M.K., Tuncbag, N. Discovering Hidden Connections in Omics Data: an Integrative Modeling Approach for Unveiling Cancer Networks, 2023, in submission <p>




