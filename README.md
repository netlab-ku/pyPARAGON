# pyPARAGON 

## Overview  

<p align="justify"><font-weight="normal"><font size=3> pyPARAGON (PAgeRAnk-flux on Graphlet-guided network for multi-Omic data integratioN) is a novel method that infers a context specific subnetwork from a given reference network, utilizing omic-hits as a seed node set and then identifies significant modules in specific subnetworks. During network inference, pyPARAGON run in three steps: 

i. Graphlet-guided network (GGN) construction, 
ii. Propagation via  Personalized PageRank (PPR) algorithm, 
iii. Edge scoring and selection via flux calculation (Figure 1A). <font> </p> 



<p align="center">    
<img " src="Concept_Figs/Network_Inference.png" width="400">
<p align="center">
<p align="justify"><font-weight="normal"><font size=> <font-weight="bold">Figure 1:</font-weight> Conceptual view of network inference in pyPARAGON <p>
 
 
 
<p align="justify" font-weight="normal"><font size=3>Graphlet module of pyPARAGON identifies, an associated region of reference network through substantially frequent non-isomorphic graphlets composed of 2-, 3-, and 4-nodes. Each graphlet includes an intermediate node that is the most interacting nodes in the graphlet. pyPARAGON collects the frequent graphlet motifs into GGN. In this way, pyPARAGON shrinks the size of the network into GGN and eliminates highly connected nodes and their unrelated interactions. <font> <p>

<p align="justify" font-weight="normal"><font size=3>Independent to the graphlet module, pyPARAGON propagates signals from omic-hits through the personalizaed PageRank algorithm in re-scoring all proteins in the reference (Figure 1B). Then, flux calculation weights each edge of the GGN, considering PageRank scores of nodes, confidence scores of edges and the number of node interactions. Finally, pyPARAGON infer the context specific network by selecting the highly scored edges (Figure 1B). <font> <p>



<p align="center"> 
<img src="Concept_Figs/Community__Analysis.png" width="400">
<p align="center"> 
<p align="justify"><font-weight="normal"><font size=> <font-weight="bold">Figure 2:</font-weight> Conceptual view of community analysis in pyPARAGON <p>



<p align="justify" font-weight="normal"><font size=3>pyPARAGON goes beyond network inference by dividing the network into functional units. Community analysis module recruites the louvain community detection method based on network topology (Figure 2). Then, pyPARAGON can identify significant modules by appying hypergeometric test for a given prior knowledge such as biological processes, pathways. <font><p>     

pyPARAGON allows researchers to precisely integrate omics data through a reference network, composed of huge prior knowledge. Researchers may model topic of interest, such as diseases, drug trials in inferred context-specific networks. Additionally, independent to network inference module,  pyPARAGON determines the substantial communities of any network associated with the biological annotations for biological interpretations such as patient stratification, survival analysis and personalized medicine.



## Installation from Git :
<p align="justify" font-weight="normal"><font size=3>
pip install git+https://github.com/metunetlab/pyPARAGON.git

                                       
                                       
                               

## Installation in linux command prompt after downloading :
<p align="justify" font-weight="normal"><font size=3>

1. Create a virtual environment for pyPARAGON
python3 -m venv pyPARAGONenv

2. Activate the pyPARAGON Environment 
source venv/bin/activate

3. Install python package for in-house setup
pip install wheel
pip install setuptools

4. Go to pyPARAGON folder, the folder including setup.py and run setup.py as using
python setup.py bdist_wheel

5. The wheel file should be stored in the "dist" folder that should be written in the step 4. run the file extended with "whl" as using
pip install /path/to/wheelfile.whl 


## Citation
<p align="justify" font-weight="normal"><font size=3>
Arici, M.K., Tuncbag, N. Discovering Hidden Connections in Omics Data: an Integrative Modeling Approach for Unveiling Cancer Networks, 2023, in submission

