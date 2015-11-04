# Sample-Research-Codes
Codes in various languages written by me for research

:: Project: Ligand Prediction 
Description: given a receptor structure, derived optimum binding features 
and use these features to search for potential binders (ligands)

main code: 
- constellationSearch.pl (find all possible patches on a candidate structure that matches the features)
- clusterGreedy_oneStep_byRMSD.pl	(cluster all matching patches by structural similarity) 

supporting modules:
- clashCheck.pm			
- extractByMatchCode.pm
- loadCoordFiles.pm

:: Project: Binding Interface Anaylsis
Description: downloaded interface analysis of a receptor-ligand complex 
from http://www.ebi.ac.uk/pdbe/pisa as an XML file;
parse file to analyse interface features

- xmlParser_pisa.pl
