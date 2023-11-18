# PRA-Pred
Protein-RNA binding affinity prediction
PRA-Pred predicts protein-RNA binding affinity using structure-based features and classes based on RNA strand, structure and function of protein. We found that RNA base-step parameters, number of atomic contacts between the RNA and protein, interaction energies, contact potentials, and hydrogen bonds are important to understand the binding affinity. PRA-Pred shows a correlation of 0.88 and a mean absolute error(MAE) of 0.71 in Self-consistency, correlation of 0.77 and a MAE of 1.02 kcal/mol in jack-knife test.
This program takes the input of protein-DNA complex and predicts the binding affinity.

Input options: User can input PDB ID of protein-RNA complex or can provide the file in PDB format. * User can enter the chain ID of protein. Multiple chain Ids can be entered with comma. (Eg. Chain ID: A,B) (Optional) User can enter the chain ID of RNA. Multiple chain Ids can be entered with comma. (Optional) User should select one of the RNA class for the prediction. *

1. Single strand
2. Double strand

User should select one of the Functional class of protein for the prediction. *
1. Enzyme
2. Regulatory
3. other 

Dependencies:

Program is basically using python3, and demands few python packages to run.
Please ensure the packages such os,re, Bio, numpy, functools,sys, time, shutil, subprocess, math, uuid, cgitb, timeit, wget, glob, urllib, pandas, warnings are installed.
Make sure to edit the naccess executable file with the path of local installation (Source: http://www.bioinf.manchester.ac.uk/naccess/).
Also needs the installation of 3vvv software, used for volume and surface area calculation (Source: https://github.com/vosslab/vossvolvox)
Kindly install the Foldx software. (Source: https://foldxsuite.crg.eu/)
Execute x3dna_setup and add the path to environment variable (http://forum.x3dna.org/downloads/3dna-download/).
Request and download the 3DNA.
Run:

Run "python3 pdpredict.py"
Enter option for PDB-ID and PDB-file
Provide the PDB ID or PDB-file with chain
Enter the RNA strand and functional class
The program will be executed, and the result will be provided in result.txt In the Form of PDB ID, Chain, Predicted affinity(kcal/mol), dissocation constant(M). You can also access it through the webserver https://web.iitm.ac.in/bioinfo2/prapred/index.html
