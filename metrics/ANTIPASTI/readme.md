This script calculates ANTIPASTI logKd binding affinity

**Note**: 
1. the environment file, antipasti-env.yml, has been modified to inlucde R-required packages
2. the preprocessing.py has been modified to accomidate flexible filename parsing and masked size
3. ANTIPASTI assume input pdb is chothia renumbered

cd ANTIPASTI
conda env create -f antipasti-env.yml
conda activate antipasti-env
pip install .
