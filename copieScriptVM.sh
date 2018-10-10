
scp -r *.py  dcbrain@wil-mig.dcbra.in:/home/dcbrain/simulation_seuils/
ssh dcbrain@wil-mig.dcbra.in; mkdir simulation_seuils/; cd simulation_seuils

mkdir -r data/lineaire_simul50Graphes_priorite_supp/degreMin/
mkdir -r data/lineaire_simul50Graphes_priorite_supp/coutMin/
mkdir -r data/lineaire_simul50Graphes_priorite_supp/aleatoire/

mkdir -r data/lineaire_simul50Graphes_priorite_ajout/degreMin/
mkdir -r data/lineaire_simul50Graphes_priorite_ajout/coutMin/
mkdir -r data/lineaire_simul50Graphes_priorite_ajout/aleatoire/

mkdir -r data/lineaire_simul50Graphes_priorite_aucune/degreMin/
mkdir -r data/lineaire_simul50Graphes_priorite_aucune/coutMin/
mkdir -r data/lineaire_simul50Graphes_priorite_aucune/aleatoire/


screen -S simulation_seuils_150G
export PATH=/home/dcbrain/anaconda3/bin:$PATH
ipython3 simulation_seuil_graphe.py
