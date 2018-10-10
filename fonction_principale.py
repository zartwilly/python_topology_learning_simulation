#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 18 14:50:54 2017

@author: willy
fonction principale 
data : velizy
"""
import re
import time
import json;
import os;

import pandas as pd
import numpy as np

import determination_matrice_adjacence as det_mat_adj

import fonctions_auxiliaires as fct_aux
import generations_mesures as mesures
import verif_correl as VerifCorrel

import clique_max as clique
import networkx as nx
import decouverte_cliques as decouvClique
import construction_DAG as cons_DAG;

#import matplotlib.pyplot as plt
#import pylab as plab
import genererMatA as geneMatA
from random import randint as valAleotoire
from pathlib import Path    # http://stackoverflow.com/questions/6004073/how-can-i-create-directories-recursively
import itertools as it
from multiprocessing import Pool


"""
lectrure matE
decouverte_cliques + correction 
orientation
distance line
"""


def fonction_principale( args_cliqs ):
    """
    lectrure matE
    decouverte_cliques + correction 
    orientation
    distance line
    
    args_cliqs = tous les parametres en relation avec la decouverte de cliques
    """
    matE = det_mat_adj.matrice_adjacence_globale( args_cliqs["chemin_dataset"],\
                                                  args_cliqs["chemin_matrices"] );
    print("====matE fini ===")
    
    ### decouverte _cliques + correction
    dico_dual_arc_sommet = dict(); dico_cliq = dict(); ens_C = []; liste_aretes = [];
    ens_C, liste_aretes, dico_cliq, som_cout_min = \
        decouvClique.decouverte_cliques( matE, dico_dual_arc_sommet, args_cliqs["seuil_U"], \
                                    args_cliqs["epsilon"],\
                                    args_cliqs["chemin_dataset"], args_cliqs["chemin_matrices"], \
                                    args_cliqs["ascendant_1"] )
    print("====decouverte cliques fini ===")
     
    ### orientation des aretes
    test_booleen = False;
    cliques, dico_orientation, dico_dual_arc_sommet = \
        cons_DAG.construire_DAG(args_cliqs,  matE, ens_C, test_booleen)
    print("==== construction DAG fini ===")
    
    ### save orientation  des aretes
    orientation_json = json.dump(dico_orientation, open(args_cliqs["chemin_matrices"]+"orientation.json", "w"))
    print("==== save  dico_orientation json  fini ===") 
    
    return ens_C, dico_orientation;
    pass

def calcul_distance_hamming_velizy(args_cliqs):
    """
    lecture du fichier matA (correspondant au reseau de velizy)
    lecture du fichier matA_predit (correspondant au reseau de velizy predit)
    calcul de la distance de hamming
    """
    
    pass
if __name__ == '__main__':
    
     start = time.time();
     nom_graphe = "Velizy"
     nom_methode = "sax" #"euclidienne"
     chemin_dataset = "data/dataReels"+nom_graphe+"/datasets/" ;
     chemin_matrices = "data/dataReels"+nom_graphe+"/matrices/" ;
     seuil_sax = 0.03 ;
     epsilon_sax = 0.01 ;
     
     # export data et creation datasets par grandeur ==> debut 
     file = "exportEdges"+nom_graphe+".csv";
     if nom_graphe == "Champlan":
         grandeurs = []; path = "data/dataReels"+nom_graphe+"/";
         det_mat_adj.create_dataset(file, path, grandeurs, nbre_rows = 100000000)
     # export data et creation datasets par grandeur==> FIN
     
     args_cliqs = dict();
     args_cliqs["seuil_U"]= 10; args_cliqs["epsilon"] = 0.75; 
     args_cliqs["ascendant_1"] = True;
     args_cliqs["chemin_dataset"] = chemin_dataset; 
     args_cliqs["chemin_matrices"] = chemin_matrices;
     
     ens_C, dico_orientation = fonction_principale( args_cliqs )