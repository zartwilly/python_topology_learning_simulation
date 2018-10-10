#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 08:28:34 2017

@author: willy
but : construction du DAG 
    * orientation des aretes du graphe non oriente
"""

import re
import json;
import os;
import time;

import pandas as pd
import numpy as np

import fonctions_auxiliaires as fct_aux
#import determination_matrice_d_adjacence as mat_adj
import generations_mesures as mesures
import verif_correl as VerifCorrel

import clique_max as clique
import networkx as nx
import decouverte_cliques as decouvClique
import simulation_50_graphes_PARALLELE as simu50
#import matplotlib.pyplot as plt
#import pylab as plab
import genererMatA as geneMatA
from random import randint as valAleotoire
from pathlib import Path    # http://stackoverflow.com/questions/6004073/how-can-i-create-directories-recursively
import threading
import multiprocessing as mp;
 

def construire_DAG(arguments, matE_calcule, ens_C,test_booleen):
    """
    A partir des cliques triees par ordre croissant, orienter le graphe
    """
    cliques = list();dico_dual_arc_sommet = dict()
    matE = pd.DataFrame(); matA = pd.DataFrame();
    # generer matE
    if (test_booleen == True):
        matE, matA, dico_dual_arc_sommet = simu50.matriceE( arguments["dimMatA"], arguments["nbre_lien"],\
                                    arguments["chemin_matrices"], arguments["chemin_dataset"], \
                                    arguments["nbre_ts"], arguments["epsilon"], \
                                    arguments["effet_joule"], arguments["test"] )
    
        cliques, aretes, dico_cliq, som_cout_min = \
        decouvClique.decouverte_cliques( matE, dico_dual_arc_sommet, arguments["seuil_U"], arguments["epsilon"], \
                                     arguments["chemin_dataset"], arguments["chemin_matrices"], \
                                     arguments["ascendant_1"])
    else:
        cliques = ens_C.copy();
        matE = matE_calcule;
        
    # construction et orientation 
    sorted_cliques = fct_aux.quicksort(cliques)
    dico = VerifCorrel.caracterisation_arc( matE.columns.tolist() )
    matA_predit = VerifCorrel.construire_matA_intelligent(dico, sorted_cliques)
    dico = VerifCorrel.orientation_arete(cliques, arguments["chemin_dataset"], matE, arguments["epsilon"])
    matA_predit.to_csv(arguments["chemin_matrices"]+"matA_predit.csv")
    list_adjacency_matA_predit = from_adjacency_to_list_matrix(matA_predit, arguments["chemin_matrices"])
    
    # representation matA et matA_predit
    if (test_booleen == True):
        aretes_matA = fct_aux.liste_arcs(matA)
        aretes_matA_predit = fct_aux.liste_arcs(matA_predit)
        G_matA = nx.Graph(); G_matA.add_edges_from(aretes_matA)
        G_matA_predit = nx.Graph(); G_matA_predit.add_edges_from(aretes_matA_predit)
    
        plt.figure(1)
        plt.subplot(211); nx.draw(G_matA, with_labels = True)
        plt.subplot(212); nx.draw(G_matA_predit, with_labels = True)
        plt.savefig(arguments["chemin_matrices"]+"matA_predit.png", format= "PNG")
        
    return cliques, dico, dico_dual_arc_sommet
    pass


def from_adjacency_to_list_matrix(matA_predit, chemin_matrices):
    cols  = matA_predit.columns.tolist()
    adjacency_list = dict()
    for node_cols in cols:
        voisins = list()
        for node_ind in cols:
            if matA_predit[node_cols][node_ind] == 1:
                voisins.append(node_ind)
        adjacency_list[node_cols] = voisins
        
    directory = chemin_matrices+str("graphe/")
    if not os.path.exists(directory):
        os.makedirs(directory)
        
    adjacency_json = json.dump(adjacency_list, open(directory+"graphe.json", "w"))
    print("****** adjacency_json = ",adjacency_json)
    
    return adjacency_json
    pass
if __name__ == '__main__':
    
    start = time.time();
    chemin_datasets = "data/datasets/"; chemin_matrices  = "data/matrices/";
    nbre_ts = 10; nbre_lien = (2,5); 
    effet_joule = 0.1; epsilon = 0.75; seuil_U = 10;
    dimMatA = 5; test = "EN TEST"; ascendant_1 = True; 
    arguments = dict()
    arguments["dimMatA"] = dimMatA; arguments["nbre_lien"] = nbre_lien; 
    arguments["chemin_matrices"] = chemin_matrices ;
    arguments["chemin_datasets"] = chemin_datasets; arguments["nbre_ts"] = nbre_ts;
    arguments["epsilon"] = epsilon; arguments["effet_joule"] = effet_joule; 
    arguments["seuil_U"] = seuil_U;
    arguments["test"] = test; arguments["ascendant_1"] = ascendant_1 ;
    #arguments[""] = ;arguments[""] = ;arguments[""] = ;arguments[""] = ;
    
    cliques, dico, dico_dual_arete_sommet = construire_DAG(arguments)
    print("cliques: ", cliques); print("dico = ", dico ); 
    print("dico_dual_arete_sommet = ",  dico_dual_arete_sommet);
    print("running time : ", time.time() - start)
        