#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 10:49:24 2018

@author: willy
"""
import re 
import os 
import time
import math
import pandas as pd
import numpy as np
import fonctions_auxiliaires as fct_aux
import generations_mesures as mesures
import decouverte_cliques as decouvClique
import genererMatA as geneMatA

#import matplotlib.pyplot as plt
from pathlib import Path    # http://stackoverflow.com/questions/6004073/how-can-i-create-directories-recursively
import itertools as it
from multiprocessing import Pool;
from pandas.parser import CParserError;
import multiprocessing as mp
from scipy.stats import truncnorm
from random import randint as valAleotoire
import random;

def liste_arcs(mat):
    """ retourne la liste des arcs ou aretes d'un graphe. """
    res = list();
    for row, col in fct_aux.range_2d(mat.columns.tolist()):
        if mat.loc[row][col] == 1 or mat.loc[col][row] == 1:
            res.append((row, col))
    return res;
    
def nommage_arcs(mat):
    """ nommage des arcs du graphe du reseau energetique.
    
    elle retourne un dictionnaire avec :
        la cle : le nom de l'arc
        la valeur : le tuple correspondant a l'arc
        
    """
    dico_arcs_sommets = dict();
    for arc in liste_arcs(mat):
        dico_arcs_sommets["_".join(arc)] = arc;
    return dico_arcs_sommets
    
def creer_reseau(chemin_datasets, chemin_matrices, args):
    """creation d'un reseau energetique.
    
    creer un graphe, ajouter des mesures de puissance sur le graphe puis
    extraire son line-graphe
    
    """
    if args["dbg"]:
        chemin_datasets = "dataTestNewCriterecorrection/datasets/";
        chemin_matrices = "dataTestNewCriterecorrection/matrices/";
        path_ = Path(chemin_datasets);
        path_.mkdir(parents=True, exist_ok=True);
        path_ = Path(chemin_matrices);
        path_.mkdir(parents=True, exist_ok=True);
        
        
    dico_graphe = {"a":["b"], "b":["c","d"], "c":["e","f"], 
                   "d":["f"], "e":["g"], "f":["h"], "g":[],
                   "h":[]};
                   
    mat = pd.DataFrame(index = dico_graphe.keys(), columns = dico_graphe.keys());
    
    for k, vals in dico_graphe.items():
        for v in vals:
            mat.loc[k,v] = 1;
    mat.fillna(value=0, inplace=True);
    
    # ajouter mesures 
    grandeurs = ["P"];    
    dico_arcs_sommets = dict()
    dico_arcs_sommets = nommage_arcs(mat)
    
#    genererMesures_all_grandeurs(matA, dico_dual_arc_sommet, liste_grandeurs, location = "data/datasets/", taille = 3, effet_joule = 0):

    mesures.genererMesures_all_grandeurs(mat, dico_arcs_sommets, grandeurs, 
                                 chemin_datasets, taille=100, effet_joule=0.1)

    
    #matrice du linegraphe du reseau de flot a determiner
    arcs = liste_arcs(mat)
    matE = mesures.creation_matE(dico_arcs_sommets, arcs)
    matE.to_csv(chemin_matrices+"matE.csv")
    return matE, mat, dico_arcs_sommets;

def simulation_nouveau_critere(args):
    """ simulation d'un nouveau critere de correction.
    
    """
    matE, mat, dico_arcs_sommets = creer_reseau(args["chemin_datasets"], 
                                                args["chemin_matrices"], args);
                                                
    
    # supprimer ou ajouter des aretes dans matE
    matE_k_alpha = matE.copy()
    rd = random.random()
    dico_k_erreurs = {"aretes_ajoutees":[], "aretes_supprimees":[]}
    if args["k_erreurs"] == 1 and rd >= args["p_correl"]:
        matE_k_alpha.loc["b_d","c_f"] = 1;
        matE_k_alpha.loc["c_f","b_d"] = 1;
        dico_k_erreurs["aretes_ajoutees"].append(("b_d","c_f"));
    elif args["k_erreurs"] == 1 and rd < args["p_correl"]:
        matE_k_alpha.loc["a_b","c_b"] = 0;
        matE_k_alpha.loc["c_b","a_b"] = 0;
        dico_k_erreurs["aretes_supprimees"].append(("a_b","b_c"));
    else:
        aretes_matE = matE_k_alpha.columns.tolist()
        for _ in range(math.ceil(args["k_erreurs"] * args["p_correl"])):         # suppression d'aretes
            irow,row = random.choice(list(enumerate(aretes_matE)));
            aretes_matE.pop(irow);
            icol,col = random.choice(list(enumerate(aretes_matE)));
            aretes_matE.pop(icol);
            
            matE_k_alpha.loc[row,col] = 0;
            matE_k_alpha.loc[col,row] = 0;
            dico_k_erreurs["aretes_supprimees"].append((row,col));
            
        aretes_matE = matE_k_alpha.columns.tolist()
        for _ in range(args["k_erreurs"] - math.ceil(args["k_erreurs"] * 
                                                    args["p_correl"])):         # ajout d'aretes
            irow,row = random.choice(list(enumerate(aretes_matE)));
            aretes_matE.pop(irow);
            icol,col = random.choice(list(enumerate(aretes_matE)));
            aretes_matE.pop(icol);
            
            matE_k_alpha.loc[row,col] = 1;
            matE_k_alpha.loc[col,row] = 1;
            dico_k_erreurs["aretes_ajoutees"].append((row,col));
    
    # algorithme de couverture
    C = list(); aretes_Ec = list();
    dico_cliq = dict(); dico_sommets_par_cliqs = dict();
    C, dico_cliq, aretes_Ec, ordre_noeuds_traites, dico_sommets_par_cliqs = \
    decouvClique.decouverte_cliques(matE_k_alpha, dico_arcs_sommets, \
                                    args["seuil_U"], args["epsilon"], \
                                    args["chemin_datasets"], 
                                    args["chemin_matrices"],
                                    args["ascendant_1"], 
                                    args["simulation"],\
                                    dict(),
                                    args)
    print("aretes_Ec={},\n C={} ,\n sommets_par_cliqs={},\n dico_k_erreurs={},\n dict_cliq={}".\
          format(len(aretes_Ec), C, dico_sommets_par_cliqs, dico_k_erreurs, 
                 dico_cliq))
    return matE_k_alpha

if __name__ == '__main__':
    
     start= time.time();
     args = {"dbg":True};
     chemin_datasets = "dataTestNewCriterecorrection/datasets/";
     chemin_matrices = "dataTestNewCriterecorrection/matrices/";
     
     bool_reseau = False;
     bool_couverture = True
     
     # test generation reseau energetique ==>OK
     if bool_reseau:
         matE, mat, dico_arcs_sommets = creer_reseau(chemin_datasets, 
                                                 chemin_matrices, args);
     
     # valeurs des parametres de decouverte_cliques
     seuil_U = 10;
     epsilon = 0.75;
     ascendant_1 = True;
     simulation = True;
     number_items_pi1_pi2 = 1;
     k_erreurs = 1; 
     p_correl = 0.5;
     args = {"dbg":True, "seuil_U":seuil_U, "epsilon":epsilon, 
             "ascendant_1":ascendant_1, "simulation":simulation,
             "number_items_pi1_pi2":number_items_pi1_pi2,
             "k_erreurs":k_erreurs, 
             "p_correl":p_correl,
             "chemin_datasets":chemin_datasets,
             "chemin_matrices":chemin_matrices}
    
     # test sur couverture en cliques ===> 
     if bool_couverture:
         matE_k_alpha = simulation_nouveau_critere(args)                                           
     