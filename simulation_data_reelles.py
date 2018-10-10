#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 13:33:04 2017

@author: willy
 simulation sur donnees reelles: champlan et velizy
"""


import re
import time
import json;
import os;
import io;

import pandas as pd
import numpy as np
import simulation_50_graphes_PARALLELE as simu50_PARALL;
import creation_datasets as create_datasets;

#import determination_matrice_adjacence as det_mat_adj

import fonctions_auxiliaires as fct_aux
import generations_mesures as mesures
import verif_correl as VerifCorrel

import clique_max as clique
import networkx as nx
import decouverte_cliques as decouvClique
#import construction_DAG as cons_DAG;

import matplotlib.pyplot as plt
#import pylab as plab
import genererMatA as geneMatA
from random import randint as valAleotoire
from pathlib import Path    # http://stackoverflow.com/questions/6004073/how-can-i-create-directories-recursively
import itertools as it
from multiprocessing import Pool
import multiprocessing as mp
import extractGraphReel as graphReel
import correlation_timeSeries as correlation;

import test_correlation_shapelet_znorm as corr_shapelet;
import distribution_correlation as distrib_correl  #  a commenter pour execution serveur

###### connectivite de graphes #######
def is_connected(gdict, vertices_encountered=set(), start_vertex=None):
    """ determines if the graph is connected """  
    vertices = list(gdict.keys()) # "list" necessary in Python 3 
    if not start_vertex:
        # chosse a vertex from graph as a starting point
        start_vertex = vertices[0]
    vertices_encountered.add(start_vertex)
    if len(vertices_encountered) != len(vertices):
        for vertex in gdict[start_vertex]:
            if vertex not in vertices_encountered:
                if is_connected(gdict, vertices_encountered, vertex):
                    return True
    else:
        return True
    return False
    
def adjacency(aretes_predit, seuil):
    nodes = set()
    for aretes in aretes_predit:
        nodes.add(aretes[0]); nodes.add(aretes[1]);
    dico_seuil = dict();
    for node in nodes:
        dico_seuil[node] = [];
    for node in nodes:
        for aretes in aretes_predit:
            if node == aretes[0]:
                dico_seuil[node].append(aretes[1])
            if node == aretes[1]:
                dico_seuil[node].append(aretes[0]);
    boolean = is_connected(dico_seuil)
    print("seuil={}, bool={}".format(seuil, boolean))
    return boolean, dico_seuil;
###### connectivite de graphes #######

def comparaison_matEReel_matEPredit(matEReel, metriques, path_root, file_matE_reel):
    """
    but: comparer les aretes differents entre matE reelles et matE predit. en un mot, 
        compter le nombre d'aretes differentes.
        
    path_root = arg_chemins["chemin_equipements"] 
    """
    dico_metrique = dict()
    for metrique in metriques:
        if matEReel.empty:
            matEReel = pd.read_csv(path_root+file_matE_reel+".csv", index_col = "Unnamed: 0")
        edges_reel = fct_aux.liste_arcs(matEReel);
        
        dico = json.load(open(path_root+"seuil_aretes_LG_"+metrique+".json","r"));
        dico_result = dict();
        for seuil, aretes_LG_predit in dico.items():
            dico = dict()
            nbre_edges_diff=np.nan; edges_diff = [];
            aretes_LG_predit = [(arete[0],arete[1]) for arete in aretes_LG_predit];
            nbre_edges_diff, edges_diff = simu50_PARALL.distance_hamming(edges_reel, aretes_LG_predit)
            dico["aretes_communes"] = len( set(edges_reel).intersection(set(aretes_LG_predit)) );
            dico["aretes_predites"] = len(aretes_LG_predit);
            dico["aretes_reel"] = len(edges_reel)
            dico["aretes_differentes"] = nbre_edges_diff;
            dico_result[seuil] = dico;
            bool_seuil, dico_adj = adjacency(aretes_LG_predit, seuil);
        print("dico_result = ", dico_result)
        dico_metrique[metrique] = dico_result
    return dico_metrique;
    pass
def nommage_aretes( matE ):
    """
    return un nom et une liste de valeurs de correlation vide  pour 
    chaque correlation (case de matE) 
    dico = {(arc1, arc2):[nomAleatoire, [divers probas]]}
    """
    arcs = matE.columns.tolist();
    cpt = 0; dico_aretes = dict();
    for i in range(len(arcs)):
        for j in range(i+1, len(arcs)):
            dico_aretes[(arcs[i],arcs[j])] = [str(cpt), []];
            cpt += 1;
    ### pour des tests ###
#    for arc1 in arcs:
#        for arc2 in arcs:
#            if arc1 != arc2 and (arc2,arc1) not in dico_aretes :
#                dico_aretes[(arc1,arc2)] = [str(cpt), []];
#                cpt += 1;
    ### pour des tests ###
    return dico_aretes;
    pass
def distribution_matrice_correlation(chemin_datasets, chemin_matrices, matE = None):
    """
    return la distribution des valeurs probabilistes formant la matrice d'adjacence
    
    Methode: lire tous les dataframes contenant les probabilites de correlations 
            par grandeur (i.e matrice_adjacence_proba)
            puis par case dans la matrice MatE, 
            calculer la probabilite en moyennant tous les correlations par grandeur
    """
    grandeurs = fct_aux.liste_grandeurs(chemin_datasets)
    dico_aretes = dict();
    if matE == None:
        matE = pd.read_csv(chemin_matrices+"matE_notLineGraph.csv", index_col = "Unnamed: 0")
        dico_aretes = nommage_aretes( matE )
        # parcourir chaque grandeur
#        for grandeur in grandeurs:
#            matE_gp = pd.read_csv(chemin_matrices+"matrice_adjacence_proba_"+grandeur+".csv", index_col = "Unnamed: 0")
#            for arete, correl in dico_aretes.items():
#                if arete[0] in matE_gp.columns.tolist() and arete[1] in matE_gp.columns.tolist():
#                    correl[1].append(matE_gp.loc[arete[0]][arete[1]])
#                    print("arete = ", correl[0], " correl = ", correl[1])
#                else:
#                    correl[1].append(0)
        ### pour des tests ###
        cptTest = 0;cptTest1 = 0;
        for grandeur in grandeurs:
            matE_gp = pd.read_csv(chemin_matrices+"matrice_adjacence_proba_"+grandeur+".csv", index_col = "Unnamed: 0")
            for arete, correl in dico_aretes.items():
                if set([arete[0],arete[1]]).issubset(matE_gp.columns):
                    correl[1].append(matE_gp.loc[arete[0]][arete[1]])
                    cptTest1 += 1;
                else:
                    cptTest += 1;
                    correl[1].append(0)
        print("cptTest = ",cptTest," cptTest1 = ",cptTest1)
        ### pour des tests ###
    else:
        dico_aretes = nommage_aretes( matE )
        # parcourir chaque grandeur
        for grandeur in grandeurs:
            matE_gp = pd.read_csv(chemin_matrices+"matrice_adjacence_proba_"+grandeur+".csv", index_col = "Unnamed: 0")
            for arete, correl in dico_aretes.items():
                if arete[0] in matE_gp.columns.tolist() and arete[1] in matE_gp.columns.tolist():
                    correl[1].append(matE_gp.loc[arete[0]][arete[1]])
                else:
                    correl[1].append(0)
    
    # moyenner les differentes correlations
    dico_aretes_correl = dict();
    for arete, correls in dico_aretes.items():
        dico_aretes_correl[arete] = [correls[0], sum(correls[1])/len(correls[1])]
        
    return dico_aretes_correl;

def distribution_correlation(chemin_datasets, chemin_matrices, matE = None):
    """
    a revoir car min_correl  == max_correl ====> bizarre
    """
    dico_aretes_correl = distribution_matrice_correlation(chemin_datasets, chemin_matrices, matE)
    values = [v[1] for k,v in dico_aretes_correl.items()]
    min_correl = min(values);
    max_correl = max(values);
    print("min:",min_correl," max ",max_correl)
    num_bins = np.linspace(min_correl, max_correl, 100); # peut etre 1000
    histo, bin_edges = np.histogram(values, num_bins)
    cdf = np.cumsum(histo*np.diff(bin_edges))
    
    #plot cdf
    plt.plot(values,cdf)
    pass

    ################# debut construction de matrice d'adjacence avec proba et valeurs binaires selon une methode de correlation ########
    ####
    ################# fin construction de matrice d'adjacence avec proba et valeurs binaires selon une methode de correlation ########
    
##### debut simulation sur donnees reelles avec divers SEUILS #######
def matrice_binaire(matE_proba, correl_seuil):
    """
    pour chaque correl (X) de matE_proba :
        * si X < correl_seuil ==> X = 0;
        * si X >= correl_seuil ==> X = 1;
    """
    print("correl_seuil={}, type={}, \nmatE_proba={}\n".format(correl_seuil,type(correl_seuil),matE_proba));
    matE_proba[matE_proba < correl_seuil] = 0;
    print("inf matE_proba={}\n".format(matE_proba))
    matE_proba[matE_proba >= correl_seuil] = 1;
    print("sup matE_proba={}\n".format(matE_proba))
    return matE_proba;

def simulation_reel(cpt_SEUIL, arg_params, arg_chemins):
    """
    SEUIL = definit sur arg_params
    """
    path_distr_chemin = arg_chemins["file"]+"/"+str(arg_params["mode_select_noeuds_1"])+"/"+\
                        "data_correl_seuil_"+str(arg_params["correl_seuil"])+"_"+\
                        arg_params["metrique_distance"]+"/distribution/";
    path_distr = Path(path_distr_chemin)
    path_distr.mkdir(parents=True, exist_ok=True);
    headers_df = ["correl_seuil", "nbre_aretes_matE", "nbre_aretes_LG", \
                  "dist_line", "nbre_aretes_diff_matE_LG", \
                  "liste_aretes_diff_matE_LG", "C_old", "C", "som_cout_min", \
                  "noeuds_corriges", "min_hamming","mean_hamming", \
                  "max_hamming","ecart_type", "max_cout", "max_permutation", \
                  "dico_som_min_permutations", "dico_dual_arc_sommet"];
    
    matE_proba = pd.read_csv(arg_chemins["chemin_matrices"]+arg_params["matE"]+\
                             "_"+arg_params["metrique_distance"]+".csv",index_col = "Unnamed: 0");
    matE = matrice_binaire(matE_proba.copy(), arg_params["correl_seuil"]);
    print("ici ok 1")
#    return 
    ## dico des correls
    dico_proba_cases = dict(); dico_dual_arc_sommet = dict(); cpt_arete = 0;
    for row, col in fct_aux.range_2d(matE_proba.columns.tolist()):
        dico_proba_cases[(row,col)] = matE_proba.loc[row][col];
        dico_dual_arc_sommet[str(cpt_arete)] = (row, col);
        cpt_arete += 1;
    print("ici ok 2")
        
    dico_permutation_cliq = dict();
    #algo corrigeant tous les noeuds a -1
    dico_permutation_cliq = \
        decouvClique.decouverte_cliques(matE.copy(), dico_dual_arc_sommet, \
                                arg_params["seuil_U"], arg_params["epsilon"], \
                                arg_chemins["chemin_datasets"], arg_chemins["chemin_matrices"],\
                                arg_params["ascendant_1"], arg_params["simulation"],\
                                dico_proba_cases,\
                                arg_params);
    print("ici ok 3")                 
    # Debut selection de la permutation de noeuds dont la distance hamming est la plus petite
    dico_sol = dict();
    dico_sol = simu50_PARALL.best_permutation(dico_permutation_cliq, matE, matE)
    # FIN selection de la permutation de noeuds dont la distance hamming est la plus petite
    print("ici ok 4")
    
    dico_som_min_permutations = dict();
    for l_noeuds_1, values in dico_permutation_cliq.items():
        if values[6] not in dico_som_min_permutations.keys():
            dico_som_min_permutations[values[6]] = [l_noeuds_1]
        else:
            dico_som_min_permutations[values[6]].append(l_noeuds_1)
    
    print("ici ok 5")
    df = pd.DataFrame( columns = headers_df);
    df.loc[cpt_SEUIL] = [arg_params["correl_seuil"], \
                         dico_sol["nbre_aretes_matE"], \
                         dico_sol["nbre_aretes_LG"], \
                         dico_sol["dist_line"], \
                         dico_sol["nbre_aretes_diff_matE_LG"], \
                         dico_sol["liste_aretes_diff_matE_LG"], \
                         dico_sol["C_old"], dico_sol["C"], \
                         dico_sol["som_cout_min"], \
                         dico_sol["noeuds_corriges"], \
                         dico_sol["min_hamming"],dico_sol["mean_hamming"], \
                         dico_sol["max_hamming"],dico_sol["ecart_type"], \
                         dico_sol["max_cout"], dico_sol["max_permutation"], \
                         dico_som_min_permutations, dico_dual_arc_sommet \
                         ]
    simu50_PARALL.save_df(df, path_distr_chemin, cpt_SEUIL, headers_df)
    print("save df pour seuil = %s " %cpt_SEUIL);
    
    f = open(path_distr_chemin+"distribution_DistLine_seuil.txt","a")
    f.write(str(cpt_SEUIL)+";"+str(dico_sol["dist_line"])+";"+ \
            str(dico_sol["som_cout_min"])+";"+\
            str(dico_sol["nbre_aretes_diff_matE_LG"])+";"+\
            str(dico_sol["nbre_aretes_matE"])+"\n")
    
    ## fichier seuil, liste aretes_LG
    PATH = arg_chemins["chemin_equipements"] +"seuil_aretes_LG_"+arg_params["metrique_distance"]+".json";
    dico_json = dict();
    dico_correl = dict(); dico_json = dict();
    dico_correl[arg_params["correl_seuil"]] = dico_sol["aretes_LG"];
    if os.path.isfile(PATH) and os.access(PATH, os.R_OK):
        print(" ok seuil ", arg_params["correl_seuil"])
        dico_json = json.load(open(PATH));
    else:
        with open(PATH, 'w') as f:
            json.dump(dico_correl, f)
    dico_json.update(dico_correl)
    with open(PATH, 'w') as f:
            json.dump(dico_json, f)
        
#        return dico_sol["C"], dico_sol["dist_line"], dico_sol["som_cout_min"];
    pass
    ##### FIN simulation sur donnees reelles avec divers SEUILS #######
if __name__ == '__main__':
    
    start= time.time();
    methode_delete_add_edges = 0 # 1:selection aretes a delete or a ajouter par liste, 0: par proba
    correl_seuil = 0.7 # seuil a partir duquel on a des correlations fausses positives (0->1) et fausses negatives (1->0)
                        # [ signification des correlations fausses {positives,negatives} a partir de MATE]
    SEUIL_PROBA = 0.8; # je pense cest pour ajouter des probas a chaque arete
    algoGreedy = False #True; # False correction tous les noeuds a -1, True: algo Greedy 
    critere_selection_pi1_pi2 = 0; # 0: moins de modif,1: ajout aretes> supp aretes, 2:ajout aretes < supp aretes,
    number_permutations_nodes_1 = 100;
    number_items_pi1_pi2 = 1
    seuil_U = 0; nbre_ts = 10000; effet_joule = 0; epsilon = 0.75;
    simulation = False; ascendant_1 = True; biais = False;
    facteur_multiplicatif = 1; exposant = 1; # 0: fct_unitaire, 1:fct_normal, 2: fct_carre, 4:fct_quadratique
    type_fct_cout = "lineaire" # ou "cloche" ou "lineaire"
    coef_fct_cout = (exposant, facteur_multiplicatif, type_fct_cout)
    mode_select_noeuds_1 = "aleatoire" #"coutMin" # "degreMin" # aleatoire
    matE = "matE";
    arg_params = {"matE":matE,\
                  "number_items_pi1_pi2": number_items_pi1_pi2,\
                  "number_permutations_nodes_1": number_permutations_nodes_1, \
                   "methode_delete_add_edges": methode_delete_add_edges, \
                   "biais": biais, \
                   "correl_seuil": correl_seuil,\
                   "seuil_U":seuil_U, "effet_joule":effet_joule,\
                   "nbre_ts":nbre_ts, "epsilon":epsilon, "ascendant_1":ascendant_1,\
                   "simulation":simulation, "algoGreedy":algoGreedy, \
                   "mode_select_noeuds_1":mode_select_noeuds_1,\
                   "coef_fct_cout":coef_fct_cout,\
                   "critere_selection_pi1_pi2":critere_selection_pi1_pi2};
    bool_creation_datasets = False;
    sub = True;
    reseau = "champlan"; fichier_mesures = "datasets_"+reseau+".json";
    rep = "champlan_newDataset"; # or champlan or velizy
    root = "/home/willy/topologyLearning/datas/"+rep+"/";
    root = root+"sous_graphe/" if sub else root;
    chemin_datasets = root+"datasets/";
    chemin_matrices = root+"matrices/";
    chemin_equipements = root+"matrices_equipements/";
    args_chemins ={"chemin_datasets":chemin_datasets,\
                   "chemin_matrices":chemin_matrices,\
                   "chemin_equipements":chemin_equipements,\
                   "file":root};    
    path_file_mesures = chemin_datasets+fichier_mesures;
    args={"path_file":path_file_mesures, "datasets":chemin_datasets};
    if bool_creation_datasets:
        create_datasets.transform_dico_to_dataFrame(args)
    
    fichier_reseau="GraphDico_directed_"+reseau+".json"; 
    path_reseauG= chemin_equipements+fichier_reseau;
    args={"path_reseauG":path_reseauG, "chemin_matrices":chemin_matrices};
    matE_reel = graphReel.matrice_adjacence_sommets_from_groovyFile(args)
    print("aretes_reseauG={} ".format(len(set(matE_reel.columns.tolist()))));
    dico=dict();
    dico["aretes_reseauG"] = set(matE_reel.columns.tolist());
    
    modifCorrel = True; matE_reel = None; 
    if modifCorrel:
        matE_reel = pd.read_csv(chemin_matrices+"matE_reelles.csv",index_col = "Unnamed: 0");
    
    methode_correl = "correlation_par_grandeur" #"correlation_par_metrique_fusion";
    metrique_distance = "lb_keogh" #"fdtw_problem"#"lb_keogh" #"sax"#"pearson";
    mode_correlation = "correlationParMorceaux" #"correlationGlissante"#"correlationParMorceaux" # ou liste_mode = ["correlationParMorceaux","correlationGlissante"]
    indice_max_correlation = 0; derivate = False; 
    interval_derivate = 10; epsilon_sax = 0.8;
    fenetre = 50;
    selected_grandeurs = []; 
    if reseau == "velizy":
        selected_grandeurs = ['AvgI1', 'AvgI']
    else:
        selected_grandeurs = ["P","I","U","S","E","U1","U2","U3","I1","I2","I3"]#["P","I"]
        selected_grandeurs = ["P"];
    args = {"methode_correl": methode_correl,"mode_correlation": mode_correlation,\
                    "indice_max_correlation":indice_max_correlation,\
                    "selected_grandeurs": selected_grandeurs,"modifCorrel": modifCorrel,\
                    "equipements": dico["aretes_reseauG"],\
                    "fenetre": fenetre,"derivate": derivate, "epsilon_sax": epsilon_sax,\
                    "interval_deriv":interval_derivate,\
                    "metrique_distance":metrique_distance}
    metrique_distances = ["lb_keogh","pearson","lcs"]
    metrique_distances = ["metrique_wil"]
    len_dim = "len2930"; args_chemins["len_dim"] = len_dim; args_chemins["grandeur"] = "P"
    matE = None;
    for metrique_distance in metrique_distances:
        args["metrique_distance"] = metrique_distance; 
        args["dbg"] = False; args["dbg_0_1"]=False; args["dbg_seuil"]=False;
        args["dbg_ajoutCorrel"] = False;
        if len_dim == "":
            matE = correlation.correlation_ts(chemin_datasets, chemin_matrices,\
                              chemin_equipements, args);
        else:
            matE = corr_shapelet.correlation(args_chemins,matE_reel.columns.tolist())
    
        correl_seuils = [0.3,0.4,0.5,0.6, 0.7, 0.8, 0.9]
#        correl_seuils = [0.6, 0.7, 0.8, 0.9]
#        correl_seuils = [0.3, 0.4, 0.5]
#        correl_seuils = [0.8];
        params = list(); N = 1;
        for correl_seuil in correl_seuils:
            arg_params["correl_seuil"] = correl_seuil;
            arg_params["metrique_distance"] = metrique_distance
            arg_params_ = arg_params.copy()
        
            params_ = list(zip( [correl_seuil]*N, [arg_params_]*N, [args_chemins]*N ))
            params.append(params_);
    
        print("params ",len(params))
        params = [par for tup in params for par in tup]
    
        #### parallele
        p = Pool(mp.cpu_count()-1) 
        p.starmap(simulation_reel, params)
        p.terminate()
    #    #### parallele
    
    # comparaison entre matE_reel et matE_predit
#    dico_metrique = comparaison_matEReel_matEPredit(matE_reel, metrique_distances, args_chemins["chemin_equipements"], "matE")  
    # distribution de line/hamming distance par rapport au seuil
    distrib_correl.distrib_line_hamming_distance_seuil({"dbg":True})
    
#    methode_correlation = "sax" # euclidienne
#    nom_graphe = "Velizy"
#    seuil_proba_to_01 = 0.03 # seuil qui change les probas en valeurs binaires (0,1)
#    epsilon_sax = 0.01;
#    
#    # algo decouverte et correction
#    methode_correction = "aleatoire" # coutMin ou degreMin ou aleatoire
#    

    print (time.time() - start)               