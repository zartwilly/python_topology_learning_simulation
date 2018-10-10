#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 14:37:35 2017

@author: willy

graphe particulier:
    graphe ou chaque sommet est etiquette a -1
"""
import pandas as pd
import numpy as np
import re
import fonctions_auxiliaires as fct_aux
import generations_mesures as mesures
import verif_correl as VerifCorrel
import time
import clique_max as clique
import networkx as nx
import decouverte_cliques as decouvClique
import correction_linegraph as corr_lineGraph
import simulation_50_graphes_PARALLELE as simu
#import matplotlib.pyplot as plt
#import pylab as plab
import genererMatA as geneMatA
from random import randint as valAleotoire
#import matplotlib.pyplot as plt
from pathlib import Path    # http://stackoverflow.com/questions/6004073/how-can-i-create-directories-recursively
import itertools as it
import multiprocessing as mp;
from multiprocessing import Pool;


def matriceE_particuliere():
    """
    a faire: generer automatiquement un graphe
    """
    """
                ["A","B","C","D","E","F","G","H","I","J","K","L"]
    dico = {"A":["D","F"],\
            "B":["J","L"],\
            "C":["G","K"],\
            "D":["A","E","G"],\
            "E":["D","F","H","I"],\
            "F":["A","E","J"],\
            "G":["D","H","C"],\
            "H":["E","G","I","K"],\
            "I":["H","E","J","L"],\
            "J":["F","I","B"],\
            "K":["C","H","L"],\
            "L":["K","I","B"]\
            }
    """
    matE = pd.read_csv("graphe_particulier/matE_particulier.csv", \
                       index_col = "nodes")
    matE.fillna( 0, inplace = True)
    return matE
    
def gamma_noeud(liste_aretes):
    """
    """
    liste_noeuds = [ [arete[0], arete[1]] for arete in liste_aretes]
    set_noeuds = set([item for sublist in liste_noeuds for item in sublist])
    dico = dict()
    for noeud in set_noeuds:
        ens = set()
        for arc in liste_aretes:
            if noeud == arc[0]:
                ens.add( arc[1] )
            if noeud == arc[1]:
                ens.add( arc[0] )
        dico[noeud] = [len(ens), ens]
    return dico   
    
def degre_noeuds_k(dico,k):
    """
    retourne la liste des noeuds dont leur degre est egale a k 
    ds notre cas k = 2
    """
    set_noeuds = set()
    for noeud, caract in dico.items():
        if caract[0] == k:
            set_noeuds.add(noeud)
    return set_noeuds
    
def G0_k(matE, k_max):
    """
    ajoute des aretes entre les noeuds de degre = 2 ===> k = 0, puis
    sur ces aretes crees on ajoute 3 noeuds (les noeuds sont adjacents entre eux)
    """
    liste_aretes = fct_aux.liste_arcs(matE)
    cpt_noeuds = len(matE.columns.tolist())
    for k in range(0,k_max+1):
        dico_gamma_noeud = gamma_noeud(liste_aretes)
        set_noeuds_degre_2 = degre_noeuds_k(dico_gamma_noeud, 2)
        aretes = list(it.combinations(set_noeuds_degre_2, 2))
        if k == k_max:
            liste_aretes.extend(aretes)
        else:
            for arete in aretes :
                noeud = str(cpt_noeuds+1)
                liste_aretes.extend([(arete[0],noeud),(noeud,arete[1])])
                cpt_noeuds += 1
            pass
    return liste_aretes
    pass
    
def test_correction_G0_k(C, liste_aretes_init, dico_cliq):
    
    liste_noeuds_1 = decouvClique.liste_noeuds_1(liste_aretes_init, dico_cliq, ascendant_1 = 1)
    ens_C, dico_cliq, liste_arcs, min_cliques_aSupprDe_ens_C = \
        corr_lineGraph.correction_noeuds(liste_noeuds_1, C, liste_aretes_init, dico_cliq)
    print("ens_C = ",ens_C)
    
    l_aretes_ensC = simu.aretes_C(ens_C)
    s_noeuds_init = { arete[0] for arete in liste_aretes_init}
    s_noeuds_init = s_noeuds_init.union({arete[1] for arete in liste_aretes_init})   
    l_noeuds_init = list(s_noeuds_init)
    
    L_G = fct_aux.transform_listeArcs_mat_adjacence( l_noeuds_init, l_aretes_ensC, oriente=False)
    g0_k = fct_aux.transform_listeArcs_mat_adjacence( l_noeuds_init, liste_aretes_init, oriente=False)
    
    dist_line = None
    liste_arc_diff = None
    dist_line, liste_arc_diff = simu.distance_hamming(g0_k, L_G)
    print("dist_line = ",dist_line)
    print("liste_arc_diff = ",liste_arc_diff)
    
#### simulation suppression k aretes sur le graphe particulier G0_k ====> DEBUT     
def simulation_G0_k(matE_G0_k, k, alpha_min, identifiant, arg_params):
    """
    but: tester la modification de correlation sur le graphe particulier G0_k 
         selon les methodes suivantes
            * le degre min avec permutation 
            * le cout min avec permutation 
            * aleatoire N = 100
    
    matE_G0_k : k etant la profondeur du graphe G0, matE_G0_k est une MATRICE D'ADJACENCE de G0_k
    
    arg_params = {"number_permutations_nodes_1": 10(30, 100), "biais": True, "algoGreedy":False, \
                  "mode_select_noeuds_1":"coutMin" or "degreMin" or "aleatoire", "number_items_pi1_pi2" = 1,\
                  "methode_deleted_add_edges": 0, "SEUIL_PROBA": 0.8, \
                  "proba_seuil": proba_seuil, \
                  "coef_fct_cout":(exposant, facteur_multiplicatif)}
    """
    
    # fichier de tracking or debug
    headers_df = ["G_cpt", "k", "alpha", "nbre_aretes_matE", "nbre_aretes_matE_k_alpha", "deleted_edges",\
                  "nbre_aretes_L_G", "nbre_aretes_diff_matE_k_alpha_LG",\
                  "dist_line", "aretes_diff_matE_k_alpha_LG",\
                  "nbre_aretes_diff_matE_LG", "hamming", "aretes_diff_matE_LG",\
                  "C","som_cout_min","noeuds_corriges",\
                  "min_hamming","mean_hamming","max_hamming","ecart_type",\
                  "max_cout","max_permutation",\
                  "dico_som_min_permutations","dico_dual_arc_sommet","ordre_noeuds_traites","C_old"]
                  
    #  creation du repertoire de calcul selon methode
    path_distr_chemin = str(arg_params["mode_select_noeuds_1"])+"_particulier/"+"data_p_"+\
                            str(arg_params["proba_seuil"])+"/distribution/"
    path_distr = Path(path_distr_chemin)
    path_distr.mkdir(parents=True, exist_ok=True)
    
    df_debug = pd.DataFrame( columns = headers_df)
    cpt_df_debug = 0;
    
    G_cpt = "G_"+str(k)+"_"+str(alpha_min);
    
    # creation repertoire contenant dataset et matrices
    # exple rep = methode_correction_nodes_1/data_p_XX/G_10 avec 10 = cpt_graphe_genere
    path = Path(str(arg_params["mode_select_noeuds_1"])+"_particulier/"+\
                            "data_p_"+str(arg_params["proba_seuil"])+"/"+G_cpt+'/datasets/');
    path.mkdir(parents=True, exist_ok=True);
    path = Path(str(arg_params["mode_select_noeuds_1"])+"_particulier/"+\
                            "data_p_"+str(arg_params["proba_seuil"])+"/"+G_cpt+'/matrices/');
    path.mkdir(parents=True, exist_ok=True);
    
                                  
    # initialisation variables
    aretes_matE = len(fct_aux.liste_arcs(matE_G0_k))
    moy_distline = 0; moy_hamming = 0; sum_distline = 0; sum_hamming = 0; correl_dl_dh = 0;
#    chemin_datasets = str(arg_params["mode_select_noeuds_1"])+"/"+\
#                          "data_p_"+str(arg_params["proba_seuil"])+"/"+G_cpt+"/datasets/";
#    chemin_matrices = str(arg_params["mode_select_noeuds_1"])+"/"+\
#                                  "data_p_"+str(arg_params["proba_seuil"])+"/"+G_cpt+"/matrices/";
    #### A EFFACER SI CA MARCHE
    
    try :
        # modification k correlations
        # TODO modifier modif_k_cases tel que deleted_edge ne se repete pas ==> UN PEU COMPLIQUE car process independant
        matE_k_alpha, dico_deleted_add_edges = simu.modif_k_cases(matE_G0_k.copy(), k, \
                                                             arg_params["methode_delete_add_edges"], 
                                                             arg_params["proba_seuil"])
        deleted_edges = list(dico_deleted_add_edges.values())
        dico_proba_cases = simu.ajouter_proba_matE(matE_k_alpha, dico_deleted_add_edges, arg_params["SEUIL_PROBA"])

        # cliques decoulant de l'algo de couverture
        liste_cliques = list(); dico_cliq = dict(); 
        ordre_noeuds_traites = [] # car liste_cliques = []
        for noeud in matE_k_alpha.columns.tolist():
            dico_cliq[noeud] = -1
        aretes_matE_alpha = fct_aux.liste_arcs(matE_k_alpha)
        dico_gamma_noeud = fct_aux.gamma_noeud(matE_k_alpha, aretes_matE_alpha) 
            
        # algo de couverture selon methodes (arg_params["mode_select_noeuds_1"])
        dico_permutations = dict();
        dico_permutations = decouvClique.solution_methode_nodes_1(dico_gamma_noeud,\
                             liste_cliques, aretes_matE_alpha, ordre_noeuds_traites, \
                             dico_cliq, dico_proba_cases, arg_params)
        
        # Debut selection de la permutation de noeuds dont la distance hamming est la plus petite
        dico_sol = dict()
        dico_sol = simu.best_permutation(dico_permutations, matE_G0_k, matE_k_alpha)
        # FIN selection de la permutation de noeuds dont la distance hamming est la plus petite
        
        
        # moyenner dist_line et hamming pour k aretes supprimes
        moy_distline = dico_sol["dist_line"]; moy_hamming = dico_sol["hamming"] ;
        if moy_hamming == 0 and moy_distline == 0:
            correl_dl_dh = 1
        else:
            correl_dl_dh = abs(moy_hamming - moy_distline)/max(moy_hamming, moy_distline)
        #print("ici")
        
        # ecrire dans un fichier pouvant etre lu pendant qu'il continue d'etre ecrit
        f = open(path_distr_chemin+"distribution_moyDistLine_moyHamming_k_"+str(k)+".txt","a")
        f.write(G_cpt+";"+str(k)+";"+str(moy_distline)+";"+str(moy_hamming)+";"+str(aretes_matE)+\
                ";"+str(correl_dl_dh)+"\n")
        f.close();
        
        
        # pour debug, log, .....
        dico_som_min_permutations = dict();
        for l_noeuds_1, values in dico_permutations.items():
            if values[6] not in dico_som_min_permutations.keys():
                dico_som_min_permutations[values[6]] = [l_noeuds_1]
            else:
                dico_som_min_permutations[values[6]].append(l_noeuds_1)
        
        dico_dual_arc_sommet = mesures.nommage_arcs( matE_G0_k )        
        df_debug.loc[len(df_debug.index)] = [G_cpt, k, alpha_min, \
                        dico_sol["nbre_aretes_matE"], dico_sol["nbre_aretes_matE_k_alpha"], \
                        deleted_edges,\
                        dico_sol["nbre_aretes_LG"], dico_sol["nbre_aretes_diff_matE_k_alpha_LG"],\
                        dico_sol["dist_line"], dico_sol["liste_aretes_diff_matE_k_alpha_LG"], \
                        dico_sol["nbre_aretes_diff_matE_LG"], dico_sol["hamming"],\
                        dico_sol["liste_aretes_diff_matE_LG"],\
                        dico_sol["C"], dico_sol["som_cout_min"], \
                        dico_sol["noeuds_corriges"], \
                        dico_sol["min_hamming"],dico_sol["mean_hamming"], \
                        dico_sol["max_hamming"],dico_sol["ecart_type"],\
                        dico_sol["max_cout"], dico_sol["max_permutation"],\
                        dico_som_min_permutations, dico_dual_arc_sommet,\
                        dico_sol["ordre_noeuds_traites"], dico_sol["C_old"]]
        if cpt_df_debug % 100 == 0:
            simu.save_df(df_debug, path_distr_chemin, identifiant, headers_df)
            df_debug = pd.DataFrame( columns = headers_df)
            print("save %s fois" %cpt_df_debug)        
    
    except Exception as e:
        print("####### EmptyDataError ", G_cpt, ": e = ", e," ####### ");
        df_debug.loc[len(df_debug.index)] = [G_cpt, k, alpha_min, \
                        "error", "error", \
                        "error",\
                        "error","error" ,\
                        "error","error" , \
                        "error","error" ,\
                        "error",\
                        "error", "error", \
                        "error", \
                        "error","error", "error",\
                        "error",\
                        "error","error" ,\
                        "error","error","error", "error"]
    
    pass

#### simulation suppression k aretes sur le graphe particulier G0_k ====> FIN     

if __name__ == '__main__':
    
    start= time.time()
    k_deep = 1#0
    matE_G0_k = matriceE_particuliere()
    noeuds_G0_k = matE_G0_k.columns.tolist()
    liste_aretes_G0_k = [];
    liste_aretes_G0_k = G0_k(matE_G0_k, k_deep)
    print("len  = ", len(liste_aretes_G0_k))
    matE_G0_1 = fct_aux.transform_list_matAdj(liste_aretes_G0_k)
    print("sommets = ", len(matE_G0_1.columns.tolist()))
    
    k_deep = 1#0
    matE_G0_k = matriceE_particuliere();
    liste_aretes_G0_k = [];
    if k_deep == 0:
        liste_aretes_G0_k = G0_k(matE_G0_k, k_deep)
    else:
        liste_aretes_G0_k = G0_k(matE_G0_k, k_deep)
        matE_G0_k = fct_aux.transform_list_matAdj(liste_aretes_G0_k)
    print("liste_aretes = ", len(liste_aretes_G0_k))
    dico_cliq = dict()
    for noeud in matE_G0_k.columns.tolist():
        dico_cliq[noeud] = -1
    
    print("len dico_cliq = ",len(dico_cliq) )
#    test_correction_G0_k( [], liste_aretes_G0_k, dico_cliq)
#    simulation_G0_k(matE_G0_k, k, 1, 1, arg_params)
#
#    
    # test simulation_G0_k
    k_max = 10 
    k_min_range = list(range(k_max)); N = len(k_min_range)
    methode_delete_add_edges = 0; # 1:selection aretes a delete or a ajouter par liste, 0: par proba
    proba_seuil = 0.5;
    SEUIL_PROBA = 0.8; # je pense cest pour ajouter des probas a chaque arete
    algoGreedy = False #True; # False correction tous les noeuds a -1, True: algo Greedy   
    biais = False; # ca ne sert a rien ICI
    puissance = 1;facteur_multiplicatif = 1; # 0: fct_unitaire, 1:fct_normal, 2: fct_carre (2,10), 4:fct_quadratique (10,4)
    coef_fct_cout = (puissance, facteur_multiplicatif)
    mode_select_noeuds_1 = "degreMin" #"coutMin" # "degreMin" # aleatoire
    arg_params = {"number_items_pi1_pi2": 1,"number_permutations_nodes_1": 10, \
                   "methode_delete_add_edges": methode_delete_add_edges, \
                   "biais": biais, \
                   "proba_seuil": proba_seuil,\
                   "SEUIL_PROBA": SEUIL_PROBA,\
                   "algoGreedy":algoGreedy, \
                   "mode_select_noeuds_1":mode_select_noeuds_1,
                   "coef_fct_cout":coef_fct_cout};
                   
##    simulation_G0_k(matE_G0_k, 1, 1, 1, arg_params)

    alpha_max = 2; alpha_range = list(range(alpha_max));                
    probas = [0.0,0.3,0.5,0.8,1.0] # [0.0,0.3]
    methodes_1 = ["degreMin","coutMin","aleatoire"]
    probas = [0.0] ; methodes_1 = ["degreMin"]
    params = list(); 
    for proba_seuil in probas:
        for mode_select_noeuds_1 in methodes_1:
            for k in range(k_max):
                arg_params["proba_seuil"] = proba_seuil;
                arg_params["mode_select_noeuds_1"] = mode_select_noeuds_1;
                arg_params_ = arg_params.copy()
    
                params_ = list(zip( [matE_G0_k]*alpha_max, [k]*alpha_max, \
                            alpha_range, [k]*alpha_max,[arg_params_]*alpha_max))
#                params_ = list(zip( [k]*alpha_max, \
#                            alpha_range, [arg_params_]*alpha_max))
                params.append(params_) 
         
     
    params = [par for tup in params for par in tup]
    print('params = ',len(params))
    p = Pool(mp.cpu_count()-1) 
    p.starmap(simulation_G0_k, params)
    p.terminate()
#    simulation_G0_k(matE_G0_k, k, alpha_min, identifiant, arg_params)
    
    k_min = "0_9";
#    g = open("tempsExecution_"+str(k_min)+".txt","w")
    ti = time.time() - start
#    g.write(str(ti))
#    g.close()
    
    print (time.time() - start)

