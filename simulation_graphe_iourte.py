#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 08:32:26 2018

@author: willy
simulation graphe iourte
"""

import re, os, time, math;
import random
import pandas as pd
import numpy as np
import fonctions_auxiliaires as fct_aux
import decouverte_cliques as decouvClique
import matplotlib.pyplot as plt
import itertools as it
import multiprocessing as mp;
from pathlib import Path    # http://stackoverflow.com/questions/6004073/how-can-i-create-directories-recursively
from multiprocessing import Pool;
import simulation_50_graphes_PARALLELE as simu50;

#CPT_DF_DEBUG
def matriceE_particuliere():
    """
    retourne la matrice d'adjacence du graphe iourte
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
    
# -------- correl aretes matE_G_K    
def aretes(matE, k0_1):
    """
    but: trouver toutes les non aretes de matE cad matE.loc[x][y] == 0
    k0_1 = {0,1}
    """
    liste_cols = matE.columns.tolist()
    tmp = list(); res_aretes = list()
    for row in liste_cols:
        for col in liste_cols:
            if row != col and (col,row) not in tmp and matE.loc[row][col] == k0_1:
                tmp.append( (col,row) ); tmp.append( (row,col) );
                res_aretes.append( (row, col) );
    return res_aretes;
    
def loi_proba(low, high):
    """
    generer une nombre entre low and high selon la loi uniforme
    """
    return random.uniform(low, high);
    
def matrice_correl_G_k(matE_G_k, p_correl, correl_seuil):
    """
    definir la matrice de correlation de G_k
    """
    infs = (0,p_correl-0.2); sups = (correl_seuil,1); dico_probas=dict()
    aretes_0 = aretes(matE_G_k, 0); aretes_1 = aretes(matE_G_k, 1);
    matE_corr_G_k = matE_G_k.copy();
    for arete in aretes_0:
        correl = loi_proba(infs[0], infs[1])
        matE_corr_G_k.loc[arete[0],arete[1]] = correl;
        matE_corr_G_k.loc[arete[1],arete[0]] = correl;
        dico_probas[arete] = correl;
    for arete in aretes_1:
        correl = loi_proba(sups[0], sups[1])
        matE_corr_G_k.loc[arete[0],arete[1]] = correl;
        matE_corr_G_k.loc[arete[1],arete[0]] = correl;
        dico_probas[arete] = correl;
    return matE_corr_G_k, dico_probas;
# -------- correl aretes matE_G_K    
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
    
def G_k(matE_G0, k_deep, p_correl, correl_seuil):
    """
    generer le graphe iourte G_k de k = k_deep profondeur
    
    ajoute des aretes entre les noeuds de degre = 2 ===> k = 0, puis
    sur ces aretes crees on ajoute 3 noeuds (les noeuds sont adjacents entre eux)
    
    A VERIFIER je ne crois pas trop
    """
    liste_aretes_G_k = fct_aux.liste_arcs(matE_G0) 
    cpt_noeuds = len(matE_G0.columns.tolist())
    for k in range(0,k_deep+1):
        dico_gamma_noeud = gamma_noeud(liste_aretes_G_k) # bizarre 
        set_noeuds_degre_2 = degre_noeuds_k(dico_gamma_noeud, 2) # bizarre
        aretes = list(it.combinations(set_noeuds_degre_2, 2)) # bizarre
        if k == k_deep:
            liste_aretes_G_k.extend(aretes)
        else:
            for arete in aretes :
                noeud = str(cpt_noeuds+1)
                liste_aretes_G_k.extend([(arete[0],noeud),(noeud,arete[1])])
                cpt_noeuds += 1;
                
    ## construire la matrice d'adjacence de G_k
    noeuds = list(set([arete for tup in liste_aretes_G_k for arete in tup]))
    matE_G_k = pd.DataFrame( index = noeuds, columns = noeuds);
    for arc in liste_aretes_G_k:
        matE_G_k.loc[ arc[0] ][arc[1]] = 1;
        matE_G_k.loc[ arc[1] ][arc[0]] = 1;
    matE_G_k.fillna(0, inplace=True);
    
    matE_corr_G_k, dico_probas = matrice_correl_G_k(matE_G_k, p_correl, correl_seuil);
    return matE_G_k, liste_aretes_G_k, dico_probas;    

# ----- best permutation iourte ------
def best_permutation_iourte(dico_permutation_cliq, matE_G_k):
    """
    selection de la permutation dont le cout et la distance de Hamming sont minimum
    dico_permutation = [ cle : une permutation
                    (A,B,C,D): [C, Ec, dico_cliq, som_cout_min, noeuds_corriges]
    ]
    return dico_sol 
    dico_sol = {"nbre_aretes_LG":,"aretes_LG":,"nbre_aretes_diff_matE_G_k_LG": ,\
                "dist_line":,"aretes_diff_matE_G_k_LG":,"C":,"som_cout_min":,\
                "noeuds_corriges":,"ordre_noeuds_traites":}
    """
    dico_sol = dict(); som_dist_line = 0
    aretes_matE_G_k = fct_aux.liste_arcs(matE_G_k);                    
    for tup_node_1, solutions in dico_permutation_cliq.items():
        aretes_LG = None;
        aretes_LG = simu50.aretes_C(solutions[0]);
        
        dist_line = None; aretes_diff_matE_G_k_LG = None;
        dist_line, aretes_diff_matE_G_k_LG = \
            simu50.distance_hamming( aretes_matE_G_k, aretes_LG )
        
                
        som_cout = solutions[6]; som_dist_line += dist_line;
#        print(" hamming=", hamming," ordre:", tup_node_1," cout:",som_cout," LG:",len(liste_aretes_LG)," MAtE:",len(liste_aretes_matE))
        if (dist_line, som_cout) not in dico_sol.keys():
            dico_sol[(dist_line, som_cout)] = [{"nbre_aretes_LG": len(aretes_LG ), \
                                "aretes_LG": aretes_LG,\
                                "nbre_aretes_diff_matE_G_k_LG": len(aretes_diff_matE_G_k_LG),\
                                "dist_line": dist_line, \
                                "aretes_diff_matE_G_k_LG": aretes_diff_matE_G_k_LG, \
                                "C": solutions[0], \
                                "som_cout_min":solutions[6], \
                                "noeuds_corriges": tup_node_1 ,\
                                "ordre_noeuds_traites": solutions[5]
                                }]
        else:
            dico_sol[(dist_line, som_cout)].append({"nbre_aretes_LG": len(aretes_LG ), \
                                "aretes_LG": aretes_LG,\
                                "nbre_aretes_diff_matE_G_k_LG": len(aretes_diff_matE_G_k_LG),\
                                "dist_line": dist_line, \
                                "aretes_diff_matE_G_k_alpha_LG": aretes_diff_matE_G_k_LG, \
                                "C": solutions[0], \
                                "som_cout_min":solutions[6], \
                                "noeuds_corriges": tup_node_1,\
                                "ordre_noeuds_traites": solutions[5] }) 
        
    # selection du min et ajout mean, max
#    print("dico_sol.keys = ", dico_sol.keys())
    print("len(dico_permutation_cliq) = ",len(dico_permutation_cliq))
    min_keys_line_cout = min(dico_sol.keys())
    max_keys_line_cout = max(dico_sol.keys())
    mean_line = som_dist_line/len(dico_permutation_cliq)
    ## ecart type --> debut
    ecart_type = 0; som = 0;
    for tup_key, list_dico in dico_sol.items():
        som += len(list_dico) * pow((tup_key[0] - mean_line),2)
    ecart_type = pow( som/len(dico_permutation_cliq), 1/2)
    ## ecart type --> fin
    
    print("----> min_line = ",min_keys_line_cout, " som_line = ", som_dist_line,\
          " mean_line= ", mean_line," len_tupl = ", len(dico_permutation_cliq) )

    dico_sol_min = dict();
    dico_sol_min = { \
            "nbre_aretes_LG": dico_sol[min_keys_line_cout][0]["nbre_aretes_LG"],\
            "aretes_LG": dico_sol[min_keys_line_cout][0]["aretes_LG"],\
            "nbre_aretes_diff_matE_G_k_LG": dico_sol[min_keys_line_cout][0]["nbre_aretes_diff_matE_G_k_LG"],\
            "dist_line": dico_sol[min_keys_line_cout][0]["dist_line"],\
            "aretes_diff_matE_G_k_LG": dico_sol[min_keys_line_cout][0]["aretes_diff_matE_G_k_LG"],\
            "C": dico_sol[min_keys_line_cout][0]["C"], \
            "som_cout_min": min_keys_line_cout[1],\
            "noeuds_corriges": dico_sol[min_keys_line_cout][0]["noeuds_corriges"],\
            "ordre_noeuds_traites": dico_sol[min_keys_line_cout][0]["ordre_noeuds_traites"],\
            "min_line": dico_sol[min_keys_line_cout][0]["dist_line"],\
            "mean_line": mean_line,\
            "max_line": max_keys_line_cout[0],\
            "max_cout": max_keys_line_cout[1],\
            "max_permutation": dico_sol[max_keys_line_cout][0]["noeuds_corriges"], \
            "ecart_type": ecart_type \
    }

    return dico_sol_min;
# ----- best permutation iourte ------

# ----- comparaison_aretes_G_k_LG ------
def comparaison_aretes_G_k_LG(aretes_matE_G_k, aretes_LG):
    """ recherche les aretes de G_k n'appartenant pas aux aretes de aretes_LG """
    cpt_aretes_G_k_notIn_LG= 0;
    for arete_G_k in aretes_matE_G_k:
        if (arete_G_k[0],arete_G_k[1]) not in aretes_LG and \
            (arete_G_k[1],arete_G_k[0]) not in aretes_LG :
            cpt_aretes_G_k_notIn_LG += 1;
    return math.ceil( (cpt_aretes_G_k_notIn_LG/len(aretes_matE_G_k))*100);

def comparaison_aretes_G_k_LG_Debug_old():
    # A EFFACER
    aretes_G_k = fct_aux.liste_arcs(matriceE_particuliere());
    matE_G_2,aretes_G_k,dico = G_k(matriceE_particuliere(), 2, 0.5, 0.8)
    aretes_LG = [];
    line_cover = [{'H','E','I'},{'L','J','I'},{'K','L'},{'F','J'},{'F','E','A','D'},\
                  {'C','A','B'},{'C','G'},{'D','H','G'}];
    for cliq in line_cover:
        aretes_LG.extend(list(it.combinations(cliq,2)));
    print("LG={}".format(aretes_LG));
    print("G_k={} , LG = {}".format(len(aretes_G_k), len(aretes_LG)))
    print("aretes_diff G_k et LG = {}".format(comparaison_aretes_G_k_LG(aretes_G_k, aretes_LG)));
def comparaison_aretes_G_k_LG_Debug():
    # A EFFACER
    aretes_G_k = fct_aux.liste_arcs(matriceE_particuliere());
    matE_G_2,aretes_G_k,dico = G_k(matriceE_particuliere(), 2, 0.5, 0.8)
    print("G_k={}".format(len(aretes_G_k)))
#    return aretes_G_k;
    aretes_LG = [];
    line_cover = [{'H','E','I'},{'L','J','I'},{'K','L'},{'F','J'},{'F','E','A','D'},\
                  {'C','A','B'},{'C','G'},{'D','H','G'}];
    line_cover = [{'F','J'},{'18','17','14'},{'H','L','K'},{'C','K'},{'E','H','I'},\
                  {'I','B','J','L'},{'15','B'},{'E','F','A'},{'A','D'},{'18','16','13'},\
                  {'15','16','17'},{'C','13','G'},{'G','D'}]
    for cliq in line_cover:
        aretes_LG.extend(list(it.combinations(cliq,2)));
    print("LG={}".format(aretes_LG));
    print("G_k={} , LG = {}".format(len(aretes_G_k), len(aretes_LG)))
    print("aretes_diff G_k et LG = {}".format(comparaison_aretes_G_k_LG(aretes_G_k, aretes_LG)));
    return aretes_G_k;
# ----- comparaison_aretes_G_k_LG ------
    
def simulation_iourte(matE_G_k, dico_proba_cases, args):
    """
    corriger le graphe G_k selon les paramatres args
    """
    # fichier de tracking or debug
    headers_df = ["G_cpt", "nbre_aretes_G_k", "nbre_aretes_LG", \
                  "dist_line", "nbre_aretes_diff_G_k_LG", \
                  "aretes_diff_G_k_LG", "C", "len(C)", "som_cout_min",\
                  "noeuds_corriges", "ordre_noeuds_traites",\
                  "min_DL", "mean_DL", "max_DL", "ecart_type",\
                  "max_cout", "max_permutation",\
                  "dico_som_min_permutations"];
            
    df_debug = pd.DataFrame( columns = headers_df);
    G_cpt = "G_"+str(args["k"]);
    
    # creation repertoire contenant distribution
    path_distr = Path(args["path_save"]+args["mode_select_noeuds_1"]+"_iourte/");
    path_distr.mkdir(parents=True, exist_ok=True);
    path_save = args["path_save"]+args["mode_select_noeuds_1"]+"_iourte/";
                                      
    # initialisation variables
    aretes_matE_G_k = fct_aux.liste_arcs(matE_G_k);
    
    try:
        ordre_noeuds_traites = [] # car liste_cliques = []
        cliques = []; # car graphe iourte;
        dico_cliq = dict();
        for noeud in matE_G_k.columns:
            dico_cliq[noeud] = -1
        dico_gamma_noeud = fct_aux.gamma_noeud(matE_G_k, aretes_matE_G_k) 
            
        # algo de correction selon methodes (arg_params["mode_select_noeuds_1"])
        dico_permutations = dict();
        dico_permutations = decouvClique.solution_methode_nodes_1(dico_gamma_noeud,\
                             cliques, aretes_matE_G_k, ordre_noeuds_traites, \
                             dico_cliq, dico_proba_cases, args);
        
        # Debut selection de la permutation de noeuds dont la distance hamming est la plus petite
        dico_sol = dict()
        dico_sol = best_permutation_iourte(dico_permutations, matE_G_k)
        # FIN selection de la permutation de noeuds dont la distance hamming est la plus petite
        
        
        # comparaison entre aretes_matE_G_k et aretes_LG
        cpt_aretes_G_k_notIn_LG = 0; #(en pourcentage)
        cpt_aretes_G_k_notIn_LG = comparaison_aretes_G_k_LG(aretes_matE_G_k, dico_sol["aretes_LG"]);

        
        # ecrire dans un fichier pouvant etre lu pendant qu'il continue d'etre ecrit
        f = open(path_save+"distribution_moyDistLine_G_k.txt","a")
        f.write(G_cpt+";"+str(args["k"])+";"+str(dico_sol["dist_line"])+";"\
                +str(len(aretes_matE_G_k))+";"+str(cpt_aretes_G_k_notIn_LG)+"\n")
        f.close();
        
        # pour debug, log, .....
        dico_som_min_permutations = dict();
        for l_noeuds_1, values in dico_permutations.items():
            if values[6] not in dico_som_min_permutations.keys():
                dico_som_min_permutations[values[6]] = [l_noeuds_1]
            else:
                dico_som_min_permutations[values[6]].append(l_noeuds_1)
                                 
        df_debug.loc[len(df_debug.index)] = [\
                        G_cpt, len(aretes_matE_G_k),\
                        dico_sol["nbre_aretes_LG"], dico_sol["dist_line"],\
                        dico_sol["nbre_aretes_diff_matE_G_k_LG"],\
                        dico_sol["aretes_diff_matE_G_k_LG"],\
                        dico_sol["C"], len(dico_sol["C"]), dico_sol["som_cout_min"],\
                        dico_sol["noeuds_corriges"], dico_sol["ordre_noeuds_traites"],\
                        dico_sol["min_line"], dico_sol["mean_line"],\
                        dico_sol["max_line"], dico_sol["ecart_type"],\
                        dico_sol["max_cout"], dico_sol["max_permutation"],\
                        dico_som_min_permutations]
#        CPT_DF_DEBUG += 1;              
        if args["k"] % 100 == 0:
            simu50.save_df(df_debug, path_save, args["k_deep"], headers_df)
            df_debug = pd.DataFrame( columns = headers_df)
            print("save {} fois".format( args["k"] ))      
    except Exception as e:
        print("####### EmptyDataError ", G_cpt, ": e = ", e," ####### ");
        df_debug.loc[len(df_debug.index)] = [G_cpt, len(aretes_matE_G_k), \
                        "error", "error", "error", "error", "error" ,\
                        "error", "error", "error", "error", "error",\
                        "error", "error", "error", "error","error", \
                        "error"];
                        
    simu50.save_df(df_debug, path_save, args["k_deep"], headers_df)
    pass    

#### ---- plot dist line ----
def equation_3k_6(row):
    return 3*row["k"]+6;
    
def plot_dist_line_vs_line_3k_6(args):
    """
    plot la courbe des dist-line en fonction de k et la droite 3k+6
    """
    file = args["path_save"]+args["mode_select_noeuds_1"]+"_iourte/distribution_moyDistLine_G_k.txt";
    df = pd.read_csv(file, names=["G_k","k","dl","nb_aretes","aretes_G_k_notIn_LG"],sep = ";");
    fig =  plt.figure(); 
    default_size = fig.get_size_inches()
    print("w =", default_size[0], " h = ",default_size[1])
    fig.set_size_inches( (default_size[0]*1.5, default_size[1]*1.5) )

    ax1 = fig.add_subplot(1,1,1);
    
    styles1 = ['bs-','r*-','--H','ro-','y^-','rs-','go-','b^-','bo-','-gD','-yp',\
               ':>','-.<','-v','-d','-h','--H','--,']    

    df["droite_3k+6"] = df.apply(lambda row: equation_3k_6(row), axis=1)
    df.plot(x="k",y=["dl","droite_3k+6","aretes_G_k_notIn_LG"], style= styles1, ax = ax1);
    ax1.set(xlabel= "G_k", ylabel= " distance line ")
    ax1.grid(True, which='both')
    plt.xticks( range(0, df["k"].max(),2) )
    plt.yticks( range(0, df["dl"].max(),6))
    path_save = args["path_save"]+args["mode_select_noeuds_1"]+"_iourte/";
    plt.savefig(path_save+"comparaison_distance_line_vs_3k_6_graphe_iourte_prior_"+\
                str(args["type_fct_cout"].split("_")[3])+"_"+\
                str(args["facteur_mult"])+".jpeg",dpi= 190);
    plt.clf();
    pass
#### ---- plot dist line ----

#### --- comparaison dl entre priorisation ajout/suppr pour w_i = [2,5,10,20,40,50]  ---
def compare_dl_priorisation(args):
    """
    creer un dataframe regroupant toute les types de priorisation suivie de leur poids w_i
    """
    df = pd.DataFrame();cols = [];
    for type_fct_cout in args["type_fct_couts"]:
        for facteur_multiplicatif in args["facteur_multiplicatifs"]:
            priorite = type_fct_cout.split("_")[3]; 
            file = args["base_rep"]+ str(priorite) +"_"+ \
                    str(facteur_multiplicatif)+"/aleatoire_iourte/distribution_moyDistLine_G_k.txt";
            dl_name = "dl_"+str(priorite)+"_"+str(facteur_multiplicatif);
            aretes_G_k_notIn_LG = "aretes_G_k_notIn_LG_"+str(priorite)+"_"+str(facteur_multiplicatif);
            df_tmp = pd.read_csv(file, names=["G_k","k",dl_name,"nb_aretes",aretes_G_k_notIn_LG],\
                                              sep = ";");
            df = pd.concat( [df, df_tmp], axis = 1);
            cols.append("dl_"+str(priorite)+"_"+str(facteur_multiplicatif));
            print("priorite = {}, facteur_multiplicatif={} termine".format(priorite,facteur_multiplicatif));
    df = df[cols];
    fig =  plt.figure(); default_size = fig.get_size_inches();
    print("w =", default_size[0], " h = ",default_size[1])
    fig.set_size_inches( (default_size[0]*1.5, default_size[1]*1.5) )

    ax1 = fig.add_subplot(1,1,1);
    styles1 = ['bs-','r*-','--H','ro-','y^-','rs-','go-','b^-','bo-','-gD','-yp',\
               ':>','-.<','-v','-d','-h','--H','--,'];
    # droite 3k+6     
    df["k"] = range(0,args["k_deep"]+1,1);     
    df["droite_3k+6"] = df.apply(lambda row: equation_3k_6(row), axis=1);
    df.plot(x="k",y=["droite_3k+6"], ax = ax1);
    # priorite suppression et ajout
    df.plot(x="k", y=cols, style= styles1, ax=ax1);
    
    # abscisses/ ordonnes
    ax1.set(xlabel= "G_k", ylabel= "distance line")
    ax1.grid(True, which='both')
    plt.xticks( range(0, df["k"].max(),2) )
    max_dl = max(df[cols].max());
    plt.yticks( range(0,max_dl,15));
    
    f, (ax2, ax3) = plt.subplots(1,2, figsize=(default_size[0]*2.5, default_size[1]*2.5));

    #### suppression 
    suppr_cols = [col for col in df.columns if "dl_supp" in col];
    df_supp = df[suppr_cols]; df_supp["k"] = range(0,args["k_deep"]+1,1);   
    df_supp.plot(x="k", y=suppr_cols, style= styles1, ax=ax2);
    df_supp["droite_3k+6"] = df_supp.apply(lambda row: equation_3k_6(row), axis=1);
    df_supp.plot(x="k",y=["droite_3k+6"], ax = ax2);
    ax2.set(xlabel= "G_k", ylabel= "distance line")
    ax2.grid(True, which='both')
    plt.xticks( range(0, df["k"].max(),2) )
    max_dl = max(df[cols].max());
    plt.yticks( range(0,max_dl,15)); 
    
    #### ajout 
    ajout_cols = [col for col in df.columns if "dl_ajout" in col];
    df_ajout = df[ajout_cols]; df_ajout["k"] = range(0,args["k_deep"]+1,1);   
    df_ajout.plot(x="k", y=ajout_cols, style= styles1, ax=ax3);
    df_ajout["droite_3k+6"] = df_ajout.apply(lambda row: equation_3k_6(row), axis=1);
    df_ajout.plot(x="k",y=["droite_3k+6"], ax = ax3);
    ax3.set(xlabel= "G_k", ylabel= "distance line")
    ax3.grid(True, which='both');
    plt.xticks( range(0, df["k"].max(),2) )
    max_dl = max(df[cols].max()); print("max_dl ={} ".format(max_dl))
    plt.yticks( range(0, max_dl,15));
    #### labels et etiquettes communes
    f.subplots_adjust(hspace=0.5);
    
    d = pd.concat([df_supp, df_ajout], axis=1);
    fig =  plt.figure(); default_size = fig.get_size_inches();
    fig.set_size_inches( (default_size[0]*1.5, default_size[1]*1.5) )
    ax1 = fig.add_subplot(1,1,1);
    d[["dl_ajout_2","dl_supp_2"]].plot(ax=ax1);
    plt.savefig(path_save+"comparaison_dl_ajout_supp_wi2_iourte_k_0_"+str(args["k_deep"])+".jpeg",dpi= 190);
    plt.clf();
    ax1 = fig.add_subplot(1,1,1);
    d[["dl_ajout_5","dl_supp_5"]].plot(ax=ax1);
    plt.savefig(path_save+"comparaison_dl_ajout_supp_wi5_iourte_k_0_"+str(args["k_deep"])+".jpeg",dpi= 190);
    plt.clf();
    ax1 = fig.add_subplot(1,1,1);
    d[["dl_ajout_10","dl_supp_10"]].plot(ax=ax1);
    plt.savefig(path_save+"comparaison_dl_ajout_supp_wi10_iourte_k_0_"+str(args["k_deep"])+".jpeg",dpi= 190);

    return df_supp, df_ajout;
    # save
    plt.savefig(args["path_save"]+"comparaison_priorite_poids_distance_line_graphe_iourte_k_0_"+args["k_deep"]+\
                ".jpeg",dpi= 190);
    pass
#### --- comparaison dl entre priorisation ajout/suppr pour k = [2,5,10,20,40,50]  ---

#### comparer dl prior 1,10,1 et dl supp 1.1.10
def comparer_dl_prior_suppr_1_10_1(args):
    df = pd.DataFrame();cols = []
    for type_fct_cout in args["type_fct_couts"]:
        for facteur_multiplicatif in args["facteur_multiplicatifs"]:
            priorite = type_fct_cout.split("_")[3]; 
            file = args["base_rep"]+ str(priorite) +"_"+ \
                    str(facteur_multiplicatif)+"/aleatoire_iourte/distribution_moyDistLine_G_k.txt";
            dl_name = "dl_"+str(priorite)+"_"+str(facteur_multiplicatif);
            aretes_G_k_notIn_LG = "aretes_G_k_notIn_LG_"+str(priorite)+"_"+str(facteur_multiplicatif);
            df_tmp = pd.read_csv(file, names=["G_k","k",dl_name,"nb_aretes",aretes_G_k_notIn_LG],\
                                              sep = ";");
            df = pd.concat( [df, df_tmp], axis = 1);
            cols.append("dl_"+str(priorite)+"_"+str(facteur_multiplicatif));
            print("priorite = {}, facteur_multiplicatif={} termine".format(priorite,facteur_multiplicatif));
    styles1 = ['bs-','r*-','--H','ro-','y^-','rs-','go-','b^-','bo-','-gD','-yp',\
               ':>','-.<','-v','-d','-h','--H','--,'];
    w_i = [1,10]; prior = ["ajout","supp"];
    cols = ["dl_ajout_1","dl_supp_1","dl_ajout_10","dl_supp_10"];
    for tu in it.combinations(it.product(w_i,prior),2):
        prior_0 = tu[0][1]; prior_1 = tu[1][1];
        w_i_0 = tu[0][0]; w_i_1 = tu[1][0];
#        if prior_0 != prior_1 :
#            if w_i_0 == 10 and w_i_1 == 10 :
#                pass
#            else:
        df_tu = df[ ["dl_"+prior_0+"_"+str(w_i_0),"dl_"+prior_1+"_"+str(w_i_1)] ];
        max_dl = max(df_tu.max());
        df_tu["k"] = range(0,args["k_deep"]+1,1); fig = plt.figure(figsize = (8,8)); 
        ax1 = fig.add_subplot(1,1,1);
        df_tu["droite_3k+6"] = df_tu.apply(lambda row: equation_3k_6(row), axis=1);
        df_tu.plot(x="k",y=["droite_3k+6"], ax = ax1);
        df_tu.plot(x="k", y=["dl_"+prior_0+"_"+str(w_i_0),"dl_"+prior_1+"_"+str(w_i_1)], style= styles1, ax=ax1);
        ax1.set(xlabel= "G_k", ylabel= "distance line")
        ax1.grid(True, which='both')
        plt.xticks( range(0, df_tu["k"].max(),2) ); 
        plt.yticks( range(0,max_dl,15)); 
        plt.savefig(args["path_save"]+"comparaison_prior_"+prior_0+"_wi_"+str(w_i_0)+\
                    "_"+prior_1+"_wi_"+ str(w_i_1) +\
                    "_distance_line_vs_3k_6_graphe_iourte.jpeg", dpi= 190);
        plt.clf();
    pass
def comparer_dl_prior_suppr_1_10_1_dbg(args):
    df = pd.DataFrame();cols = []
    fig =  plt.figure(); default_size = fig.get_size_inches();
    for type_fct_cout in args["type_fct_couts"]:
        for facteur_multiplicatif in args["facteur_multiplicatifs"]:
            priorite = type_fct_cout.split("_")[3]; 
            file = args["base_rep"]+ str(priorite) +"_"+ \
                    str(facteur_multiplicatif)+"/aleatoire_iourte/distribution_moyDistLine_G_k.txt";
            dl_name = "dl_"+str(priorite)+"_"+str(facteur_multiplicatif);
            aretes_G_k_notIn_LG = "aretes_G_k_notIn_LG_"+str(priorite)+"_"+str(facteur_multiplicatif);
            df_tmp = pd.read_csv(file, names=["G_k","k",dl_name,"nb_aretes",aretes_G_k_notIn_LG],\
                                              sep = ";");
            df = pd.concat( [df, df_tmp], axis = 1);
            cols.append("dl_"+str(priorite)+"_"+str(facteur_multiplicatif));
            print("priorite = {}, facteur_multiplicatif={} termine".format(priorite,facteur_multiplicatif));
    styles1 = ['bs-','r*-','--H','ro-','y^-','rs-','go-','b^-','bo-','-gD','-yp',\
               ':>','-.<','-v','-d','-h','--H','--,'];
    w_i = [1,10]; prior = ["ajout","supp"];
    cols = ["dl_ajout_1","dl_supp_1","dl_ajout_10","dl_supp_10"];
    for tu in it.combinations(it.product(w_i,prior),2):
        prior_0 = tu[0][1]; prior_1 = tu[1][1];
        w_i_0 = tu[0][0]; w_i_1 = tu[1][0];
#        if prior_0 != prior_1 :
#            if w_i_0 == 10 and w_i_1 == 10 :
#                pass
#            else:
        aretes_G_k_del_prior_0 = "aretes_G_k_notIn_LG_"+str(prior_0)+"_"+str(w_i_0);
        aretes_G_k_del_prior_1 = "aretes_G_k_notIn_LG_"+str(prior_1)+"_"+str(w_i_1);
        f, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(default_size[0]*2.0, default_size[1]*2.0));
        df_tu = df[ ["dl_"+prior_0+"_"+str(w_i_0),"dl_"+prior_1+"_"+str(w_i_1),\
                     aretes_G_k_del_prior_0, aretes_G_k_del_prior_1] ];
        max_dl = max(df_tu.max());
        df_tu["k"] = range(0,args["k_deep"]+1,1); fig = plt.figure(figsize = (8,8)); 
        df_tu["droite_3k+6"] = df_tu.apply(lambda row: equation_3k_6(row), axis=1);
        # ax1
        df_tu.plot(x="k",y=["droite_3k+6"], ax = ax1);
        df_tu.plot(x="k", y=["dl_"+prior_0+"_"+str(w_i_0),"dl_"+prior_1+"_"+str(w_i_1)], style= styles1, ax=ax1);
        ax1.set(xlabel= "G_k \n(a)", ylabel= "distance line", \
                xticks=range(0, df_tu["k"].max(),4),yticks=range(0,max_dl,15))
        ax1.grid(True, which='both')
        # ax2
        df_tu.plot(x="k",y=[aretes_G_k_del_prior_0, aretes_G_k_del_prior_1], style= styles1, ax=ax2);
        ax2.set(xlabel= "G_k \n(b)", ylabel= "pourcentage aretes supprimees", \
                xticks=range(0, df_tu["k"].max(),4), yticks=range(0,100,10))
        #ax3
        df_tu["ecart_"+prior_0+"_"+str(w_i_0)] = df_tu["dl_"+prior_0+"_"+str(w_i_0)] - df_tu["droite_3k+6"];
        df_tu["ecart_"+prior_1+"_"+str(w_i_1)] = df_tu["dl_"+prior_1+"_"+str(w_i_1)] - df_tu["droite_3k+6"];
        
        df_tu.plot(x="k",y=["ecart_"+prior_0+"_"+str(w_i_0), "ecart_"+prior_1+"_"+str(w_i_1)], style= styles1, ax=ax3);
        max_ecart = max(df_tu["ecart_"+prior_0+"_"+str(w_i_0)].max(), df_tu["ecart_"+prior_1+"_"+str(w_i_1)].max());
        ax3.set(xlabel= "G_k \n(c)", ylabel= "distance entre 3k+6 et DL", \
                xticks=range(0, df_tu["k"].max(),4),yticks=range(0,max_ecart,15))
              
        f.savefig(args["path_save"]+"comparaison_prior_"+prior_0+"_wi_"+str(w_i_0)+\
                    "_"+prior_1+"_wi_"+ str(w_i_1) +\
                    "_distance_line_vs_3k_6_graphe_iourte.jpeg", dpi= 190);
        plt.clf();
    pass

#### comparer dl prior 1,10,1 et dl supp 1.1.10


if __name__ == '__main__':
#    G_2 = comparaison_aretes_G_k_LG_Debug()
    start= time.time();
    k_deep = 50; 
#    ##### type_priorisations ######
#    p_correl = 0.7; correl_seuil = 0.8; # correl_seuil : valeur des correlations vrai positives
#    mode_select_noeuds_1 = "aleatoire"; number_permutations_nodes_1= 50; #100;
#    critere_selection_pi1_pi2 = 0; # 0: moins de modif,1: ajout aretes> supp aretes, 2:ajout aretes < supp aretes,
#    facteur_multiplicatif = 0; exposant = 0; # 0: fct_unitaire, 1:fct_normal, 2: fct_quadratique, 4:fct_quadruple, 5: fct_quintuple
#    type_fct_cout = "lineaire" #"cloche" # ou "lineaire"
#    
#    matE_G0 = matriceE_particuliere();
#    facteur_multiplicatifs = [1,2,5,10,20,40]; #facteur_multiplicatifs = [1,2,5,10];
#    type_fct_couts = ["lineare_iourte_priorite_supp","lineare_iourte_priorite_ajout"];
#    
#    for type_fct_cout in type_fct_couts:
#        for facteur_multiplicatif in facteur_multiplicatifs:
#            coef_fct_cout = (exposant, facteur_multiplicatif, type_fct_cout);
#            number_items_pi1_pi2 = 0.5#0.5;# 1;
#            path_save = "./graphe_particulier/graphe_particulier_prior_"+str(type_fct_cout.split("_")[3])+\
#                        "_"+str(facteur_multiplicatif)+"/"
#            args = {"k_deep": k_deep, "p_correl": p_correl, "correl_seuil": correl_seuil,\
#                    "number_items_pi1_pi2": number_items_pi1_pi2, \
#                    "number_permutations_nodes_1": number_permutations_nodes_1, \
#                    "mode_select_noeuds_1": mode_select_noeuds_1,\
#                    "coef_fct_cout": coef_fct_cout,\
#                    "critere_selection_pi1_pi2": critere_selection_pi1_pi2,\
#                    "facteur_mult": facteur_multiplicatif,\
#                    "type_fct_cout": type_fct_cout,\
#                    "path_save": path_save
#                    }
#            ## --- procedural ---
#            for k in range(0,k_deep+1):
#                matE_G_k, aretes_G_k, dico_proba_cases = \
#                G_k(matE_G0, k, args["p_correl"], args["correl_seuil"]);
#                
#                print("---> fini matE_G_k aretes ={}",format(matE_G_k.shape))
#                args["k"] = k;
#                simulation_iourte(matE_G_k, dico_proba_cases, args);
#            ## --- procedural ---
#            plot_dist_line_vs_line_3k_6(args) 
#    #####  type_priorisations ######
#    
#    ### courbe priorite poids ===> debut
    facteur_multiplicatifs = [1,2,5,10,20,40]; #facteur_multiplicatifs = [1,2,5,10];
    type_fct_couts = ["lineare_iourte_priorite_supp","lineare_iourte_priorite_ajout"];
    base_rep = "graphe_particulier/graphe_particulier_prior_"; path_save = "graphe_particulier/";
    args ={"type_fct_couts":type_fct_couts,"facteur_multiplicatifs":facteur_multiplicatifs,\
           "k_deep":k_deep,"base_rep":base_rep,"path_save":path_save};
#    compare_dl_priorisation(args);
#    ### courbe priorite poids ===> fin 

    ## priorite ajout[1,10,1]/suppr[1,1,10]
    args["facteur_multiplicatifs"] = [1,2,5,10,20,40]; #args["facteur_multiplicatifs"] = [1,2,5,10];
#    comparer_dl_prior_suppr_1_10_1(args);
    comparer_dl_prior_suppr_1_10_1_dbg(args)
    ## priorite ajout[1,10,1]/suppr[1,1,10]
    
    print (time.time() - start)
    
    
    
    
    
    
    
    
    
    
##    
###    p_correl = 0.7; correl_seuil = 0.8; # correl_seuil : valeur des correlations vrai positives
###    mode_select_noeuds_1 = "aleatoire";
###    critere_selection_pi1_pi2 = 0; # 0: moins de modif,1: ajout aretes> supp aretes, 2:ajout aretes < supp aretes,
###    number_permutations_nodes_1= 50; #100;
###    facteur_multiplicatif = 0; exposant = 0; # 0: fct_unitaire, 1:fct_normal, 2: fct_quadratique, 4:fct_quadruple, 5: fct_quintuple
###    type_fct_cout = "lineaire" #"cloche" # ou "lineaire"
###    """
###    priorite ajout aretes
###    """
####    exposant = 0; facteur_multiplicatif = 10; # 2; 5; 10
####    type_fct_cout = "lineare_iourte_priorite_ajout" #"cloche" # ou "lineaire"
###    """
###    priorite suppression aretes
###    """
###    exposant = 0; facteur_multiplicatif = 10; # 2; 5; 10 
###    type_fct_cout = "lineare_iourte_priorite_supp" #"cloche" # ou "lineaire"
###    
###    
###    coef_fct_cout = (exposant, facteur_multiplicatif, type_fct_cout)
###    number_items_pi1_pi2 = 0.5#0.5;# 1;
###    path_save = "./graphe_particulier_"+str(facteur_multiplicatif)+"/"
###    args = {"k_deep": k_deep, "p_correl": p_correl, "correl_seuil": correl_seuil,\
###            "number_items_pi1_pi2": number_items_pi1_pi2, \
###            "number_permutations_nodes_1": number_permutations_nodes_1, \
###            "mode_select_noeuds_1": mode_select_noeuds_1,\
###            "coef_fct_cout": coef_fct_cout,\
###            "critere_selection_pi1_pi2": critere_selection_pi1_pi2,\
###            "path_save": path_save
###            }
###           
###            
###    dbg = False;
###    if dbg:
###        k_deep = 10;
###    matE_G0 = matriceE_particuliere();
###    ## --- procedural ---
###    for k in range(0,k_deep+1):
###        matE_G_k, aretes_G_k, dico_proba_cases = \
###        G_k(matE_G0, k, args["p_correl"], args["correl_seuil"]);
###        
###        print("---> fini matE_G_k aretes ={}",format(matE_G_k.shape))
###        args["k"] = k;
###        simulation_iourte(matE_G_k, dico_proba_cases, args);
###    ## --- procedural ---
####
###    plot_dist_line_vs_line_3k_6(args) 
###    
###    
###    
###    ## --- parallelisation ---
###    modulo = 5; params = [];
###    for k in range(1,k_deep+1):
###        matE_G_k, aretes_G_k, dico_proba_cases = \
###        G_k(matE_G0, k, args["p_correl"], args["correl_seuil"]);
###        params.append([matE_G_k, dico_proba_cases, args]);
###        if k%modulo == 0:
####            p = Pool(mp.cpu_count()-1) 
####            p.starmap(simulation_iourte, params)
####            p.terminate();
###            print("parallel k={} params={}".format(k,len(params)))
###            params = [];
###    ## --- parallelisation ---
##
##    print (time.time() - start)
#
#
##
