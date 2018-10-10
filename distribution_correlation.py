#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 13:12:47 2017

@author: willy
"""
import math;
import os, re;
import numpy as np;
import pandas as pd;
import fonctions_auxiliaires as fct_aux;
import scipy as spy;
from scipy.stats import linregress;
from pathlib import Path;
import time;
import logging;
import collections;
import simulation_50_graphes_PARALLELE as simu50;
import simulation_data_reelles as simu_reels;
import matplotlib.pyplot as plt;
import matplotlib.font_manager as font_manager
import seaborn as sns;
import json, ujson;
import itertools as it;
from functools import reduce;
from operator import mul;

#%%
def detect_seuil(path_aleatoire):
    files = os.listdir(path_aleatoire)    
    files = [fichier for fichier in files if  re.match("^data_correl_seuil_.*",fichier)]
    files_tmp = set()
    for file in files:
        if type(file.split("_")[3]) == float:
            files_tmp.add(file.split("_")[3])
        else:
            files_tmp.add(float(file.split("_")[3].split("metrique")[0]))
 
    return files_tmp#set(files);
    pass
def distrib_line_hamming_distance_seuil(args):
    """
    distribution distance line/Hamming en fonction du seuil
    run:
        d = distrib_line_hamming_distance_seuil(args)
    """
    if args["dbg"]:
        sub = True; 
        rep = "champlan_newDataset"; # or champlan or velizy
        root = "/home/willy/topologyLearning/datas/"+rep+"/";
        root = root+"sous_graphe/" if sub else root;
        chemin_aleatoire = root+"aleatoire/"; metrique = "metrique_wil"
        chemin_matE_reel = root+"matrices/matE_reelles.csv";
        chemin_equip = root+"matrices_equipements/" 
        args={"path_aleatoire":chemin_aleatoire,"metrique":metrique,\
              "path_matE_reel":chemin_matE_reel,"path_root":chemin_equip}
    #line distance
    namefile = "distribution_DistLine_seuil.txt";
    
    path = args["path_aleatoire"]; #path_aleatoire = "/home/willy/topologyLearning/datas/champlan_newDataset/sous_graphe/aleatoire"
    seuils = detect_seuil(args["path_aleatoire"]); dico = dict();
    print("seuils={}".format(seuils))
    for seuil in seuils:
        file = path+"/data_correl_seuil_"+str(seuil)+"_metrique_wil/distribution/"+namefile
        df = pd.read_csv(file, names=["seuil","dl","som_cout_min","nbre_aretes_diff_matE_LG","nbre_aretes_matE"], sep=';')
        dico[seuil] = df.loc[0,"dl"]
    df_T_line = pd.DataFrame(dico, columns = seuils, index=[0]).T
    df_T_line.sort_index(inplace=True)
    
    # hamming distance
    matEReel = pd.read_csv(args["path_matE_reel"], index_col = "Unnamed: 0");
    edges_reel = fct_aux.liste_arcs(matEReel);
    dico = json.load(open(args["path_root"]+"seuil_aretes_LG_"+args["metrique"]+".json","r"));
    seuils = list(); dico_res = dict();
    for seuil, aretes_LG_predit in dico.items():
        seuils.append(seuil); 
        nbre_edges_diff=np.nan; edges_diff = [];
        aretes_LG_predit = [(arete[0],arete[1]) for arete in aretes_LG_predit];
        nbre_edges_diff, edges_diff = simu50.distance_hamming(edges_reel, aretes_LG_predit)
        dico_res[seuil] = nbre_edges_diff;
    df_T_ham = pd.DataFrame(dico_res, columns = seuils, index=[0]).T
    df_T_ham.sort_index(inplace=True)
    
    width = 8; height = 5.5; print("w =", width, " h = ",height)
    fig = plt.figure(figsize=(width*1.1,height*1.1));
    title = "distrib_distance_line_vs_seuil"; 
    ax1 = fig.add_subplot(1,2,1);
    df_T_line.plot(kind='bar',ax = ax1, title = title)
    ax1.set(ylabel="distance line", xlabel="seuils", title=title)
    
    title = "distrib_distance_hamming_vs_seuil"; 
    ax2 = fig.add_subplot(1,2,2);
    df_T_ham.plot(kind='bar',ax = ax2, title = title)
    ax2.set(ylabel="distance hamming", xlabel="seuils", title=title)
    
    fig.savefig(args["path_aleatoire"]+"distrib_line_hamming_distance_vs_seuil_"+args["metrique"]+".jpeg",dpi=300) 
    plt.clf();
    pass
###### distribution distance line/Hamming  selon seuil ===> Fin

###### distribution distance line/Hamming  selon seuil avant algo correction ===> debut
#%%
import test_correlation_shapelet_znorm as corr_shapelet;
def distrib_hamming_distance_seuil_avant_correction(args):
    """
    distribution des vrai positives et negatives avant execution des algos de correction
    
    run:
        args={"chemin_matrices":chemin_matrices,"chemin_datasets":chemin_datasets,\
          "metrique":metrique, "path_aleatoire":chemin_aleatoire, "dbg":True}
        distrib_hamming_distance_seuil_avant_correction(args)
    """
    if args["dbg"]:
        grandeur="P"; len_dim="len2930"; 
        args_chemins ={"len_dim":len_dim,"grandeur":grandeur, "chemin_matrices":args["chemin_matrices"], \
                       "chemin_datasets":args["chemin_datasets"]}
        args["matE_reel"] = pd.read_csv(args_chemins["chemin_matrices"]+"matE_reelles.csv",index_col = "Unnamed: 0");
        args["matE"] = corr_shapelet.correlation(args_chemins,args["matE_reel"].columns.tolist())
        
    print("matE_calcul: termine ===> debut")
    seuils = ["{:1.1f}".format(x) for x in np.linspace(0.0,1.0,11)]#[round(x,1) for x in np.linspace(0,1,11)]
    dico_ham_0, dico_ham_1 = dict(), dict();
    for seuil in seuils:
        dico_ham_0[str(seuil)] = 0; dico_ham_1[str(seuil)] = 0; 
    for row, col in fct_aux.range_2d(args["matE"].columns):
        for seuil in seuils:
            if args["matE"].loc[row, col] <= float(seuil)+0.1 and \
                args["matE"].loc[row, col] > float(seuil) and \
                row in args["matE_reel"].columns and \
                col in args["matE_reel"].columns and \
                args["matE_reel"].loc[row, col] == 1:
                dico_ham_1[str(seuil)] += 1; # dico[str(seuil)].append((row,col))
            if args["matE"].loc[row, col] <= float(seuil)+0.1 and \
                args["matE"].loc[row, col] > float(seuil) and \
                row in args["matE_reel"].columns and \
                col in args["matE_reel"].columns and \
                args["matE_reel"].loc[row, col] == 0:
                dico_ham_0[str(seuil)] += 1;
    
    print("dico_ham_0:{},dico_ham_1={}".format(dico_ham_0,dico_ham_1))
    width = 8; height = 5.5; print("w =", width, " h = ",height)
    fig = plt.figure(figsize=(width*1.1,height*1.1));
    print("df_T_ham_1: ===> debut")
    df_T_ham_1 = pd.DataFrame(dico_ham_1, columns = seuils, index=[0]).T
    df_T_ham_1.sort_index(inplace=True)
    title = "distrib_distance_hamming_vs_seuil\n_pour_correl_vrai_positives_1"; 
    ax1 = fig.add_subplot(1,2,1);
    df_T_ham_1.plot(kind='bar',ax = ax1, title = title)
    ax1.set(ylabel="distance hamming", xlabel="seuils", title=title)
    
    print("df_T_ham_0: ===> debut")
    df_T_ham_0 = pd.DataFrame(dico_ham_0, columns = seuils, index=[0]).T
    df_T_ham_0.sort_index(inplace=True)
    title = "distrib_distance_hamming_vs_seuil\n_pour_correl_vrai_negatives_0";
    ax2 = fig.add_subplot(1,2,2);
    df_T_ham_0.plot(kind='bar',ax = ax2, title = title)
    ax2.set(ylabel="distance hamming", xlabel="seuils", title=title)
    
    fig.savefig(args["path_aleatoire"]+"distrib_hamming_distance_vs_seuil_pour_correls_vrai{positives,negatives}_"+args["metrique"]+".jpeg",dpi=300) 
    plt.clf();
    pass
###### distribution distance line/Hamming  selon seuil avant algo correction ===> fin

####### test debug distribution k correlation errors  ====> debut 
def distribution_selon_errors(dico_proba_cases, dico_deleted_add_edges, seuil, error_correl):
    """
    plot, selon error_correl, la distribution des valeurs
    """
    faux_pos = "ajouter"; faux_neg = "supprimer" # {ajout, suppression} aretes
    values_distrib = list()
    if error_correl == "faux_positives":
        # faux positives
        aretes_faux_pos = dico_deleted_add_edges[faux_pos]
        values_distrib = [v for k,v in dico_proba_cases.items() if k in aretes_faux_pos]
    elif error_correl == "faux_negatives":
        # faux negatives
        aretes_faux_neg = dico_deleted_add_edges[faux_neg]
        values_distrib = [v for k,v in dico_proba_cases.items() if k in aretes_faux_neg]
    elif error_correl == "vrai_negatives":
#        key_error_correl = dico_deleted_add_edges.values();
        key_error_correl = [x for sublist in dico_deleted_add_edges.values() for x in sublist]
        values_distrib = [v for k,v in dico_proba_cases.items() if k not in key_error_correl and v < seuil]
    elif error_correl == "vrai_positives":
#        key_error_correl = dico_deleted_add_edges.values();
        key_error_correl = [x for sublist in dico_deleted_add_edges.values() for x in sublist]
        values_distrib = [v for k,v in dico_proba_cases.items() if k not in key_error_correl and v >= seuil]
        pass
    else:
        print("error correl don't exist")
        return [];
        
    return values_distrib
    pass
def test_distribution(arg_params):
    """
    generer le graphe puis son line graph LG
    sur LG retirer  k erreurs de correlation 
    et attribuer valeurs de probas
    plot distributions les 
        * correlations vrai positives 
        * correlations vrai negatives 
        * correlations faux positives 
        * correlations faux negatives 
    """
#    error_correlations = ["vrai_positives", "vrai_negatives", "faux_positives", "faux_negatives",]
    
    nbre_lien = arg_params["nbre_lien"]; #5;
    nbre_ts = 10; effet_joule = 0; epsilon = 0.75
    test = "FINI"; ascendant_1 = True; simulation = True; seuil_U = 0;
    
    path_distr_chemin = str(arg_params["mode_select_noeuds_1"])+"/"+"data_p_"+\
                            str(arg_params["p_correl"])+"/distribution/"
    path_distr = Path(path_distr_chemin)
    path_distr.mkdir(parents=True, exist_ok=True);
    
    nbre_graphe_genere = 1; cpt_graphe_genere = 0;
    while cpt_graphe_genere < nbre_graphe_genere:
        cpt_graphe_genere += 1;
        G_cpt = "G_"+str(cpt_graphe_genere)+"_"+str( arg_params["k"]);
        
        path = Path(str(arg_params["mode_select_noeuds_1"])+"/"+\
                            "data_p_"+str(arg_params["p_correl"])+"/"+G_cpt+'/datasets/');
        path.mkdir(parents=True, exist_ok=True);
        path = Path(str(arg_params["mode_select_noeuds_1"])+"/"+\
                        "data_p_"+str(arg_params["p_correl"])+"/"+G_cpt+'/matrices/');
        path.mkdir(parents=True, exist_ok=True)
        chemin_datasets = str(arg_params["mode_select_noeuds_1"])+"/"+\
                              "data_p_"+str(arg_params["p_correl"])+"/"+G_cpt+"/datasets/";
        chemin_matrices = str(arg_params["mode_select_noeuds_1"])+"/"+\
                              "data_p_"+str(arg_params["p_correl"])+"/"+G_cpt+"/matrices/";
        
        # generer matE
        matE, df_matA, dico_dual_arc_sommet = simu50.matriceE(arg_params["dimMatA"], \
                                                              nbre_lien, chemin_matrices,\
                                                              chemin_datasets, nbre_ts, epsilon, \
                                                              effet_joule, test)
        
        matE_k_alpha = None; deleted_edges = [];
        dico_proba_cases = dict();
        matE_k_alpha, dico_deleted_add_edges = simu50.modif_k_cases(matE.copy(), arg_params["k"], \
                                                             arg_params["methode_delete_add_edges"], 
                                                             arg_params["correl_seuil"])
        dico_proba_cases = simu50.ajouter_proba_matE(matE_k_alpha, dico_deleted_add_edges, \
                                              arg_params["loi_stats"],arg_params["p_correl"],\
                                              arg_params["correl_seuil"])
        print("dico_deleted_add_edges = ",dico_deleted_add_edges)
        ## parametres figure
        fig = plt.figure(1);
        fig.set_figheight(8); fig.set_figwidth(8) # width = largueur, height = longueur
        cpt_ax1 = 0; num_bins = 200;
        
        error_correlations = ["vrai_positives", "vrai_negatives", "faux_positives", "faux_negatives"]
        for error_correl in error_correlations:
            values_distrib = list()
            values_distrib = distribution_selon_errors(dico_proba_cases, \
                                                       dico_deleted_add_edges, \
                                                       arg_params["correl_seuil"], error_correl)
            if len(values_distrib) != 0:
                cpt_ax1 += 1;#cpt = num; # cpt += 1
                ax1 = fig.add_subplot(2,2,cpt_ax1);  # a chercher 
                bins = np.linspace(0,1,num_bins)
                sns.distplot(values_distrib,ax = ax1, bins = bins, kde = False)
                ax1.set(title = error_correl)
            pass
        pass
    pass
####### test debug distribution k correlation errors  ====> fin

####### distribution valeurs reelles #####=====> debut
def convert_adjlist_matrice(cliques, columns_matE_reel):
    """
    return une matrice d'adjacence a partir des cliques
    """
    matE_pred = pd.DataFrame(columns = columns_matE_reel, index = columns_matE_reel);
    for clique in cliques:
        for tuple in it.combinations(clique,2):
            if tuple[0] in columns_matE_reel and tuple[1] in columns_matE_reel:
                matE_pred.loc[tuple[0], tuple[1]] = 1;
                matE_pred.loc[tuple[1], tuple[0]] = 1;
            else:
                print("### ATTENTION: ", tuple[0]," ou ",tuple[1]," dont belong to columns_matE_reel")
    matE_pred.fillna(0.0, inplace = True)
    return matE_pred
    
def distrib_selon_errors(row, col, dico_distrib, error_correls, matEReel, matE_predit):
    """
    dire a quel error de correlation correspond l'arete (row,col)
    
    error_correls = ["vrai_positives", "vrai_negatives", "faux_positives", "faux_negatives"]
    """
    if matEReel.loc[row,col] == matE_predit.loc[row,col] and matEReel.loc[row,col] == 1:
        dico_distrib[error_correls[0]].append((row,col))
    elif matEReel.loc[row,col] == matE_predit.loc[row,col] and matEReel.loc[row,col] == 0:
        dico_distrib[error_correls[1]].append((row,col))
    elif matEReel.loc[row,col] != matE_predit.loc[row,col] and matEReel.loc[row,col] == 0 and matE_predit.loc[row,col] == 1:
        dico_distrib[error_correls[2]].append((row,col))
    elif matEReel.loc[row,col] != matE_predit.loc[row,col] and matEReel.loc[row,col] == 1 and matE_predit.loc[row,col] == 0:
        dico_distrib[error_correls[3]].append((row,col))
    else:
        print("----(",row,",",col,") dont belong to column matE_reelles----")
        print("matEReel[",row,",",col,"] = ",matEReel.loc[row,col], " != matE_predit[",row,",",col,"] = ",matE_predit.loc[row,col] )
    return dico_distrib
    
def autolabel(rects, fig):
    """
    IMPORTANT
    Attach a text label above each bar displaying its height
    """
    for rect in rects:
#        print("rect.get_x = ", rect.get_x(), " rect.get_width = ",rect.get_width() )
        height = rect.get_height()
        fig.text(rect.get_x() + rect.get_width()/2., 1.05*height,
                '%d' % int(height),ha='center', va='bottom')
def distribution_dataReelles(params):
    """
    plot, pour chaque seuil, les distributions des 
        * correlations vrai positives 
        * correlations vrai negatives 
        * correlations faux positives 
        * correlations faux negatives
    params = {"path_datasets": ,"path_matE_proba": ,"path_matE_reel": ,"path_cliques":}
    """
    matEReels_proba = pd.read_csv( params["path_matE_proba"], index_col = "Unnamed: 0")
    matEReel = pd.read_csv( params["path_matE_reel"], index_col = "Unnamed: 0")
    dico_json = json.load(open(params["path_cliques"],"r"));
    error_correls = ["vrai_positives", "vrai_negatives", "faux_positives", "faux_negatives"]
    
    for correl, cliques in dico_json.items():
#        matE_predit = matEReels_proba.copy();
#        matE_predit[matE_predit < float(correl)] = 0;
#        matE_predit[matE_predit >= float(correl)] = 1;
        ## parametres figure
        fig = plt.figure(1);
        ax = fig.add_subplot(111);
        fig.set_figheight(8); fig.set_figwidth(8) # width = largueur, height = longueur
        num_bins = 200;

        dico_distrib = {error_correls[0]:[], error_correls[1]:[], error_correls[2]:[], error_correls[3]:[]}
        matE_predit = convert_adjlist_matrice(cliques, matEReels_proba.columns.tolist())
        for row, col in fct_aux.range_2d(matEReels_proba.columns.tolist()):
            dico_distrib = distrib_selon_errors(row, col, dico_distrib, error_correls, matEReel, matE_predit)
        dico_distr_num_errors = {x: len(y) for x,y in dico_distrib.items()}
        
        rect = ax.bar(range(len(dico_distr_num_errors)),dico_distr_num_errors.values(), align='center'); 
        plt.xticks(range(len(dico_distr_num_errors)), dico_distr_num_errors.keys());
        autolabel(rect, ax)
#        for i,v in enumerate(dico_distr_num_errors.values()):
#            print("v=",v)
#            plt.text(v + 3, i + .25, str(v), color='blue', fontweight='bold')
        
        fig.savefig(params["path_save"]+"distribution_dataReels_errorCorrelations_"+\
                    str(correl)+"_"+params["metrique_distance"]+".jpeg")
        plt.clf()
    pass
####### distribution valeurs reelles #####=====> fin

##### comprendre correlation trouvee par l'algorithme. ==> debut
def tracer_courbes_equipementCorreles(nodes, params, derivates):
    # param plot
    font_prop = font_manager.FontProperties( size=11)
#    default_size = fig.get_size_inches()
#    print("w =", default_size[0], " h = ",default_size[1])
#    fig.set_size_inches( (default_size[0]*2.5, default_size[1]*2.5) )
    width = 8; height = 5.5;
    print("w =", width, " h = ",height)
    fig = plt.figure( figsize=(width*2.5, height*2.5) )
    # grandeurs
    grandeurs = fct_aux.liste_grandeurs(params["chemin_datasets"]);
    grandeurs = ["P"]
    cpt_grand = 0;
    for grandeur in grandeurs:
        df = pd.read_csv(chemin_datasets+"dataset_"+grandeur+".csv");
        df.fillna(0, inplace=True)
        if derivates:
            cols = df.columns.tolist()
            arr = df.values
            df = pd.DataFrame(np.gradient(arr,10)[0], columns= cols)
        cpt_grand += 1;
        ax1 = fig.add_subplot(3,5,cpt_grand)
        df[nodes].plot(ax = ax1)
        ax1.set(title = str(grandeur))
        ax1.legend( nodes, loc='upper center', bbox_to_anchor=(0.5, 1.00),\
                   ncol=2, fancybox=True, shadow=True, prop = font_prop);
        pass
    fig.tight_layout()
    fig.savefig(params["path_save"]+"comprendre_correl_by_mesures.jpeg",dpi=600)
    pass     
##### comprendre correlation trouvee par l'algorithme. ==> FIN

#### distribution valeurs {fausses, vrai} {negatives, positives} correlations selon methode correlation ####
def distri_selon_errors_valeurs(row, col, dico_correl, type_errors, matE_seuil, matE_reel):
    """
    return un dico dont la cle est le type error de correlation et 
                        la valeur est le nombre d'aretes.
    """
    if row in matE_seuil.columns.tolist() and col in matE_seuil.columns.tolist() and \
        row in matE_reel.columns.tolist() and col in matE_reel.columns.tolist() :
        if matE_reel.loc[row,col] == matE_seuil.loc[row,col] and matE_reel.loc[row,col] == 1:
            dico_correl[type_errors[0]] += 1;
        elif matE_reel.loc[row,col] == matE_seuil.loc[row,col] and matE_reel.loc[row,col] == 0:
            dico_correl[type_errors[1]] += 1;
        elif matE_reel.loc[row,col] != matE_seuil.loc[row,col] and matE_reel.loc[row,col] == 0 and matE_seuil.loc[row,col] == 1:
            dico_correl[type_errors[2]] += 1;
        elif matE_reel.loc[row,col] != matE_seuil.loc[row,col] and matE_reel.loc[row,col] == 1 and matE_seuil.loc[row,col] == 0:
            dico_correl[type_errors[3]] += 1;
        else:
            print("----(",row,",",col,") dont belong to column matE_reelles----")
            print("matE_reel[",row,",",col,"] = ",matE_reel.loc[row,col], " != matE_seuil[",row,",",col,"] = ",matE_seuil.loc[row,col] )
    else:
        print("ERROR col=",col," and row=",row," dont in matE_reel")
    return dico_correl;
    
def distri_selon_errors_aretes(row, col, dico_correl, type_errors, matE_seuil, matE_reel):
    """
    return un dico dont la cle est le type error de correlation et 
                        la valeur est la liste d'aretes.
    """
    if row in matE_seuil.columns.tolist() and col in matE_seuil.columns.tolist() and \
        row in matE_reel.columns.tolist() and col in matE_reel.columns.tolist() :
        if matE_reel.loc[row,col] == matE_seuil.loc[row,col] and matE_reel.loc[row,col] == 1:
            dico_correl[type_errors[0]].append( [row,col] );
        elif matE_reel.loc[row,col] == matE_seuil.loc[row,col] and matE_reel.loc[row,col] == 0:
            dico_correl[type_errors[1]].append( [row,col] );
        elif matE_reel.loc[row,col] != matE_seuil.loc[row,col] and matE_reel.loc[row,col] == 0 and matE_seuil.loc[row,col] == 1:
            dico_correl[type_errors[2]].append( [row,col] );
        elif matE_reel.loc[row,col] != matE_seuil.loc[row,col] and matE_reel.loc[row,col] == 1 and matE_seuil.loc[row,col] == 0:
            dico_correl[type_errors[3]].append( [row,col] );
        else:
            print("----(",row,",",col,") dont belong to column matE_reelles----")
            print("matE_reel[",row,",",col,"] = ",matE_reel.loc[row,col], " != matE_seuil[",row,",",col,"] = ",matE_seuil.loc[row,col] )
    else:
        print("ERROR col=",col," and row=",row," dont in matE_reel")
    return dico_correl;

def distribution_valeurs_correlation_selon_seuil_methode(args):
    """
    obtenir la distribution des valeurs de correlation selon la methode et le seuil
    
    args={"matE_proba":, "matE_reel":,"methode_correl":,}
    """
    if args["dbg_seuil"] == True:
        print("OLLALAL")
        file = "champlan";
        rep_root = "/home/willy/topologyLearning/datas/"+file+"/";

        sub_graph = "sous_graphe/" #"sous_graphe/" # "sous_graphe/" or ""
        metrique_distance = "pearson"
        chemin_datasets = rep_root+sub_graph+"datasets/";
        chemin_matrices = rep_root+sub_graph+"matrices/";
        chemin_equipements = rep_root+sub_graph+"matrices_equipements/";
        matE_proba = pd.read_csv(chemin_equipements+"matE.csv",index_col = "Unnamed: 0");
        matE_reel = pd.read_csv(chemin_equipements+"matE_reelles.csv",index_col = "Unnamed: 0");
        methode_correl = "correlationParMorceaux" # ou liste_mode = ["correlationParMorceaux","correlationGlissante"]
        path_save = chemin_equipements;
        args ={"matE_proba": matE_proba, "matE_reel": matE_reel, \
               "methode_correl":methode_correl, "path_save": path_save, \
               "metrique_distance": metrique_distance, "dbg": True};
        
    
    matE_proba = args["matE_proba"]; matE_reel = args["matE_reel"]
    if set(matE_proba.columns.tolist()) != set(matE_reel.columns.tolist()):
        print("ERROR matE_proba et matE_reel have different set columns")
        if set(matE_proba.columns.tolist()).issubset( set(matE_reel.columns.tolist())):
            print("===>ICI<===")
#        return 
    seuil_correls = np.linspace(0,1,11);
    dico = dict() # {s:dico_correl} dico_correl = {"VraiPositive":,"VraiNegative":,"FauxPosit"}
    type_errors = list()
    for seuil_correl in seuil_correls:
        matE_seuil = simu_reels.matrice_binaire(matE_proba.copy(), seuil_correl)
        dico_correl = {"VraiPositive":0,"VraiNegative":0,"FauxPositive":0,"FauxNegative":0};
        type_errors = ["VraiPositive","VraiNegative","FauxPositive","FauxNegative"] #list(dico_correl.keys())
        for row, col in fct_aux.range_2d(matE_reel.columns.tolist()):
#            print("ICI----row={},col={}".format(row,col))
            row, col = verifier_row_col(row, col, matE_seuil.columns.tolist(), matE_reel.columns.tolist())
            if row != None and col != None:
                dico_correl = distri_selon_errors_valeurs(row, col, dico_correl, type_errors, matE_seuil, matE_reel);
        dico[seuil_correl] = dico_correl;
        
    ## transform dico en dataframe
    df = pd.DataFrame(dico, columns = seuil_correls)
    df = df.T

    ### plot
    width = 8; height = 5.5;
    print("w =", width, " h = ",height)
    fig = plt.figure(figsize=(width*1.5,height*1.5))
    
    cpt_ax = 1;
    for type_error in type_errors:
        ax1 = fig.add_subplot(2,2,cpt_ax); cpt_ax += 1
        df[type_error].plot(ax = ax1, kind = "bar", title = type_error+"_error_correlation")
        ax1.set(ylabel= "nombre d'erreurs", xlabel= "valeurs de seuils")
        autolabel(ax1.patches, ax1)
    fig.tight_layout()
    fig.savefig(args["path_save"]+"error_par_methode_<"+args["methode_correl"]+","+args["metrique_distance"]+">_correlation_"+args["type_cal_matE"]+".jpeg",dpi=300);
    plt.clf();
    pass
#### distribution valeurs {fausses, vrai} {negatives, positives} correlations ####

#### distribution formalisation matE  ===> DEBUT
def distrib_formalisation_matE(args):
    """
    distribution des valeurs de correlations sur les vrai positives, negatives erreurs de correlations
    """
    matE_reel = pd.read_csv(args["chemin_matrices"]+"matE_reelles.csv", index_col = "Unnamed: 0")
#    for metrique_distance in args["metrique_distance"]:
        
    liste_Dico = ujson.load(open(args["chemin_matrices"]+"correlation_grandeurs_"+args["metrique_distance"]+".json","r"));
    cols = matE_reel.columns.tolist();
    dico_errors={ "vrai_positives":[], "vrai_negatives":[] }; 
    indice_correlations = 0;
    for dico in liste_Dico:
        arete = tuple(dico["key"]); correls = dico["value"];
        ## delete zero correlation
        correls = [x for x in correls if x != 0];
#        print("correls={}".format(correls))
        ## delete zero correlation
        if len(correls) != 0:
            if (arete[0] in cols and arete[1] in cols) and \
            (matE_reel.loc[arete[0], arete[1]] == 1  or matE_reel.loc[arete[1], arete[0]] == 1):
                max_correl = sorted(correls,reverse=True)[indice_correlations] # value max
                max_correl =  reduce(mul, correls, 1)# produit
                max_correl = np.mean(correls) # means
                dico_errors["vrai_positives"].append( max_correl )
            elif (arete[0] in cols and arete[1] in cols) and \
            (matE_reel.loc[arete[0], arete[1]] == 0  or matE_reel.loc[arete[1], arete[0]] == 0):
                max_correl = sorted(correls,reverse=True)[indice_correlations]  # value max
                max_correl =  reduce(mul, correls, 1)# produit
                max_correl = np.mean(correls) # means
                dico_errors["vrai_negatives"].append( max_correl )
            elif arete[0] not in cols and arete[1] not in cols:
                print("aretes {},{} not belong to matE_Reel".format(tuple(arete),arete[1]))
        pass   
    # plot
#        default_size = fig.get_size_inches()
#        print("w =", default_size[0], " h = ",default_size[1])
#        fig.set_size_inches( (default_size[0]*1.1, default_size[1]*1.1) );
    width = 8; height = 5.5; print("w =", width, " h = ",height)
    fig = plt.figure(figsize=(width*1.1,height*1.1)) #fig.set_size_inches( (width*1.1, height*1.1) )
    cpt_ax = 1;
    for type_error in ["vrai_positives", "vrai_negatives"]:
        ax1 = fig.add_subplot(1,2,cpt_ax); cpt_ax += 1
#        dico_errors.loc[type_error].plot(ax = ax1, kind = "bar", title = type_error+"_error_correlation")
        z = pd.DataFrame(dico_errors[type_error], columns = [type_error])
        z.fillna(0,inplace=True)
        print("z={}".format(z.values))
        bins = np.linspace(0,1,11)
        ax1.hist(z.values,bins=bins)
        ax1.set(ylabel= "nombre d'erreurs", xlabel= "valeurs de correlations", title=type_error+"_error_correlation")
        autolabel(ax1.patches, ax1)
    fig.tight_layout()
    fig.savefig(args["path_save"]+"distrib_formalisation_matE_vrai_{positives,negatives}_"+\
                args["metrique_distance"]+"_position_correl_"+str(indice_correlations)+"_"+args["type_cal_matE"]+".jpeg",dpi=300) 
    plt.clf();                      
    pass
#### distribution formalisation matE  ===> FIN

#### courbe des paires d'equipements dont leur correlation >= seuil_correl ===> debut
def plot_equipements_selon_correlations(args):
    """
    courbe des paires d'equipements dont leur correlation >= seuil_correl
    """
    if args["dbg"]:
        file = "champlan";
        rep_root = "/home/willy/topologyLearning/datas/"+file+"/";

        sub_graph = "sous_graphe/" #"sous_graphe/"  {"sous_graphe/", ""}
        metrique_distance = "metrique_wil"
        chemin_datasets = rep_root+sub_graph+"datasets/";
        chemin_matrices = rep_root+sub_graph+"matrices/";
        chemin_equipements = rep_root+sub_graph+"matrices_equipements/";
        grandeur = "P";
        matE_proba = pd.read_csv(chemin_matrices+"matrice_adjacence_proba_"+grandeur+".csv",index_col = "Unnamed: 0");
        matE_reel = pd.read_csv(chemin_equipements+"matE_reelles.csv",index_col = "Unnamed: 0");
        methode_correl = "correlationParMorceaux" # ou liste_mode = ["correlationParMorceaux","correlationGlissante"]
        path_save = chemin_equipements;
        correl_vraiNeg = 0.8; correl_vraiPos = 0.2; correl_fauxNeg = 0.7; correl_fauxPos = 0.6;
#        correl_vraiNeg = args["seuil_vrai_negatif"]; correl_vraiPos = args["seuil_vrai_positif"];
        args ={"matE_proba": matE_proba, "matE_reel": matE_reel, "grandeur":grandeur,\
               "methode_correl":methode_correl, "path_save": path_save, \
               "metrique_distance": metrique_distance, "dbg": True,\
               "chemin_datasets":chemin_datasets, "chemin_matrices":chemin_matrices,\
               "correl_vraiNeg":correl_vraiNeg, "correl_vraiPos":correl_vraiPos,\
               "correl_fauxNeg":correl_fauxNeg, "correl_fauxPos":correl_fauxPos};
    
    matE_proba = args["matE_proba"]; matE_reel = args["matE_reel"]
    if set(matE_proba.columns.tolist()) != set(matE_reel.columns.tolist()):
        print("ERROR matE_proba et matE_reel have different set columns")
        if set(matE_proba.columns.tolist()).issubset( set(matE_reel.columns.tolist())):
            print("===>ICI<===")
#        return 
    seuil_correls = np.linspace(0,1,11);
    dico = dict() # {seuil:dico_correl} dico_correl = {"VraiPositive":,"VraiNegative":,"FauxPosit"}
    type_errors = list()
    for seuil_correl in seuil_correls:
        seuil_correl = round(seuil_correl,1)
        matE_seuil = simu_reels.matrice_binaire(matE_proba.copy(), seuil_correl)
        dico_correl = {"VraiPositive":[],"VraiNegative":[],"FauxPositive":[],"FauxNegative":[]};
        type_errors = ["VraiPositive","VraiNegative","FauxPositive","FauxNegative"]#list(dico_correl.keys())
        for row, col in fct_aux.range_2d(matE_reel.columns.tolist()):
            dico_correl = distri_selon_errors_aretes(row, col, dico_correl, type_errors, matE_seuil, matE_reel);
        dico[seuil_correl] = dico_correl;
        
#    for correl, dico_correl in dico.items():
#        if round(correl,1) == args["correl_vraiNeg"]:
#            aretes = dico_correl["VraiNegative"];
#            print("vraiNegative = {}".format(len(aretes)))
#            plot_aretes(args["grandeur"], correl, aretes, args["chemin_datasets"], args["path_save"], "VraiNegative");
#        if round(correl,1) == args["correl_vraiPos"]:
#            aretes = dico_correl["VraiPositive"];
#            print("vraiPositive = {}".format(len(aretes)))
#            plot_aretes(args["grandeur"], correl, aretes[:10], args["chemin_datasets"], args["path_save"], "VraiPositive");  
#        if round(correl,1) == args["correl_fauxPos"]:
#            aretes = dico_correl["FauxPositive"];
#            print("fauxPositive = {}".format(len(aretes)))
#            plot_aretes(args["grandeur"], correl, aretes[:10], args["chemin_datasets"], args["path_save"], "FauxPositive");
#        if round(correl,1) == args["correl_fauxNeg"]:
#            aretes = dico_correl["FauxNegative"];
#            print("FauxNegative = {}".format(len(aretes)))
#            plot_aretes(args["grandeur"], correl, aretes[:10], args["chemin_datasets"], args["path_save"], "FauxNegative");
                        
    df_tmp = pd.DataFrame(dico, columns = [round(s,1)for s in seuil_correls])
    df_tmp = df_tmp.T;
    print("index = {}, columns = {}".format(df_tmp.index.tolist(), df_tmp.columns.tolist()));
    for error_correl in df_tmp.columns.tolist():
        if error_correl == "VraiPositive":
            aretes = df_tmp.loc[args["correl_vraiPos"], error_correl]
            print("vraiPositive = {}, correl={}".format(len(aretes), args["correl_vraiPos"]));
            plot_aretes(args["grandeur"],args["correl_vraiPos"], aretes[:10], \
                        args["chemin_datasets"], args["path_save"], "VraiPositive");  
        elif error_correl == "VraiNegative":
            aretes = df_tmp.loc[args["correl_vraiNeg"], error_correl]
            print("vraiNegative = {}, correl={}".format(len(aretes),args["correl_vraiNeg"]));
            plot_aretes(args["grandeur"],args["correl_vraiNeg"], aretes[:10], \
                        args["chemin_datasets"], args["path_save"], "VraiNegative");  
        elif error_correl == "FauxNegative":
            aretes = df_tmp.loc[args["correl_fauxNeg"], error_correl]
            print("fauxNegative = {}, correl={}".format(len(aretes),args["correl_fauxNeg"]));
            plot_aretes(args["grandeur"],args["correl_fauxNeg"], aretes[:10], \
                        args["chemin_datasets"], args["path_save"], "FauxNegative");
        elif error_correl == "FauxPositive":
            aretes = df_tmp.loc[args["correl_fauxPos"], error_correl]
            print("fauxPositive = {}, correl={}".format(len(aretes),args["correl_fauxPos"]));
            plot_aretes(args["grandeur"],args["correl_fauxPos"], aretes[:10], \
                        args["chemin_datasets"], args["path_save"], "FauxPositive");

    pass
def sous_plot(number_pair):
    """
    return le nombre de lignes w et columns h
    """
    if number_pair < 1:
        return None, None;
    elif number_pair >= 1 and number_pair < 4:
        return 1, number_pair;
    elif number_pair >= 4 and number_pair < 10:
        return 3, math.floor( 10/3);
    elif number_pair >= 10 and number_pair < 16:
        return 4, math.floor( 16/4);
    elif number_pair >= 16 and number_pair < 24:
        return 4, math.ceil( 24/4);
    else :
        return 5, math.ceil( number_pair/5);
def plot_aretes(grandeur, correl, pair_equips, chemin_datasets, path_save, error_correl):
    #fig.set_figheight(8); fig.set_figwidth(8) # width = largueur, height = longueur
#    default_size = fig.get_size_inches()
    width = 8; height = 5.5
    fig = plt.figure(figsize=(width*1.15,height*1.15))
    fig.set_size_inches( (width*1.15, height*1.15) )\
    if len(pair_equips) < 9 \
    else fig.set_size_inches( (width*2.75, height*2.75) ); #fig.set_size_inches( (width*2.15, height*2.15) );
    cpt_ax = 1;
    w,h = sous_plot(len(pair_equips));
    print("w={},h={}, len pair={}".format(w,h,len(pair_equips)))
    if w == None:
        print("correl={} have 0 pair_equips=aretes".format(correl))
        return;
    df = pd.read_csv(chemin_datasets+"dataset_"+grandeur+".csv");
    df = df.set_index("timestamp") if "timestamp" in df.columns.tolist() \
                                    else df.set_index("Unnamed: 0");
    df = (df - df.mean(skipna=True))/ df.std(skipna=True); # df = normalize(df.copy())
    df.fillna(0, inplace=True)
    for pair in pair_equips:
        ax = fig.add_subplot(w,h,cpt_ax); cpt_ax += 1;
        df[pair].plot(ax = ax)
        ax.set(title = pair)
    fig.tight_layout();
    fig.savefig(path_save+"plots_correl_"+str(correl)+"_error_<"+error_correl+">_grandeur_"+grandeur+".jpeg",dpi=70); 
    plt.clf();
    pass
#### courbe des paires d'equipements dont leur correlation >= seuil_correl ===> fin

#### distribution des cases de matEreelles  a 0 et a 1  ===> debut
#def distri_01(row, col, seuil, dico_correl, type_errors, matE_proba, matE_reel, excludes):
#    """ ["VraiPositive","VraiNegative"] """
###    if row in matE_proba.columns.tolist() and col in matE_proba.columns.tolist() and \
###        row in matE_reel.columns.tolist() and col in matE_reel.columns.tolist() :
##    if (row.split("->")[0]+"->"+row.split("->")[1] in matE_proba.columns.tolist() or row.split("->")[1]+"->"+row.split("->")[0] in matE_proba.columns.tolist()) and \
##       (col.split("->")[0]+"->"+col.split("->")[1] in matE_proba.columns.tolist() or col.split("->")[1]+"->"+col.split("->")[0] in matE_proba.columns.tolist()) and \
##       (row.split("->")[0]+"->"+row.split("->")[1] in matE_reel.columns.tolist() or row.split("->")[1]+"->"+row.split("->")[0] in matE_reel.columns.tolist()) and \
##       (col.split("->")[0]+"->"+col.split("->")[1] in matE_reel.columns.tolist() or col.split("->")[0]+"->"+col.split("->")[1] in matE_reel.columns.tolist()) :
###        print("seuil={} ==> mat[{},{}]={}".format(seuil,row,col,matE_proba.loc[row,col]))
##        if matE_reel.loc[row,col] == 1 and \
##            ((row,col) not in excludes or (col,row) not in excludes) and \
##            ((seuil-0.1 <= matE_proba.loc[row,col] and  matE_proba.loc[row,col] < seuil) or \
##             (seuil-0.1 <= matE_proba.loc[col,row] and  matE_proba.loc[col,row] < seuil)):
##            dico_correl[type_errors[0]].append( [row,col] );
##            excludes.append( (row,col) )
##        elif matE_reel.loc[row,col] == 0 and \
##            ((row,col) not in excludes or (col,row) not in excludes) and \
##            ((seuil-0.1 <= matE_proba.loc[row,col] and  matE_proba.loc[row,col] < seuil) or \
##             (seuil-0.1 <= matE_proba.loc[col,row] and  matE_proba.loc[col,row] < seuil) ):
##            dico_correl[type_errors[1]].append( [row,col] );
##            excludes.append( (row,col) );
##        else:
##            print("----(",row,",",col,") dont belong to column matE_reelles----")
##    else:
##        print("ERROR col=",col," and row=",row," dont in matE_reel")
##    return dico_correl, excludes;
##    pass
#
#    #### new version ==>
#    if row.split("->")[0]+"->"+row.split("->")[1] in matE_proba.columns and \
#        col.split("->")[0]+"->"+col.split("->")[1] in matE_proba.columns: 
#            arete, excludes, type_error = ajout_aretes_selon_type_error(row.split("->")[0]+"->"+row.split("->")[1], \
#                                        col.split("->")[0]+"->"+col.split("->")[1], \
#                                        seuil, excludes, matE_reel, matE_proba)
#            print("1 arete={} type={}, type_error={}".format(arete, type(arete), type_error))
#            if arete != -1 and type_error != None:
#                dico_correl[type_error].append(arete);
#    elif row.split("->")[0]+"->"+row.split("->")[1] in matE_proba.columns and \
#        col.split("->")[1]+"->"+col.split("->")[0] in matE_proba.columns:
#            arete, excludes, type_error = ajout_aretes_selon_type_error(row.split("->")[0]+"->"+row.split("->")[1], \
#                                        col.split("->")[1]+"->"+col.split("->")[0],\
#                                        seuil, excludes, matE_reel, matE_proba)
#            print("2 arete={}".format(arete))
#            if arete != -1 and type_error != None:
#                dico_correl[type_error].append(arete)
#    elif row.split("->")[1]+"->"+row.split("->")[1] in matE_proba.columns and \
#        col.split("->")[0]+"->"+col.split("->")[1] in matE_proba.columns:
#            arete, excludes, type_error = ajout_aretes_selon_type_error(row.split("->")[1]+"->"+row.split("->")[1], \
#                                        col.split("->")[0]+"->"+col.split("->")[1], \
#                                        seuil, excludes, matE_reel, matE_proba)
#            print("3 arete={}".format(arete))
#            if arete != -1 and type_error != None:
#                dico_correl[type_error].append(arete)
#    elif row.split("->")[1]+"->"+row.split("->")[0] in matE_proba.columns and \
#        col.split("->")[1]+"->"+col.split("->")[0] in matE_proba.columns:
#            arete, excludes, type_error = ajout_aretes_selon_type_error(row.split("->")[1]+"->"+row.split("->")[1], \
#                                        col.split("->")[1]+"->"+col.split("->")[0], \
#                                        seuil, excludes, matE_reel, matE_proba)
#            print("4 arete={}".format(arete))
#            if arete != -1 and type_error != None:
#                dico_correl[type_error].append(arete)
#    else:
#        print("ERROR col=",col," and row=",row," dont in matE_reel")       
#    return dico_correl, excludes;
#    #### new version ==>
def ajout_aretes_selon_type_error(row, col, seuil, excludes, matE_reel, matE_proba):
    """
    """
    arete = -1; type_error =None; excludes_tmp = excludes;
    print("row={} et col={}".format(row, col))
    if matE_reel.loc[row,col] == 1 and \
    ((row,col) not in excludes_tmp or (col,row) not in excludes_tmp) and \
    ((seuil-0.1 <= matE_proba.loc[row,col] and  matE_proba.loc[row,col] < seuil) or \
     (seuil-0.1 <= matE_proba.loc[col,row] and  matE_proba.loc[col,row] < seuil)):
        arete = (row,col); type_error = "VraiPositive";
        excludes_tmp.append( (row,col) );
    elif matE_reel.loc[row,col] == 0 and \
        ((row,col) not in excludes_tmp or (col,row) not in excludes_tmp) and \
        ((seuil-0.1 <= matE_proba.loc[row,col] and  matE_proba.loc[row,col] < seuil) or \
         (seuil-0.1 <= matE_proba.loc[col,row] and  matE_proba.loc[col,row] < seuil) ):
        arete = (row,col); type_error == "VraiNegative";
        excludes_tmp.append( (row,col) );
    return arete, excludes_tmp, type_error;

def  distri_01(row, col, seuil, dico_correl, type_errors, matE_proba, matE_reel):
    row, col = verifier_row_col(row, col, matE_proba.columns.tolist(), matE_reel.columns.tolist());
    if row!=None and col!=None:
        if matE_reel.loc[row,col] == 1 and \
            ((seuil-0.1 <= matE_proba.loc[row,col] and  matE_proba.loc[row,col] < seuil) or \
             (seuil-0.1 <= matE_proba.loc[col,row] and  matE_proba.loc[col,row] < seuil)):
                return "VraiPositive", row, col;
        elif matE_reel.loc[row,col] == 0 and \
            ((seuil-0.1 <= matE_proba.loc[row,col] and  matE_proba.loc[row,col] < seuil) or \
             (seuil-0.1 <= matE_proba.loc[col,row] and  matE_proba.loc[col,row] < seuil) ):
                return "VraiNegative", row, col;
        else:
#            print("----({},{})={}, seuil={}----".format(row,col,matE_proba.loc[row,col], seuil))
            pass
    return "error", row, col;
    
def verifier_row_col(row, col, cols_matE_proba, cols_matE_reel):
    """
    """
#    print("row={}, col={}".format(row,col))
    if row.split("->")[0]+"->"+row.split("->")[1] in cols_matE_proba and \
        col.split("->")[0]+"->"+col.split("->")[1] in cols_matE_proba and \
        row.split("->")[0]+"->"+row.split("->")[1] in cols_matE_reel and \
        col.split("->")[0]+"->"+col.split("->")[1] in cols_matE_reel:
#        print("1 verif:{},{} ".format(row.split("->")[0]+"->"+row.split("->")[1], col.split("->")[0]+"->"+col.split("->")[1]))
        return row.split("->")[0]+"->"+row.split("->")[1], col.split("->")[0]+"->"+col.split("->")[1];
    elif row.split("->")[0]+"->"+row.split("->")[1] in cols_matE_proba and \
        col.split("->")[1]+"->"+col.split("->")[0] in cols_matE_proba and \
        row.split("->")[0]+"->"+row.split("->")[1] in cols_matE_reel and \
        col.split("->")[1]+"->"+col.split("->")[0] in cols_matE_reel :
#        print("2 verif:{},{} ".format( row.split("->")[0]+"->"+row.split("->")[1], col.split("->")[1]+"->"+col.split("->")[0]))
        return row.split("->")[0]+"->"+row.split("->")[1], col.split("->")[1]+"->"+col.split("->")[0];
    elif row.split("->")[1]+"->"+row.split("->")[0] in cols_matE_proba and \
        col.split("->")[0]+"->"+col.split("->")[1] in cols_matE_proba and \
        row.split("->")[1]+"->"+row.split("->")[0] in cols_matE_reel and \
        col.split("->")[0]+"->"+col.split("->")[1] in cols_matE_reel :
#        print("3 verif:{},{} ".format( row.split("->")[1]+"->"+row.split("->")[0] ,col.split("->")[0]+"->"+col.split("->")[1]))
        return row.split("->")[1]+"->"+row.split("->")[1], col.split("->")[0]+"->"+col.split("->")[1];
    elif row.split("->")[1]+"->"+row.split("->")[0] in cols_matE_proba and \
        col.split("->")[1]+"->"+col.split("->")[0] in cols_matE_proba and \
        row.split("->")[1]+"->"+row.split("->")[0] in cols_matE_reel and \
        col.split("->")[1]+"->"+col.split("->")[0] in cols_matE_reel :
#        print("4 verif:{},{} ".format( row.split("->")[1]+"->"+row.split("->")[0] ,col.split("->")[1]+"->"+col.split("->")[0]))
        return row.split("->")[1]+"->"+row.split("->")[0], col.split("->")[1]+"->"+col.split("->")[0];
    else:
        return None, None;
def distribution_case_0_1_graphe_reelle(args):
    """
    distribution, selon la metrique de pearson, des cases a 0 et 1
    """
    if args["dbg_0_1"]:
        file = "champlan";
        rep_root = "/home/willy/topologyLearning/datas/"+file+"/";

        sub_graph = "sous_graphe/" #"sous_graphe/"  {"sous_graphe/", ""}
        metrique_distance = "metrique_wil"
        chemin_datasets = rep_root+sub_graph+"datasets/";
        chemin_matrices = rep_root+sub_graph+"matrices/";
        chemin_equipements = rep_root+sub_graph+"matrices_equipements/";
        grandeur = "P";
        matE_proba = pd.read_csv(chemin_matrices+"matrice_adjacence_proba_"+grandeur+".csv",index_col = "Unnamed: 0");
        matE_reel = pd.read_csv(chemin_equipements+"matE_reelles.csv",index_col = "Unnamed: 0");
        methode_correl = "correlationParMorceaux" # ou liste_mode = ["correlationParMorceaux","correlationGlissante"]
        path_save = chemin_equipements;
        correl_vraiNeg = 0.8; correl_vraiPos = 0.2; correl_fauxNeg = 0.7; correl_fauxPos = 0.6;
#        correl_vraiNeg = args["seuil_vrai_negatif"]; correl_vraiPos = args["seuil_vrai_positif"];
        args ={"matE_proba": matE_proba, "matE_reel": matE_reel, "grandeur":grandeur,\
               "methode_correl":methode_correl, "path_save": path_save, \
               "metrique_distance": metrique_distance, "dbg": True,\
               "chemin_datasets":chemin_datasets, "chemin_matrices":chemin_matrices,\
               "correl_vraiNeg":correl_vraiNeg, "correl_vraiPos":correl_vraiPos,\
               "correl_fauxNeg":correl_fauxNeg, "correl_fauxPos":correl_fauxPos};
    
    matE_proba = args["matE_proba"]; 
    matE_reel = args["matE_reel"];
    seuil_correls = np.linspace(0,1,11);
    seuil_correls = [round(x,1) for x in seuil_correls]
    type_errors = ["VraiPositive","VraiNegative"]; dico = dict();
    for seuil_correl in seuil_correls:
        seuil_correl = round(seuil_correl,1)
        dico_correl = {"VraiPositive":[],"VraiNegative":[]};
        for row, col in fct_aux.range_2d(matE_reel.columns.tolist()):
            type_error = "";
            type_error, row, col = distri_01(row, col, seuil_correl, dico_correl, \
                                              type_errors, matE_proba, matE_reel)
            if type_error != "error":
                dico_correl[type_error].append( (row,col) );
                
        dico[seuil_correl] = dico_correl;
    
    df_err_ = pd.DataFrame(dico, columns = [round(s,1) for s in seuil_correls])
    df_err_ = df_err_.T;
    df_err = df_err_.applymap(lambda x: len(x))
    print("df_err = {}, type={} ".format(df_err, type(df_err)))
    
    # plot
    dico_err = {"VraiPositive":1,"VraiNegative":0}
    width = 8; height = 5.5; print("w =", width, " h = ",height)
    fig = plt.figure(figsize=(width*1.1,height*1.1))
    cpt_ax = 1;
    for type_error in dico_err.keys():
        ax1 = fig.add_subplot(1,2,cpt_ax); cpt_ax += 1
        df_err.loc[:,type_error].plot(ax = ax1, kind = "bar", title = "distribution_"+str(dico_err[type_error]))
        ax1.set(ylabel= "nombre d'erreurs", xlabel= "valeurs de correlations", title="distribution_"+str(dico_err[type_error]))
        autolabel(ax1.patches, ax1)
    fig.tight_layout()
    fig.savefig(args["path_save"]+"distribution_0_1"+args["metrique_distance"]+".jpeg",dpi=100) 
    plt.clf();
    
    # comprehension distrib_0 seuil = 0.8 et distrib_1 seuil  = 2, 0.1
    error="VraiNegative"; correl = 1.0; aretes = df_err_.loc[correl, error]; print("aretes={}".format(aretes))
    plot_aretes_distr_01(args["grandeur"], correl, aretes[:8], args["chemin_datasets"], args["path_save"], error,args["date_debut"],args["date_fin"])
    
    error="VraiPositive"; correl = 0.1; aretes = df_err_.loc[correl, error]; print("aretes={}".format(aretes))
#    plot_aretes_distr_01(args["grandeur"], correl, aretes[:8], args["chemin_datasets"], args["path_save"], error,args["date_debut"],args["date_fin"])
    
def plot_aretes_distr_01(grandeur, correl, pair_equips, chemin_datasets, path_save, error, date_debut, date_fin):
    width = 8; height = 5.5; fig = plt.figure(figsize=(width*1.15,height*1.15))
    fig.set_size_inches( (width*1.15, height*1.15) )\
    if len(pair_equips) < 9 \
    else fig.set_size_inches( (width*2.75, height*2.75) ); #fig.set_size_inches( (width*2.15, height*2.15) );
    cpt_ax = 1;
    w,h = sous_plot(len(pair_equips));
    print("w={},h={}, len pair={}".format(w,h,len(pair_equips)))
    if w == None:
        print("correl={} have 0 pair_equips=aretes".format(correl))
        return;
    df = pd.read_csv(chemin_datasets+"dataset_"+grandeur+".csv");
    df = pd.read_csv(chemin_datasets+"shapelet_dataset_"+grandeur+".csv");
    df = df.set_index("timestamp") if "timestamp" in df.columns.tolist() \
                                    else df.set_index("Unnamed: 0");
    df = df.loc[date_debut:date_fin] if date_debut != 0 and date_fin != 0 else df;
    df = (df - df.mean(skipna=True))/ df.std(skipna=True); # df = normalize(df.copy())
    df.rolling(window=20).mean()
    df.fillna(0, inplace=True)
    for pair in pair_equips:
        print("pair={}".format(pair))
        ax = fig.add_subplot(w,h,cpt_ax); cpt_ax += 1;
        df[list(pair)].plot(ax = ax)
        ax.set(title = pair)
    fig.tight_layout();
    fig.savefig(path_save+"distrib01_"+error+"_correl_"+str(correl)+"_error_<"+error+">_grandeur_"+grandeur+".jpeg",dpi=70); 
    plt.clf();
    pass    
#### distribution des cases de matEreelles  a 0 et a 1  ===> FIN

if __name__ == '__main__':
    
     start= time.time();
     nbre_lien = 5; dimMatA = 20;
     methode_delete_add_edges = 0; # 1:selection aretes a delete or a ajouter par liste, 0: par proba
     p_correl = 0.5; # valeur de correlation ([0,1]) a partir duquel 
                           # *  correl = 0 est transforme en 1 et correl = 1 est transforme a 0
                           # *  on attribue des valeurs de correlations dans matE.
                           # p_correl = [0.0, 0.1, 0.2, ...., 1.0]
     correl_seuil = 0.7 # seuil a partir duquel on a des correlations fausses positives (0->1) et fausses negatives (1->0)
                        # [ signification des correlations fausses {positives,negatives} a partir de MATE]
     loi_stats = "uniforme";#"poisson";
     SEUIL_PROBA = 0.8; # je pense cest pour ajouter des probas a chaque arete # A EFFACER
     algoGreedy = False #True; # False correction tous les noeuds a -1, True: algo Greedy   
     biais = False;
     number_permutations_nodes_1= 10 #100;
     facteur_multiplicatif = 1; exposant = 1; # 0: fct_unitaire, 1:fct_normal, 2: fct_quadratique, 4:fct_quadruple, 5: fct_quintuple
     type_fct_cout = "cloche" # ou "lineaire"
     coef_fct_cout = (exposant, facteur_multiplicatif, type_fct_cout)
     mode_select_noeuds_1 = "aleatoire" #"degreMin" #"coutMin" # "degreMin" # aleatoire
     number_items_pi1_pi2 = 1; k = 50;
     arg_params = {"number_items_pi1_pi2": number_items_pi1_pi2, \
                   "nbre_lien":nbre_lien,\
                   "dimMatA":dimMatA, "k":k,\
                   "number_permutations_nodes_1": number_permutations_nodes_1, \
                   "methode_delete_add_edges": methode_delete_add_edges, \
                   "loi_stats": loi_stats,\
                   "biais": biais, \
                   "p_correl": p_correl,\
                   "correl_seuil": correl_seuil,\
                   "algoGreedy":algoGreedy, \
                   "mode_select_noeuds_1":mode_select_noeuds_1,\
                   "coef_fct_cout":coef_fct_cout};
                   
                   
#     test_distribution(arg_params)
     ti = time.time() - start    
     
     ### distribution dataReels
     rep_root = "/home/willy/topologyLearning/datasets/Champlan/matrice_equipements/"
     path_datasets = "";
     path_matE_proba = rep_root+"sub_matE.csv";
     path_matE_reel = rep_root+"sub_matE_reelles.csv";
     path_cliques = rep_root+"seuil_aretes_LG.json";
     path_save = rep_root;
     params = {"path_datasets": path_datasets,\
               "path_matE_proba": path_matE_proba,\
               "path_matE_reel": path_matE_reel,\
               "path_cliques": path_cliques,\
               "path_save": path_save}
     file = "champlan"; sub_graph = "sous_graphe/" # "sous_graphe/" or ""
#     file = "velizy"; sub_graph = "" # "sous_graphe/" or ""
     rep_root = "/home/willy/topologyLearning/datas/"+file+"/"+sub_graph;
     chemin_datasets = rep_root+"datasets/";
     chemin_matrices = rep_root+"matrices/";
     chemin_equipements = rep_root+"matrices_equipements/";
     path_matE_proba = rep_root+"matrices/"+"sub_matE.csv";
     path_matE_reel = rep_root+"matrices/"+"sub_matE_reelles.csv";
     path_cliques = rep_root+"matrices_equipements/"+"seuil_aretes_LG.json";
     path_save = rep_root+"matrices_equipements/";

     params = {"chemin_datasets":chemin_datasets,\
               "chemin_matrices":chemin_matrices,\
               "chemin_equipements":chemin_equipements,
               "path_datasets": path_datasets,\
               "path_matE_proba": path_matE_proba,\
               "path_matE_reel": path_matE_reel,\
               "path_cliques": path_cliques,\
               "path_save": path_save}
     metrique_distances = ["lb_keogh","pearson","lcs"];
     for metrique_distance in metrique_distances:
         params["metrique_distance"] = metrique_distance;
         params["path_matE_proba"] = rep_root+"matrices/"+"matE_"+metrique_distance+".csv";
         params["path_cliques"] = rep_root+"matrices_equipements/"+"seuil_aretes_LG_"+metrique_distance+".json";
#         distribution_dataReelles(params) # DECOMMENTER
     
     ### distribution valeurs {fausses, vrai} {negatives, positives} correlations ####
     matE_proba = pd.read_csv(chemin_equipements+"matE.csv",index_col = "Unnamed: 0");
     matE_reel = pd.read_csv(chemin_matrices+"matE_reelles.csv",index_col = "Unnamed: 0");
     methode_correl = "correlationParMorceaux" # ou liste_mode = ["correlationParMorceaux","correlationGlissante"]
     path_save = chemin_equipements;
     params["matE_proba"]= matE_proba; params["matE_reel"]= matE_reel;
     params["path_save"] = path_save; params["methode_correl"] = methode_correl; 
     params["dbg"] = False; 
     metrique_distances = ["lb_keogh","pearson","lcs"];
     metrique_distances = ["metrique_wil"]
     for metrique_distance in metrique_distances:
         params["metrique_distance"] = metrique_distance;
         params["matE_proba"] = pd.read_csv(chemin_matrices+"matE_"+metrique_distance+".csv",\
                                             index_col = "Unnamed: 0");
         params["matE_proba"] = pd.read_csv(chemin_matrices+"matrice_adjacence_proba_P.csv",\
                                             index_col = "Unnamed: 0"); # pour test TO DELETE AFTER
#         distribution_valeurs_correlation_selon_seuil_methode(params)
     ### distribution valeurs {fausses, vrai} {negatives, positives} correlations ####

     ### formalisation matE
     metrique_distances = ["lb_keogh","pearson","lcs"]
     metrique_distances = ["metrique_wil"]
     params["metrique_distance"] = metrique_distances;
     params["dbg_0_1"] = False;
     distrib_formalisation_matE(params)
     ### formalistaion matE
     
     ##### comprendre correlation trouvee par l'algorithme. ==> debut
     nodes = ["TGBT2","GF1","R495","CVC2"]#["TGBT1", "TGBT2"]; #["R042", "R495"]; #["DD108", "TGBT4"];
     nodes = ["TGBT2","GF1","R495"]
     derivates = True#False#True
#     tracer_courbes_equipementCorreles(nodes, params, derivates)
     ##### comprendre correlation trouvee par l'algorithme. ==> FIN
     
#     plot_equipements_selon_correlations({"dbg":True})
#     distribution_case_0_1_graphe_reelle({"dbg_0_1":True})
     print (time.time() - start)