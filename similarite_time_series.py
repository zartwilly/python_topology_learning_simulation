#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 19 18:38:20 2018

ce fichier est une particularite (amelioration) du fichier correlation.py
il est adapte pour le calcul des similarite avec la distance de Pearson
TODO : ajouter les autres distances tel que 
    * MSM, DTW, LCSS, TWE 
@author: willy
"""
REP_ROOT = "/home/willy/Documents/python_topology_learning_simulation/" 


import re, os, time, logging;
import ujson, math;
import pandas as pd
import numpy as np

from pathlib import Path;
from functools import reduce;
from operator import mul
from ast import literal_eval

import creation_datasets as create_datasets;
import extractGraphReel as graphReel;
import distances_similarite as distance_similarite;
import ajout_mesures_correlations as ajusterCorrelation;
import sax_encoding as sax_encod;
import shapelets_transformV1 as shapelet;


import fonctions_auxiliaires as fct_aux;

import courbes as courbes

####---------------------------------------------------------------------------
##              correlation par shapelets
###----------------------------------------------------------------------------
##### correlation par shapelets =====>
def transformed_dataset(chemin_dataset, grandeur, fenetre_slicing, date_debut, date_fin):
    df_gr = pd.read_csv(chemin_datasets+"dataset_"+grandeur+".csv");
    df_gr = df_gr.set_index("timestamp") if "timestamp" in df_gr.columns \
                                         else df_gr.set_index("Unnamed: 0");
    columns = df_gr.columns.tolist();
    df_gr = pd.read_csv(chemin_datasets+"shapelet_dataset_"+grandeur+".csv");
    df_gr = df_gr.set_index("timestamp") if "timestamp" in df_gr.columns \
                                         else df_gr.set_index("Unnamed: 0");
    cols_to_delete = list(df_gr.columns[df_gr.isnull().all()])
    df_gr.drop(cols_to_delete, axis = 1, inplace = True);
    if df_gr.empty:
        return df_gr;
    df_gr_norm = (df_gr - df_gr.mean(skipna=True))/df_gr.std(skipna=True);
    df_gr_norm.fillna(method='pad',inplace = True);
    df_gr_norm = df_gr_norm.rolling(window=fenetre_slicing).mean();
    df_gr_norm = df_gr_norm.loc[date_debut:date_fin]
    return df_gr_norm, columns;
    
def correlation_by_shapelets_grandeur(args):
    
    df_gr_norm,cols = transformed_dataset(args["chemin_datasets"], args["grandeur"], \
                                          args["fenetre_slicing"], args["date_debut"],\
                                          args["date_fin"]);                                 
    logging.debug("correlation_by_shapelet_grandeur : %s columns = %s",args["grandeur"], df_gr_norm.columns.tolist())
    if df_gr_norm.empty :
        dico_gp = dict();
        matE_gp = pd.DataFrame(columns=cols,index=cols)
        matE_gp.fillna(0, inplace = True);
        logging.debug("correlation_by_shapelet_grandeur : matE_(%s) is empty: all columns are empty", args["grandeur"])
        return dico_gp, matE_gp;

    label_aretes = shapelet.label_arcs(df_gr_norm.columns.tolist());
    data = [];
    for arete, label in label_aretes.items():
        data.append((df_gr_norm[arete].values, label))
    shapelet_dict = dict();
#    shapelet_dict = shapelet.extract_shapelets_v1(data, args["min_len"], args["max_len"]);
    shapelet_dict = shapelet.extract_shapelets(data, args["min_len"], \
                                               args["max_len"]);
    df_data = pd.DataFrame(shapelet_dict);
    #delete  empty columns
#    df_data.drop(list(df_data.columns[df_data.isnull().all()]), axis = 1, inplace = True);
    columns = df_gr_norm.columns;
    
#    return df_data;
    ## normalized df_data
    print("df_data: mean={}, std={}".format(df_data.mean(), df_data.std()))
    df_data = distance_similarite.normalize(df_data) #(df_data - df_data.mean()) / df_data.std();
    
    ## correlation par une  metrique choisi
    dico_gp = dict();
    matE_gp = pd.DataFrame(columns=cols,index=cols);
    for arete1, arete2 in fct_aux.range_2d(columns):
        correl = distance_similarite.distance_pearson(df_data[arete1], \
                                                      df_data[arete2]);
        matE_gp.loc[arete1, arete2] = abs(correl);
        matE_gp.loc[arete2, arete1] = abs(correl);
        dico_gp[(arete1,arete2)] = abs(correl);
        logging.debug("correlation_by_shapelet_grandeur : (%s, %s) = %s", \
                      arete1, arete2, abs(correl))
            
    if "modifCorrel" in args.keys() and args["modifCorrel"]:
        args["matE_gp"] = matE_gp;
        matE_gp = ajusterCorrelation.ajouter_mesures_correlations(args);
        dico_gp_tmp = dict();
        for arc1, arc2 in fct_aux.range_2d(matE_gp.columns):
            dico_gp_tmp[(arc1,arc2)] = matE_gp.loc[arc1, arc2]
        dico_gp = dico_gp_tmp;
    matE_gp.fillna(0, inplace = True);
    matE_gp.to_csv(args["chemin_matrices"]+"matrice_adjacence_proba_"+\
                   args["grandeur"]+".csv");
    analyse(dico_gp, args["grandeur"],args,critere=0.9);
    print("grandeur {} termine".format( args["grandeur"]) )
    logging.debug("correlation_by_shapelet_grandeur : grandeur %s termine", \
                  args["grandeur"])
    
    return dico_gp, matE_gp;
    
def transformed_dataset_v1(chemin_dataset, grandeur, fenetre_slicing, date_debut, date_fin):
    df_gr = pd.read_csv(chemin_datasets+"dataset_"+grandeur+".csv");
    df_gr = df_gr.set_index("timestamp") if "timestamp" in df_gr.columns \
                                         else df_gr.set_index("Unnamed: 0");
    cols = df_gr.columns.tolist();
    cols_to_delete = list(df_gr.columns[df_gr.isnull().all()])
    df_gr.drop(cols_to_delete, axis = 1, inplace = True);
    if df_gr.empty:
        return df_gr;
    df_gr_norm = (df_gr - df_gr.mean(skipna=True))/df_gr.std(skipna=True);
    df_gr_norm.fillna(method='pad',inplace = True);
    df_gr_norm = df_gr_norm.rolling(window=fenetre_slicing).mean();
    df_gr_norm = df_gr_norm.loc[date_debut:date_fin]
    return df_gr_norm, cols;
    
def correlation_by_shapelets_grandeur_v1(args):
    
    df_gr_norm,cols = transformed_dataset_v1(args["chemin_datasets"], \
                                             args["grandeur"], \
                                             args["fenetre_slicing"], \
                                             args["date_debut"],\
                                             args["date_fin"]);                                    
    logging.debug("correlation_by_shapelet_grandeur : %s columns = %s",\
                  args["grandeur"], df_gr_norm.columns.tolist())
    if df_gr_norm.empty :
        dico_gp = dict();
        matE_gp = pd.DataFrame(columns=cols,index=cols)
        matE_gp.fillna(0, inplace = True);
        logging.debug("correlation_by_shapelet_grandeur : matE_(%s) is empty: all columns are empty", args["grandeur"])
        return dico_gp, matE_gp;
     
    df_data = pd.read_csv(chemin_datasets+"shapelet_dataset_"+args["grandeur"]+".csv");
    df_data = df_data.set_index("timestamp") if "timestamp" in df_data.columns \
                                             else df_data.set_index("Unnamed: 0");
    #delete  empty columns
#    df_data.drop(list(df_data.columns[df_data.isnull().all()]), axis = 1, inplace = True);
    
    ## normalized df_data
#    print("df_data: mean={}, std={}".format(df_data.mean(), df_data.std()))
    df_data = distance_similarite.normalize(df_data) #(df_data - df_data.mean()) / df_data.std();
    
    ## correlation par une  metrique choisi
    dico_gp = dict();
    matE_gp = pd.DataFrame(columns=cols,index=cols);
    columns = df_gr_norm.columns;
    for arete1, arete2 in fct_aux.range_2d(columns):
        correl = distance_similarite.distance_pearson(\
                        df_data[arete1], df_data[arete2]);
        matE_gp.loc[arete1, arete2] = abs(correl);
        matE_gp.loc[arete2, arete1] = abs(correl);
        dico_gp[(arete1,arete2)] = abs(correl);
        logging.debug("correlation_by_shapelet_grandeur : (%s, %s) = %s", \
                      arete1, arete2, abs(correl))
    matE_gp.fillna(0, inplace=True)        
    if "modifCorrel" in args.keys() and args["modifCorrel"]:
        args["matE_gp"] = matE_gp;
        matE_gp = ajusterCorrelation.ajouter_mesures_correlations(args);
        dico_gp_tmp = dict();
        for arc1, arc2 in fct_aux.range_2d(matE_gp.columns):
            dico_gp_tmp[(arc1,arc2)] = matE_gp.loc[arc1, arc2]
        dico_gp = dico_gp_tmp;
    
    matE_gp.fillna(0, inplace = True);
    matE_gp.to_csv(args["chemin_matrices"]+"matrice_adjacence_proba_"+\
                   args["grandeur"]+".csv");
    analyse(dico_gp, args["grandeur"],args,critere=0.9);
    print("grandeur {} termine".format( args["grandeur"]) )
    logging.debug("correlation_by_shapelet_grandeur : grandeur %s termine", \
                  args["grandeur"])
    
    return dico_gp, matE_gp;
#### ---------------- correlation par shapelets ----------------------> FIN



####---------------------------------------------------------------------------
##              correlation par grandeur
####---------------------------------------------------------------------------
### correlation par grandeur ===> DEBUT
def remap_keys(mapping):
    return [{'key':k, 'value': v} for k, v in mapping.items()]

def calcul_correlation_slicing( ts1, ts2, metrique_distance, fenetre):
    """
    calculer la correlation selon la metrique
    """
    if metrique_distance == "pearson":
        return distance_similarite.correlation_pearson(ts1,ts2, fenetre)
    elif metrique_distance == "fdtw":
        return distance_similarite.DTW_distance_correlation(ts1, ts2, fenetre);
    elif metrique_distance == "lcs":
        return distance_similarite.LCS_correlation(ts1, ts2)
    elif metrique_distance == "distance_pearson":
        return distance_similarite.distance_pearson(ts1, ts2)
    elif metrique_distance == "metrique_pearson_damien":
        return distance_similarite.metrique_pearson_damien(ts1, ts2, fenetre)
    elif metrique_distance == "metrique_wil_histo":
        return distance_similarite.metrique_wil_histo(ts1,ts2)
    elif metrique_distance == "twed":
        return distance_similarite.twed_correlation(ts1.values, ts2.values, 0.1, 0.1)
    elif metrique_distance  == "lb_keogh":
        return distance_similarite.LB_Keogh_correlation(ts1, ts2, fenetre)
    else: 
        print("metrique_distance: {} dont exist".format(metrique_distance))
        logging.debug("metrique_distance: %s dont exist",metrique_distance)
        return 0.0;
    
    
def value_0_1(correl, correl_max):
    if np.isnan(correl)  and np.isnan(correl_max):
        return 0.0
    elif np.isnan(correl_max):
        return 0.0;
    elif correl < correl_max:
        return 1-correl/correl_max
    elif correl > correl_max:
        return 1-correl_max/correl;
    elif correl == correl_max:
        logging.debug("correl==correl_max")
        return 1
    return 0.0
    
def correlation_1e_grandeur(grandeur, args):
    """
    determiner la correlation associe a une grandeur
    """
    df_gr = pd.read_csv(args["chemin_datasets"]+"dataset_"+grandeur+".csv");
    if "timestamp" in df_gr.columns.tolist():
        df_gr = df_gr.set_index("timestamp");
    else:
        df_gr = df_gr.set_index("Unnamed: 0")
    
    logging.debug("correlation_1e_grandeur = %s columns = %s",\
                  grandeur, df_gr.columns.tolist())
    cols = df_gr.columns.tolist();
    matE_gp = pd.DataFrame(index = cols, columns = cols)
    
    #delete empty columns in df_gr
    cols_to_delete = list(df_gr.columns[df_gr.isnull().all()])
#    df_gr.drop(cols_to_delete, axis = 1, inplace = True);
    
    dico_gp = dict();
    if df_gr.empty :
        matE_gp.fillna(0, inplace = True);
        logging.debug("matE_(%s) is empty: all columns are empty", grandeur)
        return dico_gp, matE_gp;
    else:
        df_gr_norm = (df_gr - df_gr.mean(skipna=True))/df_gr.std(skipna=True); #df_gr_norm.fillna(0,inplace = True);
        df_gr_norm.fillna(method='pad',inplace = True);
        df_gr_norm = df_gr_norm.rolling(window=20).mean();
        print("*** df_gr_norm = {}".format(df_gr_norm.columns.tolist()))
        for arc1, arc2 in fct_aux.range_2d(df_gr_norm.columns.tolist()):
            correl = 0; correl_tmp = 0; 
            correl_max=0; corr_max_arc1 = 0; corr_max_arc2 = 0;
            correl_tmp = calcul_correlation_slicing(df_gr_norm[arc1], \
                        df_gr_norm[arc2], args["metrique_distance"], \
                        args["fenetre"])
            if args["metrique_distance"] not in ["pearson","distance_pearson",\
            "metrique_wil_histo","metrique_pearson_damien"]:
                corr_max_arc1 = calcul_correlation_slicing(df_gr_norm[arc1], \
                                df_gr_norm[arc1].ix[args["fenetre"]:], \
                                args["metrique_distance"], \
                                args["fenetre"])
                corr_max_arc2 = calcul_correlation_slicing(df_gr_norm[arc2], \
                                df_gr_norm[arc2].ix[args["fenetre"]:], \
                                args["metrique_distance"], \
                                args["fenetre"])
                correl_max = max(corr_max_arc1, corr_max_arc2)
                correl_tmp = value_0_1(correl_tmp, correl_max)
            correl = correl_tmp;
            if (arc1 == "GF2" or arc2 == "GF2") and \
                (arc1 == "R495" or arc2 == "R495"):
                print("correl={}, arc1={}, correl_max_arc1={}, arc2={}, \
                      correl_max_arc2={}, correl_tmp={}"\
                      .format(correl,arc1,corr_max_arc1,arc2,\
                              corr_max_arc2,correl_tmp))
            if (arc1 == "R495" or arc2 == "R495") and \
                (arc1 == "TGBT4" or arc2 == "TGBT4"):
                print("correl={}, arc1={}, arc2={}".format(correl,arc1,arc2))
            dico_gp[(arc1,arc2)] = abs(correl);
            matE_gp.loc[arc1, arc2] = abs(correl); 
            matE_gp.loc[arc2, arc1] = abs(correl);
            
            logging.debug("(%s, %s) = %s ,correl_tmp=%s,correl_max=%s matE_grand = %s ,count: arc1 = %s, arc2 = %s", \
                          arc1, arc2, correl, correl_tmp, correl_max, matE_gp.loc[arc1,arc2], \
                          df_gr_norm[arc1].count(), df_gr_norm[arc2].count())
        matE_gp.fillna(0, inplace = True);
        args["grandeur"] = grandeur;
        if "modifCorrel" in args.keys() and args["modifCorrel"]:
            args["matE_gp"] = matE_gp;
            matE_gp = ajusterCorrelation.ajouter_mesures_correlations(args);
            dico_gp_tmp = dict();
            for arc1, arc2 in fct_aux.range_2d(matE_gp.columns):
                dico_gp_tmp[(arc1,arc2)] = matE_gp.loc[arc1, arc2]
            dico_gp = dico_gp_tmp;
        matE_gp.to_csv(args["chemin_matrices"]+\
                       "matrice_adjacence_proba_"+grandeur+".csv");
        analyse(dico_gp, grandeur,args,critere=0.9);
        print("grandeur ", grandeur, " termine")
        return dico_gp, matE_gp;

def analyse(dico_gp, grandeur, args,critere ):
    dico_res = dict();
    criteres = np.linspace(0,1,11)
    for critere in criteres:
        liste_aretes = [];
        for arete, correl in dico_gp.items():
            if correl <= critere and correl >= critere-0.1:
                liste_aretes.append(arete)
        dico_res[critere] = liste_aretes;
    with open(args["chemin_equipements"]+'analyse_grandeur_'+\
              grandeur+'_'+args["metrique_distance"]+".json", 'w') as outfile:
        ujson.dump(remap_keys(dico_res), outfile, sort_keys = True, indent = 4);    
    pass        
def calcul_matE(dico_matE, args):
    """
    formalisation de matE
    """
    # TODO lire le dico_matE a partir d'un fichier, ajouter un dbg_calcul_matE
#    if args["dbg_calcul_matE"]:
#        dico_matE = ujson.load(open(args["chemin_matrices"]+"dico_matE_"+args["metrique_distance"]+".json","r"));
    dico_arcs_corrs = dict();
    for grand, dico_gr in dico_matE.items():
        for aretes, corr in dico_gr.items():
            if aretes not in dico_arcs_corrs.keys():
                dico_arcs_corrs[aretes] = [corr];
            else:
                dico_arcs_corrs[aretes].append(corr);
#    print("dico_arcs_corrs = ", dico_arcs_corrs)
    
    cols = set([ x for sub in dico_arcs_corrs.keys() for x in tuple(sub) ])
    matE = pd.DataFrame(index = cols, columns = cols)
    for aretes, correls in dico_arcs_corrs.items():
#        print("aretes = {}, tuple = {}, aretes = {}".format( type(aretes), type(tuple(aretes)), tuple(aretes)))
#        print("aretes[0] = {}, aretes[1] = {}".format( tuple(aretes)[0], tuple(aretes)[1] ) )
        if args["type_cal_matE"] == "max":
            # max correl 0
            matE.loc[aretes[0], aretes[1]] = \
                sorted(correls,reverse=True)[args["indice_max_correlation"]]
            matE.loc[aretes[1], aretes[0]] = \
                sorted(correls,reverse=True)[args["indice_max_correlation"]]
        elif args["type_cal_matE"] == "mean":
            # mean
            matE.loc[aretes[0], aretes[1]] = np.mean(correls);
            matE.loc[aretes[1], aretes[0]] = np.mean(correls);
        elif args["type_cal_matE"] == "mult":
            # produit
            matE.loc[aretes[0], aretes[1]] = reduce(mul, correls, 1);
            matE.loc[aretes[1], aretes[0]] = reduce(mul, correls, 1);
    matE.fillna(0, inplace = True);
    matE.to_csv(args["chemin_matrices"]+"matE_"+args["metrique_distance"]+".csv");
                
    with open(args["chemin_matrices"]+'correlation_grandeurs_'+\
              args["metrique_distance"]+".json", 'w') as outfile:
        ujson.dump(remap_keys(dico_arcs_corrs), outfile, sort_keys = True, \
                   indent = 4);
    return matE;

#### calcul matE mean, produit, max_correl ==> debut
def calcul_matE_mean_prod_max(dico_matE, args):
    if args["dbg_calcul_matE"]:
        dico_matE = ujson.load(open(args["chemin_matrices"]+"dico_matE_"+\
                                    args["metrique_distance"]+".json","r"));
    dico_arcs_corrs = dict();
    for grand, dico_gr in dico_matE.items():
        for aretes, corr in dico_gr.items():
#            print("{}={}".format(aretes, corr))
            if aretes not in dico_arcs_corrs.keys():
                dico_arcs_corrs[aretes] = [corr];
            else:
                dico_arcs_corrs[aretes].append(corr);
#    return dico_arcs_corrs
    cols = set([ x for sub in dico_arcs_corrs.keys() \
                for x in tuple(literal_eval(sub)) ])
    for type_cal in ["mean","mult","max"]:
        matE = pd.DataFrame(index = cols, columns = cols)
        for aretes, correls in dico_arcs_corrs.items():
            aretes = literal_eval(aretes);
#            print("{}=>{},{}".format(aretes,aretes[0],aretes[1]))
            if type_cal == "max":
                # max correl 0
                matE.loc[aretes[0], aretes[1]] = \
                    sorted(correls,reverse=True)[args["indice_max_correlation"]]
                matE.loc[aretes[1], aretes[0]] = \
                    sorted(correls,reverse=True)[args["indice_max_correlation"]]
            elif type_cal == "mean":
                # mean
                matE.loc[aretes[0], aretes[1]] = np.mean(correls);
                matE.loc[aretes[1], aretes[0]] = np.mean(correls);
            elif type_cal == "mult":
                # produit
                matE.loc[aretes[0], aretes[1]] = reduce(mul, correls, 1);
                matE.loc[aretes[1], aretes[0]] = reduce(mul, correls, 1);
        matE.fillna(0, inplace = True);
#        matE.to_csv(args["chemin_matrices"]+"matE_"+args["metrique_distance"]+".csv");
        matE_reel = pd.read_csv(args["chemin_matrices"]+\
                                "matE_reelles.csv",index_col = "Unnamed: 0");
        args["matE_proba"] = matE; args["matE_reel"] = matE_reel;
        args["path_save"] = args["chemin_equipements"];
        args["type_cal_matE"] = type_cal;
#        return matE.columns.tolist(), matE_reel.columns.tolist()
#        return matE, matE_reel
        courbes.distribution_valeurs_correlation_selon_seuil_methode(args)
        
        # distribution des valeurs de correlation selon les cases (0,1) dans le graphe reel(matE_reelles)
        courbes.distribution_case_0_1_graphe_reelle(args)
        
        # distribution formalisation matE  
        courbes.distrib_formalisation_matE(args)
        
#### calcul matE mean, produit, max_correl ==> fin
    
def correlation_par_grandeurs(args):
    """
    obtient la correlation pour chaque grandeur
    puis calcul matE
    """
    grandeurs = args["selected_grandeurs"];
    if len(args["selected_grandeurs"]) == 0:
        grandeurs = fct_aux.liste_grandeurs(args["chemin_datasets"])
    
    dico_matE = dict(); 
    for grandeur in grandeurs:
        if args["metrique_distance"] == "sax":
           dico_matE[grandeur] = sax_encod.create_matE_grandeur(grandeur, args["epsilon_sax"], \
                                   args["chemin_datasets"], args["chemin_matrices"])
            
        elif args["metrique_distance"] == "metrique_shapelets":
            args["grandeur"] = grandeur; 
            dico_gp, matE_gp = correlation_by_shapelets_grandeur_v1(args)
            dico_matE[grandeur] = dico_gp;
        else:
            dico_gp, matE_gp = correlation_1e_grandeur(grandeur, args);
            dico_matE[grandeur] = dico_gp;    
    with open(args["chemin_matrices"]+'dico_matE_'+args["metrique_distance"]+".json", 'w') as outfile:
        ujson.dump(dico_matE, outfile, sort_keys = True, indent = 4);
    matE = calcul_matE(dico_matE, args)
    return matE;

####---------------------------------------------------------------------------
##              correlation par metriqueFusion
####---------------------------------------------------------------------------
def correlation_par_metriqueFusion(chemin_mat_equipement, fenetre, mode, args):
    """
    chemin_mat_equipement = chemin de la matrice d'equipement
    mode = mode de correlation 
    return mat_R: matrice de correlation;
    """
    matrices_equipment = distance_similarite.lire_matrice_equipements(\
                            chemin_mat_equipement);
    mat_R = pd.DataFrame(index = matrices_equipment, columns = matrices_equipment)
    for nom_matA, nom_matB in fct_aux.range_2d(matrices_equipment):
#        matA = pd.read_csv(chemin_mat_equipement+"matrice_"+nom_matA+".csv", index_col = "Unnamed: 0")
#        matB = pd.read_csv(chemin_mat_equipement+"matrice_"+nom_matB+".csv", index_col = "Unnamed: 0")
        ###
        matA = pd.read_csv(chemin_mat_equipement+"matrice_"+nom_matA+".csv")
        matB = pd.read_csv(chemin_mat_equipement+"matrice_"+nom_matB+".csv")
        index_data_matA = ""; index_data_matB = "";
        if "timestamp" in matA.columns.tolist():
            matA = matA.set_index("timestamp");
            index_data_matA = "timestamp";
        else:
            matA = matA.set_index("Unnamed: 0")
            index_data_matA = "Unnamed: 0";
        if "timestamp" in matB.columns.tolist():
            matB = matB.set_index("timestamp");
            index_data_matB = "timestamp";
        else:
            matB = matB.set_index("Unnamed: 0")
            index_data_matB = "Unnamed: 0";
        ###
        # derivate: gradient ici
        if args["derivate"] == True:
            matA = pd.DataFrame(np.gradient(matA.values, args["interval_deriv"])[1], \
                                            columns=matA.columns.tolist(), index = index_data_matA)
            matB = pd.DataFrame(np.gradient(matB.values, args["interval_deriv"])[1], \
                                            columns=matB.columns.tolist(), index = index_data_matB)
        correlations = list()
        for ts_a, ts_b in distance_similarite.ts_slicing(matA, matB, fenetre, mode):
            r_ab_tmp = abs(distance_similarite.metriqueFusion_correlation(ts_a,\
                            ts_b));
            print("ts_a = {}, ts_b  = {}, r_ab_tmp = {}"\
                  .format(len(ts_a), len(ts_b), r_ab_tmp))
#            print("ts_a = {}, \nts_b  = {}".format(ts_a.index, ts_b.index) )
            correlations.append(r_ab_tmp)
        r_ab = distance_similarite.metrique_slicing(correlations)
        logging.debug("matA =  %s, matB = %s, r_ab = %s, correls=%s ", nom_matA, nom_matB, r_ab, correlations)
        
        mat_R.loc[nom_matA][nom_matB] = math.fabs(r_ab);
        mat_R.loc[nom_matB,nom_matA] = math.fabs(r_ab);
        print("(",nom_matA,",",nom_matB,"):",r_ab)
#        break
    mat_R.fillna(0, inplace = True)
    mat_R.to_csv(chemin_mat_equipement+"matE.csv")
    return mat_R;

####--------------------------------------------------------------------------- 
##              correlation time series 
####--------------------------------------------------------------------------- 
def correlation_time_series(chemin_datasets, chemin_matrices, \
                            chemin_matrice_equipements, args):
    """ 
    calcule la correlation entre tous les grandeurs selon le mode
    """
    #logging
    logging.basicConfig(filename=chemin_matrice_equipements+'equipement_by_matrice.log',\
                        filemode='w',level=logging.DEBUG, format='%(asctime)s %(message)s' )
    
    path_mat_equipment = Path(chemin_matrice_equipements)
    path_mat_equipment.mkdir(parents=True, exist_ok=True)
    
    matE = pd.DataFrame();
    if args["methode_correl"] == "correlation_par_metrique_fusion":
        bool_equip = True
        if bool_equip:
            matE = correlation_par_metriqueFusion(chemin_matrice_equipements, args["fenetre"], \
                                                  args["mode_correlation"], args)
        print(" **** creation matE===> FINI")
        logging.debug(" **** creation matE===> FINI")
    elif args["methode_correl"] == "correlation_par_grandeur":
        args["chemin_datasets"] = chemin_datasets; args["chemin_matrices"] = chemin_matrices; 
        args["chemin_equipements"] = chemin_matrice_equipements; 
        matE = correlation_par_grandeurs(args); 
    else:
        print("ERROR: methode_correl unknown".upper())
        logging.debug("ERROR => methode_correl unknown ".upper())
    
    return matE;

if __name__ == '__main__':
    
    start= time.time();
    bool_creation_datasets = False;
    sousGraphe = True;
    reseau = "champlan"; fichier_mesures = None; date_debut=None; date_fin=None; rep=None
    if reseau == "champlan":
        fichier_mesures = "datasets_"+reseau+".json"; 
        date_debut = 1359390657; date_fin = 1359477057; dbg = False; 
        rep = reseau+"_newDataset"; # or champlan 
    elif reseau == "velizy":
        fichier_mesures = "datasets_"+reseau+".json"; 
        date_debut = 1436139572; date_fin = 1436306372; dbg = False; 
        rep = reseau+"_newDataset"; # velizy
    else:
        print("ERROR reseau {} don't exist".format(reseau))
    root = REP_ROOT+"dataReels/"+rep+"/";
    root = root+"sous_graphe/" if sousGraphe else root;
    chemin_datasets = root+"datasets/"; chemin_matrices = root+"matrices/";
    chemin_equipements = root+"matrices_equipements/";
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
    modifCorrel = False # TODO A EFFACER LORS ECRITURE CHAPITRE CORRELATION
    
    selected_grandeurs = list();
    if reseau == "velizy":
        ''' 
        "avgf","avgi","avgv","avgu31","avgu12","avgu23","avgin",\
        "avgi1","avgi2","avgi3",\
        "avgqneg","avgppos","avgs","hinrow3","hinrow5","hinrow7","hinrow9",\
        "hi1row7","hi1row9","hi1row3","hi1row5",\
        "hi2row7","hi2row9","hi2row3","hi2row5",\
        "hi3row7","hi3row9","hi3row3","hi3row5",\
        "thdin","thdi1","thdi2","thdi3","thdu12","thdu23","thdu31" '''
        selected_grandeurs = ['AvgI1', 'AvgI']
        selected_grandeurs = ['avgf']; 
    else:
        selected_grandeurs = ["P","I","U","S","E","U1","U2","U3","I1","I2","I3"]#["P","I"]
        selected_grandeurs = ["P","S","U1","U2","U3","I1","I2","I3"]#["P","I"]
        selected_grandeurs = ["P"];
    # shapelets

    methode_correl = "correlation_par_grandeur" #"correlation_par_metrique_fusion";
    metrique_distance = "lb_keogh" #"fdtw_problem"#"lb_keogh" #"sax"#"pearson";
    mode_correlation = "correlationParMorceaux" #"correlationGlissante"#"correlationParMorceaux" # ou liste_mode = ["correlationParMorceaux","correlationGlissante"]
    indice_max_correlation = 0; derivate = False; interval_derivate = 10; 
    epsilon_sax = 0.8; fenetre = 200; # 50
    fenetre_slicing=20; min_len = 140; max_len = 141;#fenetre_slicing=20; min_len = 99; max_len = 100;
    # args
    args = {"methode_correl": methode_correl,"mode_correlation": mode_correlation,\
            "indice_max_correlation":indice_max_correlation,\
            "selected_grandeurs": selected_grandeurs,"modifCorrel": modifCorrel,\
            "equipements": dico["aretes_reseauG"],\
            "fenetre": fenetre,"derivate": derivate, "epsilon_sax": epsilon_sax,\
            "interval_deriv":interval_derivate,"metrique_distance":metrique_distance,\
            "date_debut":date_debut,"date_fin":date_fin,\
            "fenetre_slicing":fenetre_slicing,"min_len":min_len,"max_len":max_len,\
             "rep_root": REP_ROOT}
    metrique_distances = ["lb_keogh","pearson","lcs"]
    metrique_distances = ["distance_pearson"]
    metrique_distances = ["metrique_pearson_damien"]
    metrique_distances = ["metrique_shapelets"]
    metrique_distances = ["distance_pearson"]
    matE = None; 
    for metrique_distance in metrique_distances:
        args["metrique_distance"] = metrique_distance; 
        args["dbg"] = False; args["dbg_0_1"]=False; args["dbg_seuil"]=False;
        args["dbg_ajoutCorrel"] = False; args["dbg_calcul_matE"] = True;
        args["type_cal_matE"] = "max"; 
        matE = correlation_time_series(chemin_datasets, chemin_matrices,\
                              chemin_equipements, args)    
        
        print("ICI--- matE={}".format(matE))
        # distribution des correls {vrai, faux} {positive,negatives} selon une liste de seuil
        matE_reel = pd.read_csv(chemin_matrices+"matE_reelles.csv",index_col = "Unnamed: 0");
        args["matE_proba"] = matE; args["matE_reel"] = matE_reel;
        args["path_save"] = chemin_equipements; 
        courbes.distribution_valeurs_correlation_selon_seuil_methode(args)
        
        # distribution des valeurs de correlation selon les cases (0,1) 
        #        dans le graphe reel(matE_reelles)
#        courbes.distribution_case_0_1_graphe_reelle(args)
        
        # distribution formalisation matE  
#        courbes.distrib_formalisation_matE(args)
#

###### AEFFACER CAR UTILISATION CALCUL SIMILARITE AVEC SHAPELETS        
##    metrique1 = "distance_pearson"; metrique2 = "metrique_wil_old";
##    compare_matE_metrique(args, metrique1, metrique2)
#    
##    for metrique_distance in metrique_distances:
##        args["metrique_distance"] = metrique_distance;
##        args["dbg"] = False; args["dbg_0_1"]=False; args["dbg_seuil"]=False;args["dbg_ajoutCorrel"] = False;
##        args["chemin_equipements"]=chemin_equipements; args["chemin_datasets"]=chemin_datasets;
##        args["chemin_matrices"]=chemin_matrices;
##        args["grandeur"] = "P";
##        df_data = correlation_by_shapelets_grandeur(args)
#    args["dbg_calcul_matE"] = True;
#    args["chemin_matrices"] = chemin_matrices; args["metrique_distance"] = "metrique_shapelets";
#    args["chemin_equipements"] = chemin_equipements; args["grandeur"] = "P"; 
#    args["chemin_datasets"] = chemin_datasets;
##    d = calcul_matE_mean_prod_max(dict(), args) #TODO A DECOMMENTER