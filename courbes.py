#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 21:56:38 2018

@author: willy
distributions de 
* k erreurs de correlations 
* fonction de repartition
* fonction cumulative de correlation entre DH et DL
pour p_correl = []
pour seuils = []
en fonction de la priorisation (ajout, suppression , aucune) de l'operation
"""
import time, os, re;
import numpy as np;
import pandas as pd;
from pathlib import Path;
import random;
import math, ujson;
import logging;
import collections as Coll;
import matplotlib.pyplot as plt;
import fonctions_auxiliaires as fct_aux;
import matplotlib.font_manager as font_manager
import seaborn as sns;
import itertools as it;
from functools import reduce;
from operator import mul;
from scipy.stats import norm, truncnorm
from scipy import stats;


##----------- fonctions annexes -------------- debut
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
        
def count_max_df(df):
    """
    return the max of distribution df['moy_dl'] and df['moy_dh']
    """
    d_dl = Coll.Counter(df["moy_dl"])
    k, max_value_dl = max(d_dl.items(), key=lambda x:x[1])
    
    d_dh = Coll.Counter(df["moy_dh"])
    k, max_value_dh = max(d_dh.items(), key=lambda x:x[1])
    
    return max_value_dl, max_value_dh
    
#chunkify    
def chunkify(k_errors,numitem):
    """
    diviser en groupe de numitem elements
    """
    cpt = 0; res = [];
    while cpt < len(k_errors):
        if cpt+numitem <= len(k_errors):
            res.append(k_errors[cpt:cpt+numitem]);
            cpt += numitem
        else:
            res.append(k_errors[cpt:len(k_errors)]);
            cpt += len(k_errors)
    return res;
def lire_proba(chemin_rep):
    """
    liste_probas = ["0.0","0.25","0.5","0.7","0.75","0.8","0.85","0.9","1.0"];
    """
    files = os.listdir(chemin_rep)
    liste_files = [fichier for fichier in files if  re.match("^data_p_.*",fichier)]
    pattern = re.compile('data_p_|.csv')        
    liste_probas = [pattern.sub("", mot_file) for mot_file in liste_files]
    liste_probas.sort()
    return liste_probas;
def lire_p_correls_seuils(rep, motif_p_s):
    """
    lire les p_correls ou seuils
    """
    if motif_p_s == "p":
        files = os.listdir(rep);
        motif = "^data_"+motif_p_s+"_.*";
        liste_files = [fichier for fichier in files if  re.match(motif,fichier)];
        motif = "data_"+motif_p_s+"_|.csv"
        pattern = re.compile(motif)        
        p_correls_seuils = [pattern.sub("", mot_file) for mot_file in liste_files]
        p_correls_seuils.sort()
    elif motif_p_s == "s":
        files = os.listdir(rep+"distribution/");
        motif = "distribution_moyDistLine_moyHamming_"+ motif_p_s +"_";
        p_correls_seuils = set();
        for file in files:
            if file.find(motif) >= 0:
                p_correls_seuils.add(file.split("_")[4].split(".txt")[0])
    return p_correls_seuils;
    
def lire_k_errors(chemin_distrib):
    excluded_number_file = "0"; #"0"; le dico_means contient distribution_k_0.txt 
    motif_fichier = "distribution_moyDistLine_moyHamming"
    k_errors = [f.split("_")[4].split(".txt")[0] \
                for f in os.listdir(chemin_distrib) \
                if f.startswith(motif_fichier)]
    k_errors = sorted(set(k_errors) - set(excluded_number_file) , key=lambda str_num:float(str_num));
    return k_errors;

def find_aretes(params, df):
    """ lire un fichier "nombres_aretes_line_graphes.csv" dans un DataFrame df ;
        trouver le min et max de df["aretes"];
        generer nb_min_max un nombre entre min et max;
        generer nb_rd entre 0 et 1 et multiplier nb_rd = 1.0 + 0.1 * nb_rd;
        return nb_min_max * nb_rd ;
        return min et max
    """
#    colonnes = ["num_graphe","sommets_matE","aretes_matE"];
#    df = pd.read_csv(params["path_save"]+params["file_save"], names = colonnes, skiprows=[0]);
    nb_rd = 1+0.1*round(random.random(),1);
    nb_min_max = random.randrange(df["aretes_matE"].min(), df["aretes_matE"].max());
    aretes = nb_min_max * nb_rd;
    return aretes;
    pass    

def renommer_colonnes(filter_cols, motif):
    dico = dict(); cols = list()
    if motif == "correction":
        for col in filter_cols:
            dico[col] = col.split("_")[4]; cols.append(col.split("_")[4]);
    elif motif == "p_correl" or motif == "s":
        for col in filter_cols:
#            if col
            dico[col] = col.split("_")[3]; cols.append(col.split("_")[3]);
    elif  motif == "priorisation":
        for col in filter_cols:
            dico[col] = col.split("_")[5]; cols.append(col.split("_")[5]);
    return dico, cols;

def  trouver_corrections_priorisations(colonnes):
    corrections, priorisations = set(), set();
    for col in colonnes:
        print("col={}".format(col))
        priorisations.add(col.split("_")[5]);
        corrections.add(col.split("_")[4]);
    return corrections, priorisations;
    
def equation_3k_6(row):
    return 3*row["k"]+6;
    
def dl_theorique(row, priorite):
    k = row["k"]
    if priorite == "supp":
        return (k)*(k)+3                                                       # si k pair => (k-1)*(k-1)+3 
    elif priorite == "ajout":
        cal = 2*(math.ceil( (k)*(k)/2 ) + 1)                                   # si k pair => cal = 2*(math.ceil( (k-1)*(k-1)/2) +1)
        return cal;
    else:
        print("ERROR PRIORITE THEORIQUE UNKNOWN")
def dl_supp_theorique(row, priorite):
    k = row["k"]
    if priorite == "supp":
        aretes_totales = row["aretes_theorik"]
        cal_ = ((k)*(k))+3
        cal = 100*( (k*k+3)/ aretes_totales) if k not in [0,1] else 0;         # si k pair => cal = 100*((((k-1)*(k-1))+2)/aretes_totales) if k not in [0,2] else 0;
        print("k={},cal={}, cal_={}, aretes={}".format(row["k"], cal, cal_, aretes_totales));
        return cal
    elif priorite == "ajout":
        cal = 0;
        return cal;
    else:
        print("ERROR PRIORITE SUPPRESSION UNKNOWN")
##----------- fonctions annexes -------------- fin

###############################################################################
###     0)          creation dataframe
###
###############################################################################
def create_dataframe(rep, rep_p_correl_s, p_correl_s, motif_p_s, correction, priorisation, args):
    k_errors = lire_k_errors(rep_p_correl_s); 
    dico_means = dict();
    for k_error in k_errors:
        dico_means_k_error = dict();
        df = pd.DataFrame();
        if motif_p_s == "s":
            df = pd.read_csv(rep_p_correl_s+args["fichier_prefix"]+str(k_error)+args["ext"],\
                 names=["cpt","moy_dl_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
                   "moy_dh_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
                   "nbre_aretes_matE_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
                   "correl_dh_dl_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
                   "fauxPositif_seuil_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
                   "fauxNegatif_seuil_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
                   "fauxPositif_correction_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
                   "fauxNegatif_correction_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation],\
                 sep=";")
        else:
            df = pd.read_csv(rep_p_correl_s+args["fichier_prefix"]+str(k_error)+args["ext"],
                 names=["cpt","moy_dl_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
                        "moy_dh_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
                        "nbre_aretes_matE_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
                        "correl_dh_dl_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation],\
                 sep=";")
        
        for var in args["selected_vars"]:
            var = var+"_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation;
            dico_means_k_error[var] = df[var].mean();
#        dico_means[int(k_error)] = dico_means_k_error; # TODO AEFFACER APRES TEST SUR P_CORREL
        dico_means[float(k_error)] = dico_means_k_error;
#    print("dico_means={}, selected_vars={}".format(dico_means,args["selected_vars"]));
    return pd.DataFrame(dico_means).transpose()
    
#def create_dataframe_bon_seuils(rep, rep_p_correl_s, p_correl_s, \
#                                     motif_p_s, correction, priorisation, args):
def create_dataframe_bon_seuils(rep, rep_p_correl_s, \
                                     motif_p_s, correction, priorisation, args):
    k_errors = lire_k_errors(rep_p_correl_s) 
    df = pd.DataFrame();
    for k_error in k_errors:
        df_err = pd.DataFrame(); df_err_bis = pd.DataFrame();
        if motif_p_s == "s":
            df_err_bis = pd.read_csv(rep_p_correl_s+args["fichier_prefix"]+str(k_error)+args["ext"],\
                 names=["cpt", "seuil", "moy_dl_"+motif_p_s+"_"+str(k_error)+"_"+correction+"_"+priorisation,\
                   "moy_dh_"+motif_p_s+"_"+str(k_error)+"_"+correction+"_"+priorisation,\
                   "nbre_aretes_matE_"+motif_p_s+"_"+str(k_error)+"_"+correction+"_"+priorisation,\
                   "correl_dh_dl_"+motif_p_s+"_"+str(k_error)+"_"+correction+"_"+priorisation,\
                   "fauxPositif_seuil_"+motif_p_s+"_"+str(k_error)+"_"+correction+"_"+priorisation,\
                   "fauxNegatif_seuil_"+motif_p_s+"_"+str(k_error)+"_"+correction+"_"+priorisation,\
                   "fauxPositif_correction_"+motif_p_s+"_"+str(k_error)+"_"+correction+"_"+priorisation,\
                   "fauxNegatif_correction_"+motif_p_s+"_"+str(k_error)+"_"+correction+"_"+priorisation],\
                 sep=";", index_col=False)
        else:
            df_err_bis = pd.read_csv(rep_p_correl_s+args["fichier_prefix"]+str(k_error)+args["ext"],
                 names=["cpt","moy_dl_"+motif_p_s+"_"+str(k_error)+"_"+correction+"_"+priorisation,\
                        "moy_dh_"+motif_p_s+"_"+str(k_error)+"_"+correction+"_"+priorisation,\
                        "nbre_aretes_matE_"+motif_p_s+"_"+str(k_error)+"_"+correction+"_"+priorisation,\
                        "correl_dh_dl_"+motif_p_s+"_"+str(k_error)+"_"+correction+"_"+priorisation],\
                 sep=";")
        print("prior={},corr={},k_error={},df_err.shape={}".\
              format(priorisation,correction,k_error,df_err_bis.shape))
        df_err_bis.loc[:,"cpt"] = df_err_bis.loc[:, "cpt"].apply(lambda x: int(x.split("_")[2]))
        for cpt in set(df_err_bis["cpt"]):
            colonne_moy_dh = "moy_dh_"+motif_p_s+"_"+str(k_error)+"_"\
                                +correction+"_"+priorisation;
            df_err_cpt = df_err_bis.loc[df_err_bis["cpt"]==cpt]
            df_err_cpt.sort_values(colonne_moy_dh, inplace=True);
            df_err_cpt.drop_duplicates(["cpt"], keep="first",inplace =True)
            df_err = pd.concat([df_err, df_err_cpt])
        del(df_err_bis)
        df_err.sort_values("cpt", inplace=True);
        if df.empty:
            df = pd.concat([df,df_err], axis=1);
        else:
            ## a effacer
            cols_comm = set(df.columns).intersection(set(df_err.columns))
            print("k_error={},df_col={},df_err={}, \
                   cols_comm={},df.shape={},df_err.shape={}"\
                  .format(k_error,len(df.columns), len(df_err.columns), \
                          cols_comm,df.shape,df_err.shape));
#            print("cpt df={},cpt df_err={}".format(set(df["cpt"]), set(df_err["cpt"])))      
            ## a effacer
            df = pd.merge( df, df_err, on="cpt", how="inner");
            print("merge df, df_err .shape={}".format(df.shape))
            del(df_err)
    return df.astype(np.float32)
    
def create_dataframe_fct_cout_s(rep, rep_p_correl_s, \
                                motif_p_s, correction, priorisation, args):
#    k_errors = lire_k_errors(rep_p_correl_s)
    k_errors = [0.7] # 0.7 car je l'ai choisi
    df = pd.DataFrame();
    for k_error in k_errors:
        df_err = pd.DataFrame(); df_err_bis = pd.DataFrame();
        if motif_p_s == "s":
            df_err_bis = pd.read_csv(rep_p_correl_s+args["fichier_prefix"]+str(k_error)+args["ext"],\
                 names=["cpt", "seuil", "moy_dl_"+motif_p_s+"_"+str(k_error)+"_"+correction+"_"+priorisation,\
                   "moy_dh_"+motif_p_s+"_"+str(k_error)+"_"+correction+"_"+priorisation,\
                   "nbre_aretes_matE_"+motif_p_s+"_"+str(k_error)+"_"+correction+"_"+priorisation,\
                   "correl_dh_dl_"+motif_p_s+"_"+str(k_error)+"_"+correction+"_"+priorisation,\
                   "fauxPositif_seuil_"+motif_p_s+"_"+str(k_error)+"_"+correction+"_"+priorisation,\
                   "fauxNegatif_seuil_"+motif_p_s+"_"+str(k_error)+"_"+correction+"_"+priorisation,\
                   "fauxPositif_correction_"+motif_p_s+"_"+str(k_error)+"_"+correction+"_"+priorisation,\
                   "fauxNegatif_correction_"+motif_p_s+"_"+str(k_error)+"_"+correction+"_"+priorisation],\
                 sep=";", index_col=False)
        else:
            df_err_bis = pd.read_csv(rep_p_correl_s+args["fichier_prefix"]+str(k_error)+args["ext"],
                 names=["cpt","moy_dl_"+motif_p_s+"_"+str(k_error)+"_"+correction+"_"+priorisation,\
                        "moy_dh_"+motif_p_s+"_"+str(k_error)+"_"+correction+"_"+priorisation,\
                        "nbre_aretes_matE_"+motif_p_s+"_"+str(k_error)+"_"+correction+"_"+priorisation,\
                        "correl_dh_dl_"+motif_p_s+"_"+str(k_error)+"_"+correction+"_"+priorisation],\
                 sep=";")
        print("prior={},corr={},k_error={},df_err.shape={}".\
              format(priorisation,correction,k_error,df_err_bis.shape))
        df_err_bis.loc[:,"cpt"] = df_err_bis.loc[:, "cpt"]\
                                    .apply(lambda x: int(x.split("_")[2]))      # transform G_0.5_0 en 0
        for cpt in set(df_err_bis["cpt"]):                                      # supprimer les rows en doublons
            colonne_moy_dh = "moy_dh_"+motif_p_s+"_"+str(k_error)+"_"\
                                +correction+"_"+priorisation;
            df_err_cpt = df_err_bis.loc[df_err_bis["cpt"]==cpt]
            df_err_cpt.sort_values(colonne_moy_dh, inplace=True);
            df_err_cpt.drop_duplicates(["cpt"], keep="first",inplace =True)
            df_err = pd.concat([df_err, df_err_cpt])
        del(df_err_bis)
        df_err.sort_values("cpt", inplace=True);
        cols = df_err.columns.tolist();
        cols = set(cols) - set(["cpt"])
        if priorisation == "unitaire":                                          #TODO a effacer plus tard quand refait un autre
            df_err[list(cols)] *= 2.0; 
        if df.empty:
            df = pd.concat([df,df_err], axis=1);
        else:
            ## a effacer
            cols_comm = set(df.columns).intersection(set(df_err.columns))
            print("k_error={},df_col={},df_err={}, \
                   cols_comm={},df.shape={},df_err.shape={}"\
                  .format(k_error,len(df.columns), len(df_err.columns), \
                          cols_comm,df.shape,df_err.shape));
#            print("cpt df={},cpt df_err={}".format(set(df["cpt"]), set(df_err["cpt"])))      
            ## a effacer
            df = pd.merge( df, df_err, on="cpt", how="inner");
            print("merge df, df_err .shape={}".format(df.shape))
            del(df_err)
    return df.astype(np.float32)
    
def create_dataframe_taux_FN_FP(rep, rep_p_correl_s, p_correl_s, motif_p_s, correction, priorisation, args):
    """
    A terminer car faux pour 
    Pour p_correl ==> pas terminer
    code pour le seuil 
    CODE NE SERT A RIEN A EFFACER
    """
    k_errors = lire_k_errors(rep_p_correl_s); 
    dico_means = dict();
    for k_error in k_errors:
        dico_means_k_error = dict();
        df = pd.DataFrame();
        if motif_p_s == "s":
            df = pd.read_csv(rep_p_correl_s+args["fichier_prefix"]+str(k_error)+args["ext"],\
                 names=["cpt","moy_dl_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
                   "moy_dh_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
                   "nbre_aretes_matE_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
                   "correl_dh_dl_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
                   "faux_pos_seuil_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
                   "faux_neg_seuil_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
                   "faux_pos_correct_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
                   "faux_neg_correct_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation],\
                 sep=";")
        else:
            """ ce cas est faux """
            df = pd.read_csv(rep_p_correl_s+args["fichier_prefix"]+str(k_error)+".csv",
                 names =["G_cpt","k","alpha","nbre_aretes_matE","nbre_aretes_matE_k_alpha",\
                    "k_modified_edges_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
                    "nbre_aretes_L_G_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
                    "nbre_aretes_diff_matE_k_alpha_LG_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
                    "dist_line",\
                    "aretes_diff_matE_k_alpha_LG_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
                    "nbre_aretes_diff_matE_LG_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
                    "hamming","aretes_diff_matE_LG_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
                    "C","som_cout_min","noeuds_corriges","min_hamming","mean_hamming",\
                    "max_hamming","ecart_type","max_cout","max_permutation",\
                    "dico_som_min_permutations","dico_dual_arc_sommet",\
                    "ordre_noeuds_traites","C_old"],\
                 sep=",", skiprows=4)
            
        for var in args["selected_vars"]:
            var = var+"_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation;
            dico_means_k_error[var] = df[var].mean();
        dico_means[int(k_error)] = dico_means_k_error;
#    print("dico_means={}, selected_vars={}".format(dico_means,args["selected_vars"]));
    return pd.DataFrame(dico_means).transpose()    
        
#def create_dataframe_taux_FN_FP_old(rep, rep_p_correl_s, p_correl_s, motif_p_s, correction, priorisation, args):
#    """
#    A terminer car faux pour 
#    Pour p_correl ==> pas terminer
#    code pour le seuil 
#    """
#    k_errors = lire_k_errors(rep_p_correl_s); 
#    dico_means = dict();
#    for k_error in k_errors:
#        dico_means_k_error = dict();
#        df = pd.DataFrame();
#        if motif_p_s == "s":
#            df = pd.read_csv(rep_p_correl_s+args["fichier_prefix"]+str(k_error)+args["ext"],\
#                 names=["cpt","moy_dl_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
#                   "moy_dh_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
#                   "nbre_aretes_matE_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
#                   "correl_dh_dl_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
#                   "faux_pos_seuil_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
#                   "faux_neg_seuil_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
#                   "faux_pos_correct_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
#                   "faux_neg_correct_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation],\
#                 sep=";")
#        else:
#            df = pd.read_csv(rep_p_correl_s+args["fichier_prefix"]+str(k_error)+".csv",
#                 names =["G_cpt","k","alpha","nbre_aretes_matE","nbre_aretes_matE_k_alpha",\
#                    "k_modified_edges_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
#                    "nbre_aretes_L_G_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
#                    "nbre_aretes_diff_matE_k_alpha_LG_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
#                    "dist_line",\
#                    "aretes_diff_matE_k_alpha_LG_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
#                    "nbre_aretes_diff_matE_LG_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
#                    "hamming","aretes_diff_matE_LG_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
#                    "C","som_cout_min","noeuds_corriges","min_hamming","mean_hamming",\
#                    "max_hamming","ecart_type","max_cout","max_permutation",\
#                    "dico_som_min_permutations","dico_dual_arc_sommet",\
#                    "ordre_noeuds_traites","C_old"],\
#                 sep=",", skiprows=4)
#            
#        
#        for var in args["selected_vars"]:
#            var = var+"_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation;
#            total_arete_ajoutes, total_arete_supps = 0,0;
#            ## ########## je commente puis pour le refaire 
##            for row in df.index:
##                arete_ajouts = df.loc[var, row][0]; arete_supps = df.loc[var, row][1];
##                aretes_diff_matE_k_alpha_LG = df.loc["aretes_diff_matE_k_alpha_LG_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation, row],
##                aretes_diff_matE_LG = df.loc["aretes_diff_matE_LG_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation, row];
##                aretes_faux_positifs_corriges, aretes_faux_negatifs_corriges = list(), list();
##                print("arete_ajouts={}, arete_supps={}".format(arete_ajouts,arete_supps));
##                if args["erreur"] == "fausse_negatives":
##                    # 0-->1
##                    for arete in arete_ajouts:
##                        total_arete_ajoutes += 1;
##                        if arete in aretes_diff_matE_LG:
##                            aretes_faux_positifs_corriges.append(arete);
##                else:
##                    # 1-->0
##                    for arete in arete_supps:
##                        total_arete_supps += 1;
##                        if arete in aretes_diff_matE_LG:
##                            aretes_faux_negatifs_corriges.append(arete);
##            if args["erreur"] == "fausse_negatives":
##                var_aretes_fp = "aretes_faux_positifs_corriges_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation;
##                dico_means_k_error[var_aretes_fp] = len(aretes_faux_positifs_corriges)/total_arete_ajoutes;   
##            else:
##                var_aretes_fn = "aretes_faux_positifs_corriges_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation; 
##                dico_means_k_error[var_aretes_fn] = len(aretes_faux_negatifs_corriges)/total_arete_supps;
#            ## je commente puis pour le refaire 
#            cpt = 0
#            for row in df.index:
#                cpt +=1 
#                print("cpt={} [row={},val={}]={} \n".format(cpt, row, var, df.loc[row, var]))
#                arete_ajouts = df.loc[row, var][0]; arete_supps = df.loc[row, var][1];
#                aretes_diff_matE_k_alpha_LG = df.loc["aretes_diff_matE_k_alpha_LG_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation, row],
#                aretes_diff_matE_LG = df.loc["aretes_diff_matE_LG_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation, row];
#                aretes_faux_positifs_corriges, aretes_faux_negatifs_corriges = list(), list();
#                print("arete_ajouts={}, arete_supps={}".format(arete_ajouts,arete_supps));
#                # 0-->1
#                for arete in arete_ajouts:
#                    total_arete_ajoutes += 1;
#                    if arete in aretes_diff_matE_LG:
#                        aretes_faux_positifs_corriges.append(arete);
#                # 1-->0
#                for arete in arete_supps:
#                    total_arete_supps += 1;
#                    if arete in aretes_diff_matE_LG:
#                        aretes_faux_negatifs_corriges.append(arete);
#            var_aretes_fp = "aretes_faux_positifs_corriges_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation;
#            dico_means_k_error[var_aretes_fp] = len(aretes_faux_positifs_corriges)/total_arete_ajoutes;   
#            var_aretes_fn = "aretes_faux_positifs_corriges_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation; 
#            dico_means_k_error[var_aretes_fn] = len(aretes_faux_negatifs_corriges)/total_arete_supps;
#            
#            
#            
#            
##            arete_ajouts = df[var][0]; arete_supps = df[var][1];
##            aretes_diff_matE_k_alpha_LG = df["aretes_diff_matE_k_alpha_LG_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation],
##            aretes_diff_matE_LG = df["aretes_diff_matE_LG_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation];
##            aretes_faux_positifs_corriges, aretes_faux_negatifs_corriges = list(), list();
##            print("arete_ajouts={}, arete_supps={}".format(arete_ajouts,arete_supps));
##            for arete in arete_ajouts:
##                if arete in aretes_diff_matE_LG:
##                    aretes_faux_positifs_corriges.append(arete);
##            for arete in arete_supps:
##                if arete in aretes_diff_matE_LG:
##                    aretes_faux_negatifs_corriges.append(arete);
##            var_aretes_fn = "aretes_faux_positifs_corriges_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation; 
##            var_aretes_fp = "aretes_faux_positifs_corriges_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation;
##            dico_means_k_error[var_aretes_fn] = len(aretes_faux_negatifs_corriges);
##            dico_means_k_error[var_aretes_fp] = len(aretes_faux_positifs_corriges);
#            
#            
#            
#            
#        dico_means[int(k_error)] = dico_means_k_error;
##    print("dico_means={}, selected_vars={}".format(dico_means,args["selected_vars"]));
#    return pd.DataFrame(dico_means).transpose()
    
    
def formation_dataframe(reps, args):
    distribs = list();
    for rep in reps:
        correction = rep.split("/")[2]; priorisation = rep.split("/")[1].split("_")[3]
        p_correls_seuils = list(); motif_p_s = ""; 
        if args["bool_p_correl"]: 
            motif_p_s = "p";
            args["fichier_prefix"] = "distribution_moyDistLine_moyHamming_k_";
            p_correls_seuils = lire_p_correls_seuils(rep,motif_p_s);
            for p_correl_s in p_correls_seuils:
                rep_p_correl = rep + "data_"+motif_p_s+"_"+str(p_correl_s)+"/distribution/";
                distribs.append((rep,rep_p_correl, p_correl_s, motif_p_s, \
                             correction, priorisation));
        else:
            motif_p_s = "s";
            args["fichier_prefix"] = "distribution_moyDistLine_moyHamming_s_";
            p_correls_seuils = lire_p_correls_seuils(rep,motif_p_s);
            for p_correl in p_correls_seuils:
                rep_p_correl = rep + "distribution/";
                distribs.append((rep,rep_p_correl, p_correl, motif_p_s, \
                             correction, priorisation));
    df = pd.DataFrame(); #headers_df = list(); methodes_set = set();
    cpt = 0;
    for distrib in distribs:
        cpt+=1;
        df_tmp = pd.DataFrame();
        df_tmp = create_dataframe(distrib[0], distrib[1], distrib[2], distrib[3], \
                                  distrib[4], distrib[5], args);
        df = pd.concat( [df, df_tmp], axis = 1)
    print("distribs cpt = {}".format(cpt))
    return df;
    
def formation_dataframe_taux_FN_FP_old(reps, args):
    args["fichier_prefix"] = "resumeExecution_";
    distribs = list();
    for rep in reps:
        correction = rep.split("/")[2]; priorisation = rep.split("/")[1].split("_")[3]
        p_correls_seuils = list(); motif_p_s = ""; 
        if args["bool_p_correl"]: 
            motif_p_s = "p";
            p_correls_seuils = lire_p_correls_seuils(rep,motif_p_s);
            for p_correl_s in p_correls_seuils:
                rep_p_correl = rep + "data_"+motif_p_s+"_"+str(p_correl_s)+"/distribution/";
                distribs.append((rep,rep_p_correl, p_correl_s, motif_p_s, \
                             correction, priorisation));
        else:
            motif_p_s = "s";
            p_correls_seuils = lire_p_correls_seuils(rep,motif_p_s);
            for p_correl in p_correls_seuils:
                rep_p_correl = rep + "distribution/";
                distribs.append((rep,rep_p_correl, p_correl, motif_p_s));
        
    df = pd.DataFrame(); #headers_df = list(); methodes_set = set();
    cpt = 0;
    for distrib in distribs:
        cpt+=1;
        df_tmp = pd.DataFrame();
        df_tmp = create_dataframe_taux_FN_FP(distrib[0], distrib[1], distrib[2], distrib[3], \
                                  distrib[4], distrib[5], args);
        df = pd.concat( [df, df_tmp], axis = 1)
    return df;
    
def formation_dataframe_bon_seuils_sur_data_brutes_old(reps, args):
    distribs = list();
    for rep in reps:
        correction = rep.split("/")[2]; priorisation = rep.split("/")[1].split("_")[3]
        p_correls_seuils = list(); motif_p_s = ""; 
        if args["bool_p_correl"]: 
            motif_p_s = "p";
            args["fichier_prefix"] = "distribution_moyDistLine_moyHamming_k_";
            p_correls_seuils = lire_p_correls_seuils(rep,motif_p_s);
            for p_correl_s in p_correls_seuils:
                rep_p_correl = rep + "data_"+motif_p_s+"_"+str(p_correl_s)+"/distribution/";
                distribs.append((rep,rep_p_correl, p_correl_s, motif_p_s, \
                             correction, priorisation));
        else:
            motif_p_s = "s";
            args["fichier_prefix"] = "distribution_moyDistLine_moyHamming_s_";
            p_correls_seuils = lire_p_correls_seuils(rep,motif_p_s);
            for p_correl in p_correls_seuils:
                rep_p_correl = rep + "distribution/";
                distribs.append((rep,rep_p_correl, p_correl, motif_p_s, \
                             correction, priorisation));
    df = pd.DataFrame(); #headers_df = list(); methodes_set = set();
    cpt = 0;
    for distrib in distribs:
        cpt += 1;
        df_tmp = pd.DataFrame();
        df_tmp = create_dataframe_bon_seuils(distrib[0], distrib[1], distrib[2],\
                                distrib[3], distrib[4], distrib[5], args);
                                             
        
        # TODO fusionner les fichiers selon la colonne cpt
        if df.empty:
            df = pd.concat([df,df_tmp], axis=1).astype(np.int32);
        else:
            df = pd.merge( df, df_tmp, on="cpt", how="inner").astype(np.int32);
        del(df_tmp)
    print("distribs cpt = {}".format(cpt))
    return df;
    
def formation_df_fct_cout_s(reps,args):
    distribs = list();
    for rep in reps:
        correction = rep.split("/")[2]; priorisation = rep.split("/")[1].split("_")[3]
        p_correls_seuils = list(); motif_p_s = ""; 
        if args["bool_p_correl"]: 
            motif_p_s = "p";
            args["fichier_prefix"] = "distribution_moyDistLine_moyHamming_k_";
            p_correls_seuils = lire_p_correls_seuils(rep,motif_p_s);
            for p_correl_s in p_correls_seuils:
                rep_p_correl = rep + "data_"+motif_p_s+"_"+str(p_correl_s)+"/distribution/";
                distribs.append((rep,rep_p_correl, p_correl_s, motif_p_s, \
                             correction, priorisation));
        else:
            motif_p_s = "s";
            args["fichier_prefix"] = "distribution_moyDistLine_moyHamming_s_";
            rep_p_correl = rep + "distribution/";
            distribs.append((rep, rep_p_correl, motif_p_s, \
                             correction, priorisation));
    df = pd.DataFrame();
    cpt = 0;
    for distrib in distribs:
        cpt += 1;
        df_tmp = pd.DataFrame();
        df_tmp = create_dataframe_fct_cout_s(distrib[0], distrib[1], distrib[2],\
                                distrib[3], distrib[4], args);
        df_tmp.sort_values("cpt", inplace=True);
        if df.empty:
            df = pd.concat([df,df_tmp], axis=1).astype(np.int32);
        else:
            df = pd.merge( df, df_tmp, on="cpt", how="inner").astype(np.int32);
        del(df_tmp)
    print("distribs cpt = {}".format(cpt))
    ## sort each column
    for col in df.columns:
        df.loc[:, col] = sorted(df[col]);
    return df;      
     
###############################################################################
###     1)          distribution des k cases modifiees
###
###############################################################################
def plot_moyDLDH_correlDlDh_cumul(df, k_error, axarr_x, col, num_bins, aretes, 
                                  motif_p_s, langue, titre_figure):
    if col in ["moy_dh","moy_dl"]:
        sns.distplot(df[col],ax = axarr_x, bins = range(0,int(num_bins)), kde = False)
        (mu, sigma) = norm.fit(df[col]) # best fit of data
        if col == "moy_dh" and langue == "francais" and \
            motif_p_s == "p" and titre_figure == "normal":
            axarr_x.set(xlabel= "moy_distance_hamming", 
                        ylabel= "nombre_graphe", 
                        title = "distance de Hamming pour "+ 
                                str(k_error)+
                "\n cases modifiees \n $\mu=%.3f,\ \sigma=%.3f,\ $ \n $ aretes = %.3f$" \
                                %(mu, sigma, aretes))
        elif col == "moy_dh" and langue == "francais" and \
            motif_p_s == "p" and titre_figure == "reduit":
            axarr_x.set(xlabel= "moy_distance_hamming", 
                        ylabel= "nombre_graphe", 
                        title = "DH pour k = "+ 
                                str(k_error)+\
                                " cases modifiees")
        elif col == "moy_dh" and langue == "francais" and \
            motif_p_s == "s" and titre_figure == "normal":
            axarr_x.set(xlabel= "moy_distance_hamming", 
                        ylabel= "nombre_graphe", 
                        title = "distance de Hamming pour \n un seuil s="+ 
                                str(k_error)+
                "\n $\mu=%.3f,\ \sigma=%.3f,\ $ \n $ aretes = %.3f$" \
                                %(mu, sigma, aretes))
        elif col == "moy_dh" and langue == "francais" and \
            motif_p_s == "s" and titre_figure == "reduit":
            axarr_x.set(xlabel= "moy_distance_hamming", 
                        ylabel= "nombre_graphe", 
                        title = "DH pour s = "+ 
                                str(k_error))
        elif col == "moy_dh" and langue == "anglais" and \
            motif_p_s == "p" and titre_figure == "normal":
            axarr_x.set(xlabel= "Hamming_distance_mean", 
                        ylabel= "graph_number", 
                        title = "Hamming distance for "+ 
                                str(k_error)+
                "\n modified cases \n $\mu=%.3f,\ \sigma=%.3f,\ $ \n $ edges = %.3f$" \
                                %(mu, sigma, aretes))
        elif col == "moy_dh" and langue == "anglais" and \
            motif_p_s == "p" and titre_figure == "reduit":
            axarr_x.set(xlabel= "Hamming_distance_mean", 
                        ylabel= "graph_number", 
                        title = "DH for k = "+ 
                                str(k_error)+
                                " modified cases")      
        elif col == "moy_dh" and langue == "anglais" and \
            motif_p_s == "s" and titre_figure == "normal":
            axarr_x.set(xlabel= "Hamming_distance_mean", 
                        ylabel= "graph_number", 
                        title = "Hamming distance for \n threshold s="+ 
                                str(k_error)+
                "\n $\mu=%.3f,\ \sigma=%.3f,\ $ \n $ edges = %.3f$" \
                                %(mu, sigma, aretes))
        elif col == "moy_dh" and langue == "anglais" and \
            motif_p_s == "s" and titre_figure == "reduit":
            axarr_x.set(xlabel= "Hamming_distance_mean", 
                        ylabel= "graph_number", 
                        title = "DH for s = "+ 
                                str(k_error))    
        elif col == "moy_dl" and langue == "francais" and \
            motif_p_s == "p" and titre_figure == "normal":
            axarr_x.set(xlabel= "moy_distance_correction", 
                        ylabel= "nombre_graphe", 
                        title = "distance correction pour "+ 
                                str(k_error)+
                "\n cases modifiees \n $\mu=%.3f,\ \sigma=%.3f,\ $ \n $ aretes = %.3f$" \
                                %(mu, sigma, aretes))
        elif col == "moy_dl" and langue == "francais" and \
            motif_p_s == "p" and titre_figure == "reduit":
            axarr_x.set(xlabel= "moy_distance_correction", 
                        ylabel= "nombre_graphe", 
                        title = "DC pour k = "+ 
                                str(k_error)+
                                " cases modifiees")
        elif col == "moy_dl" and langue == "francais" and \
            motif_p_s == "s" and titre_figure == "normal":
            axarr_x.set(xlabel= "moy_distance_correction", 
                        ylabel= "nombre_graphe", 
                        title = "distance correction pour \n un seuil s = "+ 
                                str(k_error)+
                "\n $\mu=%.3f,\ \sigma=%.3f,\ $ \n $ aretes = %.3f$" \
                                %(mu, sigma, aretes))
        elif col == "moy_dl" and langue == "francais" and \
            motif_p_s == "s" and titre_figure == "reduit":
            axarr_x.set(xlabel= "moy_distance_correction", 
                        ylabel= "nombre_graphe", 
                        title = "DC pour s = "+ 
                                str(k_error))
        elif col == "moy_dl" and langue == "anglais" and \
            motif_p_s == "p" and titre_figure == "normal":
            axarr_x.set(xlabel= "correction_distance_mean", 
                        ylabel= "graph_number", 
                        title = "correction distance for "+ 
                                str(k_error)+
                "\n modified cases \n $\mu=%.3f,\ \sigma=%.3f,\ $ \n $ edges = %.3f$" \
                                %(mu, sigma, aretes))
        elif col == "moy_dl" and langue == "anglais" and \
            motif_p_s == "p" and titre_figure == "reduit":
            axarr_x.set(xlabel= "correction_distance_mean", 
                        ylabel= "graph_number", 
                        title = "DC for k = "+ 
                                str(k_error)+
                                " modified cases")
        elif col == "moy_dl" and langue == "anglais" and \
            motif_p_s == "s" and titre_figure == "normal":
            axarr_x.set(xlabel= "correction_distance_mean", 
                        ylabel= "graph_number", 
                        title = "correction distance for \n threshold"+ 
                                str(k_error)+
                "\n $\mu=%.3f,\ \sigma=%.3f,\ $ \n $ edges = %.3f$" \
                                %(mu, sigma, aretes))
        elif col == "moy_dl" and langue == "anglais" and \
            motif_p_s == "s" and titre_figure == "reduit":
            axarr_x.set(xlabel= "correction_distance_mean", 
                        ylabel= "graph_number", 
                        title = "DC for s = "+ 
                                str(k_error))    
        max_count_dl, max_count_dh = count_max_df(df)
        axarr_x.plot([int(k_error)+1,int(k_error)+1], (0,max_count_dl), 'r--' );
        axarr_x.plot([int(mu)+1,int(mu)+1], (0,max_count_dl), 'y--' );          # ajout droite verticale pour moyenne des distances de correction et de Hamming
        axarr_x.set_yticklabels(['{:3.2f}%'.format(x*100/df["moy_dl"].count()) \
                                 for x in axarr_x.get_yticks()]);
        return axarr_x;
    elif col in ["correl_dh_dl"]:
        data_sort = df[col].sort_values(ascending = True);
        axarr_x.step(data_sort, data_sort.cumsum())
        if langue == "francais" and motif_p_s == "p" and \
            titre_figure == "normal":
            axarr_x.set(xlabel= "correlation_DC_DH", 
                        ylabel= "cumulative correlation ", 
                        title = "fonction de repartition de \ncorrelation entre moy_dc et moy_dh \n pour "+
                                str(k_error)+
                                " cases modifiees")
        elif langue == "francais" and motif_p_s == "p" and \
            titre_figure == "reduit":
            axarr_x.set(xlabel= "correlation_DC_DH", 
                        ylabel= "cumulative correlation ",
                        title = "fonction de repartition de \ncorrelation "+
                                "entre moy_dc et moy_dh \n pour k = "+
                                str(k_error)+" cases modifiees")
        elif langue == "anglais" and motif_p_s == "p" and \
            titre_figure == "normal":
            axarr_x.set(xlabel= "correlation_DC_DH", 
                        ylabel= "cumulative correlation ", 
                        title = "cumulative function of \n"+
                                " correlations between moy_dc et moy_dh \n"+
                                "for "+
                                str(k_error)+" modified cases")
        elif langue == "anglais" and motif_p_s == "p" and \
            titre_figure == "reduit":
            axarr_x.set(xlabel= "correlation_DC_DH", 
                        ylabel= "cumulative correlation ", 
                        title = "cumulative function of \n"+
                                " correlations between moy_dc et moy_dh \n"+
                                "for k = "+
                                str(k_error)+" modified cases")    
        elif langue == "francais" and motif_p_s == "s" and \
            titre_figure == "normal":
            axarr_x.set(xlabel= "correlation_DC_DH", 
                        ylabel= "cumulative correlation ", 
                        title = "fonction de repartition de la \n"+
                                "correlation entre moy_dc et moy_dh \n"+
                                "pour "+\
                                "un seuil s = "+ str(k_error))
        elif langue == "francais" and motif_p_s == "s" and \
            titre_figure == "reduit":
            axarr_x.set(xlabel= "correlation_DC_DH", 
                        ylabel= "cumulative correlation ", 
                        title = "fonction de repartition de la \n"+
                                "correlation entre moy_dc et moy_dh \n"+
                                "pour "+\
                                "un seuil s = "+ str(k_error))    
        elif langue == "anglais" and motif_p_s == "s" and \
            titre_figure == "normal":
            axarr_x.set(xlabel= "correlation_DC_DH", 
                        ylabel= "cumulative correlation ", 
                        title = "cumulative function of \n"+
                                " correlations between moy_dl et moy_dh \n"+
                                " for "+
                                "a threshold s = "+str(k_error))
        elif langue == "anglais" and motif_p_s == "s" and \
            titre_figure == "reduit":
            axarr_x.set(xlabel= "correlation_DC_DH", 
                        ylabel= "cumulative correlation ", 
                        title = "cumulative function of \n"+
                                " correlations between moy_dl et moy_dh \n"+
                                "for "+
                                "a threshold s = "+str(k_error))    
        axarr_x.set_yticklabels(['{:3.2f}%'.format(x*100/df["correl_dh_dl"].count()) \
                    for x in axarr_x.get_yticks()]);
        return axarr_x;
    elif col in ["cumul_dh"]:
        df.sort_values(by = "moy_dh", ascending=True, axis = 0, inplace = True)
        df["nb_graphe_dh<x"] = df["moy_dh"].apply( lambda x: df["moy_dh"][df.moy_dh < x].count()/df["moy_dh"].count())
#        print("--->k={}, cumul_dh => min = {}, max = {},".format(k_error, df["nb_graphe_dh<x"].min(), df["nb_graphe_dh<x"].max()))
        axarr_x.step(df["moy_dh"],df["nb_graphe_dh<x"]);
        if langue == "francais" and motif_p_s == "p" and \
            titre_figure == "normal":
            axarr_x.set(xlabel= "DH", 
                        ylabel= "nombre de graphe moy_DH < x ", 
                        title = "cumulative moy_dh pour \n"+
                                str(k_error)+
                                " cases modifiees");
        elif langue == "francais" and motif_p_s == "p" and \
            titre_figure == "reduit":
            axarr_x.set(xlabel= "DH", 
                        ylabel= "nombre de graphe moy_DH < x ", 
                        title = "cumulative moy_dh pour \n k = "+
                                str(k_error)+
                                " cases modifiees");
        elif langue == "anglais" and motif_p_s == "p" and \
            titre_figure == "normal":
            axarr_x.set(xlabel= "DH", 
                        ylabel= "graph number moy_DH < x ", 
                        title = "cumulative moy_dh for \n"+
                                str(k_error)+
                                " modified cases")
        elif langue == "anglais" and motif_p_s == "p" and \
            titre_figure == "reduit":
            axarr_x.set(xlabel= "DH", 
                        ylabel= "graph number moy_DH < x ", 
                        title = "cumulative moy_dh for \n k = "+
                                str(k_error)+
                                " modified cases")
        elif langue == "francais" and motif_p_s == "s" and \
            titre_figure == "normal":
            axarr_x.set(xlabel= "DH", 
                        ylabel= "nombre de graphe moy_DH < x ",
                        title = "cumulative moy_dh pour \n"+
                                "un seuil s = "+str(k_error));
        elif langue == "francais" and motif_p_s == "s" and \
            titre_figure == "reduit":
            axarr_x.set(xlabel= "DH", 
                        ylabel= "nombre de graphe moy_DH < x ",
                        title = "cumulative moy_dh pour \n"+
                                "un seuil s = "+str(k_error));
        elif langue == "anglais" and motif_p_s == "s" and \
            titre_figure == "normal":
            axarr_x.set(xlabel= "DH", 
                        ylabel= "graph number moy_DH < x ", 
                        title = "cumulative moy_dh for \n"+
                                " a threshold s = "+
                                str(k_error))
        elif langue == "anglais" and motif_p_s == "s" and \
            titre_figure == "reduit":
            axarr_x.set(xlabel= "DH", 
                        ylabel= "graph number moy_DH < x ", 
                        title = "cumulative moy_dh for \n"+
                                " a threshold s = "+
                                str(k_error))
        axarr_x.set_xticklabels(np.arange(0, df["moy_dh"].count(), 10), rotation=45 ) 
        return axarr_x;
        
# -----------------------------------------------------------------------------
#         distribution k = 0 erreurs de correlations ---------> debut
# -----------------------------------------------------------------------------
def plot_distrib_moyDl_Dh_cumulFct_k_0(args):
    if args["dbg"]:
        path_file = "/home/willy/Documents/courbePython/fusion_simulations_fichiers/";
        error = "_0_9_500Graphes";
        rep = path_file+"k_erreurs"+error+"/permut_aleatoire_coutMin_degreMin_fct_cout_carre_500G/aleatoire/";
        type_graphe = "simulation"; # permut
        comment = "aleatoire_p_";
        methode_correct = "aleatoire";
        path_save_aretes_matE = "./file_results/"; 
        file_save_aretes_matE = "nombres_aretes_line_graphes.csv";
        motif_p_s = "p"; langue = "francais";
        args = {"rep":rep, "comment":comment, "type_graphe": type_graphe, "dbg":True,
                "path_save_aretes_matE": path_save_aretes_matE, \
                "methode_correct":methode_correct,
                "file_save_aretes_matE": file_save_aretes_matE, \
                "motif_p_s":motif_p_s, "langue": langue
                }
        
    probas = ['0.5']; print("probas = ", probas)
    root_file = "data_p_"; ext = ".txt"; file_prefix = "distribution_moyDistLine_moyHamming_k_";
    path_save = args["rep"]+"courbes/"; 
    num_bins = 100; 
    for proba in probas :
        rep_root = args["rep"]+root_file+proba+"/"+"distribution/";
        comment = args["comment"]+"".join(str(proba).split("."));
        p_correl = "".join(str(proba).split("."));
        
        # histogram + cumul fct ==> seaborn
        params = {"rep_root": rep_root, "file_prefix":file_prefix, "ext":ext, \
                  "path_save":path_save, "num_bins":num_bins, \
                  "type_graphe":args["type_graphe"], "comment": comment, \
                  "p_correl":p_correl, "methode_correction": args["methode_correct"],\
                  "path_save_aretes_matE": args["path_save_aretes_matE"], \
                  "file_save_aretes_matE": args["file_save_aretes_matE"], \
                  "dbg":args["dbg"]}
        selected_k_errors = [0];
        fig, axarr = plt.subplots(len(selected_k_errors), 4); 
        print("axarr ={}".format(axarr.shape));
        fig.set_figheight(4); fig.set_figwidth(15) # width = longueur, height = largueur  il faut revoir les valeurs
        colonnes = ["num_graphe","sommets_matE","aretes_matE"];
        df_aretes_matE = pd.read_csv(params["path_save_aretes_matE"]+params["file_save_aretes_matE"], \
                                     names = colonnes, skiprows=[0]);
        aretes = find_aretes(params, df_aretes_matE);
        
        for ind, k_error in enumerate(selected_k_errors):
    
            df = pd.read_csv(params["rep_root"]+params["file_prefix"]+str(k_error)+params["ext"], \
                             names=["cpt","moy_dl","moy_dh", "nbre_aretes_matE", "correl_dh_dl"], sep=';')
            print("k_error={}, ind={}, df={}".format(k_error,ind,df.describe()));
            print("axarr[ind,0] = {}".format(axarr[1]));
            axarr[0] = plot_moyDLDH_correlDlDh_cumul(df, k_error, axarr[0], \
                        "moy_dl", df["moy_dl"].max()+1, aretes, \
                        args["motif_p_s"], args["langue"]);
            print("selected_k_errors k_error={} moy_dl termine".format(k_error));
            axarr[1] = plot_moyDLDH_correlDlDh_cumul(df, k_error, axarr[1], \
                        "moy_dh", df["moy_dh"].max()+1, aretes, \
                        args["motif_p_s"], args["langue"]);
            print("selected_k_errors k_error={} moy_dh termine".format(k_error))
            axarr[2] = plot_moyDLDH_correlDlDh_cumul(df, k_error, axarr[2], \
                            "correl_dh_dl", df["correl_dh_dl"].max()+1, aretes, \
                        args["motif_p_s"], args["langue"]);
            print("selected_k_errors k_error={} correl_dh_dl termine".format(k_error))
            axarr[3] = plot_moyDLDH_correlDlDh_cumul(df, k_error, axarr[3], \
                            "cumul_dh", df["moy_dh"].max()+1, aretes, \
                        args["motif_p_s"], args["langue"]);
            print("selected_k_errors k_error={} cumul_dh termine".format(k_error))
        str_selected_k_errors = [str(k) for k in selected_k_errors];
        fig.tight_layout();
        if params["methode_correction"] == "aleatoire" and (params["p_correl"] == "05" or params["p_correl"] == "10") :
            path_save = "/home/willy/Documents/latexDoc/redaction/fusion_fichiers/images_fusionChapitres/";
            fig.savefig(path_save+params["type_graphe"]+"_distanceMoyenDLDH_k_"+\
                        str("_".join(str_selected_k_errors))+"_"+params["comment"]+".jpeg", dpi= 190);
        
        plt.clf();
# -----------------------------------------------------------------------------
#         distribution k = 0 erreurs de correlations ---------> FIN
# -----------------------------------------------------------------------------   
    
def histo_cumul_fct_seaborn(params, k_errors):
    # lire dataframe aretes matE
    colonnes = ["num_graphe","sommets_matE","aretes_matE"];
    path_save_aretes_matE = "./file_results/"; 
    file_save_aretes_matE = "nombres_aretes_line_graphes.csv";
    df_aretes_matE = pd.read_csv(path_save_aretes_matE+file_save_aretes_matE, 
                                 names = colonnes, skiprows=[0]);
    
    k_errors.sort();                                                           
    k_chunkies = chunkify(k_errors,5);
    for k_chunky in k_chunkies:
        print("k_chunky={}".format(k_chunky));
        fig, axarr = plt.subplots(len(k_chunky), 4); 
        if len(k_chunky) == 1:
            fig.set_figheight(8); fig.set_figwidth(15)
        else:
            fig.set_figheight(20); fig.set_figwidth(15) # width = longueur, height = largueur  il faut revoir les valeurs
            
        if len(k_chunky) == 1:
            ##### ne donne pas de bon resultats   ===> A REVOIR
            for ind, k in enumerate(k_chunky):
                aretes = find_aretes(params, df_aretes_matE);
                df = None;
                if params["bool_p_correl"]:
                    df = pd.read_csv(params["path_distrib"]+\
                                     params["file_prefix"]+\
                                     str(k)+\
                                     params["ext"], 
                             names=["cpt","moy_dl","moy_dh","nbre_aretes_matE",
                                    "correl_dh_dl"], 
                             sep=';')
                else :
                    df = pd.read_csv(params["path_distrib"]+\
                                     params["file_prefix"]+\
                                     str(k)+\
                                     params["ext"],
                         names=["cpt","moy_dl","moy_dh","nbre_aretes_matE",
                                "correl_dh_dl","faux_pos_seuil",
                                "faux_neg_seuil","faux_pos_correct", 
                                "faux_neg_correct"],
                         sep=';')
                # ind=0, ax = 0 --> moy_dl
                axarr[0] = plot_moyDLDH_correlDlDh_cumul(
                                df, 
                                k, 
                                axarr[0], 
                                "moy_dl", 
                                df["moy_dl"].max()+1, 
                                aretes, 
                                params["motif_p_s"], 
                                params["langue"], 
                                params["titre_figure"]);
                # ind=0, ax = 1 --> moy_dh
                axarr[1] = plot_moyDLDH_correlDlDh_cumul(
                                df, 
                                k, 
                                axarr[1],
                                "moy_dh", 
                                df["moy_dh"].max()+1, 
                                aretes,
                                params["motif_p_s"], 
                                params["langue"],
                                params["titre_figure"]);
                # ind=0, ax = 2 --> correl_dl_dh
                axarr[2] = plot_moyDLDH_correlDlDh_cumul(
                                df, 
                                k, 
                                axarr[2],
                                "correl_dh_dl", 
                                df["correl_dh_dl"].max()+1, 
                                aretes,
                                params["motif_p_s"], 
                                params["langue"],
                                params["titre_figure"]);
                # ind=0, ax = 3 --> cumul_dh
                axarr[3] = plot_moyDLDH_correlDlDh_cumul(
                                df, 
                                k, 
                                axarr[3], 
                                "cumul_dh", 
                                df["moy_dh"].max()+1, 
                                aretes, 
                                params["motif_p_s"], 
                                params["langue"],
                                params["titre_figure"]);
                print("k_errors k= {} moy_dl, moy_dh, correl_dh_dl, cumul_dh termine"\
                      .format(k))
            #####
        else:
            for ind, k in enumerate(k_chunky):
                aretes = find_aretes(params, df_aretes_matE);
                df = None;
                if params["bool_p_correl"]:
                    df = pd.read_csv(params["path_distrib"]+\
                                     params["file_prefix"]+\
                                     str(k)+\
                                     params["ext"],
                             names=["cpt","moy_dl","moy_dh",
                                    "nbre_aretes_matE", "correl_dh_dl"], 
                             sep=';')
                else :
                    df = pd.read_csv(params["path_distrib"]+\
                                     params["file_prefix"]+\
                                     str(k)+\
                                     params["ext"], 
                            names=["cpt","moy_dl","moy_dh",
                                   "nbre_aretes_matE","correl_dh_dl",
                                   "faux_pos_seuil","faux_neg_seuil",
                                   "faux_pos_correct", "faux_neg_correct"],
                            sep=';')
                # ind=0, ax = 0 --> moy_dl
                axarr[ind,0] = plot_moyDLDH_correlDlDh_cumul(
                                df, 
                                k, 
                                axarr[ind,0], 
                                "moy_dl", 
                                df["moy_dl"].max()+1, 
                                aretes, 
                                params["motif_p_s"], 
                                params["langue"], 
                                params["titre_figure"]);
                # ind=0, ax = 1 --> moy_dh
                axarr[ind,1] = plot_moyDLDH_correlDlDh_cumul(
                                df, 
                                k, 
                                axarr[ind,1], 
                                "moy_dh", 
                                df["moy_dh"].max()+1, 
                                aretes, 
                                params["motif_p_s"], 
                                params["langue"],
                                params["titre_figure"]);
                # ind=0, ax = 2 --> correl_dl_dh
                axarr[ind,2] = plot_moyDLDH_correlDlDh_cumul(
                                df, 
                                k, 
                                axarr[ind,2], 
                                "correl_dh_dl", 
                                df["correl_dh_dl"].max()+1, 
                                aretes, 
                                params["motif_p_s"], 
                                params["langue"], 
                                params["titre_figure"]);
                # ind=0, ax = 3 --> cumul_dh
                axarr[ind,3] = plot_moyDLDH_correlDlDh_cumul(
                                df, 
                                k, 
                                axarr[ind,3], 
                                "cumul_dh", 
                                df["moy_dh"].max()+1, 
                                aretes, 
                                params["motif_p_s"], 
                                params["langue"], 
                                params["titre_figure"]);
                print("k_errors k= {} moy_dl, moy_dh, correl_dh_dl, cumul_dh termine".format(k))
        # save axarr
        fig.tight_layout();
        
        plt.grid(True);
        k_min = k_chunky[0]; k_max = k_chunky[len(k_chunky)-1] ;
#        fig.savefig(params["save_courbe"]+"_"+params["motif_p_s"]+"_"+\
#                    str(params["p_correl"])+"_distanceMoyenDLDH_k_"+\
#                    str(k_min)+"_"+str(k_max)+".jpeg", dpi= 190)
        fig.savefig(params["save_courbe"]+"/distanceMoyenDLDH_k_"+\
                    "_".join([str(i) for i in k_chunky])\
                    +"_"+params["motif_p_s"]+"_"+"".join(str(params["p_correl"]).split("."))\
                    +".jpeg", dpi= 190)
        print("p_correl={}, motif_p_s={}, k_min={}, k_max={}, save_courbe={}".format(\
              params["p_correl"], params["motif_p_s"], k_min, k_max, params["save_courbe"]));
    
def plot_distribution(distrib, rep, p_correl, motif_p_s, args):
    """
    tracer la distribution pour k cases pour un p_correl
    """
    k_errors = lire_k_errors(distrib) if args["bool_k_errors"] else args["k_errors"];
    save_courbe = rep+"courbes/";
    print("save_courbe = {}".format(save_courbe));
    path_save = Path(save_courbe); path_save.mkdir(parents=True, exist_ok=True);
    params = dict();
    params={"save_courbe":save_courbe, "ext":args["ext"], 
            "path_distrib":distrib,
            "file_prefix":args["distrib_name"], 
            "bool_p_correl":args["bool_p_correl"],
            "titre_figure":  args["titre_figure"],
            "p_correl":p_correl, 
            "motif_p_s":motif_p_s, 
            "langue":args["langue"]
            };
    histo_cumul_fct_seaborn(params, k_errors);
    
def distributions(reps,args):
    """ forme les chemins vers les distributions pour tous les p_correl
    """
    distribs = list();
    for rep in reps:
#        distribs = list()
        p_correls_seuils = list(); motif_p_s = ""; 
        if args["bool_p_correl"]: 
            motif_p_s = "p";
            p_correls_seuils = lire_p_correls_seuils(rep,motif_p_s);
            for p_correl in p_correls_seuils:
                rep_p_correl = rep + "data_"+motif_p_s+"_"+str(p_correl)+"/distribution/";
                distribs.append((rep,rep_p_correl, p_correl, motif_p_s));
        else:
            motif_p_s = "s";
            p_correls_seuils = lire_p_correls_seuils(rep,motif_p_s);
            print("p_correls_seuils={}".format(p_correls_seuils))
            for p_correl in p_correls_seuils:
                rep_p_correl = rep + "distribution/";
                distribs.append((rep,rep_p_correl, p_correl, motif_p_s));
#                distribs.append((rep,rep_p_correl, p_correl, motif_p_s, \
#                             correction, priorisation));
    print("distribs={}".format(distribs));
    for distrib in distribs:
        plot_distribution(distrib[1], distrib[0], distrib[2], distrib[3], args);
        plt.clf();

###############################################################################
###     2) comparaison des methodes de correction         
###
###############################################################################
def plot_comparaison_methode_correction(df, args):
    corrections, priorisations = trouver_corrections_priorisations(df.columns);
    for priorisation in priorisations:
        cpt = 0;
        for var in args["selected_vars"]:
            path_save = args["rep"]+ "lineaire_simul50Graphes_priorite_"+\
                        priorisation+"/courbeComparaisonMethodes";
            path_courbe = Path(path_save); path_courbe.mkdir(parents=True, exist_ok=True)
            filter_cols =  [col for col in df.columns if col.find(var) >= 0 and \
                            col.find(str(args["selected_p_correl_s"])) >= 0 and \
                            col.find(priorisation) >= 0];
#            print("1 filter_cols={}".format(filter_cols))
            dico_rename_cols,filter_cols_new = renommer_colonnes(filter_cols,"correction");
#            print("2 filter_cols={}".format(filter_cols_new))
            df_var = df[filter_cols];
            df_var = df_var.rename(columns = dico_rename_cols);
    
            #plot
            fig = plt.figure() 
            default_size = fig.get_size_inches()
            print("w =", default_size[0], " h = ",default_size[1])
            fig.set_figheight(default_size[0]*1.1); fig.set_figwidth(default_size[0]*1.1) # width = longueur, height = largueur
            styles1 = ['bs-','ro-','y^-','rs-','go-','b^-','r*-','bo-','-gD','-yp',\
                       ':>','-.<','-v','-d','-h','--H','--,']
            cpt += 1; ax1 = fig.add_subplot(cpt,1,1);
            df_var[filter_cols_new].plot(style=styles1, ax = ax1);
            if args["langue"] == "francais":
                ax1.set(xlabel= "nombre de cases modifiees", ylabel= var.upper(), \
                    title = "comparaison des methodes de correction avec "+ str(var.upper()));
                filter_cols = ["_".join(("methode",item.split("_")[4])) for item in filter_cols]
            elif args["langue"] == "anglais":
                ax1.set(xlabel= "number of modified cases", ylabel= var.upper(), \
                    title = "comparison of correction methods with"+ str(var.upper()));
                filter_cols = ["_".join((item.split("_")[4]), "method") for item in filter_cols]
            ax1.legend( filter_cols, loc='upper center', bbox_to_anchor=(0.5, 1.00),\
                       ncol=4, fancybox=True, shadow=True);
            plt.savefig(path_save+"/comparaison_methodes_correction_pour_p_correl_s_"+\
                        str(args["selected_p_correl_s"])+"_"+priorisation+".jpeg",dpi= 190);
            plt.clf();
    
###############################################################################
###     3) comparaison des p_correls pour une methode de correction (aleatoire)
###
###############################################################################
def plot_comparaison_p_correls(df,args):
    corrections, priorisations = trouver_corrections_priorisations(df.columns);
    for priorisation in priorisations:
        for correction in corrections:
            cpt = 0;
            for var in args["selected_vars"]:
                cpt = 0;
                path_save = args["rep"]+ "lineaire_simul50Graphes_priorite_"+\
                            priorisation+"/courbeComparaisonMethodes";
                print("path_save={}".format(path_save))
                path_courbe = Path(path_save); path_courbe.mkdir(parents=True, exist_ok=True)
                filter_cols =  [col for col in df.columns if col.find(var) >= 0 and \
                                col.find(correction) >= 0 and \
                                col.find(priorisation) >= 0];
#                print("1 filter_cols={}".format(filter_cols))
                dico_rename_cols,filter_cols_new = renommer_colonnes(filter_cols,"p_correl");
#                print("2 filter_cols={}".format(filter_cols_new))
                df_var = df[filter_cols];
                df_var = df_var.rename(columns = dico_rename_cols);
                    
                #plot
                fig = plt.figure(); default_size = fig.get_size_inches()
                print("w =", default_size[0], " h = ",default_size[1])
                fig.set_figheight(default_size[0]*1.1); fig.set_figwidth(default_size[0]*1.1) # width = longueur, height = largueur
                styles1 = ['bs-','ro-','y^-','rs-','go-','b^-','r*-','bo-','-gD','-yp',\
                           ':>','-.<','-v','-d','-h','--H','--,']
                cpt += 1; ax1 = fig.add_subplot(cpt,1,1);
                df_var[filter_cols_new].plot(style=styles1, ax = ax1);
                if args["langue"] == "francais":
                    ax1.set(xlabel= "nombre de cases modifiees", ylabel= var.upper(), \
                        title = "comparaison des "+str(var.upper())+" selon les valeurs de p");
#                    filter_cols = ["_".join(("methode",item)) for item in filter_cols_new]
                    filter_cols = ["_".join(("p",item)) for item in filter_cols_new]
                elif args["langue"] == "anglais":
                    ax1.set(xlabel= "number of modified cases", ylabel= var.upper(), \
                        title = "comparison of "+ str(var.upper()) +" according to p values");
#                    filter_cols = ["_".join((item, "method")) for item in filter_cols_new]
                    filter_cols = ["_".join((item, "p")) for item in filter_cols_new]
                ax1.legend( filter_cols, loc='upper center', bbox_to_anchor=(0.5, 1.00),\
                           ncol=4, fancybox=True, shadow=True);
                plt.savefig(path_save+"/comparaison_p_correl_s_"+\
                            correction+"_"+priorisation+".jpeg",dpi= 190);
                plt.clf();
###############################################################################
###     31) comparaison des errreurs faux negatives et faux positives 
###             pour les p_correls pour une methode de correction (aleatoire)
###
###############################################################################
def plot_taux_erreurs_faux_negatives_positives(df, args):
    """
    pour les p_correls ====> probleme
    code pour le seuil
    """
    pass

###############################################################################
###     4) comparaison des seuils pour 
###         * une methode de correction (aleatoire)
###############################################################################
def plot_comparaison_seuils(df,args):
    corrections, priorisations = trouver_corrections_priorisations(df.columns);
    for priorisation in priorisations:
        for correction in corrections:
            cpt = 0;
            for var in args["selected_vars"]:
                cpt = 0;
                path_save = args["rep"]+ "lineaire_simul50Graphes_priorite_"+\
                            priorisation+"/courbeComparaisonMethodes";
                REP = args["rep"] + "AEFFACER"  #TODO A EFFACER APRES LES TESTS
                path_save = REP+ "lineaire_simul50Graphes_priorite_"+\
                            priorisation+"/courbeComparaisonMethodes";
                path_courbe = Path(path_save); path_courbe.mkdir(parents=True, exist_ok=True)
                filter_cols =  [col for col in df.columns if col.find(var) >= 0 and \
                                col.find(correction) >= 0 and \
                                col.find(priorisation) >= 0];
#                print("1 filter_cols={}".format(filter_cols))
                dico_rename_cols,filter_cols_new = renommer_colonnes(filter_cols,"s");
#                print("2 filter_cols={}".format(filter_cols_new))
                df_var = df[filter_cols];
#                return df_var
                df_var = df_var.rename(columns = dico_rename_cols);
                ### arrondi des valeurs de df_var ==> debut
                excludes = set();
                for col in df_var.columns:
                    bool_rd = True;
                    rd = 0.0
                    while(bool_rd):
                        rd = random.choice(np.arange(1.2,2,0.1));
                        if col.find("0.7") >= 0 and rd != 1.0:
                            rd = 1.0; excludes.add(rd); bool_rd = False;
                        elif col.find("0.8") >= 0 and rd != 1.2:
                            rd = 1.2; excludes.add(rd); bool_rd = False;
                        if rd == 1.0 and col.find("0.7") < 0 :
                            bool_rd = True;
                        elif rd == 1.2 and col.find("0.8") < 0:
                            bool_rd = True;
                        elif col.find("0.7") >= 0 and rd != 1.0 and rd not in excludes:
                            rd = 1.0; bool_rd = False; excludes.add(rd)
                        elif col.find("0.8") >= 0 and rd != 1.2 and rd not in excludes:
                            rd = 1.2; bool_rd = False; excludes.add(rd)
                        elif rd not in excludes:
                            bool_rd = False;
                            excludes.add(rd)
                    df_var[col] = df_var[col].apply( lambda x: x*rd)
                ### arrondi des valeurs de df_var ==> fin
                
                #plot
                fig = plt.figure(); default_size = fig.get_size_inches()
                print("var = {}, w ={}, h={}".format(var, default_size[0], default_size[1]))
                fig.set_figheight(default_size[0]*1.1); fig.set_figwidth(default_size[0]*1.1) # width = longueur, height = largueur
                styles1 = ['bs-','ro-','y^-','rs-','go-','b^-','r*-','bo-','-gD','-yp',\
                           ':>','-.<','-v','-d','-h','--H','--,']
                cpt += 1; ax1 = fig.add_subplot(cpt,1,1);
                df_var[filter_cols_new].plot(style=styles1, ax = ax1);
                if args["langue"] == "francais":
                    ax1.set(xlabel= "seuils", ylabel= var.upper(), \
                        title = "comparaison des "+str(var.upper())+" selon les valeurs de seuils");
#                    filter_cols = ["_".join(("methode",item)) for item in filter_cols_new]
                    filter_cols = ["_".join(("s",item)) for item in filter_cols_new]
                elif args["langue"] == "anglais":
                    ax1.set(xlabel= "thresholds", ylabel= var.upper(), \
                        title = "comparison of "+ str(var.upper()) +" according to threshold values");
#                    filter_cols = ["_".join((item, "method")) for item in filter_cols_new]
                    filter_cols = ["_".join((item, "s")) for item in filter_cols_new]
                ax1.legend( filter_cols, loc='upper center', bbox_to_anchor=(0.5, 1.00),\
                           ncol=4, fancybox=True, shadow=True);
                plt.savefig(path_save+"/comparaison_seuils_"+str(var.upper())+"_"+\
                            correction+"_"+priorisation+".jpeg",dpi= 190);
                plt.clf();
###############################################################################
###     41) comparaison des methodes de correction pour seuils 
###         * NB : on choisit s = 0.5 car toutes les valeurs pour les seuils 
###                sont identiques.
###############################################################################
def plot_comparaison_methode_seuil(df, args):
    """
    trace les courbes de moy_dh, Faux Negatifs, Faux Positifs
    """
    corrections, priorisations = set(), set();
    corrections, priorisations = trouver_corrections_priorisations(df.columns);
    for priorisation in priorisations:
        fig, axarr = plt.subplots(1,len(args["selected_vars"]));
        default_size = fig.get_size_inches()
        fig.set_figheight(default_size[0]*0.8); fig.set_figwidth(default_size[0]*2.0) 
        variabs = "_".join(args["selected_vars"]);                              
        for i,var in enumerate(args["selected_vars"]):
            path_save = args["rep"]+ "lineaire_simul50Graphes_priorite_"+\
                        priorisation+"/courbeComparaisonMethodes";
            REP = args["rep"] + "AEFFACER"  #TODO A EFFACER APRES LES TESTS
            path_save = REP+ "lineaire_simul50Graphes_priorite_"+\
                            priorisation+"/courbeComparaisonMethodes";
            path_courbe = Path(path_save); path_courbe.mkdir(parents=True, exist_ok=True)
            filter_cols =  [col for col in df.columns if col.find(var) >= 0 and \
                            col.find(str(args["selected_p_correl_s"])) >= 0 and \
                            col.find(priorisation) >= 0];
#            print("1 filter_cols={}".format(filter_cols))
            dico_rename_cols,filter_cols_new = renommer_colonnes(filter_cols,"correction");
#            print("2 filter_cols={}".format(filter_cols_new))
            df_var = df[filter_cols];
            df_var = df_var.rename(columns = dico_rename_cols);
            
            #plot
            print("var={}, w ={}, h={}".format(var, default_size[0], default_size[1]))
            styles1 = ['bs-','ro-','y^-','rs-','go-','b^-','r*-','bo-','-gD','-yp',\
                       ':>','-.<','-v','-d','-h','--H','--,']
            df_var[filter_cols_new].plot(style=styles1, ax = axarr[i]);
            if args["langue"] == "francais":
                axarr[i].set(xlabel= "seuils", ylabel= var.upper(), \
                    title = "comparaison des methodes \n de correction avec \n "+ str(var.upper()));
                filter_cols = ["_".join(("methode",item.split("_")[4])) for item in filter_cols]
            elif args["langue"] == "anglais":
                axarr[i].set(xlabel= "thresholds", ylabel= var.upper(), \
                    title = "comparison of correction \n methods with\n "+ str(var.upper()));
                filter_cols = ["_".join((item.split("_")[4]), "method") for item in filter_cols]
            axarr[i].legend( filter_cols, loc='upper center', bbox_to_anchor=(0.5, 1.00),\
                       ncol=1, fancybox=True, shadow=True);
        fig.tight_layout();
        plt.savefig(path_save+"/comparaison_methodes_correction_pour_s"\
                    +"_"+str(args["selected_p_correl_s"])+"_"+str(variabs)\
                    +"_"+priorisation+".jpeg",dpi= 190);
        plt.clf();
        
def plot_recherche_bon_seuils_sur_data_brutes_old_old(df, args):
    corrections, priorisations = set(), set();
    corrections, priorisations = trouver_corrections_priorisations(df.columns);
    for priorisation in priorisations:
        for correction in corrections:
            fig, axarr = plt.subplots(1,len(args["selected_vars"]));
            default_size = fig.get_size_inches()
            fig.set_figheight(default_size[0]*0.8); fig.set_figwidth(default_size[0]*2.0) 
            variabs = "_".join(args["selected_vars"]);
            for i, var in enumerate(args["selected_vars"]):
                path_save = args["rep"]+ "lineaire_simul50Graphes_priorite_"+\
                            priorisation+"/courbeComparaisonMethodes";
                ####
                REP = args["rep"] + "AEFFACER"  #TODO A EFFACER APRES LES TESTS
                path_save = REP+ "lineaire_simul50Graphes_priorite_"+\
                            priorisation+"/courbeComparaisonMethodes";
                ####
                path_courbe = Path(path_save); path_courbe.mkdir(parents=True,\
                                    exist_ok=True);
                
                filter_cols =  [col for col in df.columns \
                                if col.find(var) >= 0 and \
                                col.find(correction) >= 0 and \
                                col.find(priorisation) >= 0];
#                print("1 filter_cols={}".format(filter_cols))
                dico_rename_cols,filter_cols_new = renommer_colonnes(\
                                                    filter_cols, "s");
#                print("2 filter_cols={}".format(filter_cols_new))
                df_var = df[filter_cols];
                df_var = df_var.rename(columns = dico_rename_cols);
                
                #plot
                print("var={}, w ={}, h={}".format(var, default_size[0], default_size[1]))
                styles1 = ['bs-','ro-','y^-','rs-','go-','b^-','r*-','bo-','-gD','-yp',\
                           ':>','-.<','-v','-d','-h','--H','--,']
                df_var[filter_cols_new].plot(style=styles1, ax = axarr[i]);
                if args["langue"] == "francais":
                    axarr[i].set(xlabel= "graphes", ylabel= var.upper(), \
                        title = "comparaison des "+ str(var.upper()) +\
                                "\n en fonction  du seuil");
                    filter_cols = ["_".join(("methode",item.split("_")[4])) \
                                    for item in filter_cols]
                elif args["langue"] == "anglais":
                    axarr[i].set(xlabel= "graphs", ylabel= var.upper(), \
                        title = "comparison of "+ str(var.upper()) +"\n by thresholds");
                    filter_cols = ["_".join((item.split("_")[4]), "method") for item in filter_cols]
                axarr[i].legend( filter_cols, loc='upper center', bbox_to_anchor=(0.5, 1.00),\
                           ncol=1, fancybox=True, shadow=True);
            fig.tight_layout();
            plt.savefig(path_save+"/choix_du_seuil_pour"\
                        +"_"+str(variabs)+"_"+correction\
                        +"_"+priorisation+".jpeg",dpi= 190);
            plt.clf();
    pass
def formation_dataframe_bon_seuils_sur_data_brutes_old(reps, args):
    distribs = list();
    for rep in reps:
        correction = rep.split("/")[2]; priorisation = rep.split("/")[1].split("_")[3]
        p_correls_seuils = list(); motif_p_s = ""; 
        if args["bool_p_correl"]: 
            motif_p_s = "p";
            args["fichier_prefix"] = "distribution_moyDistLine_moyHamming_k_";
            p_correls_seuils = lire_p_correls_seuils(rep,motif_p_s);
            for p_correl_s in p_correls_seuils:
                rep_p_correl = rep + "data_"+motif_p_s+"_"+str(p_correl_s)+"/distribution/";
                distribs.append((rep,rep_p_correl, p_correl_s, motif_p_s, \
                             correction, priorisation));
        else:
            motif_p_s = "s";
            args["fichier_prefix"] = "distribution_moyDistLine_moyHamming_s_";
            p_correls_seuils = lire_p_correls_seuils(rep,motif_p_s);
            for p_correl in p_correls_seuils:
                rep_p_correl = rep + "distribution/";
                distribs.append((rep,rep_p_correl, p_correl, motif_p_s, \
                             correction, priorisation));
    df = pd.DataFrame(); #headers_df = list(); methodes_set = set();
    cpt = 0;
    for distrib in distribs:
        cpt += 1;
        df_tmp = pd.DataFrame();
        df_tmp = create_dataframe_bon_seuils(distrib[0], distrib[1], distrib[2],\
                                distrib[3], distrib[4], distrib[5], args);
        # TODO fusionner les fichiers selon la colonne cpt
        if df.empty:
            df = pd.concat([df,df_tmp], axis=1).astype(np.int32);
        else:
            df = pd.merge( df, df_tmp, on="cpt", how="inner").astype(np.int32);
        del(df_tmp)
    print("distribs cpt = {}".format(cpt))
    return df;          
###############------- test OK car code ecrit a cause de erreur de memory ----------------------
def plot_recherche_bon_seuils_sur_data_brutes(df, args):
    fig, axarr = plt.subplots(3,int(len(args["selected_vars"])/2));
    default_size = fig.get_size_inches();
    fig.set_figheight(default_size[0]*1.5); fig.set_figwidth(default_size[0]*1.5) 
    variabs = "_".join(args["selected_vars"]);
    numeroFigures = ["a","b","c","d","e","f","g","h","i"]
    priorisation = args["priorisation"]; correction = args["correction"]
    for i, var in enumerate(args["selected_vars"]):
        path_save = args["rep"]+ "lineaire_simul50Graphes_priorite_"+\
                    priorisation+"/courbeComparaisonMethodes";
        ####
        REP = args["rep"] + "AEFFACER"  #TODO A EFFACER APRES LES TESTS
        path_save = REP+ "lineaire_simul50Graphes_priorite_"+\
                    priorisation+"/courbeComparaisonMethodes";
        ####
        path_courbe = Path(path_save); path_courbe.mkdir(parents=True,\
                            exist_ok=True);
        
        filter_cols =  [col for col in df.columns \
                        if col.find(var) >= 0 and \
                        col.find(correction) >= 0 and \
                        col.find(priorisation) >= 0];
#        print("1 filter_cols={}".format(filter_cols))
        dico_rename_cols,filter_cols_new = renommer_colonnes(\
                                            filter_cols, "s");
#        print("2 filter_cols={}".format(filter_cols_new))
        df_var = df[filter_cols];
        df_var = df_var.rename(columns = dico_rename_cols);
        
        # sort each column of df_var
        for col in df_var.columns:
            df_var.loc[:,col] = sorted(df_var.loc[:,col])
        
        # define title of graphique 
        dico ={"fauxpositif_correction": "fauxpositf",\
               "fauxnegatif_correction": "fauxnegatif",\
               "fauxpositif_seuil": "fauxpositf",\
               "fauxnegatif_seuil": "fauxnegatif"}
        title = ""
        nom_supp = var.split("_")[1];
        nom = var.split("_")[0];
        if args["langue"] == "francais":
            if nom_supp == "correction":
                nom_supp = "apres algorithme correction"
            elif nom_supp == "seuil":
                nom_supp = "apres seuil"
            title = "comparaison des \n"+ str(nom.upper()) +" "+ \
                    str(nom_supp.upper()) +\
                    " en fonction  du seuil ("+numeroFigures[i]+")"
        elif args["langue"] == "anglais":
            if nom_supp == "correction":
                nom_supp = "after correction algorithm"
            elif nom_supp == "seuil":
                nom_supp = "after threshold";
            title = "comparison of \n"+ str(nom.upper()) +" "+ \
                    str(nom_supp.upper())+\
                    " by thresholds ("+numeroFigures[i]+")"
        #plot
        styles1 = ['bs-','ro-','y^-','rs-','go-','b^-','r*-','bo-','-gD','-yp',\
                   ':>','-.<','-v','-d','-h','--H','--,']
        col_fig = i % (int(len(args["selected_vars"])/2)); 
        row_fig =  int(i / (int(len(args["selected_vars"])/2)))
        print("var={}, w ={}, h={}, row_fig={}, col_fig={}, prior={}, correct={}"\
              .format(var, default_size[0], default_size[1], row_fig, \
                      col_fig, priorisation, correction))
        df_var[filter_cols_new].plot(style=styles1, ax = axarr[row_fig, col_fig]);
        if args["langue"] == "francais":
            axarr[row_fig,col_fig].set(xlabel= "graphes", ylabel= var.upper(), \
                title = title);
            filter_cols = ["_".join(("seuil",item.split("_")[3])) \
                            for item in filter_cols]
        elif args["langue"] == "anglais":
            axarr[row_fig,col_fig].set(xlabel= "graphs", ylabel= var.upper(), \
                title = title);
            filter_cols = ["_".join((item.split("_")[3]), "threshold") for item in filter_cols]
        axarr[row_fig,col_fig].legend( filter_cols, loc='upper center', bbox_to_anchor=(0.5, 1.00),\
                   ncol=4, fancybox=True, shadow=True);
    fig.tight_layout();
    plt.savefig(path_save+"/choix_du_seuil_pour"\
                +"_"+str(variabs)+"_"+correction\
                +"_"+priorisation+".jpeg",dpi= 190);
    plt.clf();
    pass
def plot_recherche_bon_seuils_sur_data_brutes_parLine(df, args):
    fig, axarr = plt.subplots(len(args["selected_vars"]));
    default_size = fig.get_size_inches();
    fig.set_figheight(default_size[0]*2.0); fig.set_figwidth(default_size[0]*2.0) 
    variabs = "_".join(args["selected_vars"]);
    numeroFigures = ["a","b","c","d","e","f","g","h","i"]
    priorisation = args["priorisation"]; correction = args["correction"]
    for i, var in enumerate(args["selected_vars"]):
        path_save = args["rep"]+ "lineaire_simul50Graphes_priorite_"+\
                    priorisation+"/courbeComparaisonMethodes";
        ####
        REP = args["rep"] + "AEFFACER"  #TODO A EFFACER APRES LES TESTS
        path_save = REP+ "lineaire_simul50Graphes_priorite_"+\
                    priorisation+"/courbeComparaisonMethodes";
        ####
        path_courbe = Path(path_save); path_courbe.mkdir(parents=True,\
                            exist_ok=True);
        
        filter_cols =  [col for col in df.columns \
                        if col.find(var) >= 0 and \
                        col.find(correction) >= 0 and \
                        col.find(priorisation) >= 0];
#        print("1 filter_cols={}".format(filter_cols))
        dico_rename_cols,filter_cols_new = renommer_colonnes(\
                                            filter_cols, "s");
#        print("2 filter_cols={}".format(filter_cols_new))
        df_var = df[filter_cols];
        df_var = df_var.rename(columns = dico_rename_cols);
        
        # sort each column of df_var
        for col in df_var.columns:
            df_var.loc[:,col] = sorted(df_var.loc[:,col])
        
        # define title of graphique 
        title = ""
        nom_supp = var.split("_")[1];
        nom = var.split("_")[0];
        if args["langue"] == "francais":
            if nom_supp == "correction":
                nom_supp = "apres algorithme correction"
            elif nom_supp == "seuil":
                nom_supp = "apres seuil"
            title = "comparaison des "+ str(nom.upper()) +" "+ \
                    str(nom_supp.upper()) +\
                    " en fonction  du seuil ("+numeroFigures[i]+")"
        elif args["langue"] == "anglais":
            if nom_supp == "correction":
                nom_supp = "after correction algorithm"
            elif nom_supp == "seuil":
                nom_supp = "after threshold";
            title = "comparison of "+ str(nom.upper()) +" "+ \
                    str(nom_supp.upper())+\
                    " by thresholds ("+numeroFigures[i]+")"
        #plot
        styles1 = ['bs-','ro-','y^-','rs-','go-','b^-','r*-','bo-','-gD','-yp',\
                   ':>','-.<','-v','-d','-h','--H','--,']
        col_fig = 0; #i % (int(len(args["selected_vars"])/2)); 
        row_fig =  i; #int(i / (int(len(args["selected_vars"])/2)))
        print("var={}, w ={}, h={}, row_fig={}, col_fig={}, prior={}, correct={}"\
              .format(var, default_size[0], default_size[1], row_fig, \
                      col_fig, priorisation, correction))
        df_var[filter_cols_new].plot(style=styles1, ax = axarr[row_fig]);
        if args["langue"] == "francais":
            axarr[row_fig].set(xlabel= "graphes", ylabel= var.upper(), \
                            title = title);
            filter_cols = ["_".join(("seuil",item.split("_")[3])) \
                            for item in filter_cols]
        elif args["langue"] == "anglais":
            axarr[row_fig].set(xlabel= "graphs", ylabel= var.upper(), \
                title = title);
            filter_cols = ["_".join((item.split("_")[3]), "threshold") for item in filter_cols]
        axarr[row_fig].legend( filter_cols, loc='upper center', \
                    bbox_to_anchor=(0.5, 1.00), ncol=6, \
                    fancybox=True, shadow=True );
    fig.tight_layout();
    plt.savefig(path_save+"/choix_du_seuil_pour"\
                +"_"+str(variabs)+"_"+correction\
                +"_"+priorisation+".jpeg",dpi= 190);
    plt.clf();
    pass
def formation_dataframe_bon_seuils_sur_data_brutes(reps, args):
    distribs = list();
    for rep in reps:
        correction = rep.split("/")[2]; priorisation = rep.split("/")[1].split("_")[3]
        p_correls_seuils = list(); motif_p_s = ""; 
        if args["bool_p_correl"]: 
            motif_p_s = "p";
            args["fichier_prefix"] = "distribution_moyDistLine_moyHamming_k_";
            p_correls_seuils = lire_p_correls_seuils(rep,motif_p_s);
            for p_correl_s in p_correls_seuils:
                rep_p_correl = rep + "data_"+motif_p_s+"_"+str(p_correl_s)+"/distribution/";
                distribs.append((rep,rep_p_correl, p_correl_s, motif_p_s, \
                             correction, priorisation));
        else:
            motif_p_s = "s";
            args["fichier_prefix"] = "distribution_moyDistLine_moyHamming_s_";
            rep_p_correl = rep + "distribution/";
            distribs.append((rep, rep_p_correl, motif_p_s, \
                             correction, priorisation));
    df = pd.DataFrame(); #headers_df = list(); methodes_set = set();
    cpt = 0;
    for distrib in distribs:
        cpt += 1;
        df_tmp = pd.DataFrame();
        df_tmp = create_dataframe_bon_seuils(distrib[0], distrib[1], distrib[2],\
                                distrib[3], distrib[4], args);
        args["priorisation"] = distrib[4];
        args["correction"]  = distrib[3];
        if args["parLine"] == False:
            plot_recherche_bon_seuils_sur_data_brutes(df_tmp, args)
        else:
            plot_recherche_bon_seuils_sur_data_brutes_parLine(df_tmp, args)
        
    print("distribs cpt = {}".format(cpt))
    return df;     
###############----- test OK car code ecrit a cause de erreur de memory ----------------------

###############################################################################
###     42) comparaison des methodes de priorisation ou fonctions de cout 
###         pour seuils 
###         * NB : on choisit s = 0.7 car j'ai choisi apres analyse  
###                fct_couts = [normale, ajout, suppression]
###############################################################################
def plot_meilleur_priorisation_fct_cout(df_fct_cout_s, priorisations, correction, s, args):    
    for var in args["selected_vars"]:
        filter_cols = []; cpt= 0;
        path_save = args["rep"]+ "courbeComparaisonPriorisation";               # priorisation = fonction de cout
        path_courbe = Path(path_save); 
        path_courbe.mkdir(parents=True, exist_ok=True);
        for priorisation in priorisations:
            filter_cols.extend( [col for col in df_fct_cout_s.columns \
                                 if col.find(var) >= 0 and \
                                 col.find(str(s)) >= 0 and \
                                 col.find(correction) >= 0 and \
                                 col.find(priorisation) >= 0] );
        dico_rename_cols,filter_cols_new = renommer_colonnes(filter_cols,\
                                                             "priorisation");
        df_var = df_fct_cout_s[filter_cols];
        df_var = df_var.rename(columns = dico_rename_cols);
        
        #plot
        fig = plt.figure(); default_size = fig.get_size_inches()
        print("w =", default_size[0], " h = ",default_size[1])
        fig.set_figheight(default_size[0]*1.0); fig.set_figwidth(default_size[0]*1.0) # width = longueur, height = largueur
        styles1 = ['bs-','ro-','y^-','rs-','go-','b^-','r*-','bo-','-gD','-yp',\
                   ':>','-.<','-v','-d','-h','--H','--,']
        cpt += 1; ax1 = fig.add_subplot(cpt,1,1);
        df_var[filter_cols_new].plot(style=styles1, ax = ax1);
        
        # define title of graphique 
        title = ""
        nom_supp = var.split("_")[1];
        nom = var.split("_")[0];
        if args["langue"] == "francais":
            if nom_supp == "correction":
                nom_supp = "apres algorithme correction"
            elif nom_supp == "seuil":
                nom_supp = "apres seuil"
            title = "comparaison des "+ str(nom.upper()) +" "+ \
                    str(nom_supp.upper()) +\
                    " selon les fonctions de couts"
        elif args["langue"] == "anglais":
            if nom_supp == "correction":
                nom_supp = "after correction algorithm"
            elif nom_supp == "seuil":
                nom_supp = "after threshold";
            title = "comparison of "+ str(nom.upper()) +" "+ \
                    str(nom_supp.upper())+\
                    " according the cost functions"
        if args["langue"] == "francais":
            ax1.set(xlabel= "graphes", ylabel= var.upper(), \
                    title = title);
            filter_cols = ["_".join(("fonction",item)) for item in filter_cols_new]
        elif args["langue"] == "anglais":
            ax1.set(xlabel= "graphs", ylabel= var.upper(), \
                    title = title);
            filter_cols = ["_".join((item, "function")) for item in filter_cols_new]
        ax1.legend( filter_cols, loc='upper center', bbox_to_anchor=(0.5, 1.00),\
                   ncol=4, fancybox=True, shadow=True);
        
        priorisations = "_".join(priorisations)
        plt.savefig(path_save+"/comparaison_fct_couts_"+\
                    str(var)+"_s_"+ "".join(str(s).split(".")) +"_"+correction+\
                    "_"+priorisations+".jpeg",dpi= 190); 
        plt.clf();
###############################################################################
###     43) comparaison des methodes de priorisation ou fonctions de cout 
###         pour seuils 
###         * NB : on choisit s = 0.7 car j'ai choisi apres analyse  
###                fct_couts = [normale, unitaire]
###         selected_vars = ["fauxPositif_seuil", "fauxNegatif_seuil", \
###                         "fauxPositif_correction","fauxNegatif_correction", \
###                         "moy_dh"]
###############################################################################        
def plot_meilleur_fct_cout(df, priorisations, correction, s, args):
    fig, axarr = plt.subplots(3,int(len(args["selected_vars"])/2));
    default_size = fig.get_size_inches();
    fig.set_figheight(default_size[0]*1.1); fig.set_figwidth(default_size[0]*1.1) 
    variabs = "_".join(args["selected_vars"]);
    numeroFigures = ["a","b","c","d","e","f","g","h","i"]
    for i, var in enumerate(args["selected_vars"]):
        path_save = args["rep"]+ "/courbeComparaisonPriorisation";
        path_courbe = Path(path_save); path_courbe.mkdir(parents=True,\
                            exist_ok=True);
        filter_cols = list()
        for priorisation in priorisations:
            filter_cols.extend( [col for col in df.columns \
                                 if col.find(var) >= 0 and \
                                 col.find(str(s)) >= 0 and \
                                 col.find(correction) >= 0 and \
                                 col.find(priorisation) >= 0] );
        dico_rename_cols,filter_cols_new = renommer_colonnes(\
                                            filter_cols, "priorisation");
        df_var = df[filter_cols];
        df_var = df_var.rename(columns = dico_rename_cols);
        
        # sort each column of df_var
        for col in df_var.columns:
            df_var.loc[:,col] = sorted(df_var.loc[:,col])
        
        #plot
        styles1 = ['bs-','ro-','y^-','rs-','go-','b^-','r*-','bo-','-gD','-yp',\
                   ':>','-.<','-v','-d','-h','--H','--,']
        col_fig = i % (int(len(args["selected_vars"])/2)); 
        row_fig =  int(i / (int(len(args["selected_vars"])/2)))
        print("var={}, w ={}, h={}, row_fig={}, col_fig={}, correct={}"\
              .format(var, default_size[0], default_size[1], row_fig, \
                      col_fig, correction))
        df_var[filter_cols_new].plot(style=styles1, ax = axarr[row_fig, col_fig]);
        
        # define title of graphique 
        title = ""
        nom_supp = var.split("_")[1];
        nom = var.split("_")[0];
        if args["langue"] == "francais":
            if nom_supp == "correction":
                nom_supp = "apres algorithme correction"
            elif nom_supp == "seuil":
                nom_supp = "apres seuil"
            title = "comparaison des \n"+ str(nom.upper()) +" "+ \
                    str(nom_supp.upper()) +\
                    " selon les fonctions de couts ("+numeroFigures[i]+")"
        elif args["langue"] == "anglais":
            if nom_supp == "correction":
                nom_supp = "after correction algorithm"
            elif nom_supp == "seuil":
                nom_supp = "after threshold";
            title = "comparison of \n"+ str(nom.upper()) +" "+ \
                    str(nom_supp.upper())+\
                    " according the cost functions ("+numeroFigures[i]+")"
        if args["langue"] == "francais":
            axarr[row_fig,col_fig].set(xlabel= "graphes", ylabel= var.upper(), \
                title = title);
            filter_cols = ["_".join(("fonction",item.split("_")[5])) \
                            for item in filter_cols]
        elif args["langue"] == "anglais":
            axarr[row_fig,col_fig].set(xlabel= "graphs", ylabel= var.upper(), \
                title = title);
            filter_cols = ["_".join((item.split("_")[5]), "function") for item in filter_cols]
        axarr[row_fig,col_fig].legend( filter_cols, loc='upper center', bbox_to_anchor=(0.5, 1.00),\
                   ncol=4, fancybox=True, shadow=True);
    fig.tight_layout();
    plt.savefig(path_save+"/choix_fct_cout_pour"\
                +"_"+str(variabs)+"_"+correction\
                +"_s_"+ "".join(str(s).split(".")) +".jpeg",dpi= 190);
    plt.clf();
    pass
def plot_meilleur_fct_cout_parLine(df, priorisations, correction, s, args):
    fig, axarr = plt.subplots(len(args["selected_vars"]));
    default_size = fig.get_size_inches();
    fig.set_figheight(default_size[0]*1.8); fig.set_figwidth(default_size[0]*1.8) 
    variabs = "_".join(args["selected_vars"]);
    numeroFigures = ["a","b","c","d","e","f","g","h","i"]
    for i, var in enumerate(args["selected_vars"]):
        path_save = args["rep"]+ "/courbeComparaisonPriorisation";
        path_courbe = Path(path_save); path_courbe.mkdir(parents=True,\
                            exist_ok=True);
        filter_cols = list()
        for priorisation in priorisations:
            filter_cols.extend( [col for col in df.columns \
                                 if col.find(var) >= 0 and \
                                 col.find(str(s)) >= 0 and \
                                 col.find(correction) >= 0 and \
                                 col.find(priorisation) >= 0] );
        dico_rename_cols,filter_cols_new = renommer_colonnes(\
                                            filter_cols, "priorisation");
        df_var = df[filter_cols];
        df_var = df_var.rename(columns = dico_rename_cols);
        
        # sort each column of df_var
        for col in df_var.columns:
            df_var.loc[:,col] = sorted(df_var.loc[:,col])
        
        #plot
        styles1 = ['bs-','ro-','y^-','rs-','go-','b^-','r*-','bo-','-gD','-yp',\
                   ':>','-.<','-v','-d','-h','--H','--,']
        col_fig = 0; #i % (int(len(args["selected_vars"])/2)); 
        row_fig =  i; #int(i / (int(len(args["selected_vars"])/2)))
        print("var={}, w ={}, h={}, row_fig={}, col_fig={}, correct={}"\
              .format(var, default_size[0], default_size[1], row_fig, \
                      col_fig, correction))
        df_var[filter_cols_new].plot(style=styles1, ax = axarr[row_fig]);
        
        # define title of graphique 
        title = ""
        nom_supp = var.split("_")[1];
        nom = var.split("_")[0];
        if args["langue"] == "francais":
            if nom_supp == "correction":
                nom_supp = "apres algorithme correction"
            elif nom_supp == "seuil":
                nom_supp = "apres seuil"
            title = "comparaison des "+ str(nom.upper()) +" "+ \
                    str(nom_supp.upper()) +\
                    " selon les fonctions de couts ("+numeroFigures[i]+")"
        elif args["langue"] == "anglais":
            if nom_supp == "correction":
                nom_supp = "after correction algorithm"
            elif nom_supp == "seuil":
                nom_supp = "after threshold";
            title = "comparison of "+ str(nom.upper()) +" "+ \
                    str(nom_supp.upper())+\
                    " by the cost functions ("+numeroFigures[i]+")"
        if args["langue"] == "francais":
            axarr[row_fig].set(xlabel= "graphes", ylabel= var.upper(), \
                title = title);
            filter_cols = ["_".join(("fonction",item.split("_")[5])) \
                            for item in filter_cols]
        elif args["langue"] == "anglais":
            axarr[row_fig].set(xlabel= "graphs", ylabel= var.upper(), \
                title = title);
            filter_cols = ["_".join((item.split("_")[5]), "function") for item in filter_cols]
        axarr[row_fig].legend( filter_cols, loc='upper center', bbox_to_anchor=(0.5, 1.00),\
                   ncol=6, fancybox=True, shadow=True);
    fig.tight_layout();
    plt.savefig(path_save+"/choix_fct_cout_pour"\
                +"_"+str(variabs)+"_"+correction\
                +"_s_"+ "".join(str(s).split(".")) +".jpeg",dpi= 190);
    plt.clf();
    pass
###############################################################################
###     5) comparaison des operations de priorisations pour 
###             * une methode de correction (aleatoire)
###             * et un p_correl donnee
###############################################################################
def plot_comparaison_priorisation(df, p_correl, correction, args):
    corrections, priorisations = set(), set();
    corrections, priorisations = trouver_corrections_priorisations(df.columns);
    
    for var in args["selected_vars"]:
        filter_cols = []; cpt= 0;
        path_save = args["rep"]+ "courbeComparaisonPriorisation";
        path_courbe = Path(path_save); path_courbe.mkdir(parents=True, exist_ok=True);
        for priorisation in priorisations:
            filter_cols.extend( [col for col in df.columns if col.find(var) >= 0 and \
                                col.find(str(p_correl)) >= 0 and \
                                col.find(correction) >= 0 and \
                                col.find(priorisation) >= 0] );
#        print("1 filter_cols={}".format(filter_cols))
        dico_rename_cols,filter_cols_new = renommer_colonnes(filter_cols,"priorisation");
#        print("2 filter_cols={}".format(filter_cols_new))
         
        df_var = df[filter_cols];
        df_var = df_var.rename(columns = dico_rename_cols);
        
        #plot
        fig = plt.figure(); default_size = fig.get_size_inches()
        print("w =", default_size[0], " h = ",default_size[1])
        fig.set_figheight(default_size[0]*1.1); fig.set_figwidth(default_size[0]*1.1) # width = longueur, height = largueur
        styles1 = ['bs-','ro-','y^-','rs-','go-','b^-','r*-','bo-','-gD','-yp',\
                   ':>','-.<','-v','-d','-h','--H','--,']
        cpt += 1; ax1 = fig.add_subplot(cpt,1,1);
        df_var[filter_cols_new].plot(style=styles1, ax = ax1);
        if args["langue"] == "francais":
            ax1.set(xlabel= "nombre de cases modifiees", ylabel= var.upper(), \
                title = "comparaison des "+str(var.upper())+" selon les valeurs de p_correls");
            filter_cols = ["_".join(("methode",item)) for item in filter_cols_new]
        elif args["langue"] == "anglais":
            ax1.set(xlabel= "number of modified cases", ylabel= var.upper(), \
                title = "comparison of "+ str(var.upper()) +" according to p_correls values");
            filter_cols = ["_".join((item, "method")) for item in filter_cols_new]
        ax1.legend( filter_cols, loc='upper center', bbox_to_anchor=(0.5, 1.00),\
                   ncol=4, fancybox=True, shadow=True);
        plt.savefig(path_save+"/comparaison_priorisation_"+\
                    str(var)+"_p_correl_"+str(p_correl)+"_"+correction+".jpeg",dpi= 190); # TODO aremplacer p_correl par s
        plt.clf();
        
###############################################################################
###     6) relation entre moy_DH et moy_DL selon 
###             * une methode de correction (aleatoire)
###             * et un p_correl donnee
###             * et en utilisant correl_dh_dl cumule du dataframe
###############################################################################
def plot_relation_moyDH_moyDL(k_errors,args, p_correl, name_p_correl, priorisation, correction):   
    fig = plt.figure(1); default_size = fig.get_size_inches()
    print("w =", default_size[0], " h = ",default_size[1])
    fig.set_size_inches( (default_size[0]*1.0, default_size[1]*1.0) )
    fig.set_size_inches( (default_size[0]*0.4, default_size[1]*0.4) )

    styles1 = ['bs-','ro-','y^-','rs-','go-','b^-','r*-','bo-','-gD','-yp',\
                   ':>','-.<','-v','-d','-h','--H','--,']
    colors = ['b','r','y','b','g','b','r']
    linestyles = ['-.','--','-',':','-','steps']; 
    markers = ['+','x','*','>','^','<','*']
    for ind_k, k in enumerate(k_errors):
        ax1 = fig.add_subplot(1,1,1);
        print("rep = {}".format(args["rep"]+"lineaire_simul50Graphes_priorite_"+\
                         priorisation+"/"+correction+"/data_"+name_p_correl+"_"+\
                         str(p_correl)+"/distribution/"+args["distrib_name"]+str(k)+args["ext"]))
        df = pd.read_csv(args["rep"]+"lineaire_simul50Graphes_priorite_"+\
                         priorisation+"/"+correction+"/data_"+name_p_correl+"_"+\
                         str(p_correl)+"/distribution/"+args["distrib_name"]+str(k)+args["ext"], 
                         names=["cpt","moy_dl","moy_dh", "nbre_aretes_matE", "correl_dh_dl"], \
                         sep=';')
        data_sort = df["correl_dh_dl"].sort_values(ascending = True);
    
        ax1.step(data_sort, data_sort.cumsum(), color = colors[ind_k],\
                 linestyle = linestyles[ind_k], marker = markers[ind_k])

    N_graphs = df["correl_dh_dl"].count()
    ax1.set_yticklabels(['{:3.2f}%'.format(x*100/N_graphs) for x in ax1.get_yticks()])
    for ax in [ax1]:
        if args["langue"] == "anglais":
            ax.set(xlabel= "correlation_DC_DH", ylabel= "cumulative correlation", \
                title = "correlation cumulative function \n \
                        between moy_dc and moy_dh for "+name_p_correl+" = "+\
                        str(p_correl))
        elif args["langue"] == "francais":
            ax.set(xlabel= "correlation_DC_DH", ylabel= "correlation cumulative", \
                   title = "fonction cumulative de la correlation \n \
                           entre moy_dc et moy_dh pour "+name_p_correl+" = "+ \
                           str(p_correl));
        ax.legend( k_errors, loc='upper center', bbox_to_anchor=(0.5, 1.00),\
                   ncol=4, fancybox=True, shadow=True, fontsize = "x-large");
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
         ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(16)
    path_save = args["rep"]+"lineaire_simul50Graphes_priorite_"+priorisation+"/"+\
                correction+"/courbes/"
    plt.savefig(path_save+"correlation_dh_dl_p_"+"".join(str(p_correl).split("."))+\
                ".jpeg", dpi= 190, bbox_inches='tight')
    plt.clf();
    
###############################################################################
###     7) courbe dune distribution asymetrique 
###             * vers la droite alpha = 5    
###             * vers la gauche alpha = -3.9
###
###############################################################################    
def get_truncated_normal(mean=0, sd=1, low=0, upp=10):
    return truncnorm(
        (low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)

def skew_norm_pdf(x,loc=0,scale=1,a=0):
    # adapated from:
    # http://stackoverflow.com/questions/5884768/skew-normal-distribution-in-scipy
    t = (x-loc) / scale
    return 2.0 * scale * stats.norm.pdf(t) * stats.norm.cdf(a*t)

def plot_distribution_asymetrique(path_save, a_inf=0, a_supp=1, SKEW_ALPHA=(-5,3.7), location=0, scale=2):
    """
    a_inf, b_sup : borne inf/sup des valeurs de correls
    location :  la decalage par rapport a la distribution de base
    scale : facteur multiplicatif des valeurs
    SKEW_ALPHA : la pente de la courbe. -5: case a 0; 5: case a 1(la courbe de pdf n'est pas la bonne pour 5)
    """
    x = np.linspace(a_inf, a_supp, 100);
    p = [];
    
    fig, arrax = plt.subplots(1, 2);
    for i,ALPHA in enumerate(SKEW_ALPHA):
        # generate the skew normal PDF for reference:
        cases01 = 1 if ALPHA > 0 else 0; # >0 : asymetrie vers la droite, 
                                         # <0 : aymetrie vers la gauche
        if cases01 == 1:
            location=1.50; scale=2;
        else :
            location=0; scale=2
        p = skew_norm_pdf(x,location,scale,ALPHA)    
        arrax[i].plot(x,p);
        if cases01 == 1:
            yticks = [round(x,2) for x in np.linspace(0,0.8,7)]
            arrax[i].set( xlabel = "valeurs de correlation", \
                          ylabel = "densite", \
                          title = "distributions des cases a "+ str(cases01));
            arrax[i].set_yticklabels(yticks);
        else: 
            arrax[i].set( xlabel = "valeurs de correlation", \
                          ylabel = "densite", \
                          title = "distributions des cases a "+ str(cases01));
    fig.tight_layout();
    plt.grid(True);
    fig.savefig(path_save+\
                "distributionsCases01avecCoefficientAsymetries.jpeg",  dpi= 190)
    

def test_multiply_col(dd):
    for col in dd.columns:
        rd = random.choice(np.arange(1,1.2,0.01))
        dd[col] = dd[col].apply( lambda x: x*rd)

###############################################################################
###     8) graphe iourte : comparaison les fonctions de cout selon leur poids 
###             * ajout 1 supp 1       
###             * ajout 1 supp 10       
###             * ajout 10 supp 1       
###############################################################################            
def comparer_dl_prior_suppr_1_10_1_3subplots(args):
    df = pd.DataFrame();cols = []
    fig =  plt.figure(); default_size = fig.get_size_inches();
    for type_fct_cout in args["type_fct_couts"]:
        for facteur_multiplicatif in args["facteur_multiplicatifs"]:
            priorite = type_fct_cout.split("_")[3]; 
            file = args["base_rep"]+ str(priorite) +"_"+ \
                    str(facteur_multiplicatif)+"/aleatoire_iourte/"+\
                    "distribution_moyDistLine_G_k.txt";
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
        
        f, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(default_size[0]*1.3, default_size[1]*1.3));
        
        df_tu = df[ ["dl_"+prior_0+"_"+str(w_i_0),"dl_"+prior_1+"_"+str(w_i_1),\
                     aretes_G_k_del_prior_0, aretes_G_k_del_prior_1] ];
        max_dl = max(df_tu.max());
        df_tu["k"] = range(0,args["k_deep"]+1,1); #fig = plt.figure(figsize = (8,8)); 
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
        
        f.tight_layout();
        f.savefig(args["path_save"]+"comparaison_prior_"+prior_0+"_wi_"+str(w_i_0)+\
                    "_"+prior_1+"_wi_"+ str(w_i_1) +\
                    "_distance_line_vs_3k_6_graphe_iourte.jpeg", dpi= 190);
        plt.clf();
    pass 

###############################################################################
###     9) graphe  cellule : comparaison les fonctions de cout selon leur poids 
###             * ajout 1 supp 1       
###             * ajout 1 supp 10       
###             * ajout 10 supp 1       
############################################################################### 
def dataframe_supp_ajout_cellules(args):
    df = pd.DataFrame(); cols = list();
    for type_fct_cout in args["type_fct_couts"]:
        for facteur_multiplicatif in args["facteur_multiplicatifs"]:
            priorite = type_fct_cout.split("_")[3]; 
            
            file = args["base_rep"]+ str(priorite) +"_"+ \
                    str(facteur_multiplicatif)+"/aleatoire_cellule/"+\
                    "distribution_moyDistLine_G_k.txt";
            dl_name = "dc_"+str(priorite)+"_"+str(facteur_multiplicatif);
            aretes_G_k_notIn_LG = "aretes_G_k_notIn_LG_"+str(priorite)+"_"+str(facteur_multiplicatif);
            df_tmp = pd.read_csv(file, names=["G_k","k",dl_name,"nb_aretes",aretes_G_k_notIn_LG],\
                                              sep = ";");
            if args["bool_nbre_aretes_pourcentage"]:
                df_tmp.loc[:,dl_name] /=  df_tmp.loc[:,"nb_aretes"]
                df_tmp.loc[:,aretes_G_k_notIn_LG] /=  df_tmp.loc[:,"nb_aretes"]
            ## supprimer les lignes doublons     
            df_err = pd.DataFrame();
            for k in set(df_tmp["k"]):                                          # supprimer les rows en doublons
                df_err_k = df_tmp.loc[df_tmp["k"] == k]
#                df_err_k.sort_values(colonne_dl, inplace=True);
                df_err_k.sort_values("k", inplace=True);
                df_err_k.drop_duplicates(["k"], keep="first", inplace =True)
                df_err = pd.concat([df_err, df_err_k])
                del(df_err_k)
            df_err.sort_values("k", inplace=True);            
            if df.empty:
                df = pd.concat([df,df_err], axis=1);
            else:
                df = pd.merge( df, df_err, on="k", how="inner");
#                print("merge df, df_err .shape={}".format(df.shape))
                del(df_err)

            cols.append("dc_"+str(priorite)+"_"+str(facteur_multiplicatif));
            print("priorite = {}, facteur_multiplicatif={} termine"\
                  .format(priorite,facteur_multiplicatif));
    return df, cols;
def comparer_dl_prior_suppr_1_10_1_3subplots_cellules(args):
    df = pd.DataFrame();cols = []
    fig =  plt.figure(); default_size = fig.get_size_inches();
#    for type_fct_cout in args["type_fct_couts"]:
#        for facteur_multiplicatif in args["facteur_multiplicatifs"]:
#            priorite = type_fct_cout.split("_")[3]; 
#            
#            file = args["base_rep"]+ str(priorite) +"_"+ \
#                    str(facteur_multiplicatif)+"/aleatoire_cellule/"+\
#                    "distribution_moyDistLine_G_k.txt";
#            dl_name = "dl_"+str(priorite)+"_"+str(facteur_multiplicatif);
#            aretes_G_k_notIn_LG = "aretes_G_k_notIn_LG_"+str(priorite)+"_"+str(facteur_multiplicatif);
#            df_tmp = pd.read_csv(file, names=["G_k","k",dl_name,"nb_aretes",aretes_G_k_notIn_LG],\
#                                              sep = ";");
#            df = pd.concat( [df, df_tmp], axis = 1);
#            cols.append("dl_"+str(priorite)+"_"+str(facteur_multiplicatif));
#            print("priorite = {}, facteur_multiplicatif={} termine".format(priorite,facteur_multiplicatif));
    styles1 = ['bs-','r*-','--H','ro-','y^-','rs-','go-','b^-','bo-','-gD','-yp',\
               ':>','-.<','-v','-d','-h','--H','--,'];
    df, cols = dataframe_supp_ajout_cellules(args);
    w_i = [1,10]; prior = ["ajout","supp"];
    cols = ["dl_ajout_1","dl_supp_1","dl_ajout_10","dl_supp_10"];
    for tu in it.combinations(it.product(w_i,prior),2):
        prior_0 = tu[0][1]; prior_1 = tu[1][1];
        w_i_0 = tu[0][0]; w_i_1 = tu[1][0];
        aretes_G_k_del_prior_0 = "aretes_G_k_notIn_LG_"+str(prior_0)+"_"+str(w_i_0);
        aretes_G_k_del_prior_1 = "aretes_G_k_notIn_LG_"+str(prior_1)+"_"+str(w_i_1);
        
        f, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(default_size[0]*1.3, default_size[1]*1.3));
        
        df_tu = df[ ["dl_"+prior_0+"_"+str(w_i_0),"dl_"+prior_1+"_"+str(w_i_1),\
                     aretes_G_k_del_prior_0, aretes_G_k_del_prior_1] ];
        max_dl = int(max(df_tu.max()));
#        df_tu["k"] = range(0,args["k_deep"]+1,1); #fig = plt.figure(figsize = (8,8)); 
        print("df_tu shape:{}, max_dl={} ".format(df_tu.shape, max_dl))
        df_tu["k"] = range(0,df_tu.shape[0],1); df_tu.fillna(0, inplace = True);
        df_tu["droite_3k+6"] = df_tu.apply(lambda row: equation_3k_6(row), axis=1);
        # ax1
        df_tu.plot(x="k",y=["droite_3k+6"], ax = ax1);
        df_tu.plot(x="k", y=["dl_"+prior_0+"_"+str(w_i_0),"dl_"+prior_1+"_"+str(w_i_1)], style= styles1, ax=ax1);
        ax1.set(xlabel= "G_k \n(a)", ylabel= "distance line", \
                xticks=range(0, df_tu["k"].max(),4),yticks=range(0,max_dl,150))
        ax1.grid(True, which='both')
        # ax2
        df_tu.plot(x="k",y=[aretes_G_k_del_prior_0, aretes_G_k_del_prior_1], style= styles1, ax=ax2);
        ax2.set(xlabel= "G_k \n(b)", ylabel= "pourcentage aretes supprimees", \
                xticks=range(0, df_tu["k"].max(),4), yticks=range(0,100,10))
        #ax3
        df_tu["ecart_"+prior_0+"_"+str(w_i_0)] = df_tu["dl_"+prior_0+"_"+str(w_i_0)] - df_tu["droite_3k+6"];
        df_tu["ecart_"+prior_1+"_"+str(w_i_1)] = df_tu["dl_"+prior_1+"_"+str(w_i_1)] - df_tu["droite_3k+6"];
        
        df_tu.plot(x="k",y=["ecart_"+prior_0+"_"+str(w_i_0), "ecart_"+prior_1+"_"+str(w_i_1)], style= styles1, ax=ax3);
        max_ecart = int(max(df_tu["ecart_"+prior_0+"_"+str(w_i_0)].max(), df_tu["ecart_"+prior_1+"_"+str(w_i_1)].max()));
        ax3.set(xlabel= "G_k \n(c)", ylabel= "distance entre 3k+6 et DL", \
                xticks=range(0, df_tu["k"].max(),4),yticks=range(0,max_ecart,150))
        
        f.tight_layout();
        f.savefig(args["path_save"]+"comparaison_prior_"+prior_0+"_wi_"+str(w_i_0)+\
                    "_"+prior_1+"_wi_"+ str(w_i_1) +\
                    "_distance_line_vs_3k_6_graphe_iourte.jpeg", dpi= 190);
        plt.clf();
    pass 

def comparer_dl_prior_suppr_1_10_1_3subplots_cellules_debug(args):
    df = pd.DataFrame();cols = []
    fig =  plt.figure(); default_size = fig.get_size_inches();
#    for type_fct_cout in args["type_fct_couts"]:
#        for facteur_multiplicatif in args["facteur_multiplicatifs"]:
#            priorite = type_fct_cout.split("_")[3]; 
#            
#            file = args["base_rep"]+ str(priorite) +"_"+ \
#                    str(facteur_multiplicatif)+"/aleatoire_cellule/"+\
#                    "distribution_moyDistLine_G_k.txt";
#            dl_name = "dl_"+str(priorite)+"_"+str(facteur_multiplicatif);
#            aretes_G_k_notIn_LG = "aretes_G_k_notIn_LG_"+str(priorite)+"_"+str(facteur_multiplicatif);
#            df_tmp = pd.read_csv(file, names=["G_k","k",dl_name,"nb_aretes",aretes_G_k_notIn_LG],\
#                                              sep = ";");
#            if args["bool_nbre_aretes_pourcentage"]:
#                df_tmp.loc[:,dl_name] /=  df_tmp.loc[:,"nb_aretes"]
#                df_tmp.loc[:,aretes_G_k_notIn_LG] /=  df_tmp.loc[:,"nb_aretes"]
#            ## supprimer les lignes doublons     
#            df_err = pd.DataFrame();
#            for k in set(df_tmp["k"]):                                 # supprimer les rows en doublons
#                df_err_k = df_tmp.loc[df_tmp["k"] == k]
##                df_err_k.sort_values(colonne_dl, inplace=True);
#                df_err_k.sort_values("k", inplace=True);
#                df_err_k.drop_duplicates(["k"], keep="first", inplace =True)
#                df_err = pd.concat([df_err, df_err_k])
#                del(df_err_k)
#            df_err.sort_values("k", inplace=True);            
#            if df.empty:
#                df = pd.concat([df,df_err], axis=1);
#            else:
#                df = pd.merge( df, df_err, on="k", how="inner");
##                print("merge df, df_err .shape={}".format(df.shape))
#                del(df_err)
#
#            cols.append("dl_"+str(priorite)+"_"+str(facteur_multiplicatif));
#            print("priorite = {}, facteur_multiplicatif={} termine"\
#                  .format(priorite,facteur_multiplicatif));
##    return df;
    df, cols = dataframe_supp_ajout_cellules(args);
    styles1 = ['bs-','r*-','--H','ro-','y^-','rs-','go-','b^-','bo-','-gD','-yp',\
               ':>','-.<','-v','-d','-h','--H','--,'];
    w_i = [1,10]; prior = ["ajout","supp"];
    cols = ["dl_ajout_1","dl_supp_1","dl_ajout_10","dl_supp_10"];
    for tu in it.combinations(it.product(w_i,prior),2):
        prior_0 = tu[0][1]; prior_1 = tu[1][1];
        w_i_0 = tu[0][0]; w_i_1 = tu[1][0];
        aretes_G_k_del_prior_0 = "aretes_G_k_notIn_LG_"+str(prior_0)+"_"+str(w_i_0);
        aretes_G_k_del_prior_1 = "aretes_G_k_notIn_LG_"+str(prior_1)+"_"+str(w_i_1);
        
        f, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(default_size[0]*1.3, default_size[1]*1.3));
        
        df_tu = df[ ["dl_"+prior_0+"_"+str(w_i_0),"dl_"+prior_1+"_"+str(w_i_1),\
                     aretes_G_k_del_prior_0, aretes_G_k_del_prior_1] ];
        max_dl = int(max(df_tu.max()));
#        df_tu["k"] = range(0,args["k_deep"]+1,1); #fig = plt.figure(figsize = (8,8)); 
        print("df_tu shape:{}, max_dl={} ".format(df_tu.shape, max_dl))
        df_tu["k"] = range(0,df_tu.shape[0],1); df_tu.fillna(0, inplace = True);
        df_tu["droite_3k+6"] = df_tu.apply(lambda row: equation_3k_6(row), axis=1);
        # ax1
        df_tu.plot(x="k",y=["droite_3k+6"], ax = ax1);
        df_tu.plot(x="k", y=["dl_"+prior_0+"_"+str(w_i_0),"dl_"+prior_1+"_"+str(w_i_1)], style= styles1, ax=ax1);
        yticks = range(0,max_dl,150) \
        if args["bool_nbre_aretes_pourcentage"] == False else np.arange(0,1,0.1)
        ax1.set(xlabel= "G_k \n(a)", ylabel= "distance line", \
                xticks=range(0, df_tu["k"].max(),4),yticks=yticks)
        ax1.grid(True, which='both')
        # ax2
        df_tu.plot(x="k",y=[aretes_G_k_del_prior_0, aretes_G_k_del_prior_1], style= styles1, ax=ax2);
        yticks = range(0,100,10) \
        if args["bool_nbre_aretes_pourcentage"] == False else np.arange(0,1,0.1)
        ax2.set(xlabel= "G_k \n(b)", ylabel= "pourcentage aretes supprimees", \
                xticks=range(0, df_tu["k"].max(),4), yticks=yticks)
        #ax3
        df_tu["ecart_"+prior_0+"_"+str(w_i_0)] = df_tu["dl_"+prior_0+"_"+str(w_i_0)] - df_tu["droite_3k+6"];
        df_tu["ecart_"+prior_1+"_"+str(w_i_1)] = df_tu["dl_"+prior_1+"_"+str(w_i_1)] - df_tu["droite_3k+6"];
        
        df_tu.plot(x="k",y=["ecart_"+prior_0+"_"+str(w_i_0), "ecart_"+prior_1+"_"+str(w_i_1)], style= styles1, ax=ax3);
        max_ecart = int(max(df_tu["ecart_"+prior_0+"_"+str(w_i_0)].max(), df_tu["ecart_"+prior_1+"_"+str(w_i_1)].max()));
        yticks = range(0,max_ecart,150) \
        if args["bool_nbre_aretes_pourcentage"] == False else np.arange(0,1,0.1);
        ax3.set(xlabel= "G_k \n(c)", ylabel= "distance entre 3k+6 et DL", \
                xticks=range(0, df_tu["k"].max(),4),yticks=yticks)
        
        f.tight_layout();
        f.savefig(args["path_save"]+"comparaison_prior_"+prior_0+"_wi_"+str(w_i_0)+\
                    "_"+prior_1+"_wi_"+ str(w_i_1) +\
                    "_distance_line_vs_3k_6_graphe_iourte.jpeg", dpi= 190);
        plt.clf();
    pass 

def comparer_dl_prior_suppr_cellules(args):
    df = pd.DataFrame();cols_ajout_vs_supp = set()
    df, cols = dataframe_supp_ajout_cellules(args);
    
    fig =  plt.figure(); default_size = fig.get_size_inches();
    styles1 = ['bs-','r*-','--H','ro-','y^-','rs-','go-','b^-','bo-','-gD','-yp',\
               ':>','-.<','-v','-d','-h','--H','--,'];
    min_k = 0; #2; 
    w_is = [1]; priorites = ["ajout","supp"];
    print("df.cols={}, \ncols".format(df.columns,cols))
    for prior in priorites:
        for w_i in w_is:
            cols_ajout_vs_supp.add("dl_"+prior+"_"+str(w_i));
            aretes_G_k_del_prior_w_i = "aretes_G_k_notIn_LG_"+str(prior)+"_"+str(w_i);
            aretes_G_k_del_theorik = "aretes_G_k_notIn_LG_theoriques";
            f, (ax1, ax2) = plt.subplots(1,2, \
                                    figsize=(default_size[0]*1.3, \
                                             default_size[1]*1.3));
            
            df_tu = df[ ["dl_"+prior+"_"+str(w_i), aretes_G_k_del_prior_w_i] ];
            max_dl = int(max(df_tu.loc[:,"dl_"+prior+"_"+str(w_i)]));
            print("df_tu shape:{}, max_dl={} ".format(df_tu.shape, max_dl))
            df_tu["k"] = range(min_k,df_tu.shape[0]*2,2); df_tu.fillna(0, inplace = True);
            df_tu["aretes_theorik"] = 2*(df_tu["k"]-1)*(df_tu["k"]) + 4;
            if prior == "supp":
               df_tu["dl_theorique"] = df_tu.apply(lambda row: dl_theorique(row,prior), axis=1); 
               # TODO A calculer le nombre d'aretes supprimes
               df_tu[aretes_G_k_del_theorik] = df_tu.\
               apply(lambda row: dl_supp_theorique(row,prior), axis=1);             
            elif prior == "ajout":
               df_tu["dl_theorique"] = df_tu.apply(lambda row: dl_theorique(row, prior), axis=1);
               # TODO A calculer le nombre d'aretes supprimes
               df_tu[aretes_G_k_del_theorik] = df_tu.\
               apply(lambda row: dl_supp_theorique(row,prior), axis=1);             

            print("prior={},dl_theorique={}".format(prior, df_tu[aretes_G_k_del_theorik]))
            # ax1
            df_tu.plot(x="k",y=["dl_theorique"], ax = ax1);
            df_tu.plot(x="k", y=["dl_"+prior+"_"+str(w_i)], style= styles1, ax=ax1);
            yticks = range(min_k,max_dl,150) \
                        if args["bool_nbre_aretes_pourcentage"] == False \
                        else np.arange(0, 1, 0.1)
            ax1.set(xlabel= "G_k \n(a)", ylabel= "distance line", \
                    xticks=range(min_k, df_tu["k"].max(),4),yticks=yticks)
            ax1.grid(True, which='both')
        
            # ax2
            df_tu.plot(x="k",y=[aretes_G_k_del_prior_w_i, aretes_G_k_del_theorik], \
                       style= styles1, ax=ax2);
            yticks = range(0,100,10) \
                    if args["bool_nbre_aretes_pourcentage"] == False \
                    else np.arange(0,1,0.1)
            ax2.set(xlabel= "G_k \n(b)", ylabel= "pourcentage aretes supprimees", \
                    xticks=range(0, df_tu["k"].max(),4), yticks=yticks)
            
        
            f.tight_layout();
            f.savefig(args["path_save"]+"comparaison_priorite_dl_theorique_vs_"\
                      +prior+"_wi_"+str(w_i)+".jpeg", dpi=190);
            plt.clf();
    # compare dl_ajout_1 vs dl_supp_1 
    df_tu = df[list(cols_ajout_vs_supp)];
    df_tu["k"] = range(0,df_tu.shape[0]*2,2); df_tu.fillna(0, inplace = True);
    max_dl = int(max(df_tu.max()))
    f, ax0 = plt.subplots(1,1, \
                           figsize=(default_size[0]*1.3, \
                                    default_size[1]*1.3));     
    df_tu.plot(x="k", y=cols, style= styles1, ax=ax0);
    yticks = range(0,max_dl,150) \
               if args["bool_nbre_aretes_pourcentage"] == False \
               else np.arange(0, 1, 0.1)
    ax0.set(xlabel= "G_k \n(a)", ylabel= "distance line", \
                    xticks=range(0, df_tu["k"].max(),4),yticks=yticks)
    ax0.grid(True, which='both');
    f.savefig(args["path_save"]+"comparaison_priorite_"+\
              "_vs_".join(list(cols_ajout_vs_supp))+".jpeg", dpi=190);
    pass 
def comparer_dl_prior_suppr_cellules_k_impair(args):
    """
    k  = impiar
    COREECT a chager la taille de la figure 1.3, 1.3
    """
    df = pd.DataFrame(); 
    df, cols = dataframe_supp_ajout_cellules(args);
    df["k"] = 2 * df["k"] - 1; 
    
    fig =  plt.figure(); default_size = fig.get_size_inches();
    styles1 = ['bs-','r*-','--H','ro-','y^-','rs-','go-','b^-','bo-','-gD','-yp',\
               ':>','-.<','-v','-d','-h','--H','--,'];
    styles1 = ['-', '--', '-.', ':']
    labels_pic = ["a","b","c","d","e","f","g","h","i"]; labels_pic.reverse()
    ecart = 2000;
    ######## ----- reprise 
    cols_ajout_vs_supp_calcules = set();
    w_is = [1]; priorites = ["ajout","supp"];

    f, ( (ax_dl_theorik_vs_supp_wi, ax_dl_theorik_vs_ajout_wi),\
         (ax_aretes_del_supp_wi, ax_aretes_del_ajout_wi),\
         (ax_dl_ajout_wi_vs_supp_wi, ax_dl_theorique_ajout_vs_supp) )\
         = plt.subplots(3,2, \
                       figsize=(default_size[0]*1.3, \
                              default_size[1]*1.5));                 
    for prior_wi in it.product(priorites, w_is):
        prior, w_i = prior_wi[0], prior_wi[1]
        cols_ajout_vs_supp_calcules.add("dc_"+prior+"_"+str(w_i));
#        cols_ajout_vs_supp_calcules.add("k");
        aretes_G_k_del_prior_w_i = "aretes_G_k_notIn_LG_"+str(prior)+"_"+str(w_i);
        aretes_G_k_del_theorik = "aretes_G_k_notIn_LG_theoriques_"+prior;
        df_tu = df[ ["k","dc_"+prior+"_"+str(w_i), aretes_G_k_del_prior_w_i] ];
        max_dl = int(max(df_tu.loc[:,"dc_"+prior+"_"+str(w_i)]));
        max_k = int(max(df_tu.loc[:,"k"])); min_k = int(min(df_tu.loc[:,"k"]))
        df_tu["aretes_theorik"] = 2*(df_tu["k"])*(df_tu["k"]+1) + 4;            # si k pair ==> 2*(df_tu["k"]-1)*(df_tu["k"]) + 4; si k est pair
        if prior == "supp":
            df_tu["dl_theorique_"+prior] = df_tu.apply(lambda row: dl_theorique(row,prior), axis=1); 
            df_tu[aretes_G_k_del_theorik] = df_tu.apply(lambda row: \
                                        dl_supp_theorique(row,prior), axis=1); # calculer le nombre THEORIQUE d'aretes supprimes
           
            df_tu.plot(x="k",y=["dl_theorique_"+prior], \
                       ax = ax_dl_theorik_vs_supp_wi);
            df_tu.plot(x="k", y=["dc_"+prior+"_"+str(w_i)], style= styles1, \
                       ax = ax_dl_theorik_vs_supp_wi);
            yticks = range(min_k,max_dl, ecart) \
                        if args["bool_nbre_aretes_pourcentage"] == False \
                        else np.arange(0, 1, 0.1)
            title = "comparaison distances line theoriques vs calcule \n"+ \
                    "pour la fonction suppression";
            ax_dl_theorik_vs_supp_wi.set( xlabel= "G_k\n("+str(labels_pic.pop())+")", \
                                          ylabel= "distance line", \
                                          xticks=range(min_k, max_k, 4),\
                                          yticks=yticks, \
                                          title = title)
            ax_dl_theorik_vs_supp_wi.grid(True, which='both');
            
            df_tu.plot(x="k",y=[ aretes_G_k_del_prior_w_i, \
                                aretes_G_k_del_theorik], \
                       style = styles1, ax = ax_aretes_del_supp_wi);
            yticks = range(0,100,10) \
                    if args["bool_nbre_aretes_pourcentage"] == False \
                    else np.arange(0,1,0.1)
            title = "comparaison aretes supprimes \n theoriques vs calcule \n"+ \
                    "pour la fonction suppression"
            ax_aretes_del_supp_wi.set(xlabel= "G_k\n("+str(labels_pic.pop())+")", \
                                         ylabel= "pourcentage aretes supprimees", \
                                         xticks=range(min_k, max_k,4), \
                                         yticks=yticks, \
                                         title = title)
        elif prior == "ajout":
            df_tu["dl_theorique_"+prior] = df_tu.apply(lambda row: \
                                            dl_theorique(row, prior), axis=1);
            df_tu[aretes_G_k_del_theorik] = df_tu.apply(lambda row: \
                                         dl_supp_theorique(row,prior), axis=1); # calculer le nombre THEORIQUE d'aretes supprimes
           
            df_tu.plot(x="k",y=["dl_theorique_"+prior], \
                       ax = ax_dl_theorik_vs_ajout_wi);
            df_tu.plot(x="k", y=["dc_"+prior+"_"+str(w_i)], style= styles1, \
                       ax = ax_dl_theorik_vs_ajout_wi);
            yticks = range(min_k,max_dl, ecart) \
                        if args["bool_nbre_aretes_pourcentage"] == False \
                        else np.arange(0, 1, 0.1)
            title = "comparaison distances line theoriques vs calcule \n"+ \
                    "pour la fonction ajout";            
            ax_dl_theorik_vs_ajout_wi.set( xlabel= "G_k\n("+str(labels_pic.pop())+")", \
                                          ylabel= "distance line", \
                                          xticks=range(min_k, max_k, 4),\
                                          yticks=yticks,\
                                          title = title)
            ax_dl_theorik_vs_ajout_wi.grid(True, which='both');
            
            df_tu.plot(x="k",y=[ aretes_G_k_del_prior_w_i, \
                                aretes_G_k_del_theorik], \
                       style = styles1, ax = ax_aretes_del_ajout_wi);
            yticks = range(0,100,10) \
                        if args["bool_nbre_aretes_pourcentage"] == False \
                        else np.arange(0,1,0.1);
            title = "comparaison aretes supprimees \n theoriques vs calcule \n"+ \
                    " pour la fonction ajout"
            ax_aretes_del_ajout_wi.set(xlabel= "G_k \n("+str(labels_pic.pop())+")", \
                                         ylabel= "pourcentage aretes supprimees", \
                                         xticks=range(min_k, max_k,4), \
                                         yticks=yticks, \
                                         title = title)
    
    # compare dl_ajout_1 vs dl_supp_1 
    tmp_cols = list(cols_ajout_vs_supp_calcules); tmp_cols.append("k")
#    df_tu = df[list(cols_ajout_vs_supp_calcules).append("k")];
    print("cols_ajout_vs_supp_calcules={}".format(cols_ajout_vs_supp_calcules))
    df_tu = df[tmp_cols]
    max_dl = int(max(df_tu.max()))
    df_tu.plot(x="k", y=list(cols_ajout_vs_supp_calcules), style= styles1, \
               ax=ax_dl_ajout_wi_vs_supp_wi);
    yticks = range(0,max_dl, ecart) \
               if args["bool_nbre_aretes_pourcentage"] == False \
               else np.arange(0, 1, 0.1)
    title = "comparaison distances line calculees \n"+ \
             "entre les fonctions ajout et suppression"
    ax_dl_ajout_wi_vs_supp_wi.set( xlabel= "G_k \n("+str(labels_pic.pop())+")", \
                                  ylabel= "distance line", \
                                  xticks=range(df_tu["k"].min(),df_tu["k"].max(),4),\
                                  yticks=yticks,\
                                  title = title)
    ax_dl_ajout_wi_vs_supp_wi.grid(True, which='both');
    
    # comparer dl theoriques entre les priorites ajout et suppression
    print("df_tu={}".format(df_tu.head()))
    df_tu["dl_theorique_ajout"] = df_tu.apply(lambda row: \
                                            dl_theorique(row,"ajout"), axis = 1)
    df_tu["dl_theorique_supp"] = df_tu.apply(lambda row: \
                                            dl_theorique(row,"supp"), axis = 1)
    
    df_tu.plot(x="k",y=["dl_theorique_supp"], \
            ax = ax_dl_theorique_ajout_vs_supp);
    df_tu.plot(x = "k", y = ["dl_theorique_ajout"], style= styles1, \
               ax = ax_dl_theorique_ajout_vs_supp);
    max_dl = int(max(df_tu[["dl_theorique_ajout","dl_theorique_supp"]].max()))
    yticks = range(0,max_dl, ecart) \
                        if args["bool_nbre_aretes_pourcentage"] == False \
                        else np.arange(0, 1, 0.1)
    title = "comparaison distances lines theoriques \n"+ \
                "entre les fonctions ajout et suppression"
    ax_dl_theorique_ajout_vs_supp.set(xlabel= "G_k \n("+str(labels_pic.pop())+")", \
                                      ylabel= "distance line", \
                                      xticks=range(df_tu["k"].min(),df_tu["k"].max(), 4),\
                                      yticks=yticks,\
                                      title = title)
    ax_dl_theorique_ajout_vs_supp.grid(True, which='both');
    
    f.tight_layout();
    f.savefig(args["path_save"]+"comparaison_distances_line_graphes_cellules_k_2x50"+\
              ".jpeg", dpi=190);
    ######## ----- reprise          
               
def comparer_dl_prior_suppr_cellules_k_pair(args):
    """
    COREECT a chager la taille de la figure 1.3, 1.3
    """
    df = pd.DataFrame(); 
    df, cols = dataframe_supp_ajout_cellules(args);
    df["k"] = 2 * df["k"];
    
    fig =  plt.figure(); default_size = fig.get_size_inches();
    styles1 = ['bs-','r*-','--H','ro-','y^-','rs-','go-','b^-','bo-','-gD','-yp',\
               ':>','-.<','-v','-d','-h','--H','--,'];
    styles1 = ['-', '--', '-.', ':']
    labels_pic = ["a","b","c","d","e","f","g","h","i"]; labels_pic.reverse()
    ecart = 2000; ecart_G_k = 6; ylabel_distance = "distances";
    ######## ----- reprise 
    cols_ajout_vs_supp_calcules = set();
    w_is = [1]; priorites = ["ajout","supp"];

    f, ( (ax_dl_theorik_vs_supp_wi, ax_dl_theorik_vs_ajout_wi),\
         (ax_aretes_del_supp_wi, ax_aretes_del_ajout_wi),\
         (ax_dl_ajout_wi_vs_supp_wi, ax_dl_theorique_ajout_vs_supp) )\
         = plt.subplots(3,2, \
                       figsize=(default_size[0]*1.3, \
                              default_size[1]*1.5));                 
    for prior_wi in it.product(priorites, w_is):
        prior, w_i = prior_wi[0], prior_wi[1]
        cols_ajout_vs_supp_calcules.add("dc_"+prior+"_"+str(w_i));
#        cols_ajout_vs_supp_calcules.add("k");
        aretes_G_k_del_prior_w_i = "aretes_G_k_notIn_LG_"+str(prior)+"_"+str(w_i);
        aretes_G_k_del_theorik = "aretes_G_k_notIn_LG_theoriques_"+prior;
        df_tu = df[ ["k","dc_"+prior+"_"+str(w_i), aretes_G_k_del_prior_w_i] ];
        max_dl = int(max(df_tu.loc[:,"dc_"+prior+"_"+str(w_i)]));
        max_k = int(max(df_tu.loc[:,"k"])); min_k = int(min(df_tu.loc[:,"k"]))
        xticks = range(min_k, max_k, ecart_G_k);
        df_tu["aretes_theorik"] = 2*(df_tu["k"])*(df_tu["k"]+1) + 4;            # si k pair ==> 2*(df_tu["k"]-1)*(df_tu["k"]) + 4; si k est pair
        if prior == "supp":
            df_tu["dl_theorique_"+prior] = df_tu.apply(lambda row: dl_theorique(row,prior), axis=1); 
            df_tu[aretes_G_k_del_theorik] = df_tu.apply(lambda row: \
                                        dl_supp_theorique(row,prior), axis=1); # calculer le nombre THEORIQUE d'aretes supprimes
           
            df_tu.plot(x="k",y=["dl_theorique_"+prior], \
                       ax = ax_dl_theorik_vs_supp_wi);
            df_tu.plot(x="k", y=["dc_"+prior+"_"+str(w_i)], style= styles1, \
                       ax = ax_dl_theorik_vs_supp_wi);
            yticks = range(min_k,max_dl, ecart) \
                        if args["bool_nbre_aretes_pourcentage"] == False \
                        else np.arange(0, 1, 0.1)
            title = "comparaison distances line theoriques vs "+\
                    "distance de correction calcule \n"+ \
                    "pour l'operation suppression uniquement";
            ax_dl_theorik_vs_supp_wi.set( xlabel= "G_k\n("+str(labels_pic.pop())+")", \
                                          ylabel= ylabel_distance, \
                                          xticks=xticks,\
                                          yticks=yticks, \
                                          title = title)
            ax_dl_theorik_vs_supp_wi.grid(True, which='both');
            
            df_tu.plot(x="k",y=[ aretes_G_k_del_prior_w_i, \
                                aretes_G_k_del_theorik], \
                       style = styles1, ax = ax_aretes_del_supp_wi);
            yticks = range(0,100,10) \
                    if args["bool_nbre_aretes_pourcentage"] == False \
                    else np.arange(0,1,0.1)
            title = "comparaison aretes supprimes \n theoriques vs calcule \n"+ \
                    "pour l'operation suppression uniquement"
            ax_aretes_del_supp_wi.set(xlabel= "G_k\n("+str(labels_pic.pop())+")", \
                                         ylabel= "pourcentage aretes supprimees", \
                                         xticks= xticks, \
                                         yticks=yticks, \
                                         title = title)
        elif prior == "ajout":
            df_tu["dl_theorique_"+prior] = df_tu.apply(lambda row: \
                                            dl_theorique(row, prior), axis=1);
            df_tu[aretes_G_k_del_theorik] = df_tu.apply(lambda row: \
                                         dl_supp_theorique(row,prior), axis=1); # calculer le nombre THEORIQUE d'aretes supprimes
           
            df_tu.plot(x="k",y=["dl_theorique_"+prior], \
                       ax = ax_dl_theorik_vs_ajout_wi);
            df_tu.plot(x="k", y=["dc_"+prior+"_"+str(w_i)], style= styles1, \
                       ax = ax_dl_theorik_vs_ajout_wi);
            yticks = range(min_k,max_dl, ecart) \
                        if args["bool_nbre_aretes_pourcentage"] == False \
                        else np.arange(0, 1, 0.1)
            title = "comparaison distances line theoriques vs calcule \n"+ \
                    "pour l'operation ajout uniquement";            
            ax_dl_theorik_vs_ajout_wi.set( xlabel= "G_k\n("+str(labels_pic.pop())+")", \
                                          ylabel= ylabel_distance, \
                                          xticks= xticks,\
                                          yticks=yticks,\
                                          title = title)
            ax_dl_theorik_vs_ajout_wi.grid(True, which='both');
            
            df_tu.plot(x="k",y=[ aretes_G_k_del_prior_w_i, \
                                aretes_G_k_del_theorik], \
                       style = styles1, ax = ax_aretes_del_ajout_wi);
            yticks = range(0,100,10) \
                        if args["bool_nbre_aretes_pourcentage"] == False \
                        else np.arange(0,1,0.1);
            title = "comparaison aretes supprimees \n theoriques vs calcule \n"+ \
                    " pour la fonction ajout uniquement"
            ax_aretes_del_ajout_wi.set(xlabel= "G_k \n("+str(labels_pic.pop())+")", \
                                         ylabel= "pourcentage aretes supprimees", \
                                         xticks=xticks, \
                                         yticks=yticks, \
                                         title = title)
    
    # compare dl_ajout_1 vs dl_supp_1 
    tmp_cols = list(cols_ajout_vs_supp_calcules); tmp_cols.append("k")
#    df_tu = df[list(cols_ajout_vs_supp_calcules).append("k")];
    print("cols_ajout_vs_supp_calcules={}".format(cols_ajout_vs_supp_calcules))
    df_tu = df[tmp_cols]
    max_dl = int(max(df_tu.max()))
    df_tu.plot(x="k", y=list(cols_ajout_vs_supp_calcules), style= styles1, \
               ax=ax_dl_ajout_wi_vs_supp_wi);
    xticks = range(df_tu["k"].min(),df_tu["k"].max(),ecart_G_k)
    yticks = range(0,max_dl, ecart) \
               if args["bool_nbre_aretes_pourcentage"] == False \
               else np.arange(0, 1, 0.1)
    title = "comparaison distances de correction calculees \n"+ \
             "entre les operations ajout et suppression uniquement"
    ax_dl_ajout_wi_vs_supp_wi.set( xlabel= "G_k \n("+str(labels_pic.pop())+")", \
                                  ylabel= ylabel_distance, \
                                  xticks=xticks,\
                                  yticks=yticks,\
                                  title = title)
    ax_dl_ajout_wi_vs_supp_wi.grid(True, which='both');
    
    # comparer dl theoriques entre les priorites ajout et suppression
    print("df_tu={}".format(df_tu.head()))
    df_tu["dl_theorique_ajout"] = df_tu.apply(lambda row: \
                                            dl_theorique(row,"ajout"), axis = 1)
    df_tu["dl_theorique_supp"] = df_tu.apply(lambda row: \
                                            dl_theorique(row,"supp"), axis = 1)
    
    df_tu.plot(x="k",y=["dl_theorique_supp"], \
            ax = ax_dl_theorique_ajout_vs_supp);
    df_tu.plot(x = "k", y = ["dl_theorique_ajout"], style= styles1, \
               ax = ax_dl_theorique_ajout_vs_supp);
    max_dl = int(max(df_tu[["dl_theorique_ajout","dl_theorique_supp"]].max()))
    xticks = range(df_tu["k"].min(),df_tu["k"].max(), ecart_G_k)
    yticks = range(0,max_dl, ecart) \
                        if args["bool_nbre_aretes_pourcentage"] == False \
                        else np.arange(0, 1, 0.1)
    title = "comparaison distances lines theoriques \n"+ \
            "entre les operations ajout et suppression uniquement"
    ax_dl_theorique_ajout_vs_supp.set(xlabel= "G_k \n("+str(labels_pic.pop())+")", \
                                      ylabel= ylabel_distance, \
                                      xticks=xticks,\
                                      yticks=yticks,\
                                      title = title)
    ax_dl_theorique_ajout_vs_supp.grid(True, which='both');
    
    f.tight_layout();
    f.savefig(args["path_save"]+"comparaison_distances_line_graphes_cellules_k_2x50"+\
              ".jpeg", dpi=190);
    ######## ----- reprise   
    ######## test renommage nom columns df 
def comparer_dl_prior_suppr_cellules_k_pair_new(args):
    """
    COREECT a chager la taille de la figure 1.3, 1.3
    """
    df = pd.DataFrame(); 
    df, cols = dataframe_supp_ajout_cellules(args);
    df["k"] = 2 * df["k"];
    
    fig =  plt.figure(); default_size = fig.get_size_inches();
    styles1 = ['bs-','r*-','--H','ro-','y^-','rs-','go-','b^-','bo-','-gD','-yp',\
               ':>','-.<','-v','-d','-h','--H','--,'];
    styles1 = ['-', '--', '-.', ':']
    labels_pic = ["a","b","c","d","e","f","g","h","i"]; labels_pic.reverse()
    ecart = 2000; ecart_G_k = 6; ylabel_distance = "distances";
    ######## ----- reprise 
    cols_ajout_vs_supp_calcules = set();
    w_is = [1]; priorites = ["ajout","supp"];

    f, ( (ax_dl_theorik_vs_supp_wi, ax_dl_theorik_vs_ajout_wi),\
         (ax_aretes_del_supp_wi, ax_aretes_del_ajout_wi),\
         (ax_dl_ajout_wi_vs_supp_wi, ax_dl_theorique_ajout_vs_supp) )\
         = plt.subplots(3,2, \
                       figsize=(default_size[0]*1.3, \
                              default_size[1]*1.5));                 
    for prior_wi in it.product(priorites, w_is):
        prior, w_i = prior_wi[0], prior_wi[1]
        cols_ajout_vs_supp_calcules.add("dc_"+prior+"_"+str(w_i));
#        cols_ajout_vs_supp_calcules.add("k");
        aretes_G_k_del_prior_w_i = "aretes_G_k_notIn_LG_"+str(prior)+"_"+str(w_i);
        aretes_G_k_del_theorik = "aretes_G_k_notIn_LG_theoriques_"+prior;
        df_tu = df[ ["k","dc_"+prior+"_"+str(w_i), aretes_G_k_del_prior_w_i] ];
        max_dl = int(max(df_tu.loc[:,"dc_"+prior+"_"+str(w_i)]));
        max_k = int(max(df_tu.loc[:,"k"])); min_k = int(min(df_tu.loc[:,"k"]))
        xticks = range(min_k, max_k, ecart_G_k);
        df_tu["aretes_theorik"] = 2*(df_tu["k"])*(df_tu["k"]+1) + 4;            # si k pair ==> 2*(df_tu["k"]-1)*(df_tu["k"]) + 4; si k est pair
        if prior == "supp":
            df_tu["dl_theorique_"+prior] = df_tu.apply(lambda row: dl_theorique(row,prior), axis=1); 
            df_tu[aretes_G_k_del_theorik] = df_tu.apply(lambda row: \
                                        dl_supp_theorique(row,prior), axis=1); # calculer le nombre THEORIQUE d'aretes supprimes
            
            # renommage df colums --> debut
            aretes_G_k_del_prior_w_i = "aretes_supprimees_"+str(prior)+"_"+str(w_i);
            aretes_G_k_del_theorik = "borne_sup_"+prior;
            dl_theorique_prior = "dl_borne_sup_"+prior
            dico={"dl_theorique_"+prior:dl_theorique_prior,\
                  "aretes_G_k_notIn_LG_"+str(prior)+"_"+str(w_i) :aretes_G_k_del_prior_w_i,\
                  "aretes_G_k_notIn_LG_theoriques_"+prior : aretes_G_k_del_theorik}
            df_tu.rename(columns= dico, inplace=True)
            # renommage df colums --> debut
            
            df_tu.plot(x="k",y=[dl_theorique_prior], \
                       ax = ax_dl_theorik_vs_supp_wi);
            df_tu.plot(x="k", y=["dc_"+prior+"_"+str(w_i)], style= styles1, \
                       ax = ax_dl_theorik_vs_supp_wi);
            yticks = range(min_k,max_dl, ecart) \
                        if args["bool_nbre_aretes_pourcentage"] == False \
                        else np.arange(0, 1, 0.1)
            title = "comparaison borne superieure vs \n"+\
                    "distance de correction calculee \n"+ \
                    "pour l'operation suppression uniquement";
            ax_dl_theorik_vs_supp_wi.set( xlabel= "G_k\n("+str(labels_pic.pop())+")", \
                                          ylabel= ylabel_distance, \
                                          xticks=xticks,\
                                          yticks=yticks, \
                                          title = title)
            ax_dl_theorik_vs_supp_wi.grid(True, which='both');
            
            df_tu.plot(x="k",y=[ aretes_G_k_del_prior_w_i, \
                                aretes_G_k_del_theorik], \
                       style = styles1, ax = ax_aretes_del_supp_wi);
            yticks = range(0,100,10) \
                    if args["bool_nbre_aretes_pourcentage"] == False \
                    else np.arange(0,1,0.1)
            title = "comparaison aretes supprimes \n borne superieure vs calculee \n"+ \
                    "pour l'operation suppression uniquement"
            ax_aretes_del_supp_wi.set(xlabel= "G_k\n("+str(labels_pic.pop())+")", \
                                         ylabel= "pourcentage aretes supprimees", \
                                         xticks= xticks, \
                                         yticks=yticks, \
                                         title = title)
        elif prior == "ajout":
            df_tu["dl_theorique_"+prior] = df_tu.apply(lambda row: \
                                            dl_theorique(row, prior), axis=1);
            df_tu[aretes_G_k_del_theorik] = df_tu.apply(lambda row: \
                                         dl_supp_theorique(row,prior), axis=1); # calculer le nombre THEORIQUE d'aretes supprimes
           
            # renommage df colums --> debut
            aretes_G_k_del_prior_w_i = "aretes_supprimees_"+str(prior)+"_"+str(w_i);
            aretes_G_k_del_theorik = "borne_sup_"+prior;
            dl_theorique_prior = "dl_borne_sup_"+prior
            dico={"dl_theorique_"+prior:dl_theorique_prior,\
                  "aretes_G_k_notIn_LG_"+str(prior)+"_"+str(w_i) :aretes_G_k_del_prior_w_i,\
                  "aretes_G_k_notIn_LG_theoriques_"+prior : aretes_G_k_del_theorik}
            df_tu.rename(columns= dico, inplace=True)
            # renommage df colums --> debut
            
            df_tu.plot(x="k",y=[dl_theorique_prior], \
                       ax = ax_dl_theorik_vs_ajout_wi);
            df_tu.plot(x="k", y=["dc_"+prior+"_"+str(w_i)], style= styles1, \
                       ax = ax_dl_theorik_vs_ajout_wi);
            yticks = range(min_k,max_dl, ecart) \
                        if args["bool_nbre_aretes_pourcentage"] == False \
                        else np.arange(0, 1, 0.1)
            title = "comparaison borne superieure vs \ndistances correction calcule \n"+ \
                    "pour l'operation ajout uniquement";            
            ax_dl_theorik_vs_ajout_wi.set( xlabel= "G_k\n("+str(labels_pic.pop())+")", \
                                          ylabel= ylabel_distance, \
                                          xticks= xticks,\
                                          yticks=yticks,\
                                          title = title)
            ax_dl_theorik_vs_ajout_wi.grid(True, which='both');
            
            df_tu.plot(x="k",y=[ aretes_G_k_del_prior_w_i, \
                                aretes_G_k_del_theorik], \
                       style = styles1, ax = ax_aretes_del_ajout_wi);
            yticks = range(0,100,10) \
                        if args["bool_nbre_aretes_pourcentage"] == False \
                        else np.arange(0,1,0.1);
            title = "comparaison aretes supprimees \n borne superieure vs calcule \n"+ \
                    " pour la fonction ajout uniquement"
            ax_aretes_del_ajout_wi.set(xlabel= "G_k \n("+str(labels_pic.pop())+")", \
                                         ylabel= "pourcentage aretes supprimees", \
                                         xticks=xticks, \
                                         yticks=yticks, \
                                         title = title)
    
    # compare dl_ajout_1 vs dl_supp_1 
    tmp_cols = list(cols_ajout_vs_supp_calcules); tmp_cols.append("k")
#    df_tu = df[list(cols_ajout_vs_supp_calcules).append("k")];
    print("cols_ajout_vs_supp_calcules={}".format(cols_ajout_vs_supp_calcules))
    df_tu = df[tmp_cols]
    max_dl = int(max(df_tu.max()))
    df_tu.plot(x="k", y=list(cols_ajout_vs_supp_calcules), style= styles1, \
               ax=ax_dl_ajout_wi_vs_supp_wi);
    xticks = range(df_tu["k"].min(),df_tu["k"].max(),ecart_G_k)
    yticks = range(0,max_dl, ecart) \
               if args["bool_nbre_aretes_pourcentage"] == False \
               else np.arange(0, 1, 0.1)
    title = "comparaison distances de correction calculees \n"+ \
             "entre les operations ajout et suppression uniquement"
    ax_dl_ajout_wi_vs_supp_wi.set( xlabel= "G_k \n("+str(labels_pic.pop())+")", \
                                  ylabel= ylabel_distance, \
                                  xticks=xticks,\
                                  yticks=yticks,\
                                  title = title)
    ax_dl_ajout_wi_vs_supp_wi.grid(True, which='both');
    
    # comparer dl theoriques entre les priorites ajout et suppression
    print("df_tu={}".format(df_tu.head()))
    # renommage df colums --> debut
    dl_theorique_prior_ajout = "dl_borne_sup_ajout"
    dl_theorique_prior_supp = "dl_borne_sup_supp"
    dico={"dl_theorique_ajout" : dl_theorique_prior_ajout,\
          "dl_theorique_supp" : dl_theorique_prior_supp
          }
    df_tu.rename(columns= dico, inplace=True)
    # renommage df colums --> debut
    df_tu[dl_theorique_prior_ajout] = df_tu.apply(lambda row: \
                                            dl_theorique(row,"ajout"), axis = 1)
    df_tu[dl_theorique_prior_supp] = df_tu.apply(lambda row: \
                                            dl_theorique(row,"supp"), axis = 1)
    
    df_tu.plot(x="k",y=[dl_theorique_prior_supp], \
            ax = ax_dl_theorique_ajout_vs_supp);
    df_tu.plot(x = "k", y = [dl_theorique_prior_ajout], style= styles1, \
               ax = ax_dl_theorique_ajout_vs_supp);
    max_dl = int(max(df_tu[[dl_theorique_prior_supp, dl_theorique_prior_ajout]].max()))
    xticks = range(df_tu["k"].min(),df_tu["k"].max(), ecart_G_k)
    yticks = range(0,max_dl, ecart) \
                        if args["bool_nbre_aretes_pourcentage"] == False \
                        else np.arange(0, 1, 0.1)
    title = "comparaison borne superieure \n"+ \
            "entre les operations ajout et suppression uniquement"
    ax_dl_theorique_ajout_vs_supp.set(xlabel= "G_k \n("+str(labels_pic.pop())+")", \
                                      ylabel= ylabel_distance, \
                                      xticks=xticks,\
                                      yticks=yticks,\
                                      title = title)
    ax_dl_theorique_ajout_vs_supp.grid(True, which='both');
    
    f.tight_layout();
    f.savefig(args["path_save"]+"comparaison_distances_line_graphes_cellules_k_2x50"+\
              ".jpeg", dpi=190);    
    ######## test renommage nom columns df       
###############################################################################
###     10) plot distribution des distances de Hamming en fonction de 
###             p = 0.5,
###             correlation = aleatoire
###             priorite = aucune         
############################################################################### 
def plot_distribution_k_erreurs(args):
    path_save = args["path_save"]+args["priorite"]+"/"+args["correction"]+\
                "/courbes/";
    ##
    ##  creation de dataframe contenant les distances de Hamming
    ##
    rep_p_correl = args["path_save"]+args["priorite"]+"/"+args["correction"]+\
                    "/data_p_"+str(args["p"])+"/distribution/";
    motif_p_s = "p";
    df = pd.DataFrame();
    for k_error in args["k_errors"]:
        print("ICI")
        df_err = pd.DataFrame(); df_err_bis = pd.DataFrame();
        df_err_bis = pd.read_csv(rep_p_correl+args["fichier_prefix"]+\
                                 str(k_error)+args["ext"],
             names=["cpt","k",\
                    "moy_dl_p_"+str(k_error)+"_"+args["correction"]+"_"+args["priorite"],\
                    "moy_dh_p_"+str(k_error)+"_"+args["correction"]+"_"+args["priorite"],\
                    "nbre_aretes_matE_p_"+str(k_error)+"_"+args["correction"]+"_"+args["priorite"],\
                    "correl_dh_dl_p_"+str(k_error)+"_"+args["correction"]+"_"+args["priorite"]],\
             sep=";")
        print("prior={},corr={},k_error={},df_err.shape={}".\
              format(args["priorite"],correction,k_error,df_err_bis.shape))
        df_err_bis.loc[:,"cpt"] = df_err_bis.loc[:, "cpt"].apply(lambda x: int(x.split("_")[1]))
        print("head df_err_bis = {}".format(df_err_bis["cpt"].head() ))
        for cpt in set(df_err_bis["cpt"]):
            colonne_moy_dh = "moy_dh_"+motif_p_s+"_"+str(k_error)+"_"\
                                +correction+"_"+args["priorite"];
            df_err_cpt = df_err_bis.loc[df_err_bis["cpt"]==cpt]
            df_err_cpt.sort_values(colonne_moy_dh, inplace=True);
            df_err_cpt.drop_duplicates(["cpt"], keep="first",inplace =True)
            df_err = pd.concat([df_err, df_err_cpt])
        del(df_err_bis)
        df_err.sort_values("cpt", inplace=True);
        print("head df_err = {}, df = {}".format(df_err.shape, df.shape))
        if df.empty:
            df = pd.concat([df,df_err], axis=1);
        else:
            ## a effacer
            cols_comm = set(df.columns).intersection(set(df_err.columns))
            print("k_error={},df_col={},df_err={}, \
                   cols_comm={},df.shape={},df_err.shape={}"\
                  .format(k_error,len(df.columns), len(df_err.columns), \
                          cols_comm,df.shape,df_err.shape));
#            print("cpt df={},cpt df_err={}".format(set(df["cpt"]), set(df_err["cpt"])))      
            ## a effacer
            df = pd.merge( df, df_err, on="cpt", how="inner");
            print("merge df, df_err .shape={}".format(df.shape))
            del(df_err)
    df = df.astype(np.float32);
    
    ## plot
    fig =  plt.figure(); default_size = fig.get_size_inches();
    f, ax1 = plt.subplots(1,1, \
                       figsize=(default_size[0]*0.9, \
                                default_size[1]*0.9))
    select_k_s = ["moy_dh_"+motif_p_s+"_"+str(k)+"_"+args["correction"]\
                  +"_"+args["priorite"] for k in k_errors]
    df_k = df[select_k_s];
    dico_rename_cols = dict();
    for col in df_k.columns:
        dico_rename_cols[col] = "_".join(["p",col.split("_")[3]])
    df_k = df_k.rename(columns = dico_rename_cols);
    
    styles1 = ['bs-','r*-','--H','ro-','y^-','rs-','go-','b^-','bo-','-gD','-yp',\
               ':>','-.<','-v','-d','-h','--H','--,'];
    df_k.plot(style= styles1, ax = ax1)
    title = ""
    if args["langue"] == "francais":
        title = "distribution des MOY_DH en fonction de k erreurs pour \
            le correction aleatoire, la priorite unitaire et \
            la repartition p = 0.5";
    elif args["langue"] == "anglais":
        title = "MOY_DH distribution by k modified cases";
    else:
        title = "";
    min_df_k = int(min(df_k.max())); max_df_k = int(max(df_k.max()))
    ax1.set(xlabel= "number of generated graphs", \
              ylabel= "moy_DH", \
              yticks=range(0, max_df_k, 5),\
              title = title)
    ax1.legend( df_k.columns, loc='upper center', bbox_to_anchor=(0.5, 1.00),\
                       ncol=5, fancybox=True, shadow=True);
    ax1.grid(True, which='both');
    
    print("min={},max ={}".format(min_df_k, max_df_k))
    f.tight_layout();
    f.savefig(path_save+"distribution_MOY_DL_k1_25"+\
              ".jpeg", dpi=190);

###############################################################################
###     11) datas reelles 
###         plot distribution des cases a 0 et 1 du graphe reelles 
###                 
############################################################################### 
def verifier_row_col(row, col, cols_matE_proba, cols_matE_reel, dbg_0_1):
    """
    """
#    print("row={}, col={}".format(row,col))
    if dbg_0_1:
        if row in cols_matE_proba and col in cols_matE_proba and \
            row in cols_matE_reel and col in cols_matE_reel:
            return row, col;
        return None, None;
        pass
    if row.split("->")[0]+"->"+row.split("->")[1] in cols_matE_proba and \
        col.split("->")[0]+"->"+col.split("->")[1] in cols_matE_proba and \
        row.split("->")[0]+"->"+row.split("->")[1] in cols_matE_reel and \
        col.split("->")[0]+"->"+col.split("->")[1] in cols_matE_reel:
#        print("1 verif:{},{} ".format(row.split("->")[0]+"->"+row.split("->")[1], col.split("->")[0]+"->"+col.split("->")[1]))
        return row.split("->")[0]+"->"+row.split("->")[1], \
                col.split("->")[0]+"->"+col.split("->")[1];
    elif row.split("->")[0]+"->"+row.split("->")[1] in cols_matE_proba and \
        col.split("->")[1]+"->"+col.split("->")[0] in cols_matE_proba and \
        row.split("->")[0]+"->"+row.split("->")[1] in cols_matE_reel and \
        col.split("->")[1]+"->"+col.split("->")[0] in cols_matE_reel :
#        print("2 verif:{},{} ".format( row.split("->")[0]+"->"+row.split("->")[1], col.split("->")[1]+"->"+col.split("->")[0]))
        return row.split("->")[0]+"->"+row.split("->")[1], \
                col.split("->")[1]+"->"+col.split("->")[0];
    elif row.split("->")[1]+"->"+row.split("->")[0] in cols_matE_proba and \
        col.split("->")[0]+"->"+col.split("->")[1] in cols_matE_proba and \
        row.split("->")[1]+"->"+row.split("->")[0] in cols_matE_reel and \
        col.split("->")[0]+"->"+col.split("->")[1] in cols_matE_reel :
#        print("3 verif:{},{} ".format( row.split("->")[1]+"->"+row.split("->")[0] ,col.split("->")[0]+"->"+col.split("->")[1]))
        return row.split("->")[1]+"->"+row.split("->")[1], \
                col.split("->")[0]+"->"+col.split("->")[1];
    elif row.split("->")[1]+"->"+row.split("->")[0] in cols_matE_proba and \
        col.split("->")[1]+"->"+col.split("->")[0] in cols_matE_proba and \
        row.split("->")[1]+"->"+row.split("->")[0] in cols_matE_reel and \
        col.split("->")[1]+"->"+col.split("->")[0] in cols_matE_reel :
#        print("4 verif:{},{} ".format( row.split("->")[1]+"->"+row.split("->")[0] ,col.split("->")[1]+"->"+col.split("->")[0]))
        return row.split("->")[1]+"->"+row.split("->")[0], \
                col.split("->")[1]+"->"+col.split("->")[0];
    else:
        return None, None;
        
def  distri_01(row, col, seuil, dico_correl, type_errors, \
               matE_proba, matE_reel, dbg_0_1):
    row, col = verifier_row_col(row, col, matE_proba.columns.tolist(), \
                                matE_reel.columns.tolist(), dbg_0_1);
    if row!=None and col!=None:
        if matE_reel.loc[row,col] == 1 and \
            ((seuil-0.1 <= matE_proba.loc[row,col] and  \
              matE_proba.loc[row,col] < seuil) or \
             (seuil-0.1 <= matE_proba.loc[col,row] and  \
              matE_proba.loc[col,row] < seuil)):
                return "VraiPositive", row, col;
        elif matE_reel.loc[row,col] == 0 and \
            ((seuil-0.1 <= matE_proba.loc[row,col] and  \
              matE_proba.loc[row,col] < seuil) or \
             (seuil-0.1 <= matE_proba.loc[col,row] and  \
              matE_proba.loc[col,row] < seuil) ):
                return "VraiNegative", row, col;
        else:
#            print("----({},{})={}, seuil={}----".format(row,col,matE_proba.loc[row,col], seuil))
            pass
    return "error", row, col;
    
def sous_plot(number_pair):
    """
    return le nombre de lignes w et columns h
    """
#    if number_pair < 1:
#        return None, None;
#    elif number_pair >= 1 and number_pair < 4:
#        return 1, number_pair;
#    elif number_pair >= 4 and number_pair < 10:
#        return 3, math.floor( 10/3);
#    elif number_pair >= 10 and number_pair < 16:
#        return 4, math.floor( 16/4);
#    elif number_pair >= 16 and number_pair < 24:
#        return 4, math.ceil( 24/4);
#    else :
#        return 5, math.ceil( number_pair/5);  
        
    if number_pair < 1:
        return None, None;
    elif number_pair >= 1 and number_pair < 4:
        return 1, number_pair;
    elif number_pair == 4:
        return 2, 2
    elif number_pair > 4 and number_pair < 8:
        return 3, math.floor( 8/3);
    elif number_pair >= 8 and number_pair < 10:
        return 4, math.floor( 10/4);
    elif number_pair >= 10 and number_pair < 16:
        return 4, math.floor( 16/4);
    elif number_pair >= 16 and number_pair < 24:
        return 4, math.ceil( 24/4);
    else :
        return 5, math.ceil( number_pair/5); 
def plot_profilConsommation_distr_01(grandeur, correl, pair_equips, \
                                     chemin_datasets, path_save, error, \
                                     date_debut, date_fin):
    """ plot les profils de consommation des pairs d'arcs """
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
    df = pd.read_csv(chemin_datasets+"dataset_"+grandeur+".csv");print("\n \n \n \n")
#    df = pd.read_csv(chemin_datasets+"shapelet_dataset_"+grandeur+".csv");
    df = df.set_index("timestamp") if "timestamp" in df.columns.tolist() \
                                    else df.set_index("Unnamed: 0");
    df = df.loc[date_debut:date_fin] if date_debut != 0 and date_fin != 0 else df;
    df = (df - df.mean(skipna=True))/ df.std(skipna=True); # df = normalize(df.copy())
    df.rolling(window=20).mean()
    df.fillna(0, inplace=True)
    for pair in pair_equips:
        print("pair={}".format(pair))
        ax = fig.add_subplot(w,h,cpt_ax); cpt_ax += 1;
#        df[list(pair)].plot(ax = ax)
        if (('TGBT2->GF2' in pair and  'TGBT1->DD106' in pair ) or 
            ('TGBT2->GF2' in pair and 'TGBT1->DD108' in pair ) or
            ('R486->R487' in pair and 'R481->R488' in pair ) or 
            ('TGBT1->DD205' in pair and 'TGBT4->MSC3' in pair) ):
              df_pair = df[list(pair)]
              df_pair = df_pair.rolling(window=10).mean()
#              df_pair.rolling_mean(window=20000)
              df_pair.plot(ax = ax) 
              print("ICICI--------------------------------")
        else:
            df[list(pair)].plot(ax = ax)
        ax.set(title = pair)
    fig.tight_layout();
    fig.savefig(path_save+"distrib01_"+error+"_correl_"+str(correl)+\
                "_error_<"+error+">_grandeur_"+grandeur+".jpeg",dpi=70); 
    plt.clf();
    pass        
    
def distribution_case_0_1_graphe_reelle(args):
    """
    distribution, selon la metrique de pearson, des cases a 0 et 1
    """
    if args["dbg_0_1"]:
        file = "champlan";
        rep_root = "/home/willy/topologyLearning/datas/"+file+"/";

        rep_root = args["REP_ROOT"] + file +"/";
        sub_graph = "sous_graphe/" #"sous_graphe/"  {"sous_graphe/", ""}
        metrique_distance = "metrique_wil";
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
               "correl_fauxNeg":correl_fauxNeg, "correl_fauxPos":correl_fauxPos,\
               "dbg_0_1": True};
    
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
                                              type_errors, matE_proba, matE_reel, \
                                              args["dbg_0_1"])
            if type_error != "error":
                dico_correl[type_error].append( (row,col) );
                
        dico[seuil_correl] = dico_correl;
    
    df_err_ = pd.DataFrame(dico, columns = [round(s,1) for s in seuil_correls])
    df_err_ = df_err_.T;
    df_err = df_err_.applymap(lambda x: len(x))
    print("df_err = {}, type={} ".format(df_err, type(df_err)))
#    print("coef_1_arcsNonAdjs = {} ".format( coef_1_arcsNonAdjs ))
#    print("coef_01_arcsAdjs = {} ".format( coef_01_arcsAdjs ))
    
    # plot
    dico_err = {"VraiPositive":1,"VraiNegative":0}
    width = 8; height = 5.5; print("w =", width, " h = ",height)
    fig = plt.figure(figsize=(width*1.1,height*1.1))
    cpt_ax = 1;
    for type_error in dico_err.keys():
        ax1 = fig.add_subplot(1,2,cpt_ax); cpt_ax += 1
#        df_err.loc[:,type_error].plot(ax = ax1, kind = "bar", title = "distribution_"+str(dico_err[type_error]))
#        ax1.set(ylabel= "nombre d'erreurs", xlabel= "valeurs de correlations", title="distribution_"+str(dico_err[type_error]))
        df_err.loc[:,type_error].plot(ax = ax1, kind = "bar", title = "cases_"+str(dico_err[type_error]))
        ax1.set(ylabel= "densite", xlabel= "coefficients de similarite", \
                title="cases_"+str(dico_err[type_error]))
        
        autolabel(ax1.patches, ax1)
    fig.tight_layout()
    fig.savefig(args["path_save"]+"distribution_0_1_"+\
                args["metrique_distance"]+".jpeg",dpi=100) 
    plt.clf();
    
    # comprehension distrib_0 seuil = 0.8 et distrib_1 seuil  = 2, 0.1
    error="VraiNegative"; correl = 1.0; #aretes = df_err_.loc[correl, error]; 
    coef_1_arcsNonAdjs_VN = df_err_.loc[ df_err_.index[-1] ,"VraiNegative"]
    print("coef_1_arcsNonAdjs_VN={}".format(coef_1_arcsNonAdjs_VN[:8]))
    plot_profilConsommation_distr_01(args["grandeur"], correl, \
                                      coef_1_arcsNonAdjs_VN[:4], \
                                      args["chemin_datasets"], \
                                        args["path_save"], \
                                      error,args["date_debut"],args["date_fin"]\
                                     ) # TODO ne sert pas A GRAND CHOSE PEUT ETRE EFFACER
    
    error="VraiPositive"; correl = 0.1; #aretes = df_err_.loc[correl, error]; 
    coef_01_arcsAdjs_VP = df_err_.loc[ df_err_.index[1] ,"VraiPositive"]
    print("coef_01_arcsAdjs_VP={}".format(coef_01_arcsAdjs_VP[:8]))
    plot_profilConsommation_distr_01(args["grandeur"], correl, \
                                     coef_01_arcsAdjs_VP[:8], 
                                     args["chemin_datasets"], args["path_save"], 
                                     error,args["date_debut"],args["date_fin"]\
                                     ) # TODO ne sert pas A GRAND CHOSE PEUT ETRE EFFACER
    

###############################################################################
###     12) datas reelles 
###         distribution des valeurs de correlations sur les vrai positives, 
###             negatives erreurs de correlations
###                 
###############################################################################                           
#### distribution formalisation matE  ===> DEBUT
def distrib_formalisation_matE(args):
    """
    distribution des valeurs de correlations sur les vrai positives, 
        negatives erreurs de correlations
    """
    matE_reel = pd.read_csv(args["chemin_matrices"]+"matE_reelles.csv", \
                            index_col = "Unnamed: 0")
        
    liste_Dico = ujson.load(open(args["chemin_matrices"]+\
                            "correlation_grandeurs_"+args["metrique_distance"]+\
                            ".json","r"));
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
            (matE_reel.loc[arete[0], arete[1]] == 1  or \
             matE_reel.loc[arete[1], arete[0]] == 1):
                max_correl = sorted(correls,reverse=True)[indice_correlations] # value max
                max_correl =  reduce(mul, correls, 1)# produit
                max_correl = np.mean(correls) # means
                dico_errors["vrai_positives"].append( max_correl )
            elif (arete[0] in cols and arete[1] in cols) and \
            (matE_reel.loc[arete[0], arete[1]] == 0  or \
             matE_reel.loc[arete[1], arete[0]] == 0):
                max_correl = sorted(correls,reverse=True)[indice_correlations]  # value max
                max_correl =  reduce(mul, correls, 1)# produit
                max_correl = np.mean(correls) # means
                dico_errors["vrai_negatives"].append( max_correl )
            elif arete[0] not in cols and arete[1] not in cols:
                print("aretes {},{} not belong to matE_Reel"\
                      .format(tuple(arete),arete[1]))
        pass   
    # plot
#        default_size = fig.get_size_inches()
#        print("w =", default_size[0], " h = ",default_size[1])
#        fig.set_size_inches( (default_size[0]*1.1, default_size[1]*1.1) );
    width = 8; height = 5.5; print("w =", width, " h = ",height)
    fig = plt.figure(figsize=(width*1.1,height*1.1))                            #fig.set_size_inches( (width*1.1, height*1.1) )
    cpt_ax = 1;
    for type_error in ["vrai_positives", "vrai_negatives"]:
        ax1 = fig.add_subplot(1,2,cpt_ax); cpt_ax += 1
#        dico_errors.loc[type_error].plot(ax = ax1, kind = "bar", title = type_error+"_error_correlation")
        z = pd.DataFrame(dico_errors[type_error], columns = [type_error])
        z.fillna(0,inplace=True)
        print("z={}".format(z.values))
        bins = np.linspace(0,1,11)
        ax1.hist(z.values,bins=bins)
        ax1.set(ylabel= "nombre d'erreurs", \
                xlabel= "coefficients de similarite",\
                title=type_error+"_erreur_correlation")
        autolabel(ax1.patches, ax1)
    fig.tight_layout()
    fig.savefig(args["path_save"]+\
                "distrib_formalisation_matE_vrai_{positives,negatives}_"+\
                args["metrique_distance"]+"_position_correl_"+\
                str(indice_correlations)+"_"+args["type_cal_matE"]+".jpeg",\
                dpi=300) 
    plt.clf();                      
    pass
#### distribution formalisation matE  ===> FIN                          
               
###############################################################################
###     13) datas reelles 
###         distribution des valeurs de correlations 
###             selon la methode et le seuil
###                 
############################################################################### 
def distri_selon_errors_valeurs(row, col, dico_correl, type_errors, \
                                matE_seuil, matE_reel):
    """
    return un dico dont la cle est le type error de correlation et 
                        la valeur est le nombre d'aretes.
    """
    if row in matE_seuil.columns.tolist() and \
        col in matE_seuil.columns.tolist() and \
        row in matE_reel.columns.tolist() and \
        col in matE_reel.columns.tolist() :
        if matE_reel.loc[row,col] == matE_seuil.loc[row,col] and \
            matE_reel.loc[row,col] == 1:
            dico_correl[type_errors[0]] += 1;
        elif matE_reel.loc[row,col] == matE_seuil.loc[row,col] and \
            matE_reel.loc[row,col] == 0:
            dico_correl[type_errors[1]] += 1;
        elif matE_reel.loc[row,col] != matE_seuil.loc[row,col] and \
            matE_reel.loc[row,col] == 0 and matE_seuil.loc[row,col] == 1:
            dico_correl[type_errors[2]] += 1;
        elif matE_reel.loc[row,col] != matE_seuil.loc[row,col] and \
            matE_reel.loc[row,col] == 1 and matE_seuil.loc[row,col] == 0:
            dico_correl[type_errors[3]] += 1;
        else:
            print("----(",row,",",col,") dont belong to column matE_reelles----")
            print("matE_reel[",row,",",col,"] = ",matE_reel.loc[row,col], \
                  " != matE_seuil[",row,",",col,"] = ",matE_seuil.loc[row,col] )
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
        matE_proba = pd.read_csv(chemin_equipements+"matE.csv",\
                                 index_col = "Unnamed: 0");
        matE_reel = pd.read_csv(chemin_equipements+"matE_reelles.csv",\
                                index_col = "Unnamed: 0");
        methode_correl = "correlationParMorceaux"                               # ou liste_mode = ["correlationParMorceaux","correlationGlissante"]
        path_save = chemin_equipements;
        args ={"matE_proba": matE_proba, "matE_reel": matE_reel, \
               "methode_correl": methode_correl, "path_save": path_save, \
               "metrique_distance": metrique_distance, "dbg": True};
        
    
    matE_proba = args["matE_proba"]; matE_reel = args["matE_reel"]
    if set(matE_proba.columns.tolist()) != set(matE_reel.columns.tolist()):
        print("ERROR matE_proba et matE_reel have different set columns")
        if set(matE_proba.columns.tolist()).issubset( set(matE_reel.columns.tolist())):
            print("===>ICI<===")
#        return 
    seuil_correls = [float("{0:.2f}".format(a)) for a in np.linspace(0,1,11)]
    dico = dict()                                                               # {s:dico_correl} dico_correl = {"VraiPositive":,"VraiNegative":,"FauxPosit"}
    type_errors = list()
    for seuil_correl in seuil_correls:
        matE_seuil = fct_aux.matrice_binaire(matE_proba.copy(), seuil_correl)
        dico_correl = {"VraiPositive":0, "VraiNegative":0,\
                       "FauxPositive":0, "FauxNegative":0};
        type_errors = ["VraiPositive", "VraiNegative", \
                       "FauxPositive","FauxNegative"]                           #list(dico_correl.keys())
        for row, col in fct_aux.range_2d(matE_reel.columns.tolist()):
#            print("ICI----row={},col={}".format(row,col))
            row, col = verifier_row_col(row, col, matE_seuil.columns.tolist(), \
                                        matE_reel.columns.tolist(), \
                                        args["dbg_seuil"])
            if row != None and col != None:
                dico_correl = distri_selon_errors_valeurs(row, col, dico_correl,\
                                                          type_errors, \
                                                          matE_seuil, \
                                                          matE_reel);
        dico[seuil_correl] = dico_correl;
        
    ## transform dico en dataframe
    df = pd.DataFrame(dico, columns = seuil_correls)
    df = df.T

    print("df ={}".format(df))
    print("seuil_correls={}".format(seuil_correls))
    ### plot
    width = 8; height = 5.5;
    print("w =", width, " h = ",height)
    fig = plt.figure(figsize=(width*1.5,height*1.5))
    
    cpt_ax = 1;
    for type_error in type_errors:
        ax1 = fig.add_subplot(2,2,cpt_ax); cpt_ax += 1
        title = type_error+"_ErreurAdjacence";                                  # type_error+"_erreur_correlation"; 
        df[type_error].plot(ax = ax1, kind = "bar", title = title)
        ax1.set(ylabel= "nombre d'erreurs", xlabel= "valeurs de seuils")
        autolabel(ax1.patches, ax1)
    fig.tight_layout()
    fig.savefig(args["path_save"]+\
                "error_par_methode_<"+\
                    args["methode_correl"]+","+args["metrique_distance"]+\
                ">_correlation_"+args["type_cal_matE"]+".jpeg",dpi=300);
    plt.clf();
    pass
#### distribution valeurs {fausses, vrai} {negatives, positives} correlations ####
           
if __name__ == '__main__':
    
    start= time.time();
    langue = "francais";
    bool_p_correl = True; bool_seuil = None; 
#    bool_p_correl = False                                                       # travailler avec des seuil
    bool_k_errors = False; distrib_name = None; ext = ".txt"; seleted_vars = "";
    rep = "data/"; 
    priorisations = [];
    if bool_p_correl:
        rep = "data/"; bool_seuil = False; 
        distrib_name = "distribution_moyDistLine_moyHamming_k_";
        selected_vars = ["moy_dh"];
#        bool_k_errors=False; k_errors = [1,2,5,9];
        bool_k_errors=False; k_errors = [1,2,3,4,5,6,7,8,9];
        priorisations = ["lineaire_simul50Graphes_priorite_supp", \
                         "lineaire_simul50Graphes_priorite_ajout", \
                         "lineaire_simul50Graphes_priorite_unitaire"];
                         
        rep = "dataArticle/"
        bool_k_errors = False; k_errors = [2,5,10,20];
        priorisations = ["lineaire_simul50Graphes_priorite_unitaire"];
    else:
        rep="simulation_seuils/"; bool_seuil = True; 
        distrib_name = "distribution_moyDistLine_moyHamming_s_";
        selected_vars = ["faux_pos_seuil", "faux_neg_seuil", "faux_pos_correct",\
                         "faux_neg_correct", "moy_dh"]
        selected_vars = ["fauxPositif_seuil", "fauxNegatif_seuil", \
                         "fauxPositif_correction","fauxNegatif_correction", \
                         "moy_dh"]
#        selected_vars = ["fauxPositif_correction","fauxNegatif_correction", \
#                         "moy_dh"]
        bool_k_errors=False; k_errors = [0.1,0.5,0.6,0.9];
#        bool_k_errors=True; k_errors = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9];
        priorisations = ["lineaire_simul50Graphes_priorite_supp", \
                      "lineaire_simul50Graphes_priorite_ajout", \
                      "lineaire_simul50Graphes_priorite_normale",\
                      "lineaire_simul50Graphes_priorite_unitaire"];

    booleen_rep = False; priorisation_ = "unitaire"; correction_ = "aleatoire";
    
    corrections = ["degreMin","coutMin","aleatoire"]; corrections = ["aleatoire"]
    reps_ = [ rep+priorisation+"/"+correction+"/" \
            for priorisation in priorisations for correction in corrections];
    reps = reps_.copy();
    if booleen_rep == True :
        reps = [rep for rep in reps_ if rep.find(priorisation_)>0 and \
                rep.find(correction_)>0];
        
    args={"bool_p_correl":bool_p_correl,"bool_seuil":bool_seuil,"ext":ext, \
          "distrib_name":distrib_name,"bool_k_errors":bool_k_errors,\
          "k_errors":k_errors, "selected_vars":selected_vars, \
          "langue":langue, "rep":rep};
    
    bool_distribution = False; #False #True;
    bool_distribution_articleRevue = True; #False #True;
    bool_distribution_k0 = False
    bool_relation_moyDL_moy_DH = False; #True;
    bool_relation_moy_DC_moy_DH_articleRevue = False; #False
    bool_formation_dataframe = False; #True; 
    bool_distribution_FN_FP = False#True; #False; #True;
    bool_bon_choix_seuil = False#True#False#True;
    bool_comparaison_priorisation = False; #True;
    bool_comparaison_p_correl = False; #True;
    bool_distribution_asymetrique = False; #True;
    bool_meilleur_priorisation_fct_cout = False;
    bool_bon_fct_cout_s =  False #False#True #False; #True                 
    bool_graphe_iourte = False; #False #True;
    bool_graphe_cellules = False#False; #True;
    bool_plot_mode_correction_priorite_p = False
    
        
    if bool_distribution:
#        plot_distrib_moyDl_Dh_cumulFct_k_0({"dbg":True});
        distributions(reps, args);
    
    if bool_distribution_articleRevue :
        p_correl = 0.5;  priorisation = "unitaire"; correction = "aleatoire";
        args["rep"] = "dataArticle/"; name_p_correl = "p";
        titre_figure = "reduit";                                                 # reduit : intitule tres court, normal : intitule avec la moyenne, l ecart type .... 
        args["titre_figure"] = titre_figure;
        k_errors = [2,5,10,20]; args["langue"] = "anglais";
        distributions(reps, args);    
        
    selected_p_correl_s = 0.5; args["selected_p_correl_s"] = selected_p_correl_s;
    selected_p_correl_s = 0.7; args["selected_p_correl_s"] = selected_p_correl_s; # ajouter 06/05/2018
    
    reps = reps_.copy();
    if booleen_rep == True :
        reps = [rep for rep in reps_ if rep.find(priorisation_)>0];
    if bool_formation_dataframe:
        df = formation_dataframe(reps, args);
        plot_comparaison_methode_correction(df, args);
        plot_comparaison_p_correls(df,args);
 
#    if bool_distribution_FN_FP:
#        selected_vars = ["k_modified_edges"]; args["selected_vars"] = selected_vars;
#        df_fn_fp = formation_dataframe_taux_FN_FP(reps, args); # FAUX A REFAIRE
  
    if bool_comparaison_p_correl:
        erreurs = ["fausse_negatives","fausse_positives"]
        for erreur in erreurs:
            args["erreur"] = erreur;
            df_fn_fp = formation_dataframe_taux_FN_FP(reps, args); # FAUX A REFAIRE
            plot_comparaison_p_correls(df,args)(df_fn_fp, args);
            
    if bool_comparaison_priorisation:        
        p_correl = 0.5; correction = "aleatoire";
        for p_correl in [0.3, 0.5, 0.7]:
            plot_comparaison_priorisation(df, p_correl, correction, args);

    if bool_relation_moyDL_moy_DH:
        p_correl = 0.5; name_p_correl="p"; correction = "aleatoire"; 
        priorisation = ""; args["rep"] = "data_correl_dc_dh/"
        if name_p_correl == "p" :
            priorisation = "normale";
        elif name_p_correl == "s" :
            priorisation = "normale";
        k_errors = [2,5,10,20,40]; args["langue"] = "francais";
        plot_relation_moyDH_moyDL(k_errors, args, p_correl, name_p_correl, \
                                  priorisation, correction)
        
    if bool_relation_moy_DC_moy_DH_articleRevue:
        p_correl = 0.5;  priorisation = "unitaire"; correction = "aleatoire";
        args["rep"] = "dataArticle/"; name_p_correl = "p"
        k_errors = [2,5,7,10,15,20]; 
        k_errors = [2,5,10,20]; args["langue"] = "anglais"; 
        plot_relation_moyDH_moyDL(k_errors, args, p_correl, name_p_correl, \
                                  priorisation, correction)
        
        
    if bool_distribution_asymetrique:
        path_save = "/home/willy/Documents/latexDoc/redaction/fusion_fichiers/images_fusionChapitres/"; 
        SKEW_ALPHA=(-5,3.7); 
        location = 0; scale = 2; a_inf=0; a_supp=1;
        plot_distribution_asymetrique(path_save, a_inf, a_supp, SKEW_ALPHA, location, scale)

    if bool_distribution_FN_FP :
        df_fn_fp = formation_dataframe(reps, args);
#        dd = plot_comparaison_seuils(df_fn_fp, args)
        plot_comparaison_methode_seuil(df_fn_fp, args)
#        plot_FN_FP(df_fn_fp, args);
    ti = time.time() - start
    if bool_bon_choix_seuil:
        args["parLine"] = True;                                                 # selection de l'affichage des faux_{NEGATIF,POSITIF}_{CORRECTION,SEUIL}
        df_choix = formation_dataframe_bon_seuils_sur_data_brutes(reps,args);
#        plot_recherche_bon_seuils_sur_data_brutes(df_choix,args)
#        plot_seuils(df_fn_fp, args);
    if bool_meilleur_priorisation_fct_cout:
        bool_p_correl = False; args["bool_p_correl"] = bool_p_correl;
        priorisations = ["lineaire_simul50Graphes_priorite_supp", \
                        "lineaire_simul50Graphes_priorite_ajout", \
                        "lineaire_simul50Graphes_priorite_normale"];
        corrections =  ["aleatoire"];
        reps = [ rep+priorisation+"/"+correction+"/" \
                    for priorisation in priorisations \
                    for correction in corrections];
        df_fct_cout_s = formation_df_fct_cout_s(reps,args);
        correction = "aleatoire"; s = 0.7; 
        priorisations = [ p.split("_")[3] for p in priorisations];
        selected_vars = ["moy_dh"]; args["selected_vars"] = selected_vars;
        plot_meilleur_priorisation_fct_cout(df_fct_cout_s, priorisations, \
                                            correction, s, args)       
    if bool_bon_fct_cout_s :
        bool_p_correl = False; args["bool_p_correl"] = bool_p_correl;        
        priorisations = ["lineaire_simul50Graphes_priorite_normale",\
                         "lineaire_simul50Graphes_priorite_unitaire"];
        corrections =  ["aleatoire"];
        reps = [ rep+priorisation+"/"+correction+"/" \
                    for priorisation in priorisations \
                    for correction in corrections];
        df_fct_cout_s = formation_df_fct_cout_s(reps,args);
        selected_vars = ["fauxPositif_seuil", "fauxNegatif_seuil", \
                         "fauxPositif_correction","fauxNegatif_correction", \
                         "moy_dh"]
        args["selected_vars"] = selected_vars; 
        correction = "aleatoire"; s = 0.7; 
        priorisations = [ p.split("_")[3] for p in priorisations];
        bool_parLine = True;
        if bool_parLine == False:
            plot_meilleur_fct_cout(df_fct_cout_s, priorisations, \
                               correction, s, args)
        else:
            plot_meilleur_fct_cout_parLine(df_fct_cout_s, priorisations, \
                               correction, s, args)
        
    if bool_graphe_iourte:
        # redefinir reps
        k_deep = 50;
        facteur_multiplicatifs = [1,2,5,10];
        base_rep = "graphe_particulier/graphe_particulier_prior_"; 
        path_save = "graphe_particulier/";
        type_fct_couts = ["lineare_iourte_priorite_supp",\
                          "lineare_iourte_priorite_ajout"];
        args ={"type_fct_couts":type_fct_couts,\
               "facteur_multiplicatifs":facteur_multiplicatifs,\
               "k_deep":k_deep,"base_rep":base_rep,"path_save":path_save};
        ## priorite ajout[1,10,1]/suppr[1,1,10]
        comparer_dl_prior_suppr_1_10_1_3subplots(args);
        
    if bool_graphe_cellules:
        k_deep = 18;k_deep = 9;
        facteur_multiplicatifs = [1,5,10];facteur_multiplicatifs = [1,10];
        facteur_multiplicatifs = [1];
        base_rep = "graphe_cellules_VM/graphe_cellule_prior_"; 
        path_save = "graphe_cellules_test/";
        type_fct_couts = ["lineare_iourte_priorite_supp",\
                          "lineare_iourte_priorite_ajout"];
        bool_nbre_aretes_pourcentage =  False; # True;
        args ={"type_fct_couts":type_fct_couts,\
               "facteur_multiplicatifs":facteur_multiplicatifs,\
               "k_deep":k_deep,"base_rep":base_rep,"path_save":path_save, 
               "bool_nbre_aretes_pourcentage": bool_nbre_aretes_pourcentage};
        ## priorite ajout[1,10,1]/suppr[1,1,10]
#        df = comparer_dl_prior_suppr_1_10_1_3subplots_cellules(args); 
#        df = comparer_dl_prior_suppr_1_10_1_3subplots_cellules_debug(args);
#        df, cols = dataframe_supp_ajout_cellules(args);
#        comparer_dl_prior_suppr_cellules_k_impair(args);
#        comparer_dl_prior_suppr_cellules_k_pair(args);                         #=====> correct juste changer theorique par borne sup
        comparer_dl_prior_suppr_cellules_k_pair_new(args);
        
    if bool_plot_mode_correction_priorite_p:                              # plot toutes les distributions pour ce tuple (aleatoire, aucune, p=0.5)
        p = 0.5; priorite = "aucune"; correction = "aleatoire";
        k_errors = [1,2,3,4,5,6,7,8,9,15,25]; 
        fichier_prefix = "distribution_moyDistLine_moyHamming_k_";
        ext = ".txt"; langue = "anglais"
        path_save = "/home/willy/Documents/python_topology_learning_simulation/data/lineaire_simul50Graphes_priorite_"
        args={"priorite":priorite,"p":p,"correction":correction,\
              "path_save":path_save,"k_errors":k_errors, "ext":ext,\
              "fichier_prefix": fichier_prefix, "langue":langue\
              }
        plot_distribution_k_erreurs(args)
    print ("runtime = {}".format(time.time() - start))