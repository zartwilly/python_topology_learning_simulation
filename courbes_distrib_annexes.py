#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  5 22:58:43 2018

@author: willy
courbes distributions annexes
fct_cout unitaire, ajout, suppression
"""

import time, os, re;
import numpy as np;
import pandas as pd;
from pathlib import Path;
import random;
import math;
import logging;
import collections as Coll;
import matplotlib.pyplot as plt;
import matplotlib.font_manager as font_manager
import seaborn as sns;
import itertools as it;
from functools import reduce;
from operator import mul;
from scipy.stats import norm, truncnorm
from scipy import stats;


##----------- fonctions annexes -------------- debut
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
    elif motif == "p":
        for col in filter_cols:
            dico[col] = "_".join([col.split("_")[4],col.split("_")[5]]); 
            cols.append(col.split("_")[4]+"_"+col.split("_")[5]);
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
##----------- fonctions annexes -------------- fin

###############################################################################
###     0)          creation dataframe
###
###############################################################################
def create_dataframe(rep, rep_p_correl_s, p_correl_s, motif_p_s, correction, priorisation, args):
    k_errors = [];
    if args["bool_excludes_sup_10cases"]:
        k_errors = [1,2,3,4,5,6,7,8,9]
    else:
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
        dico_means[float(k_error)] = dico_means_k_error;
    return pd.DataFrame(dico_means).transpose()
    
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

def formation_dataframe(reps, args):
    distribs = list();
    for rep in reps:
        correction = rep.split("/")[6]; priorisation = rep.split("/")[5]
        p_correls_seuils = list(); motif_p_s = ""; 
        if args["bool_p"]: 
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
    
###############################################################################
###     1)          distribution des k cases modifiees
###
###############################################################################
def plot_moyDLDH_correlDlDh_cumul(df, k_error, axarr_x, col, num_bins, aretes, motif_p_s, langue):
    if col in ["moy_dh","moy_dl"]:
        sns.distplot(df[col],ax = axarr_x, bins = range(0,int(num_bins)), kde = False)
        (mu, sigma) = norm.fit(df[col]) # best fit of data
        if col == "moy_dh" and langue == "francais" and motif_p_s == "p":
            axarr_x.set(xlabel= "moy_distance_line", ylabel= "nombre_graphe", \
                title = "distance de Hamming pour "+ str(k_error)+\
                "\n cases modifiees \n $\mu=%.3f,\ \sigma=%.3f,\ $ \n $ aretes = %.3f$" \
                %(mu, sigma, aretes))
        elif col == "moy_dh" and langue == "francais" and motif_p_s == "s":
            axarr_x.set(xlabel= "moy_distance_line", ylabel= "nombre_graphe", \
                title = "distance de Hamming pour \n un seuil s="+ str(k_error)+\
                "\n $\mu=%.3f,\ \sigma=%.3f,\ $ \n $ aretes = %.3f$" \
                %(mu, sigma, aretes))
        elif col == "moy_dh" and langue == "anglais" and motif_p_s == "p":
            axarr_x.set(xlabel= "moy_distance_line", ylabel= "number_graphe", \
                title = "Hamming distance for "+ str(k_error)+\
                "\n modified cases \n $\mu=%.3f,\ \sigma=%.3f,\ $ \n $ edges = %.3f$" \
                %(mu, sigma, aretes))
        elif col == "moy_dh" and langue == "anglais" and motif_p_s == "s":
            axarr_x.set(xlabel= "moy_distance_line", ylabel= "number_graphe", \
                title = "Hamming distance for \n threshold s="+ str(k_error)+\
                "\n $\mu=%.3f,\ \sigma=%.3f,\ $ \n $ edges = %.3f$" \
                %(mu, sigma, aretes))
        elif col == "moy_dl" and langue == "francais" and motif_p_s == "p":
            axarr_x.set(xlabel= "moy_distance_line", ylabel= "nombre_graphe", \
                title = "distance line pour "+ str(k_error)+\
                "\n cases modifiees \n $\mu=%.3f,\ \sigma=%.3f,\ $ \n $ aretes = %.3f$" \
                %(mu, sigma, aretes))
        elif col == "moy_dl" and langue == "francais" and motif_p_s == "s":
            axarr_x.set(xlabel= "moy_distance_line", ylabel= "nombre_graphe", \
                title = "distance line pour \n un seuil s = "+ str(k_error)+\
                "\n $\mu=%.3f,\ \sigma=%.3f,\ $ \n $ aretes = %.3f$" \
                %(mu, sigma, aretes))
        elif col == "moy_dl" and langue == "anglais" and motif_p_s == "p":
            axarr_x.set(xlabel= "moy_distance_line", ylabel= "number_graphe", \
                title = "line distance for "+ str(k_error)+\
                "\n modified cases \n $\mu=%.3f,\ \sigma=%.3f,\ $ \n $ edges = %.3f$" \
                %(mu, sigma, aretes))
        elif col == "moy_dl" and langue == "anglais" and motif_p_s == "s":
            axarr_x.set(xlabel= "moy_distance_line", ylabel= "number_graphe", \
                title = "line distance for \n threshold"+ str(k_error)+\
                "\n $\mu=%.3f,\ \sigma=%.3f,\ $ \n $ edges = %.3f$" \
                %(mu, sigma, aretes))
        max_count_dl, max_count_dh = count_max_df(df)
        axarr_x.plot([int(k_error)+1,int(k_error)+1], (0,max_count_dl), 'r--' );
        axarr_x.set_yticklabels(['{:3.2f}%'.format(x*100/df["moy_dl"].count()) \
                                 for x in axarr_x.get_yticks()]);
        return axarr_x;
    elif col in ["correl_dh_dl"]:
        data_sort = df[col].sort_values(ascending = True);
        axarr_x.step(data_sort, data_sort.cumsum())
        if langue == "francais" and motif_p_s == "p":
            axarr_x.set(xlabel= "correlation_DL_DH", ylabel= "cumulative correlation ", \
                title = "fonction de repartition de \ncorrelation entre moy_dl et moy_dh \n pour "+\
                str(k_error)+" cases modifiees")
        elif  langue == "anglais" and motif_p_s == "p":
            axarr_x.set(xlabel= "correlation_DL_DH", ylabel= "cumulative correlation ", \
                title = "cumulative function of \n correlations between moy_dl et moy_dh \n for "+\
                str(k_error)+" modified cases")
        elif langue == "francais" and motif_p_s == "s":
            axarr_x.set(xlabel= "correlation_DL_DH", ylabel= "cumulative correlation ", \
                title = "fonction de repartition de la \ncorrelation entre moy_dl et moy_dh \n pour "+\
                        " un seuil s = "+ str(k_error))
        elif langue == "anglais" and motif_p_s == "s":
            axarr_x.set(xlabel= "correlation_DL_DH", ylabel= "cumulative correlation ", \
                title = "cumulative function of \n correlations between moy_dl et moy_dh \n for "+\
                        "a threshold s = "+str(k_error))
        axarr_x.set_yticklabels(['{:3.2f}%'.format(x*100/df["correl_dh_dl"].count()) \
                    for x in axarr_x.get_yticks()]);
        return axarr_x;
    elif col in ["cumul_dh"]:
        df.sort_values(by = "moy_dh", ascending=True, axis = 0, inplace = True)
        df["nb_graphe_dh<x"] = df["moy_dh"].apply( lambda x: df["moy_dh"][df.moy_dh < x].count()/df["moy_dh"].count())
#        print("--->k={}, cumul_dh => min = {}, max = {},".format(k_error, df["nb_graphe_dh<x"].min(), df["nb_graphe_dh<x"].max()))
        axarr_x.step(df["moy_dh"],df["nb_graphe_dh<x"]);
        if langue == "francais" and motif_p_s == "p":
            axarr_x.set(xlabel= "DH", ylabel= "number graph moy_DH < x ", \
                title = "cumulative moy_dh pour \n"+str(k_error)+" cases modifiees");
        elif langue == "anglais" and motif_p_s == "p":
            axarr_x.set(xlabel= "DH", ylabel= "number graph moy_DH < x ", \
                title = "cumulative moy_dh for \n"+str(k_error)+"  modified cases")
        elif langue == "francais" and motif_p_s == "s":
            axarr_x.set(xlabel= "DH", ylabel= "number graph moy_DH < x ", \
                title = "cumulative moy_dh pour \n un seuil s = "+str(k_error));
        elif langue == "anglais" and motif_p_s == "s":
            axarr_x.set(xlabel= "DH", ylabel= "number graph moy_DH < x ", \
                title = "cumulative moy_dh for \n a threshold s = "+str(k_error))
        axarr_x.set_xticklabels(np.arange(0, df["moy_dh"].count(), 10), rotation=45 ) 
        return axarr_x;
    
def histo_cumul_fct_seaborn(params, k_errors):
    # lire dataframe aretes matE
    colonnes = ["num_graphe","sommets_matE","aretes_matE"];
    path_save_aretes_matE = "./file_results/"; file_save_aretes_matE = "nombres_aretes_line_graphes.csv";
    df_aretes_matE = pd.read_csv(path_save_aretes_matE+file_save_aretes_matE, names = colonnes, skiprows=[0]);
    
    k_errors.sort();                                                           
    k_chunkies = chunkify(k_errors,5);
    for k_chunky in k_chunkies:
        print("k_chunky={}".format(k_chunky));
        fig, axarr = plt.subplots(len(k_chunky), 4); 
        fig.set_figheight(20); fig.set_figwidth(15) # width = longueur, height = largueur  il faut revoir les valeurs
        for ind, k in enumerate(k_chunky):
            aretes = find_aretes(params, df_aretes_matE);
            df = None;
            if params["bool_p"]:
                df = pd.read_csv(params["path_distrib"]+params["file_prefix"]+str(k)+params["ext"], \
                         names=["cpt","moy_dl","moy_dh","nbre_aretes_matE", "correl_dh_dl"], sep=';')
            else :
                df = pd.read_csv(params["path_distrib"]+params["file_prefix"]+str(k)+params["ext"], \
                     names=["cpt","moy_dl","moy_dh","nbre_aretes_matE","correl_dh_dl",\
                            "faux_pos_seuil","faux_neg_seuil","faux_pos_correct", "faux_neg_correct"],\
                     sep=';')
            # ind=0, ax = 0 --> moy_dl
            axarr[ind,0] = plot_moyDLDH_correlDlDh_cumul(df, k, axarr[ind,0], \
                            "moy_dl", df["moy_dl"].max()+1, aretes, \
                            params["motif_p_s"], params["langue"]);
            # ind=0, ax = 1 --> moy_dh
            axarr[ind,1] = plot_moyDLDH_correlDlDh_cumul(df, k, axarr[ind,1], \
                            "moy_dh", df["moy_dh"].max()+1, aretes, \
                            params["motif_p_s"], params["langue"]);
            # ind=0, ax = 2 --> correl_dl_dh
            axarr[ind,2] = plot_moyDLDH_correlDlDh_cumul(df, k, axarr[ind,2], \
                            "correl_dh_dl", df["correl_dh_dl"].max()+1, aretes,\
                            params["motif_p_s"], params["langue"]);
            # ind=0, ax = 3 --> cumul_dh
            axarr[ind,3] = plot_moyDLDH_correlDlDh_cumul(df, k, axarr[ind,3], \
                            "cumul_dh", df["moy_dh"].max()+1, aretes, \
                            params["motif_p_s"], params["langue"]);
            print("k_errors k= {} moy_dl, moy_dh, correl_dh_dl, cumul_dh termine".format(k))
        # save axarr
        fig.tight_layout();
        
        plt.grid(True);
        k_min = k_chunky[0]; k_max = k_chunky[len(k_chunky)-1] ;
#        fig.savefig(params["save_courbe"]+"_"+params["motif_p_s"]+"_"+\
#                    str(params["p_correl"])+"_distanceMoyenDLDH_k_"+\
#                    str(k_min)+"_"+str(k_max)+".jpeg", dpi= 190)
        fig.savefig(params["save_courbe"]+"/distanceMoyenDLDH_k_"+str(k_chunky)\
                    +"_"+params["motif_p_s"]+"_"+str(params["p"])+".jpeg", dpi= 190)
        print("p_correl={}, motif_p_s={}, k_min={}, k_max={}, save_courbe={}".format(\
              params["p"], params["motif_p_s"], k_min, k_max, params["save_courbe"]));
    del(df_aretes_matE);
    
def plot_distribution(distrib, rep, p, motif_p_s, args):
    """
    tracer la distribution pour k cases pour un p_correl
    """
    k_errors = lire_k_errors(distrib) if args["bool_k_errors"] else args["k_errors"];
    k_errors = [1,2,3,4,5,6,7,8,9] if args["bool_k_errors"] else args["k_errors"];
    save_courbe = rep+"courbes/";
    print("save_courbe = {}".format(save_courbe));
    path_save = Path(save_courbe); path_save.mkdir(parents=True, exist_ok=True);
    params = dict();
    params={"save_courbe":save_courbe, "ext":args["ext"], "path_distrib":distrib,\
            "file_prefix":args["distrib_name"], "bool_p":args["bool_p"],\
            "p":p, "motif_p_s":motif_p_s, "langue":args["langue"]};
    histo_cumul_fct_seaborn(params, k_errors);
    
def distributions(reps,args):
    """ forme les chemins vers les distributions pour tous les p
    """
    distribs = list();
    for rep in reps:
#        distribs = list()
        p_correls_seuils = list(); motif_p_s = ""; 
        if args["bool_p"]: 
            motif_p_s = "p";
            p_seuils = lire_p_correls_seuils(rep,motif_p_s);
            for p in p_seuils:
                rep_p_correl = rep+"data_"+motif_p_s+"_"+str(p)+"/distribution/";
                distribs.append((rep,rep_p_correl, p , motif_p_s));
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
    priorisations = ["unitaire"];
    for priorisation in priorisations:
        cpt = 0;
        for var in args["selected_vars"]:
            path_save = args["rep"]+priorisation+"/courbeComparaisonMethodes";
            path_courbe = Path(path_save); path_courbe.mkdir(parents=True, exist_ok=True)
            filter_cols =  [col for col in df.columns if col.find(var) >= 0 and \
                            col.find(str(args["selected_p_s"])) >= 0 and \
                            col.find(priorisation) >= 0];
#            print("1 filter_cols={}".format(filter_cols))
            dico_rename_cols,filter_cols_new = renommer_colonnes(filter_cols,"p");
            print("2 filter_cols_new={}".format(filter_cols_new))
            df_var = df[filter_cols];
            df_var = df_var.rename(columns = dico_rename_cols);
    
            #plot
            fig = plt.figure() 
            default_size = fig.get_size_inches()
            print("w =", default_size[0], " h = ",default_size[1])
            fig.set_figheight(default_size[0]*1.0); fig.set_figwidth(default_size[0]*1.0) # width = longueur, height = largueur
            styles1 = ['bs-','ro-','y^-','rs-','go-','b^-','r*-','bo-','-gD','-yp',\
                       ':>','-.<','-v','-d','-h','--H','--,']
            cpt += 1; ax1 = fig.add_subplot(cpt,1,1);
            df_var[filter_cols_new].plot(style=styles1, ax = ax1);
            if args["langue"] == "francais":
                ax1.set(xlabel= "nombre de cases modifiees", ylabel= var.upper(), \
                    title = "comparaison des modes de correction avec "+ str(var.upper()));
            elif args["langue"] == "anglais":
                ax1.set(xlabel= "number of modified cases", ylabel= var.upper(), \
                    title = "comparison of correction modes with"+ str(var.upper()));
            ax1.legend( filter_cols_new, loc='upper center', bbox_to_anchor=(0.5, 1.00),\
                       ncol=3, fancybox=True, shadow=True);
            plt.savefig(path_save+"/comparaison_methodes_correction_pour_p_"+\
                        str(args["selected_p_s"])+"_"+priorisation+".jpeg",dpi= 190);
            plt.clf();
          
###############################################################################
###     3) comparaison des p pour une methode de correction (aleatoire)
###
###############################################################################
def plot_comparaison_p(df,args):
    corrections, priorisations = trouver_corrections_priorisations(df.columns);
    priorisations = ["unitaire"];
    corrections = ["aleatoire_sansRemise"];
    for priorisation in priorisations:
        for correction in corrections:
            cpt = 0;
            for var in args["selected_vars"]:
                cpt = 0;
                path_save = args["rep"]+priorisation+"/courbeComparaisonMethodes";
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
                plt.savefig(path_save+"/comparaison_p_s_"+\
                            correction+"_"+priorisation+".jpeg",dpi= 190);
                plt.clf();
                
###############################################################################
###     6) relation entre moy_DH et moy_DL selon 
###             * une methode de correction (aleatoire_sansRemise)
###             * et un p donnee
###             * et en utilisant correl_dh_dl cumule du dataframe
###############################################################################
def plot_relation_moyDH_moyDL(k_errors, args, p_correl, name_p_correl, \
                              priorisation, correction):   
    fig = plt.figure(1); default_size = fig.get_size_inches()
    print("w =", default_size[0], " h = ",default_size[1])
    fig.set_size_inches( (default_size[0]*1.0, default_size[1]*1.0) )
    styles1 = ['bs-','ro-','y^-','rs-','go-','b^-','r*-','bo-','-gD','-yp',\
                   ':>','-.<','-v','-d','-h','--H','--,']
    colors = ['b','r','y','b','g','b','r']
    styles = ['-.','--','-',':','-.']
    for ind_k, k in enumerate(k_errors):
        ax1 = fig.add_subplot(1,1,1);
        df = pd.read_csv(args["rep"]+priorisation+"/"+correction+\
                         "/data_"+name_p_correl+"_"+str(p_correl)+"/distribution/"+\
                         args["distrib_name"]+str(k)+args["ext"], 
                         names=["cpt","moy_dl","moy_dh", "nbre_aretes_matE", "correl_dh_dl"], \
                         sep=';')
        data_sort = df["correl_dh_dl"].sort_values(ascending = True);
    
        ax1.step(data_sort, data_sort.cumsum(),color = colors[ind_k],linestyle= styles[ind_k])

    N_graphs = df["correl_dh_dl"].count()
    ax1.set_yticklabels(['{:3.2f}%'.format(x*100/N_graphs) for x in ax1.get_yticks()])
    p_correl = 0.5
    for ax in [ax1]:
        ax1.set(xlabel= "correlation_DL_DH", ylabel= "correlation cumulative", \
                title = "fonction cumulative de la correlation \n \
    entre moy_dl et moy_dh pour "+name_p_correl+" = "+ str(p_correl));
        ax1.legend( k_errors, loc='upper center', bbox_to_anchor=(0.5, 1.00),\
                   ncol=4, fancybox=True, shadow=True);
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
         ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(14)
    path_save = args["rep"]+priorisation+"/"+correction+"/courbes/";
    plt.savefig(path_save+"correlation_dh_dl_p_"+"".join(str(p_correl).split("."))+\
                ".jpeg", dpi= 190, bbox_inches='tight')
    plt.clf();
    
    
###############################################################################
###     7) main
###############################################################################    
if __name__ == '__main__':
    
    start= time.time();
    langue = "francais";
    bool_p = True; bool_seuil = None; bool_excludes_sup_10cases = True
    bool_k_errors=False; distrib_name = None; ext = ".txt"; seleted_vars = "";
    rep = "data/"; 
    if bool_p:
        rep = "/home/willy/Documents/courbePython05052018/"; bool_seuil = False; 
        distrib_name = "distribution_moyDistLine_moyHamming_k_";
        selected_vars = ["moy_dh"];
        bool_k_errors=False; k_errors = [1,2,5,9];
        bool_k_errors=True; k_errors = [1,2,3,4,5,6,7,8,9];
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
    booleen_rep = False; priorisation_ = "unitaire"; correction_ = "aleatoire_sansRemise";
    priorisations = ["unitaire", "ajout", "suppression"]; priorisations = ["unitaire"]
    corrections =  ["degreMin_sansRemise","coutMin_sansRemise",\
                    "degreMin_avecRemise","coutMin_avecRemise",\
                    "aleatoire_sansRemise"];
    reps_ = [ rep+priorisation+"/"+correction+"/" \
            for priorisation in priorisations for correction in corrections];
    reps = reps_.copy();
    if booleen_rep == True :
        reps = [rep for rep in reps_ if rep.find(priorisation_)>0 and rep.find(correction_)>0];
        
    args={"bool_p":bool_p,"bool_seuil":bool_seuil,"ext":ext, \
          "distrib_name":distrib_name,"bool_k_errors":bool_k_errors,\
          "k_errors":k_errors, "selected_vars":selected_vars, \
          "langue":langue, "rep":rep, \
          "bool_excludes_sup_10cases":bool_excludes_sup_10cases};
    
    bool_distribution = True; #False #True;
    bool_formation_dataframe = False; #True; 
    bool_comparaison_p = False; #True;
    bool_relation_moyDL_moy_DH = False; #True, False;
    
    
    if bool_distribution:
        distributions(reps, args);
    
    selected_p_s = 0.5; args["selected_p_s"] = selected_p_s;
    reps = reps_.copy();
    if booleen_rep == True :
        reps = [rep for rep in reps_ if rep.find(priorisation_)>0];
    if bool_formation_dataframe:
        df = formation_dataframe(reps, args);
        plot_comparaison_methode_correction(df, args);
        plot_comparaison_p(df,args);
 
    if bool_relation_moyDL_moy_DH:
        p = 0.7; name_p ="p"; correction = "aleatoire_sansRemise"; priorisation = "unitaire";
        k_errors = [2,5,10,20,40];
        plot_relation_moyDH_moyDL(k_errors, args, p, name_p, priorisation, correction)
    
    
    print ("runtime = {}".format(time.time() - start))