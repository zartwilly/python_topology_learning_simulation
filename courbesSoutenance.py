#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 14:22:58 2018

@author: willy
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

#---------------------- distribution horizontale   =  debut -------------------
def create_dataframe(args):
    dico_means = dict();
    priorisation = args["priorite"]; motif_p_s = args["motif_p_s"];
    p_correl_s = args["p_value"]; correction = args["correction"];
    rep_p_correl_s = args["path_save"]+args["priorite"]+"/"+args["correction"]+\
                    "/data_p_"+str(args["p_value"])+"/distribution/";
    for k_error in args["k_errors"]:
        dico_means_k_error = dict();
        df = pd.DataFrame();
        if motif_p_s == "s":
            df = pd.read_csv(rep_p_correl_s+args["fichier_prefix"]+str(k_error)+args["ext"],\
                 names=["cpt","moy_dc_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
                   "moy_dh_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
                   "nbre_aretes_matE_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
                   "correl_dh_dc_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
                   "fauxPositif_seuil_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
                   "fauxNegatif_seuil_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
                   "fauxPositif_correction_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
                   "fauxNegatif_correction_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation],\
                 sep=";")
        elif motif_p_s == "p":
            df = pd.read_csv(rep_p_correl_s+args["fichier_prefix"]+str(k_error)+args["ext"],
                 names=["cpt","moy_dc_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
                        "moy_dh_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
                        "nbre_aretes_matE_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation,\
                        "correl_dh_dc_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation],\
                 sep=";")
        else: 
            print("motif n existe pas")
            return pd.Dataframe();
        for var in args["selected_vars"]:
            var = var+"_"+motif_p_s+"_"+str(p_correl_s)+"_"+correction+"_"+priorisation;
            dico_means_k_error[var] = df[var].mean();
        dico_means[float(k_error)] = dico_means_k_error;
#    print("dico_means={}, selected_vars={}".format(dico_means,args["selected_vars"]));
    return pd.DataFrame(dico_means).transpose()

from collections import Counter;
def count_max_df(df):
    """
    return the max of distribution df['moy_dl'] and df['moy_dh']
    """
    d_dl = Counter(df["moy_dc"])
    k, max_value_dl = max(d_dl.items(), key=lambda x:x[1])
    
    d_dh = Counter(df["moy_dh"])
    k, max_value_dh = max(d_dh.items(), key=lambda x:x[1])
    
    return max_value_dl, max_value_dh
    pass

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
    nb_rd = 1+1.6*round(random.random(),1);
    nb_min_max = random.randrange(df["aretes_matE"].min(), df["aretes_matE"].max());
    aretes = nb_min_max * nb_rd;
    return aretes;
    pass    

def distribution_horizontale(args):
    """
    mettre les pourcentages sur l'axe y
    nbreFileNotDisplay: nbre de fichier a ne pas afficher
    
    https://stackoverflow.com/questions/31357611/format-y-axis-as-percent
    
    PRobleme : affiche trop petit 1700*1900 pixels
    """
    number_files = [2,5,10,20];
    nbreFileNotDisplay = 0;
    comment = "";
    num_bins = args["num_bins"];
    rep = args["path_save"]+args["correction"]+\
                    "/data_p_"+str(args["p_value"])+"/distribution/";
    w = 4; h = 1;                                                             # width = largueur, height = longueur
    fig = plt.figure( figsize=(w,h) );             
    cpt_ax1 = 0;
    for num in number_files:
        print("num = ", num)
        num = int(num)
        cpt_ax1 += 1;#cpt = num; # cpt += 1
        
        # ax1
        ax1 = fig.add_subplot(2,len(number_files),cpt_ax1);
        df = pd.read_csv(rep+args["fichier_prefix"] +str(num)+args["ext"], \
                         names=["cpt","moy_dc","moy_dh", "nbre_aretes_matE", "correl_dh_dl"], \
                         sep=';')
        N_graphs = df["moy_dc"].count()
        
        # best fit of data
        (mu, sigma) = norm.fit(df["moy_dc"])
        num_bins = df["moy_dc"].max()+1
        bins = range(0,int(num_bins)); bins = range(0, 100)
        print("---> bins = ", bins, " min = ",df["moy_dc"].min(), \
              " max = ",df["moy_dc"].max())
        
        max_count_dl, max_count_dh = count_max_df(df)
        
        sns.distplot(df["moy_dc"], ax = ax1, bins = bins, kde = False)
        ax1.set(xlabel= "moy_distance_correction", ylabel= "nombre_graphe", \
                title = "distance de correction pour \n "+ str(num)+\
                " cases modifiees \n $\mu=%.3f,\ \sigma=%.3f$, " %(mu, sigma)+ \
                " \n $aretes = %.3f$" %(df["nbre_aretes_matE"].mean))
        ax1.plot([num+1,num+1], (0,max_count_dl), 'r--' )
        ax1.set_yticklabels(['{:3.2f}%'.format(x*100/N_graphs) \
                             for x in ax1.get_yticks()])
        
        # ax2
        cpt_ax2 = cpt_ax1 +len(number_files);  #cpt = num+len(number_files); # cpt +=1 ;
        ax2 = fig.add_subplot(2,len(number_files),cpt_ax2);
        N_graphs = df["moy_dh"].count()
        # best fit of data
        (mu, sigma) = norm.fit(df["moy_dh"])
        
        num_bins = df["moy_dh"].max()+1
        bins = range(0 ,int(num_bins)); bins = range(0, 100)

        sns.distplot(df["moy_dh"], ax = ax2, bins = bins, kde = False, color = 'red')
        ax2.set(xlabel= "moy_distance_hamming", ylabel= "nombre_graphe", \
                title = "distance de Hamming pour \n "+ str(num)+ \
                " cases modifiees \n $\mu=%.3f,\ \sigma=%.3f$, " %(mu, sigma) + \
                " \n $aretes = %.3f$" %(df["nbre_aretes_matE"].mean()))
#        ax2.set_xticklabels(bins, rotation=90)
        ax2.plot([num+1,num+1], (0,max_count_dh), 'r--' )
        ax2.set_yticklabels(['{:3.2f}%'.format(x*100/N_graphs) \
                             for x in ax2.get_yticks()])
    
        for ax in [ax1,ax2]:
            for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
                item.set_fontsize(8)
                
#    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.grid(True)
    comment += "_horizontale";
    plt.savefig(args["path_save"]+args["correction"]+"/courbes/"+\
                "distributionHorizontale_k_0_"+str(number_files[len(number_files)-1])+\
                "_"+comment+".jpeg", \
                dpi= 190)
    pass        

def distribution_horizontale_new(args):
    """
    mettre les pourcentages sur l'axe y
    nbreFileNotDisplay: nbre de fichier a ne pas afficher
    
    https://stackoverflow.com/questions/31357611/format-y-axis-as-percent
    """
    number_files = [2,5,10,20];
    rep = args["path_save"]+args["correction"]+\
                    "/data_p_"+str(args["p_value"])+"/distribution/";
    fig =  plt.figure(); default_size = fig.get_size_inches();                
    f, ax_arrs = plt.subplots(2, 4, figsize=(default_size[0]*2.2, \
                                             default_size[1]*1.5), \
                                );
    cpt1 = 0; cpt2 = 1; tab_bins = [20, 40, 80,100, 100, 100]
    for ind, k in enumerate(number_files) :
        df = pd.read_csv(rep+args["fichier_prefix"] +str(k)+args["ext"], \
                         names=["cpt","moy_dc","moy_dh", "aretes_matE", "correl_dh_dl"], \
                         sep=';')
        N_graphs = df["moy_dc"].count();
        
        aretes = find_aretes(args, df)
        bins = range(0, (ind+1)*args["num_bins"]);
        bins = range(0, tab_bins[ind]);
        max_count_dl, max_count_dh = count_max_df(df)
        
        # plot ax1
        (mu, sigma) = norm.fit(df["moy_dc"]);                                   # best fit of data
        sns.distplot(df["moy_dc"], ax = ax_arrs[cpt1, ind], bins = bins, kde = False)
        ax_arrs[cpt1, ind].set(xlabel= "moy_distance_correction", \
                ylabel= "nombre_graphe", \
                title = "distance de correction pour \n "+ str(k)+ \
                " cases modifiees")
        ax_arrs[cpt1, ind].plot([k+1,k+1], (0,max_count_dl), 'r--' )
        ax_arrs[cpt1, ind].set_yticklabels(['{:3.2f}%'.format(x*100/N_graphs) \
                             for x in ax_arrs[cpt1, ind].get_yticks()])
        
        #plot ax2
        (mu, sigma) = norm.fit(df["moy_dh"]);                                   # best fit of data
        sns.distplot(df["moy_dh"], ax = ax_arrs[cpt2, ind], bins = bins, kde = False)
        ax_arrs[cpt2, ind].set(xlabel= "moy_distance_correction", \
                ylabel= "nombre_graphe", \
                title = "distance Hamming pour \n "+ str(k)+\
                " cases modifiees")
        ax_arrs[cpt2, ind].plot([k+1,k+1], (0,max_count_dl), 'r--' )
        ax_arrs[cpt2, ind].set_yticklabels(['{:3.2f}%'.format(x*100/N_graphs) \
                             for x in ax_arrs[cpt2, ind].get_yticks()])
    
#    fig = ax_arrs[0,0].figure ;
#    fig.text(0.5,0.04, "Some very long and even longer xlabel", ha="center", va="center")
#    fig.text(0.05,0.5, "Some quite extensive ylabel", ha="center", va="center", rotation=90)

        
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.grid(True)
    plt.savefig(args["path_save"]+args["correction"]+"/courbes/"+\
                "distributionMoyDCDHp05k1251020.jpeg",\
                dpi= 250) #190                                    
#---------------------- distribution horizontale   =  fin ---------------------
                               
#---------------------- test distribution horizontale   =  debut -------------------
def distribution_horizontale_par_blocs(args):
    """
    mettre les pourcentages sur l'axe y
    nbreFileNotDisplay: nbre de fichier a ne pas afficher
    
    https://stackoverflow.com/questions/31357611/format-y-axis-as-percent
    """
    number_files = [2,5,10,20]; fact_w = 1; fact_h = 1;
    rep = args["path_save"]+args["correction"]+\
                    "/data_p_"+str(args["p_value"])+"/distribution/";
    fig =  plt.figure(); default_size = fig.get_size_inches();                
    f25, ax_arrs_k25 = plt.subplots(2, 2, figsize=(default_size[0]*fact_w, \
                                             default_size[1]*fact_h), \
                                );
    
    cpt1 = 0; cpt2 = 1; tab_bins = [20, 40, 80,100, 100, 100]
    for ind, k in enumerate(number_files[:len(number_files)//2]) :
        df = pd.read_csv(rep+args["fichier_prefix"] +str(k)+args["ext"], \
                         names=["cpt","moy_dc","moy_dh", "aretes_matE", "correl_dh_dl"], \
                         sep=';')
        N_graphs = df["moy_dc"].count();
        
#        bins = range(0, (ind+1)*args["num_bins"]);
        bins = range(0, tab_bins[ind]);
        max_count_dl, max_count_dh = count_max_df(df)
        
        # plot ax1
        (mu, sigma) = norm.fit(df["moy_dc"]);                                   # best fit of data
        sns.distplot(df["moy_dc"], ax = ax_arrs_k25[cpt1, ind], \
                     bins = bins, kde = False)
        ax_arrs_k25[cpt1, ind].set(xlabel= "moy_distance_correction", \
                ylabel= "nombre_graphe", \
                title = "distance de correction pour \n "+ str(k)+ \
                " cases modifiees")
        ax_arrs_k25[cpt1, ind].plot([k+1,k+1], (0,max_count_dl), 'r--' )
        ax_arrs_k25[cpt1, ind].set_yticklabels(['{:3.2f}%'.format(x*100/N_graphs) \
                             for x in ax_arrs_k25[cpt1, ind].get_yticks()])
        
        #plot ax2
        (mu, sigma) = norm.fit(df["moy_dh"]);                                   # best fit of data
        sns.distplot(df["moy_dh"], ax = ax_arrs_k25[cpt2, ind], \
                     bins = bins, kde = False)
        ax_arrs_k25[cpt2, ind].set(xlabel= "moy_distance_correction", \
                ylabel= "nombre_graphe", \
                title = "distance Hamming pour \n "+ str(k)+\
                " cases modifiees")
        ax_arrs_k25[cpt2, ind].plot([k+1,k+1], (0,max_count_dl), 'r--' )
        ax_arrs_k25[cpt2, ind].set_yticklabels(['{:3.2f}%'.format(x*100/N_graphs) \
                             for x in ax_arrs_k25[cpt2, ind].get_yticks()])
    
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.grid(True)
    plt.savefig(args["path_save"]+args["correction"]+"/courbes/"+\
                "distributionMoyDCDHp05k25.jpeg",\
                dpi= 250)
        
    fig =  plt.figure(); default_size = fig.get_size_inches();        
    f1020, ax_arrs_k1020 = plt.subplots(2, 2, figsize=(default_size[0]*fact_w, \
                                             default_size[1]*fact_h), \
                                );    
    for ind, k in enumerate(number_files[len(number_files)//2:]) :
        df = pd.read_csv(rep+args["fichier_prefix"] +str(k)+args["ext"], \
                         names=["cpt","moy_dc","moy_dh", "aretes_matE", "correl_dh_dl"], \
                         sep=';')
        N_graphs = df["moy_dc"].count();
        
        bins = range(0, tab_bins[len(number_files)//2:][ind]);
        max_count_dl, max_count_dh = count_max_df(df)
        
       # plot ax1
        (mu, sigma) = norm.fit(df["moy_dc"]);                                   # best fit of data
        sns.distplot(df["moy_dc"], ax = ax_arrs_k1020[cpt1, ind], \
                     bins = bins, kde = False)
        ax_arrs_k1020[cpt1, ind].set(xlabel= "moy_distance_correction", \
                ylabel= "nombre_graphe", \
                title = "distance de correction pour \n "+ str(k)+ \
                " cases modifiees")
        ax_arrs_k1020[cpt1, ind].plot([k+1,k+1], (0,max_count_dl), 'r--' )
        ax_arrs_k1020[cpt1, ind].set_yticklabels(['{:3.2f}%'.format(x*100/N_graphs) \
                                 for x in ax_arrs_k1020[cpt1, ind].get_yticks()]) 
        
        #plot ax2
        (mu, sigma) = norm.fit(df["moy_dh"]);                                   # best fit of data
        sns.distplot(df["moy_dh"], ax = ax_arrs_k1020[cpt2, ind], \
                     bins = bins, kde = False)
        ax_arrs_k1020[cpt2, ind].set(xlabel= "moy_distance_correction", \
                ylabel= "nombre_graphe", \
                title = "distance Hamming pour \n "+ str(k)+\
                " cases modifiees")
        ax_arrs_k1020[cpt2, ind].plot([k+1,k+1], (0,max_count_dl), 'r--' )
        ax_arrs_k1020[cpt2, ind].set_yticklabels(['{:3.2f}%'.format(x*100/N_graphs) \
                                 for x in ax_arrs_k1020[cpt2, ind].get_yticks()])
#    fig = ax_arrs[0,0].figure ;
#    fig.text(0.5,0.04, "Some very long and even longer xlabel", ha="center", va="center")
#    fig.text(0.05,0.5, "Some quite extensive ylabel", ha="center", va="center", rotation=90)

        
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.grid(True)
    plt.savefig(args["path_save"]+args["correction"]+"/courbes/"+\
                "distributionMoyDCDHp05k1020.jpeg",\
                dpi= 250) #190                      
#---------------------- test distribution horizontale   =  fin ---------------------


#---------------------- profil consommation  =  debut -------------------------
def profil_consommation(params):
    """
    generer deux listes de 200 items a_i et a_k
    multiplier par 1.5 la liste a_i element par element
    plot les trois listes
    """
    N = 50; min_val = 120; max_val = 200;
    df = pd.DataFrame( data = np.random.rand(N, 3), columns = ["a_i","a_j","a_k"] );
#    df = pd.DataFrame( data = np.random.uniform(low=min_val, high=max_val, size=(N, 3)),\
#                      columns = ["a_i","a_j","a_k"] );
    df["a_j"] = 0.3 + df["a_i"] *( 1 + np.random.random(N,)*df["a_i"]);
    df["a_k"] = df["a_k"] ; 
    df["a_i"] = df["a_i"] + 0.6;
    df["a_j"] = df["a_j"] + 0.6;
    
    fig, ax = plt.subplots()
    df.plot( ax = ax);
    ax.set_yticklabels(['{:3.1f}'.format(x*100) \
                                 for x in ax.get_yticks()]);
    ax.legend( df.columns, loc='upper center', bbox_to_anchor=(0.5, 1.00),\
                       ncol=3, fancybox=True, shadow=True);                   
    plt.savefig(params["path_save"]+\
                "profilDeConsommationSeriesTemporelles.eps",\
                dpi= 190)   
    
#---------------------- profil consommation   =  fin --------------------------


if __name__ == '__main__':
    
    start= time.time();
    
    motif_p_s = "p"; vars_k = ["moy_dc","moy_dh"];
    correction = "aleatoire"; priorite = "unitaire"; p_value = 0.5;
    fichier_prefix = "distribution_moyDistLine_moyHamming_k_"; ext = ".txt"
    k_errors = [2,5,10,20]; num_bins = 10#50,100;
    path_save = "/home/willy/Documents/python_topology_learning_simulation/dataArticle/lineaire_simul50Graphes_priorite_unitaire/"; 
    params = {"motif_p_s":motif_p_s, "correction": correction,\
              "priorite":priorite,"p_value": p_value,\
              "vars_k":vars_k,"k_errors":k_errors,\
              "fichier_prefix":fichier_prefix, "ext": ext,\
              "path_save":path_save, "num_bins":num_bins,\
              }
    
    bool_distrib_horiz = True #False;
    bool_profil_consom = False;
    
    if bool_distrib_horiz :
#        histo_seaborn_horizontale(params)
#        distributions_horizontale(params);
#        distribution_horizontale_new(params)
        distribution_horizontale_par_blocs(params)
    if bool_profil_consom :
        params["path_save"] = "/home/willy/Documents/latexDoc/redaction/presentationSoutenanceThese/slide_v1/images/"
        profil_consommation(params)
    