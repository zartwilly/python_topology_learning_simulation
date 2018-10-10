#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  2 10:35:19 2018

@author: willy
test codes
"""

import time, os, re;
import numpy as np;
import pandas as pd;
from pathlib import Path;
import random;
import math;
import logging;
import collections as Coll;
import fonctions_auxiliaires as fct_aux
import matplotlib.pyplot as plt;
import matplotlib.font_manager as font_manager
import seaborn as sns;
import itertools as it;
from functools import reduce;
from operator import mul;
from scipy.stats import norm, truncnorm
from scipy import stats;

#################################
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
##################################

def create_dataframe_bon_seuils(rep, rep_p_correl_s, p_correl_s, motif_p_s, correction, priorisation, args):
    k_errors = lire_k_errors(rep_p_correl_s); 
    df = pd.DataFrame();
    for k_error in k_errors:
        df_err = pd.DataFrame();
        if motif_p_s == "s":
            df_err = pd.read_csv(rep_p_correl_s+args["fichier_prefix"]+str(k_error)+args["ext"],\
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
            df_err = pd.read_csv(rep_p_correl_s+args["fichier_prefix"]+str(k_error)+args["ext"],
                 names=["cpt","moy_dl_"+motif_p_s+"_"+str(k_error)+"_"+correction+"_"+priorisation,\
                        "moy_dh_"+motif_p_s+"_"+str(k_error)+"_"+correction+"_"+priorisation,\
                        "nbre_aretes_matE_"+motif_p_s+"_"+str(k_error)+"_"+correction+"_"+priorisation,\
                        "correl_dh_dl_"+motif_p_s+"_"+str(k_error)+"_"+correction+"_"+priorisation],\
                 sep=";")
    
        df_err.loc[:,"cpt"] = df_err.loc[:, "cpt"].apply(lambda x: int(x.split("_")[2]))
        df_err.sort_values("cpt", inplace=True);
        if df.empty:
            df = pd.concat([df,df_err], axis=1);
        else:
            ## a effacer
            cols_comm = set(df.columns).intersection(set(df_err.columns))
            print("k_error={},df_col={},df_err={}, cols_comm={},df_shape={}"\
                  .format(k_error,len(df.columns), len(df_err.columns), cols_comm,df.shape));
#            print("cpt df={},cpt df_err={}".format(set(df["cpt"]), set(df_err["cpt"])))      
            ## a effacer
            df = pd.merge( df, df_err, on="cpt", how="outer");
#            df = pd.merge( df, df_err, left_on="cpt",right_on="cpt");
#            df = pd.concat([df,df_err], axis=1);
            del(df_err)
    return df.astype(np.int32)



def plot_recherche_bon_seuils_sur_data_brutes(df, args):
    fig, axarr = plt.subplots(3,int(len(args["selected_vars"])/2));
    default_size = fig.get_size_inches();
    fig.set_figheight(default_size[0]*2.0); fig.set_figwidth(default_size[0]*2.0) 
    variabs = "_".join(args["selected_vars"]);
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
                title = "comparaison des \n"+ str(var.upper()) +\
                        "\n en fonction  du seuil");
            filter_cols = ["_".join(("seuil",item.split("_")[3])) \
                            for item in filter_cols]
        elif args["langue"] == "anglais":
            axarr[row_fig,col_fig].set(xlabel= "graphs", ylabel= var.upper(), \
                title = "comparison of \n"+ str(var.upper()) +"\n by thresholds");
            filter_cols = ["_".join((item.split("_")[3]), "threshold") for item in filter_cols]
        axarr[row_fig,col_fig].legend( filter_cols, loc='upper center', bbox_to_anchor=(0.5, 1.00),\
                   ncol=1, fancybox=True, shadow=True);
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
        args["priorisation"] = distrib[5];
        args["correction"]  = distrib[4];
        print("rep={}, rep_p_correl_s={}, p_correl_s={}, motif_p_s={}, \
              correction={}, priorisation={}".format(distrib[0], distrib[1], \
              distrib[2],distrib[3], distrib[4], distrib[5],))
        return df_tmp                                   
        plot_recherche_bon_seuils_sur_data_brutes(df_tmp, args)
#        return 
        
    print("distribs cpt = {}".format(cpt))
    return df; 

    
def create_dataframe_bon_seuils_test(rep, rep_p_correl_s, p_correl_s, \
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
#            print("cpt={}, df_err_cpt.shape={}".format(cpt,df_err_cpt.shape))
            df_err = pd.concat([df_err, df_err_cpt])
        del(df_err_bis)
        print("df_err.shape={}".format(df_err.shape))
        df_err.sort_values("cpt", inplace=True);
        if df.empty:
            df = pd.concat([df,df_err], axis=1);
#            return df.astype(np.float32)
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
#            df = pd.merge( df, df_err, left_on="cpt",right_on="cpt");
#            df = pd.concat([df,df_err], axis=1);
            print("merge df, df_err .shape={}".format(df.shape))
            del(df_err)
    return df.astype(np.float32)


#------------------------------------------------------------------------------
#                   graphe cellule G_nn 
#                       nn = n*n est le nombre de cellules 
#------------------------------------------------------------------------------
import itertools as it;
def aretes(matE, k0_1):
    """
    but: trouver toutes les non aretes de matE cad matE.loc[x][y] == 0
    k0_1 = {0,1}
    """
    liste_cols = matE.columns.tolist()
    tmp = list(); res_aretes = list(); dico_proba_cases = dict()
    for row in liste_cols:
        for col in liste_cols:
            if row != col and (col,row) not in tmp and matE.loc[row][col] == k0_1:
                tmp.append( (col,row) ); tmp.append( (row,col) );
                res_aretes.append( (row, col) );
                dico_proba_cases[(row, col)] = 0.5                              # valeur constante car l'exposant  = 0 dans correction_linegraph.py (fonction de cout)
    return res_aretes, dico_proba_cases;
def matrice_G_nn(N,M):
    """ 
    creer le dataframe M_G_nn avec pour colonnes et index les sommets de G_nn
    """
    # liste de sommets de G_nn
    sommets = list(it.product(range(N), range(M), repeat = 1));
    # creation dataframe avec colonnes = liste de sommets G_nn
    M_G_nn = pd.DataFrame(0, columns = sommets, index = sommets);
    return M_G_nn;

def G_nn_debug(N, M):
    # creation Matrice M_G_nn
    M_G_nn = matrice_G_nn(N,M)
    # affectation valeurs dans M_G_nn
    row, col, cpt = 0, 0, 0;
    dico_sommets = dict(); aretes_G_nn_row_col = list();
    for sommet in M_G_nn.columns:
        row = sommet[0]; col = sommet[1];
        cpt += 1;dico_sommets[(row,col)] = str(cpt);
        if col+1 < M:
#            print("col+1={}".format(col+1))
            M_G_nn.at[(row,col),(row,col+1)] = 1;
            M_G_nn.at[(row,col+1),(row,col)] = 1;
            aretes_G_nn_row_col.append((row,col))
        if col-1 >= 0:
#            print("col-1={}".format(col-1))
            M_G_nn.at[(row,col),(row,col-1)] = 1;
            M_G_nn.at[(row,col-1),(row,col)] = 1;
            aretes_G_nn_row_col.append((row,col))
        if row+1 < N:
#            print("row+1={}".format(row+1))
            M_G_nn.at[(row,col),(row+1,col)] = 1;
            M_G_nn.at[(row+1,col),(row,col)] = 1;
            aretes_G_nn_row_col.append((row,col))
        if row-1 >= 0:
#            print("row-1={}".format(row-1))
            M_G_nn.at[(row,col),(row-1,col)] = 1;
            M_G_nn.at[(row-1,col),(row,col)] = 1;
            aretes_G_nn_row_col.append((row,col))
    # renommer index et colonnes
    M_G_nn.rename(columns = dico_sommets, inplace = True);
    M_G_nn.rename(index = dico_sommets, inplace = True);
    aretes_G_nn = list();
    aretes_G_nn, dico_proba_cases = aretes(M_G_nn, k0_1=1);                                        #aretes_G_nn= aretes obtenus a partir du renommage des sommets (row, col)
#    return M_G_nn, dico_sommets, aretes_G_nn, dico_proba_cases;
    print("fini")
    return M_G_nn;
    
#def G_nn(N, M):
#    """ definir une matrice de N*M cases 
#        idealement N = M
#        row = {0,...,N-1}
#        col = {0,...,M-1}
#    """ 
#    # creation Matrice M_G_nn
#    M_G_nn = matrice_G_nn(N,M)
#    # affectation valeurs dans M_G_nn
#    row, col, cpt = 0, 0,0;
#    dico_sommets = dict();
#    for row in range(N):
#        for col in range(M):
#            cpt += 1;
#            dico_sommets[(row,col)] = str(cpt);
#            if row+1 < N:
#                M_G_nn.at[(row,col),(row+1,col)] = 1;
#                M_G_nn.at[(row+1,col),(row,col)] = 1;
#            if col+1 < M:
#                M_G_nn.at[(row,col),(row,col+1)] = 1;
#                M_G_nn.at[(row,col+1),(row,col)] = 1;
#    # renommer index et colonnes
##    M_G_nn.rename(columns = dico_sommets, inplace = True);
##    M_G_nn.rename(index = dico_sommets, inplace = True);
#    return M_G_nn, dico_sommets;
    
def G_nn(N, M):
    """ definir une matrice de N*M cases 
        idealement N = M
        row = {0,...,N-1}
        col = {0,...,M-1}
    """ 
    # creation Matrice M_G_nn
    M_G_nn = matrice_G_nn(N,M)
    # affectation valeurs dans M_G_nn
    row, col, cpt = 0, 0, 0;
    dico_sommets = dict(); 
    aretes_G_nn_row_col = list();                                                # aretes_G_nn_row_col = liste aretes en fonction des row et cols; 
    for row in range(N):
        for col in range(M):
            cpt += 1;
            dico_sommets[(row,col)] = str(cpt);
            if row+1 < N:
                M_G_nn.at[(row,col),(row+1,col)] = 1;
                M_G_nn.at[(row+1,col),(row,col)] = 1;
                aretes_G_nn_row_col.append((row,col))
            if col+1 < M:
                M_G_nn.at[(row,col),(row,col+1)] = 1;
                M_G_nn.at[(row,col+1),(row,col)] = 1;
                aretes_G_nn_row_col.append((row,col))
    # renommer index et colonnes
    M_G_nn.rename(columns = dico_sommets, inplace = True);
    M_G_nn.rename(index = dico_sommets, inplace = True);
    aretes_G_nn = list();
    aretes_G_nn, dico_proba_cases = aretes(M_G_nn, k0_1=1);                                        #aretes_G_nn= aretes obtenus a partir du renommage des sommets (row, col)
    return M_G_nn
    return M_G_nn, dico_sommets, aretes_G_nn, dico_proba_cases;

def compare_df(m_g_nn, m_g_nn_dbg):
#    united_data = pd.concat([m_g_nn, m_g_nn_dbg])
#    # group the data by the whole row to find duplicates
#    united_data_grouped = united_data.groupby(list(united_data.columns))
#    # detect the row indices of unique rows
#    uniq_data_idx = [x[0] for x in united_data_grouped.indices.values() if len(x) == 1]
#    # extract those unique values
#    uniq_data = united_data.iloc[uniq_data_idx] 
#    print("uniq_data = {}".format(uniq_data))
    
    for row in m_g_nn.columns:
        for col in m_g_nn.columns:
#            print("***row={}, col ={}".format(row, col))
            if m_g_nn.at[row,col] != m_g_nn_dbg.at[row,col]:
                print("row={}, col ={}, val_g_nn={}, val_g_nn_dbg={}"\
                      .format(row, col, m_g_nn.at[row,col], m_g_nn_dbg.at[row,col]))
                

#------------------------------------------------------------------------------        
#--------------------- new test creation df  graphe cellule--------------------    
#           creer un liste d'adjacence puis
#           une transformer la liste d'adjacence en une matrice d'adjacence
#------------------------------------------------------------------------------
def liste_adjacence_G_nn(N,M):
    sommets = list(it.product(range(N), range(M), repeat = 1));
    dico_nom_sommet = dict(); cpt = 0; dico_adjacence=dict();
    aretes_G_nn_row_col = set();
    for sommet in it.product(range(N), range(M), repeat = 1):
        row = sommet[0]; col = sommet[1];
        cpt += 1; dico_nom_sommet[(row,col)] = cpt;
        dico_adjacence[(row,col)] = [];
        if col+1 < M:
            dico_adjacence[(row,col)].append((row,col+1));
            aretes_G_nn_row_col.update([(row,col),(row,col+1)])
        if col-1 >= 0:
            dico_adjacence[(row,col)].append((row,col-1));
            aretes_G_nn_row_col.update([(row,col),(row,col-1)])
        if row+1 < N:
            dico_adjacence[(row,col)].append((row+1,col));
            aretes_G_nn_row_col.update([(row,col),(row+1,col)])
        if row-1 >= 0:
            dico_adjacence[(row,col)].append((row-1,col));
            aretes_G_nn_row_col.update([(row,col),(row-1,col)])
        if row == 0 and col == 0:
            dico_adjacence[(row,col)].append((0,M-1));                         # M_G_nn.at[(0,0),(0,M-1)] = 1
            dico_adjacence[(row,col)].append((N-1,0));                         # M_G_nn.at[(N-1,0),(0,0)] = 1;
            aretes_G_nn_row_col.update([(row,col),(0,M-1),(N-1,0)])
        if row == N-1 and col == 0:
            dico_adjacence[(row,col)].append((N-1,M-1));                       # M_G_nn.at[(N-1,0),(N-1,M-1)] = 1; 
            dico_adjacence[(row,col)].append((0,0));                           # M_G_nn.at[(N-1,0),(0,0)] = 1;
            aretes_G_nn_row_col.update([(row,col),(N-1,M-1),(0,0)])
        if row == 0 and col == M-1:
            dico_adjacence[(row,col)].append((N-1,M-1));                       # M_G_nn.at[(0,M-1),(N-1,M-1)] = 1; 
            dico_adjacence[(row,col)].append((0,0));                           # M_G_nn.at[(0,M-1),(0,0)] = 1;
            aretes_G_nn_row_col.update([(row,col),(N-1,M-1),(N-1,0)])
        if row == N-1 and col == M-1:
            dico_adjacence[(row,col)].append((0,M-1));                         # M_G_nn.at[(N-1,M-1),(0,M-1)] = 1; 
            dico_adjacence[(row,col)].append((N-1,0));                         # M_G_nn.at[(N-1,M-1),(N-1,0)] = 1;
            aretes_G_nn_row_col.update([(row,col),(0,M-1),(N-1,0)])
    return dico_adjacence;
        
        
#     M_G_nn.at[(0,0),(0,M-1)] = 1 ; M_G_nn.at[(0,M-1),(0,0)] = 1
#     M_G_nn.at[(N-1,0),(0,0)] = 1 ; M_G_nn.at[(0,0),(N-1,0)] = 1; 
#     M_G_nn.at[(N-1,0),(N-1,M-1)] = 1; M_G_nn.at[(N-1,M-1),(N-1,0)] = 1;
#     M_G_nn.at[(N-1,0),(0,0)] = 1 ; M_G_nn.at[(0,0),(N-1,0)] = 1 
#     M_G_nn.at[(0,M-1),(N-1,M-1)] = 1 ; M_G_nn.at[(N-1,M-1),(0,M-1)] = 1 
#     M_G_nn.at[(0,M-1),(0,0)] = 1 ; M_G_nn.at[(0,0),(0,M-1)] = 1  
#     M_G_nn.at[(N-1,M-1),(0,M-1)] = 1 ; M_G_nn.at[(0,M-1),(N-1,M-1)] = 1;  
#     M_G_nn.at[(N-1,M-1),(N-1,0)] = 1 ; M_G_nn.at[(N-1,0),(N-1,M-1)] = 1;
        
    pass

def G_nn_new_avec_Cpt(N, M):
    sommets = list(it.product(range(N), range(M), repeat = 1));
    # creation dataframe avec colonnes = liste de sommets G_nn
    M_G_nn = pd.DataFrame(0, columns = sommets, index = sommets);
    
    dico_sommets = dict(); cpt = 0; 
    aretes_G_nn_row_col = set(); aretes_G_nn = set(); dico_proba_cases = dict();
    
    for sommet in sommets:
        row = sommet[0]; col = sommet[1];
        cpt += 1; dico_sommets[(row,col)] = str(cpt); 
#        dico_proba_cases[(row, col)] = 0.5; 
        if col+1 < M:
            M_G_nn.at[(row,col),(row,col+1)] = 1;
            M_G_nn.at[(row,col+1),(row,col)] = 1;
            aretes_G_nn_row_col.update([(row,col),(row,col+1)])
            if (row,col+1) not in dico_sommets.keys() and \
                (col+1,row) not in dico_sommets.keys():
                cpt += 1; dico_sommets[(row,col+1)] = str(cpt);
            elif (row,col+1) in dico_sommets.keys():
                aretes_G_nn.add((dico_sommets[(row,col)],dico_sommets[(row,col+1)]))
            elif (col+1,row) in dico_sommets.keys():
                aretes_G_nn.add((dico_sommets[(row,col)],dico_sommets[(col+1,row)]))
#            if (dico_sommets[(row,col+1)], dico_sommets[(row,col)]) not in aretes_G_nn:
#            aretes_G_nn.add((dico_sommets[(row,col)],dico_sommets[(row,col+1)]))
        if col-1 >= 0:
            M_G_nn.at[(row,col),(row,col-1)] = 1;
            M_G_nn.at[(row,col-1),(row,col)] = 1;
            aretes_G_nn_row_col.update([(row,col),(row,col-1)])
            if (row,col-1) not in dico_sommets.keys() and \
                (col-1,row) not in dico_sommets.keys() :
                cpt += 1; dico_sommets[(row,col-1)] = str(cpt);
            elif (row,col-1) in dico_sommets.keys() :
                aretes_G_nn.add((dico_sommets[(row,col)],dico_sommets[(row,col-1)]))
            elif (col-1,row) in dico_sommets.keys() :
                aretes_G_nn.add((dico_sommets[(row,col)],dico_sommets[(col-1,row)]))
#            if (dico_sommets[(row,col-1)],dico_sommets[(row,col)]) not in aretes_G_nn:
#            aretes_G_nn.add((dico_sommets[(row,col)],dico_sommets[(row,col-1)]))
        if row+1 < N:
            M_G_nn.at[(row,col),(row+1,col)] = 1;
            M_G_nn.at[(row+1,col),(row,col)] = 1;
            aretes_G_nn_row_col.update([(row,col),(row+1,col)])
            if (row+1,col) not in dico_sommets.keys() and \
                (col,row+1) not in dico_sommets.keys():
                cpt += 1; dico_sommets[(row+1,col)] = str(cpt);
                aretes_G_nn.add((dico_sommets[(row,col)],dico_sommets[(row+1,col)]))
            elif (row+1,col) in dico_sommets.keys() :
                aretes_G_nn.add((dico_sommets[(row,col)],dico_sommets[(row+1,col)]))
            elif (col,row+1) in dico_sommets.keys() :
                aretes_G_nn.add((dico_sommets[(row,col)],dico_sommets[(col,row+1)]))
        if row-1 >= 0:
            M_G_nn.at[(row,col),(row-1,col)] = 1;
            M_G_nn.at[(row-1,col),(row,col)] = 1;
            aretes_G_nn_row_col.update([(row,col),(row-1,col)])
            if (row-1,col) not in dico_sommets.keys() and \
                (col,row-1) not in dico_sommets.keys():
                cpt += 1; dico_sommets[(row-1,col)] = str(cpt);
                aretes_G_nn.add((dico_sommets[(row,col)],dico_sommets[(row-1,col)]))
            elif (row-1,col) in dico_sommets.keys() :
                aretes_G_nn.add((dico_sommets[(row,col)],dico_sommets[(row-1,col)]))
            elif (col,row-1) in dico_sommets.keys():
                aretes_G_nn.add((dico_sommets[(row,col)],dico_sommets[(col,row-1)]))            
    
    M_G_nn.at[(0,0),(0,M-1)] = 1; M_G_nn.at[(0,M-1),(0,0)] = 1; 
    M_G_nn.at[(N-1,0),(0,0)] = 1; M_G_nn.at[(0,0),(N-1,0)] = 1;
    aretes_G_nn_row_col.update([(0,0),(0,M-1),(N-1,0)])
    if (0,M-1) not in dico_sommets.keys() and (M-1,0) not in dico_sommets.keys():
        cpt += 1; dico_sommets[(0,M-1)] = str(cpt);
    if (N-1,0) not in dico_sommets.keys() and (0, N-1) not in dico_sommets.keys():
        cpt += 1; dico_sommets[(N-1,0)] = str(cpt);
    if (dico_sommets[(0,M-1)],dico_sommets[(0,0)]) not in aretes_G_nn and \
        (dico_sommets[(N-1,0)],dico_sommets[(0,0)]) not in aretes_G_nn:
        aretes_G_nn.update([(dico_sommets[(0,0)],dico_sommets[(0,M-1)]),\
                        (dico_sommets[(0,0)],dico_sommets[(N-1,0)])])
    
    M_G_nn.at[(N-1,0),(N-1,M-1)] = 1; M_G_nn.at[(N-1,M-1),(N-1,0)] = 1; 
    M_G_nn.at[(N-1,0),(0,0)] = 1; M_G_nn.at[(0,0),(N-1,0)] = 1;
    aretes_G_nn_row_col.update([(N-1,0),(N-1,M-1),(0,0)])
    if (N-1,M-1) not in dico_sommets.keys() and \
        (M-1,N-1) not in dico_sommets.keys():
        cpt += 1; dico_sommets[(N-1,M-1)] = str(cpt);
    if (N-1,0) not in dico_sommets.keys() and \
        (0,N-1) not in dico_sommets.keys():
        cpt += 1; dico_sommets[(N-1,0)] = str(cpt);
    if (dico_sommets[(N-1,M-1)],dico_sommets[(N-1,0)]) not in aretes_G_nn and \
        (dico_sommets[(0,0)],dico_sommets[(N-1,0)]) not in aretes_G_nn:
        aretes_G_nn.update([(dico_sommets[(N-1,0)],dico_sommets[(N-1,M-1)]),\
                        (dico_sommets[(N-1,0)],dico_sommets[(0,0)])])

    M_G_nn.at[(0,M-1),(N-1,M-1)] = 1; M_G_nn.at[(N-1,M-1),(0,M-1)] = 1;
    M_G_nn.at[(0,M-1),(0,0)] = 1; M_G_nn.at[(0,0),(0,M-1)] = 1;
    aretes_G_nn_row_col.update([(0,M-1),(N-1,M-1),(N-1,0)])
    if (N-1,M-1) not in dico_sommets.keys() and \
        (M-1,N-1) not in dico_sommets.keys():
        cpt += 1; dico_sommets[(N-1,M-1)] = str(cpt);
    if (0,M-1) not in dico_sommets.keys() and \
        (M-1,0) not in dico_sommets.keys():
        cpt += 1; dico_sommets[(0,M-1)] = str(cpt);
    if (dico_sommets[(N-1,M-1)],dico_sommets[(0,M-1)]) not in aretes_G_nn and \
       (dico_sommets[(0,0)],dico_sommets[(0,M-1)]) not in aretes_G_nn :
        aretes_G_nn.update([(dico_sommets[(0,M-1)],dico_sommets[(N-1,M-1)]),\
                        (dico_sommets[(0,M-1)],dico_sommets[(0,0)])])
    
    M_G_nn.at[(N-1,M-1),(0,M-1)] = 1; M_G_nn.at[(0,M-1),(N-1,M-1)] = 1; 
    M_G_nn.at[(N-1,M-1),(N-1,0)] = 1; M_G_nn.at[(N-1,0),(N-1,M-1)] = 1; 
    aretes_G_nn_row_col.update([(N-1,M-1),(0,M-1),(N-1,0)])
    if (N-1,M-1) not in dico_sommets.keys() and \
        (M-1,N-1) not in dico_sommets.keys():
        cpt += 1; dico_sommets[(N-1,M-1)] = str(cpt);
    if (0,M-1) not in dico_sommets.keys() and \
        (M-1,0) not in dico_sommets.keys():
        cpt += 1; dico_sommets[(0,M-1)] = str(cpt);
    if (N-1,0) not in dico_sommets.keys() and \
        (0,N-1) not in dico_sommets.keys():
        cpt += 1; dico_sommets[(N-1,0)] = str(cpt);
    if (dico_sommets[(M-1,N-1)],dico_sommets[(M-1,0)]) not in aretes_G_nn and \
        (dico_sommets[(N-1,0)],dico_sommets[(N-1,M-1)]) not in aretes_G_nn:
        aretes_G_nn.update([(dico_sommets[(N-1,M-1)],dico_sommets[(0,M-1)]),\
                        (dico_sommets[(N-1,M-1)],dico_sommets[(N-1,0)])])
    
    # renommer index et colonnes
    M_G_nn.rename(columns = dico_sommets, inplace = True);
    M_G_nn.rename(index = dico_sommets, inplace = True);
    
    print("aretes_G_nn")
    # aretes G_nn
    aretes_G_nn = set()
    for row, col in fct_aux.range_2d(M_G_nn):
        if M_G_nn.at[row, col] == 1 and (row,col) not in aretes_G_nn:
            aretes_G_nn.add((row,col));
        dico_proba_cases[(row, col)] = 0.5;
    print("M_G_nn=\n{}\n, dico_sommets={}\n, aretes_G_nn = {}\n, dico_proba_cases={}"\
          .format(M_G_nn, dico_sommets, aretes_G_nn, dico_proba_cases));
    return M_G_nn, dico_sommets, aretes_G_nn, dico_proba_cases;
    pass

def G_nn_new_sans_Cpt(N, M):
    sommets = [ "_".join([str(tu[0]),str(tu[1])]) \
               for tu in it.product(range(N), range(M), repeat = 1)];
    # creation dataframe avec colonnes = liste de sommets G_nn
    M_G_nn = pd.DataFrame(0, columns = sommets, index = sommets);
    
    dico_sommets, dico_proba_cases = dict(), dict();
    aretes_G_nn = set();
    for sommet in sommets:
        row = int(sommet.split("_")[0]); col = int(sommet.split("_")[1]);
        dico_sommets[(row,col)] = sommet;
        dico_proba_cases[sommet] = 0.5; 
        if col+1 < M:
            row_col1 = "_".join([str(row),str(col+1)])                         # row_col1 = (row,col+1)
            M_G_nn.at[sommet,row_col1] = 1;
            M_G_nn.at[row_col1,sommet] = 1;
            dico_sommets[(row,col+1)] = row_col1;
            if (row_col1,sommet) not in aretes_G_nn:
                aretes_G_nn.update([(sommet,row_col1)])
        if col-1 >= 0:
            row_col_1 = "_".join([str(row),str(col-1)])                        # row_col1 = (row,col-1)
            M_G_nn.at[sommet,row_col_1] = 1;
            M_G_nn.at[row_col_1,sommet] = 1;
            dico_sommets[(row,col-1)] = row_col_1;
            if (row_col_1,sommet) not in aretes_G_nn:
                aretes_G_nn.update([(sommet,row_col_1)])            
        if row+1 < N:
            row1_col = "_".join([str(row+1),str(col)])                         # row1_col = (row+1,col)
            M_G_nn.at[sommet,row1_col] = 1;
            M_G_nn.at[row1_col,sommet] = 1;
            dico_sommets[(row+1,col)] = row1_col;
            if (row1_col,sommet) not in aretes_G_nn:
                aretes_G_nn.update([(sommet,row1_col)])
        if row-1 >= 0:
            row_1_col = "_".join([str(row-1),str(col)])                        # row_1_col = (row-1,col)
            M_G_nn.at[sommet,row_1_col] = 1;
            M_G_nn.at[row_1_col,sommet] = 1;
            dico_sommets[(row-1,col)] = row_1_col;
            if (row_1_col,sommet) not in aretes_G_nn:
                aretes_G_nn.update([(sommet,row_1_col)])
    
    row_col_0_0 = "_".join([str(0),str(0)])                                    # row_col_0_0 = (0,0)
    row_col_0_M_1 = "_".join([str(0),str(M-1)])                                # row_col_0_M_1 = (0,M-1)
    row_col_0_N_1 = "_".join([str(0),str(N-1)])                                # row_col_0_N_1 = (0,N-1)
    M_G_nn.at[row_col_0_0,row_col_0_M_1] = 1; 
    M_G_nn.at[row_col_0_M_1,row_col_0_0] = 1; 
    M_G_nn.at[row_col_0_N_1,row_col_0_0] = 1; 
    M_G_nn.at[row_col_0_0,row_col_0_N_1] = 1;
    if (row_col_0_M_1,row_col_0_0) not in aretes_G_nn:
        aretes_G_nn.update([ (row_col_0_0, row_col_0_M_1)]);
    if (row_col_0_N_1,row_col_0_0) not in aretes_G_nn:
        aretes_G_nn.update([(row_col_0_0,row_col_0_N_1)])
    
    row_col_N_1_0 = "_".join([str(N-1),str(0)])                                # row_col_N_1_0 = (N-1,0)
    row_col_N_1_M_1 = "_".join([str(N-1),str(M-1)])                            # row_col_N_1_M_1 = (N-1,M-1)
    row_col_N_1_0 = "_".join([str(N-1),str(0)])                                # row_col_N_1_0 = (N-1,0)
    M_G_nn.at[row_col_N_1_0,row_col_N_1_M_1] = 1; 
    M_G_nn.at[row_col_N_1_M_1,row_col_N_1_0] = 1; 
    M_G_nn.at[row_col_N_1_0,row_col_0_0] = 1; 
    M_G_nn.at[row_col_0_0,row_col_N_1_0] = 1;
    if (row_col_N_1_M_1,row_col_N_1_0) not in aretes_G_nn:
        aretes_G_nn.update([(row_col_N_1_0,row_col_N_1_M_1)])
    if (row_col_0_0,row_col_N_1_0) not in aretes_G_nn:
        aretes_G_nn.update([(row_col_N_1_0,row_col_0_0)])
    
    M_G_nn.at[row_col_0_M_1,row_col_N_1_M_1] = 1; 
    M_G_nn.at[row_col_N_1_M_1,row_col_0_M_1] = 1;
    M_G_nn.at[row_col_0_M_1,row_col_0_0] = 1; 
    M_G_nn.at[row_col_0_0,row_col_0_M_1] = 1;
    if (row_col_N_1_M_1,row_col_0_M_1) not in aretes_G_nn:
        aretes_G_nn.update([(row_col_0_M_1,row_col_N_1_M_1)])
    if (row_col_0_M_1,row_col_0_0) not in aretes_G_nn:
        aretes_G_nn.update([(row_col_0_0,row_col_0_M_1)])
        
    M_G_nn.at[row_col_N_1_M_1, row_col_0_M_1] = 1; 
    M_G_nn.at[row_col_0_M_1, row_col_N_1_M_1] = 1; 
    M_G_nn.at[row_col_N_1_M_1, row_col_0_N_1] = 1; 
    M_G_nn.at[row_col_N_1_0, row_col_N_1_M_1] = 1; 
    if (row_col_0_M_1, row_col_N_1_M_1) not in aretes_G_nn:
        aretes_G_nn.update([(row_col_N_1_M_1,row_col_0_M_1)])
    if (row_col_N_1_0, row_col_N_1_M_1) not in aretes_G_nn:
        aretes_G_nn.update([(row_col_N_1_M_1, row_col_0_N_1)])
        
    print("aretes_G_nn")
#    print("M_G_nn=\n{}\n, dico_sommets={}:{}\n, aretes_G_nn ={}: {}\n, \
#          dico_proba_cases={}:{}"\
#          .format(M_G_nn, dico_sommets, len(dico_sommets), aretes_G_nn, \
#                  len(aretes_G_nn), len(dico_proba_cases), dico_proba_cases));
    print("M_G_nn=\n{}\n, dico_sommets={}\n, aretes_G_nn={}:{}\n, \
          dico_proba_cases={}"\
          .format(M_G_nn.shape, len(dico_sommets), \
                  len(aretes_G_nn), aretes_G_nn,len(dico_proba_cases)));
    # networkx
#    G = nx.Graph(M_G_nn.values);                                              # plot graph with networkX
#    nx.draw(G, pos=nx.spring_layout(G), with_labels=True);
#    print("aretes_G_nn={},M_G_nn={},dico_sommets={} ".format(aretes_G_nn, M_G_nn.shape,dico_sommets ))
    return M_G_nn, dico_sommets, aretes_G_nn, dico_proba_cases;
    pass      

#------------------------------------------------------------------------------
#--------------------- new test creation df  graphe cellule-------------    
#------------------------------------------------------------------------------        

if __name__ == '__main__':
    
    start= time.time();
    ### definition Matrice
    N,M= 10,10; N,M= 30,30; N,M= 100,100; N, M = 20,20 ;N,M = 4,4 ;
#    m_g_nn = G_nn(N, M);
    print("fin g_nn = {}".format( time.time() - start))
    start= time.time();
#    m_g_nn_dbg = G_nn_debug(N, M);
    print("fin g_nn_dbg = {}".format( time.time() - start))
#    compare_df(m_g_nn, m_g_nn_dbg)
    
    start= time.time();
    m_g_nn_new = G_nn_new_sans_Cpt(N, M);
    print("fin g_nn_new_sans_cpt = {}".format( time.time() - start))
    
#    compare_df(m_g_nn_new[0], m_g_nn_dbg)
    
    
    
    langue = "francais";
    bool_p_correl = False; bool_seuil = None; 
    bool_k_errors=False; distrib_name = None; ext = ".txt"; seleted_vars = "";
    rep = "data/"; 
    if bool_p_correl:
        rep = "data/"; bool_seuil = False; 
        distrib_name = "distribution_moyDistLine_moyHamming_k_";
        selected_vars = ["moy_dh"];
        bool_k_errors=False; k_errors = [1,2,5,9];
#        bool_k_errors=True; k_errors = [1,2,3,4,5,6,7,8,9];
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
    booleen_rep = False; priorisation_ = "aucune"; correction_ = "aleatoire";
    priorisations = ["lineaire_simul50Graphes_priorite_supp", \
                      "lineaire_simul50Graphes_priorite_ajout", \
                      "lineaire_simul50Graphes_priorite_aucune"];
    corrections =  ["degreMin","coutMin","aleatoire"];
    reps_ = [ rep+priorisation+"/"+correction+"/" \
            for priorisation in priorisations for correction in corrections];
    reps = reps_.copy();
    if booleen_rep == True :
        reps = [rep for rep in reps_ if rep.find(priorisation_)>0 and rep.find(correction_)>0];
        
    args={"bool_p_correl":bool_p_correl,"bool_seuil":bool_seuil,"ext":ext, \
          "distrib_name":distrib_name,"bool_k_errors":bool_k_errors,\
          "k_errors":k_errors, "selected_vars":selected_vars, \
          "langue":langue, "rep":rep};
    
    bool_bon_choix = True;    
    selected_p_correl_s = 0.5; args["selected_p_correl_s"] = selected_p_correl_s;

    if bool_bon_choix:
#        df_choix = formation_dataframe_bon_seuils_sur_data_brutes(reps,args);
#        plot_recherche_bon_seuils_sur_data_brutes(df_choix,args)
#        plot_seuils(df_fn_fp, args);
        pass


###### test merge df_a df_b by k_errors= [0.1, 0.2]
    rep="simulation_seuils/lineaire_simul50Graphes_priorite_supp/degreMin/"; 
    rep_p_correl_s="simulation_seuils/lineaire_simul50Graphes_priorite_supp/degreMin/distribution/";
    p_correl_s=0.6; motif_p_s="s"; correction="degreMin"; priorisation = "supp"
    args["motif_p_s"] = "s";
    args["fichier_prefix"] = "distribution_moyDistLine_moyHamming_s_";
#    df_choix = create_dataframe_bon_seuils_test(rep, rep_p_correl_s, p_correl_s, \
#                                     motif_p_s, correction, priorisation, args);
###### test merge df_a df_b by k_errors= [0.1, 0.2]
