#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 19 19:24:51 2018

@author: willy
"""


import math;
import os, re;
import numpy as np;
import pandas as pd;
import fonctions_auxiliaires as fct_aux;
import ajout_mesures_correlations as ajusterCorrelation;
import creation_datasets as create_datasets;
import scipy as spy;
from scipy.stats import linregress;
from pathlib import Path;
import time;
import logging;
import collections;
import extractGraphReel as graphReel;
import sax_encoding as sax_encod;
import ujson;
from functools import reduce;
from operator import mul;
import shapelets_transformV1 as shapelet;
from ast import literal_eval

#############################################################################
###     correlation selon la metrique                               #########
###     X = trace( (A-A_bar).transpose * (B-B_bar)
###    r_ab = \frac{X}{square_root(X) * square_root(X) }
###                                                                 #########
####################### === debut === #######################################
def lire_matrice_equipements(chemin_mat_equipement):
    """
    return la liste des equipements, se trouvant dans le repertoire "chemin_mat_equipement"
    """
    files = os.listdir(chemin_mat_equipement)
    liste_files = [fichier for fichier in files if  re.match("^matrice_.*",fichier)]
 
    pattern = re.compile('matrice_|.csv')        
    equipements = [pattern.sub("", mot_file) for mot_file in liste_files]
    return equipements;
    
def metrique_slicing(correls_pearson):
    """
    creer un dico dont les cles sont (0.0,0.1),(0.1,0.2),...,(0.9,1),(1,)
    et les valeurs sont les correlations dans cet intervalle.
    """
    correls_pearson = [0 if np.isnan(x) else x for x in correls_pearson]
    hist, bin_edges = np.histogram(correls_pearson, bins = 10)
    if len(set(correls_pearson)) == 1:
        min_correl = correls_pearson[0]
    elif len(set(hist)) == 1:
        min_correl = bin_edges[0]
    else:
        max_ind = max(enumerate(hist), key = lambda x: x[1])[0]
        min_correl = bin_edges[max_ind+1]
    return min_correl


def conjugee_matrice(matA):
    """
    calculer le conjuge de matA (matA barre)
    appliquer mean a chaque df column (apply mean to each df column)
    """
    matA_bar = pd.DataFrame();
    for col in matA.columns.tolist():
        mean_df_col = matA[col].mean()
        matA_bar[col] = matA[col].apply(lambda value: value - mean_df_col)
#        print("col",col,"mean = ", mean_df_col, " mean matA_bar", matA_bar[col].mean())
    return matA_bar;

def diff_A_A_bar(matA, matA_bar):
    """
    A - A_bar
    """
    return matA - matA_bar

def metriqueFusion_correlation(matA, matB):
    """
    X = trace( (A-A_bar).transpose * (B-B_bar)
    r_ab = \frac{X}{square_root(X) * square_root(X) }
    """
#    print("matA = ", matA.describe())
    matA_bar = conjugee_matrice(matA);
    matB_bar = conjugee_matrice(matB);
#    print("matA {} matA_bar {} matB {} matB_bar {}".format(matA.isnull().all().any(),\
#           matA_bar.isnull().all().any() \
#          ,matB.isnull().all().any(),matB_bar.isnull().all().any()) )
    mat_A_B_res = np.matrix( (matA - diff_A_A_bar(matA,matA_bar)).transpose()) \
                            * np.matrix((matB - diff_A_A_bar(matB,matB_bar)))
#    print("mat_A_B_res = ", pd.DataFrame(mat_A_B_res).isnull().all().any())
#    print("mat_A_B_res = {}".format(mat_A_B_res))
    trace_mat_A_B = np.trace(mat_A_B_res);
#    print("trace_mat_A_B =", trace_mat_A_B)
    
    mat_A_A_res = np.matrix( (matA - diff_A_A_bar(matA,matA_bar)).transpose()) \
                            * np.matrix((matA - diff_A_A_bar(matA,matA_bar)))
    trace_mat_A_A_bar = np.trace(mat_A_A_res)
    
    mat_B_B_res = np.matrix( (matB - diff_A_A_bar(matB,matB_bar)).transpose()) \
                    * np.matrix( (matB - diff_A_A_bar(matB,matB_bar)) )
    trace_mat_B_B_bar = np.trace(mat_B_B_res)
    print("trace_mat_A_B={}, trace_mat_B_B_bar={}, trace_mat_A_A_bar={}"\
          .format(trace_mat_A_B,trace_mat_B_B_bar, trace_mat_A_A_bar))
#    print("mat_A_B_res={},\n mat_A_A_res={},\n mat_B_B_res={}"\
#           .format(mat_A_B_res,mat_A_A_res,mat_B_B_res))
    r_ab = trace_mat_A_B / (math.sqrt(trace_mat_A_A_bar) * math.sqrt(trace_mat_B_B_bar))
    return r_ab
    
def ts_slicing(ts1, ts2, n, mode):
    """
    permet de recuperer les n elements des 2 times series
    n = la fenetre de glissements
    mode = mode de correlation exple: (correlationParMorceaux ou correlationGlissante)
    NB: * ts1 et ts2 ont la meme taille
        * si ts1 et ts2 n'ont pas les memes tailles alors
            - methode1: on complete le plus petit ts avec des nan ==> agrandir index 1
            - methode2: on reduit le plus grand ts a la taille du plus petit ts ===> reduire index 2
    """
    size_ts = max(ts1.shape[0], ts2.shape[0]);
    if mode == "correlationParMorceaux":
        for cpt_init in range(0, size_ts, n):
            ts1_f = ts1.iloc[cpt_init: cpt_init+n];
            ts2_f = ts2.iloc[cpt_init: cpt_init+n];
            if ts1_f.shape[0] - ts2_f.shape[0] > 0:
                print("ici1")
                print("ts1_a = ",  len(ts1_f));print("ts_2 = ",  len(ts2_f))
                ###
                ts2_f = ts2_f.reindex(index=range(ts1_f.shape[0])) # ==> BON 1
#                n = ts1_f.shape[0] - ts2_f.shape[0]               # ==> 2
#                ts1_f.drop( ts1_f.tail(n).index,inplace=True)     # ==> 2
                ###
                pass
            elif ts1_f.shape[0] - ts2_f.shape[0] < 0:
                print("ici2")
                ###
                ts1_f = ts1_f.reindex(index=range(ts2_f.shape[0])) # ==> BON
#                n = ts2_f.shape[0] - ts1_f.shape[0]
#                ts2_f.drop( ts2_f.tail(n).index,inplace=True)
                ###
                pass
            ts1_f.fillna(method="pad", inplace = True); 
            ts2_f.fillna(method="pad", inplace=True);
            yield ts1_f, ts2_f;
    elif mode == "correlationGlissante":
        for cpt_init in range(0, size_ts-n):
            ts1_f = ts1.iloc[cpt_init: cpt_init+n];
            ts2_f = ts2.iloc[cpt_init: cpt_init+n];
            if ts1_f.shape[0] - ts2_f.shape[0] > 0:
                ts2_f = ts2_f.reindex(index=range(ts1_f.shape[0]))
            elif ts1_f.shape[0] - ts2_f.shape[0] < 0:
                ts2_f = ts2_f.reindex(index=range(ts1_f.shape[0]))
            ts1_f.fillna(method="pad", inplace = True); 
            ts2_f.fillna(method="pad", inplace=True);
            yield ts1_f, ts2_f;
    else:
        print("ERROR : methode unknown".upper())


def normalize(df):
    result = df.copy()
#    print("shape result = {}".format(result.shape) )
    for feature_name in df.columns.tolist():
        mean = df[feature_name].mean(skipna=True);
        std = df[feature_name].std(skipna=True);
        if df[feature_name].isnull().all() != True and std !=0 and \
            np.isnan(std) and np.isnan(mean):
            result[feature_name] = (df[feature_name] - mean) / (std)
#    print("result %s ", result.isnull().all().any())
#   result = (df - df.mean(skipna=True))/df.std(skipna=True) 
    return result;

def correlation_pearson(ts1, ts2, fenetre):
    """
    calculer le pearson par fenetre glissante
    ts1 et ts2 ont la meme taille
    """
    liste_pearson = list();
#    ts1 = pd.Series(ts1); ts2 = pd.Series(ts2);
    size_ts = ts1.count() - fenetre +1;
    cpt_ini = 0
    while size_ts > 0:
        ts1_f = ts1.iloc[cpt_ini:cpt_ini+fenetre];
        ts2_f = ts2.iloc[cpt_ini:cpt_ini+fenetre];
        ts1_f.fillna(method="pad", inplace = True);
        ts2_f.fillna(method="pad", inplace = True);
        valeur_corr = linregress(ts1_f,ts2_f).rvalue;
        liste_pearson.append(valeur_corr);
        cpt_ini += 1;
        size_ts -= 1; 
    result = metrique_slicing(liste_pearson)
    return result;    
    
def range_2d_ts(ts1, ts2):
    min_len = min(len(ts1),len(ts2))
    for i in range(0, min_len):
        yield ts1[i], ts2[i];
        
def distance_pearson(ts1,ts2):
    """
    utilisons la relation existante entre distance euclidienne et la coeff de correlation
    """
    """  """
    distances = list();
    val_na = np.inf #-1200 # -1000
    where_are_NaNs = np.isnan(ts1); ts1[where_are_NaNs] = val_na;
    where_are_NaNs = np.isnan(ts2); ts2[where_are_NaNs] = val_na;
    ts1 = [x for x in ts1 if x != val_na ]
    ts2 = [x for x in ts2 if x != val_na ]
    for ts1_a, ts2_b in range_2d_ts(ts1, ts2):
        distances.append(lpNorm(ts1_a, ts2_b,2));
        
    n = min(len(ts1),len(ts2)); correl = 0;
    correl = 1 - 1/(2*n) * sum(distances) if n != 0 else 0;
#    print("ts1={}, ts2={}, n={}, sum={}, correl={} "\
#           .format(ts1_var, ts2_var,n,sum(distances),correl))
    return correl
    pass
def metrique_wil_histo(ts1,ts2):
    """
    BAD
    faire la lpNorm entre x(t) et y(t) for t in [1,dim(ts)]
    difference entre les lpNorm(i) et lpNorm(i+1) => e_i = lpNorm_i - lpNorm_i+1
    distribution des E={e_i}
    recuperer la distribution maximale x 
    diviser x par len(E)
    """
    val_na = np.nan # -1000
    where_are_NaNs = np.isnan(ts1); ts1[where_are_NaNs] = val_na;
    where_are_NaNs = np.isnan(ts2); ts2[where_are_NaNs] = val_na;
    ts1 = [x for x in ts1 if x != val_na ]
    ts2 = [x for x in ts2 if x != val_na ]
    E = list()
    for t1_i, t2_i in range_2d_ts(ts1, ts2):
        E.append( round(lpNorm(t1_i,t2_i,1),2) );
    max_E = max(E); min_E = min(E);
    bins = np.linspace(min_E, max_E,11);
    val_batons, intervals = np.histogram(E, bins); 
#    print("val_batons ={}, E={}, bins={}".format(val_batons, len(E),intervals))
    correl = max(val_batons)/(len(E))
    return correl;
    
def LCS_correlation(X ,Y, epsilon=0.1):
    # find the length of the strings
    X = X.values; Y = Y.values
    m = len(X)
    n = len(Y)
 
    # declaring the array for storing the dp values
    L = [[None]*(n+1) for i in range(m+1)]
 
    """Following steps build L[m+1][n+1] in bottom up fashion
    Note: L[i][j] contains length of LCS of X[0..i-1]
    and Y[0..j-1]"""
    for i in range(m+1):
        for j in range(n+1):
            if i == 0 or j == 0 :
                L[i][j] = 0
            elif abs(X[i-1] - Y[j-1]) < epsilon:
                L[i][j] = L[i-1][j-1]+1
            else:
                L[i][j] = max(L[i-1][j] , L[i][j-1])
 
    # L[m][n] contains the length of LCS of X[0..n-1] & Y[0..m-1]
    return L[m][n]/min(m,n)

def lpNorm(ai, bj, exposant=2):
#    return pow(abs(ai-bj), exposant) if ai>=bj else  pow(abs(bj-ai), exposant);
    return pow(abs(ai-bj), exposant)
    
def twed_correlation(ts1, ts2, lam, mu):
    """
    lam = 0.1 => une penalite
    mu = 0 => no stiffness DtW search
    mu = inf => stronger stiffness ED
    """
    where_are_NaNs = np.isnan(ts1); ts1[where_are_NaNs] = 0;
    where_are_NaNs = np.isnan(ts2); ts2[where_are_NaNs] = 0;
    
    D = np.zeros((len(ts1), len(ts2)));
    D[1,0] = ts1[0] * ts1[0];
    D[0,1] = ts2[0] * ts2[0];
    for i in range(2, len(ts1)):
        D[i,0] = D[i-1,0] + lpNorm(ts1[i-1], ts1[i]);
    for i in range(2, len(ts2)):
        D[0,i] = D[0,i-1] + lpNorm(ts1[i-1], ts2[i]);
        
    for i in range(1, len(ts1)):
        for j in range(1, len(ts2)):
            #insertion
            if i>0 and j>1:
                insertion = D[i-1,j-1] + 2*lam * abs(i-j) + lpNorm(ts1[i-1], ts2[j-1])+ \
                            lpNorm(ts1[i], ts2[j])
            else:
                insertion = D[i-1,j-1] + lam * abs(i-j) + lpNorm(ts1[i-1], ts2[j-1]);
                    
            #delete
            if i>1 :
                delete = D[i-1,j] + lpNorm(ts1[i], ts1[i-1]) + lam + mu;
            else:
                delete = D[i-1,j] + ts1[i-1] * ts1[i-1] + lam;
                
            #match
            if j>1 :
                match = D[i,j-1] + lpNorm(ts2[j], ts2[j-1]) + lam + mu;
            else:
                match = D[i,j-1] + ts2[j-1] * ts2[j-1] + lam;
                
            D[i,j] = min(insertion, delete, match)
            
    return D[len(ts1)-1, len(ts2)-1]
    
def LB_Keogh_correlation(s1,s2,r):
    LB_sum=0
    for ind,i in enumerate(s1):

        lower_bound=min(s2[(ind-r if ind-r>=0 else 0):(ind+r)])
        upper_bound=max(s2[(ind-r if ind-r>=0 else 0):(ind+r)])

        if i>upper_bound:
            LB_sum=LB_sum+(i-upper_bound)**2
        elif i<lower_bound:
            LB_sum=LB_sum+(i-lower_bound)**2

    return math.sqrt(LB_sum)
    
def DTW_distance_correlation(s1, s2,w):
    DTW={}

    w = max(w, abs(len(s1)-len(s2)))

    for i in range(-1,len(s1)):
        for j in range(-1,len(s2)):
            DTW[(i, j)] = float('inf')
    DTW[(-1, -1)] = 0

    for i in range(len(s1)):
        for j in range(max(0, i-w), min(len(s2), i+w)):
            dist= (s1[i]-s2[j])**2
            DTW[(i, j)] = dist + min(DTW[(i-1, j)],DTW[(i, j-1)], DTW[(i-1, j-1)])
            
#    TODO: utiliser une metrique pour la mettre entre [0 et 1]

    return math.sqrt(DTW[len(s1)-1, len(s2)-1])
    
def metrique_pearson_damien(ts1, ts2, fenetre):
    """
    correlation glissante sur la fenetre
    """
    # rows with NaN only
    ts1.dropna(axis=0, how='all', inplace=True);
    ts2.dropna(axis=0, how='all', inplace=True);
    rolling_pearson_corr = ts1.rolling(window=fenetre, center=True).corr(other=ts2)
    if np.isnan(rolling_pearson_corr.mean()) or np.isinf(rolling_pearson_corr.mean()):
        print("{}".format(rolling_pearson_corr.mean() ))
        return 0;
    return abs(rolling_pearson_corr.mean());