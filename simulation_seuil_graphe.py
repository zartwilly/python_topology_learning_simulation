#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 19:17:55 2018

@author: willy
simulation seuils
"""

import pandas as pd
import numpy as np
import re, os, time
import fonctions_auxiliaires as fct_aux
import generations_mesures as mesures
import decouverte_cliques as decouvClique
#import tmp_decouverteClique_1006 as decouvClique
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import pylab as plab
import seaborn as sns
import genererMatA as geneMatA
from random import randint as valAleotoire
import random;
from pathlib import Path    # http://stackoverflow.com/questions/6004073/how-can-i-create-directories-recursively
import itertools as it
from multiprocessing import Pool;
from pandas.parser import CParserError;
import multiprocessing as mp
from scipy.stats import truncnorm;
from scipy import stats;
from bisect import bisect


import simulation_50_graphes_PARALLELE as simu50;


#------------------------------------------------------------------------------
#       generer n nombres entre 0 et 1 selon la distribution asymmetrique vers la droite
#       diviser intervalle entre 10 segments inters = [0,1,11], seg_x =[0,0.1[
#       trouver n_seg_x = le nombre de proba entre chaque segment seg_x  tel que sum_{x=[0,10]}(n_seg_x)=1 
#       generer 100 nombres aleatoires x_i uniformement entre 0 et 100 
#       si x_i = [0,n_seg_x] alors attribuer a la case 0 la proba (par exple 0.1)
#   ----- debut -----------------------
#------------------------------------------------------------------------------
def get_truncated_normal(mean=0, sd=1, low=0, upp=10):
    return truncnorm(
        (low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)

def skew_norm_pdf(x,loc=0,scale=1,a=0):
    # adapated from:
    # http://stackoverflow.com/questions/5884768/skew-normal-distribution-in-scipy
    t = (x-loc) / scale
    return 2.0 * scale * stats.norm.pdf(t) * stats.norm.cdf(a*t)
def weighted_choice(choices):
    values, weights = zip(*choices)
    total = 0; cum_weights = [];
    for w in weights:
        total += w; cum_weights.append(total)
    x = random.random() * total; i = bisect(cum_weights, x)
    return values[i]
def generate_proba_skew_normal(a_inf, a_supp, SKEW_ALPHA, location=0, scale=2):
    """
    a_inf, b_sup : borne inf/sup des valeurs de correls
    location :  la decalage par rapport a la distribution de base
    scale : facteur multiplicatif des valeurs
    SKEW_ALPHA : la pente de la courbe. -5: case a 0; 5: case a 1(la courbe de pdf n'est pas la bonne pour 5)
    """
    x = np.linspace(a_inf, a_supp, 100);
    p = []; dico_density_intervalsCorrels=dict()
    
    # generate the skew normal PDF for reference:
    p = skew_norm_pdf(x,location,scale,SKEW_ALPHA)
    # n.b. note that alpha (SKEW_ALPHA) is a parameter that controls skew, but the 'skewness'
    # as measured will be different. see the wikipedia page:
    # https://en.wikipedia.org/wiki/Skew_normal_distribution
    counts, bins = None, None;
    if SKEW_ALPHA < 0:
        counts, bins = np.histogram(p, bins=10, range=(0, 1));
    else :
        counts, bins = np.histogram(p, bins=10, range= (0, 1));
        bins = bins[::-1];
    intervals = [bins[i:i+2] for i in range(len(bins))]; intervals.pop();
    cpt_densite_inf = 0; cpt_densite_sup = 0; 
    for i in range(0,len(bins)-1):
        k = (bins[i:i+2][0], bins[i:i+2][1]);#dico_density_intervalsCorrels[k] = counts[i];
        cpt_densite_inf = cpt_densite_sup;  cpt_densite_sup += counts[i];
        dico_density_intervalsCorrels[(cpt_densite_inf, cpt_densite_sup)] = (k,counts[i]/counts.sum());
    return dico_density_intervalsCorrels;

def generate_correlation_value(cases_0_1, dico_value_correls, dico_density_intervalsCorrels,step_choiced_tup=0.001,):
    """ cases_0_1 = [(row1,col1), (row2,col2), ...]"""
    for case_0_1 in cases_0_1:
        choiced_tup = weighted_choice(dico_density_intervalsCorrels.values());
        val = random.choice(np.arange(choiced_tup[0],choiced_tup[1], step_choiced_tup));
        dico_value_correls[case_0_1] = val;
    return dico_value_correls;
    
def selectionner_cases_0_1(M_C):
    cases_0, cases_1 = list(),list();
    cpt_cases = 0;
    for row,col in fct_aux.range_2d(M_C.columns):
        cpt_cases += 1;
        if row != col and M_C.loc[row,col] == 0:
            cases_0.append((row,col));
        elif row != col and M_C.loc[row,col] == 1:
            cases_1.append((row,col));
    print("cpt_cases ={}, cases_0_1={}, cases_0={}, cases_1={}".format(\
          cpt_cases, len(cases_0)+len(cases_1), round(len(cases_0)/(len(cases_0)+len(cases_1)),3),\
          round(len(cases_1)/(len(cases_0)+len(cases_1)),3)));
    return  cases_0, cases_1;
def generer_correlations(matE,chemin_mat,SKEW_ALPHA_0_1=(-5,5), step_choiced_tup=0.001):
    """
    choisir les cases a 0 et les cases a 1 puis
    pour chaque cases_0/cases_1, generer les values de correlations
    SKEW_ALPHA_0_1[0] = -5 #cases_0 ; SKEW_ALPHA_0_1[1] = 5;#cases_1
    """
    cases_0, cases_1 = selectionner_cases_0_1(matE);
    dico_value_correls = dict();
    dico_density_intervalsCorrels = generate_proba_skew_normal(\
                                    a_inf=0, a_supp=1, SKEW_ALPHA=SKEW_ALPHA_0_1[0], \
                                    location=0, scale=2);
    dico_value_correls = generate_correlation_value(cases_0, dico_value_correls, \
                         dico_density_intervalsCorrels,step_choiced_tup=0.001);
    dico_density_intervalsCorrels = generate_proba_skew_normal(\
                                    a_inf=0, a_supp=1, SKEW_ALPHA=SKEW_ALPHA_0_1[1], \
                                    location=1.50, scale=2);
    dico_value_correls = generate_correlation_value(cases_1, dico_value_correls, \
                         dico_density_intervalsCorrels,step_choiced_tup=-0.001);
    # plot distribution
    f, ax = plt.subplots(figsize=(12,4));
    hist, bins = np.histogram(list(dico_value_correls.values()));
    ax.bar(bins[:-1], hist.astype(np.float32) / hist.sum(), width=(bins[1]-bins[0]), color='blue')
    cases01 = len(cases_0)+len(cases_1); n0=round(len(cases_0)/cases01,3);n1=round(len(cases_1)/cases01,3);
    label = "distribution cases 0 et 1: cases0={}, cases1={}, cases01={}".format(n0,n1, cases01)
    ax.set_title(label);
    f.savefig(chemin_mat+"distribution01.jpeg");
    plt.clf();
                           
    M_C = matE.copy();
    for arete, correl in dico_value_correls.items():
         M_C.loc[arete[0],arete[1]] = correl;
         M_C.loc[arete[1],arete[0]] = correl;                                 
    return M_C, dico_value_correls; 
#------------------------------------------------------------------------------
#   ----- fin -----------------------
#------------------------------------------------------------------------------
def matrice_binaire_seuil(M_C, matE_LG, seuil):
    cols = M_C.columns.tolist(); #cols_LG = matE_LG.columns.tolist();
#    print("cols_M_C={}, cols_LG={}".format(cols, cols_LG))
    matE = M_C.copy(); faux_pos = []; faux_neg = [];
    for row,col in fct_aux.range_2d(cols):
        if M_C.loc[row,col] > seuil:
            matE.loc[row,col] = 1;
            matE.loc[col,row] = 1;
        else:
            matE.loc[row,col] = 0;
            matE.loc[col,row] = 0;
        if matE_LG.loc[row,col] == 0 and matE.loc[row,col] == 1:
            faux_pos.append((row,col));
        elif matE_LG.loc[row,col] == 1 and matE.loc[row,col] == 0:
            faux_neg.append((row,col));
    return matE, faux_pos, faux_neg;
    
def calculer_faux_pos_neg_correction(aretes_LG, matE):
    faux_pos_correct, faux_neg_correct = list(), list();
    for row,col in fct_aux.range_2d(matE.columns):
        if (row,col) in aretes_LG or (col,row) in aretes_LG:
            if matE.loc[row, col] == 0:
                faux_pos_correct.append((row,col));
        elif matE.loc[row, col] == 1:
            faux_neg_correct.append((row,col));
    return faux_pos_correct, faux_neg_correct;
    
def save_df(df, path_distr_chemin, identifiant, headers):
    """
    merge df to the .csv loaded dataframe 
    """
    if not os.path.exists(path_distr_chemin+"resumeExecution_"+str(identifiant)+".csv"):
        df.colums = headers
        df.to_csv(path_distr_chemin+"resumeExecution_"+str(identifiant)+".csv", sep=",")
    else:
        df_csv = pd.read_csv(path_distr_chemin+"resumeExecution_"+str(identifiant)+".csv", \
                         sep=",", names = headers)
        pd.concat([df_csv, df])\
                  .to_csv(path_distr_chemin+"resumeExecution_"+str(identifiant)+".csv", sep=",")
    

def simulation_seuil(args, seuil):
    """
    """    
#    print("****seuil = {}, chemin_distrib={}".format(seuil, args["chemin_distrib"]));
    # fichier de tracking or debug
    headers_df = ["G_cpt","seuil","num_graphe","nbre_aretes_matE_s","nbre_aretes_LG",\
                  "nbre_aretes_diff_matE_s_LG","dist_line","liste_aretes_diff_matE_s_LG","nbre_aretes_diff_matE_LG","hamming",\
                  "liste_aretes_diff_matE_LG","C","som_cout_min","noeuds_corriges", "min_hamming",\
                  "mean_hamming","max_hamming","ecart_type","max_cout","max_permutation",\
                  "dico_som_min_permutations", "dico_arc_sommet", "ordre_noeuds_traites","C_old",\
                  "faux_pos_seuil","faux_pos_apres_correction", "faux_neg_seuil",\
                  "faux_neg_apres_correction" ]
            
    df_debug_s = pd.DataFrame( columns = headers_df);
    G_s = "G_"+str(seuil)+"_"+str(args["num_graphe"]);
                                      
    # initialisation variables    
    moy_distline = 0; moy_hamming = 0; sum_distline = 0; sum_hamming = 0; correl_dl_dh = 0;
    faux_pos_correct, faux_neg_correct, faux_pos_seuil, faux_neg_seuil = list(),list(),list(),list();
    # application du seuil
    
#    print("(**) M_C={}, matE={}".format(args["M_C"].columns.tolist(), args["matE"].columns.tolist()))
    matE_s, faux_pos_seuil, faux_neg_seuil = matrice_binaire_seuil(args["M_C"].copy(), \
                                                       args["matE"].copy(), seuil);
#    print("MC = \n{}, matE_s={}\n".format(args["M_C"], matE_s));                                 
    aretes_matE_s = len(fct_aux.liste_arcs(matE_s));
    try:
        dico_permutation_cliq = dict();
        #algo corrigeant tous les noeuds a -1
#        print("1");  
        args["correl_seuil"] = seuil;
        dico_permutation_cliq = \
            decouvClique.decouverte_cliques(matE_s, args["dico_arc_sommet"], \
                                    args["seuil_U"], args["epsilon"], \
                                    args["chemin_datasets"], args["chemin_matrices"],\
                                    args["ascendant_1"], args["simulation"],\
                                    args["dico_proba_cases"],\
                                    args);
#        print("2");
        # Debut selection de la permutation de noeuds dont la distance hamming est la plus petite
        dico_sol = dict()
        dico_sol = simu50.best_permutation(dico_permutation_cliq, args["matE"], matE_s);
        # FIN selection de la permutation de noeuds dont la distance hamming est la plus petite
        dico_som_min_permutations = dict();
        for l_noeuds_1, values in dico_permutation_cliq.items():
            if values[6] not in dico_som_min_permutations.keys():
                dico_som_min_permutations[values[6]] = [l_noeuds_1]
            else:
                dico_som_min_permutations[values[6]].append(l_noeuds_1);
        faux_pos_correct, faux_neg_correct = calculer_faux_pos_neg_correction(\
                                            dico_sol["aretes_LG"],args["matE"]);
        #--- ajouter le debug par seuil ---> DEBUT
        cpt_s = args["occurence_seuil"];
        df_debug_s.loc[int(cpt_s)] = [G_s, seuil, args["num_graphe"],\
                    dico_sol["nbre_aretes_matE"], \
                    dico_sol["nbre_aretes_LG"], dico_sol["liste_aretes_diff_matE_k_alpha_LG"],\
                    dico_sol["dist_line"], dico_sol["liste_aretes_diff_matE_k_alpha_LG"], \
                    dico_sol["nbre_aretes_diff_matE_LG"], dico_sol["hamming"],\
                    dico_sol["liste_aretes_diff_matE_LG"],\
                    dico_sol["C"], dico_sol["som_cout_min"], \
                    dico_sol["noeuds_corriges"], \
                    dico_sol["min_hamming"],dico_sol["mean_hamming"], \
                    dico_sol["max_hamming"],dico_sol["ecart_type"],\
                    dico_sol["max_cout"], dico_sol["max_permutation"],\
                    dico_som_min_permutations, args["dico_arc_sommet"],\
                    dico_sol["ordre_noeuds_traites"], dico_sol["C_old"],\
                    faux_pos_seuil,faux_pos_correct, faux_neg_seuil,faux_neg_correct]
        
        #--- ajouter le debug par seuil ---> FIN
        sum_distline += dico_sol["dist_line"] # dico_sol["dist_line"]/dico_sol["nbre_aretes_matE_k_alpha"] 
        sum_hamming += dico_sol["hamming"] # dico_sol["hamming"]/dico_sol["nbre_aretes_matE"]
        print("{},s={},num_graphe={}, noeuds_1={},fct_cout={} --> TERMINE".\
              format(G_s, seuil, args["num_graphe"], args["mode_select_noeuds_1"], args["coef_fct_cout"][2])) 
              
    except Exception as e:
        print("####### EmptyDataError fct_cout={} noeuds_1={}, s={}, num_graphe={}, {} : e = {} #######".format(
              args["coef_fct_cout"][2],arg_params["mode_select_noeuds_1"], \
              seuil, args["num_graphe"], G_s, e));
        # a ajouter le debug seuil
        cpt_s = args["occurence_seuil"]
        df_debug_s[cpt_s] = [G_s, seuil, args["num_graphe"],\
                    "error", "error","error","error", "error",\
                    "error", "error","error","error", "error", \
                    "error", "error","error","error", "error",\
                    "error", "error","error", args["dico_arc_sommet"],"error", \
                    "error"]
    moy_distline = sum_distline;
    moy_hamming = sum_hamming;
    if moy_hamming == 0 and moy_distline == 0:
        correl_dl_dh = 1
    else:
        correl_dl_dh = abs(moy_hamming - moy_distline)/max(moy_hamming, moy_distline)
            
    # ecrire dans un fichier pouvant etre lu pendant qu'il continue d'etre ecrit
    f = open(args["chemin_distrib"]+"distribution_moyDistLine_moyHamming_s_"+str(seuil)+".txt","a")
    f.write(G_s+";"+str(seuil)+";"+str(moy_distline)+";"+str(moy_hamming)+";"+str(aretes_matE_s)+\
            ";"+str(correl_dl_dh)+";"+str(len(faux_pos_seuil))+";"+str(len(faux_neg_seuil))+";"+\
            str(len(faux_pos_correct))+";"+str(len(faux_neg_correct))+"\n")
    f.close();
    save_df(df_debug_s, args["chemin_distrib"], seuil, headers_df)
    
 
#-------------------------------------
#     generer N graphes ---> debut

def generer_N_graphes(dimMatA, nbre_graphes,args):
    graphes = list()
#    path_Gs = Path(args["rep"]); path_Gs.mkdir(parents=True, exist_ok=True);
    for nbre_graphe in range(nbre_graphes):
        for type_fct_cout in args["type_fct_couts"]:
            for mode_select_noeuds_1 in args["mode_select_noeuds_1s"]:
                matE = None;
                chemin_mat = args["rep"]+"/"+type_fct_cout+"/"+mode_select_noeuds_1+"/"+"G_"+str(nbre_graphe)+"/matrices/"
                chemin_data = args["rep"]+"/"+type_fct_cout+"/"+mode_select_noeuds_1+"/"+"G_"+str(nbre_graphe)+"/datasets/"
                chemin_distrib = args["rep"]+"/"+type_fct_cout+"/"+mode_select_noeuds_1+"/distribution/";
                path_Gs = Path(chemin_mat); path_Gs.mkdir(parents=True, exist_ok=True);
                path_Gs = Path(chemin_data); path_Gs.mkdir(parents=True, exist_ok=True);
                path_Gs = Path(chemin_distrib); path_Gs.mkdir(parents=True, exist_ok=True);
                matE, matA, dico_arc_sommet = simu50.matriceE(dimMatA, args["nbre_lien"], chemin_mat,\
                                                           chemin_data, args["nbre_ts"], args["epsilon"], \
                                                           args["effet_joule"], "Pas test");
                """
                attribuer des valeurs de correlations pour chaque case de matE
                et mettre dans M_C
                """
                M_C = matE.copy();
                M_C, dico_proba_cases = generer_correlations(M_C,chemin_mat,args["SKEW_ALPHA_0_1"])
                graphes.append( (M_C, matE, matA, dico_arc_sommet,nbre_graphe, \
                                 chemin_mat, chemin_data, chemin_distrib, dico_proba_cases) );
#                print("*** M_C={}, matE={}".format(M_C.columns.tolist(), matE.columns.tolist()))
    return graphes;
#     generer N graphes ---> fin
#-------------------------------------

#-------------------------------------
##     representation seuil ---> debut
def lire_seuils(path_files):
    """ ----distribution_moyDistLine_moyHamming_s_0.1 ----"""
    files = os.listdir(path_files)
    liste_files = [fichier for fichier in files if 
                   re.match("^distribution_moyDistLine_moyHamming_s_.*",fichier)]
 
    seuils = [file.split("_")[-1].split(".txt")[0] for file in liste_files];
    return seuils;
def create_dataframe(args, seuils):
    """create un dataframe contenant les moy_DH ou moy_DL de chaque seuil donnee """
    df = pd.DataFrame(); rep = args["path_files"];
    for seuil in seuils:
        df_seuil = pd.read_csv(rep+args["root_files"]+seuil+args["ext"], \
                         names=["cpt","moy_dl"+"_s_"+str(seuil),\
                                "moy_dh"+"_s_"+str(seuil), \
                                "nbre_aretes_matE"+"_s_"+str(seuil), \
                                "correl_dh_dl"+"_s_"+str(seuil)], \
                         sep=';');
#        df_seuil["moy_dh"+"_s_"+str(seuil)] = round( df_seuil["moy_dh"+"_s_"+str(seuil)] / df_seuil["nbre_aretes_matE"+"_s_"+str(seuil)], 3)
#        df_seuil["moy_dl"+"_s_"+str(seuil)] = round( df_seuil["moy_dl"+"_s_"+str(seuil)] / df_seuil["nbre_aretes_matE"+"_s_"+str(seuil)], 3)
#        df_seuil = df_seuil.sort_values(by=["moy_dl"+"_s_"+str(seuil),"moy_dh"+"_s_"+str(seuil)]).\
#                    reset_index(drop=True);
        df_seuil = df_seuil.sort_values(by=["moy_dh"+"_s_"+str(seuil),"moy_dl"+"_s_"+str(seuil)]).\
                    reset_index(drop=True);            
        #df_seuil.reset_index(drop=True,inplace=True);                               
        df = pd.concat( [df, df_seuil], axis = 1);
    return df;
    
def plot_seuil(df, path_save, selected_var):
    if not os.path.exists(path_save):
        path_distr = Path(path_save);
        path_distr.mkdir(parents=True, exist_ok=True);
    
    cols = df.columns.tolist();print("cols={}".format(cols))
#    df = df.sort_values(by = cols);
    fig = plt.figure(); default_size = fig.get_size_inches()
    print("w =", default_size[0], " h = ",default_size[1])
    fig.set_size_inches( (default_size[0]*1.5, default_size[1]*1.5) )
    ax1 = fig.add_subplot(1,1,1);
        
    styles1 = ['bs-','ro-','y^-','rs-','go-','b^-','r*-','bo-','-gD','-yp',\
                  ':>','-.<','-v','-d','-h','--H','--,']
    df.plot(style=styles1, ax = ax1)
        
    ax1.set(xlabel= "graphes", ylabel= selected_var.upper(), \
            title = "choix du seuil en fonction des "+str(selected_var.upper() ));
    font_prop = font_manager.FontProperties( size=14)
    ax1.legend( df.columns.tolist(), loc='upper center', bbox_to_anchor=(0.5, 1.00),\
                ncol=2, fancybox=True, shadow=True, prop = font_prop);
                   
    ### modif font size
    for ax in [ax1]:
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
         ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(12)
            
    nom_fichier = "choixSeuils_by_"+ selected_var.upper();         
    plt.savefig(path_save+"/"+nom_fichier+".jpeg",dpi= 190, bbox_inches='tight')
    plt.clf();

def plot_choix_seuil_selon_variables(args):
    seuils = lire_seuils(args["path_files"]); print("seuils={}".format(seuils));
    df = create_dataframe(args, seuils);
    for var in args["selected_vars"]:
        selected_vars = [ var+"_s_"+str(s) for s in seuils];
        plot_seuil(df[selected_vars], args["path_save"], var);
    return df
##     representation seuil ---> fin
#-------------------------------------


if __name__ == '__main__':
    start= time.time();
    dimMatA = 20; nbre_graphes = 10#150#10#50;
    nbre_lien = (2,5); nbre_ts = 10; epsilon = 0.75; effet_joule = 0.1;
    test = "FINI";ascendant_1 = True; simulation = True; seuil_U = 0;
    
    critere_selection_pi1_pi2 = 2; # 0: moins de modif,1: ajout aretes> supp aretes, 2:ajout aretes < supp aretes,
    number_permutations_nodes_1= 10 #100;
    facteur_multiplicatif = 1; exposant = 1; # 0: fct_unitaire, 1:fct_normal, 2: fct_quadratique, 4:fct_quadruple, 5: fct_quintuple
    type_fct_cout = "cloche" # ou "lineaire" ou "cloche"
    coef_fct_cout = (exposant, facteur_multiplicatif, type_fct_cout)
    mode_select_noeuds_1 = "aleatoire" #"degreMin" #"coutMin" # "degreMin" # aleatoire
    number_items_pi1_pi2 = 0.5#1;

    arg_params = {"nbre_lien":nbre_lien, "nbre_ts":nbre_ts, "epsilon":epsilon, \
                  "effet_joule":effet_joule, "test": test, "ascendant_1":ascendant_1,\
                  "simulation": simulation, "seuil_U":seuil_U,\
                  "number_items_pi1_pi2": number_items_pi1_pi2,\
                   "number_permutations_nodes_1": number_permutations_nodes_1, \
                   "mode_select_noeuds_1":mode_select_noeuds_1,\
                   "coef_fct_cout":coef_fct_cout,\
                   "critere_selection_pi1_pi2":critere_selection_pi1_pi2};
    
    rep = "simulationSeuils_fct_cout_unitaire";
    type_fct_couts = ["lineaire","cloche"];
    type_fct_couts = ["lineaire_simul50Graphes_priorite_supp", "lineaire_simul50Graphes_priorite_ajout", "lineaire_simul50Graphes_priorite_aucune"];
    mode_select_noeuds_1s = ["degreMin","coutMin","aleatoire"]; 
    SKEW_ALPHA_0_1 = (-5,3.7)#(-5,5); #SKEW_ALPHA_0_1 = (-7,5);
    correl_seuils = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]; N_seuils = len(correl_seuils);
    
    mode_select_noeuds_1s = ["aleatoire"];
#    correl_seuils = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]; N_seuils = len(correl_seuils);
    correl_seuils = [0.5]; N_seuils = len(correl_seuils);
    mode_select_noeuds_1s = ["aleatoire"]; type_fct_couts = ["lineaire"];
    
    #### cas fonction de cout unitaire OUBLIE
    correl_seuils = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]; N_seuils = len(correl_seuils);
    mode_select_noeuds_1s = ["aleatoire"]; type_fct_couts = ["lineaire"];
    facteur_multiplicatif = 1; exposant = 0; type_fct_cout = "lineaire";
    coef_fct_cout = (exposant, facteur_multiplicatif, type_fct_cout)
    arg_params["coef_fct_cout"] = coef_fct_cout;
    dimMatA = 8; nbre_graphes = 150#10#50;
    #### cas fonction de cout unitaire OUBLIE
    
    args_graphes = {"rep":rep,"nbre_lien":nbre_lien,"nbre_ts":nbre_ts,\
                    "epsilon":epsilon, "effet_joule":effet_joule, "SKEW_ALPHA_0_1":SKEW_ALPHA_0_1,\
                    "type_fct_couts":type_fct_couts,"mode_select_noeuds_1s":mode_select_noeuds_1s}
    graphes = [];
    graphes = generer_N_graphes(dimMatA, nbre_graphes, args_graphes);
                   
    
    params = list(); 
    for type_fct_cout in type_fct_couts:
        for mode_select_noeuds_1 in mode_select_noeuds_1s:
            occurence_seuil = 0;
            for graphe in graphes:
                arg_params["mode_select_noeuds_1"] = mode_select_noeuds_1;
                arg_params["coef_fct_cout"] = (arg_params["coef_fct_cout"][0], \
                                        arg_params["coef_fct_cout"][1],\
                                        type_fct_cout);                
                arg_params["M_C"] = graphe[0]; arg_params["matE"] = graphe[1];
                arg_params["dico_arc_sommet"] = graphe[3];
                arg_params["num_graphe"] = graphe[4];
                arg_params["chemin_matrices"] = graphe[5];
                arg_params["chemin_datasets"] = graphe[6];
                arg_params["chemin_distrib"] = graphe[7];
                arg_params["dico_proba_cases"] = graphe[8];
                occurence_seuil += 1
                arg_params["occurence_seuil"] = occurence_seuil;
                
                arg_params_ = arg_params.copy();
            
                params_ = list(zip( [arg_params_]*N_seuils, correl_seuils))
                params.append(params_)
    params = [par for tup in params for par in tup]
    print("params ",len(params))
    
    p = Pool(mp.cpu_count()-1) 
    p.starmap(simulation_seuil, params); p.terminate()
    s_min_max = "01_09"; g = open("tempsExecution_seuil"+str(s_min_max)+".txt","w");
#    #### A EFFACER
##    for param in params:
##        simulation_seuil(param[0], param[1])
##    s_min_max = "01_09"; g = open("tempsExecution_seuil"+str(s_min_max)+".txt","w");
#    #### A EFFACER
#    ti = time.time() - start; g.write(str(ti)); g.close()
    
    ##### plot choix seuils ---> DEBUT ###
#    df = None;
#    for type_fct_cout in type_fct_couts:
#        for mode_select_noeuds_1 in mode_select_noeuds_1s:
#            path_files = rep+"/"+type_fct_cout+"/"+mode_select_noeuds_1+"/distribution/";
#            root_files = "distribution_moyDistLine_moyHamming_s_"; ext = ".txt";
#            path_save = rep+"/"+type_fct_cout+"/"+mode_select_noeuds_1+"/courbes/";
#            selected_vars = ["moy_dl", "moy_dh"];
#            args_={"path_files":path_files, "root_files":root_files,"ext":ext,\
#                   "selected_vars":selected_vars,"path_save":path_save};
#            df = plot_choix_seuil_selon_variables(args_)
#    ##### plot choix seuils ---> FIN #####
    
    print("runtime = {}".format(time.time() - start))
    
    
    