#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 10:49:24 2018

@author: willy
"""
import re 
import os 
import time
import math
import random
import logging
import pandas as pd
import numpy as np
import fonctions_auxiliaires as fct_aux
import generations_mesures as mesures
import decouverte_cliques as decouvClique
import genererMatA as geneMatA
import algorithme_correction as algoCorrection

#import matplotlib.pyplot as plt
from pathlib import Path    # http://stackoverflow.com/questions/6004073/how-can-i-create-directories-recursively
import itertools as it
from multiprocessing import Pool;
from pandas.parser import CParserError;
import multiprocessing as mp
from scipy.stats import truncnorm
from random import randint as valAleotoire



def distance_hamming(aretes_matE_k_alpha,  aretes_matE_corrige):
    """ calculer le nombre de cases modifies entre G_k_alpha et G 
    
    aretes_matE_corrige = {frozenset(a,b), ... }
    """
    aretes_matE_k_alpha = set(map(frozenset, aretes_matE_k_alpha))
    aretes_matE_corrige = set(map(frozenset, aretes_matE_corrige))
    aretes_diff = aretes_matE_k_alpha.union(aretes_matE_corrige) - \
                  aretes_matE_k_alpha.intersection(aretes_matE_corrige);
                  
    return len(aretes_diff), aretes_diff;

def comparer_cliques(C, C_old):
    """ retourne les cliques differentes et identiques entre C et C_old
    """
    C_old = set(map(frozenset, C_old))
    cliques_identiques = set();
    cliques_differentes = set();
            
    cliques_identiques = C.intersection(C_old);
    cliques_differentes = C.union(C_old)- C.intersection(C_old);
    
    return cliques_identiques, cliques_differentes;
    
def is_locked(filepath, df, G_k):
    """Checks if a file is locked by opening it in append mode.
    If no exception thrown, then the file is not locked.
    """
    locked = None
    file_object = None
    
    try:
        print("Trying to open {} resumeExecution.".format(G_k))
        buffer_size = 8
        # Opening file in append mode and read the first 8 characters.
        file_object = open(filepath, 'a', buffer_size)
        if file_object:
            df_resExec = pd.read_csv(filepath, sep = ",", index_col=0).\
                            reset_index();
            # merger df_resExec et df en gardant les index (fusionner leur index)
            df_resExec = pd.merge(df_resExec, df, on="index", how="outer");
            df_resExec.to_csv(filepath, sep=',', index=False);
            locked = False;
    except IOError as message:
        print("resumeExecution_{}.csv is not locked ({}).".format( \
                  G_k.split("_")[2], message ))
        locked = True;
    finally:
        if file_object:
            file_object.close();
            print("resumeExecution_{}.csv  closed.".format( \
                  G_k.split("_")[2]))
    
    return locked;
    
def sauver_df_resume(df, name_save_df, G_k):
    """ sauvegarder le dataframe contenant la colonne G_numeroGraphe_k dans le dataframe generale 
    resumeExecution.csv
    
    df : contient une seule colonne "G_numeroGrapke_k"
    """
    # open file
    # verifier si ce fichier nest pas ouvert par un autre fichier
    #   si oui attendre
    #   sinon enregistrer.
    my_file = Path(name_save_df);
    
    temps_attente = 0.010;                                                      # attente de 10 ms
    if my_file.is_file():
        while is_locked(name_save_df, df, G_k):
            print("reseumeExecution_{} is currently in use. Waiting {} milliseconds.".\
                  format((G_k.split("_")[2], temps_attente)))
            time.sleep(temps_attente);
    else:
        df.to_csv(name_save_df, sep=',', index=False);
###############################################################################
#               generation graphes de flots ---> debut
###############################################################################  
def nommage_arcs(mat):
    """ nommage des arcs du graphe du reseau energetique.
    
    elle retourne un dictionnaire avec :
        la cle : le nom de l'arc
        la valeur : le tuple correspondant a l'arc
        
    """
    dico_arcs_sommets = dict();
    for arc in fct_aux.liste_arcs(mat):
        dico_arcs_sommets["_".join(arc)] = arc;
    return dico_arcs_sommets
    
def creer_reseau(chemin_datasets, chemin_matrices, args):
    """creation d'un reseau energetique.
    
    creer un graphe, ajouter des mesures de puissance sur le graphe puis
    extraire son line-graphe
    
    """
    if args["dbg"]:
        chemin_datasets = "dataNewCriterecorrectionGrapheConnu/datasets/";
        chemin_matrices = "dataNewCriterecorrectionGrapheConnu/matrices/";
        path_ = Path(chemin_datasets);
        path_.mkdir(parents=True, exist_ok=True);
        path_ = Path(chemin_matrices);
        path_.mkdir(parents=True, exist_ok=True);
        
    logger = logging.getLogger('creer_reseau')
    logger.debug("creation de mat et matE")
    dico_graphe = {"a":["b"], "b":["c","d"], "c":["e","f"], 
                   "d":["f"], "e":["g"], "f":["h"], "g":[],
                   "h":[]};
                   
    mat = pd.DataFrame(index = dico_graphe.keys(), columns = dico_graphe.keys());
    
    for k, vals in dico_graphe.items():
        for v in vals:
            mat.loc[k,v] = 1;
    mat.fillna(value=0, inplace=True);
    
    # ajouter mesures 
    grandeurs = ["P"];    
    dico_arcs_sommets = dict()
    dico_arcs_sommets = nommage_arcs(mat)
    
#    genererMesures_all_grandeurs(matA, dico_dual_arc_sommet, liste_grandeurs, location = "data/datasets/", taille = 3, effet_joule = 0):

    mesures.genererMesures_all_grandeurs(mat, dico_arcs_sommets, grandeurs, 
                                 chemin_datasets, taille=100, effet_joule=0.1)

    
    #matrice du linegraphe du reseau de flot a determiner
    arcs = fct_aux.liste_arcs(mat)
    matE = mesures.creation_matE(dico_arcs_sommets, arcs)
    matE.to_csv(chemin_matrices+"matE.csv")
    logger.debug("mat cree : arcs={}, sommets={}, matE cree : aretes={}, ".
                 format(len(arcs), len(dico_graphe.keys()), 
                            fct_aux.liste_arcs(matE)))
    return matE, mat, dico_arcs_sommets;
    
def generer_reseau(dim_mat, nbre_lien,
                   chemin_datasets,
                   chemin_matrices,
                   nbre_ts, epsilon, effet_joule):
    """ generer un reseau de sommets = dim_mat avec ses mesures de flots.
    
        dim_mat : l'ordre du graphe
        chemin_datasets : chemin pour sauvegarder les mesures de 
                        grandeurs physiques
        chemin_matrices : chemin pour sauvegarder les matrices de notre 
                        reseau cad matrice d'adjacence du linegraphe (matE) 
                        et son graphe racine (mat)
        epsilon : 0.75
        seuil : valeur par defaut definissant une adjacence entre 2 sommets
        nbre_ts : nombre de time series 
    """
    
    #generer reseau de flots (A DETERMINER) avec mesures 
    matA = None;
    matA = geneMatA.genererMatriceA(dim_mat, nbre_lien)
    matA.to_csv(chemin_matrices+"mat_generer.csv")
    dico_arcs_sommets = mesures.nommage_arcs( mat )
    mesures.create_datasets(mat, dico_arcs_sommets, 
                            chemin_datasets, nbre_ts, effet_joule) 
    
    #matrice du linegraphe du reseau de flot a determiner
    arcs = fct_aux.liste_arcs(mat)
    matE = mesures.creation_matE(dico_arcs_sommets, arcs)
    matE.to_csv(chemin_matrices+"matE_generer.csv")
    return matE, mat, dico_arcs_sommets;
###############################################################################
#               generation graphes de flots ---> fin
############################################################################### 

###############################################################################
#               ajout de valeur de probabilites aux cases de matE---> debut
###############################################################################
def get_truncated_normal(mean=0, sd=1, low=0, upp=10):
    return truncnorm((low - mean) / sd, 
                     (upp - mean) / sd, 
                     loc=mean, scale=sd)
    
def loi_de_probalibilite(debut_proba, fin_proba, 
                         nbre_correlations, correl_seuil):
    """
    generer des valeurs de correlations entre debut_proba et fin_proba.
    """
    mean = (fin_proba + debut_proba)/2; 
    sd = abs((fin_proba - debut_proba)/(2 * 4)) \
            if fin_proba != debut_proba else abs((fin_proba)/(2 * 4))
    correlations  = get_truncated_normal(mean * 100, 
                                         sd * 100, 
                                         low = 1, 
                                         upp = nbre_correlations * 100);
    return np.random.choice(correlations.rvs(nbre_correlations) / 100);
    
def ajouter_probabilite(dico_proba_cases, cases, loi_proba, 
                        p_correl, val_limite_proba_1, 
                        case_a_selectionner):
    """ ajouter  des valeurs de probabilites a chaque case de cases 
        selon la type de cases 'case_a_selectionner' et 
                la loi de proba 'loi_proba'.
        
        val_limite_proba_1 :  limite inferieure de toutes les cases a 1.
    """
    if loi_proba == "poisson":
        if case_a_selectionner == "0_0":
            for case in cases:
                dico_proba_cases[case] = round(
                                        random.choice(
                                        np.union1d(np.linspace(0, 0.117, 100), \
                                        np.linspace(0.35, 0.5, 100))),
                                                      3)
        elif case_a_selectionner == "1_1":
            for case in cases:
                debut_proba = val_limite_proba_1;                               # 0.8
                fin_proba = 1;
                dico_proba_cases[case] = loi_de_probalibilite(
                                                debut_proba, 
                                                fin_proba,
                                                val_limite_proba_1,
                                                nbre_correlations = 250);
        elif case_a_selectionner == "1_0":
            for case in cases:
                dico_proba_cases[case] = round(
                                        random.choice(
                                        np.union1d(np.linspace(0.6, 0.644, 100), \
                                        np.linspace(0.72, 0.79, 100))),
                                                3)
        
        elif case_a_selectionner == "0_1":
            for case in cases:
                debut_proba = p_correl;
                fin_proba = val_limite_proba_1 - 0.001;                         # 0.79
                dico_proba_cases[case] = loi_de_probalibilite(
                                                debut_proba, 
                                                fin_proba,
                                                val_limite_proba_1,
                                                nbre_correlations = 250);
    elif loi_proba == "gaussienne":
        if case_a_selectionner == "0_0":
            for case in cases:
                debut_proba = 0;
                fin_proba = val_limite_proba_1 - p_correl - 0.01;               # 0.3
                dico_proba_cases[case] = loi_de_probalibilite(
                                                debut_proba, 
                                                fin_proba,
                                                val_limite_proba_1,
                                                nbre_correlations = 250);
        elif case_a_selectionner == "1_1":
            for case in cases:
                debut_proba = val_limite_proba_1;                               # 0.8
                fin_proba = 1;
                dico_proba_cases[case] = loi_de_probalibilite(
                                                debut_proba, 
                                                fin_proba,
                                                val_limite_proba_1,
                                                nbre_correlations = 250);
        elif case_a_selectionner == "1_0":
            for case in cases:
                debut_proba = val_limite_proba_1 - p_correl;                    # 0.3
                fin_proba = p_correl;
                dico_proba_cases[case] = loi_de_probalibilite(
                                                debut_proba, 
                                                fin_proba,
                                                val_limite_proba_1,
                                                nbre_correlations = 250);
        
        elif case_a_selectionner == "0_1":
            for case in cases:
                debut_proba = p_correl;
                fin_proba = val_limite_proba_1 - 0.001;                         # 0.79
                dico_proba_cases[case] = loi_de_probalibilite(
                                                debut_proba, 
                                                fin_proba,
                                                val_limite_proba_1,
                                                nbre_correlations = 250);
    return dico_proba_cases;
###############################################################################
#               ajout de valeur de probabilites aux cases de matE---> fin
###############################################################################

###############################################################################
#               modification de k cases de matE---> debut
############################################################################### 
def modifier_k_cases(matE, p_correl, k, loi_proba, val_limite_proba_1):
    """ modification de k cases selon la proba p_correl.
    
    si p_correl = 0, toutes les cases sont a 1
    si p_correl = 1, toutes les cases sont a 0
    
    NE PAS OUBLIER D ATTRIBUER LES VALEURS DE PROBAS POUR CHAQUE CASE DE MATE
    """
    
    dico_k_erreurs = dict();
    dico_k_erreurs["aretes_ajoutees"] = []; 
    dico_k_erreurs["aretes_supprimees"] = [];
    
    cases_0 = fct_aux.liste_nonArcs(matE, 0);
    cases_1 = fct_aux.liste_arcs(matE);
    dico_proba_cases = dict();
    dico_proba_cases = ajouter_probabilite(dico_proba_cases, 
                                           cases_0, 
                                           loi_proba,
                                           p_correl,
                                           val_limite_proba_1, 
                                           case_a_selectionner = "0_0");
    dico_proba_cases = ajouter_probabilite(dico_proba_cases, 
                                           cases_1, 
                                           loi_proba,
                                           p_correl,
                                           val_limite_proba_1,
                                           case_a_selectionner = "1_1")                                       
    for k_ in range(0,k):
        proba = random.random();
    
        m0 = len(cases_0);
        m1 = pow(matE.shape[0],2) - m0;
        if proba < p_correl:
            # case 0 -> 1
            p0 = 1/m0;
            bool = True;
            while bool:
                arete = random.choice(cases_0)
                if random.random() > p0 and \
                    arete not in dico_k_erreurs["aretes_ajoutees"] and \
                    arete not in dico_k_erreurs["aretes_supprimees"] :
                    matE.loc[arete[0], arete[1]] = 1; 
                    matE.loc[arete[1], arete[0]] = 1;
                    dico_k_erreurs["aretes_ajoutees"].append(arete);
                    cases_1.append(arete);
                    dico_proba_cases = ajouter_probabilite(dico_proba_cases, 
                                                           [arete], 
                                                           loi_proba,
                                                           p_correl,
                                                           val_limite_proba_1,
                                                           "0_1");
                    bool = False;
                    break;
        else:
            # case 1 -> 0
            p1 = 1/m1;
            bool = True;
            while bool:
                arete = random.choice(cases_1)
                if random.random() > p1 and \
                    arete not in dico_k_erreurs["aretes_ajoutees"] and \
                    arete not in dico_k_erreurs["aretes_supprimees"]:
                    matE.loc[arete[0], arete[1]] = 0; 
                    matE.loc[arete[1], arete[0]] = 0;
                    dico_k_erreurs["aretes_supprimees"].append(arete);
                    cases_0.append(arete);
                    dico_proba_cases = ajouter_probabilite(dico_proba_cases, 
                                                           [arete], 
                                                           loi_proba,
                                                           p_correl,
                                                           val_limite_proba_1,
                                                           "1_0");
                    bool = False;
                    break;
    return matE, dico_k_erreurs, dico_proba_cases;    
    pass

###############################################################################
#               modification de k cases de matE---> fin
############################################################################### 

###############################################################################
#               simulation nouveau critere sur un graphe connu ---> debut
############################################################################### 
def simulation_nouveau_critere(args):
    """ simulation d'un nouveau critere de correction.
    
    """
    logging.basicConfig(format='%(asctime)s - %(message)s', 
                        level=logging.DEBUG,
                        datefmt='%d-%b-%y %H:%M:%S',
                        filename=args["log_file"],
                        filemode="w")
    logger = logging.getLogger('***** simulation_nouveau_critere')
    matE, mat, dico_arcs_sommets = creer_reseau(args["chemin_datasets"], 
                                                args["chemin_matrices"], args);             
    
    
    # supprimer ou ajouter des aretes dans matE
    logger.debug("***** Suppression et ajout aretes dans matE")
    matE_k_alpha = matE.copy()
    rd = random.random()
    dico_k_erreurs = {"aretes_ajoutees":[], "aretes_supprimees":[]}
    if args["k_erreurs"] == 1 and rd >= args["p_correl"]:
        matE_k_alpha.loc["b_d","c_f"] = 1;
        matE_k_alpha.loc["c_f","b_d"] = 1;
        dico_k_erreurs["aretes_ajoutees"].append(("b_d","c_f"));
        logger.debug(" * Aretes ajoutees : {}".format( ("b_d","c_f") ))
    elif args["k_erreurs"] == 1 and rd < args["p_correl"]:
        matE_k_alpha.loc["a_b","c_b"] = 0;
        matE_k_alpha.loc["c_b","a_b"] = 0;
        dico_k_erreurs["aretes_supprimees"].append(("a_b","b_c"));
        logger.debug(" * Aretes supprimees : {}".format( ("a_b","b_c") ))
    else:
        aretes_mat = matE_k_alpha.columns.tolist()
        for _ in range(math.ceil(args["k_erreurs"] * args["p_correl"])):         # suppression d'aretes
            irow,row = random.choice(list(enumerate(aretes_mat)));
            aretes_mat.pop(irow);
            icol,col = random.choice(list(enumerate(aretes_mat)));
            aretes_mat.pop(icol);
            
            matE_k_alpha.loc[row,col] = 0;
            matE_k_alpha.loc[col,row] = 0;
            dico_k_erreurs["aretes_supprimees"].append((row,col));
            logger.debug(" * Aretes supprimees : ({},{})".format(row,col))
        aretes_mat = matE_k_alpha.columns.tolist()
        for _ in range(args["k_erreurs"] - math.ceil(args["k_erreurs"] * 
                                                    args["p_correl"])):         # ajout d'aretes
            irow,row = random.choice(list(enumerate(aretes_mat)));
            aretes_mat.pop(irow);
            icol,col = random.choice(list(enumerate(aretes_mat)));
            aretes_mat.pop(icol);
            
            matE_k_alpha.loc[row,col] = 1;
            matE_k_alpha.loc[col,row] = 1;
            dico_k_erreurs["aretes_ajoutees"].append((row,col));
            logger.debug(" * Aretes ajoutees : ({},{})".format(row,col))
    
    # algorithme de couverture
    logger.debug("***** Algorithme de couverture")
    C = list(); aretes_Ec = list();
    dico_cliq = dict(); dico_sommets_par_cliqs = dict();
    dico_gamma_sommets = dict()
    C, dico_cliq, aretes_Ec, ordre_noeuds_traites, \
    dico_sommets_par_cliqs, dico_gamma_sommets = \
    decouvClique.decouverte_cliques_new(matE_k_alpha, dico_arcs_sommets, \
                                    args["seuil_U"], args["epsilon"], \
                                    args["chemin_datasets"], 
                                    args["chemin_matrices"],
                                    args["ascendant_1"], 
                                    args["simulation"],\
                                    dict(),
                                    args)
    print("aretes_Ec={},\n C={} ,\n sommets_par_cliqs={},\n dico_k_erreurs={},\n dict_cliq={}".\
          format(len(aretes_Ec), C, dico_sommets_par_cliqs, dico_k_erreurs, 
                 dico_cliq))
    logger.debug("aretes_Ec={}".format(len(aretes_Ec)));
    logger.debug("C={}".format(len(C))); 
    logger.debug("sommets_par_cliques={}".format(dico_sommets_par_cliqs));
    sommets_1 = {sommet for sommet, etat in dico_cliq if dico_cliq[sommet]==-1}
    logger.debug("sommets_1={}".format(sommets_1))
    
    ### correction
    if len(aretes_Ec) != 0:
        logger.debug("***** Algorithme de correction ")
        C_old = C.copy(); 
        args["C"] = C.copy();
        args["dico_sommets_par_cliqs"] = dico_sommets_par_cliqs;
        args["dico_cliq"] = dico_cliq;
        args["aretes_Ec"] = fct_aux.liste_arcs(matE_k_alpha);
        args["dico_gamma_sommets"] = dico_gamma_sommets;
        dico_solution = algoCorrection.correction_graphe_correlation(args);
        return dico_solution;
    else:
        logger.debug("***** Pas de Correction *****")
        cout_correction = 0; noeuds_corriges = list()
        return {cout_correction:[C, dico_cliq, ordre_noeuds_traites,
                                 dico_sommets_par_cliqs, noeuds_corriges,
                                 C_old]}
#    return matE_k_alpha                                                        # commenter a cause du type de retour.
###############################################################################
#               simulation nouveau critere sur un graphe connu ---> fin
############################################################################### 

###############################################################################
#               simulation de graphes en parallele --> debut
###############################################################################
def simulation_parallele(mat, matE, k, alpha, dico_arcs_sommets, 
                         numero_graphe, args):
    """ simulation d'un reseau avec le nombre d erreurs 
        variant de 0 a k.
    
        * graphe deja genere
        * modification de k cases -> matE_k_alpha
        * algorithme de de decouverte et correction sur G_k_alpha
        * calcul de DH(G_k, G_k_alpha) --> DC (distance de correction)
        * calcul de DH(G, G_k_alpha) ---> DH (distance de Hamming)
        * enregistrement dans un fichier 
            * colonnes -> numero_graphe, k, DC, DH, 
                          nombre_aretes, nombre_cliques_trouve, 
                          nombre_sommets_mat, etc 
    """
    
    logging.basicConfig(format='%(asctime)s - %(message)s', 
                        level=logging.DEBUG,
                        datefmt='%d-%b-%y %H:%M:%S',
                        filename=args["log_file"],
                        filemode="w")
    logger = logging.getLogger('***** simulation_parallele_graphes_generes');
    
    path_distribution = args["dir_base"]+ \
                        args["critere_selection_compression"]+ "/" + \
                        args["mode_correction"]+ "/" + \
                        "data_p_"+str(args["p_correl"]) + \
                        "/distribution";                                        # repertoire des distribution cree
    path_distr = Path(path_distribution);    
    if not path_distr.isdir() :
        path_distr.mkdir(parents=True, exist_ok=True)
        
    G_k = "G_"+str(args["numero_graphe"])+"_"+str(args["k"]);
    aretes_init_matE = fct_aux.liste_arcs(matE.columns.tolist());
    nbre_sommets_matE = len(dico_arcs_sommets.keys())                           # les sommets de matE sont les arcs de mat
    
    chemin_datasets = args["dir_base"]+ \
                        args["critere_selection_compression"]+ "/" + \
                        args["mode_correction"]+ "/" + \
                        "data_p_"+str(args["p_correl"]) + "/" + \
                        G_k + "/" + \
                        args["chemin_datasets"];                                # creation repertoire chemin_datasets
    path_datasets = Path(chemin_datasets);
    if not path_datasets.is_dir():
        path_datasets.mkdir(parents=True, exist_ok=True)
                    
    chemin_matrices = args["dir_base"]+ \
                        args["critere_selection_compression"]+ "/" + \
                        args["mode_correction"]+ "/" + \
                        "data_p_"+str(args["p_correl"]) + "/" + \
                        G_k + "/" + \
                        args["chemin_matrices"];                                # creation repertoire chemin_matrices
    path_matrices = Path(chemin_matrices);
    if not path_matrices.is_dir():
        path_matrices.mkdir(parents=True, exist_ok=True)
        
    matE.to_csv(chemin_matrices+"matE_generer.csv")
    mat.to_csv(chemin_matrices+"mat_generer.csv")

    
    moy_correction = 0; moy_hamming = 0; 
    sum_correction = 0; sum_hamming = 0; 
    correl_dc_dh = 0;
    for alpha in args["alpha"]:
        try :
            print("G_k={}, k={}, alpha={}".format(G_k,k,alpha))
            matE_k_alpha = None;
            dico_k_erreurs = {"aretes_ajoutees":[], "aretes_supprimees":[]};
            dico_proba_cases = dict();
            matE_k_alpha, dico_k_erreurs, dico_proba_cases = \
                    modifier_k_cases(matE, 
                                     args["p_correl"], 
                                     args["k"], 
                                     args["loi_proba"], 
                                     args["val_limite_proba_1"])
            matE_k_alpha.to_csv(chemin_matrices +
                                "matE_" +
                                str(k) + "_" +
                                str(alpha) +
                                ".csv")
            
            # algorithme de couverture
            logger.debug("***** Algorithme de couverture :")
            C = list(); aretes_Ec = list();
            dico_cliq = dict(); dico_sommets_par_cliqs = dict();
            dico_gamma_sommets = dict()
            C, dico_cliq, aretes_Ec, ordre_noeuds_traites, \
            dico_sommets_par_cliqs, dico_gamma_sommets = \
            decouvClique.decouverte_cliques_new(matE_k_alpha, \
                                                dico_arcs_sommets, \
                                                args["seuil_U"], 
                                                args["epsilon"], \
                                                chemin_datasets, 
                                                chemin_matrices,
                                                args["ascendant_1"], 
                                                args["simulation"],\
                                                dict(),
                                                args)
            print("aretes_Ec={},\n C={} ,\n sommets_par_cliqs={},\n dico_k_erreurs={},\n dict_cliq={}".\
                  format(len(aretes_Ec), C, dico_sommets_par_cliqs, dico_k_erreurs, 
                         dico_cliq))
            logger.debug("aretes_Ec={}".format(len(aretes_Ec)));
            logger.debug("C={}".format(len(C))); 
            logger.debug("sommets_par_cliques={}".format(dico_sommets_par_cliqs));
            sommets_1 = {sommet for sommet, etat in dico_cliq if dico_cliq[sommet]==-1}
            logger.debug("sommets_1={}".format(sommets_1))
            
            ### correction
            dico_sommets_par_cliqs_new = dico_sommets_par_cliqs.copy();
            dico_solution = dict(); 
            args_res = dict();
            if len(aretes_Ec) != 0:
                logger.debug("***** Algorithme de correction ")
                C_old = C.copy(); 
                args["C"] = C.copy();
                args["dico_sommets_par_cliqs"] = dico_sommets_par_cliqs;
                args["dico_cliq"] = dico_cliq;
                args["aretes_Ec"] = fct_aux.liste_arcs(matE_k_alpha);
                args["dico_gamma_sommets"] = dico_gamma_sommets;
                args_res, dico_solution = \
                        algoCorrection.correction_graphe_correlation(args);
            else:
                logger.debug("***** Pas de Correction *****")
            
            cliques_identiques_C_C_old = list();   
            cliques_differentes_C_C_old = list();
            dc = list();
            dh = list();
            if args_res:
                dc, aretes_diff_dc = distance_hamming(
                            set(fct_aux.liste_arcs(matE_k_alpha.columns.tolist())), 
                            args_res["aretes_Ec"])
                dh, aretes_diff_dh = distance_hamming(set(aretes_init_matE),
                                                      args_res["aretes_Ec"])
                sum_correction += dc;
                sum_hamming += dh;
                cliques_identiques_C_C_old, cliques_differentes_C_C_old = \
                    comparer_cliques(args_res["C"], C.copy())                  # C.copy() = C_old
                dico_sommets_par_cliqs_new = args_res["dico_sommets_par_cliqs_new"]
            df_dico = dict();
            df_dico["G_k"] = G_k; df_dico["k"] = k; df_dico["alpha"] = alpha;
            df_dico["nbre_sommets_matE"] = nbre_sommets_matE;
            df_dico["aretes_init_matE"] = aretes_init_matE; 
            df_dico["aretes_ajoutees"] = dico_k_erreurs["aretes_ajoutees"]; 
            df_dico["aretes_supprimees"] = dico_k_erreurs["aretes_supprimees"]; 
            df_dico["dc"] = dc; df_dico["dh"] = dh; 
            df_dico["C_old"] = len(C_old); 
            df_dico["C"] = len(args["C"]);
            df_dico["cliques_identiques_C_C_old"] = len(cliques_identiques_C_C_old);
            df_dico["cliques_differentes_C_C_old"] = len(cliques_differentes_C_C_old); 
            #df_dico[""] = ; df_dico[""] = ;
            for cpt_sommet, value in dico_solution.items():
                df_dico["etape_"+str(cpt_sommet[0])+"_sommet_1"] = \
                                    cpt_sommet[1]
                df_dico["etape_"+str(cpt_sommet[0])+"_sommets_corriges"] = \
                                    len(value["sommets_corriges"])
                                    
                df_dico["etape_"+str(cpt_sommet[0])+"_nbre_aretes_ajoutes_p1"]=\
                                    len(value["cout_T"]["aretes_ajoutes_p1"])
                df_dico["etape_"+str(cpt_sommet[0])+"_aretes_ajoutes_p1"] = \
                                    value["cout_T"]["aretes_ajoutes_p1"]
                df_dico["etape_"+str(cpt_sommet[0])+"_aretes_p1"] = \
                                list(it.combinations(value["compression_p1"],2))
                                    
                df_dico["etape_"+str(cpt_sommet[0])+"_aretes_ajoutes_p2"] = \
                                    value["cout_T"]["aretes_ajoutes_p2"]
                df_dico["etape_"+str(cpt_sommet[0])+"_nbre_aretes_ajoutes_p2"]=\
                                    len(value["cout_T"]["aretes_ajoutes_p2"])
                df_dico["etape_"+str(cpt_sommet[0])+"_aretes_p2"] = \
                                list(it.combinations(value["compression_p2"],2))
                                    
                df_dico["etape_"+str(cpt_sommet[0])+"_aretes_supprimes"] = \
                                    value["cout_T"]["aretes_supprimes"]
                df_dico["etape_"+str(cpt_sommet[0])+"_nbre_aretes_supprimes"] = \
                                    len(value["cout_T"]["aretes_supprimes"])
                                    
                df_dico["etape_"+str(cpt_sommet[0])+"_min_c1"] = \
                                    value["cout_T"]["min_c1"]
                df_dico["etape_"+str(cpt_sommet[0])+"_min_c2"] = \
                                    value["cout_T"]["min_c2"]
            # mettre un for pour dico_sommets_par_cliqs_new
            for sommet, cliques in dico_sommets_par_cliqs_new.items():
                df_dico[str(sommet)] = len(cliques);
            # convertir df_dico en dataframe
            df = pd.DataFrame.from_dict(df_dico, orient="index");
            df.columns = [G_k];
            # save dataframe
            name_save_df = args["dir_base"]+ \
                            args["critere_selection_compression"]+ "/" + \
                            args["mode_correction"]+ "/" + \
                            "data_p_"+str(args["p_correl"]) + \
                            "/distribution" + "/" + \
                            "resumeExecution_"+ \
                            str(k) + \
                            ".csv";
            sauver_df_resume(df, name_save_df, G_k);   
        except Exception as e:
            df_dico["G_k"] = G_k; df_dico["k"] = k; df_dico["alpha"] = alpha;
            df_dico["nbre_sommets_matE"] = nbre_sommets_matE;
            df_dico["aretes_init_matE"] = aretes_init_matE; 
            df_dico["aretes_ajoutees"] = "error";
            df_dico["aretes_supprimees"] = "error"; 
            df_dico["aretes_matE_k_alpha"] = "error";
            df_dico["dc"] = "error"; df_dico["dh"] = "error"; 
            df_dico["C_old"] = "error"; 
            df_dico["C"] = "error";
            df_dico["cliques_identiques_C_C_old"] = "error";
            df_dico["cliques_differentes_C_C_old"] = "error";
    # moyenner dist_line et hamming pour k aretes supprimes
    moy_correction = sum_correction / alpha_max;
    moy_hamming = sum_hamming / alpha_max;
    if moy_hamming == 0 and moy_correction == 0:
        correl_dc_dh = 1
    else:
        correl_dc_dh = abs(moy_hamming - moy_correction) / \
                        max(moy_hamming, moy_correction)
    
    # ecrire dans un fichier pouvant etre lu pendant qu'il continue d'etre ecrit
    f = open(path_distribution + 
             "distribution_moyDistLine_moyHamming_k_" +
             str(k) + 
             ".txt","a")
    f.write(G_k + ";" +\
            str(k) + ";" + \
            str(moy_correction) + ";" + \
            str(moy_hamming) + ";" + \
            str(aretes_init_matE) + ";" + \
            str(correl_dc_dh) + "\n")
    f.close();        


###############################################################################
#               simulation de graphes en parallele --> fin
###############################################################################



if __name__ == '__main__':
    
     start= time.time();
     args = {"dbg":True};
     log_file = "DEBUG_CRITERE_C2C1.log";
     log_simulation = "DEBUG_simulation_parallele";
     
     bool_reseau = False;
     bool_couverture_graphe_connu = False;
     bool_simulation = False;
     bool_parallele = False;
     bool_test_critere_correction = True;
     
     # valeurs des parametres de decouverte_cliques
     seuil_U = 10; epsilon = 0.75; effet_joule = 0;
     nbre_ts = 10; nbre_lien = (2,5);
     ascendant_1 = True;
     simulation = True;
     k_erreurs = 1; 
     p_correl = 0.5;
     loi_proba = "gaussienne";
     val_limite_proba_1 = 0.8
     number_items_pi1_pi2 = 1;
     mode_correction = "aleatoire_sans_remise";                                 # "aleatoire_sans_remise", degre_min_sans_remise, cout_min_sans_remise, aleatoire_avec_remise", degre_min_avec_remise, cout_min_avec_remise, 
     critere_selection_compression = "voisins_corriges"                         # "voisins_corriges" (C2), "nombre_aretes_corriges" (C1), "voisins_nombre_aretes_corriges" (C2 puis C1)
     args = {"dbg":True, "seuil_U":seuil_U, "epsilon":epsilon, 
             "ascendant_1":ascendant_1, "simulation":simulation,
             "k_erreurs":k_erreurs, 
             "p_correl":p_correl,
             "loi_proba": loi_proba, 
             "val_limite_proba_1": val_limite_proba_1,
             "number_items_pi1_pi2":number_items_pi1_pi2,
             "mode_correction":mode_correction,
             "critere_selection_compression":critere_selection_compression,
             "log_file":log_file, 
             "log_simulation":log_simulation}
    
     # generation reseau energetique ==>OK
     dim_mat = 5;
     nbre_graphes = 10;
     graphes = list();
     if bool_reseau:
         chemin_datasets = "dataNewCriterecorrectionGrapheConnu/datasets/";
         chemin_matrices = "dataNewCriterecorrectionGrapheConnu/matrices/";
         dir_base = "dataNewCriterecorrectionGrapheConnu/";
         args["chemin_datasets"] = chemin_datasets;
         args["chemin_matrices"] = chemin_matrices;
         args["dir_base"] = dir_base;
         matE, mat, dico_arcs_sommets = creer_reseau(chemin_datasets, 
                                                 chemin_matrices, args);
     else:
         chemin_datasets = "datasets/";
         chemin_matrices = "matrices/";
         args["chemin_datasets"] = chemin_datasets;
         args["chemin_matrices"] = chemin_matrices;
         dir_base = "dataNewCriterecorrection/";
         args["dir_base"] = dir_base;
         for numero_graphe in range(nbre_graphes):
             matE, mat, dico_arcs_sommets = generer_reseau(dim_mat, nbre_lien,
                                                           chemin_datasets,
                                                           chemin_matrices, 
                                                           nbre_ts, epsilon, 
                                                           effet_joule)
             graphes.append((matE, mat, dico_arcs_sommets, numero_graphe))
                                                         
     
     # test sur couverture en cliques ===> 
     if bool_couverture_graphe_connu:
         dico_solution = simulation_nouveau_critere(args);                                           
     
     if bool_simulation:
         alpha_max = 1;
         k_min = 0; k_max = 5; step_range = 1;
         k_range = range(k_min, k_max, step_range);
         p_correl_max = 1;
         p_correl_min = 0;
         step_range_p = 0.1;
         
         
         modes_correction = list();
         criteres_selection_compression = list();
         if bool_test_critere_correction:
             modes_correction = ["aleatoire_sans_remise"]
             criteres_selection_compression = ["voisins_corriges"]
             p_correls = [0.5];
         else:
             modes_correction = ["aleatoire_sans_remise", 
                                 "degre_min_sans_remise", 
                                 "cout_min_sans_remise", 
                                 "aleatoire_avec_remise", 
                                 "degre_min_avec_remise", 
                                 "cout_min_avec_remise"]
             criteres_selection_compression = ["voisins_corriges", 
                                               "nombre_aretes_corriges", 
                                               "voisins_nombre_aretes_corriges"]
             p_correls = range(p_correl_min, p_correl_max, step_range_p);                                 
         params = list();
         for tupl in it.product(p_correls, 
                                modes_correction, 
                                criteres_selection_compression,
                                graphes,
                                k_range,
                                [alpha_max]):
             args["p_correl"] = tupl[0];
             args["mode_correction"] = tupl[1];
             args["critere_selection_compression"] = tupl[2];
             matE = tupl[3][0];
             mat = tupl[3][1];
             dico_arcs_sommets = tupl[3][2];
             numero_graphe = tupl[3][3];
             k = tupl[4]
             alpha = tupl[5]; 
             params.append( (mat, matE, k, alpha, 
                             dico_arcs_sommets, numero_graphe, 
                             args.copy()) );
             
         print("params = {}".format(len(params)))
         
         #parallelisation avec multiprocessing
         if bool_parallele:
             p = Pool(mp.cpu_count()-1) 
             p.starmap(simulation_parallele, params)
             p.terminate()
         else:
             arguments = params[0];
             simulation_parallele(arguments[0],
                                  arguments[1],
                                  arguments[2],
                                  arguments[3],
                                  arguments[4],
                                  arguments[5],
                                  arguments[6])
         #### parallele
         
         g = open("tempsExecution_SIMULATION_k"+
                  "_".join(map(str, k_range))+
                  ".txt","w")
         ti = time.time() - start
         g.write(str(ti))
         g.close()