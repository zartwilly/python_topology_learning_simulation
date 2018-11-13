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
    
    path_datasets = Path(args["chemin_datasets"]);
    if not path_datasets.is_dir():
        path_datasets.mkdir(parents=True, exist_ok=True)
        
    path_matrices = Path(args["chemin_matrices"]);
    if not path_matrices.is_dir():
        path_matrices.mkdir(parents=True, exist_ok=True)
    
   
    path_distribution = args["dir_base"]+ \
                        args["critere_selection_compression"]+ "/" + \
                        args["mode_correction"]+ "/" + \
                        "data_p_"+str(args["p_correl"]) + \
                        "/distribution";                                        # repertoire des distribution cree
    path_distr = Path(path_distribution);    
    if not path_distr.isdir() :
        path_distr.mkdir(parents=True, exist_ok=True)
        
    
    G_k = "G_"+str(args["numero_graphe"])+str(args["k"]);
    
    

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
     bool_test_critere_correction = True;
     
     # valeurs des parametres de decouverte_cliques
     seuil_U = 10; epsilon = 0.75; effet_joule = 0;
     nbre_ts = 10; nbre_lien = (2,5);
     ascendant_1 = True;
     simulation = True;
     k_erreurs = 1; 
     p_correl = 0.5;
     number_items_pi1_pi2 = 1;
     mode_correction = "aleatoire_sans_remise";                                 # "aleatoire_sans_remise", degre_min_sans_remise, cout_min_sans_remise, aleatoire_avec_remise", degre_min_avec_remise, cout_min_avec_remise, 
     critere_selection_compression = "voisins_corriges"                         # "voisins_corriges" (C2), "nombre_aretes_corriges" (C1), "voisins_nombre_aretes_corriges" (C2 puis C1)
     args = {"dbg":True, "seuil_U":seuil_U, "epsilon":epsilon, 
             "ascendant_1":ascendant_1, "simulation":simulation,
             "k_erreurs":k_erreurs, 
             "p_correl":p_correl,
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
         chemin_datasets = "dataNewCriterecorrection/datasets/";
         chemin_matrices = "dataNewCriterecorrection/matrices/";
         dir_base = "dataNewCriterecorrection/";
         args["chemin_datasets"] = chemin_datasets;
         args["chemin_matrices"] = chemin_matrices;
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
         p = Pool(mp.cpu_count()-1) 
         p.starmap(simulation_parallele, params)
         p.terminate()
         #### parallele
         
         g = open("tempsExecution_SIMULATION_k"+
                  "_".join(map(str, k_range))+
                  ".txt","w")
         ti = time.time() - start
         g.write(str(ti))
         g.close()