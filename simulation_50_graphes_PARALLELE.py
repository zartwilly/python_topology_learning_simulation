#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 10:09:53 2017

export PATH=/home/dcbrain/anaconda3/bin:$PATH
scp  *.py dcbrain@wil.dcbra.in:/home/dcbrain/simulation_refactoring_withTryExcept/

@author: willy
"""
import re, os, time
import pandas as pd
import numpy as np
import fonctions_auxiliaires as fct_aux
import generations_mesures as mesures
import decouverte_cliques as decouvClique
import genererMatA as geneMatA

#import matplotlib.pyplot as plt
from pathlib import Path    # http://stackoverflow.com/questions/6004073/how-can-i-create-directories-recursively
import itertools as it
from multiprocessing import Pool;
from pandas.parser import CParserError;
import multiprocessing as mp
from scipy.stats import truncnorm
from random import randint as valAleotoire
import random;


def distance_hamming(edge_list_matE, edge_list_matE_LG): 
    """
    but: comparer les 2 listes d'aretes et compter le nombre d'aretes differentes entre eux.
    Ce nombre est le nombre de HAMMING
    
    nbre_hamming : nbre d'elements differents entre matE et matE_
    liste_cpt : liste d'aretes ou arcs  n'appartenant pas aux 2 matrices matE et matE_
    """
    liste_arc_diff = list(); set_arc_diff = set();
    for arete in edge_list_matE:
        """
        case de 1 a 0 dans matE_LG
        aretes ABSENTE(s) de matE_LG mais PRESENCE(s) dans matE
        """
        if (arete[0], arete[1]) not in edge_list_matE_LG and \
            (arete[1], arete[0]) not in edge_list_matE_LG:
                liste_arc_diff.append( (arete[0], arete[1]) )
                set_arc_diff.add( (arete[0], arete[1]) )
    for arete in edge_list_matE_LG:
        """
        case de 1 a 0 dans matE
        aretes ABSENTE(s) de matE mais PRESENCE(s) dans matE_LG
        """
        if (arete[0], arete[1]) not in edge_list_matE and \
            (arete[1], arete[0]) not in edge_list_matE:
                liste_arc_diff.append( (arete[0], arete[1]) )
                set_arc_diff.add( (arete[0], arete[1]) )
    
#    return len(liste_arc_diff), liste_arc_diff;
    return len(set_arc_diff), set_arc_diff;

def select_items_df( nbre_1, matE, list_excluded_rows_cols):
    '''
    but selectionner aleatoirement la ligne et colonnes d'un elt dont la valeur est egale a nbre_1
        ds matE
    Particularite: NE PAS CHOISIR 2 FOIS LE MEME ELEMENT
    return row, col, list_excluded_rows_cols
    '''
    bol = True
    col_matE = matE.columns.tolist()
    len_col_matE = len(col_matE)
    while bol :
        ind_row = np.random.randint(len_col_matE)
        ind_col = np.random.randint(len_col_matE)
        row = col_matE[ind_row]
        col = col_matE[ind_col]
        if matE.loc[row][col] == nbre_1 and (row,col) not in list_excluded_rows_cols:
            bol = False
            list_excluded_rows_cols.append( (row,col) )
            list_excluded_rows_cols.append( (col,row) )
            return row, col, list_excluded_rows_cols

def supprimer_k_aretes(matE, k):
    """
    but: supprimer k aretes de tel sorte qu'il est un nbre proportionnel 
         de 0 et 1 dans la matE resultante (matE_resul)
    n1: nbre de 1 dans matE
    n0: nbre de 0 dans matE
    P: proba de reference pour le nbre de 1 dans matE
    """
    deleted_edges = list()
    if k == 0:
        return matE, deleted_edges 
        
    n1 = 0
    n0 = 0
    list_excluded_rows_cols = list() # [(row,col),(col,row),....]
    col_matE = matE.columns.tolist()
    for row in col_matE:
        for col in col_matE:
            if matE.loc[row][col] == 1 :
                n1 += 1
            else:
                n0 += 1
    p = 0
    p = 1/n1 if n1 > n0 else 1/n0
    bol = True
    cpt_arete_supp = 0
    while bol:
        ps = np.random.rand()
        if ps > p : # auparavant cetait ecrit ps< p
            # changement de 1 en 0
            #print("p= ",round(p,3)," ps:",round(ps,3), " ",cpt_arete_supp+1," arete SUPPRIME")
            nbre_1 = 1
            row, col, list_excluded_rows_cols = select_items_df( nbre_1, matE, list_excluded_rows_cols)
            matE.loc[row][col] = 0
            matE.loc[col][row] = 0
            deleted_edges.append( (row, col) )
            cpt_arete_supp += 1
        bol = False if cpt_arete_supp == k else True 
    return matE, deleted_edges

###### modif k cases dans  matE ======> debut
def delete_aretes(j,liste_aretes_0_1,dico_deleted_add_edges, operation):
    """
    suppression de l'arete a l'indice j puis ajout dans le dico_deleted_add_edges
    soit avec la cle "ajouter" ou "supprimer"
    """
#    if operation == 0:
#        # case 0 -> 1: case from 0 to 1
#        if "ajouter" not in dico_deleted_add_edges:
#            dico_deleted_add_edges["ajouter"] = [liste_aretes_0_1[j]]
#        else:
#            dico_deleted_add_edges["ajouter"].append(liste_aretes_0_1[j])
#        del liste_aretes_0_1[j]   
#    elif operation == 1:
#        # case 1 -> 0: case from 1 to 0
#        if "supprimer" not in dico_deleted_add_edges:
#            dico_deleted_add_edges["supprimer"] = [liste_aretes_0_1[j]]
#        else:
#            dico_deleted_add_edges["supprimer"].append(liste_aretes_0_1[j])
#        del liste_aretes_0_1[j] 
    if operation == 0:
        # case 0 -> 1: case from 0 to 1
        dico_deleted_add_edges["ajouter"].append(liste_aretes_0_1[j])
        del liste_aretes_0_1[j]   
    elif operation == 1:
        # case 1 -> 0: case from 1 to 0
        dico_deleted_add_edges["supprimer"].append(liste_aretes_0_1[j])
        del liste_aretes_0_1[j] 
    
    return liste_aretes_0_1, dico_deleted_add_edges   
    
def ajout_suppression_arete(proba, liste_aretes_0, liste_aretes_1, k_cases):
    """
    selection de cases et modification
    """
    dico_deleted_add_edges = dict();
    dico_deleted_add_edges["ajouter"] = []; dico_deleted_add_edges["supprimer"] = [];
    for k in range(0, k_cases):
        if proba < 0.5:
            # case 0 -> 1
            i = random.choice(range(0, len(liste_aretes_0)))
            liste_aretes_1.append(liste_aretes_0[i])
            liste_aretes_0, dico_deleted_add_edges = delete_aretes(i, liste_aretes_0,\
                                                                   dico_deleted_add_edges, 0)
        else:
            # case 1 -> 0
            j = random.choice(range(0, len(liste_aretes_1)))
            liste_aretes_0.append(liste_aretes_1[j])
            liste_aretes_1, dico_deleted_add_edges = delete_aretes(j, liste_aretes_1,\
                                                                   dico_deleted_add_edges, 1)
    return liste_aretes_0, liste_aretes_1, dico_deleted_add_edges
    pass
        
def modif_k_cases(matE, k_cases, methode, p_correl_seuil):
    """
    changer k cases dans matE tel que les cases
        * 0 passent a 1 ( FAUX POSITIFS )
        * 1 passent a 0 ( FAUX NEGATIFS )
    le but est de reduire le biais cree par la suppression d aretes dans matE
    methode 0: avec proba
    methode 1: avec liste
    p_correl_seuil = 0.5, 0.75, 1 # si p_correl_seuil = 0 => on ajoute que des aretes
                                  # si p_correl_seuil = 1 => on supprime que des aretes
                                  # si p_correl_seuil = 0.5 => on ajoute et supprime des aretes
    explications:
        0           faux neg     faux pos      1
        ---------------------|----------------------------> p_correl
        0->0        1->0     |    0->1       1->1
                             |
         suppresion aretes   |     ajout aretes
                             |
                        p_correl_seuil
    
    test methode:
        dimMat = 5;nbre_lien = 5;nbre_ts = 10
        effet_joule = 0;epsilon = 0.75;test = "FINI";
        chemin_datasets = "data/G1/datasets/";chemin_matrices = "data/G1/matrices/";
        methode = 0
        matE, df_matA, dico_arc = simu50.matriceE(dimMat, nbre_lien, chemin_matrices,\
                                                  chemin_datasets, nbre_ts, epsilon, \
                                                  effet_joule, test)
        matE_, dico_deleted_add_edges= modif_k_cases(matE.copy(), 2, methode)
        DH, l_DH = simu50.distance_hamming(fct_aux.liste_arcs(matE),fct_aux.liste_arcs(matE_))
        print("DH = ", DH," l_dh = ", l_DH," dico_deleted_add_edges = ", dico_deleted_add_edges.values())
    """
    if methode == 0:
        # methode avec biais
        dico_deleted_add_edges = dict();
        dico_deleted_add_edges["ajouter"] = []; 
        dico_deleted_add_edges["supprimer"] = [];
        
        for k in range(0, k_cases):
            proba = random.random();
        
            liste_aretes_0 = fct_aux.liste_nonArcs(matE, 0);
            liste_aretes_1 = fct_aux.liste_arcs(matE);
            
            m0 = len(liste_aretes_0);
            m1 = pow(matE.shape[0],2) - m0;    
                     
            if proba < p_correl_seuil:
                # case 0 -> 1
                p0 = 1/m0;
                bool = True;
                while bool:
                    arete = random.choice(liste_aretes_0)
                    if random.random() > p0 and \
                        arete not in dico_deleted_add_edges["ajouter"] and \
                        arete not in dico_deleted_add_edges["supprimer"] :
                        matE.loc[arete[0]][arete[1]] = 1; matE.loc[arete[1]][arete[0]] = 1;
                        dico_deleted_add_edges["ajouter"].append(arete) 
                        bool = False
                        break;
#                for arete in liste_aretes_0 :
#                    if random.random() > p0 and arete not in dico_deleted_add_edges["ajouter"] \
#                        and arete not in dico_deleted_add_edges["supprimer"] :
#                        matE.loc[arete[0]][arete[1]] = 1; matE.loc[arete[1]][arete[0]] = 1;
#                        dico_deleted_add_edges["ajouter"].append(arete)   
#                        break;
            else:
                # case 1 -> 0
                p1 = 1/m1;
                bool = True;
                while bool:
                    arete = random.choice(liste_aretes_1)
                    if random.random() > p1 and \
                        arete not in dico_deleted_add_edges["ajouter"] and \
                        arete not in dico_deleted_add_edges["supprimer"]:
                        matE.loc[arete[0]][arete[1]] = 0; matE.loc[arete[1]][arete[0]] = 0;
                        dico_deleted_add_edges["supprimer"].append(arete)   
                        bool = False;
                        break;   
#                for arete in liste_aretes_1: 
#                    if random.random() > p1 and arete not in dico_deleted_add_edges["ajouter"] \
#                        and arete not in dico_deleted_add_edges["supprimer"]:
#                        matE.loc[arete[0]][arete[1]] = 0; matE.loc[arete[1]][arete[0]] = 0;
#                        dico_deleted_add_edges["supprimer"].append(arete)   
#                        break;
        return matE, dico_deleted_add_edges;
    else:
        # methode avec liste
        liste_aretes_0 = fct_aux.liste_nonArcs(matE, 0);
        liste_aretes_1 = fct_aux.liste_arcs(matE);
        proba = random.random();
        dico_deleted_add_edges = dict();
        liste_aretes_0, liste_aretes_1, dico_deleted_add_edges = \
        ajout_suppression_arete(proba, liste_aretes_0, liste_aretes_1, k_cases)
        
        for key, aretes in dico_deleted_add_edges.items():
            if key == "ajouter":
                for arete in aretes:
                    matE.loc[arete[0]][arete[1]] = 1; matE.loc[arete[1]][arete[0]] = 1;
            elif key == "supprimer":
                for arete in aretes:
                    matE.loc[arete[0]][arete[1]] = 0; matE.loc[arete[1]][arete[0]] = 0; 
                    
        return matE, dico_deleted_add_edges;
###### modif k cases dans  matE ======> fin
    
##### ajout proba de chaque case de matE ====> debut
def loi_proba(x,s):
    return pow(x-s, 2);
def loi_de_probalibilite_old(debut_prob, fin_proba, nbre_correlations, correl_seuil):
    """
    genere des valeurs compris [debut_prob, fin_proba] et leur proba associe
    """
    correlations = np.linspace(debut_prob, fin_proba, nbre_correlations)
    proba_correlations =  [loi_proba(x,correl_seuil) for x in correlations]
    # Normalising to 1.0
    proba_correlations /= np.sum(proba_correlations)
    return round(np.random.choice(correlations, 1, p=proba_correlations)[0], 3)
    
def get_truncated_normal(mean=0, sd=1, low=0, upp=10):
    return truncnorm(
        (low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)
    
def loi_de_probalibilite(debut_proba, fin_proba, nbre_correlations, correl_seuil):
    """
    generer des valeurs de correlations entre debut_proba et fin_proba.
    """
    mean = (fin_proba + debut_proba)/2; 
    sd = abs((fin_proba - debut_proba)/(2*4)) if fin_proba != debut_proba else abs((fin_proba)/(2*4))
    correlations  = get_truncated_normal(mean*100, sd*100, low=1, upp = nbre_correlations*100);
    return np.random.choice(correlations.rvs(nbre_correlations)/100);

###############################################################################
#######################  A REVOIR    ##########################################
def ajouter_proba_matE_old_old(matE, dico_deleted_add_edges, loi_stats, p_correl, correl_seuil=0.8):
    """
    le but est de definir des probabilites pour chaque case de matE pour simuler 
    la matrice de correlation obtenue apres calcul. les probas (correls) sont les suivantes:
        * case 0 -----> 0: proba entre 0 - 0.5
        * case 1 -----> 0: proba entre 0.6 - 0.79 ( X-0.2 - X-0.1 ) ==> faux negatifs
        * case 0 -----> 1: proba a 0.8 = X                          ==> faux positifs
        * case 1 -----> 1: proba entre 0.8 - 1
        
        * case 0 -----> 0: proba entre 0 - 0.5 (0, correl_seuil)                   ==> vrai negatifs
        * case 1 -----> 0: proba entre 0.6 - 0.79 ( p_correl-0.2 - p_correl-0.01 ) ==> faux negatifs
        * case 0 -----> 1: proba entre p_correl et correl_seuil -0.001             ==> faux positifs
        * case 1 -----> 1: proba entre 0.8 - 1 (correl_seuil+0.001, 1)             ==> vrai positifs
    NB: X = correl_seuil
    dico_deleted_add_edges = {"ajouter":[(a,b),...],"supprimer":[(c,d),...]}
    loi_stats = loi uniforme, de poisson
    """
    cases_traites = list(); dico_proba_case = dict(); nbre_correlations = 250
    for row in matE.columns.tolist():
        for col in matE.columns.tolist():
            if col != row and (row,col) not in cases_traites and (row,col) not in cases_traites:
#                print("ici  ==> (row,col)=",(row,col))
#                print("ici matE.loc[row][col]=", matE.loc[row][col])
#                print("dico_deleted_add_edges = ", dico_deleted_add_edges)
                if col != row and matE.loc[row][col] == 0 :
                    if (row,col) in dico_deleted_add_edges["supprimer"] or \
                        (col,row) in dico_deleted_add_edges["supprimer"]:
                        # proba entre 0.6 et 0.79 (entre X-0.2 -- X-0.01) ===> faux negatifs (distributions uniforme)
#                        print("ici 1 ==> (row,col)=",(row,col))
                        if loi_stats == "poisson":
                            dico_proba_case[(row,col)] = round(random.choice(np.union1d(np.linspace(0.6,0.644,100), \
                                                            np.linspace(0.72,0.79,100))),3)
                        else:
                            correlation = loi_de_probalibilite(p_correl-0.2, p_correl-0.01, nbre_correlations, correl_seuil)
#                            dico_proba_case[(row,col)] = round(random.uniform(correl_seuil-0.2, correl_seuil-0.01), 3);
                            dico_proba_case[(row,col)] = correlation;
                        cases_traites.append( (row,col) )
                    else:
                        # proba entre 0 - 0.5 ( 0 -- X-0.3)  ===> vrai negatifs (distributions asymetrique selon fonction loi_proba)
#                        print("ici 2 ==> (row,col)=",(row,col))
                        if loi_stats == "poisson":
                            dico_proba_case[(row,col)] = round(random.choice(np.union1d(np.linspace(0,0.117,100), \
                                                            np.linspace(0.35,0.5,100))),3)
                        else:
                            correlation = loi_de_probalibilite(0, correl_seuil, nbre_correlations, correl_seuil)
#                            dico_proba_case[(row,col)] = round(random.uniform(0, correl_seuil-0.3), 3);
                            dico_proba_case[(row,col)] =  correlation;              
                        cases_traites.append( (row,col) )
                elif col != row and matE.loc[row][col] == 1 :
                    if (row,col) in dico_deleted_add_edges["ajouter"] or \
                        (col,row) in dico_deleted_add_edges["ajouter"]:
                        # proba entre 0.8 = X ===> faux positifs (distributions uniforme selon fonction loi_proba)
#                        print("ici 3 ==> (row,col)=",(row,col))
                        correlation = loi_de_probalibilite(p_correl, correl_seuil -0.001, nbre_correlations, correl_seuil)
#                        dico_proba_case[(row,col)] = correl_seuil;
                        dico_proba_case[(row,col)] = correlation
                        cases_traites.append( (row,col) );
                    else:
                        # proba entre 0.8001 - 1 ( X -- 1 ) ===> vrai positifs (distributions asymetrique selon fonction loi_proba)
#                        print("ici 4 ==> (row,col)=",(row,col))
                        correlation = loi_de_probalibilite(correl_seuil+0.001, 1, nbre_correlations, correl_seuil)
#                        dico_proba_case[(row,col)] = round(random.uniform(correl_seuil+0.001, 1), 3);
                        dico_proba_case[(row,col)] = correlation
                        cases_traites.append( (row,col) );
    return dico_proba_case;               
def ajouter_proba_matE_old(matE, dico_deleted_add_edges, loi_stats, p_correl, correl_seuil=0.8):
    """
    le but est de definir des probabilites pour chaque case de matE pour simuler 
    la matrice de correlation obtenue apres calcul. les probas (correls) sont les suivantes:
        * case 0 -----> 0: proba entre 0 - 0.5
        * case 1 -----> 0: proba entre 0.6 - 0.79 ( X-0.2 - X-0.1 ) ==> faux negatifs
        * case 0 -----> 1: proba a 0.8 = X                          ==> faux positifs
        * case 1 -----> 1: proba entre 0.8 - 1
        
        * case 0 -----> 0: proba entre 0 - 0.5 (0, correl_seuil-p_correl-0.01)           ==> vrai negatifs
        * case 1 -----> 0: proba entre 0.6 - 0.79 ( correl_seuil-p_correl, correl_seuil) ==> faux negatifs
        * case 0 -----> 1: proba entre p_correl et correl_seuil-0.001 (p_correl,correl_seuil-0.01)==> faux positifs
        * case 1 -----> 1: proba entre 0.8 - 1 (correl_seuil, 1)             ==> vrai positifs
    NB: X = correl_seuil
    dico_deleted_add_edges = {"ajouter":[(a,b),...],"supprimer":[(c,d),...]}
    loi_stats = loi uniforme, de poisson
    """
    cases_traites = list(); dico_proba_case = dict(); nbre_correlations = 250
    for row in matE.columns.tolist():
        for col in matE.columns.tolist():
            if col != row and (row,col) not in cases_traites and (row,col) not in cases_traites:
#                print("ici  ==> (row,col)=",(row,col))
#                print("ici matE.loc[row][col]=", matE.loc[row][col])
#                print("dico_deleted_add_edges = ", dico_deleted_add_edges)
                if col != row and matE.loc[row][col] == 0 :
                    if (row,col) in dico_deleted_add_edges["supprimer"] or \
                        (col,row) in dico_deleted_add_edges["supprimer"]:
                        # proba entre 0.6 et 0.79 (entre X-0.2 -- X-0.01) ===> faux negatifs (distributions uniforme)
#                        print("ici 1 ==> (row,col)=",(row,col))
                        if loi_stats == "poisson":
                            dico_proba_case[(row,col)] = round(random.choice(np.union1d(np.linspace(0.6,0.644,100), \
                                                            np.linspace(0.72,0.79,100))),3)
                        else:
                            correlation = loi_de_probalibilite(correl_seuil-p_correl, p_correl, \
                                                               nbre_correlations, correl_seuil)
#                            dico_proba_case[(row,col)] = round(random.uniform(correl_seuil-0.2, correl_seuil-0.01), 3);
                            dico_proba_case[(row,col)] = correlation;
                        cases_traites.append( (row,col) )
                    else:
                        # proba entre 0 - 0.5 ( 0 -- X-0.3)  ===> vrai negatifs (distributions asymetrique selon fonction loi_proba)
#                        print("ici 2 ==> (row,col)=",(row,col))
                        if loi_stats == "poisson":
                            dico_proba_case[(row,col)] = round(random.choice(np.union1d(np.linspace(0,0.117,100), \
                                                            np.linspace(0.35,0.5,100))),3)
                        else:
                            correlation = loi_de_probalibilite(0, correl_seuil-p_correl-0.01, \
                                                               nbre_correlations, correl_seuil);
#                            dico_proba_case[(row,col)] = round(random.uniform(0, correl_seuil-0.3), 3);
                            dico_proba_case[(row,col)] =  correlation;              
                        cases_traites.append( (row,col) )
                elif col != row and matE.loc[row][col] == 1 :
                    if (row,col) in dico_deleted_add_edges["ajouter"] or \
                        (col,row) in dico_deleted_add_edges["ajouter"]:
                        # proba entre 0.8 = X ===> faux positifs (distributions uniforme selon fonction loi_proba)
#                        print("ici 3 ==> (row,col)=",(row,col))
                        correlation = loi_de_probalibilite(p_correl, correl_seuil -0.001, nbre_correlations, correl_seuil)
#                        dico_proba_case[(row,col)] = correl_seuil;
                        dico_proba_case[(row,col)] = correlation
                        cases_traites.append( (row,col) );
                    else:
                        # proba entre 0.8001 - 1 ( X -- 1 ) ===> vrai positifs (distributions asymetrique selon fonction loi_proba)
#                        print("ici 4 ==> (row,col)=",(row,col))
                        correlation = loi_de_probalibilite(correl_seuil, 1, nbre_correlations, correl_seuil)
#                        dico_proba_case[(row,col)] = round(random.uniform(correl_seuil+0.001, 1), 3);
                        dico_proba_case[(row,col)] = correlation
                        cases_traites.append( (row,col) );
    return dico_proba_case;   
def ajouter_proba_matE(matE, dico_deleted_add_edges, loi_stats, p_correl, correl_seuil=0.8):
    """
    le but est de definir des probabilites pour chaque case de matE pour simuler 
    la matrice de correlation obtenue apres calcul. les probas (correls) sont les suivantes:
        
        * case 0 -----> 0: proba entre 0 - 0.5 (0, correl_seuil-p_correl-0.01)           ==> vrai negatifs
        * case 1 -----> 0: proba entre 0.6 - 0.79 ( correl_seuil-p_correl, correl_seuil) ==> faux negatifs
        * case 0 -----> 1: proba entre p_correl et correl_seuil-0.001 (p_correl,correl_seuil-0.01)==> faux positifs
        * case 1 -----> 1: proba entre 0.8 - 1 (correl_seuil, 1)             ==> vrai positifs
    NB: X = correl_seuil
    dico_deleted_add_edges = {"ajouter":[(a,b),...],"supprimer":[(c,d),...]}
    loi_stats = loi uniforme, de poisson
    """
    cases_traites = list(); dico_proba_case = dict(); nbre_correlations = 250
    for row in matE.columns.tolist():
        for col in matE.columns.tolist():
            if col != row and (row,col) not in cases_traites and (row,col) not in cases_traites:
                if col != row and matE.loc[row][col] == 0 :
                    if (row,col) in dico_deleted_add_edges["supprimer"] or \
                        (col,row) in dico_deleted_add_edges["supprimer"]:
                        # proba entre 0.6 et 0.79 (entre X-0.2 -- X-0.01) ===> faux negatifs (distributions uniforme)
                        if loi_stats == "poisson":
                            dico_proba_case[(row,col)] = round(random.choice(np.union1d(np.linspace(0.6,0.644,100), \
                                                            np.linspace(0.72,0.79,100))),3)
                        else:
                            correlation = loi_de_probalibilite(correl_seuil-p_correl, p_correl, \
                                                               nbre_correlations, correl_seuil)
#                            dico_proba_case[(row,col)] = round(random.uniform(correl_seuil-0.2, correl_seuil-0.01), 3);
                            dico_proba_case[(row,col)] = correlation;
                        cases_traites.append( (row,col) )
                    else:
                        # proba entre 0 - 0.5 ( 0 -- X-0.3)  ===> vrai negatifs (distributions asymetrique selon fonction loi_proba)
                        if loi_stats == "poisson":
                            dico_proba_case[(row,col)] = round(random.choice(np.union1d(np.linspace(0,0.117,100), \
                                                            np.linspace(0.35,0.5,100))),3)
                        else:
                            correlation = loi_de_probalibilite(0, correl_seuil-p_correl-0.01, \
                                                               nbre_correlations, correl_seuil);
#                            dico_proba_case[(row,col)] = round(random.uniform(0, correl_seuil-0.3), 3);
                            dico_proba_case[(row,col)] =  correlation;              
                        cases_traites.append( (row,col) )
                elif col != row and matE.loc[row][col] == 1 :
                    if (row,col) in dico_deleted_add_edges["ajouter"] or \
                        (col,row) in dico_deleted_add_edges["ajouter"]:
                        # proba entre 0.8 = X ===> faux positifs (distributions uniforme selon fonction loi_proba)
                        correlation = loi_de_probalibilite(p_correl, correl_seuil-0.001,\
                                                           nbre_correlations, correl_seuil)
                        dico_proba_case[(row,col)] = correlation
                        cases_traites.append( (row,col) );
                    else:
                        # proba entre 0.8001 - 1 ( X -- 1 ) ===> vrai positifs (distributions asymetrique selon fonction loi_proba)
                        correlation = loi_de_probalibilite( correl_seuil, 1, \
                                                           nbre_correlations, correl_seuil)
#                        dico_proba_case[(row,col)] = round(random.uniform(correl_seuil+0.001, 1), 3);
                        dico_proba_case[(row,col)] = correlation
                        cases_traites.append( (row,col) );
    return dico_proba_case;               
##### ajout proba de chaque case de matE ====> fin            
#######################  A REVOIR    ##########################################
###############################################################################


def matriceE(dimMat, nbre_lien, chemin_matrices, chemin_datasets, nbre_ts, epsilon, effet_joule, test):
    """
    cette fonction permet de :
        * generer aleatoirement une matrice matA
        * generer les mesures du graphe matA en faisant une propagation de flots 
            descendantes( des sources vers les puits)
        * generer la matrice matE (qui doit etre un linegraph)
        
    IMPORTANT :============> epsilon = 0.75 <===============
    dimMat : dimension de la matrice matA
    seuil : valeur par defaut definissant une adjacence entre 2 sommets
    chemin_datasets = "data/datasets/"
    
    nbre_ts = nombre de time series 
    """
    #generer reseau de flots (A DETERMINER) avec mesures 
    df_matA = None
    if test == "EN TEST":
        df_matA = pd.read_csv("datas/data_test/matrices/df_matA_generer.csv", \
                              index_col = "nodes");
    else:
        df_matA = geneMatA.genererMatriceA(dimMat, nbre_lien)
    df_matA.to_csv(chemin_matrices+"df_matA_generer.csv")
    dico_dual_arc_sommet = mesures.nommage_arcs( df_matA )
    #print("dico_dual_arc_sommet = ", dico_dual_arc_sommet)
    mesures.create_datasets(df_matA, dico_dual_arc_sommet,chemin_datasets, nbre_ts, effet_joule ) 
    
    #matrice du linegraphe du reseau de flot a determiner
    listeArcs = fct_aux.liste_arcs(df_matA)
    matE = mesures.creation_matE(dico_dual_arc_sommet, listeArcs)
    matE.to_csv(chemin_matrices+"matE_test.csv")
    return matE, df_matA, dico_dual_arc_sommet

######### test simulation k-G-alpha conservation parametres debut  ============>  
def aretes_C(C):
    liste_aretes = []
    for Cu in C:
        liste_aretes.extend( list(it.combinations(Cu,2)) )
    return liste_aretes
    
def save_df(df, path_distr_chemin, identifiant, headers):
    """
    merge df to the .csv loaded dataframe 
    """
    if not os.path.exists(path_distr_chemin+"resumeExecution_"+str(identifiant)+".csv"):
        df.to_csv(path_distr_chemin+"resumeExecution_"+str(identifiant)+".csv", sep=",")
        df.colums = headers
        return
        
    df_csv = pd.read_csv(path_distr_chemin+"resumeExecution_"+str(identifiant)+".csv", \
                         sep=",", names = headers)
    pd.concat([df_csv, df])\
        .to_csv(path_distr_chemin+"resumeExecution_"+str(identifiant)+".csv", sep=",")
    
def best_permutation(dico_permutation_cliq, matE, matE_k_alpha):
    """
    selection de la permutation dont le cout et la distance de Hamming sont minimum
    dico_permutation = [ cle : une permutation
                    (A,B,C,D): [C, Ec, dico_cliq, som_cout_min, noeuds_corriges]
    ]
    return dico_sol 
    dico_sol = ["nbre_aretes_matE": , "nbre_aretes_matE_k_alpha": ,\
                "nbre_aretes_LG": , "nbre_aretes_diff_matE_k_alpha_LG": ,\
                "dist_line": , "liste_aretes_diff_matE_k_alpha_LG": , \
                "nbre_aretes_diff_matE_LG": , "hamming": , "liste_aretes_diff_matE_LG": ,\
                "C": , "som_cout_min": , "noeuds_corriges": , "ordre_noeuds_traites":]
    """
    dico_sol = dict(); som_hamming = 0;
    liste_aretes_matE = fct_aux.liste_arcs(matE);
    liste_aretes_matE_k_alpha = fct_aux.liste_arcs(matE_k_alpha);                    
    for tup_node_1, solutions in dico_permutation_cliq.items():
        liste_aretes_LG = None;
        liste_aretes_LG = aretes_C(solutions[0]);
        dist_line = None; liste_aretes_diff_matE_k_alpha_LG = None;
        dist_line, liste_aretes_diff_matE_k_alpha_LG = \
            distance_hamming(liste_aretes_matE_k_alpha, liste_aretes_LG )
        
        hamming = None; liste_aretes_diff_matE_LG = None;
        hamming, liste_aretes_diff_matE_LG = distance_hamming(liste_aretes_matE, liste_aretes_LG)
        
        som_cout = solutions[6];
        som_hamming += hamming;
#        print(" hamming=", hamming," ordre:", tup_node_1," cout:",som_cout," LG:",len(liste_aretes_LG)," MAtE:",len(liste_aretes_matE))
        if (hamming, som_cout) not in dico_sol.keys():
            dico_sol[(hamming, som_cout)] = [{"nbre_aretes_LG": len(liste_aretes_LG ), \
                                "aretes_LG": liste_aretes_LG,\
                                "nbre_aretes_diff_matE_k_alpha_LG": len(liste_aretes_diff_matE_k_alpha_LG),\
                                "dist_line": dist_line, \
                                "liste_aretes_diff_matE_k_alpha_LG": liste_aretes_diff_matE_k_alpha_LG, \
                                "nbre_aretes_diff_matE_LG": len(liste_aretes_diff_matE_LG), \
                                "hamming": hamming, \
                                "liste_aretes_diff_matE_LG": liste_aretes_diff_matE_LG,\
                                "C": solutions[0], "C_old":solutions[7],\
                                "som_cout_min":solutions[6], \
                                "noeuds_corriges": tup_node_1 ,\
                                "ordre_noeuds_traites": solutions[5]
                                }]
        else:
            dico_sol[(hamming, som_cout)].append({"nbre_aretes_LG": len(liste_aretes_LG ), \
                                "aretes_LG": liste_aretes_LG,\
                                "nbre_aretes_diff_matE_k_alpha_LG": len(liste_aretes_diff_matE_k_alpha_LG),\
                                "dist_line": dist_line, \
                                "liste_aretes_diff_matE_k_alpha_LG": liste_aretes_diff_matE_k_alpha_LG, \
                                "nbre_aretes_diff_matE_LG": len(liste_aretes_diff_matE_LG), \
                                "hamming": hamming, \
                                "liste_aretes_diff_matE_LG": liste_aretes_diff_matE_LG,\
                                "C": solutions[0],  "C_old":solutions[7],\
                                "som_cout_min":solutions[6], \
                                "noeuds_corriges": tup_node_1,\
                                "ordre_noeuds_traites": solutions[5] }) 
        
    # selection du min et ajout mean, max
#    print("dico_sol.keys = ", dico_sol.keys())
    print("len(dico_permutation_cliq) = ",len(dico_permutation_cliq))
    min_keys_hamming_cout = min(dico_sol.keys())
    max_keys_hamming_cout = max(dico_sol.keys())
    mean_hamming = som_hamming/len(dico_permutation_cliq)
    ## ecart type --> debut
    ecart_type = 0; som = 0;
    for tup_key, list_dico in dico_sol.items():
        som += len(list_dico) * pow((tup_key[0] - mean_hamming),2)
    ecart_type = pow( som/len(dico_permutation_cliq), 1/2)
    ## ecart type --> fin
    
    print("----> min_hamming = ",min_keys_hamming_cout, " som_ham = ", som_hamming," mean_Ham= ", mean_hamming," len_tupl = ", len(dico_permutation_cliq) )

    
    dico_sol_min = dict();
    dico_sol_min = { \
            "nbre_aretes_LG": dico_sol[min_keys_hamming_cout][0]["nbre_aretes_LG"],\
            "aretes_LG": dico_sol[min_keys_hamming_cout][0]["aretes_LG"],\
            "nbre_aretes_diff_matE_k_alpha_LG": dico_sol[min_keys_hamming_cout][0]["nbre_aretes_diff_matE_k_alpha_LG"],\
            "dist_line": dico_sol[min_keys_hamming_cout][0]["dist_line"],\
            "liste_aretes_diff_matE_k_alpha_LG": dico_sol[min_keys_hamming_cout][0]["liste_aretes_diff_matE_k_alpha_LG"],\
            "nbre_aretes_diff_matE_LG": dico_sol[min_keys_hamming_cout][0]["nbre_aretes_diff_matE_LG"],\
            "hamming": dico_sol[min_keys_hamming_cout][0]["hamming"],\
            "liste_aretes_diff_matE_LG": dico_sol[min_keys_hamming_cout][0]["liste_aretes_diff_matE_LG"],\
            "C": dico_sol[min_keys_hamming_cout][0]["C"], \
            "C_old":dico_sol[min_keys_hamming_cout][0]["C_old"],\
            "som_cout_min": min_keys_hamming_cout[1],\
            "noeuds_corriges": dico_sol[min_keys_hamming_cout][0]["noeuds_corriges"],\
            "ordre_noeuds_traites": dico_sol[min_keys_hamming_cout][0]["ordre_noeuds_traites"],\
            "nbre_aretes_matE": len(liste_aretes_matE),\
            "nbre_aretes_matE_k_alpha": len(liste_aretes_matE_k_alpha),\
            "min_hamming": dico_sol[min_keys_hamming_cout][0]["hamming"],\
            "mean_hamming": mean_hamming,\
            "max_hamming": max_keys_hamming_cout[0],\
            "max_cout": max_keys_hamming_cout[1],\
            "max_permutation": dico_sol[max_keys_hamming_cout][0]["noeuds_corriges"], \
            "ecart_type": ecart_type \
    }

    return dico_sol_min;
    
######################## test simulation k-G-alpha conservation parametres  avec try except DEBUT  ============> 
def simulation_p_correl(nbre_graphe, dimMatA, k_min, k_max, alpha_max, identifiant, arg_params):
    """
    tester BON
    but: tester les valeurs de dist-line entre 2 matrices de correlations
            1ere matrice de correlation ayant k aretes supprimees et 
            2eme matrice de correlation etant sa matrice reconstruite
         puis comparer la matrice de correlation correcte et la matrice reconstruite
         
    nbre_graphe: nbre de graphes generees
    kmax: nombre maximum d'aretes a supprimer k=[1..kmax]
    k_min: nombre minimum d'aretes a supprimer en general k = 0
    alpha_max: nombre maximum de matrices ayant k aretes supprimees
    cpt_graphe_genere: nbre courant de graphe etant genere
    dimMatA: dimension de la matrice A 
    ascendant_1 : definit l'ordre dans lequel les noeuds sont selectionnees.dans la correction de matE
                sil est a True alors cest du plus petit au plus grand
                sil est a True alors cest du plus grand au plus petit
                
    matE_k_aplha: matrice matE supprime de k aretes a alpha-ieme fois
    L_G : la Matrice de correlation corrige a k aretes supprimes et au rang alpha. il est note dans le cahier M'_alpha
    number_permutations_nodes_1: nombre de permutations de sommets -1. 
                                une permutation de noeuds est une selection uniforme de noeuds a -1 dans la liste list_noeuds_1
                                
    arg_params = {"number_permutations_nodes_1": 10(30, 100), "biais": True, "algoGreedy":True, \
                  "mode_select_noeuds_1":"coutMin" or "degreMin", "number_items_pi1_pi2" = 1,\
                  methode_deleted_add_edges = 0, SEUIL_PROBA = 0.8, "proba_seuil": proba_seuil}
    arg_params = {"number_permutations_nodes_1": 10(30, 100), "biais": True, "algoGreedy":True, \
                  "mode_select_noeuds_1":"coutMin" or "degreMin", "number_items_pi1_pi2" = 1,\
                  methode_deleted_add_edges = 0, "p_correl": p_correl, \
                  "correl_seuil":correl_seuil}
    """ 
    nbre_lien = (2,5);
    nbre_lien = arg_params["nbre_lien"]; #5;
    nbre_ts = 10
    effet_joule = 0#0.1
    epsilon = 0.75
    test = "FINI";
    ascendant_1 = True
    simulation = True
    seuil_U = 0;

    headers_df = ["G_cpt", "k", "alpha", "nbre_aretes_matE", "nbre_aretes_matE_k_alpha", "k_modified_edges",\
                  "nbre_aretes_L_G", "nbre_aretes_diff_matE_k_alpha_LG",\
                  "dist_line", "aretes_diff_matE_k_alpha_LG",\
                  "nbre_aretes_diff_matE_LG", "hamming", "aretes_diff_matE_LG",\
                  "C","som_cout_min","noeuds_corriges",\
                  "min_hamming","mean_hamming","max_hamming","ecart_type",\
                  "max_cout","max_permutation",\
                  "dico_som_min_permutations","dico_dual_arc_sommet","ordre_noeuds_traites","C_old"]
        
    path_graphes_chemin = 'representationGraphique/graphes/'
    path_graphes = Path(path_graphes_chemin)
    path_graphes.mkdir(parents=True, exist_ok=True)

    path_distr_chemin = arg_params["rep"]+str(arg_params["coef_fct_cout"][2])+"/"+\
                        str(arg_params["mode_select_noeuds_1"])+"/"+\
                        "data_p_"+str(arg_params["p_correl"])+"/distribution/";
    path_distr = Path(path_distr_chemin)
    path_distr.mkdir(parents=True, exist_ok=True)
    
    df = pd.DataFrame( columns = headers_df)
    cpt_df = 0

    #for k in range(1, k_max+1):
    for k in range(k_min, k_max):
        cpt_graphe_genere = 0;
        while cpt_graphe_genere < nbre_graphe:
            #print("cpt_graphe_genere= ",cpt_graphe_genere, " k=",k);
            cpt_graphe_genere += 1
            G_cpt = "G_"+str(cpt_graphe_genere)+"_"+str(k);
    
            # creer repertoire data contenant matrices et datasets
            # rep =  contient le chemin vers le repertoire correspondant du graphe n genere
            # par exple rep = methode_correction_nodes_1/data_p_XX/G_10 avec 10 = cpt_graphe_genere
            chemin_datasets = arg_params["rep"]+str(arg_params["coef_fct_cout"][2])+"/"+\
                              str(arg_params["mode_select_noeuds_1"])+"/"+\
                              "data_p_"+str(arg_params["p_correl"])+"/"+G_cpt+"/datasets/";
            chemin_matrices = arg_params["rep"]+str(arg_params["coef_fct_cout"][2])+"/"+\
                              str(arg_params["mode_select_noeuds_1"])+"/"+\
                              "data_p_"+str(arg_params["p_correl"])+"/"+G_cpt+"/matrices/";
            path = Path(chemin_datasets);
            path.mkdir(parents=True, exist_ok=True);
            path = Path(chemin_matrices);
            path.mkdir(parents=True, exist_ok=True);
    
            # generer matE
            matE, df_matA, dico_dual_arc_sommet = matriceE(dimMatA, nbre_lien, chemin_matrices,\
                                                           chemin_datasets, nbre_ts, epsilon, \
                                                           effet_joule, test)
            
            aretes_matE = len(fct_aux.liste_arcs(matE))
            moy_distline = 0; moy_hamming = 0; sum_distline = 0; sum_hamming = 0; correl_dl_dh = 0;
            for alpha in range(alpha_max):
                print(G_cpt,"k = ",k," alpha = ",alpha,"DEBUT ")
                # suppression k arteres alpha_max fois
                try :
                    matE_k_alpha = None; k_modified_edges = []; dico_proba_cases = dict();
                    matE_k_alpha, dico_deleted_add_edges = modif_k_cases(matE.copy(), k, \
                                                                         arg_params["methode_delete_add_edges"], 
                                                                         arg_params["p_correl"])
                    k_modified_edges = list(dico_deleted_add_edges.values())
                    dico_proba_cases = ajouter_proba_matE(matE_k_alpha, dico_deleted_add_edges, \
                                                              arg_params["loi_stats"], arg_params["p_correl"],\
                                                              arg_params["correl_seuil"])
                                            
                    matE_k_alpha.to_csv(chemin_matrices+"matE_"+str(k)+"_"+str(alpha)+".csv")
                    # decouverte puis correction matE_k_alpha
                    #{0:[C, dico_cliq, E0, min_cliques_aSupprDe_ens_C, noeuds_corriges, ordre_noeuds_traites, som_cout_min]}
                    dico_permutation_cliq = dict()
                    if arg_params["algoGreedy"] == False:
                        #algo corrigeant tous les noeuds a -1
                        dico_permutation_cliq = \
                            decouvClique.decouverte_cliques(matE_k_alpha, dico_dual_arc_sommet, \
                                                    seuil_U, epsilon, \
                                                    chemin_datasets , chemin_matrices,\
                                                    ascendant_1, simulation,\
                                                    dico_proba_cases,\
                                                    arg_params)
                    else:                        
                        #algo corrigeant un seul noeud a -1 puis modif de matE et enfin algo de couverture + correction
                        dico_permutation_cliq = \
                            decouvClique.decouverte_cliques_corriger_noeuds_new(\
                                                    matE_k_alpha.copy(), dico_dual_arc_sommet, \
                                                    seuil_U, epsilon, \
                                                    chemin_datasets , chemin_matrices,\
                                                    ascendant_1, simulation,\
                                                    dico_proba_cases,\
                                                    arg_params)
                    # Debut selection de la permutation de noeuds dont la distance hamming est la plus petite
                    dico_sol = dict()
                    dico_sol = best_permutation(dico_permutation_cliq, matE, matE_k_alpha)
                    # FIN selection de la permutation de noeuds dont la distance hamming est la plus petite
                    dico_som_min_permutations = dict();
                    for l_noeuds_1, values in dico_permutation_cliq.items():
                        if values[6] not in dico_som_min_permutations.keys():
                            dico_som_min_permutations[values[6]] = [l_noeuds_1]
                        else:
                            dico_som_min_permutations[values[6]].append(l_noeuds_1)
                
                    cpt_df += 1;
                    df.loc[cpt_graphe_genere] = [G_cpt, k, alpha, \
                                    dico_sol["nbre_aretes_matE"], dico_sol["nbre_aretes_matE_k_alpha"], \
                                    k_modified_edges,\
                                    dico_sol["nbre_aretes_LG"], dico_sol["nbre_aretes_diff_matE_k_alpha_LG"],\
                                    dico_sol["dist_line"], dico_sol["liste_aretes_diff_matE_k_alpha_LG"], \
                                    dico_sol["nbre_aretes_diff_matE_LG"], dico_sol["hamming"],\
                                    dico_sol["liste_aretes_diff_matE_LG"],\
                                    dico_sol["C"], dico_sol["som_cout_min"], \
                                    dico_sol["noeuds_corriges"], \
                                    dico_sol["min_hamming"],dico_sol["mean_hamming"], \
                                    dico_sol["max_hamming"],dico_sol["ecart_type"],\
                                    dico_sol["max_cout"], dico_sol["max_permutation"],\
                                    dico_som_min_permutations, dico_dual_arc_sommet,\
                                    dico_sol["ordre_noeuds_traites"], dico_sol["C_old"]]
                    if cpt_df % 100 == 0:
                        save_df(df, path_distr_chemin, identifiant, headers_df)
                        df = pd.DataFrame( columns = headers_df)
                        print("save %s fois" %cpt_df)
                    # sommer distline et hamming
                    # je divise l'erreur trouve(DL,DH) 
                    # par le nbre d'aretes de matE_k_alpha(DL) ou
                    # par le nbre d'aretes de matE
                    sum_distline += dico_sol["dist_line"] # dico_sol["dist_line"]/dico_sol["nbre_aretes_matE_k_alpha"] 
                    sum_hamming += dico_sol["hamming"] # dico_sol["hamming"]/dico_sol["nbre_aretes_matE"]
                    print(G_cpt,"k = ",k," alpha = ",alpha," ===> FIN termine")
                    print("{}, k={}, {}, p_correl={}, {}, s={}".format(G_cpt, k, arg_params["coef_fct_cout"][2], \
                          arg_params["p_correl"], arg_params["mode_select_noeuds_1"],\
                          arg_params["correl_seuil"] ))
                except Exception as e:
                    print("####### EmptyDataError {} {}, s={}, p={}, k={}, {}: e = {} #######".format(
                          arg_params["coef_fct_cout"][2],arg_params["mode_select_noeuds_1"], \
                          arg_params["correl_seuil"], arg_params["p_correl"], k, G_cpt, e));
                    cpt_df += 1;
                    df.loc[cpt_graphe_genere] = [G_cpt, k, alpha, \
                                    "error", "error", \
                                    "error",\
                                    "error","error" ,\
                                    "error","error" , \
                                    "error","error" ,\
                                    "error",\
                                    "error", "error", \
                                    "error", \
                                    "error","error", "error",\
                                    "error",\
                                    "error","error" ,\
                                    "error","error","error", "error"]
            # moyenner dist_line et hamming pour k aretes supprimes
            moy_distline = sum_distline / alpha_max;
            moy_hamming = sum_hamming / alpha_max;
            if moy_hamming == 0 and moy_distline == 0:
                correl_dl_dh = 1
            else:
                correl_dl_dh = abs(moy_hamming - moy_distline)/max(moy_hamming, moy_distline)
            
            # ecrire dans un fichier pouvant etre lu pendant qu'il continue d'etre ecrit
            f = open(path_distr_chemin+"distribution_moyDistLine_moyHamming_k_"+str(k)+".txt","a")
            f.write(G_cpt+";"+str(k)+";"+str(moy_distline)+";"+str(moy_hamming)+";"+str(aretes_matE)+\
                    ";"+str(correl_dl_dh)+"\n")
            f.close();
    print("k_min = ",k_min," termine")
    save_df(df, path_distr_chemin, identifiant, headers_df)
######################## test simulation k-G-alpha conservation parametres  avec try except FIN  ============> 

if __name__ == '__main__':
    
     start= time.time();
     
     ### test parallelisation\
     rep = "data/";
     nbre_lien = 5;
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
     critere_selection_pi1_pi2 = 2; # 0: moins de modif,1: ajout aretes> supp aretes, 2:ajout aretes < supp aretes,
     number_permutations_nodes_1= 10 #100;
     facteur_multiplicatif = 1; exposant = 1; # 0: fct_unitaire, 1:fct_normal, 2: fct_quadratique, 4:fct_quadruple, 5: fct_quintuple
     type_fct_cout = "cloche" # ou "lineaire"
     coef_fct_cout = (exposant, facteur_multiplicatif, type_fct_cout)
     mode_select_noeuds_1 = "aleatoire" #"degreMin" #"coutMin" # "degreMin" # aleatoire
     number_items_pi1_pi2 = 1;
     arg_params = {"number_items_pi1_pi2": number_items_pi1_pi2, "nbre_lien":nbre_lien,\
                   "number_permutations_nodes_1": number_permutations_nodes_1, \
                   "methode_delete_add_edges": methode_delete_add_edges, \
                   "loi_stats": loi_stats,"biais": biais, \
                   "p_correl": p_correl, "correl_seuil": correl_seuil,\
                   "algoGreedy":algoGreedy, \
                   "mode_select_noeuds_1":mode_select_noeuds_1,\
                   "coef_fct_cout":coef_fct_cout,\
                   "critere_selection_pi1_pi2":critere_selection_pi1_pi2, "rep":rep};
#     simulation_test(nbre_graphe, dimMatA, k_min, k_min+1, alpha_max,k_min, arg_params)

     #### parallele
     k=10; k = 10; 
     dimMatA = 6#10#5 #10 #15 
     nbre_graphe = 150#30#10#300#5#200
     alpha_max = 1;
     k_min_range = list(range(0,k)) # k=6;list(range(0, k*10, 10)); #[10,20,30,40,50]
     k_max_range = list(range(1,k+1)) # k=6;list(range(1, k*10, 10)); #[11,21,31,41,50]
     
     dimMatA = 8 #5 #10 #15 
     nbre_graphe = 100
     k_min_range = [15,25] #[15, 25, 30] 
     k_max_range = [16,26]#[16, 26, 31]
     N = len(k_max_range);
    
     ##### p_correls seuils fct_cout ===> debut
     critere_selection_pi1_pi2 = 0; 
     arg_params["critere_selection_pi1_pi2s"] = critere_selection_pi1_pi2;
     type_fct_couts = ["lineaire","cloche"]; 
     type_fct_couts = ["lineaire_simul50Graphes_priorite_supp", "lineaire_simul50Graphes_priorite_ajout", "lineaire_simul50Graphes_priorite_aucune"];
     type_fct_couts = ["lineaire_simul50Graphes_priorite_aucune"]
     p_correls = [0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9]; p_correls = [0.5];
     correl_seuils = [0.1, 0.2, 0.3, 0.4, 0.5 ,0.6, 0.7, 0.8, 0.9]; correl_seuils = [0.6];
     mode_select_noeuds_1s = ["degreMin","coutMin","aleatoire"]; mode_select_noeuds_1s = ["aleatoire"]; 
     facteur_multiplicatif = 1; exposant = 0;
     params = list();
     for type_fct_cout in type_fct_couts:
         arg_params["coef_fct_cout"] = (exposant,facteur_multiplicatif,type_fct_cout);
         for p_correl in p_correls:
             arg_params["p_correl"] = p_correl;
             for mode_select_noeuds_1 in mode_select_noeuds_1s:
                 arg_params["mode_select_noeuds_1"] = mode_select_noeuds_1;
                 arg_params_ = arg_params.copy();
                 
                 params_ = list(zip( [nbre_graphe]*N, [dimMatA]*N, k_min_range, k_max_range, \
                                    [alpha_max]*N, k_min_range, [arg_params_]*N))
                 params.append(params_)
                 
     params = [par for tup in params for par in tup]
     print("params = {}".format(len(params)))
     ##### p_correls seuils fct_cout ===> fin 
     print("params ",len(params))
     p = Pool(mp.cpu_count()-1) 
     p.starmap(simulation_p_correl, params)
     p.terminate()
     #### parallele
     k_min = "0_9";
     g = open("tempsExecution_"+str(k_min)+".txt","w")
     ti = time.time() - start
     g.write(str(ti))
     g.close()
     
     print (time.time() - start)
     
     #print("FIN")

     
