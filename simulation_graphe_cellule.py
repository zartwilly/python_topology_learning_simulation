#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 19 00:39:00 2018

@author: willy
"""
import time;
import re, os, random, math;
import itertools as it;
import pandas as pd;
import numpy as np;
import fonctions_auxiliaires as fct_aux
import decouverte_cliques as decouvClique
import multiprocessing as mp;
from pathlib import Path    # http://stackoverflow.com/questions/6004073/how-can-i-create-directories-recursively
from multiprocessing import Pool;
import simulation_50_graphes_PARALLELE as simu50;
from collections import Counter

import matplotlib.pyplot as plt
import networkx as nx;

#------------------------------------------------------------------------------
#                   fonctions communes aux codes
#------------------------------------------------------------------------------
def isOdd( number ):
	return number % 2 != 0
 
def isEven( number ):
	return number % 2 == 0
	
def countOdd( numbers ):
	return len( 
		list(filter( lambda x: isOdd( x ), numbers ))
	)
def countEven( numbers ):
	return len( 
		list(filter( lambda x: isEven( x ), numbers ))
	)

def aretes_C(C):
    aretes = list()
    for clique in C:
        if len(clique)>1:
            aretes.extend(list(it.combinations(clique,2)))
    return aretes;
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

def gamma_noeud(liste_aretes):
    """
    """
    liste_noeuds = [ [arete[0], arete[1]] for arete in liste_aretes]
    set_noeuds = set([item for sublist in liste_noeuds for item in sublist])
    dico = dict()
    for noeud in set_noeuds:
        ens = set()
        for arc in liste_aretes:
            if noeud == arc[0]:
                ens.add( arc[1] )
            if noeud == arc[1]:
                ens.add( arc[0] )
        dico[noeud] = [len(ens), ens]
    return dico   
    
# ----- best permutation iourte ------
def best_permutation_cellule(dico_permutation_cliq, matE_G_k):
    """
    selection de la permutation dont le cout et la distance de Hamming sont minimum
    dico_permutation = [ cle : une permutation
                    (A,B,C,D): [C, Ec, dico_cliq, som_cout_min, noeuds_corriges]
    ]
    return dico_sol 
    dico_sol = {"nbre_aretes_LG":,"aretes_LG":,"nbre_aretes_diff_matE_G_k_LG": ,\
                "dist_line":,"aretes_diff_matE_G_k_LG":,"C":,"som_cout_min":,\
                "noeuds_corriges":,"ordre_noeuds_traites":}
    """
    dico_sol = dict(); som_dist_line = 0
    aretes_matE_G_k = fct_aux.liste_arcs(matE_G_k);                    
    for tup_node_1, solutions in dico_permutation_cliq.items():
        aretes_LG = None;
        aretes_LG = simu50.aretes_C(solutions[0]);
        
        dist_line = None; aretes_diff_matE_G_k_LG = None;
        dist_line, aretes_diff_matE_G_k_LG = \
            simu50.distance_hamming( aretes_matE_G_k, aretes_LG )
        
                
        som_cout = solutions[6]; som_dist_line += dist_line;
#        print(" hamming=", hamming," ordre:", tup_node_1," cout:",som_cout," LG:",len(liste_aretes_LG)," MAtE:",len(liste_aretes_matE))
        if (dist_line, som_cout) not in dico_sol.keys():
            dico_sol[(dist_line, som_cout)] = [{"nbre_aretes_LG": len(aretes_LG ), \
                                "aretes_LG": aretes_LG,\
                                "nbre_aretes_diff_matE_G_k_LG": len(aretes_diff_matE_G_k_LG),\
                                "dist_line": dist_line, \
                                "aretes_diff_matE_G_k_LG": aretes_diff_matE_G_k_LG, \
                                "C": solutions[0], \
                                "som_cout_min":solutions[6], \
                                "noeuds_corriges": tup_node_1 ,\
                                "ordre_noeuds_traites": solutions[5]
                                }]
        else:
            dico_sol[(dist_line, som_cout)].append({"nbre_aretes_LG": len(aretes_LG ), \
                                "aretes_LG": aretes_LG,\
                                "nbre_aretes_diff_matE_G_k_LG": len(aretes_diff_matE_G_k_LG),\
                                "dist_line": dist_line, \
                                "aretes_diff_matE_G_k_alpha_LG": aretes_diff_matE_G_k_LG, \
                                "C": solutions[0], \
                                "som_cout_min":solutions[6], \
                                "noeuds_corriges": tup_node_1,\
                                "ordre_noeuds_traites": solutions[5] }) 
        
    # selection du min et ajout mean, max
#    print("dico_sol.keys = ", dico_sol.keys())
    print("len(dico_permutation_cliq) = ",len(dico_permutation_cliq))
    min_keys_line_cout = min(dico_sol.keys())
    max_keys_line_cout = max(dico_sol.keys())
    mean_line = som_dist_line/len(dico_permutation_cliq)
    ## ecart type --> debut
    ecart_type = 0; som = 0;
    for tup_key, list_dico in dico_sol.items():
        som += len(list_dico) * pow((tup_key[0] - mean_line),2)
    ecart_type = pow( som/len(dico_permutation_cliq), 1/2)
    ## ecart type --> fin
    
    print("----> min_line = ",min_keys_line_cout, " som_line = ", som_dist_line,\
          " mean_line= ", mean_line," len_tupl = ", len(dico_permutation_cliq) )

    dico_sol_min = dict();
    dico_sol_min = { \
            "nbre_aretes_LG": dico_sol[min_keys_line_cout][0]["nbre_aretes_LG"],\
            "aretes_LG": dico_sol[min_keys_line_cout][0]["aretes_LG"],\
            "nbre_aretes_diff_matE_G_k_LG": dico_sol[min_keys_line_cout][0]["nbre_aretes_diff_matE_G_k_LG"],\
            "dist_line": dico_sol[min_keys_line_cout][0]["dist_line"],\
            "aretes_diff_matE_G_k_LG": dico_sol[min_keys_line_cout][0]["aretes_diff_matE_G_k_LG"],\
            "C": dico_sol[min_keys_line_cout][0]["C"], \
            "som_cout_min": min_keys_line_cout[1],\
            "noeuds_corriges": dico_sol[min_keys_line_cout][0]["noeuds_corriges"],\
            "ordre_noeuds_traites": dico_sol[min_keys_line_cout][0]["ordre_noeuds_traites"],\
            "min_line": dico_sol[min_keys_line_cout][0]["dist_line"],\
            "mean_line": mean_line,\
            "max_line": max_keys_line_cout[0],\
            "max_cout": max_keys_line_cout[1],\
            "max_permutation": dico_sol[max_keys_line_cout][0]["noeuds_corriges"], \
            "ecart_type": ecart_type \
    }

    return dico_sol_min;
# ----- best permutation iourte ------

# ----- comparaison_aretes_G_k_LG ------
def comparaison_aretes_G_k_LG(aretes_matE_G_k, aretes_LG):
    """ recherche les aretes de G_k n'appartenant pas aux aretes de aretes_LG """
    cpt_aretes_G_k_notIn_LG= 0;
    for arete_G_k in aretes_matE_G_k:
        if (arete_G_k[0],arete_G_k[1]) not in aretes_LG and \
            (arete_G_k[1],arete_G_k[0]) not in aretes_LG :
            cpt_aretes_G_k_notIn_LG += 1;
    return math.ceil( (cpt_aretes_G_k_notIn_LG/len(aretes_matE_G_k))*100);
# ----- comparaison_aretes_G_k_LG ------
    
#------------------------------------------------------------------------------
#                   graphe cellule G_nn 
#                       nn = n*n est le nombre de cellules 
#------------------------------------------------------------------------------
def matrice_G_nn(N,M):
    """ 
    creer le dataframe M_G_nn avec pour colonnes et index les sommets de G_nn
    """
    # liste de sommets de G_nn
    sommets = list(it.product(range(N), range(M), repeat = 1));
    # creation dataframe avec colonnes = liste de sommets G_nn
    M_G_nn = pd.DataFrame(0, columns = sommets, index = sommets);
    return M_G_nn;
    
def G_nn_old(N, M):
    """ definir une matrice de N*M cases 
        idealement N = M
        row = {0,...,N-1}
        col = {0,...,M-1}
    trop lent pour des N,M > 20
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
                
    M_G_nn.at[(0,0),(0,M-1)] = 1 ; M_G_nn.at[(0,M-1),(0,0)] = 1
    M_G_nn.at[(N-1,0),(0,0)] = 1 ; M_G_nn.at[(0,0),(N-1,0)] = 1; 
    M_G_nn.at[(N-1,0),(N-1,M-1)] = 1; M_G_nn.at[(N-1,M-1),(N-1,0)] = 1;
    M_G_nn.at[(N-1,0),(0,0)] = 1 ; M_G_nn.at[(0,0),(N-1,0)] = 1 
    M_G_nn.at[(0,M-1),(N-1,M-1)] = 1 ; M_G_nn.at[(N-1,M-1),(0,M-1)] = 1 
    M_G_nn.at[(0,M-1),(0,0)] = 1 ; M_G_nn.at[(0,0),(0,M-1)] = 1  
    M_G_nn.at[(N-1,M-1),(0,M-1)] = 1 ; M_G_nn.at[(0,M-1),(N-1,M-1)] = 1;  
    M_G_nn.at[(N-1,M-1),(N-1,0)] = 1 ; M_G_nn.at[(N-1,0),(N-1,M-1)] = 1;
    # renommer index et colonnes
    M_G_nn.rename(columns = dico_sommets, inplace = True);
    M_G_nn.rename(index = dico_sommets, inplace = True);
    aretes_G_nn = list();
    aretes_G_nn, dico_proba_cases = aretes(M_G_nn, k0_1=1);                                        #aretes_G_nn= aretes obtenus a partir du renommage des sommets (row, col)
    return M_G_nn, dico_sommets, aretes_G_nn, dico_proba_cases;
    
def G_nn_debug(N, M):
    """ trop lent pour des N,M > 20"""
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
            aretes_G_nn_row_col.append((row,col))                               # TODO a revoir le (row,col)
            
    M_G_nn.at[(0,0),(0,M-1)] = 1 ; M_G_nn.at[(0,M-1),(0,0)] = 1
    M_G_nn.at[(N-1,0),(0,0)] = 1 ; M_G_nn.at[(0,0),(N-1,0)] = 1; 
    M_G_nn.at[(N-1,0),(N-1,M-1)] = 1; M_G_nn.at[(N-1,M-1),(N-1,0)] = 1;
    M_G_nn.at[(N-1,0),(0,0)] = 1 ; M_G_nn.at[(0,0),(N-1,0)] = 1 
    M_G_nn.at[(0,M-1),(N-1,M-1)] = 1 ; M_G_nn.at[(N-1,M-1),(0,M-1)] = 1 
    M_G_nn.at[(0,M-1),(0,0)] = 1 ; M_G_nn.at[(0,0),(0,M-1)] = 1  
    M_G_nn.at[(N-1,M-1),(0,M-1)] = 1 ; M_G_nn.at[(0,M-1),(N-1,M-1)] = 1;  
    M_G_nn.at[(N-1,M-1),(N-1,0)] = 1 ; M_G_nn.at[(N-1,0),(N-1,M-1)] = 1;
    
    # renommer index et colonnes
    M_G_nn.rename(columns = dico_sommets, inplace = True);
    M_G_nn.rename(index = dico_sommets, inplace = True);
    aretes_G_nn = list();
    aretes_G_nn, dico_proba_cases = aretes(M_G_nn, k0_1=1);                                        #aretes_G_nn= aretes obtenus a partir du renommage des sommets (row, col)
    return M_G_nn, dico_sommets, aretes_G_nn, dico_proba_cases;
    
def G_nn_old_probleme_avec_cpt_aretes_Gnn_k(N, M):
    """ definir une matrice de N*M cases 
        idealement N = M
        row = {0,...,N-1}
        col = {0,...,M-1}
    trop lent pour des N,M > 20
    """ 
    sommets = list(it.product(range(N), range(M), repeat = 1));
    # creation dataframe avec colonnes = liste de sommets G_nn
    M_G_nn = pd.DataFrame(0, columns = sommets, index = sommets);
    
    dico_sommets = dict(); cpt = 0; dico_proba_cases = dict();
    aretes_G_nn_row_col = set(); aretes_G_nn = set(); dico_proba_cases = dict();
    
    for sommet in sommets:
        row = sommet[0]; col = sommet[1];
        cpt += 1; dico_sommets[(row,col)] = str(cpt); 
        dico_proba_cases[(row, col)] = 0.5; 
        if col+1 < M:
            M_G_nn.at[(row,col),(row,col+1)] = 1;
            M_G_nn.at[(row,col+1),(row,col)] = 1;
            aretes_G_nn_row_col.update([(row,col),(row,col+1)])
            if (row,col+1) not in dico_sommets.keys():
                cpt += 1; dico_sommets[(row,col+1)] = str(cpt);
            aretes_G_nn.add((dico_sommets[(row,col)],dico_sommets[(row,col+1)]))
        if col-1 >= 0:
            M_G_nn.at[(row,col),(row,col-1)] = 1;
            M_G_nn.at[(row,col-1),(row,col)] = 1;
            aretes_G_nn_row_col.update([(row,col),(row,col-1)])
            if (row,col-1) not in dico_sommets.keys():
                cpt += 1; dico_sommets[(row,col-1)] = str(cpt);
            aretes_G_nn.add((dico_sommets[(row,col)],dico_sommets[(row,col-1)]))
        if row+1 < N:
            M_G_nn.at[(row,col),(row+1,col)] = 1;
            M_G_nn.at[(row+1,col),(row,col)] = 1;
            aretes_G_nn_row_col.update([(row,col),(row+1,col)])
            if (row+1,col) not in dico_sommets.keys():
                cpt += 1; dico_sommets[(row+1,col)] = str(cpt);
            aretes_G_nn.add((dico_sommets[(row,col)],dico_sommets[(row+1,col)]))
        if row-1 >= 0:
            M_G_nn.at[(row,col),(row-1,col)] = 1;
            M_G_nn.at[(row-1,col),(row,col)] = 1;
            aretes_G_nn_row_col.update([(row,col),(row-1,col)])
            if (row-1,col) not in dico_sommets.keys():
                cpt += 1; dico_sommets[(row-1,col)] = str(cpt);
            aretes_G_nn.add((dico_sommets[(row,col)],dico_sommets[(row-1,col)]))
    
    M_G_nn.at[(0,0),(0,M-1)] = 1; M_G_nn.at[(0,M-1),(0,0)] = 1; 
    M_G_nn.at[(N-1,0),(0,0)] = 1; M_G_nn.at[(0,0),(N-1,0)] = 1;
    aretes_G_nn_row_col.update([(0,0),(0,M-1),(N-1,0)])
    if (0,M-1) not in dico_sommets.keys() :
        cpt += 1; dico_sommets[(0,M-1)] = str(cpt);
    if (N-1,0) not in dico_sommets.keys() :
        cpt += 1; dico_sommets[(N-1,0)] = str(cpt);
    aretes_G_nn.update([(dico_sommets[(0,0)],dico_sommets[(0,M-1)]),\
                        (dico_sommets[(0,0)],dico_sommets[(N-1,0)])])
    
    M_G_nn.at[(N-1,0),(N-1,M-1)] = 1; M_G_nn.at[(N-1,M-1),(N-1,0)] = 1; 
    M_G_nn.at[(N-1,0),(0,0)] = 1; M_G_nn.at[(0,0),(N-1,0)] = 1;
    aretes_G_nn_row_col.update([(N-1,0),(N-1,M-1),(0,0)])
    if (N-1,M-1) not in dico_sommets.keys() :
        cpt += 1; dico_sommets[(N-1,M-1)] = str(cpt);
    if (N-1,0) not in dico_sommets.keys() :
        cpt += 1; dico_sommets[(N-1,0)] = str(cpt);
    aretes_G_nn.update([(dico_sommets[(N-1,0)],dico_sommets[(N-1,M-1)]),\
                        (dico_sommets[(N-1,0)],dico_sommets[(0,0)])])

    M_G_nn.at[(0,M-1),(N-1,M-1)] = 1; M_G_nn.at[(N-1,M-1),(0,M-1)] = 1;
    M_G_nn.at[(0,M-1),(0,0)] = 1; M_G_nn.at[(0,0),(0,M-1)] = 1;
    aretes_G_nn_row_col.update([(0,M-1),(N-1,M-1),(N-1,0)])
    if (N-1,M-1) not in dico_sommets.keys() :
        cpt += 1; dico_sommets[(N-1,M-1)] = str(cpt);
    if (0,M-1) not in dico_sommets.keys() :
        cpt += 1; dico_sommets[(0,M-1)] = str(cpt);
    aretes_G_nn.update([(dico_sommets[(0,M-1)],dico_sommets[(N-1,M-1)]),\
                        (dico_sommets[(0,M-1)],dico_sommets[(0,0)])])
    
    M_G_nn.at[(N-1,M-1),(0,M-1)] = 1; M_G_nn.at[(0,M-1),(N-1,M-1)] = 1; 
    M_G_nn.at[(N-1,M-1),(N-1,0)] = 1; M_G_nn.at[(N-1,0),(N-1,M-1)] = 1; 
    aretes_G_nn_row_col.update([(N-1,M-1),(0,M-1),(N-1,0)])
    if (N-1,M-1) not in dico_sommets.keys() :
        cpt += 1; dico_sommets[(N-1,M-1)] = str(cpt);
    if (0,M-1) not in dico_sommets.keys() :
        cpt += 1; dico_sommets[(0,M-1)] = str(cpt);
    if (N-1,0) not in dico_sommets.keys() :
        cpt += 1; dico_sommets[(N-1,0)] = str(cpt);
    aretes_G_nn.update([(dico_sommets[(N-1,M-1)],dico_sommets[(0,M-1)]),\
                        (dico_sommets[(N-1,M-1)],dico_sommets[(N-1,0)])])
    
#    G = nx.Graph(M_G_nn.values);                                              # plot graph with networkX
#    nx.draw(G, pos=nx.spring_layout(G), with_labels=True);
#    print("aretes_G_nn={},M_G_nn={},dico_sommets={} ".format(aretes_G_nn, M_G_nn.shape,dico_sommets ))
    
    # aretes dans M_G_nn
#    aretes_G_nn = 
    
    # renommer index et colonnes
    M_G_nn.rename(columns = dico_sommets, inplace = True);
    M_G_nn.rename(index = dico_sommets, inplace = True);
    print("M_G_nn={}".format(M_G_nn.columns))

    return M_G_nn, dico_sommets, aretes_G_nn, dico_proba_cases;

def G_nn_commentaire(N, M):
    sommets = list(it.product(range(N), range(M), repeat = 1));
    # creation dataframe avec colonnes = liste de sommets G_nn
    M_G_nn = pd.DataFrame(0, columns = sommets, index = sommets);
    
    dico_sommets = dict(); cpt = 0; 
    aretes_G_nn_row_col = set(); aretes_G_nn = set(); dico_proba_cases = dict();
    
    for sommet in sommets:
        row = sommet[0]; col = sommet[1];
        cpt += 1; dico_sommets[(row,col)] = str(cpt); 
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
    
    # networkx
#    G = nx.Graph(M_G_nn.values);                                              # plot graph with networkX
#    nx.draw(G, pos=nx.spring_layout(G), with_labels=True);
#    print("aretes_G_nn={},M_G_nn={},dico_sommets={} ".format(aretes_G_nn, M_G_nn.shape,dico_sommets ))

        
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

def G_nn(N, M):
    """
    definir une matrice de N*M cases 
        idealement N = M = nombre pair
        row = {0,...,N-1}
        col = {0,...,M-1}
    """
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
        
#    print("aretes_G_nn")
#    print("M_G_nn=\n{}\n, dico_sommets={}\n, aretes_G_nn={}\n, \
#          dico_proba_cases={}"\
#          .format(M_G_nn.shape, len(dico_sommets), \
#                  len(aretes_G_nn), len(dico_proba_cases)));
    # networkx
#    G = nx.Graph(M_G_nn.values);                                              # plot graph with networkX
#    nx.draw(G, pos=nx.spring_layout(G), with_labels=True);
#    print("aretes_G_nn={},M_G_nn={},dico_sommets={} ".format(aretes_G_nn, M_G_nn.shape,dico_sommets ))
    return M_G_nn, dico_sommets, aretes_G_nn, dico_proba_cases;
    pass 
#------------------------------------------------------------------------------
#                   generation k graphes cellule G_nn 
#                       k est le nombre de profondeur de G_nn
#                       k = {0,1,2,...,50} 
#------------------------------------------------------------------------------
def generer_k_graphes(k_min,k_max, rep):
    graphes = []; 
    rep = rep+"graphes_k/"
    for k in range(k_min,k_max+1):                                              # k_min = 1
        N = 2*k;
        start = time.time()
        M_G_nn, dico_sommets_row_col, aretes_G_nn, dico_proba_cases = G_nn(N, N);
        
        path_distr = Path(rep); path_distr.mkdir(parents=True, exist_ok=True);
#        M_G_nn.to_csv(rep+"graphe_"+str(k)+".csv", sep = ";");
        dico_graphes = dict();
        dico_graphes = {"k":k,\
                        "M_G_nn":M_G_nn, \
                        "aretes_G_nn_k":list(aretes_G_nn), \
                        "dico_proba_cases":dico_proba_cases,\
                        "dico_sommets_row_col":dico_sommets_row_col}
        graphes.append(dico_graphes);
        print("k={} runtime={} termine".format(k, time.time() - start))
    return graphes;
    
#------------------------------------------------------------------------------
#                   simulation d'un graphe cellule
#------------------------------------------------------------------------------
def simulation_graphe_cellule(args):
    headers_df = ["G_cpt", "nbre_aretes_G_k", "nbre_aretes_LG", \
                  "dist_line", "nbre_aretes_diff_G_k_LG", \
                  "aretes_diff_G_k_LG", "orientation_aretes_corrigees","count_cliques",\
                  "C", "len(C)", "som_cout_min",\
                  "noeuds_corriges", "ordre_noeuds_traites",\
                  "min_DL", "mean_DL", "max_DL", "ecart_type",\
                  "max_cout", "max_permutation",\
                  "dico_som_min_permutations"];
            
    df_debug = pd.DataFrame( columns = headers_df);
    G_cpt = "G_"+str(args["k"]);
    
    # creation repertoire contenant distribution
    path_distr = Path(args["path_save"]+args["mode_select_noeuds_1"]+"_cellule/");
    path_distr.mkdir(parents=True, exist_ok=True);
    path_save = args["path_save"]+args["mode_select_noeuds_1"]+"_cellule/";
    
    try:
        ordre_noeuds_traites = [] # car liste_cliques = []
        cliques = []; # car graphe iourte;
        dico_cliq = dict();
        for noeud in args["M_G_nn"].columns:
            dico_cliq[noeud] = -1
        dico_gamma_noeud = fct_aux.gamma_noeud(\
                                    args["M_G_nn"], args["aretes_G_nn_k"]) 
        
        # algo de correction selon methodes (arg_params["mode_select_noeuds_1"])
        dico_permutations = dict();
        dico_permutations = decouvClique.solution_methode_nodes_1(dico_gamma_noeud,\
                             cliques, args["aretes_G_nn_k"], ordre_noeuds_traites, \
                             dico_cliq, args["dico_proba_cases"], args);
        
        # Debut selection de la permutation de noeuds dont la distance hamming est la plus petite
        dico_sol = dict()
        dico_sol = best_permutation_cellule(dico_permutations, args["M_G_nn"])
        # FIN selection de la permutation de noeuds dont la distance hamming est la plus petite
        
        
        # comparaison entre aretes_matE_G_k et aretes_LG
        cpt_aretes_G_k_notIn_LG = 0; #(en pourcentage)
        cpt_aretes_G_k_notIn_LG = comparaison_aretes_G_k_LG(\
                                 args["aretes_G_nn_k"], dico_sol["aretes_LG"]);

        
        # ecrire dans un fichier pouvant etre lu pendant qu'il continue d'etre ecrit
        f = open(path_save+"distribution_moyDistLine_G_k.txt","a")
        f.write(G_cpt+";"+str(args["k"])+";"+str(dico_sol["dist_line"])+";"\
                +str(len(args["aretes_G_nn_k"]))+";"+str(cpt_aretes_G_k_notIn_LG)+"\n")
        f.close();
        
        #----------- pour debug, log, ..... -------------
        dico_som_min_permutations = dict();
        for l_noeuds_1, values in dico_permutations.items():
            if values[6] not in dico_som_min_permutations.keys():
                dico_som_min_permutations[values[6]] = [l_noeuds_1]
            else:
                dico_som_min_permutations[values[6]].append(l_noeuds_1)
        
        # combien de cliques a 4, 3 2 et 1
        count_cliques = Counter([len(C) for C in dico_sol["C"]])
        
        # comment sont ajoute/supprimer les aretes (verticalement/diagonalement/horizontalement)
        notes_aretes_diff = []; orientation_aretes_corrigees = ""
        for arete_diff in list(dico_sol["aretes_diff_matE_G_k_LG"])[:150]:
            sommet_gauche = arete_diff[0]; sommet_droit = arete_diff[1]
            print("sommet_gauche={}, sommet_droit={}".format(sommet_gauche, sommet_droit))
            if arete_diff[0] in args["dico_sommets_row_col"].keys():
                sommet_gauche = int(args["dico_sommets_row_col"][arete_diff[0]].split("_")[0])
            if arete_diff[1] in args["dico_sommets_row_col"].keys():
                sommet_droit = int(args["dico_sommets_row_col"][arete_diff[1]].split("_")[1])
            som = 0;
            som = int(sommet_gauche.split("_")[0])+int(sommet_gauche.split("_")[1])+\
                  int(sommet_droit.split("_")[0])+int(sommet_droit.split("_")[1]);
            notes_aretes_diff.append(som);
        if countEven(notes_aretes_diff) > countOdd(notes_aretes_diff):
            orientation_aretes_corrigees = "PAIR"
        elif countEven(notes_aretes_diff) < countOdd(notes_aretes_diff):
            orientation_aretes_corrigees = "IMPAIR"
        else:
            orientation_aretes_corrigees = "NONE";
            
        df_debug.loc[len(df_debug.index)] = [\
                        G_cpt, len(args["aretes_G_nn_k"]),\
                        dico_sol["nbre_aretes_LG"], dico_sol["dist_line"],\
                        dico_sol["nbre_aretes_diff_matE_G_k_LG"],\
                        dico_sol["aretes_diff_matE_G_k_LG"],\
                        orientation_aretes_corrigees,\
                        count_cliques,\
                        dico_sol["C"], len(dico_sol["C"]), dico_sol["som_cout_min"],\
                        dico_sol["noeuds_corriges"], dico_sol["ordre_noeuds_traites"],\
                        dico_sol["min_line"], dico_sol["mean_line"],\
                        dico_sol["max_line"], dico_sol["ecart_type"],\
                        dico_sol["max_cout"], dico_sol["max_permutation"],\
                        dico_som_min_permutations]
        if args["k"] % 100 == 0:
            simu50.save_df(df_debug, path_save, args["k"], headers_df)
            df_debug = pd.DataFrame( columns = headers_df)
            print("save {} fois".format( args["k"] )) 
        if args["debug"]:
            edges_C = aretes_C(dico_sol["C"]); #mylabels = get_labels()
            print("C={}, \n edges_C={}, \n aretes_G_nn_k={}, \n dico_sommets_row_col={}"\
                  .format(dico_sol["C"], edges_C, args["aretes_G_nn_k"], \
                          args["dico_sommets_row_col"]))
            G = nx.Graph();
#            G.add_edges_from(edges_C); 
            G.add_edges_from(dico_sol["aretes_LG"])
            nx.draw(G, node_size=500, with_labels=True);
            plt.savefig(path_save+"Graph_"+str(args["k"])+".png", format="PNG")
            plt.clf();
        #----------- pour debug, log, ..... -------------
    except Exception as e:
        print("####### EmptyDataError ", G_cpt, ": e = ", e," ####### ");
        df_debug.loc[len(df_debug.index)] = [G_cpt, len(args["aretes_G_nn_k"]), \
                        "error", "error", "error", "error", "error" ,\
                        "error", "error", "error", "error", "error",\
                        "error", "error", "error", "error","error", \
                        "error", "error", "error"];
                        
    simu50.save_df(df_debug, path_save, args["k"], headers_df)
    


if __name__ == '__main__':
    
    start= time.time();
    ### definition Matrice
#    m = G_nn(N=4, M=4);
    
    rep = "./graphe_cellules_test/"; 
    debug = True; debug = False; debug = True
    correl_seuil = 0.5                                                          # ne sert a rien mais urilise ds la fonction de cout
    mode_select_noeuds_1 = "aleatoire";
    critere_selection_pi1_pi2 = 0;                                              # 0: moins de modif,1: ajout aretes> supp aretes, 2:ajout aretes < supp aretes,
    number_permutations_nodes_1 = 25 #50;
    facteur_multiplicatifs = [1,2,5,10,20,40]; facteur_multiplicatifs = [1,2,5,10];
    facteur_multiplicatifs = [1,5,10]; facteur_multiplicatifs = [1,10];
    type_fct_couts = ["lineare_iourte_priorite_supp","lineare_iourte_priorite_ajout"];
    
    # construction de 50 graphes cellules k = {1,2,3,...,50}, k_min = 1
    k_max = 50; k_max = 25; k_min = 1; k_max = 50;  
    k_min = 2; k_max = 2;                                                      # A EFFACER
    graphes = generer_k_graphes(k_min, k_max, rep);
    params = list();
    for graphe in graphes:
        for type_fct_cout in type_fct_couts:
            for facteur_multiplicatif in facteur_multiplicatifs:
                exposant = 0; number_items_pi1_pi2 = 0.5#0.5;# 1;
                coef_fct_cout = (exposant, facteur_multiplicatif, type_fct_cout);
                path_save = rep+"graphe_cellule_prior_"+\
                            str(type_fct_cout.split("_")[3])+"_"+\
                            str(facteur_multiplicatif)+"/"
                args = {"k": graphe["k"],\
                    "M_G_nn": graphe["M_G_nn"], \
                    "aretes_G_nn_k":graphe["aretes_G_nn_k"], \
                    "dico_proba_cases": graphe["dico_proba_cases"],\
                    "dico_sommets_row_col": graphe["dico_sommets_row_col"],
                    "number_items_pi1_pi2": number_items_pi1_pi2, \
                    "number_permutations_nodes_1": number_permutations_nodes_1, \
                    "mode_select_noeuds_1": mode_select_noeuds_1,\
                    "coef_fct_cout": coef_fct_cout,\
                    "critere_selection_pi1_pi2": critere_selection_pi1_pi2,\
                    "facteur_mult": facteur_multiplicatif,\
                    "type_fct_cout": type_fct_cout,\
                    "correl_seuil":correl_seuil,\
                    "path_save": path_save, \
                    "debug":debug
                    }
                params.append([args])
                
    print("params = {}".format(len(params)))

    p = Pool(mp.cpu_count()-1) 
    p.starmap(simulation_graphe_cellule, params)
    p.terminate() 
           
    g = open("runtime_graphe_cellule_"+str(k_max)+".txt","w")
    ti = time.time() - start
    g.write(str(ti))
    g.close()
     
    print ("runtine = {}".format(time.time() - start))