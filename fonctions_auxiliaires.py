#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 08:25:51 2016

ce fichier contient toutes les fonctions secondaires utilises dans tous les fichiers .py

@author: willy
"""

import pandas as pd
import numpy as np
import os
import re
from sklearn.metrics import mean_squared_error
from fastdtw import fastdtw
import itertools as it;
import json;

def is_file_exist(file_Path):
    """
    return True if file exist false otherwise
    """
    if os.path.isfile(file_Path) and os.access(file_Path, os.R_OK):
        return True;
    else:
        return False;
    
def range_2d(columns):
    treated = list();
    for row in columns:
        for col in columns:
            if row != col and (row,col) not in treated and (col,row) not in treated:
                treated.append( (row,col) )
                yield row, col;
                
def liste_grandeurs(chemin):
    files = os.listdir(chemin)
    liste_files = [fichier for fichier in files if  re.match("^dataset_.*",fichier)]
 
    pattern = re.compile('dataset_|.csv')        
    liste_grandeur = [pattern.sub("", mot_file) for mot_file in liste_files]
    return liste_grandeur

def remove_values_from_list(the_list, val):
   return [value for value in the_list if value != val]

def transform_nestedList_to_simpleList(liste_resultats):
    liste_res = list()
    for elt in liste_resultats:
        if type(elt) == list:
            tmp = [item for item in elt ]
            liste_res.extend (tmp)
        elif type(elt) == tuple:
            liste_res.append(elt)
    return liste_res
    
def export_from_json_to_dico(path_file):
    """
    exporter le json en dico;
    """
    return json.load(open(path_file,"r"));
    
def ajouter_correlation_dico(dico_aretes, arete, val_distance):
    """
    ajouter la valeur "val_distance" dans le dico "dico_aretes"
    """
    if arete in dico_aretes.keys():
        dico_aretes[arete].append(val_distance);
    else:
        dico_aretes[arete] = [val_distance];
    return dico_aretes        
        
def generer_tuple(liste):
    size = len(liste)
    list_tuple = list()
    if size == 0:
        return []
    else:
        for i in range(size):
            for j in range(i+1, size):
                list_tuple.append( (liste[i],liste[j]) )
    return list_tuple

def liste_arcs_old(matE):
    liste_cols = matE.columns.tolist()
    tmp = list()
    res = list()
    for row in liste_cols:
        for col in liste_cols:
            if row != col and matE.loc[row][col] == 1:
                if (col,row) not in tmp:
                    tmp.append( (col,row) )
                    tmp.append( (row,col) )
                    res.append( (row, col) )
    return res
def liste_arcs(mat):
    """ retourne la liste des arcs ou aretes d'un graphe. """
    res = list();
    for row, col in range_2d(mat.columns.tolist()):
        if mat.loc[row][col] == 1 or mat.loc[col][row] == 1:
            res.append((row, col))
    return res;
    
def liste_nonArcs(matE, k0):
    """
    but: trouver toutes les non aretes de matE cad matE.loc[x][y] == 0
    k0 = 0
    """
    liste_cols = matE.columns.tolist()
    tmp = list()
    res = list()
    for row in liste_cols:
        for col in liste_cols:
            if row != col and matE.loc[row][col] == k0:
                if (col,row) not in tmp:
                    tmp.append( (col,row) )
                    tmp.append( (row,col) )
                    res.append( (row, col) )
    return res
    pass

##### transforme liste d'aretes en une matrice d'adjacence
def transform_listeArcs_mat_adjacence(liste_noeuds_matIni, listeArcs, oriente):
    """
    but transformer une liste d'aretes en une matrice d'adjacence en utilisant 
        la liste des noeuds de la matrice initiale
        
    NB: liste_noeuds_matIni: liste des noeuds de la matrice initiale
        oriente est un booleen qui informe si ma matrice est oriente(directed) ou pas 
    """
    matE = pd.DataFrame( index = liste_noeuds_matIni, columns = liste_noeuds_matIni)
    ##print ("col matE ==== ", matE.columns.tolist())
    ##print ("ind matE ==== ", matE.index.tolist())
    
    if oriente == True:
        # on a un DAG
        for arc in listeArcs:
            matE.loc[ arc[0] ][arc[1]] = 1
    else:
        for arc in listeArcs:
            matE.loc[ arc[0] ][arc[1]] = 1
            matE.loc[ arc[1] ][arc[0]] = 1
    matE.fillna(0, inplace=True)
    return matE
##### transforme liste d'aretes en une matrice d'adjacence
def transform_list_matAdj(liste_aretes):
    """
    transforme liste_aretes par une matrice d'adjacence
    """
    l_noeuds = list(set([arete for tup in liste_aretes for arete in tup]))
    matE = pd.DataFrame( index = l_noeuds, columns = l_noeuds);
    for arc in liste_aretes:
        matE.loc[ arc[0] ][arc[1]] = 1;
        matE.loc[ arc[1] ][arc[0]] = 1;
    matE.fillna(0, inplace=True)
    return matE    
    pass

#### A EFFACER CAR PLUS utiliser
def voisins(liste_arcs, noeud):
    """
    recherche pour chaque arc si noeud est une extremite de cet arc
    """        
    liste_voisins = list()
    for arc in liste_arcs:
        if noeud == arc[0]:
            liste_voisins.append( arc[1] )
        if noeud == arc[1]:
            liste_voisins.append( arc[0] )
    return liste_voisins
#### A EFFACER CAR PLUS utiliser
def gamma_noeud_(matE):
    """
    but: determine le nombre de voisin d'un noeud et la liste des noeuds de son voisinage
    return dico
    dico = {"2":[3,{"1","3","4"},....]}
    """
    dico = dict()
    l_noeuds = matE.columns.tolist()
    for noeud in l_noeuds:
        ens = set([xj for xj in matE.columns.tolist() if matE.loc[noeud][xj] == 1 ])
        dico[noeud] = [len(ens), ens]
    return dico
def gamma_noeud(matE, liste_aretes):
    """
    but: determine le nombre de voisin d'un noeud et la liste des noeuds de son voisinage
    return dico
    """
    dico = dict()
    l_noeuds = matE.columns.tolist()
    for noeud in l_noeuds:
        ens = set()
        for arc in liste_aretes:
            if noeud == arc[0]:
                ens.add( arc[1] )
            if noeud == arc[1]:
                ens.add( arc[0] )
        dico[noeud] = [len(ens), ens]
    return dico   
def gamma(liste_arcs, noeud):
    """
    recherche le voisinage de "noeud"
    cad pour chaque arc si noeud est une extremite de cet arc
    """        
    liste_voisins = list()
    for arc in liste_arcs:
        if noeud == arc[0]:
            liste_voisins.append( arc[1] )
        if noeud == arc[1]:
            liste_voisins.append( arc[0] )
    return liste_voisins

def degre_noeud(liste_arcs, noeud):
    """
    retourne le nbre d arcs ayant un noeud en commun 
    """
    cpt = 0
    for arc in liste_arcs:
        if noeud == arc[0] or noeud == arc[1]:
           cpt += 1 
    return cpt
    
def listeAdjacence(matE, arc):
        
    liste = list()
    for elt in matE.columns.tolist():
        if elt != arc and matE.loc[arc][elt] == 1:
            liste.append(elt)
    return liste            
#    liste = list(set( [ei for ei in range(numeroEdge+1, matE[numeroEdge].shape[0]) \
#            if matE[numeroEdge][ei] == 1  ]))
#    return liste
def listeAdjacence_new(matE, arc, grandeur, split_boolean):
    liste = list()
    for elt in matE.columns.tolist():
        if elt != arc and matE.loc[arc][elt] == 1 and split_boolean == False:
            liste.append(elt)
        if elt != arc and matE.loc[arc][elt] == 1 and split_boolean == True:
            liste.append(elt.split(grandeur)[0])
    return liste 
    pass

def listeAdjacenceORACLE(matE, arc, dico):
        
    liste = list()
    for elt in matE.columns.tolist():
        if elt != arc and matE.loc[arc][elt] == 1:
            liste.append( dico[elt])
            if dico[arc] not in liste:
                liste.append(dico[arc])
    return liste      
    
def quicksort(arr):
    if len(arr) <= 1:
        return arr
    pivot = len( arr[len(arr) // 2] )
    left = [x for x in arr if len(x) < pivot]
    middle = [x for x in arr if len(x) == pivot]
    right = [x for x in arr if len(x) > pivot]
    return quicksort(left) + middle + quicksort(right)

###### calcull de similarite entre 2 arcs debut
def znorm(x):
    return (x - x.mean())/x.std()
        
def distance_euclidienne(x,y):
    if type(x) == list :
        x = np.asarray(x)
    if type(y) == list :
        y = np.asarray(y)
    return np.sqrt(np.sum(np.square(x-y)))


def distance_rmse(y_cible , y_pred):
    return mean_squared_error(y_cible, y_pred) ** 0.5

def calcul_distance(z_tuple0, z_tuple1, nom_distance):    
    if nom_distance == "euclidienne":
        return distance_euclidienne(z_tuple0, z_tuple1)
    elif nom_distance == "rmse":
        return distance_rmse(z_tuple0, z_tuple1)
    elif nom_distance == "dtw":
        return fastdtw(z_tuple0, z_tuple1, dist=distance_euclidienne)[0]

###### calcull de similarite entre 2 arcs fin    


###### coefficient de similarite debut
def strIntersection(s1, s2):
    """
    recherche et retourne l'element en commun dans les 2 strings
    """
    out = ""
    for c in s1:
        if c in s2 and not c in out:
            out += c
    return out
    
def trouver_arcs_entrant_sortant( noeud, liste_arcs_predits, dico_predit):
    liste_entrant_sortant = list()
    for arc in liste_arcs_predits:
        if noeud in arc:
            for k, v in dico_predit.items():
                if v == arc:
                    liste_entrant_sortant.append(k)
    
    return liste_entrant_sortant
    
###### coefficient de similarite  fin
   
###### generer combi pour C1_S1
def set_e(e, l):
    s = set();
    for item in range(len(e)):
        s = set().union(s).union(l[e[item]])
    return s
def x_(a):
    #a = [1,2,3,4]
    groupL = []
    for i in range(1,len(a)+1):
        groupL.extend(list(it.combinations(a,i)))
    L = []
    for e in list( it.combinations(groupL,2) ) :
        if len(set(e[0]).intersection(e[1])) == 0 and \
                len(set(e[0]).union(e[1]).intersection(a)) == len(a) :
            L.append(e)
    return L
def y1(l, C1_tuple1):
    # l = [{'A','B','C'},{'E'},{'F'}]
    # return  [({'A', 'D', 'F'}, {'E', 'G'}, []),({'E'}, {'A', 'D', 'F', 'G'}, []), ...]
    a = [i for i in range(len(l))]
    L = x_(a)
    L.append(( tuple(a),tuple()))
    r_l = []
    for e in L:
        set_e0 = set_e(e[0], l)
        set_e1 = set_e(e[1], l)
        r_l.append( (set(set_e0), set(set_e1), C1_tuple1) )
    return r_l
    
def generer_combi(C1_tuple, S1_):
    #S1_ = [{'E'}, {'G'}]
    #C1_tuple = ({'A', 'F', 'D'}, [{'A', 'F', 'D'}, set()])
    
    l = list(); r_l = list()
    if len(C1_tuple) == 0 :
        l.extend(S1_)
        r_l = y1(l, [])      
    else:
        l.extend(C1_tuple[0]);l.extend(S1_)
        r_l = y1(l, C1_tuple[1])
    #print("r_l = ", r_l)
    return r_l
    pass
###### generer combi pour C1_S1



def matrice_binaire(matE_proba, correl_seuil):
    """
    pour chaque correl (X) de matE_proba :
        * si X < correl_seuil ==> X = 0;
        * si X >= correl_seuil ==> X = 1;
    """
#    print("correl_seuil={}, type={}, \nmatE_proba={}\n"\
#          .format(correl_seuil,type(correl_seuil),matE_proba));
    matE_proba[matE_proba < correl_seuil] = 0;
#    print("inf matE_proba={}\n".format(matE_proba))
    matE_proba[matE_proba >= correl_seuil] = 1;
#    print("sup matE_proba={}\n".format(matE_proba))
    return matE_proba;