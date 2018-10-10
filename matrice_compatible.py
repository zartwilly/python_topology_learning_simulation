# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 21:16:53 2016

@author: willy
but creer une matrice compatible 'matE_compa' de matE 
        tel que l'intersection de l'ensemble d'arcs adjacents de 2 arcs non adjacents est null cad
        si 2 arcs sont non adjacents alors leur ensemble d'arc adjacents est disjoint
        pour le faire, on a:
        * creation de matE_compa avec des 0
        * creation de tuples avec (ai,aj) ==(aj,ai). cela evite de parcourir tous la matrice
        * pour un tuple donne, on recupere l'ensemble d'arcs adjacents de 1er elt et 2nd elt du tuple
        * on fait un xor entre les 2 ensembles 
        * sil est egale a True on ajoute 1 dans matE_compa entre le 1er elt et le 2nd elt du tuple
            sinon on met 0
fonction principal : creer_matrice_compatible
"""

import numpy as np
import pandas as pd
import time
import sys 
sys.setrecursionlimit(10000) 


def xor(vect_ai, vect_aj):
    '''
        ===> code tester bon <====
    but ou exclusif entre vect_ai et vect_aj
        * faire un OR entre  vect_ai[cpt] et vect_aj[cpt] qu on met dans result
        * compte le nombre de 0 dans result 
            * si cpt_0 = len(result) return true => cela veut dire kon peut mettre un 1 entre les arcs ai et aj
        * compte le nombre de 1 dans result, vect_ai et vect_aj
            * si cpt_1_result == cpt_1_vect_ai + cpt_1_vect_aj return true => mettre un 1 entre les arcs ai et aj
            * si cpt_1_result != cpt_1_vect_ai + cpt_1_vect_aj return False => mettre un 0 entre les arcs ai et aj
    '''
    result = []
    for cpt in range( len(vect_ai) ):
        result.append( vect_ai[cpt] ^ vect_aj[cpt] )

    cpt_0_result = 0
    for k in result:
        if k == 0:
            cpt_0_result += 1
    if cpt_0_result == len(result):
        # tous les items de result sont = a 0
        return True
            
    cpt_1_vect_ai = 0 
    for i in vect_ai:
        if i == 1:
            cpt_1_vect_ai += 1
    cpt_1_vect_aj = 0
    for j in vect_aj:
        if j == 1:
            cpt_1_vect_aj += 1
    cpt_1_result = 0
    for k in result:
        if k == 1:
            cpt_1_result += 1
    if cpt_1_result == cpt_1_vect_ai + cpt_1_vect_aj:
        # |vect_ai XOR vect_aj | = |vect_ai| + |vect_aj|
        return True
    else:
        # |vect_ai XOR vect_aj | != |vect_ai| + |vect_aj|
        return False

def creer_tuples(matE):
    '''
        ===> code tester bon <====
    but : (a,b) tel que matA[a][b] == {0,1} ou matA[b][a] == {0,1}
    '''
    liste_tmp = []
    liste = []
    liste_cols = matE.columns.tolist()
    for row in liste_cols:
        for col in liste_cols:
            if col != row :
                if (col, row) not in liste_tmp:
                    liste_tmp.append( (row, col) )
                    liste.append( (row, col) )
    return liste
    
def creer_matrice_compatible(matE):
    '''
        ==> code tester bon <===
    but creer une matrice compatible 'matE_compa' de matE 
        tel que l'intersection de l'ensemble d'arcs adjacents de 2 arcs non adjacents est null cad
        si 2 arcs sont non adjacents alors leur ensemble d'arc adjacents est disjoint
        pour le faire, on a:
        * creation de matE_compa avec des 0
        * creation de tuples avec (ai,aj) == (aj,ai). cela evite de parcourir tous la matrice
        * pour un tuple donne, on recupere l'ensemble d'arcs adjacents de 1er elt et 2nd elt du tuple
        * on fait un xor entre les 2 ensembles 
        * sil est egale a True on ajoute 1 dans matE_compa entre le 1er elt et le 2nd elt du tuple
            sinon on met 0
    '''
    liste_cols = matE.columns.tolist()
    matE_compa = pd.DataFrame(index = liste_cols, columns = liste_cols)
    liste_tuple = creer_tuples(matE)
    for elt_tuple in liste_tuple:
        vec_ai = matE.loc[elt_tuple[0]]
        vec_aj = matE.loc[elt_tuple[1]]
        booleen = xor(vec_ai, vec_aj)
        if booleen :
            matE_compa.loc[ elt_tuple[0] ][ elt_tuple[1] ] = 1
            matE_compa.loc[ elt_tuple[1] ][ elt_tuple[0] ] = 1
        else:
            matE_compa.loc[ elt_tuple[0] ][ elt_tuple[1] ] = 0
            matE_compa.loc[ elt_tuple[1] ][ elt_tuple[0] ] = 0
    
    matE_compa.fillna(0, inplace = True) 
    return matE_compa.astype(int)
    
def creer_matrice_proba_compatible(matE_compa, matE_proba):
    '''
    matE_proba_compa est la matrice de probabilite associee a la matrice matE_compa. 
    en d'autres termes la proba entre 2 arcs est ajoutee a cette matrice ssi il existe une adjacence dans matE_compa
    matE_proba est la matrice de probabilites associe a une grandeur
    
    retourne matE_proba_compa
    '''

    liste_cols = matE_compa.columns.tolist()
    matE_proba_compa = pd.DataFrame(index = liste_cols, columns = liste_cols)
    liste_tuple = creer_tuples(matE_compa)
    for elt_tuple in liste_tuple:
        arc_ai = elt_tuple[0]
        arc_aj = elt_tuple[1]
        proba = matE_proba.loc[arc_ai][arc_aj]
        matE_proba_compa.loc[arc_ai][arc_aj] = proba
    
    matE_proba_compa.fillna(0, inplace = True)
    return matE_proba_compa.astype(float)
    
    