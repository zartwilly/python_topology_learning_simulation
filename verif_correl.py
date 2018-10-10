#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 11:57:48 2017

@author: willy
Verif_correl
"""

import pandas as pd;
import numpy as np;
import itertools as it;
import re;
import fonctions_auxiliaires as fct_aux;
import time;

def caracterisation_arc(liste_colonnes):
    dico = dict()
    for arc in liste_colonnes:
        dico[arc] = {"nom":arc, "marque":None, "entree":list(), "sortie":list()}
    return dico
    
def ajouter_arc_marque_correct_tuples(list_arc_marke, liste_tuples, dico):
    '''
    tester BON
    ce script ajoute les arcs marques dans les combinaisons de tuples generees
    soit en entree lorsque l'arc est marque en sortie
    soit en sortie lorsque l'arc est marque en entree
    '''
    list_EnSortie = list()
    list_EnEntree = list()
    listes_tuples_complet = list()
    for arc in list_arc_marke:
        if dico[arc]['marque'] == "entree":
            list_EnSortie.append(arc)
        elif dico[arc]["marque"] == "sortie":
            list_EnEntree.append(arc)
    if len(liste_tuples) == 0:
        list_EnEntree_ = list_EnEntree.copy()
        list_EnSortie_ = list_EnSortie.copy()
        
        tupleE = list(set(list_EnEntree_))
        tupleS = list(set(list_EnSortie_))
        
        tuple_final = (tupleE, tupleS)
        listes_tuples_complet.append( tuple_final )
    else:
        for tuple_ in liste_tuples:
            tupleE = tuple_[0]
            tupleS = tuple_[1]
            list_EnEntree_ = list_EnEntree.copy()
            list_EnSortie_ = list_EnSortie.copy()
            if tupleE == None:
                tupleE = list_EnEntree_
            else:
                list_EnEntree_.extend(tupleE)
                tupleE = list(set(list_EnEntree_))
        
            if tupleS == None:
                tupleS = list_EnSortie_
            else:
                list_EnSortie_.extend(tupleS)
                tupleS = list(set(list_EnSortie_))
        
            tuple_final = (tupleE, tupleS)
            listes_tuples_complet.append( tuple_final )
    
    #print("liste_tuples_complet ", listes_tuples_complet)
    return listes_tuples_complet    
    
def combinaison(liste):
    '''
    tester BON
    liste = [1,2,3]
    return ([1,2],[3]), ([1,3],[2]), ....
    '''
    liste_tuple = list()
    taille = len(liste)
    if len(liste) < 1:
        #print (" error liste doit etre superieure ou egale a 1")
        pass
    else:
        if (None, liste[:]) not in liste_tuple and (liste[:], None) not in liste_tuple:
            liste_tuple.append( ([None], liste[:]) )
            liste_tuple.append( (liste[:], [None]) )
        
        '''
        complexite en 4*n-2
        '''
        for i in range(2,taille-1 ):
            tuple1 = (liste[:i], liste[i:]); tuple2 = (liste[i:], liste[:i]);
            if tuple1 not in liste_tuple:
                liste_tuple.append( tuple1 )
            if tuple2 not in liste_tuple:
                liste_tuple.append( tuple2 )
        '''
        complexite en n*n
        '''
        tmp = liste.copy(); s_tmp_ = set(tmp);
        for i in s_tmp_:
            s_tmp = set(tmp); s_tmp.remove(i);
            tuple1 = ([i], list(s_tmp)) if len(s_tmp) != 0 else ([i], [None]);
            tuple2 = (list(s_tmp),[i]) if len(s_tmp) != 0 else ([None], [i]);
            if tuple1 not in liste_tuple:
                liste_tuple.append( tuple1 )
            if tuple2 not in liste_tuple:
                liste_tuple.append( tuple2 )           
    return liste_tuple

def combinaison_generale(liste, dico):
    '''
    tester BON
    liste = [arc1, arc2, ...] : liste des arcs adjacents
    avec arc1 = dict("marque":{entre, sortie, None}, "entree": [], "sortie":[] )
    '''
    #print ('liste = ', liste)
    list_combinaison = list()
    list_arc_marke = list()
    liste_tuples = list()
    listes_tu = list()
    for arc in liste:
        if dico[arc]["marque"] != None :
            list_arc_marke.append( arc )
        else:
            list_combinaison.append( arc )
    
    liste_tuples = combinaison(list_combinaison)
    #print ("liste_tuples = ", liste_tuples)
    listes_tu = ajouter_arc_marque_correct_tuples(list_arc_marke, liste_tuples, dico)
    #print("listes_tu= ",listes_tu)
    return listes_tu

#def merger_dataframe(liste_grandeurs, chemin_dataset = "data/datasets/"):
#    '''
#    tester bon
#    forme un new dataset "df_fusion" avec certains columns de df_tmp
#    '''
#    df_fusion = pd.DataFrame()
#    for grandeur in liste_grandeurs:
#        df_tmp = pd.read_csv(chemin_dataset+'dataset_'+grandeur+'.csv')
#        df_fusion = pd.concat([df_fusion, df_tmp], axis=1)
#    return df_fusion
def check_grandeur_name_columns(cols,grandeur):
    """
    on verifie si grandeur est contenu dans les noms de grandeur 
    """
    for col in cols:
        if col.find("_"+grandeur) != -1:
            return True;
    return False;
    pass
def merger_dataframe(liste_grandeurs, chemin_dataset = "data/datasets/"):
    '''
    tester bon
    forme un new dataset "df_fusion" avec certains columns de df_tmp
    '''
    df_fusion = pd.DataFrame()
    for grandeur in liste_grandeurs:
        df_tmp = pd.read_csv(chemin_dataset+'dataset_'+grandeur+'.csv')
        cols = df_tmp.columns.tolist();
        dico_oldName_newName = dict();
        bool = check_grandeur_name_columns(cols, grandeur)
        if bool == False:
            for col in cols:
                dico_oldName_newName[col] = col+"_"+grandeur
            df_tmp.rename(columns= dico_oldName_newName, inplace = True)
        df_fusion = pd.concat([df_fusion, df_tmp], axis=1)
    return df_fusion
    
def compare_row(row, seuil):
    combinaisons = it.combinations(row.index.tolist(),2)
    liste_comparaison_val = []
    for col1, col2 in combinaisons:
        liste_comparaison_val.append( np.abs(row[col1] -row[col2]) )
    if round(np.max(liste_comparaison_val), 3) > seuil:
        return False
    else:
        return True

def VERIF_CORREL_grandeur(df_fusion, tuple_, grandeur, liste_grandeurs_enumeree, seuil_U = 10, epsilon = 0.75):
    """
    selon les elements du tuple tu_E[0], tu_S[1], selectionner les columnes de df_fusion
        * sommer( moyenne) les lignes par lignes 
        * mettre le resultat des sommes dans un dataset df_som (resp df_mean)
        * faire la comparaison element par element ==> (element1 - element2) < epsilon
        
    seuil_U: seuil de tension entre les arcs entrants et sortants d'un noeud. en general a 10
            seuil_U = 10
    """
    #### type de grandeur debut
    p = re.compile(".*P[pos,eng]*|.*I[1-3,12,23,13,21,32,31,n]{1}")
    liste_grandeurs_sommables = set()
    liste_grandeurs_non_sommables = set()
    for grandeur_ in liste_grandeurs_enumeree:
#        print("grandeur_: ",grandeur_);print(" p.match(grandeur_):",p.match(grandeur_))
        if p.match(grandeur_) != None:
            liste_grandeurs_sommables.add(grandeur_)
        else:
            liste_grandeurs_non_sommables.add(grandeur_)
    liste_grandeurs_non_sommables = list(liste_grandeurs_non_sommables)
    liste_grandeurs_sommables = list(liste_grandeurs_sommables)
    #### type de grandeur FIN
    
    liste_tuple_res  = list()
    tupleE = tuple_[0]
    tupleS = tuple_[1]
#    print('tuple_ =', tuple_)
    if None in tupleE:
        tupleE.remove(None)
    if None in tupleS:
        tupleS.remove(None)
    liste_arcs_ = tupleE + tupleS
#    print("liste_arcs_ :", liste_arcs_)
    liste_arcs = [x+"_"+grandeur for x in liste_arcs_]
    tupleE = [x+"_"+grandeur for x in tupleE]
    tupleS = [x+"_"+grandeur for x in tupleS]
    
    df_E_S = df_fusion[liste_arcs]
    
    if len(tupleE) != 0 and len(tupleS) != 0:
        if grandeur in liste_grandeurs_sommables:
            df_E_S['som_E'] = df_E_S[tupleE].sum(axis = 1)
            df_E_S['som_S'] = df_E_S[tupleS].sum(axis = 1)
            df_E_S['ecart_som'] = np.abs(df_E_S['som_S'] / df_E_S['som_E'])
            if ( df_E_S['ecart_som'].min() > epsilon and df_E_S['ecart_som'].min() <= 1 ):
                liste_tuple_res.append(tuple_)
        if grandeur in liste_grandeurs_non_sommables:
            df_E_S['mean_E'] = df_E_S[tupleE].mean(axis = 1)
            df_E_S['mean_S'] = df_E_S[tupleS].mean(axis = 1)
            df_E_S['ecart_mean'] = np.abs(df_E_S['mean_S'] / df_E_S['mean_E'])
            #print ( df_E_S[ ['ecart_mean','mean_E','mean_S'] ].head())
            
            if ( df_E_S['ecart_mean'].min() > epsilon and df_E_S['ecart_mean'].min() <= 1):
                liste_tuple_res.append(tuple_)
            
    elif (len(tupleE) != 0 and len(tupleS) == 0) :
        if grandeur in liste_grandeurs_non_sommables:
            res = df_E_S[tupleE].apply(lambda row: compare_row(row,seuil_U), axis = 1)
            if res.mean() > epsilon:
                liste_tuple_res.append(tuple_)
                
    elif (len(tupleE) == 0 and len(tupleS) != 0):
        if grandeur in liste_grandeurs_non_sommables:
            res = df_E_S[tupleS].apply(lambda row: compare_row(row,seuil_U), axis = 1)
            if res.mean() > epsilon:
                liste_tuple_res.append(tuple_)
         
    #### retour solution
    if len(liste_tuple_res) == 0:
        return False, list()
    elif len(liste_tuple_res) == 1:
        return True,liste_tuple_res[0]         

def VERIF_CORREL(l_tup_cu, df_fusion, seuil_U = 10, epsilon = 0.75 ,chemin_dataset = "data/datasets/"):
    """
    trouve le tuple verifiant les regles pour le maximum de grandeurs
    retrouve le tuple, un dico et un score
    """
#    print("**** l_tup_cu = ",l_tup_cu)
    liste_grandeurs_ = fct_aux.liste_grandeurs(chemin_dataset)
    liste_resultats = list()
    for grandeur in liste_grandeurs_:
        """ GROS PROBLEME ICI : => erreur de parametres PUIS ne calcul pas la correlation pour chaque tuple"""
#        resultat_oracle_grandeur = VERIF_CORREL_grandeur(df_fusion, grandeur, \
#                                                         liste_grandeurs_, l_tup_cu, \
#                                                         seuil_U, epsilon)
##                                                         seuil_U, epsilon, chemin_dataset)
#        liste_resultats.append(resultat_oracle_grandeur)
        #######
        for tup_cu in l_tup_cu:
            resultat_oracle_grandeur = VERIF_CORREL_grandeur(df_fusion, tup_cu, \
                                                             grandeur, liste_grandeurs_, \
                                                             seuil_U, epsilon)
            liste_resultats.append(resultat_oracle_grandeur)
        #######
        
    l_resultats_flatten = fct_aux.transform_nestedList_to_simpleList(liste_resultats)
    
    return l_resultats_flatten
    
def meilleur_score(l_resultats_flatten):
    """
    pour chaque tuple compte le nombre de fois qu'il apparait dans cette liste et ajouter dans un dico
    pour chaque valeur compter le nombre de cle ayant cette cle
    puis former la nouvelle valeur de cette cle 
    enfin retourner le 1er element du dico ayant la valeur max
    retourne tuple_max, n, k
    n: le nbre de bipartitions ayant verifier les regles au moins une grandeur
    k: le nombre de bipartitions ayant n
    """
    dico = dict(); tuple_ = list();
    if len(l_resultats_flatten):
        return tuple_, 0, 0;
        
    for tuple_ in l_resultats_flatten:
        tu0 = tuple_[0]; tu1 = tuple_[1];
        
        #print("____tu0 = ", tu0, "____tu1 = ", tu1)
        # creer cle
        key0 = None; key1 = None; key = None
        
        if len(tu0) != 0:
            sort_tu0 = sorted(tu0)
            key0 = ".".join(sort_tu0) 
        if len(tu1) != 0:
            sort_tu1 = sorted(tu1)
            key1 = ".".join(sort_tu1) 
        if key0 == None:
            key = "_"+key1
        elif key1 == None:
            key = key0+"_"
        else:
            key = key0+"_"+key1
            
        #ajout cle dans dico
        if key in dico:
            dico[key] += 1
        else:
            dico[key] = 1
    #recherche cle max
    dico_finale = dict()
    for k, v in dico.items():
        if v not in dico_finale:
            liste__= list()
            liste__.append(k)
            dico_finale[v] = liste__
        else:
            liste__= list()
            liste__= dico_finale[v]
            liste__.append(k)
            dico_finale[v] = liste__
    k_max = max( dico_finale.keys())
    #print("k_max=", k_max); print("dico_finale[k_max]= ", dico_finale[k_max])
    #print("dico_finale = ", dico_finale )
    valeur_k_max = dico_finale[k_max][0]

    # split valeur_k_max en 2 tuples
    s = valeur_k_max.split("_")
    tu0_max = s[0].split(".")
    tu1_max = s[1].split(".")
    tuple_max = (tu0_max,tu1_max)
    return tuple_max, k_max, len(valeur_k_max);
        
def noeud_concourant(cu, arguments_MAJ):
    """
    l_tup_cu = [([1,2],[3]), ([1,3],[2]), .... ]
    dico_scores ={ 4: [([1,2],[3]),...], 3:[([1,3],[2]),....]}
    """
    l_tup_cu = combinaison(list(cu));
#    print("cu: ",cu," l_tup_cu: ",l_tup_cu)
    resultats_flatten = VERIF_CORREL(l_tup_cu, arguments_MAJ["df_fusion"], \
                                     arguments_MAJ["seuil_U"], arguments_MAJ["epsilon"],\
                                     arguments_MAJ["chemin_dataset"])
    tuple_cu, k, n = meilleur_score(resultats_flatten)
    cu = set().union(*tuple_cu)
    aretes_cu1 = list();
    for noeud in cu:
        aretes_cu1.append( set(arguments_MAJ["dico_sommet_arete"][noeud]) )
    return set().union(*aretes_cu1).intersection(*aretes_cu1);
    
def choisir_best_cliques(subsets, arguments_MAJ):
    l_tup_subsets = list()
#    print("**** subsets: ",subsets)
    for cu in subsets:
        l_tup_subsets.extend(combinaison(list(cu)))
#    print("------ choisir_best_cliques l_tup_subsets: ", len(l_tup_subsets))
    resultats_flatten = list()
    resultats_flatten = VERIF_CORREL(l_tup_subsets, arguments_MAJ["df_fusion"], \
                                     arguments_MAJ["seuil_U"], arguments_MAJ["epsilon"],\
                                     arguments_MAJ["chemin_dataset"])
    tuple_cu, k, n = meilleur_score(resultats_flatten);
    cu = set().union(*tuple_cu);
    return cu;

def score_tuple(tup_cu, arguments_MAJ):
    """
    calcule le score de la tup_cu
    """
    pass
def meilleur_tuple(dico_scores):
    pass    
def noeud_concourant_old(cu, arguments_MAJ):
    """
    l_tup_cu = [([1,2],[3]), ([1,3],[2]), .... ]
    dico_scores ={ 4: [([1,2],[3]),...], 3:[([1,3],[2]),....]}
    """
    l_tup_cu = combinaison(list(cu));
    dico_scores = dict();
    for tup_cu in l_tup_cu:
        score_tup_cu = score_tuple(tup_cu, arguments_MAJ) 
        if score_tup_cu not in dico_scores.keys():
            dico_scores[score_tup_cu] = [tup_cu]
        else:
            l_same_scores = dico_scores[score_tup_cu];
            l_same_scores.append(tup_cu)
            dico_scores[score_tup_cu] = l_same_scores
    tuple_cu = meilleur_tuple(dico_scores)
    cu = set().union(*tuple_cu)
    aretes_cu1 = list();
    for noeud in cu:
        aretes_cu1.append( set(arguments_MAJ["dico_sommet_arete"][noeud]) )
    return set().union(*aretes_cu1).intersection(*aretes_cu1);


#### orientation des aretes DEBUT
def orientation_arete(liste_cliques, chemin_dataset, chemin_matE, epsilon):
    """
    but: orienter  les aretes du graphe non oriente trouve
    liste_cliques: liste des cliques par ordre croissant Cu etant une clique centre sur le noeud u
    
    1. fusionner les colonnes de toutes les grandeurs
    2. pour chaque clique generer les combinaisons
        2.1. appel du meta oracle ou VERIF-CORREL
    """
    pass
def construire_matA_intelligent(dico_arcs, sorted_cliques ):
    """
    tester ==> BON
    process :
        1 quel sont les cliques auquel appartient chaque arc
        2 # si cette arc appartient a 2 cliques 
          # mets une arete entre elles
          # sinon 
            # creer un noeud ext_arc 
            # mettre un arete entre ext_arc et clique
        3 transforme listes de couples en une liste
        4 construire la matrice d'adjacence des noeuds
        
    NB : sorted_cliques est rangee par ordre croissant 
    """
#### orientation des aretes FIN

if __name__ == '__main__':
    c = {1,2,3,4};
    l = []
    for i in range(len(c)):
        l.append(list(c)[i:])
    print("l=",l)