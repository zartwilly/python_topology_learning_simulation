#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 11:28:51 2016

@author: willy
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Created on Wed Oct 19 13:24:21 2016

@author: willy

ce fichier genere les mesures a partir du graphe fournir dans un fichier .csv
"""
import pandas as pd
import numpy as np
import clique_max as clique
#import oracle as ORACLE
import fonctions_auxiliaires as fct_aux
import time
import os
import generations_mesures as gene_mes
import json

def liste_arcs(matE):
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

def nommage_arcs(matA):
    """
    le but est de nommer les arcs entre 2 sommets de la maniere suivante:
        dico_dual_arc_sommet= {"A":("u","v"),"B":("u","w"),"C":("w","v"),"D":("z","w")}
    """
    dico_dual_arc_sommet = dict()
    dimensionAlphabet = 1000
    
    alphabet = [str(i) for i in range(dimensionAlphabet)]
    liste_arcs_ = liste_arcs(matA)
    for cpt in range(len(liste_arcs_)):
        dico_dual_arc_sommet[ alphabet[cpt] ] = liste_arcs_[cpt]
    return dico_dual_arc_sommet

def recherche_tuple(dico, tuple_):
    for k, v in dico.items():
        if v == tuple_:
            return k
    
def creation_matE_bis(dico, listeArcs):
    matE = pd.DataFrame(index = dico.keys(), columns = dico.keys())
    for lettre , tuple_ in dico.items():
        for arc_i in range( len(listeArcs) ):
            if lettre != listeArcs[arc_i]: 
                tu0 = tuple_[0]
                tu1 = tuple_[1]
                tu0_arc_i = listeArcs[arc_i][0]
                tu1_arc_i = listeArcs[arc_i][1]
                if tu0 == tu0_arc_i or tu0 == tu1_arc_i or \
                tu1 == tu0_arc_i or tu1 == tu1_arc_i:
                    lettre_trouve = recherche_tuple(dico, (tu0_arc_i, tu1_arc_i))
                    matE.loc[lettre][lettre_trouve] = 1
                    matE.loc[lettre_trouve][lettre] = 1
    matE.fillna( 0, inplace = True)
    for k in dico.keys():
        matE.loc[k][k] = 0
    return matE

def creation_matE(dico, listeArcs):
    matE = pd.DataFrame(index = dico.keys(), columns = dico.keys())
    for lettre , tuple_ in dico.items():
        for arc_i in range( len(listeArcs) ):
            tu0 = tuple_[0]
            tu1 = tuple_[1]
            tu0_arc_i = listeArcs[arc_i][0]
            tu1_arc_i = listeArcs[arc_i][1]
            if tu0 == tu0_arc_i or tu0 == tu1_arc_i or \
                tu1 == tu0_arc_i or tu1 == tu1_arc_i:
                lettre_trouve = recherche_tuple(dico, (tu0_arc_i, tu1_arc_i))
                if lettre != lettre_trouve:
                    matE.loc[lettre][lettre_trouve] = 1
                    matE.loc[lettre_trouve][lettre] = 1
    matE.fillna( 0, inplace = True)
    for k in dico.keys():
        matE.loc[k][k] = 0
        
    matE.index.rename("aretes",inplace=True)
    matE.loc["aretes"] = [str(i) for i in matE.index]
    return matE
        
def liste_predecesseurs(matE):
    '''
    but : definir la liste des predecesseurs de chaque noeud et le mettre ds un dico
            retour dico
    '''
    liste_col = matE.columns.tolist()
    dico = dict()
    for col in liste_col:
        liste_row = list()
        for row in liste_col:
            if matE.loc[row][col] == 1:
                liste_row.append(row)
                
        dico[col] = liste_row
    return dico

def liste_successeurs(matE):
    '''
    but : definir la liste des successeurs de chaque noeud et le mettre ds un dico
            retour dico
    '''
    dico = dict()
    liste_columns = matE.columns.tolist()
    for row in liste_columns:
        liste_col = list()
        for col in liste_columns:
            if matE.loc[row][col] == 1:
                liste_col.append(col)
                
        dico[row] = liste_col
    return dico

def genererMesures_ascendant_old(matA, dico_dual_arc_sommet, grandeur, taille = 3, effet_joule = 0.2):
    '''  tester bon ==> 
            ===> ERREUR lorsque 2 noeuds du meme niveau hierachiques ont une arete alors la propagation Puissance est faux

         propagation du flux des puits a la source
    
    taille : la dimension de chaque serie de mesures = 3000
    effet_joule: pourcentage de perte due a la resistance des cables. ici correspond a la difference entre l'entre et la sortie du noeud
    
    dico_pred : dictionnaire des predecesseurs avec cle : le noeud et valeurs: ses predecesseurs
    dico_succ : dictionnaire des successeurs avec cle : le noeud et valeurs: ses successeurs
    liste_grandeurs_sommables : liste des grandeurs dont les valeurs instantanees peuvent etre sommees (sum_entrant = sum_sortant => conservation de flux)
    intervalle_grandeur : dictionnaire contenant les valeurs inf et sup de chaque grandeur
    dico_grandeur : dictionnaire contenant les series de mesures pour chaque noeud( cle)
    liste_feuille : liste des noeuds n'ayant aucuns successeurs
    list_succ    : la liste des successeurs a un noeud
    
    liste_pourcentage = np.random.dirichlet(np.ones(nbre_pred),size=1)[0] : generation de nombre dont la somme est = a 1 et on recupere le 1er tableau ( [0] )
    '''
    dico_pred = liste_predecesseurs(matA)
    dico_succ = liste_successeurs(matA)
    liste_grandeurs_sommables = ["I12", "I23", "I31", "P", "S"]
    intervalle_grandeur = {"I12": (150, 200), "I23": (150, 200), "I31": (150, 200), "P": (33000, 62500),\
                           "U12": (220, 250), "U23": (220, 250), "U31": (220, 250), "S": (33000, 62500)}
    
    dico_grandeur = dict()
    dico_arc = dict()
    dico_grandeur_nom = dict()
    dico_arc_nom = dict()
    
    liste_cols = matA.columns.tolist()
    for noeud in liste_cols:
        if grandeur in liste_grandeurs_sommables:
            dico_grandeur[noeud] = np.zeros(taille)
        else: 
            dico_grandeur[noeud] = np.ones(taille)
    
    liste_arcs_ = liste_arcs(matA)
    for arc in liste_arcs_:
        if grandeur in liste_grandeurs_sommables:
            dico_arc[arc] = np.zeros(taille)
        else: 
            dico_arc[arc] = np.ones(taille)
            
    #recherche des puits de ce graphe
    liste_feuille = list()
    for cle, val in dico_succ.items():
        if len(val) == 0:
            liste_feuille.append(cle)
    ##print ("liste_feuille = ", liste_feuille)
        
    # generation de valeurs aux noeuds sources
    if len(liste_feuille) != 0:
        inf = intervalle_grandeur[grandeur][0] 
        sup = intervalle_grandeur[grandeur][1] 
        for feuille in liste_feuille:
            if  grandeur in liste_grandeurs_sommables:
                dico_grandeur[ feuille] = dico_grandeur[ feuille ] + np.random.uniform( inf, sup, taille)
            else:
                dico_grandeur[ feuille ] = dico_grandeur[ feuille ] * np.random.uniform( inf, sup, taille)
    ##print ("dico_grandeur = ", dico_grandeur)
    
    bool = True 
    while bool :
        noeud_feuille = liste_feuille.pop() # retirer les elts de liste_feuilles comme une file
        
        list_pred = dico_pred[noeud_feuille]
        ##print(" noeud_feuille ",noeud_feuille, "list_pred = ", list_pred)
        if len(list_pred) == 0:
            ##print("noeud source ", noeud_feuille)
            ##print (" ON PASSE AU NOEUD FEUILLE SUIVANT ")
            pass
        else:
            nbre_pred = len(list_pred)
            
            cpt = 0
            while len(list_pred) != 0 :
                pred = list_pred.pop()
                
                if grandeur in liste_grandeurs_sommables:
                    dico_grandeur[pred] += dico_grandeur[noeud_feuille] * (1/nbre_pred)*(1 - effet_joule)
                    arc_ = (noeud_feuille, pred)
                    dico_arc[arc_] += dico_grandeur[noeud_feuille] * (1/nbre_pred)*(1 - effet_joule)
                else:
                    dico_grandeur[pred] = dico_grandeur[ noeud_feuille ] *( 1 - effet_joule/nbre_pred )
                    arc_ = ( noeud_feuille, pred)
                    dico_arc[arc_] = dico_grandeur[noeud_feuille]*( 1 - effet_joule/nbre_pred )
                ##print ("dico_grandeur[", pred ,"] = ", dico_grandeur[pred])
                liste_feuille.insert(0, pred) # ajout de noeud successeur dans la queue de la file
                cpt += 1
        
        if len( liste_feuille ) == 0 :
            bool = False

    for cle, serie in dico_grandeur.items():
        nom_cle =  cle +"_"+grandeur
        dico_grandeur_nom[ nom_cle ] = serie
    for arc, tuple_ in dico_dual_arc_sommet.items():
        tuple_inv = (tuple_[1], tuple_[0])
        if tuple_ in dico_arc.keys():
            nom_cle =  arc +"_"+grandeur
            dico_arc_nom[ nom_cle ] = dico_arc[tuple_]
        if tuple_inv in dico_arc.keys():
            nom_cle =  arc +"_"+grandeur
            dico_arc_nom[ nom_cle ] = dico_arc[tuple_inv] 

    return dico_arc_nom

def classer_feuilles_pred(liste_feuille, dico_pred):
    """
    but : classer les feuilles de sorte que 
            * la feuille ayant le plus petit nbre de predecesseurs soit a l'indice 0 de liste_feuille_return
            * la feuille ayant le plus grand nbre de predecesseurs soit a l'indice n de liste_feuille_return
            
    test function:
        d_pred= {"A":2, "B":3, "C":4, "D":1, "E": 0}
        l_f = ["B","C","A","D","E"]
        l_f = ["D","E","B","C","A"]
        classer_feuilles_pred_bis(l_f, d_pred)
    """ 
    liste_feuille_return = list() # liste ordonne de mamniere croissante en fct du nombre de predecesseurs
    for feuille in liste_feuille:
        if len(liste_feuille_return) == 0:
            liste_feuille_return.append(feuille)
        else:      
            size_lfr = len(liste_feuille_return)         
            indice_avant = None                 # indice de l'elt de liste_feuille_return dont le nbre de predecesseurs est < au nbre de predecesseurs de feuille
            for cpt_liste in range(size_lfr):
                if dico_pred[feuille] < dico_pred[ liste_feuille_return[cpt_liste] ]:
                    indice_avant = cpt_liste
                    break
            if indice_avant != None:            # si indice_avant est non nulle alors on ajoute la feuille entre indice_avant et indice_avant+1 sinon on ajoute a la fin de liste_feuille_return
                liste_feuille_return = liste_feuille_return[:indice_avant] + [feuille] + liste_feuille_return[indice_avant:]
            else: 
                liste_feuille_return.append(feuille)

    ##print ("liste_feuille_return = ", liste_feuille_return)
    return liste_feuille_return


def genererMesures_ascendant(matA, dico_dual_arc_sommet, grandeur, taille = 3, effet_joule = 0.2):
    '''  tester bon ==> 
            ===> ERREUR lorsque 2 noeuds du meme niveau hierachiques ont une arete alors la propagation Puissance est faux

         propagation du flux des puits a la source
    
    taille : la dimension de chaque serie de mesures = 3000
    effet_joule: pourcentage de perte due a la resistance des cables. ici correspond a la difference entre l'entre et la sortie du noeud
    
    dico_pred : dictionnaire des predecesseurs avec cle : le noeud et valeurs: ses predecesseurs
    dico_succ : dictionnaire des successeurs avec cle : le noeud et valeurs: ses successeurs
    liste_grandeurs_sommables : liste des grandeurs dont les valeurs instantanees peuvent etre sommees (sum_entrant = sum_sortant => conservation de flux)
    intervalle_grandeur : dictionnaire contenant les valeurs inf et sup de chaque grandeur
    dico_grandeur : dictionnaire contenant les series de mesures pour chaque noeud( cle)
    liste_feuille : liste des noeuds n'ayant aucuns successeurs
    list_succ    : la liste des successeurs a un noeud
    
    liste_pourcentage = np.random.dirichlet(np.ones(nbre_pred),size=1)[0] : generation de nombre dont la somme est = a 1 et on recupere le 1er tableau ( [0] )
    '''
    dico_pred = liste_predecesseurs(matA)
    dico_succ = liste_successeurs(matA)
    liste_grandeurs_sommables = ["I12", "I23", "I31", "P", "S"]
    intervalle_grandeur = {"I12": (150, 200), "I23": (150, 200), "I31": (150, 200), "P": (33000, 62500),\
                           "U12": (220, 250), "U23": (220, 250), "U31": (220, 250), "S": (33000, 62500)}
    
    dico_grandeur = dict()
    dico_arc = dict()
    dico_grandeur_nom = dict()
    dico_arc_nom = dict()
    
    liste_cols = matA.columns.tolist()
    for noeud in liste_cols:
        if grandeur in liste_grandeurs_sommables:
            dico_grandeur[noeud] = np.zeros(taille)
        else: 
            dico_grandeur[noeud] = np.ones(taille)
    
    liste_arcs_ = liste_arcs(matA)
    for arc in liste_arcs_:
        if grandeur in liste_grandeurs_sommables:
            dico_arc[arc] = np.zeros(taille)
        else: 
            dico_arc[arc] = np.ones(taille)
            
    #recherche des puits de ce graphe et classer par ordre croissant 
    # par rapport au nbre de predecesseurs que chaque feuille possede
    # par exemple: soient A et B des feuilles ayant comme nbre de predecesseurs 1 et 2 resp.
    # alors liste_feuille = [A,B]
    liste_feuille = list()
    for cle, val in dico_succ.items():
        if len(val) == 0:
            liste_feuille.append(cle)
    liste_feuille = classer_feuilles_pred(liste_feuille, dico_pred)
    ##print ("liste_feuille = ", liste_feuille)
        
    # generation de valeurs aux noeuds sources
    if len(liste_feuille) != 0:
        inf = intervalle_grandeur[grandeur][0] 
        sup = intervalle_grandeur[grandeur][1] 
        for feuille in liste_feuille:
            if  grandeur in liste_grandeurs_sommables:
                dico_grandeur[ feuille] = dico_grandeur[ feuille ] + np.random.uniform( inf, sup, taille)
            else:
                dico_grandeur[ feuille ] = dico_grandeur[ feuille ] * np.random.uniform( inf, sup, taille)
    ##print ("dico_grandeur = ", dico_grandeur)
    
    bool = True 
    while bool :
        noeud_feuille = liste_feuille.pop() # retirer les elts de liste_feuilles comme une file
        
        list_pred = dico_pred[noeud_feuille]
        ##print(" noeud_feuille ",noeud_feuille, "list_pred = ", list_pred)
        if len(list_pred) == 0:
            #print("noeud source ", noeud_feuille)
            #print (" ON PASSE AU NOEUD FEUILLE SUIVANT ")
            pass
        else:
            nbre_pred = len(list_pred)
            
            cpt = 0
            list_pred = classer_feuilles_pred(list_pred, dico_pred)
            while len(list_pred) != 0 :
                pred = list_pred.pop()
                
                if grandeur in liste_grandeurs_sommables:
                    dico_grandeur[pred] += dico_grandeur[noeud_feuille] * (1/nbre_pred)*(1 - effet_joule)
                    arc_ = (noeud_feuille, pred)
                    dico_arc[arc_] += dico_grandeur[noeud_feuille] * (1/nbre_pred)*(1 - effet_joule)
                else:
                    dico_grandeur[pred] = dico_grandeur[ noeud_feuille ] *( 1 - effet_joule/nbre_pred )
                    arc_ = ( noeud_feuille, pred)
                    dico_arc[arc_] = dico_grandeur[noeud_feuille]*( 1 - effet_joule/nbre_pred )
                ##print ("dico_grandeur[", pred ,"] = ", dico_grandeur[pred])
                liste_feuille.insert(0, pred) # ajout de noeud successeur dans la queue de la file
                cpt += 1
        
        if len( liste_feuille ) == 0 :
            bool = False

    for cle, serie in dico_grandeur.items():
        nom_cle =  cle +"_"+grandeur
        dico_grandeur_nom[ nom_cle ] = serie
    for arc, tuple_ in dico_dual_arc_sommet.items():
        tuple_inv = (tuple_[1], tuple_[0])
        if tuple_ in dico_arc.keys():
            nom_cle =  arc +"_"+grandeur
            dico_arc_nom[ nom_cle ] = dico_arc[tuple_]
        if tuple_inv in dico_arc.keys():
            nom_cle =  arc +"_"+grandeur
            dico_arc_nom[ nom_cle ] = dico_arc[tuple_inv] 

    return dico_arc_nom


def genererMesures_descendant(matA, dico_dual_arc_sommet, grandeur, taille = 3, effet_joule = 0.2):
    '''  tester bon
     propagation du flux de la source aux puits
    
    taille : la dimension de chaque serie de mesures = 3000
    effet_joule: pourcentage de perte due a la resistance des cables. ici correspond a la difference entre l'entre et la sortie du noeud
    
    dico_pred : dictionnaire des predecesseurs avec cle : le noeud et valeurs: ses predecesseurs
    dico_succ : dictionnaire des successeurs avec cle : le noeud et valeurs: ses successeurs
    liste_grandeurs_sommables : liste des grandeurs dont les valeurs instantanees peuvent etre sommees (sum_entrant = sum_sortant => conservation de flux)
    intervalle_grandeur : dictionnaire contenant les valeurs inf et sup de chaque grandeur
    dico_grandeur : dictionnaire contenant les series de mesures pour chaque noeud( cle)
    liste_feuille : liste des noeuds n'ayant aucuns successeurs
    list_succ    : la liste des successeurs a un noeud
    
    liste_pourcentage = np.random.dirichlet(np.ones(nbre_pred),size=1)[0] : generation de nombre dont la somme est = a 1 et on recupere le 1er tableau ( [0] )
    '''
    dico_pred = liste_predecesseurs(matA)
    dico_succ = liste_successeurs(matA)
    liste_grandeurs_sommables = ["I12", "I23", "I31", "P", "S"]
    intervalle_grandeur = {"I12": (150, 200), "I23": (150, 200), "I31": (150, 200), "P": (33000, 62500),\
                           "U12": (220, 250), "U23": (220, 250), "U31": (220, 250), "S": (33000, 62500)}
    
    dico_grandeur = dict()
    dico_arc = dict()
    dico_grandeur_nom = dict()
    dico_arc_nom = dict()
    
    liste_cols = matA.columns.tolist()
    for noeud in liste_cols:
        if grandeur in liste_grandeurs_sommables:
            dico_grandeur[noeud] = np.zeros(taille)
        else: 
            dico_grandeur[noeud] = np.ones(taille)
    
    liste_arcs_ = liste_arcs(matA)
    for arc in liste_arcs_:
        if grandeur in liste_grandeurs_sommables:
            dico_arc[arc] = np.zeros(taille)
        else: 
            dico_arc[arc] = np.ones(taille)
            
    #recherche des sources de ce graphe
    liste_source = list()

    for cle, val in dico_pred.items():
        if len(val) == 0:
            liste_source.append(cle)
    ##print ("liste_source = ", liste_source)
        
    # generation de valeurs aux noeuds sources
    if len(liste_source) != 0:
        inf = intervalle_grandeur[grandeur][0] 
        sup = intervalle_grandeur[grandeur][1] 
        for source in liste_source:
            if  grandeur in liste_grandeurs_sommables:
                dico_grandeur[ source ] = dico_grandeur[ source ] + np.random.uniform( inf, sup, taille)
            else:
                dico_grandeur[ source ] = dico_grandeur[ source ] * np.random.uniform( inf, sup, taille)
#    #print ("dico_grandeur = ", dico_grandeur)
    
    bool = True 
    while bool :
        noeud_source = liste_source.pop() # retirer les elts de liste_feuilles comme une file
        
        list_succ = dico_succ[noeud_source]
#        #print(" noeud_source ",noeud_source, "list_succ = ", list_succ)
        if len(list_succ) == 0:
            ##print("noeud source ", noeud_source)
            ##print (" ON PASSE AU NOEUD FEUILLE SUIVANT ")
            pass
        else:
            nbre_succ = len(list_succ)
            
            cpt = 0
            while len(list_succ) != 0 :
                succ = list_succ.pop()
                
                if grandeur in liste_grandeurs_sommables:
                    dico_grandeur[succ] += dico_grandeur[noeud_source] * (1/nbre_succ)*(1 - effet_joule)
                    arc_ = (noeud_source, succ)
                    dico_arc[arc_] += dico_grandeur[noeud_source] * (1/nbre_succ)*(1 - effet_joule)
                else:
                    dico_grandeur[succ] = dico_grandeur[noeud_source] *( 1 - effet_joule/nbre_succ )
                    arc_ = (noeud_source, succ)
                    dico_arc[arc_] = dico_grandeur[noeud_source]*( 1 - effet_joule/nbre_succ )
#                #print ("dico_grandeur[", succ ,"] = ", dico_grandeur[succ])
                liste_source.insert(0, succ) # ajout de noeud successeur dans la queue de la file
                cpt += 1
        
        if len( liste_source ) == 0 :
            bool = False

    for cle, serie in dico_grandeur.items():
        nom_cle =  cle +"_"+grandeur
        dico_grandeur_nom[ nom_cle ] = serie
    for arc, tuple_ in dico_dual_arc_sommet.items():
        tuple_inv = (tuple_[1], tuple_[0])
        if tuple_ in dico_arc.keys():
            nom_cle =  arc +"_"+grandeur
            dico_arc_nom[ nom_cle ] = dico_arc[tuple_]
        if tuple_inv in dico_arc.keys():
            nom_cle =  arc +"_"+grandeur
            dico_arc_nom[ nom_cle ] = dico_arc[tuple_inv] 

    return dico_arc_nom

def genererMesures_all_grandeurs(matA, dico_dual_arc_sommet, liste_grandeurs, location = "data/datasets/", taille = 3, effet_joule = 0):
    '''
        generer pour chaque noeud des series de mesures associe aux differents grandeurs,
        puis la mettre ds un dico dont les cles sont les noms des arcs suivi de la grandeur like "XXX_U12", 
        et l'enregistrer dans des datasets en fonctions des grandeurs like "dataset_U12"
        
    
    taille : la dimension de chaque serie de mesures
    '''
    #effet_joule = 0 #0.2
    dico = dict()
    for grandeur in liste_grandeurs:
        dico[grandeur] = genererMesures_descendant(matA, dico_dual_arc_sommet, grandeur, taille, effet_joule)
        
        
        df = pd.DataFrame.from_dict(dico[grandeur], orient='columns', dtype=None)
        df.to_csv(location+"dataset_"+grandeur+".csv", index = False) 
        
        #print("grandeur "+ grandeur +" ...... termine")


def create_datasets(df_matA, dico_dual_arc_sommet, chemin_datasets, nbre_ts, effet_joule):
    """
    but :
        generer le dataset de mesures pour chaque grandeur
    particularite:
        df_matA est le graphe genere aleatoirement 
        
    nbre_ts : nombre de time series pour chaque grandeur
    """
    
    liste_grandeurs = ["I12", "I23", "I31", "P", "S", "U12", "U23", "U31"]
    genererMesures_all_grandeurs(df_matA, dico_dual_arc_sommet, liste_grandeurs, chemin_datasets, nbre_ts, effet_joule)

####### test generation mesures et datasets de mesures pour divers grandeurs ==> debut
def genererMesures_all_grandeurs_new(matA, 
                                     dico_dual_arc_sommet, 
                                     liste_grandeurs, 
                                     taille = 3, effet_joule = 0) :
    '''
        generer pour chaque noeud des series de mesures associe aux differents grandeurs,
        puis la mettre ds un dico dont les cles sont les noms des arcs suivi de la grandeur like "XXX_U12", 
        et l'enregistrer dans des datasets en fonctions des grandeurs like "dataset_U12"
        
    
    taille : la dimension de chaque serie de mesures
    '''
    #effet_joule = 0 #0.2
    dico = dict(); datasets = list();
    for grandeur in liste_grandeurs:
        dico[grandeur] = genererMesures_descendant(matA, dico_dual_arc_sommet, grandeur, taille, effet_joule)
        
        
        df = pd.DataFrame.from_dict(dico[grandeur], orient='columns', dtype=None)
        datasets.append((grandeur, df));
#        df.to_csv(location+"dataset_"+grandeur+".csv", index = False) 
    return datasets        

def create_datasets_new(df_matA, 
                        dico_dual_arc_sommet, 
                        nbre_ts, effet_joule):
    """
    but :
        generer le dataset de mesures pour chaque grandeur
    particularite:
        df_matA est le graphe genere aleatoirement 
        
    nbre_ts : nombre de time series pour chaque grandeur
    """
    
    liste_grandeurs = ["I12", "I23", "I31", "P", "S", "U12", "U23", "U31"]
    datasets = list()
    datasets = genererMesures_all_grandeurs_new(
                                            df_matA, 
                                            dico_dual_arc_sommet, 
                                            liste_grandeurs, 
                                            nbre_ts, effet_joule)
    return datasets;
####### test generation mesures et datasets de mesures pour divers grandeurs ==> fin
    
def grouper_liste(liste):
    dico = dict()
    for tup in liste:
        key = list()
        for l in tup:
            key.append("".join(l))
        key = "_".join(key)
        if key in dico:
            dico[key] += 1
        else:
            dico[key] = 1
    maxi = 0
    maxi_key = None
    for k, v in dico.items():
        if v > maxi:
            maxi = v
            maxi_key = k
    #print ("maxi_key = ", maxi_key)
    return dico
        
if __name__ == "__main__" :
    
    start= time.time()
     
    chemin_datasets = "data/datasets/"
    nbre_ts = 10
    effet_joule = 0.1
    

    matA = pd.read_csv("df_matA_generer.csv", index_col="nodes")
#    matA.set_index("nodes", inplace = True)
#    matA.index = [str(i) for i in matA.index]
    dico = nommage_arcs( matA )
    create_datasets(matA, dico, chemin_datasets, nbre_ts, effet_joule ) 
    #print("dico\n", dico)
    listeArcs = fct_aux.liste_arcs(matA)
    #print("listeArcs\n", listeArcs)
    matE = creation_matE(dico, listeArcs)
    
    #print(" matE = \n", matE)
    #print (time.time() - start) 
