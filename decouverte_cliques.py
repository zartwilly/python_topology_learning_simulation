#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 18:12:05 2017
Created on Sat Jun 10 08:15:38 2017
Created on Mon Apr 17 09:19:11 2017

@author: willy
"""

from scipy.stats import norm
from scipy import stats, integrate

import pandas as pd
import numpy as np
import clique_max as clique
import verif_correl as VerifCorrel
import fonctions_auxiliaires as fct_aux
import generations_mesures as gene_mes
import correction_linegraph as corr_lineGraph

import networkx as nx
import time
import itertools as it
import random;
import math;
##### A EFFACER
import os
#import matplotlib.pyplot as plt
import simulation_50_graphes_PARALLELE as simu50
from pathlib import Path
##### A EFFACER

################## ambiguite et graphes doubles DEBUT #########################
def is_isomorphe_ambiguite(liste_matP):
    """
        A ECRIRE L ISOMORPHISME ENTRE 2 GRAPHES
    """
    """ correspond a la figure 2
    A = nx.Graph()
    B = nx.Graph()
    C = nx.Graph()
    D = nx.Graph()
    edgeA = [(1,2),(1,3),(2,3)]
    edgeB = [(1,2),(1,3),(2,3),(2,4),(4,3)]
    edgeC = [(1,2),(1,3),(1,4),(2,3),(3,4),(2,5),(3,5),(5,4)]
    edgeD = [(1,2),(1,4),(1,3),(1,5),(2,3),(2,4),(2,6),(6,4),(6,5),(6,3),(4,5),(5,3)]
    """
    A = nx.Graph()
    B = nx.Graph()
    C = nx.Graph()
    C_ = nx.Graph()
    D = nx.Graph()
    E = nx.Graph()
    E_ = nx.Graph()
    F = nx.Graph()
    edgeA = [(1,2),(2,3)]
    edgeB = [(1,2),(1,3),(2,3)]
    edgeD = [(1,2),(1,3),(2,3),(1,4)]
    edgeC = [(1,2),(1,3),(1,4,),(3,2)] # semble etre identique a edgeD
    edgeC_ = [(1,2),(1,3),(2,3),(1,4),(3,4)]
    #edgeE = [(1,2),(1,3),(2,3),(1,4),(1,5),(4,5)]
    #edgeE_ = [(1,2),(1,3),(2,3),(1,4),(1,5),(4,5),(2,5),(3,4)]
    #edgeF = [(1,2),(1,5),(2,5),(1,3),(1,4),(3,4)]
    
    A.add_edges_from(edgeA)
    B.add_edges_from(edgeB)
    C.add_edges_from(edgeC)
    C_.add_edges_from(edgeC_)
    D.add_edges_from(edgeD)
    #E.add_edges_from(edgeE)
    #E_.add_edges_from(edgeE_)
    #F.add_edges_from(edgeF)
    
    ### liste arcs matP et conversion nx.Graph()
    matP_ = nx.Graph()
    matP_.add_edges_from(liste_matP)
    ###
    
    ##print("isomatP, A", nx.is_isomorphic(matP_,A) == True)
    ##print("isomatP, B", nx.is_isomorphic(matP_,B) == True)
    ##print("isomatP, C", nx.is_isomorphic(matP_,C) == True)
    ##print("isomatP, C_", nx.is_isomorphic(matP_,C_) == True)
    ##print("isomatP, D", nx.is_isomorphic(matP_,D) == True)
    ##print("isomatP, E", nx.is_isomorphic(matP_,E) == True)
    ##print("isomatP, E_", nx.is_isomorphic(matP_,E_) == True)
    ##print("isomatP, F", nx.is_isomorphic(matP_,F) == True)
    
    if  nx.is_isomorphic(matP_,A) == True or \
        nx.is_isomorphic(matP_,B) == True or \
        nx.is_isomorphic(matP_,C) == True or \
        nx.is_isomorphic(matP_,C_) == True or \
        nx.is_isomorphic(matP_,D) == True :
        #nx.is_isomorphic(matP_,E) == True or \
        #nx.is_isomorphic(matP_,E_) == True or \
        #nx.is_isomorphic(matP_,F) == True  :
            return True
        
    return False
    

def is_isomorphe_graphe_double(liste_aretes_Ec):
    """
    but: determiner l'isomorphisme entre matP et les graphes doubles predefinis dans cette fonction
    retourne :
        True sil ya isomorphisme 
        False si il n'ya pas d'isomorphisme
    """
    A = nx.Graph()
    B = nx.Graph()
    C = nx.Graph()
    D = nx.Graph()
    edgeA = [(1,2),(1,3),(2,3)]
    edgeB = [(1,2),(1,3),(2,3),(4,2),(4,3)]
    edgeC = [(3,1),(3,2),(3,4),(3,5),(1,2),(1,4),(5,4),(5,2)]
    edgeD = [(1,2),(1,3),(1,4),(1,5),(6,2),(6,3),(6,4),(6,5),(2,3),(2,5),(4,5)]
    
    A.add_edges_from(edgeA)
    B.add_edges_from(edgeB)
    C.add_edges_from(edgeC)
    D.add_edges_from(edgeD)
    ### liste arcs matP et conversion nx.Graph()
    matP = nx.Graph()
    matP.add_edges_from(liste_aretes_Ec)
    ###
    
    if nx.is_isomorphic(matP,A) == True or \
        nx.is_isomorphic(matP,B) == True or \
        nx.is_isomorphic(matP,C) == True or \
        nx.is_isomorphic(matP,D) == True :
            return True
        
    return False
################## ambiguite et graphes doubles FIN ###########################

#####  fonctions secondaires --- DEBUT #################
def aretes_a_supprimer(liste_arcs, liste_Cu):
    for Cu in liste_Cu:
        #print(">>>0   Cu ", Cu)
        liste_tuples = fct_aux.generer_tuple(Cu)
        liste_arcs1 =  len(liste_arcs)
        #print(">>>1   Cu ", Cu)
        for tuple_ in liste_tuples:
            tuple_inverse = tuple( (tuple_[1], tuple_[0]) )
            if tuple_ in liste_arcs:
                liste_arcs.remove(tuple_)
            elif tuple_inverse in liste_arcs:
                liste_arcs.remove(tuple_inverse)
                #print("====>Suppression liste_arcs  AVANT = ",liste_arcs1," APRES = ", len(liste_arcs))  
    return liste_arcs  
    
def verifier_noeuds_belongs_2_cliques(ens_C, dico_cliq):
    """
    but: verifier si tous les noeuds de ma dico_cliq sont couvert par 1 ou 2 cliques max.
        si cest verifie alors return True
        sinon return False
    """
    for noeud in dico_cliq.keys():
        cpt = 0 
        for C in ens_C:
            if noeud in C:
                cpt += 1
        if cpt > 2:
            return False
    return True    
    
def liste_noeuds_1(liste_arcs, dico_cliq, ascendant_1):
    """
    but: retourne les listes des noeuds a -1 par ordre 
            croissant(ascendant) ou decroissant(descendant) de leur degre( voisinage)
    NB: ascendant_1 : definit l'ordre dans lequel les noeuds sont selectionnees
    """
    liste_1_sorted = list() #[c,b,a,e,d] ==> asc or [d,e,a,b,c] ==> desc
    dico_1 = dict()
    sorted_list_tuple = list() ### [('c', 1), ('b', 1), ('a', 3), ('e', 6), ('d', 10)]
    for noeud in dico_cliq.keys():
        if dico_cliq[noeud] == -1:
            degre_noeud = 0
            degre_noeud = fct_aux.degre_noeud(liste_arcs, noeud)
            dico_1[noeud] = degre_noeud

    if ascendant_1 == True:
        sorted_list_tuple =  sorted(dico_1.items(), key=lambda x: x[1], reverse=False)
    else:
        sorted_list_tuple =  sorted(dico_1.items(), key=lambda x: x[1], reverse=True) 
        
    for tup in sorted_list_tuple:
        liste_1_sorted.append(tup[0])
    
    return liste_1_sorted 
    
def liste_permutations_nodes_1(dico_gamma_noeud,dico_cliq, number_permutations): 
    nodes_1 = [node for node, label in dico_cliq.items() if label == -1]
    number_permutations = math.factorial(len(nodes_1)) \
    if math.factorial(len(nodes_1)) < number_permutations \
    else number_permutations;
    
    boolean = True; l_permutations = list()
    while boolean :
        tmp_nodes_1 = nodes_1.copy();
        une_permutation = list()
        while( len(tmp_nodes_1) !=0 ):
            node = random.choice(tmp_nodes_1)
            une_permutation.append(node)
            tmp_nodes_1.remove(node)
        if tuple(une_permutation) not in l_permutations:
           l_permutations.append(tuple(une_permutation));
           number_permutations += -1;
        boolean = False if number_permutations < 1 else True;
    
#    # ajout nodes_1 ordonnes soit par ordre croissant ou decroissant
#    l_tuples = list()
#    for node in nodes_1:
#        l_tuples.append( (node, dico_gamma_noeud[node][0]) )
#        
#    if croissant:
#        l_tuples = sorted(l_tuples, key = lambda x: x[1])
#    else :
#        l_tuples = sorted(l_tuples, key = lambda x: x[1], reverse = True)
#    
#    nodes_1_ordre_decroissant = [x[0] for x in l_tuples ]
##    print("**** nodes_1_ordre_decroissant= ", nodes_1_ordre_decroissant)
#    l_permutations.append(tuple(nodes_1_ordre_decroissant))
    return l_permutations;
    pass
    
def permutation_degre_min_max(dico_gamma_noeud, dico_cliq, croissant = False):
    """
    but retourne une permutaion dans les elements sont ordonnees par le degre. 
        ici on choisit le DEGRE MINIMUM
    """
    noeuds_1 = [x for x, v in dico_cliq.items() if v == -1]
    l_tuples = list()
    for node in noeuds_1:
        l_tuples.append( (node, dico_gamma_noeud[node][0]) )
    
    # ajout nodes_1 ordonnes soit par ordre croissant ou decroissant
    if croissant:
        l_tuples = sorted(l_tuples, key = lambda x: x[1])
    else :
        l_tuples = sorted(l_tuples, key = lambda x: x[1], reverse = True)
    
    nodes_1_ordre_decroissant = [x[0] for x in l_tuples ]
    return [tuple(nodes_1_ordre_decroissant)]
    pass
#####  fonctions secondaires --- FIN #################

###### degre minimum, sommet_tagge02_voisinage_clique, sommet_tagge0_voisinage_clique, voisinage_noeud -----> DEBUT ###########################################
def sommet_tagge0_voisinage_clique(dico_cliq, listeArcsEc, dico_gamma_noeud):  
    """ 
    fonction qui recherche tous les sommets dont cliq = {0,2} et Gamma forme une clique
    retourne :
        liste_sommets_cliq02_gamma : liste des sommets tel que cliq = {0,2} et Gamma forme une clique
        liste_sommets_NotCliq02_NotGamma : liste des sommets tel que cliq != {0,2} ou Gamma ne forme pas une clique

    dico_gamma_noeud: dico de degre de tous les noeuds et le set des noeuds voisins a "noeud"
    
    pour un test: 
    dico_gamma_noeud = {7:[2,{4,5}],1:[4,{0,2,4,3}],2:[3,{0,1,5}],3:[3,{1,4,6}],4:[4,{1,3,5,7}],5:[3,{2,4,7}],6:[1,{3}],0:[2,{1,2}]}
    dico_cliq = {7:0,1:0,2:0,3:0,4:0,5:0,6:0,0:0}
    listeArcsEc = [(1,0),(1,2),(0,2),(2,5),(4,5),(5,7),(4,7),(4,1),(4,3),(1,3),(3,6)]
    l0 = [1, 4, 2, 3, 5, 0, 7, 6]
    """
    liste_sommets_cliq02_gamma = list()
    liste_sommets_NotCliq02_NotGamma = list()
    
    list_noeuds_02 = [ noeud for noeud,val in dico_cliq.items() if val == 0 or val == 2]
    list_noeuds_diff_02 = [ noeud for noeud,val in dico_cliq.items() if val != 0 and val != 2]
    list_noeuds_02.sort(key = lambda noeud: dico_gamma_noeud[noeud], reverse = True)
    #print("+++++++++++ list_noeuds_02 : ", list_noeuds_02)
    while len(list_noeuds_02) != 0:
        #list_noeuds_0.sort(key = lambda noeud: dico_gamma_noeud[noeud], reverse = True)
        noeud = list_noeuds_02.pop()
        set_vois_noeud = dico_gamma_noeud[noeud][1]
        booleen = clique.is_clique(listeArcsEc, set_vois_noeud)
        if booleen == True:
            liste_sommets_cliq02_gamma.append(noeud)
            # pas besoin de supprimer tous les noeuds Voisins de "noeuds"
            # list_noeuds_0 = [ node for node in list_noeuds_0 if node not in set_vois_noeud  ]
        else:
            liste_sommets_NotCliq02_NotGamma.append(noeud)
            
    liste_sommets_NotCliq02_NotGamma.extend( list_noeuds_diff_02 )
    return set(liste_sommets_cliq02_gamma), set(liste_sommets_NotCliq02_NotGamma)
    
def sommet_tagge02_voisinage_clique(dico_cliq, listeArcsEc, dico_gamma_noeud):
    """
    liste des sommets labellise a 0 ou 2 et dont le voisinage forme une clique
    retourne 
        liste_sommets_cli02_gamma : liste des sommets tel que cliq = 0 ou 2 et Gamma forme une clique
        
    dico_gamma_noeud[noeud][1] = set de noeuds voisins a "noeud"
    """
    liste_sommets_cliq02_gamma = list()
    liste_sommets_cliq02_gamma = set([ noeud for noeud,val in dico_cliq.items() if (val == 0 or val == 2) and \
                          clique.is_clique(listeArcsEc, dico_gamma_noeud[noeud][1] )])
    
    return liste_sommets_cliq02_gamma
    
def degre_minimum(dico_gamma_noeud, l_sommet_notTag0_NotGama, dico_cliq, code):
    """
    but: recherche le noeud de degre minimum dans le graphe
    
    code in ["0","02"]
    min_deg : degre minimum
    res_noeud : noeud dont le degre est minimum
    dico_gamma_noeud = {7:[1,{5}],1:[3,{0,4,3}],2:[2,{0,5}],3:[3,{1,4,6}],4:[3,{1,3,5}],5:[3,{2,4,7}],6:[1,{3}],0:[2,{1,2}]}
        ===> dico_gamma_noeud[noeud=7][0] = 1
    """
    l_sommet_notTag0_NotGama.sort(key = lambda noeud: dico_gamma_noeud[noeud], reverse = True)
    min_deg = dico_gamma_noeud[l_sommet_notTag0_NotGama[0]][0]
    res_noeud = None
    for noeud in l_sommet_notTag0_NotGama:
        ##print("debug noeud=",noeud," min_deg= ",min_deg," degreNoeud= ", dico_gamma_noeud[noeud][0]," dico_cliq= ", dico_cliq[noeud]," res_noeud=", res_noeud)
        if min_deg >= dico_gamma_noeud[noeud][0] and dico_cliq[noeud] == 0 and code == "0":
            min_deg = dico_gamma_noeud[noeud][0]
            res_noeud = noeud
        if min_deg >= dico_gamma_noeud[noeud][0] and dico_cliq[noeud] in [0,2] and code == "02":
            min_deg = dico_gamma_noeud[noeud][0]
            res_noeud = noeud
    
    if res_noeud != None:  
        #print("debug res_noeud=",res_noeud," min_deg= ",min_deg," degreNoeud= ", dico_gamma_noeud[res_noeud][0]," dico_cliq= ", dico_cliq[res_noeud])
        pass
    else:
        #print("debug res_noeud=",res_noeud, " min_deg= 0")
        pass
    return res_noeud
    
def voisinage_noeud(listeArcs, liste_noeuds):
    """
    but : retourne le degree de tous les noeuds du graphe
    NB: dico = gamma_graph
    """
    dico = dict()
    for noeud in liste_noeuds:
        dico[noeud] = len(fct_aux.voisins(listeArcs, noeud) )
    return dico
###### degre minimum, sommet_tagge02_voisinage_clique, sommet_tagge0_voisinage_clique, voisinage_noeud ------> FIN ###########################################
    
########## reecriture MAJ DEBUT ###############################
def intersection_aretes_graphe( cu1, arguments_MAJ):
    """
    But:
        cas : simulation
        faire l intersection de toutes les aretes issues des noeuds de cu1 
            pour savoir si elles sont concourantes en un noeud du 
            graphe initial ou de base
        cas: dataReels
            A ECRIRE
            Recherche le noeud concourant des aretes issues des noeuds de cu1
    """
    if arguments_MAJ["simulation"] == True:
        aretes_cu1 = list();
        for noeud in cu1:
            aretes_cu1.append( set(arguments_MAJ["dico_sommet_arete"][noeud]) )
        #print("ICI aretes_cu1 = ", aretes_cu1)
        return set().union(*aretes_cu1).intersection(*aretes_cu1);
    else:
        s_noeud_inter = VerifCorrel.noeud_concourant(cu1, arguments_MAJ)
        return s_noeud_inter

def choisir_cliques(noeud_u, cu1, s_noeuds_inter, arguments_MAJ):
    """
    But:
        quel est le sous-ensemble de cu1 dont les aretes correspondantes concourent 
        en un noeud dans le graphe initial ou de base
        methodes:
            cas : simulation
            - les sous ensembles de taille len(cu1)-1 dont s_noeud_inter n'est pas inclu
            - trouver le sous-ensemble dont les aretes correspondantes forment un noeud 
                dans le graphe initial 
            
            cas: data reels
            - les sous ensembles de taille len(cu1)-1 dont s_noeud_inter n'est pas inclu
            - trouver le sous-ensemble dont les aretes correspondantes donne le score ou 
                proba la plus elevee. Ce score "eleve" indique la presence d'un noeud  
        
    """
    if arguments_MAJ["simulation"] == True:
        if len(cu1) == 3 and len(s_noeuds_inter) == 0:
            noeud_inter = intersection_aretes_graphe( cu1, arguments_MAJ)
            #print("ICI cu1=", cu1," noeud_inter=",noeud_inter," len=", len(noeud_inter))
            #if len(noeud_inter) != 0 and noeud_inter == {noeud_u}:
            if len(noeud_inter) != 0 and noeud_u in cu1:
                return cu1
            else:
                subsets = [set(c) for c in it.combinations(cu1, len(cu1)-1) ]
                for cliq in subsets:
                    noeud_inter = list();
                    noeud_inter = intersection_aretes_graphe( cliq, arguments_MAJ)
                    #print("noeud_inter: ", noeud_inter, " noeud_u: ", noeud_u, "cliq: ",cliq)
                    #if len(noeud_inter) != 0 and noeud_inter == {noeud_u}:
                    if len(noeud_inter) != 0 and noeud_u in cliq:
                        #print("cliq: ", cliq)
                        return cliq
    
        subsets = [set(c) for c in it.combinations(cu1, len(cu1)-1) if not s_noeuds_inter.issubset(c)]
        #print("subsets=", subsets, " s_noeuds_inter=", s_noeuds_inter)
        for cliq in subsets:
            noeud_inter = list();
            noeud_inter = intersection_aretes_graphe( cliq, arguments_MAJ)
            #print("noeud_inter: ", noeud_inter, " noeud_u: ", noeud_u, "cliq: ",cliq)
            #if len(noeud_inter) != 0 and noeud_inter == {noeud_u}:
            if len(noeud_inter) != 0 and noeud_u in cliq:
                #print("cliq: ", cliq)
                return cliq
    else:
        #### a REECRIRE car semble faux
        if len(cu1) == 3: 
            noeud_inter = intersection_aretes_graphe( cu1, arguments_MAJ)
            #print("noeud_inter: ", noeud_inter, " noeud_u: ", noeud_u)
            if len(noeud_inter) != 0 and noeud_inter == {noeud_u}:
                return cu1
        #### a REECRIRE car semble faux
        
#        print("choisir_cliques: cu1: ",cu1," noeud_u: ", noeud_u ," s_noeuds_inter = ",s_noeuds_inter)
        subsets = list();
        if len(s_noeuds_inter) <1:
            subsets = [set(c) for c in it.combinations(cu1, len(cu1)-1) if s_noeuds_inter.issubset(c)]
        else:
            subsets = [set(c) for c in it.combinations(cu1, len(cu1)-1) if not s_noeuds_inter.issubset(c)]
#        print("choisir_cliques: subsets:",subsets)
        cliq = VerifCorrel.choisir_best_cliques(subsets, arguments_MAJ)
        return cliq;
    return None
    
def MAJ_simulation(noeud_u, cu1, cu2, ver_u, matE, arguments_MAJ):
    """
    resout la MAJ pour data generes et pour data reels
    """
    s_noeuds_inter = cu1.intersection(cu2);
    if len(s_noeuds_inter) > 2 :
        return None, None, 0;
    elif len(s_noeuds_inter) < 1:
        # cas len(cu1) == 3
        cu1 = choisir_cliques(noeud_u,cu1, s_noeuds_inter, arguments_MAJ)
        if cu1 != None:
            return cu1, cu2, 1;
    elif len(s_noeuds_inter) == 1:
        noeud_inter_cu1 = intersection_aretes_graphe( cu1, arguments_MAJ);
        noeud_inter_cu2 = intersection_aretes_graphe( cu2, arguments_MAJ);
        if len(noeud_inter_cu1) == 1 and len(noeud_inter_cu2) == 1:
            return cu1, cu2, 1;
    elif len(s_noeuds_inter) == 2:
        noeud_inter_cu1 = intersection_aretes_graphe( cu1, arguments_MAJ)
        noeud_inter_cu2 = intersection_aretes_graphe( cu2, arguments_MAJ)
        if len(noeud_inter_cu1) == 0 and len(noeud_inter_cu2) == 0:
            cu1 = None; cu2 = None;
            #print("MAJ BIZARRE noeuds_inter = 2")
        elif len(noeud_inter_cu1) == 1 and len(noeud_inter_cu2) == 0:
            cu2 = choisir_cliques(noeud_u, cu2, s_noeuds_inter, arguments_MAJ)
        elif len(noeud_inter_cu1) == 0 and len(noeud_inter_cu2) == 1:
            cu1 = choisir_cliques(noeud_u, cu1, s_noeuds_inter, arguments_MAJ)
    
    if cu1 != None and cu2 != None:
        return cu1, cu2, 1
    else:
        return None, None, 0

def MAJ(noeud_u, tup_cliq, ver_u, matE, arguments_MAJ):
    """
    pas besoin de rechercher des isomorphismes
    """
    cu1 = set(); cu2 = set();
    if len(tup_cliq[0]) >= len(tup_cliq[1]):
        cu1 = tup_cliq[0]; cu2 = tup_cliq[1];
    else:
        cu1 = tup_cliq[1]; cu2 = tup_cliq[0];
    
    cu1, cu2, ver_u = MAJ_simulation(noeud_u, cu1, cu2, ver_u, matE, arguments_MAJ)
    return cu1, cu2, ver_u;

def MAJ_1clique(noeud_u, cu, matE, arguments_MAJ):
    if arguments_MAJ["simulation"] == True:
        """
        * On recherche les sous-ensembles dont 'noeud_u' est un element
        * et on cree des couples d ensembles (cu1,cu2)
        * pour chaque (cu1,cu2) donnee, on verifie si les aretes des noeuds de cu1 (resp cu2) 
            concourent en un point.
            - si cest le cas alors on return (cu1,cu2)
        * si aucune couple ne donne une bonne partition alors on retourne None, None, 0
        """
        l_paire_subsets = []
        if len(cu) == 2:
            arete_paire0 = set(arguments_MAJ["dico_sommet_arete"][list(cu)[0]]) 
            arete_paire1 = set(arguments_MAJ["dico_sommet_arete"][list(cu)[1]]) 
            print("paire0={},paire1={}".format(arete_paire0, arete_paire1))
            if len(arete_paire0.intersection(arete_paire1)) != 0:
                return cu, set(), 1;
        elif len(cu) > 2:
            subsets = [set(c) for c in it.combinations(cu, len(cu)-1) if {noeud_u}.issubset(c)]
            l_paire_subsets = [paire for paire in it.combinations(subsets, 2)]
            print("subsets={}".format(subsets))
            for paire in l_paire_subsets:
                aretes_paire0 = list();
                for noeud in paire[0]:
                    aretes_paire0.append( set(arguments_MAJ["dico_sommet_arete"][noeud]) )
                aretes_paire1 = list();
                for noeud in paire[1]:
                    aretes_paire1.append( set(arguments_MAJ["dico_sommet_arete"][noeud]) )
                noeud0_inter = set().union(*aretes_paire0).intersection(*aretes_paire0)
                noeud1_inter = set().union(*aretes_paire1).intersection(*aretes_paire1)
                if len(noeud0_inter) != 0 and len(noeud1_inter) != 0:
                    return paire[0], paire[1], 1;
        return None, None, 0;
    else:
        """
        A tester
        * On recherche les sous-ensembles dont 'noeud_u' est un element
        * utilisation de choisir_best_cliques pour determiner la clique cu1
        * ensuite suppression de cu1 dans la liste des sous-ensembles
        * utilisation encore de choisir_best_cliques pour determiner la clique cu2
        """
        if len(cu) < 3:
            #print("**** noeud_u:", noeud_u," Cu:",cu)
            return cu, set(),1;
#        print("MAJ_1clique: cu: ",cu," noeud_u: ", noeud_u )
        subsets = [set(c) for c in it.combinations(cu, len(cu)-1) if {noeud_u}.issubset(c)]
        cu1 = VerifCorrel.choisir_best_cliques(subsets, arguments_MAJ)
        #print("MAJ_1clique cu1 =",cu1)
        if cu1 != None and len(cu1)>0:
            subsets.remove(cu1)
        cu2 = VerifCorrel.choisir_best_cliques(subsets, arguments_MAJ)
        #print("noeud_u:", noeud_u," cu1:",cu1," cu2:",cu2, " Cu:",cu)
        return cu1, cu2, 1;
        
def MAJ_ambiguite(noeud_u, l_tup_cliq_coherentes, ver_u, matE, arguments_MAJ):
    """
    return la meilleur clique parmi les listes des cliques coherentes s'il peut la trouver
            dans le cas contraire, on laisse correction de cliques s'en charger.
    """
    for tup_cliq_coh in l_tup_cliq_coherentes:
        cu1 = None; cu2 = None; ver_u = 0;
        cu1, cu2, ver_u = MAJ_simulation(noeud_u, tup_cliq_coh[0], tup_cliq_coh[1], \
                                  ver_u, matE, arguments_MAJ)
        if cu1 != None and cu2 != None:
            return cu1, cu2, 1

    return None, None, 0;
########## reecriture MAJ FIN #################################

########## reecriture PARTITION DEBUT ###############################
def PARTITION(noeud_u, s_gamma_u, liste_aretes_Ec, dico_gamma_noeud, dico_cliq, matE, arguments_MAJ):
    """
    But: recherche la partition de s_gamma_u en 2 cliques cu1, cu2 tel que:
            * cu1.union(cu2)  = s_gamma_u
            * cu1.intersection(cu2) = {noeud_u} or {noeud_u, autre_noeud}
            * v in cu1 (resp cu2) a au max 1 voisin dans cu2(resp cu1) quand cliq(v) = 0 or
               v in cu1 (resp cu2) a au max 0 voisin dans cu2(resp cu1) quand cliq(v) = 2
               (coherence cliques cu1, cu2)
    """
    
    l_cliques = clique.find_clique(matE, s_gamma_u, [])
    l_cliques = [set(c) for c in l_cliques if noeud_u in c]
    #print("l_cliques  = ", l_cliques)
    if len(l_cliques) == 0:
        return None, None;
    if len(l_cliques) == 1:
        #print("ERROR PARTITION l_cliques")
        ver_u = 0;
        cu1, cu2, ver_u = MAJ_1clique(noeud_u, l_cliques[0], matE, arguments_MAJ);
#        return None, None; # commenter pour debug sur graphe test
        return cu1, cu2;

    
    # DEBUT ==== recherche bonne paire de cliques cad union(cu1,cu2) = s_gamma_u et coherente(cu1,cu2)
    l_tup_cliq_coherentes = [];
    for tup_cu1_cu2 in it.combinations(l_cliques, 2):
        noeud_inter_cu1_cu2 = tup_cu1_cu2[0].intersection(tup_cu1_cu2[1])
        if tup_cu1_cu2[0].union(tup_cu1_cu2[1]) == s_gamma_u and \
        is_cliques_coherente(noeud_u, tup_cu1_cu2, noeud_inter_cu1_cu2, dico_gamma_noeud, dico_cliq):
            l_tup_cliq_coherentes.append(tup_cu1_cu2)
#    print("l_tup_cliques_coherentes = ", l_tup_cliq_coherentes)
    # FIN ==== recherche bonne paire de cliques cad union(cu1,cu2) = s_gamma_u et coherente(cu1,cu2)
    
    if len(l_tup_cliq_coherentes) == 0:
#        print("Error tuple cliques coherentes")
        return None, None ;
    
    for tup_cliq in l_tup_cliq_coherentes:
        for i in range( len(tup_cliq) ):
            # recherche ambiguite ---> debut
            liste_aretes_P = list()
            for n1, n2 in it.combinations(tup_cliq[i],2):
                liste_aretes_P.append((n1,n2))
            if is_isomorphe_ambiguite(liste_aretes_P):
                ver_u = 0;
                #print(" = MAJ ambiguite tup_cliq : ", tup_cliq);
                #cu1, cu2, ver_u = MAJ(noeud_u, tup_cliq, ver_u, matE, arguments_MAJ)
                cu1, cu2, ver_u = MAJ_ambiguite(noeud_u, l_tup_cliq_coherentes, ver_u, matE, arguments_MAJ)
                return cu1, cu2;
            # recherche ambiguite ---> fin
            ver_u = 0;
            cu1, cu2, ver_u = MAJ(noeud_u, tup_cliq, ver_u, matE, arguments_MAJ)
            if ver_u == 1:
                return cu1, cu2;
    return None, None;
    
def is_cliques_coherente(noeud_u, tup_cu1_cu2, noeud_inter_cu1_cu2, dico_gamma_noeud, dico_cliq):
    """
    on verifie la regle de la clique coherente):
            * chaque sommet de Cu1 a au plus 1 voisin dans Cu2 si dico_cliq[noeud_u] = 0 et
            * chaque sommet de Cu1 a 0 voisin dans Cu2 si dico_cliq[noeud_u] = 2
    exple:
        dico_gamma_noeud = {"2":[3,{"1","3","4"}],....}
    """
    """cu1 = tup_cu1_cu2[0]; cu2 = tup_cu1_cu2[1]; """
    #print("tup_cu1_cu2 : ", tup_cu1_cu2, " noeud_inter_cu1_cu2 = ", noeud_inter_cu1_cu2)
    for noeud in tup_cu1_cu2[0] - noeud_inter_cu1_cu2:
        vois_noeud = dico_gamma_noeud[noeud][1];
        if dico_cliq[noeud] == 0 and len(vois_noeud.intersection(tup_cu1_cu2[1])-set(noeud_u)) > 1:
            #print("0 noeud: ", noeud, " incoherent 0 tuple=(", tup_cu1_cu2[0], tup_cu1_cu2[1],")","vois_noeud=",vois_noeud," inter_Vnode_tupC1C2S :", vois_noeud.intersection(tup_cu1_cu2[1]))
            return False
        if dico_cliq[noeud] == 2 and len(vois_noeud.intersection(tup_cu1_cu2[1])-set(noeud_u)) > 0:
            #print("1 noeud: ", noeud, " incoherent 2 tuple=", tup_cu1_cu2[0])
            return False
    for noeud in tup_cu1_cu2[1] - noeud_inter_cu1_cu2:
        vois_noeud = dico_gamma_noeud[noeud][1];
        if dico_cliq[noeud] == 0 and len(vois_noeud.intersection(tup_cu1_cu2[0])-set(noeud_u)) > 1:
            #print("2 noeud: ", noeud, " incoherent 0 tuple=", tup_cu1_cu2[1])
            return False
        if dico_cliq[noeud] == 2 and len(vois_noeud.intersection(tup_cu1_cu2[0])-set(noeud_u)) > 0:
            #print("3 noeud: ", noeud, " incoherent 2 tuple=", tup_cu1_cu2[1])
            return False
    return True
########## reecriture PARTITION FIN #################################

########## reecriture decouverte de cliques DEBUT ###############################
def couverture_en_cliques(dico_cliq, dico_gamma_noeud, liste_aretes_Ec, matE, 
                          dico_ver, arguments_MAJ):
    """
    algorithme de couverture en cliques.
    
    fonction qui determine le partitionnement en cliques de tous les sommets 
    du graphe en une ou deux cliques.
    
    Parametres nommes : 
    * dico_cliq -- un dictionnaire contenant les etats (0,1,2,3) de 
                        chaque sommet
    * dico_gamma_noeud -- dictionnaire contenant le voisinage de chaque sommet
    * liste_aretes_Ec -- liste des aretes du graphe 
    * matE -- matrice d'adjacence du graphe
    * dico_ver -- dictionnaire contenant NE SERT A RIEN
    * arguments_MAJ -- dictionnaire contenant les infos sur le reseau 
                        energetique cad la correspondante entre aretes et 
                        sommets, le seuil et les pertes par effetsjoules choisi
                        sur les mesures, les mesures dans un dataframe, etc 
        
    Elle retourne :
    * C -- la liste des cliques decouvertes dans le graphe
    * dico_cliq -- un dictionnaire contenant les etats (0,1,2,3) de 
                        chaque sommet
    * liste_aretes_Ec -- liste des aretes du graphe a la fin de la 
                                couverture. 
                                Si elle est vide, le graphe est un line-graphe.
    * ordre_noeuds_traites -- l'ordre dans lequel les sommets sont traites.
    * dico_sommets_par_cliqs -- dictionnaire contenant les cliques couvrant 
                                un sommet
    
    """
    C = list();
    ordre_noeuds_traites = list();
    dico_sommets_par_cliqs = dict();
    
    # initialisation du dictionnaire dico_sommets_par_cliqs
    for sommet in dico_cliq.keys():
        dico_sommets_par_cliqs[sommet] = list();
    
    ### while 
    while 0 in dico_cliq.values():
        
        l_sommets_tag02_gama, l_sommets_notTag02_NotGama = \
                      sommet_tagge0_voisinage_clique(dico_cliq, 
                                                     liste_aretes_Ec, 
                                                     dico_gamma_noeud)
#            print ("01 l_sommet_tag02_gama:",l_sommets_tag02_gama)
#            print ("01 l_sommet_notTag02_NotGama: ",l_sommets_notTag02_NotGama)              
        if len(l_sommets_tag02_gama) == 0:   
        #if len(l_sommets_notTag02_NotGama) != 0:  +=====>>> FAUX  FAUX 
            noeud_u = degre_minimum(dico_gamma_noeud, 
                                    list(l_sommets_notTag02_NotGama), 
                                    dico_cliq, code="0")                        # noeud de degre min dans l_sommet_notTag02_NotGama
#                print ("00 noeud_u :", noeud_u )
            if noeud_u is None:
                l_sommets_notTag02_NotGama = []
                continue;
            ordre_noeuds_traites.append(noeud_u)

            l_sommets_notTag02_NotGama.remove( noeud_u )
            gamma_u = dico_gamma_noeud[noeud_u][1]; gamma_u.add(noeud_u)
            #print ("01 noeud_u :", noeud_u ); #print ("01 gamma_u :", gamma_u)
            
            Cu1, Cu2 = PARTITION(noeud_u, gamma_u, liste_aretes_Ec, 
                                 dico_gamma_noeud, dico_cliq, 
                                 matE, arguments_MAJ)
            #print("01 noeud_u :", noeud_u ," PARTITION Cu1 = ", Cu1," Cu2 = ",Cu2,' gamma_u = ',gamma_u)
            
            if Cu1 is not None or Cu2 is not None:
                dico_cliq[noeud_u] = 1; dico_ver[noeud_u] = 1;
                #print("01 Cu1: ",Cu1," Cu2: ",Cu2)
                C.append(frozenset(Cu1)) if set(Cu1) not in C and \
                                            len(Cu1) != 0 else None
                C.append(frozenset(Cu2)) if set(Cu2) not in C and \
                                            len(Cu2) != 0 else None
                liste_aretes_Ec = aretes_a_supprimer(liste_aretes_Ec, 
                                                     [list(Cu1), list(Cu2)] )
            
                # recalcul gamma de chaque noeud ===> gamma_graph
                dico_gamma_noeud = fct_aux.gamma_noeud(matE, liste_aretes_Ec)
                for z in gamma_u-{noeud_u}:
                    if len(dico_gamma_noeud[z][1]) == 0 and dico_cliq[z] != -1:
                        dico_cliq[z] = 1;                                       #print("noeud_u :", noeud_u ," z= ", z, "1")
                    elif dico_cliq[z] == 0:                                    #or dico_cliq[z] == 2:
                         dico_cliq[z] = 2;                                      #print("noeud_u :", noeud_u ," z= ", z, "2")
                    else:
                         dico_cliq[z] = -1;                                     #print("01 traitement_Erreur1 noeud=",z)
                                  
            else:
                # il n'existe pas 2 partitions de U_gammaU
                # Cu1 == None ou Cu2 == None
                #print("01 il n'existe pas 2 partitions de U_gammaU noeud_1=", noeud_u, " a -1");
                dico_cliq[noeud_u] = -1;
          # fin if len(l_sommet_tag02_gama) == 0
        else:
#                print("****C = ",C)
            dico_gamma_noeud = fct_aux.gamma_noeud(matE, liste_aretes_Ec)
            while len(l_sommets_tag02_gama) != 0:
            #if len(l_sommets_tag02_gama) != 0:
                
                noeud_u = degre_minimum(dico_gamma_noeud, 
                                        list(l_sommets_tag02_gama), 
                                        dico_cliq, code="02")                   # noeud de degre min dans l_sommets_tag02_gama
                if noeud_u is None:
                    l_sommets_tag02_gama = []; continue;
                ordre_noeuds_traites.append(noeud_u)
                    
                l_sommets_tag02_gama.remove(noeud_u)
                
                gamma_u = dico_gamma_noeud[noeud_u][1]
                #gamma_u_comp = set(fct_aux.gamma(liste_aretes_Ec, noeud_u))
                #print("compare gamma_u-u_comp ", gamma_u == gamma_u_comp) # A EFFACER
                gamma_u.add(noeud_u)
                Cu = set(); Cu = gamma_u.copy();
#                    print("02 noeud_u: ",noeud_u, " gamma_u: ", gamma_u," Cu: ",Cu)
                
                if len(Cu) == 3:
                    Cu, Cu2,dico_ver[noeud_u] = MAJ(noeud_u, 
                                                    [set(Cu), set([])],
                                                    dico_ver[noeud_u], matE,
                                                    arguments_MAJ)
#                        print("ICI noeud_u: ",noeud_u, "Cu: ",Cu, " Cu2: ", Cu2)
                #print("sansPArtition noeud_u: ",noeud_u, "Cu: ",Cu, " gamma_u = ", gamma_u)

                dico_cliq[noeud_u] = 1; dico_ver[noeud_u] = 1;
                if len(Cu) != 1 and Cu not in C:
                    C.append(frozenset(Cu));
#                        print("ICI noeud_u: ",noeud_u, " add clique ",Cu)
                    
                #print(" -----1 liste_aretes_Ec: ", len(liste_aretes_Ec), " Cu: ",Cu )
                liste_aretes_Ec = aretes_a_supprimer(liste_aretes_Ec, 
                                                     [list(Cu)])
                #print(" -----2 liste_aretes_Ec: ", len(liste_aretes_Ec) )
                # recalcul gamma de chaque noeud ===> gamma_graph
                dico_gamma_noeud = fct_aux.gamma_noeud(matE, liste_aretes_Ec)
                
                for z in gamma_u-{noeud_u}:
                    #print("dico_cliq[",z,"]=", dico_cliq[z])
                    if len(dico_gamma_noeud[z][1]) == 0 and dico_cliq[z] != -1:
                        dico_cliq[z] = 1
                    elif dico_cliq[z] == 0 : # or dico_cliq[z] == 2:
                        dico_cliq[z] = 2
                    else:
                        dico_cliq[z] = -1
                        #print("02 traitement_Erreur z =", z)
                          
                l_sommets_tag02_gama = sommet_tagge02_voisinage_clique(
                                        dico_cliq, liste_aretes_Ec, 
                                        dico_gamma_noeud)
               
    ### while fin
    
    ## determiner couverture par sommets
    dico_sommets_par_cliqs = fct_aux.couverture_par_sommets(C)
            
            
    return C, dico_cliq, liste_aretes_Ec, ordre_noeuds_traites, \
            dico_sommets_par_cliqs;

def decouverte_cliques(matE, dico_sommet_arete, seuil_U=10, epsilon=0.75,
                        chemin_dataset="data/datasets/", 
                        chemin_matrices="data/matrices/",
                        ascendant_1=True, simulation=True, 
                        dico_proba_cases=dict(),
                        arg_params=dict()):
    """
    valeurs Parametres Par defaut: 
        seuil_U= 10, epsilon = 0.75, \
        chemin_dataset = "data/datasets/", chemin_matrices = "data/matrices/",\
        ascendant_1 = True, simulation = True, number_items_pi1_pi2 = 1, \
        number_permutations_nodes_1 = 10,
        
    BUT: recherche une couverture en cliques du graphe matE.
         Si matE est correct, aucun sommet n'est labellise par -1 cad dico_cliq.values != -1.
             on retourne la couverture la liste de toutes les cliques decouvertes
             variables returned :
                 C, E0, dico_cliq, som_cout_min
         si matE est incorrect, certains sommets sont labellises a -1
             on conserve les cliques decouvertes dans ens_C_i
             puis on execute la fonction correction_noeuds
             variables returned :
                 ens_C, liste_arcs, dico_cliq, som_cout_min
                 
    definition:
        C : ensemble de cliques. est la couverture en cliques
        ens_C: ensemble de cliques dont certains cliques sont labellises a -1
        liste_aretes = liste_arcs: ensemble de tous les aretes de matE
        E0 : liste aretes de matE. egale a liste_aretes. E0 ne subit aucune modification.
                elle est envoyee a la fonction "correction_noeuds"
        dico_cliq: dictionnaire des noeuds labellises, chaque noeud a trois labels: 0,1,2,-1
                    0: etat initiale. n'est couvert par aucune clique
                    1: est couvert par 2 cliques au maximum
                    2: est couvert par 1 clique. peut etre couvert par une autre clique
                    -1: est couvert par plus de 2 cliques.
        som_cout_min: en cas de correction de noeuds, on compte le nombre d'aretes supprimees 
                        et ajoutees.
        dico_sommet_arete: dico de dualite entre les sommets et leur aretes. 
                            sommet est associe a un linegraph, 
                            aretes est associe au reseau de flots non orientes
        epsilon: valeur a partir de laquel une clique correspond a un sommet.
                    est utilisee dans MAJ.
        seuil_U: ne sert pas actu
        chemin_dataset: chemin du repertoire contenant les datasets de chaque grandeur
        chemin_dataset: chemin du repertoire contenant les matrices d'adjacence (matE, matA)
        ascendant_1: definit l'ordre dans lequel les noeuds a -1 sont selectionnees.
                    s'il est a True alors cest du plus petit au plus grand
                    s'il est a True alors cest du plus grand au plus petit
        arg_params = {"number_permutations_nodes_1": 10, "biais": True, "algoGreedy":True, \
                  "mode_select_noeuds_1":"coutMin" or "degreMin" or "aleatoire", "number_items_pi1_pi2" = 1,\
                  "methode_delete_add_edges": 0, "coef_fct_cout":(1,1)}
    """
    #initialisation cliq et ver et C
    dico_cliq = dict(); dico_ver = dict(); C = list(); noeuds_corriges = list();
    dico_sommets_par_cliqs = dict();
    ordre_noeuds_traites = list();
    for sommet in matE.columns.tolist():   # nbre de noeuds dans le graphe
        dico_cliq[sommet] = 0; dico_ver[sommet] = 0
    
    # fusion des datasets 
    liste_grandeurs = fct_aux.liste_grandeurs(chemin_dataset)
    df_fusion = VerifCorrel.merger_dataframe(liste_grandeurs, chemin_dataset) 
    df_fusion.fillna(0, inplace = True)
    arguments_MAJ = {"dico_sommet_arete": dico_sommet_arete, 
                     "df_fusion": df_fusion, 
                     "seuil_U": seuil_U, "epsilon":epsilon, 
                     "chemin_dataset": chemin_dataset,
                     "simulation": simulation, "grandeurs": liste_grandeurs}
#    print("df_fusion: ", df_fusion.describe());
#    print("liste_grandeurs: ", liste_grandeurs)
#    print(" chemin_dataset: ",chemin_dataset);
#    print("chemin_matrices: ", chemin_matrices)                 
    # copy E0 <- Ec
    liste_aretes_Ec = fct_aux.liste_arcs(matE)
    dico_gamma_noeud = fct_aux.gamma_noeud(matE, liste_aretes_Ec) # {"2":[3,{"1","3","4"}],....}
    
    E0 = liste_aretes_Ec.copy()
    
    if is_isomorphe_graphe_double(liste_aretes_Ec) :
        """
        DEMANDER A DOMINIK CE QU"EST CE un GRAPHE DOUBLE ==> trouver
        """
        #print("le traiter avec Verif_correl ou ORACLE")
        return [], [], None, 0
    else:
        C, dico_cliq, liste_aretes_Ec, ordre_noeuds_traites, \
        dico_sommets_par_cliqs = couverture_en_cliques(dico_cliq.copy(), 
                                                     dico_gamma_noeud.copy(), 
                                                     liste_aretes_Ec.copy(), 
                                                     matE.copy(), 
                                                     dico_ver, 
                                                     arguments_MAJ.copy())
#    return C, dico_cliq, liste_aretes_Ec,\
#            ordre_noeuds_traites, dico_sommets_par_cliqs
#    return {0:[C, dico_cliq, E0, [], [], [], -11110,C]} # a DELETE
    print("Avant Correction liste_aretes_Ec = ", len(liste_aretes_Ec))#, "@@@@ C = ",C)
    som_cout_min = 0; 
    C_old = C.copy();
    if len(liste_aretes_Ec) != 0:
        #return C, liste_aretes_Ec, dico_cliq, 0
        #return [], [], dico_cliq #######====================############################### POUr TEST  trouver un matE qui n'est pas un linegraph
        matE.to_csv(chemin_matrices+'matE_notLineGraph.csv')
        dico_gamma_noeud = fct_aux.gamma_noeud(matE, E0) # {"2":[3,{"1","3","4"}],....}

        dico_permutations = dict(); 
        dico_permutations = solution_methode_nodes_1(dico_gamma_noeud, C.copy(), \
                             E0.copy(), ordre_noeuds_traites, dico_cliq, \
                             dico_proba_cases, arg_params)

        return dico_permutations
    else:
        return {0:[C, dico_cliq, E0, [], noeuds_corriges, ordre_noeuds_traites, som_cout_min, C_old]};             
    pass
#ens_C_i, dico_cliq, liste_aretes_E_C_i, min_cliques_aSupprDe_ens_C, noeuds_traites, som_cout_min;  
#####

########################## decouverte cliques nouveau #########################
def decouverte_cliques_new(matE, dico_sommet_arete, seuil_U=10, epsilon=0.75,
                        chemin_dataset="data/datasets/", 
                        chemin_matrices="data/matrices/",
                        ascendant_1=True, simulation=True, 
                        dico_proba_cases=dict(),
                        arg_params=dict()):
    """
    valeurs Parametres Par defaut: 
        seuil_U= 10, epsilon = 0.75, \
        chemin_dataset = "data/datasets/", chemin_matrices = "data/matrices/",\
        ascendant_1 = True, simulation = True, number_items_pi1_pi2 = 1, \
        number_permutations_nodes_1 = 10,
        
    BUT: recherche une couverture en cliques du graphe matE.
         Si matE est correct, aucun sommet n'est labellise par -1 cad dico_cliq.values != -1.
             on retourne la couverture la liste de toutes les cliques decouvertes
             variables returned :
                 C, E0, dico_cliq, som_cout_min
         si matE est incorrect, certains sommets sont labellises a -1
             on conserve les cliques decouvertes dans ens_C_i
             puis on execute la fonction correction_noeuds
             variables returned :
                 ens_C, liste_arcs, dico_cliq, som_cout_min
                 
    definition:
        C : ensemble de cliques. est la couverture en cliques
        ens_C: ensemble de cliques dont certains cliques sont labellises a -1
        liste_aretes = liste_arcs: ensemble de tous les aretes de matE
        E0 : liste aretes de matE. egale a liste_aretes. E0 ne subit aucune modification.
                elle est envoyee a la fonction "correction_noeuds"
        dico_cliq: dictionnaire des noeuds labellises, chaque noeud a trois labels: 0,1,2,-1
                    0: etat initiale. n'est couvert par aucune clique
                    1: est couvert par 2 cliques au maximum
                    2: est couvert par 1 clique. peut etre couvert par une autre clique
                    -1: est couvert par plus de 2 cliques.
        som_cout_min: en cas de correction de noeuds, on compte le nombre d'aretes supprimees 
                        et ajoutees.
        dico_sommet_arete: dico de dualite entre les sommets et leur aretes. 
                            sommet est associe a un linegraph, 
                            aretes est associe au reseau de flots non orientes
        epsilon: valeur a partir de laquel une clique correspond a un sommet.
                    est utilisee dans MAJ.
        seuil_U: ne sert pas actu
        chemin_dataset: chemin du repertoire contenant les datasets de chaque grandeur
        chemin_dataset: chemin du repertoire contenant les matrices d'adjacence (matE, matA)
        ascendant_1: definit l'ordre dans lequel les noeuds a -1 sont selectionnees.
                    s'il est a True alors cest du plus petit au plus grand
                    s'il est a True alors cest du plus grand au plus petit
        arg_params = {"number_permutations_nodes_1": 10, "biais": True, "algoGreedy":True, \
                  "mode_select_noeuds_1":"coutMin" or "degreMin" or "aleatoire", "number_items_pi1_pi2" = 1,\
                  "methode_delete_add_edges": 0, "coef_fct_cout":(1,1)}
    """
    #initialisation cliq et ver et C
    dico_cliq = dict(); dico_ver = dict(); C = list();
    dico_sommets_par_cliqs = dict();
    ordre_noeuds_traites = list();
    for sommet in matE.columns.tolist():   # nbre de noeuds dans le graphe
        dico_cliq[sommet] = 0; dico_ver[sommet] = 0
    
    # fusion des datasets 
    liste_grandeurs = fct_aux.liste_grandeurs(chemin_dataset)
    df_fusion = VerifCorrel.merger_dataframe(liste_grandeurs, chemin_dataset) 
    df_fusion.fillna(0, inplace = True)
    arguments_MAJ = {"dico_sommet_arete": dico_sommet_arete, 
                     "df_fusion": df_fusion, 
                     "seuil_U": seuil_U, "epsilon":epsilon, 
                     "chemin_dataset": chemin_dataset,
                     "simulation": simulation, "grandeurs": liste_grandeurs}
#    print("df_fusion: ", df_fusion.describe());
#    print("liste_grandeurs: ", liste_grandeurs)
#    print(" chemin_dataset: ",chemin_dataset);
#    print("chemin_matrices: ", chemin_matrices)                 
    # copy E0 <- Ec
    liste_aretes_Ec = fct_aux.liste_arcs(matE)
    dico_gamma_noeud = fct_aux.gamma_noeud(matE, liste_aretes_Ec) # {"2":[3,{"1","3","4"}],....}
    
    E0 = liste_aretes_Ec.copy()
    
    if is_isomorphe_graphe_double(liste_aretes_Ec) :
        """
        DEMANDER A DOMINIK CE QU"EST CE un GRAPHE DOUBLE ==> trouver
        """
        #print("le traiter avec Verif_correl ou ORACLE")
        return [], [], None, 0
    else:
        C, dico_cliq, liste_aretes_Ec, ordre_noeuds_traites, \
        dico_sommets_par_cliqs = couverture_en_cliques(dico_cliq.copy(), 
                                                     dico_gamma_noeud.copy(), 
                                                     liste_aretes_Ec.copy(), 
                                                     matE.copy(), 
                                                     dico_ver, 
                                                     arguments_MAJ.copy())
    return C, dico_cliq, liste_aretes_Ec,\
            ordre_noeuds_traites, dico_sommets_par_cliqs
#ens_C_i, dico_cliq, liste_aretes_E_C_i, min_cliques_aSupprDe_ens_C, noeuds_traites, som_cout_min;  

########################## decouverte cliques nouveau #########################

def solution_methode_nodes_1(dico_gamma_noeud,cliques, aretes_E0, ordre_noeuds_traites, \
                             dico_cliq, dico_proba_cases, arg_params):
    """
    return un dico_permutations correspondant a tous les solutions selon une methode
    """
    dico_permutations = dict()
#    print("00")
    if arg_params["mode_select_noeuds_1"] == "aleatoire":
#        print("10")
        l_permutations_1 = liste_permutations_nodes_1(\
                                                dico_gamma_noeud,dico_cliq, \
                                                arg_params["number_permutations_nodes_1"])
#        print("11")
        for une_permutation in l_permutations_1:
            dico_permutations[une_permutation] = \
            corr_lineGraph.correction_noeuds(une_permutation, cliques.copy(), aretes_E0.copy(), \
                                         dico_cliq, ordre_noeuds_traites, arg_params["number_items_pi1_pi2"], \
                                         dico_proba_cases, arg_params["coef_fct_cout"], \
                                         arg_params["correl_seuil"], arg_params["critere_selection_pi1_pi2"])
#        print("12")
    elif arg_params["mode_select_noeuds_1"] == "degreMin":
#        print("20")
        l_permut_degre_min = permutation_degre_min_max(dico_gamma_noeud, dico_cliq, False)
#        print("21")
        for permut_degre_min in l_permut_degre_min:
            dico_permutations[permut_degre_min] = \
                corr_lineGraph.correction_noeuds(permut_degre_min, cliques.copy(), aretes_E0.copy(), \
                                         dico_cliq, ordre_noeuds_traites, arg_params["number_items_pi1_pi2"], \
                                         dico_proba_cases, arg_params["coef_fct_cout"],\
                                         arg_params["correl_seuil"], arg_params["critere_selection_pi1_pi2"])
#        print("22")
    elif arg_params["mode_select_noeuds_1"] == "coutMin":
#        print("30")
        noeuds_1 = [x for x, v in dico_cliq.items() if v == -1];
#        print("31")
        solution = corr_lineGraph.corriger_noeuds_1(\
                            noeuds_1, cliques.copy(), aretes_E0.copy(),\
                            dico_cliq, ordre_noeuds_traites, \
                            arg_params["number_items_pi1_pi2"], dico_proba_cases,\
                            arg_params["coef_fct_cout"], arg_params["correl_seuil"],\
                            arg_params["critere_selection_pi1_pi2"])
        dico_permutations[solution[4]] = solution
#        print("32")
    else:
        print("..........METHOD DOESNT EXIST: methode de correction n existe pas...........")
    return dico_permutations
#####    

def decouverte_cliques_corriger_noeuds( matE, dico_sommet_arete, seuil_U= 10, epsilon = 0.75, \
                        chemin_dataset = "data/datasets/", chemin_matrices = "data/matrices/",\
                        ascendant_1 = True, simulation = True, arg_params = dict() ):
    """
    BUT: recherche une couverture en cliques du graphe matE.
         Si matE est correct, aucun sommet n'est labellise par -1 cad dico_cliq.values != -1.
             on retourne la couverture la liste de toutes les cliques decouvertes
             variables returned :
                 C, E0, dico_cliq, som_cout_min
         si matE est incorrect, certains sommets sont labellises a -1
             on conserve les cliques decouvertes dans ens_C_i
             puis on execute la fonction correction_noeuds
             variables returned :
                 ens_C, liste_arcs, dico_cliq, som_cout_min
                 
    definition:
        C : ensemble de cliques. est la couverture en cliques
        ens_C: ensemble de cliques dont certains cliques sont labellises a -1
        liste_aretes = liste_arcs: ensemble de tous les aretes de matE
        E0 : liste aretes de matE. egale a liste_aretes. E0 ne subit aucune modification.
                elle est envoyee a la fonction "correction_noeuds"
        dico_cliq: dictionnaire des noeuds labellises, chaque noeud a trois labels: 0,1,2,-1
                    0: etat initiale. n'est couvert par aucune clique
                    1: est couvert par 2 cliques au maximum
                    2: est couvert par 1 clique. peut etre couvert par une autre clique
                    -1: est couvert par plus de 2 cliques.
        som_cout_min: en cas de correction de noeuds, on compte le nombre d'aretes supprimees 
                        et ajoutees.
        dico_sommet_arete: dico de dualite entre les sommets et leur aretes. 
                            sommet est associe a un linegraph, 
                            aretes est associe au reseau de flots non orientes
        epsilon: valeur a partir de laquel une clique correspond a un sommet.
                    est utilisee dans MAJ.
        seuil_U: ne sert pas actu
        chemin_dataset: chemin du repertoire contenant les datasets de chaque grandeur
        chemin_dataset: chemin du repertoire contenant les matrices d'adjacence (matE, matA)
        ascendant_1: definit l'ordre dans lequel les noeuds a -1 sont selectionnees.
                    s'il est a True alors cest du plus petit au plus grand
                    s'il est a True alors cest du plus grand au plus petit
                    
        arg_params = {"number_permutations_nodes_1": 10(30, 100), "biais": True, "algoGreedy":True, \
                  "mode_select_noeuds_1":"coutMin" or "degreMin" or "aleatoire", \
                  "number_items_pi1_pi2" = 1, methode_deleted_add_edges = 0, \
                  "correl_seuil": correl_seuil}
    """
    boolean = True; som_cout_min = 0; noeuds_corriges = list(); number_noeud_1_min = pow(10,9);
    while (boolean):
        #initialisation cliq et ver et C
        dico_cliq = dict(); dico_ver = dict(); C = list();
        ordre_noeuds_traites = list();
        for sommet in matE.columns.tolist():   # nbre de noeuds dans le graphe
            dico_cliq[sommet] = 0; dico_ver[sommet] = 0
    
        # fusion des datasets 
        liste_grandeurs = fct_aux.liste_grandeurs(chemin_dataset)
        df_fusion = VerifCorrel.merger_dataframe(liste_grandeurs, chemin_dataset) 
        arguments_MAJ = {"dico_sommet_arete": dico_sommet_arete, "df_fusion": df_fusion, \
                         "seuil_U": seuil_U, "epsilon":epsilon, "chemin_dataset": chemin_dataset,\
                         "simulation": simulation, "grandeurs": liste_grandeurs}
                     
        # copy E0 <- Ec
        liste_aretes_Ec = fct_aux.liste_arcs(matE)
        dico_gamma_noeud = fct_aux.gamma_noeud(matE, liste_aretes_Ec) # {"2":[3,{"1","3","4"}],....}
    
        E0 = liste_aretes_Ec.copy()
    
        if is_isomorphe_graphe_double(liste_aretes_Ec) :
            """
            DEMANDER A DOMINIK CE QU"EST CE un GRAPHE DOUBLE ==> trouver
            """
            #print("le traiter avec Verif_correl ou ORACLE")
            return [], [], None, 0
        else:
            C, dico_cliq, liste_aretes_Ec, ordre_noeuds_traites = \
            couverture_en_cliques(dico_cliq.copy(), dico_gamma_noeud.copy(), liste_aretes_Ec.copy(), \
                                  matE.copy(), dico_ver, arguments_MAJ.copy())

        #print("Avant Correction liste_aretes_Ec = ", len(liste_aretes_Ec), "noeuds-1 = ",dico_cliq.values())
        number_noeuds_1 = sum( x == -1 for x in dico_cliq.values() )

        if -1 in dico_cliq.values() and number_noeud_1_min > number_noeuds_1:
            boolean = True; number_noeud_1_min = number_noeuds_1;
                
            matE.to_csv(chemin_matrices+'matE_notLineGraph.csv')
            #print("**** IMPORTANT: C: ", C, " dico_cliq: ", dico_cliq)
            list_noeuds_1 = liste_noeuds_1( E0, dico_cliq, ascendant_1);
            som_cout_z = 0;
            noeud_z, matE, cliques_a_ajouter, som_cout_z = \
            corr_lineGraph.corriger_noeud(list_noeuds_1, C, dico_cliq, \
                                          E0, matE, arg_params["number_items_pi1_pi2"],\
                                          arg_params["coef_fct_cout"], \
                                          arg_params["correl_seuil"], \
                                          arg_params["critere_selection_pi1_pi2"])
            som_cout_min += som_cout_z;
            noeuds_corriges.append(noeud_z)
            #print("1&&& boolean = True")
        elif -1 in dico_cliq.values() and number_noeud_1_min <= number_noeuds_1:
            boolean = False;
            matE.to_csv(chemin_matrices+'matE_notLineGraph.csv')
            list_noeuds_1 = liste_noeuds_1( E0, dico_cliq, ascendant_1);
            C, dico_cliq, aretes_Ec, cliques_aSupprDe_C, noeuds_corriges, som_cout_min = \
                    corr_lineGraph.correction_noeuds(list_noeuds_1, C, E0, dico_cliq, \
                                                     arg_params["number_items_pi1_pi2"],\
                                                     arg_params["coef_fct_cout"], \
                                                     arg_params["correl_seuil"],\
                                                     arg_params["critere_selection_pi1_pi2"])
            #print("2&&& boolean = false, number_noeud_1_min: ", number_noeud_1_min," <= number_noeuds_1: ",number_noeuds_1)
        else:
            boolean = False;
            #print("3&&& boolean = false")
                
    #print("C= ",C,"\n som_cout_min= ",som_cout_min,"\n noeuds_corriges= ",noeuds_corriges)           
    return C, E0, dico_cliq, som_cout_min, noeuds_corriges             
                    
#######^^^^^^^^^ new version 
def decouverte_cliques_corriger_noeuds_new( matE, dico_sommet_arete, seuil_U= 10, epsilon = 0.75, \
                        chemin_dataset = "data/datasets/", chemin_matrices = "data/matrices/",\
                        ascendant_1 = True, simulation = True, \
                        dico_proba_cases = dict(), arg_params = dict() ):
    """
    BUT: recherche une couverture en cliques du graphe matE.
         Si matE est correct, aucun sommet n'est labellise par -1 cad dico_cliq.values != -1.
             on retourne la couverture la liste de toutes les cliques decouvertes
             variables returned :
                 C, E0, dico_cliq, som_cout_min
         si matE est incorrect, certains sommets sont labellises a -1
             on conserve les cliques decouvertes dans ens_C_i
             puis on execute la fonction correction_noeuds
             variables returned :
                 ens_C, liste_arcs, dico_cliq, som_cout_min
                 
    definition:
        C : ensemble de cliques. est la couverture en cliques
        ens_C: ensemble de cliques dont certains cliques sont labellises a -1
        liste_aretes = liste_arcs: ensemble de tous les aretes de matE
        E0 : liste aretes de matE. egale a liste_aretes. E0 ne subit aucune modification.
                elle est envoyee a la fonction "correction_noeuds"
        dico_cliq: dictionnaire des noeuds labellises, chaque noeud a trois labels: 0,1,2,-1
                    0: etat initiale. n'est couvert par aucune clique
                    1: est couvert par 2 cliques au maximum
                    2: est couvert par 1 clique. peut etre couvert par une autre clique
                    -1: est couvert par plus de 2 cliques.
        som_cout_min: en cas de correction de noeuds, on compte le nombre d'aretes supprimees 
                        et ajoutees.
        dico_sommet_arete: dico de dualite entre les sommets et leur aretes. 
                            sommet est associe a un linegraph, 
                            aretes est associe au reseau de flots non orientes
        epsilon: valeur a partir de laquel une clique correspond a un sommet.
                    est utilisee dans MAJ.
        seuil_U: ne sert pas actu
        chemin_dataset: chemin du repertoire contenant les datasets de chaque grandeur
        chemin_dataset: chemin du repertoire contenant les matrices d'adjacence (matE, matA)
        ascendant_1: definit l'ordre dans lequel les noeuds a -1 sont selectionnees.
                    s'il est a True alors cest du plus petit au plus grand
                    s'il est a True alors cest du plus grand au plus petit
        arg_params = {"number_permutations_nodes_1": 10, "biais": True, "algoGreedy":True, \
                  "mode_select_noeuds_1":"coutMin" or "degreMin", "number_items_pi1_pi2" = 1,\
                  "methode_delete_add_edges": 0, "proba_seuil": proba_seuil}
    """
    boolean = True; som_cout_min = 0; noeuds_corriges = list(); number_noeud_1_min = pow(10,9);
    dico_permutations = dict();  cpt_gnrle = 0;
    while (boolean ):
        cpt_gnrle += 1;
        
        #initialisation cliq et ver et C
        dico_cliq = dict(); dico_ver = dict(); C = list();
        ordre_noeuds_traites = list();
        for sommet in matE.columns.tolist():   # nbre de noeuds dans le graphe
            dico_cliq[sommet] = 0; dico_ver[sommet] = 0;
    
        # fusion des datasets 
        liste_grandeurs = fct_aux.liste_grandeurs(chemin_dataset)
        df_fusion = VerifCorrel.merger_dataframe(liste_grandeurs, chemin_dataset) 
        arguments_MAJ = {"dico_sommet_arete": dico_sommet_arete, "df_fusion": df_fusion, \
                         "seuil_U": seuil_U, "epsilon":epsilon, "chemin_dataset": chemin_dataset,\
                         "simulation": simulation, "grandeurs": liste_grandeurs}
                     
        # copy E0 <- Ec
        liste_aretes_Ec = fct_aux.liste_arcs(matE)
        dico_gamma_noeud = fct_aux.gamma_noeud(matE, liste_aretes_Ec) # {"2":[3,{"1","3","4"}],....}
    
        E0 = liste_aretes_Ec.copy()
    
        if is_isomorphe_graphe_double(liste_aretes_Ec) :
            """
            DEMANDER A DOMINIK CE QU"EST CE un GRAPHE DOUBLE ==> trouver
            """
            #print("le traiter avec Verif_correl ou ORACLE")
            return [], [], None, 0
        else:
            C, dico_cliq, liste_aretes_Ec, ordre_noeuds_traites = \
            couverture_en_cliques(dico_cliq.copy(), dico_gamma_noeud.copy(), liste_aretes_Ec.copy(), \
                                  matE.copy(), dico_ver, arguments_MAJ.copy())

        #print("Avant Correction liste_aretes_Ec = ", len(liste_aretes_Ec), "noeuds-1 = ",dico_cliq.values())
        C_old = C.copy();
        number_noeuds_1 = sum( x == -1 for x in dico_cliq.values() )

        #print("number_noeud_1_min = ", number_noeud_1_min, ", number_noeuds_1 = ", number_noeuds_1)
        if -1 in dico_cliq.values() and number_noeud_1_min >= number_noeuds_1 and cpt_gnrle < 5:
            boolean = True; number_noeud_1_min = number_noeuds_1;
                
            matE.to_csv(chemin_matrices+'matE_notLineGraph.csv')
            list_noeuds_1 = liste_noeuds_1( E0, dico_cliq, ascendant_1);
            som_cout_z = 0;
            
            # selectionner le noeud_z de cout minimum ou de degre minimum
            critere_selection_noeuds_1 = ""
            if arg_params["mode_select_noeuds_1"] == "coutMin":
                critere_selection_noeuds_1 = "cout_minimum";
            elif arg_params["mode_select_noeuds_1"] == "degreMin":
                critere_selection_noeuds_1 = "degre_minimum";
            elif arg_params["mode_select_noeuds_1"] == "aleatoire":
                critere_selection_noeuds_1 = "aleatoire";
            else:
                critere_selection_noeuds_1 = ""
            ##### IMPORTANT a revoir correction + ajout choix cout min + degre min
            noeud_z, matE, cliques_a_ajouter, som_cout_z = \
            corr_lineGraph.corriger_noeud(list_noeuds_1, C, dico_cliq, \
                                              E0, matE.copy(), \
                                              arg_params["number_items_pi1_pi2"], \
                                              critere_selection_noeuds_1, \
                                              dico_proba_cases, arg_params["coef_fct_cout"], \
                                              arg_params["correl_seuil"],\
                                              arg_params["critere_selection_pi1_pi2"])
            som_cout_min += som_cout_z;
            noeuds_corriges.append(noeud_z)
            #print("noeud_z = ", noeud_z)
            #print("1&&& boolean = True")
        elif -1 in dico_cliq.values() and number_noeud_1_min < number_noeuds_1:
            boolean = False;
            matE.to_csv(chemin_matrices+'matE_notLineGraph.csv')
            list_noeuds_1 = liste_noeuds_1( E0, dico_cliq, ascendant_1);
            
            if arg_params["mode_select_noeuds_1"] == "coutMin":
                noeuds_1 = [x for x, v in dico_cliq.items() if v == -1]
                dico_permutations = dict(); 
                solution = corr_lineGraph.corriger_noeuds_1(\
                            noeuds_1, C, E0,\
                            dico_cliq, ordre_noeuds_traites, \
                            arg_params["number_items_pi1_pi2"], \
                            dico_proba_cases, arg_params["coef_fct_cout"], \
                            arg_params["correl_seuil"], \
                            arg_params["critere_selection_pi1_pi2"])
                dico_permutations[solution[4]] = solution;
            elif arg_params["mode_select_noeuds_1"] == "aleatoire":
                l_permutations_1 = liste_permutations_nodes_1(dico_gamma_noeud,dico_cliq, \
                                                          arg_params["number_permutations_nodes_1"])
                for une_permutation in l_permutations_1:
                    dico_permutations[une_permutation] = \
                    corr_lineGraph.correction_noeuds(une_permutation, C.copy(), E0.copy(), \
                                             dico_cliq, ordre_noeuds_traites, \
                                             arg_params["number_items_pi1_pi2"],\
                                             dico_proba_cases, arg_params["coef_fct_cout"],\
                                             arg_params["correl_seuil"],\
                                             arg_params["critere_selection_pi1_pi2"])
            elif arg_params["mode_select_noeuds_1"] == "degreMin":
                l_permut_degre_min = permutation_degre_min_max(dico_gamma_noeud, dico_cliq, False)
                for permut_degre_min in l_permut_degre_min:
                    dico_permutations[permut_degre_min] = \
                        corr_lineGraph.correction_noeuds(permut_degre_min, C.copy(), E0.copy(), \
                                             dico_cliq, ordre_noeuds_traites, arg_params["number_items_pi1_pi2"], \
                                             dico_proba_cases, arg_params["coef_fct_cout"],\
                                             arg_params["correl_seuil"], \
                                             arg_params["critere_selection_pi1_pi2"])
            else:
                dico_permutations = dict();
            
        else:
            boolean = False;
            dico_permutations = dict(); som_cout_min = 0;
            dico_permutations = {0:[C, dico_cliq, E0, [], noeuds_corriges, ordre_noeuds_traites, som_cout_min, C_old]}
            #print("3&&& boolean = false")
                
    #print("C= ",C,"\n som_cout_min= ",som_cout_min,"\n noeuds_corriges= ",noeuds_corriges)  
    return  dico_permutations;        
#######^^^^^^^^^ new version 
########## reecriture decouverte de cliques FIN #################################

################### debug decouverte_clique ###################
def test_decouverte_clique(dimMat, nb_lien, epsilon = 0.75, \
                           chemin_matrices= "data/matrices/", chemin_dataset = "data/datasets/", \
                           ascendant_1 = True, simulation = True, number_items_pi1_pi2 = 1, \
                           number_permutations_nodes_1 = 10 ):
    """
    but tester la fonction decouverte_cliques
    """
    test = "EN TEST" ; #test = "FINI"
    effet_joule = 0;
    seuil_U= 10;
    nbre_ts = 10
                                
    matE, matA, dico_dual_arc_sommet = simu50.matriceE(dimMat, nb_lien, chemin_matrices, \
                                                       chemin_dataset, nbre_ts, epsilon, \
                                                       effet_joule, test)
    print(" dico_dual_arc_sommet= ",dico_dual_arc_sommet)
    ### aretes a supprimer ==> voir data/G1
    print("----(avant)= mat[4,2]=", matE.loc['4']['2'])
    matE.loc['4']['2'] = 0; matE.loc['2']['4'] = 0;
    print("----(apres)= mat[2,4]=", matE.loc['2']['4'])
    ### aretes a supprimer ==> voir data/G1
    
    AAA = decouverte_cliques( matE, dico_dual_arc_sommet, seuil_U, epsilon, chemin_dataset, \
                             chemin_matrices, ascendant_1, simulation, number_items_pi1_pi2)
    for key in AAA.keys():
        print("C = ",AAA[key][0])
        print("dico_cliq = ",AAA[key][1])
        edge_list_matE =  fct_aux.liste_arcs(matE)
        # representation graphes
#        path_graphes_chemin = chemin_matrices +'representationGraphique/'
#        path_graphes = Path(path_graphes_chemin)
#        path_graphes.mkdir(parents=True, exist_ok=True)
#        #### matE
#        edge_list_matE =  fct_aux.liste_arcs(matE)
#        H_matE = nx.Graph(edge_list_matE); nx.draw(H_matE,with_labels=True);
#        plt.savefig(path_graphes_chemin+"matE.png"); plt.clf();
#        #### matA
#        edge_list_matA =  fct_aux.liste_arcs(matA)
#        H_matA = nx.Graph(edge_list_matA); nx.draw(H_matA,with_labels=True);
#        plt.savefig(path_graphes_chemin+"matA.png"); plt.clf();
#        print(" TEST dico_dual_arc_sommet = ", dico_dual_arc_sommet)

        # DH entre matE et LG
        edge_list_LG = simu50.aretes_C( AAA[key][0] )
        nbre_hamming, liste_arc_diff = simu50.distance_hamming(edge_list_matE, edge_list_LG)
        print("key :", key, " nbre_hamming: ", nbre_hamming); 
        print("key :", key, " liste_arc_diff: ",liste_arc_diff)
    
def test_decouverte_clique_corriger_noeuds(dimMat, nb_lien, epsilon = 0.75, \
                                           chemin_matrices= "data/matrices/", \
                                           chemin_dataset = "data/datasets/", \
                                           ascendant_1 = True, \
                                           nbre_aretes_Adelete = 1, \
                                           number_items_pi1_pi2= 1):
    """
    A MODIFIER
    but tester la fonction decouverte_cliques
    """
    test = "EN TEST" ; test = "FINI"; dimMat = 15
    simulation = True;
    effet_joule = 0
    seuil_U= 10;
    nbre_ts = 10
                                
    matE, matA, dico_dual_arc_sommet = simu50.matriceE(dimMat, nb_lien, chemin_matrices, \
                                                       chemin_dataset, nbre_ts, epsilon, \
                                                       effet_joule, test)

    ### suppression aretes dans matE ---> debut
    liste_aretes = fct_aux.liste_arcs(matE)
    for i in range(nbre_aretes_Adelete):
        random_nombre = random.randint(0,len(liste_aretes)-1);
        index_delete = liste_aretes[random_nombre]
        #print("----arete(s) a supprimer= ", index_delete)
        #print("----(avant)= mat[",index_delete[0],",",index_delete[1],"]=", matE.loc[index_delete[0]][index_delete[1]])
        matE.loc[index_delete[0]][index_delete[1]] = 0;
        matE.loc[index_delete[1]][index_delete[0]] = 0;
        #print("----(apres)= mat[",index_delete[0],",",index_delete[1],"]=", matE.loc[index_delete[0]][index_delete[1]])
    ### suppression aretes dans matE ---> fin
        
    AAA = decouverte_cliques_corriger_noeuds( matE, dico_dual_arc_sommet, seuil_U, epsilon, chemin_dataset, \
                             chemin_matrices, ascendant_1, simulation, number_items_pi1_pi2)
    
    
    print("C = ",AAA[0], " som_cout_min = ", AAA[3], "\n noeuds_corriges = ", AAA[4])
    
    
    # representation graphes
#    path_graphes_chemin = chemin_matrices +'representationGraphique/'
#    path_graphes = Path(path_graphes_chemin)
#    path_graphes.mkdir(parents=True, exist_ok=True)
#    #### matE
#    edge_list_matE =  fct_aux.liste_arcs(matE)
#    H_matE = nx.Graph(edge_list_matE); nx.draw(H_matE,with_labels=True);
#    plt.savefig(path_graphes_chemin+"matE.png"); plt.clf();
#    #### matA
#    edge_list_matA =  fct_aux.liste_arcs(matA)
#    H_matA = nx.Graph(edge_list_matA); nx.draw(H_matA,with_labels=True);
#    plt.savefig(path_graphes_chemin+"matA.png"); plt.clf();
#    #print(" TEST dico_dual_arc_sommet = ", dico_dual_arc_sommet)
    
    # DH entre matE et LG
    edge_list_LG = simu50.aretes_C( AAA[0] );
    edge_list_matE =  fct_aux.liste_arcs(matE);
    nbre_hamming, liste_arc_diff = simu50.distance_hamming(edge_list_matE, edge_list_LG)
    print("nbre_hamming: ", nbre_hamming); 
    print("liste_arc_diff: ",liste_arc_diff)

from pathlib import Path
def test_decouvClique_dataReal():
    """
    tester la decouverte de cliques sur des donnees reelles. 
    On appelle verif_correl pour verifier la loi de conservation
    """
    seuil_U= 10; epsilon = 0.75; nbre_ts = 10; effet_joule = 0;nbre_lien = 5; dimMat = 5;
    ascendant_1 = True; simulation = True;
    chemin_datasets = "data/datasets/"; chemin_matrices = "data/matrices/";
    dico_proba_cases= dict(); dico_sommet_arete = dict();
    methode_delete_add_edges = 0; proba_seuil = 0.5; SEUIL_PROBA = 0.8; 
    algoGreedy = False; biais = False;
    facteur_multiplicatif = 1; exposant = 1; coef_fct_cout = (exposant, facteur_multiplicatif)
    mode_select_noeuds_1 = "degreMin";
    arg_params = dict(); 
    arg_params = {"number_items_pi1_pi2": 1,"number_permutations_nodes_1": 10, \
                   "methode_delete_add_edges": methode_delete_add_edges, \
                   "biais": biais, \
                   "proba_seuil": proba_seuil,\
                   "SEUIL_PROBA": SEUIL_PROBA,\
                   "algoGreedy":algoGreedy, \
                   "mode_select_noeuds_1":mode_select_noeuds_1,
                   "coef_fct_cout":coef_fct_cout};
    ### path ===> debut
    path_chemin_matrices = Path( chemin_matrices)
    path_chemin_matrices.mkdir(parents=True, exist_ok=True)
    path_chemin_datasets = Path( chemin_datasets )
    path_chemin_datasets.mkdir(parents=True, exist_ok=True)
    ### path ===> FIN
    #### matE ====> debut
    test = "EN TEST" ; test = "FINI"; 
    matE, matA, dico_sommet_arete = simu50.matriceE(dimMat, nbre_lien, \
                                                    chemin_matrices, chemin_datasets, \
                                                    nbre_ts, epsilon, effet_joule, test)
    
    dico_proba_cases = simu50.ajouter_proba_matE(matE, {"ajouter":[],"supprimer":[]}, arg_params["SEUIL_PROBA"])
    matE.to_csv(chemin_matrices+"matE.csv")
    #### matE =====> fin
    
    dico_permutations = dict();
    dico_permutations = decouverte_cliques( matE, dico_sommet_arete, seuil_U, epsilon, \
                        chemin_dataset, chemin_matrices,\
                        ascendant_1, simulation, dico_proba_cases,\
                        arg_params)
    
    #### best solution ===> debut
    dico_sol = dict()
    dico_sol = simu50.best_permutation(dico_permutations, matE, matE);
    #### best solution ===> fin
    
    print("C : ", dico_sol['C'])
    pass
################### debug decouverte_clique ###################

#### algo couverture new version ####### ===> debut
def degre_minimum_new(noeuds_notProcess,dico_gamma_noeud, dico_cliq, code):
    """
    but: recherche le noeud de degre minimum dans le graphe
    
    code in ["0","02","03"];
    min_deg : degre minimum
    res_noeud : noeud dont le degre est minimum
    sommets_notTag0_NotGama: liste sommets tel que cliq != {0,2} ou Gamma ne forme pas une clique
        
    dico_gamma_noeud = {7:[1,{5}],1:[3,{0,4,3}],2:[2,{0,5}],3:[3,{1,4,6}],4:[3,{1,3,5}],5:[3,{2,4,7}],6:[1,{3}],0:[2,{1,2}]}
        ===> dico_gamma_noeud[noeud=7][0] = 1;
    
    """
    min_deg = dico_gamma_noeud[noeuds_notProcess[0]][0]
    res_noeud = noeuds_notProcess[0];
#    print("0 res_noeud={} min_deg={}".format(res_noeud,min_deg))
    for noeud in set(noeuds_notProcess)-set(res_noeud):
#        print("debug noeud=",noeud," min_deg= ",min_deg," degreNoeud= ", dico_gamma_noeud[noeud][0]," dico_cliq= ", dico_cliq[noeud]," res_noeud=", res_noeud)
        if min_deg >= dico_gamma_noeud[noeud][0] and dico_cliq[noeud] == 0 and code == "0":
            min_deg = dico_gamma_noeud[noeud][0]
            res_noeud = noeud
        if min_deg >= dico_gamma_noeud[noeud][0] and dico_cliq[noeud] in [0,2] and code == "02":
            min_deg = dico_gamma_noeud[noeud][0]
            res_noeud = noeud
        if min_deg >= dico_gamma_noeud[noeud][0] and dico_cliq[noeud] in [0,3] and code == "03":
            min_deg = dico_gamma_noeud[noeud][0]
            res_noeud = noeud
#    print("1 res_noeud={}, min_deg={}".format(res_noeud,min_deg))
    return res_noeud
    
def aretes_cliques(cliques):
    """
    return edges of each clique belonging to cliques
    """
    aretes = set()
    for cliq in cliques:
        aretes = aretes.union( set(it.combinations(cliq,2)) )
    return list(aretes);
    
def cover_new_version(dico_cliq, dico_gamma_noeud, liste_aretes_Ec, matE, dico_ver, arguments_MAJ):
    """
    arguments_MAJ = {"dico_sommet_arete": dico_sommet_arete, "df_fusion": df_fusion, \
                     "seuil_U": seuil_U, "epsilon":epsilon, "chemin_dataset": chemin_dataset,\
                     "simulation": simulation, "grandeurs": liste_grandeurs}
    """
    C = list(); ordre_noeuds_traites = list();noeuds_notProcess = list(dico_cliq.keys())
    while 0 in dico_cliq.values() or 3 in dico_cliq.values():
        noeud_u = degre_minimum_new(noeuds_notProcess,dico_gamma_noeud, dico_cliq, code="03")
#        if noeud_u == None:
#            sommets_notTag03_NotGama = [];
#            continue;
        ordre_noeuds_traites.append(noeud_u);
        noeuds_notProcess.remove( noeud_u )
        gamma_u = dico_gamma_noeud[noeud_u][1]; gamma_u.add(noeud_u)
        #print ("01 noeud_u :", noeud_u ); #print ("01 gamma_u :", gamma_u)
        
        Cu1, Cu2 = PARTITION(noeud_u, gamma_u, liste_aretes_Ec, dico_gamma_noeud,\
                             dico_cliq, matE, arguments_MAJ);
        print("01 noeud_u :", noeud_u ," PARTITION Cu1 = ", Cu1," Cu2 = ",Cu2,' gamma_u = ',gamma_u)
        
        if (Cu1 != None and Cu2 != None) or \
            (Cu1 != None and len(Cu2)==0 and dico_cliq[noeud_u]==3) or \
            (Cu2 != None and len(Cu1)==0 and dico_cliq[noeud_u]==3):
            if dico_cliq[noeud_u] == 0 and len(Cu2) != 0:
                dico_cliq[noeud_u] = 3;
            elif dico_cliq[noeud_u] == 0 and len(Cu2) == 0:
                dico_cliq[noeud_u] = 1;
            else:
                dico_cliq[noeud_u] = 2;
            # ajout Cu1 et cu2 a C ==> debut
            dico_ver[noeud_u] = 1;
            if len(Cu1) > 1 and len(Cu2) > 1:
                C.append( Cu1 ); C.append( Cu2 );
            elif len(Cu1) > 1 and len(Cu2) < 1:
                C.append( Cu1 );
            elif len(Cu1) < 1 and len(Cu2) < 1:
                print("<=====Cu1={}, Cu2={} BIZARRE =====>".format(Cu1,Cu2))
            # ajout Cu1 et cu2 a C ==> fin
            aretes_noeud_u = aretes_cliques([Cu1, Cu2]) if len(Cu2)>1 else aretes_cliques([Cu1]) 
            for w in gamma_u:
                E_epsi_w = [arete for arete in liste_aretes_Ec \
                            if (arete[0],arete[1]) not in aretes_noeud_u and \
                               (arete[1],arete[0]) not in aretes_noeud_u]
                alpha_w = len([arete for arete in E_epsi_w if arete[0] == w or arete[1] == w]);
                if alpha_w > 0 and dico_cliq[w] == 0:
                    dico_cliq[w] = 3;
                elif alpha_w > 0 and dico_cliq[w] == 3:
                    dico_cliq[w] = -1;
                elif dico_cliq[w] == 0:
                    dico_cliq[w] = 1;
                elif dico_cliq[w] == 3:
                    dico_cliq[w] = 2;
                else:
                    print("<====== Noeud {} BIZARRE ======>".format(w))
            # mise a jour variables
            liste_aretes_Ec = aretes_a_supprimer(liste_aretes_Ec, [list(Cu1), list(Cu2)]);
            dico_gamma_noeud = fct_aux.gamma_noeud(matE, liste_aretes_Ec);
        else:
            dico_cliq[noeud_u] = -1;
            print("<====== Noeud {} -1 ======>".format(noeud_u))
    ### while fin
    return C, dico_cliq, liste_aretes_Ec, ordre_noeuds_traites;
def debug_couverture_new():
#    ---
    test = "EN TEST"; effet_joule = 0; seuil_U = 10; nbre_ts = 10; dimMat = 0; nb_lien = 0;
    chemin_matrices = "/home/willy/topologyLearning/datas/data_test/matrices/"
    chemin_datasets = "/home/willy/topologyLearning/datas/data_test/datasets/"
    matE, matA, dico_dual_arc_sommet = simu50.matriceE(dimMat, nb_lien, chemin_matrices, \
                                                       chemin_datasets, nbre_ts, epsilon, \
                                                       effet_joule, test)
    print(" dico_dual_arc_sommet= ",dico_dual_arc_sommet)
    
    simulation = True; 
    liste_grandeurs = fct_aux.liste_grandeurs(chemin_datasets)
    df_fusion = VerifCorrel.merger_dataframe(liste_grandeurs, chemin_datasets) 
    df_fusion.fillna(0, inplace = True);
    arguments_MAJ = {"dico_sommet_arete": dico_dual_arc_sommet, "df_fusion": df_fusion, \
                     "seuil_U": seuil_U, "epsilon":epsilon, "chemin_dataset": chemin_datasets,\
                     "simulation": simulation, "grandeurs": liste_grandeurs}
    dico_cliq = dict(); dico_ver = dict(); C = list(); noeuds_traites = list();
    for sommet in matE.columns.tolist():   # nbre de noeuds dans le graphe
        dico_cliq[sommet] = 0; dico_ver[sommet] = 0
    
    liste_aretes_Ec = fct_aux.liste_arcs(matE)
    dico_gamma_noeud = fct_aux.gamma_noeud(matE, liste_aretes_Ec) # {"2":[3,{"1","3","4"}],....}
    
    C, dico_cliq, liste_aretes_Ec, noeuds_traites = \
    cover_new_version(dico_cliq, dico_gamma_noeud, liste_aretes_Ec.copy(), \
                      matE, dico_ver, arguments_MAJ)
    
    print("len={}, C={}".format(len(C),C))
    pass
#### algo couverture new version ####### ===> fin
if __name__ == '__main__':
    
    start= time.time()
    
    nbre_ts = 10; effet_joule = 0.1
    epsilon = 0.75; nb_lien = (2,5)
    dimMat = 5; test = "EN TEST"
    ascendant_1 = True; 
    chemin_matrices  = "data/matrices/"; 
    chemin_dataset = "data/datasets/";
    nbre_aretes_a_supprimer = 2;
    number_items_pi1_pi2 = 1;
#    test_decouverte_clique(dimMat, nb_lien, epsilon, chemin_matrices, \
#                           chemin_dataset, ascendant_1, number_items_pi1_pi2)

#    test_decouverte_clique_corriger_noeuds(dimMat, nb_lien, epsilon, chemin_matrices, \
#                                           chemin_dataset, ascendant_1, \
#                                           nbre_aretes_a_supprimer, number_items_pi1_pi2)

#    test_decouvClique_dataReal()
    debug_couverture_new();
    print ("running time: ",time.time() - start)
    