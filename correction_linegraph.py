#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 11:31:09 2017

@author: willy
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 14:50:58 2017

@author: willy
"""
import pandas as pd
import numpy as np
import clique_max as clique
import fonctions_auxiliaires as fct_aux
import itertools as it
import time
import math;
import graphe_particulier as graph_part
#import tmp_decouverteClique_1006 as decouvClique
import decouverte_cliques as decouvClique
import logging;



def C_z(noeud_z, ens_C_i):
    """
    recherche toutes les cliques contenant le noeud "noeud_z"
    c_z: liste des cliques contenant noeud_z et dont 
            tous les noeuds ne contiennent que noeud_z pris 2 a 2
    ens_C_i: ensemble de cliques obtenu a l'etat i-1 et nous sommes a l'etat i
    """
#    #print("ens_C_i ", ens_C_i)
    node_z_on_cliques = [C for C in ens_C_i if noeud_z in C and len(C) > 2]
    c_z = list()
    #print("node_z_on_cliques= ", node_z_on_cliques)
    while len(node_z_on_cliques):
        C1 = node_z_on_cliques.pop()
        booleen = True
        for C in node_z_on_cliques:
            if C1.intersection(C) != {noeud_z}:
                #print("C1=",C1," C=",C)
                booleen = False; break;
        if booleen == True:
            c_z.append(C1)
    return c_z
    pass

def aretes_ens_C_i(ens_C_i):
    """
    liste des aretes de chaque clique contenue dans ens_C_i
    aretes = list des aretes 
    NB: on a que les aretes dans une direction (1,2)
        pour les tests, il faut verifier l'arete (1,2) et (2,1)
    """
    aretes = list()
    for cliq in ens_C_i:
        aretes.extend( list(it.combinations(cliq,2)) )
    return aretes
    pass

def cliques_contractables(noeud_z, C_z, ens_C_i, liste_aretes_E_C_i, liste_aretes_ens_C_i):
    """
    retourne la liste de couples de cliques contractables
     ======> format de return [({C1,C2},[C1, C2]), ({C3,C4},[C3, C4]),....]
     
    liste_aretes_E_C_i: liste des aretes ds le graphe a l'etape i cad a la fin de l'etape i+1 
    """
    l_cliqs_contractables = list()
#    #print("C_z = ",C_z)
    #for tuple_cliq in it.combinations(ens_C_i, 2): =====> A Effacer
    for tuple_cliq in it.combinations(C_z, 2):
        #tu0 = tuple_cliq[0]; tu1 = tuple_cliq[1];
        booleen = True
        for arete in it.product( tuple_cliq[0] - {noeud_z} , tuple_cliq[1] - {noeud_z} ):
            if ( (arete[0],arete[1]) in liste_aretes_E_C_i or (arete[1],arete[0]) in liste_aretes_E_C_i) \
                and ( (arete[0],arete[1]) in liste_aretes_ens_C_i or (arete[1],arete[0]) in liste_aretes_ens_C_i )\
                and set(arete) not in ens_C_i:
                booleen = False; #print("**** not contractable: ", tuple_cliq, " arete_problem: ", arete)
                break;
        if booleen == True:
            l_cliqs_contractables.append( (tuple_cliq[0].union(tuple_cliq[1]), [tuple_cliq[0],tuple_cliq[1]] ))
    
    #ajout des cliques contractables dont unee clique est l'ensemble vide
    for tuple_cliq in it.combinations(C_z, 1):
        l_cliqs_contractables.append( (tuple_cliq[0], [tuple_cliq[0]]) )
        
    return l_cliqs_contractables
    pass

def S_z(noeud_z, gamma_z, ens_C_i, liste_aretes_ens_C_i):
    """
    gamma_z tel que 
        - {noeud_z, voisin_z} in ens_C_i ou
        - (noeud_z, voisin_z) not in liste_aretes_ens_C_i and (voisin_z, noeud_z) not in liste_aretes_ens_C_i
    return la liste des voisin_z
    """
    s_z = list()
    for voisin_z in gamma_z:
        if {noeud_z, voisin_z} in ens_C_i or \
            ((noeud_z, voisin_z) not in liste_aretes_ens_C_i and \
             (voisin_z, noeud_z) not in liste_aretes_ens_C_i):
            s_z.append(voisin_z)
    return s_z

def voisine_z(noeud_z, C, C_z, S_z, n = 2):
    """
    clique C n'appartenant pas a C_z et card(C \cup S_Z) > n
    return un booleen tel que 
        si False alors la clique C n'est pas voisine de noeud_z
        si True alors la clique C est voisine de noeud_z
    """
    if C in C_z:
        return False
    #if len( C.intersection(S_z) ) > n-1:   # == >= 1
    if len( C.intersection(S_z) ) > n-2:   # == >= 1
        return True

def D_z_c( C, gamma_z, C_z):
    """
    determine si une clique c_z contenant "noeud_z" est dependance avec C
            si oui retourne les listes des cliques dependante avec C
    return
        list_clique: liste de cliques dependantes de C 
        
    C_z: liste des cliques contenant noeud_z = z 
    gamma_z: un ens de voisins de noeud_z = z cad il est un set
    
    NB : C est une clique voisine de noeud_z
    """
    list_clique = list()
    list_clique = [cliq_z for cliq_z in C_z if len(C.intersection(cliq_z.intersection(gamma_z))) != 0]
    return list_clique
def ens_cliques_contractables( C, liste_cliques, noeud_z, C_z, ens_C_i, liste_aretes_ens_C_i):
    """
    verifie si une paire de cliques de liste_cliques est contractable.
    Si cette paire est contractable alors
        * on ajoute [ paire[0], paire[1], C] a l'ensemble des listes augmentantes

    """
    if len(liste_cliques) == 0:
        return [(C, set())]
    if len(liste_cliques) == 1:
#        #print("-------> liste_cliques = ", len(liste_cliques) )
        return [( liste_cliques.pop(), C)]
    else:
        res_liste_cliques = list();
        for tu_contr_possib in it.combinations(liste_cliques, 2):
            set_inter = set(tu_contr_possib[0]).intersection(set(tu_contr_possib[1]));
#            #print("+++++++> liste_cliques =", liste_cliques ," tu_contr_possib[0] = ", tu_contr_possib[0], " tu_contr_possib[1] = ", tu_contr_possib[1], " set_inter = ",set_inter)
            set_inter = set_inter - {noeud_z};
#            #print("+++++++> set_inter = ",set_inter)
            tu0 = set(tu_contr_possib[0]) - set_inter ; tu1 = set(tu_contr_possib[1]) - set_inter;
        
#            #print("+++++++> tu0 = ", tu0, " tu1 = ", tu1)
            booleen = True;
            for arete in it.product( tu0 , tu1 ):
                if ( (arete[0],arete[1]) in liste_aretes_ens_C_i \
                    or (arete[1],arete[0]) in liste_aretes_ens_C_i) \
                    and set(arete) not in ens_C_i:
                        booleen = False; break;
            if booleen == True:
                res_liste_cliques.append( (tu_contr_possib[0], tu_contr_possib[1], C) )
                
        return res_liste_cliques   
def augmentation_z(noeud_z, gamma_z, C_z, S_z, ens_C_i, liste_aretes_ens_C_i, n=2):
    """
    format l_augm_z: [ ({1,2,3,z},[{1,2,3}]), ... ]
    """
    if noeud_z == "R08_2":
        print("noeud_z=",noeud_z,";gamma_z=",gamma_z,";C_z=",C_z,";S_z=",S_z,";ens_C_i=",ens_C_i,";liste_aretes_ens_C_i=",liste_aretes_ens_C_i)
#    logging.debug("noeud_z=%s; gamma_z=%s; C_z=%s; S_z=%s; ens_C_i=%s; liste_aretes_ens_C_i= %s",\
#                   noeud_z, gamma_z, C_z, S_z, ens_C_i, liste_aretes_ens_C_i)
    l_augm_z = list()
    for C in ens_C_i:
        liste_cliques = list()
        if voisine_z(noeud_z, C, C_z, S_z, n) == True:
            liste_cliques = D_z_c( C, gamma_z, C_z)
            #print("----> voisin(z) = ", C," DzC = ", liste_cliques);
            liste_cliques = ens_cliques_contractables( C, liste_cliques, \
                                                       noeud_z, C_z, ens_C_i, liste_aretes_ens_C_i)
            while len(liste_cliques) != 0:
                cliq_contratables = liste_cliques.pop();
                l_augm_z.append( (set().union( *cliq_contratables).union({noeud_z}), list(cliq_contratables)) )
#    for tu in l_augm_z:
#        #print ("-----> tu ", tu)
    return l_augm_z

###================= reecriture tuple_4_pi1_pi2_pis ============================
"""
format l_4uplet_pi1_pi2_pis_cliqsASupp:
    (pi1, pi2, pis, cliqsASupp)
format l_pi1_pi2:
    (pi1, pi2, [cliqsASupp])
    combinaison 2 par 2 pour un C1 et S1 donne
format pi1 === je ne sais pas trop
"""
def all_subsets(ss):
    #print("all ss ", ss)
    return it.chain(*map(lambda x: it.combinations(ss, x), range(1, len(ss)+1)))

def combinaison_possible_C1_S1(C1, s1):
    """
    p1 = c1[0]; p2= s1.pop(), p3 = s1.pop(), .....
    resultat p1p2p3, p1p2, p1p3, p2p3
    
    NB: TOUTES LES COMBINAISONS POSSIBLE DE GROUPE DE 1, 2 ,... LEN(N) 
    AVEC N LA TAILLE DE LA LISTE DE [C1, s1] 
    """
#    #print("C1 = ", C1)
#    #print("s1 = ", s1)
    l = list(); l_r = list()
    l.append(C1[0]); l.extend(s1);
    for subset in all_subsets(l):
        if C1[0] in subset:
            l_r.append( (set().union(*subset),C1[1]) )
    return l_r
    
def C1_S1(noeud_z, C1, S_z, ens_C_i):
    """
    sommets v de S_z n'appartenant a aucune clique de C1 tel que:
        v in S_z, x in C1, !C' in ens_C_i dont len(C') >2 et {v,x} convert par C' ===> 1
        
    C1: clique contractable de la forme ({C1,C2}, [C1,C2])
    return liste de sous ensemble verifiant 1
    
    exple: 
        C1 = ({1,2,3},[{1,2,3}]);s1= [{4},{5}];
        ===> return [ ({1, 2, 3}, [{1, 2, 3}]), ({1, 2, 3, 4}, [{1, 2, 3}]),
         ({1, 2, 3, 5}, [{1, 2, 3}]), ({1, 2, 3, 4, 5}, [{1, 2, 3}])]
    """
    #TODO mettre une entree dans fichier log "run fonction S1 ===> debut"
    l_sous_ensemble_C1_S1 = list()
    if len(C1) == 0:
        for i in range(1,len(S_z)+1):
            for tu in it.combinations(S_z, i):
                l_sous_ensemble_C1_S1.append( (set().union( *[{noeud_z},tu] ),[]) )   #[{'z'},{1,2,3}]
    else:
        S1_tmp = list(); s1 = list()
        S1_tmp = [s for s in S_z if s not in C1[0]] 
        ens_C_i_sup2 = [Cliq for Cliq in ens_C_i if len(Cliq) > 2]
        l_aretes_cliques = aretes_ens_C_i(ens_C_i_sup2)
        
        while len(S1_tmp) != 0:
            booleen = True; noeud_S1_tmp = S1_tmp.pop()
            for arete in it.product({noeud_S1_tmp},C1[0]):
                if (arete[0], arete[1]) in l_aretes_cliques or (arete[1], arete[0]) in l_aretes_cliques:
                    booleen = False; 
                    # TODO mettre une entree dans fichier log
                    break;
            if booleen == True:
                s1.append({noeud_S1_tmp})
                # TODO ici une entre dans le log
        
        # ici faire combinaison de C1[0] et S1 A DEPLACER
        l = combinaison_possible_C1_S1(C1, s1);
        l_sous_ensemble_C1_S1.extend(l)
    #TODO mettre une entree dans fichier log "run fonction S1 ===> fin"
    return l_sous_ensemble_C1_S1
        
def set_C1_S1(noeud_z, cliqs_contract, S_z, ens_C_i):
    """
    return tous les C1 \cup S1
    format return 
    [({'E', 'H', 'I'}, [{'E', 'H', 'I'}]), 
     ({'F', 'E'}, []), 
     ({'D', 'E', 'G'}, [{'D', 'E', 'G'}]) ]
    """
    l_C1_S1 = list()
    if len(cliqs_contract) == 0:
        l_C1_S1.extend( C1_S1(noeud_z, (), S_z, ens_C_i) )
    for C1 in cliqs_contract:
        #print("cliq contractables = ", C1)
        sub_set_S1 = C1_S1(noeud_z, C1, S_z, ens_C_i);
        l_C1_S1.extend(sub_set_S1)
#    #print("l_C1_S1 = ", l_C1_S1)
    
    # a ajouter toutes combinaison possibles dans l_C1_S1
    return l_C1_S1

def pi1_pi2(noeud_z, l_augm_z, l_C1_S1, number_items_pi1_pi2):
    """
    format l_C1_S1 : [ ({1, 2, 3}, [{1, 2, 3}]), ({1, 2, 3, 4}, [{1, 2, 3}]), ....]
    format l_augm_z: [ ({1,2,3,z},[{1,2,3}]), ... ]
    * on fusionne les 2 listes l_C1_S1 et l_augm_z, on a fusion_liste_C1S1_augmz
    * ensuite on fait une combinaison par 2 de la liste fusion_liste_C1S1_augmz tel que 
        - dans chaque tuple, l intersection des 2 items du tuple == {noeud_z} 
        
    number_items_pi1_pi2: "pourcentage" d'elements de l_pi1_pi2 a recuperer ( <=1 )
    """
    fusion_liste_C1S1_augmz = l_C1_S1 + l_augm_z;
    #print("noeud traite: ", noeud_z," l_C1_S1: ",len(l_C1_S1)," l_augm_z: ",len(l_augm_z)," fusion_l_C1S1_augm_z: ", len(fusion_liste_C1S1_augmz))
    l_pi1_pi2 = list()
#    for tuple_possible in it.combinations(fusion_liste_C1S1_augmz, 2):
#        if tuple_possible[0][0].intersection(tuple_possible[1][0]) == {noeud_z}:
#            #tmp_pi1_pi2 = (tuple_possible[0][0],  tuple_possible[1][0], tuple_possible[0][1]+tuple_possible[1][1])
#            pi1 = tuple_possible[0][0]; 
#            pi2 = tuple_possible[1][0];
#            l_cont = tuple_possible[0][1]+tuple_possible[1][1];
#            l_pi1_pi2.append( (pi1,pi2, l_cont) )
#    return l_pi1_pi2
    
    # test debut
    if len(fusion_liste_C1S1_augmz) == 1:
        tuple_possible = fusion_liste_C1S1_augmz.pop()
        # tuple_possible = ({'D', 'F', 'A', 'E'}, [{'D', 'F', 'A', 'E'}])
        pi1 = tuple_possible[0]
        l_cont = tuple_possible[1]
        l_pi1_pi2.append( (pi1, set(), l_cont) )
        return l_pi1_pi2
    #for tuple_possible in it.combinations(fusion_liste_C1S1_augmz, 2):
    number_items = math.ceil(len(fusion_liste_C1S1_augmz) * number_items_pi1_pi2)
    for tuple_possible in it.islice(it.combinations(fusion_liste_C1S1_augmz, 2), number_items):
        ##print("tuple_possible= ", tuple_possible)
        if tuple_possible[0][0].intersection(tuple_possible[1][0]) == {noeud_z}:
            #tmp_pi1_pi2 = (tuple_possible[0][0],  tuple_possible[1][0], tuple_possible[0][1]+tuple_possible[1][1])
            pi1 = tuple_possible[0][0]; 
            pi2 = tuple_possible[1][0];
            l_cont = tuple_possible[0][1]+tuple_possible[1][1];
            l_pi1_pi2.append( (pi1, pi2, l_cont) )
        else:
            pi1 = tuple_possible[0][0].intersection(tuple_possible[1][0]);
            pi2 = set().union( tuple_possible[0][0].union( tuple_possible[1][0] ) - pi1).union({noeud_z})
            l_cont = [pi1]
            l_pi1_pi2.append( (pi1, pi2, l_cont) )
    #print("l_pi1_pi2: ", len(l_pi1_pi2))
    return l_pi1_pi2
    # test fin
    
def pis(noeud_z, gamma_z, l_pi1_pi2, ens_C_i):
    """
    retourne le n-uplet (pi1,pi2,pis,[cliques a supprimer])
    l_pi1_pi2: liste des pi1_pi2 de la forme 
                [({"a","b"}, {"c","d"},[{'A','B','C'}]), ... ]
    """
    l_4uplet_pi1_pi2_pis_cliqsASupp = list()
#    #print("===== l_pi1_pi2 = ", l_pi1_pi2)
    for p1_p2 in l_pi1_pi2:
#        #print("p1_p2 = ", p1_p2)
        pi1 = p1_p2[0] # pi1 de la forme pi1 = ({1,2,3},[{1,2},{2,3}])
        pi2 = p1_p2[1]
        gammaZ_pi1 = set(); gammaZ_pi2 = set(); Y = set(); pi_s = set()
        gammaZ_pi1 = set(pi1).intersection(gamma_z)
        gammaZ_pi2 = set(pi2).intersection(gamma_z)
        Y = gammaZ_pi1.union(gammaZ_pi2)
        pi_s = set(gamma_z) - Y
        
        l_cliques_ASupp = list()
        # pour eviter les doublons de cliste a supprimer et aussi les cliques a suppirmer 
        # n'etant pas dans ens_C_i
        l_cliques_ASupp = list()
        for cliq_a_suppr in p1_p2[2]:
            if cliq_a_suppr in ens_C_i and cliq_a_suppr not in l_cliques_ASupp:
                l_cliques_ASupp.append(cliq_a_suppr)
        l_4uplet_pi1_pi2_pis_cliqsASupp.append( (pi1, pi2, pi_s, l_cliques_ASupp) )
        
    ##print("l_4uplet_pi1_pi2_pis_cliqsASupp = ", l_4uplet_pi1_pi2_pis_cliqsASupp)
    return l_4uplet_pi1_pi2_pis_cliqsASupp
###================= reecriture tuple_4_pi1_pi2_pis ============================     

#### definition de la fonction de cout d'une arete ====> debut
def fct_cout(arete, dico_proba_cases, coef_fct_cout, label, correl_seuil):
    """
    le cout d'une arete est fonction du seuil de correlation et 
    de la valeur de correlation.
    En effet le cout est tres eleve (proche de 1) quand la valeur de correlation tend vers 0 ou 1.
    Par compte, la valeur de correlation proche de correl_seuil implique un cout tres faible.
        
    FACTEUR_MULTIPLICATIF = 1;
    EXPOSANT = 1
    coef_fct_cout = (exposant = 1, facteur_multiplicatif = 1)
    """
    FACTEUR_MULTIPLICATIF = coef_fct_cout[1]; EXPOSANT = coef_fct_cout[0]; 
    type_fct_cout = coef_fct_cout[2]
    cout = 0;
    #### A EFFACER
#    if (arete[0], arete[1]) in dico_proba_cases:
#        print("e={} OK_01".format((arete[0], arete[1])))
#    elif (arete[1], arete[0]) in dico_proba_cases:
#        print("e={} OK_10".format((arete[1], arete[0])))
#    else:
#        print("e={} not in dico_proba_cases".format((arete[1], arete[0])))
    #### A EFFACER
    if (arete[0], arete[1]) in dico_proba_cases:
        if type_fct_cout == "cloche":
#            cout = round(abs( pow(pow(dico_proba_cases[(arete[0],arete[1])],1) - correl_seuil, 3) ), 3) 
#            F_c = | 4\cdot((p-s) - (s-0.5))^2 |^{1.5}  
            cout = pow( abs( 4*pow( (dico_proba_cases[(arete[0],arete[1])] - correl_seuil) - (correl_seuil-0.5) ,2) ), 1.5)
        elif type_fct_cout == "lineaire":
            if label == "suppression_arete":
                cout = FACTEUR_MULTIPLICATIF * pow(dico_proba_cases[(arete[0],arete[1])], EXPOSANT);
            elif label == "ajout_arete":
                cout = FACTEUR_MULTIPLICATIF * pow(1 - dico_proba_cases[(arete[0],arete[1])], EXPOSANT);
        elif type_fct_cout == "lineare_iourte_priorite_supp":
            # car suppression coute moins cher.
            if label == "suppression_arete":
#                cout = FACTEUR_MULTIPLICATIF * pow(dico_proba_cases[(arete[0],arete[1])], EXPOSANT);
                cout = 1 * pow(dico_proba_cases[(arete[0],arete[1])], EXPOSANT); 
            elif label == "ajout_arete":
#                cout = 1 * pow(1 - dico_proba_cases[(arete[0],arete[1])], EXPOSANT);
                cout = FACTEUR_MULTIPLICATIF * pow(dico_proba_cases[(arete[0],arete[1])], EXPOSANT);
        elif type_fct_cout == "lineare_iourte_priorite_ajout": 
            # ajout coute moins cher
            if label == "suppression_arete":
#                cout = 1 * pow(dico_proba_cases[(arete[0],arete[1])], EXPOSANT);
                cout = FACTEUR_MULTIPLICATIF * pow(dico_proba_cases[(arete[0],arete[1])], EXPOSANT);
            elif label == "ajout_arete":
#                cout = FACTEUR_MULTIPLICATIF * pow(1 - dico_proba_cases[(arete[0],arete[1])], EXPOSANT);
                cout = 1 * pow(dico_proba_cases[(arete[0],arete[1])], EXPOSANT);
        elif type_fct_cout == "lineaire_simul50Graphes_priorite_ajout": 
            # ajout coute moins cher
            if label == "suppression_arete":
                cout = FACTEUR_MULTIPLICATIF * pow(dico_proba_cases[(arete[0],arete[1])], EXPOSANT);
            elif label == "ajout_arete":
                cout = 1 * pow(dico_proba_cases[(arete[0],arete[1])], EXPOSANT); 
        elif type_fct_cout == "lineaire_simul50Graphes_priorite_supp":
            # car suppression coute moins cher.
            if label == "suppression_arete":
                cout = 1 * pow(dico_proba_cases[(arete[0],arete[1])], EXPOSANT); 
            elif label == "ajout_arete":
                cout = FACTEUR_MULTIPLICATIF * pow(dico_proba_cases[(arete[0],arete[1])], EXPOSANT);
        elif type_fct_cout == "lineaire_simul50Graphes_priorite_aucune":
            # car suppression et ajout coutent 1.
            cout = 1 * pow(dico_proba_cases[(arete[0],arete[1])], EXPOSANT); 
            
    elif (arete[1], arete[0]) in dico_proba_cases:
        if type_fct_cout == "cloche":
            cout = round(abs( pow(pow(dico_proba_cases[(arete[1],arete[0])],1) - correl_seuil, 3) ), 3) 
        elif type_fct_cout == "lineaire":
            if label == "suppression_arete":
                cout = FACTEUR_MULTIPLICATIF * pow(dico_proba_cases[(arete[1],arete[0])], EXPOSANT);
            elif label == "ajout_arete":
                cout = FACTEUR_MULTIPLICATIF * pow(1 - dico_proba_cases[(arete[1],arete[0])], EXPOSANT);
        elif type_fct_cout == "lineare_iourte_priorite_supp":
            if label == "suppression_arete":
                cout = 1 * pow(dico_proba_cases[(arete[1],arete[0])], EXPOSANT);
            elif label == "ajout_arete":
                cout = FACTEUR_MULTIPLICATIF * pow(dico_proba_cases[(arete[1],arete[0])], EXPOSANT);
        elif type_fct_cout == "lineare_iourte_priorite_ajout":
            if label == "suppression_arete":
                cout = FACTEUR_MULTIPLICATIF * pow(dico_proba_cases[(arete[1],arete[0])], EXPOSANT);
            elif label == "ajout_arete":
                cout = 1 * pow(dico_proba_cases[(arete[1],arete[0])], EXPOSANT);
        elif type_fct_cout == "lineaire_simul50Graphes_priorite_ajout": 
            # ajout coute moins cher
            if label == "suppression_arete":
                cout = FACTEUR_MULTIPLICATIF * pow(dico_proba_cases[(arete[1],arete[0])], EXPOSANT);
            elif label == "ajout_arete":
                cout = 1 * pow(dico_proba_cases[(arete[0],arete[1])], EXPOSANT); 
        elif type_fct_cout == "lineaire_simul50Graphes_priorite_supp":
            # car suppression coute moins cher.
            if label == "suppression_arete":
                cout = 1 * pow(dico_proba_cases[(arete[1],arete[0])], EXPOSANT); 
            elif label == "ajout_arete":
                cout = FACTEUR_MULTIPLICATIF * pow(dico_proba_cases[(arete[1],arete[0])], EXPOSANT);
        elif type_fct_cout == "lineaire_simul50Graphes_priorite_aucune":
            # car suppression et ajout coutent 1.
            cout = 1 * pow(dico_proba_cases[(arete[0],arete[1])], EXPOSANT); 
    return cout
#### definition de la fonction de cout d'une arete ====> fin
def cout_pi1_pi2_pis(noeud_z, tuple_pi1_pi2_pis, liste_aretes_E_C_i, \
                     liste_aretes_E0, dico_proba_cases, coef_fct_cout, \
                     correl_seuil):
    """
    but: cherche le cout d'un triplet 
    ps1, ps2, pis sont des ensembles
    liste_aretes_aAjouter: liste des aretes a ajouter a liste_arcs pour obtenir un linegraph
    liste_aretes_aSuppr: liste des aretes a supprimer a liste_arcs pour obtenir un linegraph
    cliques_aSupprDe_ens_C : liste des cliques a supprimer de ens_C
    
    NB: dico_proba_cases: probabilite de chaque case de matE. son prototype est 
        dico_proba_cases[(a,b)] = 0..1
    """
    liste_aretes_aAjouter = list(); liste_aretes_aSuppr = list(); 
    pi1 = tuple_pi1_pi2_pis[0]; pi2 = tuple_pi1_pi2_pis[1]
    pis = tuple_pi1_pi2_pis[2]; cliques_aSupprDe_ens_C_i = tuple_pi1_pi2_pis[3];
    
    #pi1 = aretes a ajouter
    cout_pi1  = 0 
    liste_aretes_pi1 = aretes_ens_C_i([pi1]);
#    print("YYY 1")
    for arete in liste_aretes_pi1:
        if (arete[0],arete[1]) not in liste_aretes_E_C_i and \
            (arete[1],arete[0]) not in liste_aretes_E_C_i:
            liste_aretes_aAjouter.append( (arete[0],arete[1]) )
            if (arete[0],arete[1]) not in liste_aretes_E0 or \
                (arete[1],arete[0]) not in liste_aretes_E0:
#                    cout_pi1 += 1; 
                    label_ajout = "ajout_arete";
#                    print("ZZZ 11")
                    cout_pi1 += fct_cout(arete, dico_proba_cases, coef_fct_cout,\
                                         label_ajout, correl_seuil)
#                    print("ZZZ 12")
    #pi2 = aretes a ajouter
    cout_pi2  = 0 
    liste_aretes_pi2 = aretes_ens_C_i([pi2]);
#    print("YYY 2")
    for arete in liste_aretes_pi2:
        if (arete[0],arete[1]) not in liste_aretes_E_C_i and \
            (arete[1],arete[0]) not in liste_aretes_E_C_i:
            liste_aretes_aAjouter.append( (arete[0],arete[1]) )
            if (arete[0],arete[1]) not in liste_aretes_E0 or \
                (arete[1],arete[0]) not in liste_aretes_E0:
#                    cout_pi2 += 1
                    label_ajout = "ajout_arete";
#                    print("ZZZ 21")
                    cout_pi2 += fct_cout(arete, dico_proba_cases, coef_fct_cout,\
                                         label_ajout, correl_seuil)
#                    print("ZZZ 22")
    #pis = aretes a supprimer
    cout_pis = 0
#    print("YYY 3")
    for noeud_pis in pis:
        if (noeud_z, noeud_pis) in liste_aretes_E_C_i or (noeud_pis, noeud_z) in liste_aretes_E_C_i:
            liste_aretes_aSuppr.append( (noeud_z, noeud_pis) )
#            cout_pis += 1;
            label_supp = "suppression_arete";
#            print("ZZZ 31")
            cout_pis += fct_cout((noeud_z, noeud_pis), dico_proba_cases, coef_fct_cout, \
                                 label_supp, correl_seuil)
#            print("ZZZ 32")
            
    cout_noeud_z = 0
    cout_noeud_z = cout_pi1 + cout_pi2 + cout_pis
    #cout_noeud_z = cout_pi1 + cout_pi2 - cout_pis
    return cout_noeud_z, liste_aretes_aAjouter, liste_aretes_aSuppr, cliques_aSupprDe_ens_C_i 

#def critere_selection_min_debug(dico_cout_5uplets):
#    """
#    le critere de selection du minimum est le maximum de la difference entre les aretes_aAjouter et les aretes_aSupp
#        max_diff = max( len(aretes_aAjouter) - len(aretes_aSupp) )
#    ON PRIVILEGIE CELUI QUI  FAIT LE MOINS DE MODIFICATION POSSIBLE
#    ON PRIVILEGE L'AJOUT D'ARETES par rapport a la SUPPRESSION 
#    A REECRIRE ===> NON car BON
#    """
#    min_liste_arete_aAjouter = list(); min_liste_aretes_aSuppr = list()
#    min_cliques_aSupprDe_ens_C = list(); min_tuple_solution = None; min_cout = 0
#    
#    min_cout = min(dico_cout_5uplets.keys()) if len(dico_cout_5uplets) != 0 else None
#    
##    #print("critere_selection_min===dico_cout_5uplets : ", dico_cout_5uplets, "critere_selection_min===len(dico_cout_5uplets) : ", len(dico_cout_5uplets), "critere_selection_min=== min_cout= ",min_cout)
#    if min_cout == None :
#        return min_cout, min_tuple_solution, min_liste_arete_aAjouter, min_liste_aretes_aSuppr, min_cliques_aSupprDe_ens_C 
#
#    if len(dico_cout_5uplets[min_cout]) == 1:
##        #print("dico_cout_5uplets[min_cout]: ", dico_cout_5uplets[min_cout])
#        min_tuple_solution = dico_cout_5uplets[min_cout][0][1]
#        min_liste_arete_aAjouter = dico_cout_5uplets[min_cout][0][2]
#        min_liste_aretes_aSuppr = dico_cout_5uplets[min_cout][0][3]
#        min_cliques_aSupprDe_ens_C = dico_cout_5uplets[min_cout][0][4]
#        
#    elif len(dico_cout_5uplets[min_cout]) > 1:
#        max_diff = 0; max_diff = -1000
#        for first_5uplet in dico_cout_5uplets[min_cout]:
#            diff = np.abs(len(first_5uplet[2]) - len(first_5uplet[3]))
##            #print("diff=",diff," max_diff=",max_diff)
#            if max_diff <= diff:
#                max_diff = diff
#                min_tuple_solution = first_5uplet[1]
#                min_liste_arete_aAjouter = first_5uplet[2]
#                min_liste_aretes_aSuppr = first_5uplet[3]
#                min_cliques_aSupprDe_ens_C = first_5uplet[4]
##                #print("max_diff: ", max_diff, " first_5uplet: ",first_5uplet)
#    return min_cout, min_tuple_solution, min_liste_arete_aAjouter, min_liste_aretes_aSuppr, min_cliques_aSupprDe_ens_C 
##### debug critere => selection du type de critere ===> debut
#%%
def critere_selection_min(dico_cout_5uplets, critere):
    """
    si critere == 0:
        ON PRIVILEGIE CELUI QUI  FAIT LE MOINS DE MODIFICATION POSSIBLE : 
        max_diff = max( len(aretes_aAjouter) - len(aretes_aSupp) ) 
    si critere == 1:
        ON PRIVILEGE CELUI QUI  FAIT LE MOINS DE MODIFICATION POSSIBLE et dont L'AJOUT D'ARETES >=  len(aretes_aAjouter)
        max_diff = max( len(aretes_aAjouter) - len(aretes_aSupp) ) and max_ajout <= len(aretes_aAjouter)
    si critere == 2:
        ON PRIVILEGE CELUI QUI  FAIT LE MOINS DE MODIFICATION POSSIBLE et dont La SUPPRESSION D'ARETES < len(aretes_aSupp)
        max_diff = max( len(aretes_aAjouter) - len(aretes_aSupp) ) and max_supp <= len(aretes_aSupp)
    """
    min_liste_arete_aAjouter = list(); min_liste_aretes_aSuppr = list()
    min_cliques_aSupprDe_ens_C = list(); min_tuple_solution = None; min_cout = 0
    
    min_cout = min(dico_cout_5uplets.keys()) if len(dico_cout_5uplets) != 0 else None
    
#    #print("critere_selection_min===dico_cout_5uplets : ", dico_cout_5uplets, "critere_selection_min===len(dico_cout_5uplets) : ", len(dico_cout_5uplets), "critere_selection_min=== min_cout= ",min_cout)
    if min_cout == None :
        return min_cout, min_tuple_solution, min_liste_arete_aAjouter, min_liste_aretes_aSuppr, min_cliques_aSupprDe_ens_C 

    if len(dico_cout_5uplets[min_cout]) == 1:
#        #print("dico_cout_5uplets[min_cout]: ", dico_cout_5uplets[min_cout])
        min_tuple_solution = dico_cout_5uplets[min_cout][0][1]
        min_liste_arete_aAjouter = dico_cout_5uplets[min_cout][0][2]
        min_liste_aretes_aSuppr = dico_cout_5uplets[min_cout][0][3]
        min_cliques_aSupprDe_ens_C = dico_cout_5uplets[min_cout][0][4]
        
    elif len(dico_cout_5uplets[min_cout]) > 1:
        max_diff = 0; max_diff = -1000; max_ajout = -1000; max_supp = -1000;
        for first_5uplet in dico_cout_5uplets[min_cout]:
            diff = np.abs(len(first_5uplet[2]) - len(first_5uplet[3]))
#            #print("diff=",diff," max_diff=",max_diff)
            if critere == 1:
                if max_diff <= diff and max_ajout <= len(first_5uplet[2]):
                    max_diff = diff; max_ajout = len(first_5uplet[2]);
                    min_tuple_solution = first_5uplet[1]
                    min_liste_arete_aAjouter = first_5uplet[2]
                    min_liste_aretes_aSuppr = first_5uplet[3]
                    min_cliques_aSupprDe_ens_C = first_5uplet[4]
            elif critere == 2:
                if max_diff <= diff and max_supp <= len(first_5uplet[3]):
                    max_diff = diff; max_supp = len(first_5uplet[3]);
                    min_tuple_solution = first_5uplet[1]
                    min_liste_arete_aAjouter = first_5uplet[2]
                    min_liste_aretes_aSuppr = first_5uplet[3]
                    min_cliques_aSupprDe_ens_C = first_5uplet[4]
            else :
                # critere == 0
                if max_diff <= diff:
                    max_diff = diff
                    min_tuple_solution = first_5uplet[1]
                    min_liste_arete_aAjouter = first_5uplet[2]
                    min_liste_aretes_aSuppr = first_5uplet[3]
                    min_cliques_aSupprDe_ens_C = first_5uplet[4]
#                #print("max_diff: ", max_diff, " first_5uplet: ",first_5uplet)
    return min_cout, min_tuple_solution, min_liste_arete_aAjouter, min_liste_aretes_aSuppr, min_cliques_aSupprDe_ens_C 
##### debug critere => selection du type de critere ===> fin
#%%
def cout_min(noeud_z, liste_tuple_pi1_pi2_pis, liste_aretes_E_C_i, liste_aretes_E0,\
             dico_proba_cases, coef_fct_cout, correl_seuil, critere_selection_pi1_pi2):
    """
    recherche le minimum des couts selon le poids ou cout des aretes ajoutes ou supprimes.
    le cout d'une arete est une fonction de cout fct_cout
    fct_cout((a,b)) = abs( pow(pow(p_correl,1) - correl_seuil, 3) ) 
    avec p_correl, la correlation entre a et b
    """
    min_cout = pow(10,9);
#    print("XXX 1")
    dico_cout_5uplets = dict()
    for tuple_pi1_pi2_pis in liste_tuple_pi1_pi2_pis:
#        print("XXX 2")
        min_cout_pi1_pi2_pis, liste_aretes_aAjouter, \
        liste_aretes_aSuppr, cliques_aSupprDe_ens_C = \
        cout_pi1_pi2_pis(noeud_z, tuple_pi1_pi2_pis, liste_aretes_E_C_i, \
                         liste_aretes_E0, dico_proba_cases, coef_fct_cout, \
                         correl_seuil)
#        print("XXX min_cout_pi1_pi2_pis = ", min_cout_pi1_pi2_pis, " tuple_pi1_pi2_pis = ",tuple_pi1_pi2_pis)
        pi1 = tuple_pi1_pi2_pis[0]; pi2 = tuple_pi1_pi2_pis[1]; pis = tuple_pi1_pi2_pis[2]
        if min_cout_pi1_pi2_pis in dico_cout_5uplets.keys():
            dico_cout_5uplets[min_cout_pi1_pi2_pis].append( (min_cout_pi1_pi2_pis, \
                                                             (pi1,pi2,pis),\
                                                             liste_aretes_aAjouter,\
                                                             liste_aretes_aSuppr, \
                                                             cliques_aSupprDe_ens_C ) \
                                                          )
        else:
            dico_cout_5uplets[min_cout_pi1_pi2_pis] =  [(min_cout_pi1_pi2_pis, \
                                                         (pi1,pi2,pis),\
                                                         liste_aretes_aAjouter,\
                                                         liste_aretes_aSuppr, \
                                                         cliques_aSupprDe_ens_C ) \
                                                        ]
#    print ("XXX dico_cout_5uplets = ",dico_cout_5uplets)                                                   
    min_liste_arete_aAjouter = list(); min_liste_aretes_aSuppr = list()
    min_cliques_aSupprDe_ens_C = list(); min_tuple_solution = None
    min_cout, min_tuple_solution, min_liste_arete_aAjouter, min_liste_aretes_aSuppr, min_cliques_aSupprDe_ens_C = \
    critere_selection_min(dico_cout_5uplets, critere_selection_pi1_pi2)
    
    return min_cout, min_tuple_solution, min_liste_arete_aAjouter, min_liste_aretes_aSuppr, min_cliques_aSupprDe_ens_C

def supprimer_cliq_couvert_min_solution(min_tuple_solution_0, ens_C):
    """
    A REVOIR var 
    {1,2,3}.issubset({1,2}) ====> FALSE
    """
    tmp_ens_C = ens_C.copy()
    for cliq in tmp_ens_C:
        if cliq.issubset(min_tuple_solution_0):
            ens_C.remove(cliq)
    return ens_C  

def correction_noeuds(liste_noeuds_1, ens_C_i, liste_aretes_E_C_i, dico_cliq, \
                      ordre_noeuds_traites, number_items_pi1_pi2, \
                      dico_proba_cases, coef_fct_cout, correl_seuil, \
                      critere_selection_pi1_pi2):
    """
    but: corriger tous les noeuds labelises a -1 
    
    liste_noeuds_1 :  liste des noeuds labelises a -1
    ens_C: ensemble de cliques obtenu par l'algo de decouverte de cliques(devrant couvrir tout le graphe)
    liste_aretes: liste des aretes du linegraph initial
    dico_cliq: dictionnaire des labels de chaque noeud du graphe,
    dico_proba_cases: probabilite de chaque case de matE,
    correl_seuil: seuil de correlation pour lequel les correlations sont transforme en 0 ou 1
    """
    ##print("liste_noeuds_1 = ",liste_noeuds_1)
    """
    liste_noeuds_1 =  ['A','D','F', 'G', 'K', 'J', 'L', 'B', 'C', 'H', 'I', 'E']
    avec cet ordre on obtient somme_cout = 10 
    """
    logging.basicConfig(filename='DEBUG_corrlineGraph.log',\
                        filemode='w',level=logging.DEBUG, format='%(asctime)s %(message)s' )
    som_cout_min = 0; noeuds_traites = list(); min_cliques_aSupprDe_ens_C = list()
    
    liste_aretes_E0 = liste_aretes_E_C_i.copy();
    ens_C_i_old = ens_C_i.copy();
#    print("liste_noeuds_1 = ", liste_noeuds_1, " E0:",len(liste_aretes_E_C_i))
    for noeud_z in liste_noeuds_1:
        #print("###### =====> noeud_z = ", noeud_z)
        #print("taille liste_aretes_E_C_i = ", len(liste_aretes_E_C_i))
        
        noeuds_traites.append(noeud_z);

        gamma_z = fct_aux.voisins(liste_aretes_E_C_i, noeud_z)
        #c_z = C_z( noeud_z, liste_aretes_E_C_i)
        c_z = C_z( noeud_z, ens_C_i)
        liste_aretes_ens_C_i = aretes_ens_C_i(ens_C_i)
        l_cliq_contr = cliques_contractables(noeud_z, c_z, ens_C_i, liste_aretes_E_C_i, liste_aretes_ens_C_i)
        s_z = S_z(noeud_z, gamma_z, ens_C_i, liste_aretes_ens_C_i)
        
        l_augm_z = augmentation_z(noeud_z, gamma_z, c_z, s_z, ens_C_i, liste_aretes_E_C_i, 2)
        l_C1_S1 = set_C1_S1(noeud_z, l_cliq_contr, s_z, ens_C_i)
        l_pi1_pi2 = pi1_pi2(noeud_z, l_augm_z, l_C1_S1,number_items_pi1_pi2)
        l_4uplet_pi1_pi2_pis_cliqsASupp = pis(noeud_z, gamma_z, l_pi1_pi2, ens_C_i)
        
        #print("correction_noeuds gamma_z : ", gamma_z)
        if noeud_z == "R08_2":
            print("----- noeud_z:",noeud_z,"---DEBUT---")
            print("gamma_z: ", gamma_z)
            print("c_z : ", c_z)
            print("s_z : ", s_z)
            print("l_cliq_contr : ", l_cliq_contr)
            #print("list_voisine_z : ", list_voisin_z)
            print("list_augm_z: ",l_augm_z)
            print("----- noeud_z:",noeud_z,"---FIN---")
#        logging.debug("----- noeud_z: %s ---DEBUT---",noeud_z)
#        logging.debug("gamma_z: %s ",gamma_z)
#        logging.debug("c_z: %s ",c_z)
#        logging.debug("s_z: %s ",s_z)
#        logging.debug("l_cliq_contr: %s ",l_cliq_contr)
#        logging.debug("l_pi1_pi2: %s ",l_pi1_pi2)
        
#        print("100");
        min_tuple_solution = None; min_cout = 0;
        min_cout, min_tuple_solution, min_liste_arete_Aajouter, \
        min_liste_aretes_aSuppr, min_cliques_aSupprDe_ens_C = \
        cout_min(noeud_z, l_4uplet_pi1_pi2_pis_cliqsASupp, liste_aretes_E_C_i, \
                 liste_aretes_E0, dico_proba_cases, coef_fct_cout, correl_seuil,\
                 critere_selection_pi1_pi2)
#        print("101");
        
        som_cout_min += min_cout if min_cout != None else 0;
        
        if min_tuple_solution != None:
            # supprimer aretes dans liste_arcs_E_C_i
#            #print("aretes_a_Supprimer = ", min_liste_aretes_aSuppr)
#            logging.debug("min_tuple_solution: %s ", min_liste_aretes_aSuppr)
#            print("1011");
            for arete in min_liste_aretes_aSuppr:
                if (arete[0],arete[1]) in liste_aretes_E_C_i:
                    liste_aretes_E_C_i.remove( (arete[0],arete[1]) )
                    ##print("ICI ", (arete[0],arete[1]))
                    for Clik in ens_C_i:
                        if arete[0] in Clik and arete[1] in Clik:
#                            #print("1 cliques_a_Supprimer = ", Clik)
                            ens_C_i.remove(Clik)
                if (arete[1],arete[0]) in liste_aretes_E_C_i:
                    liste_aretes_E_C_i.remove( (arete[1],arete[0]) )
                    ##print("ICI ", (arete[1],arete[0]))
                    for Clik in ens_C_i:
                        if arete[0] in Clik and arete[1] in Clik:
#                            #print("2 cliques_a_Supprimer = ", Clik)
                            ens_C_i.remove(Clik)
            
#            logging.debug("min_tuple_solution: arete_Aajouter: len: %s %s ", len(min_liste_arete_Aajouter),min_liste_arete_Aajouter)
#            logging.debug("min_tuple_solution: AVANT liste_aretes_E_C_i: %s ",len(liste_aretes_E_C_i))
            liste_aretes_E_C_i.extend(min_liste_arete_Aajouter)
#            logging.debug("min_tuple_solution: APRES liste_aretes_E_C_i: %s ",len(liste_aretes_E_C_i))
            #print("aretes_a_ajouter = ", min_liste_arete_Aajouter)
            dico_cliq[noeud_z] = 1
#            print("1012");
        
            # supprimer les cliques de min_cliques_aSupprDe_ens_C dans ens_C ===> Ne pas oublier
#            #print("Avant ens_C_i: ", ens_C_i)
            for clik_a_supprimer in min_cliques_aSupprDe_ens_C :
                if  clik_a_supprimer in ens_C_i:
                    ens_C_i.remove(clik_a_supprimer)
                    #print("3 cliques_a_Supprimer = ", clik_a_supprimer)
#            #print("Apres ens_C_i: ", ens_C_i)
#            print("1013");
            
            #ajouter nouvelle clique dans ens_c
            if len(min_tuple_solution[0]) != 0:
#                print("10131");
#                logging.debug("min_tuple_solution[0]:  %s ", min_tuple_solution[0])
#                logging.debug("min_tuple_solution[0]: ens_C_i: len %s ", len(ens_C_i))
                ens_C_i = supprimer_cliq_couvert_min_solution(min_tuple_solution[0], ens_C_i)
                ens_C_i.append(min_tuple_solution[0])
#                logging.debug("min_tuple_solution[0]: ens_C_i: len %s ", len(ens_C_i))
            if len(min_tuple_solution[1]) != 0:
#                print("10132");
#                logging.debug("min_tuple_solution[1]:  %s ", min_tuple_solution[1])
#                logging.debug("min_tuple_solution[1]: ens_C_i: len %s ", len(ens_C_i))
                ens_C_i = supprimer_cliq_couvert_min_solution(min_tuple_solution[1], ens_C_i)
                ens_C_i.append(min_tuple_solution[1])
#                logging.debug("min_tuple_solution[1]: ens_C_i: len %s ", len(ens_C_i))
            #print("Ajout ens_C_i: ", ens_C_i)
            #print("******** min_cout : ", min_cout)
            #print("******** min_tuple : ", min_tuple_solution)
#        logging.debug("----- noeud_z: %s ---FIN---",noeud_z)
            
    #print("som_cout_min: ",som_cout_min)       
    return ens_C_i, dico_cliq, liste_aretes_E_C_i, min_cliques_aSupprDe_ens_C, \
            noeuds_traites, ordre_noeuds_traites, som_cout_min, ens_C_i_old;  
            
            
############### glouton  debut #############################################################################
def calcul_cout_correction(noeud_z, ens_C_i, liste_aretes_E_C_i, liste_aretes_E0,\
                           dico_cliq, number_items_pi1_pi2, dico_proba_cases, \
                           coef_fct_cout, critere_selection_pi1_pi2):
    """
    but: calcul le cout de la correction d'un noeud
    return  som_cout = min_cout: cout correction du noeud, 
            min_tuple_solution: les 2 cliques qui correspondent a la correction du noeud "noeud_z", contiennent aussi le noeud "noeud_z", 
            min_liste_arete_Aajouter: liste des aretes a ajouter pdt la correction de noeud "noeud_z", 
            min_liste_aretes_aSuppr: liste des aretes a supprimer pdt la correction de noeud "noeud_z", 
            min_cliques_aSupprDe_ens_C: liste des cliques a supprimer de l'ensemble de cliques "ens_C_i" pdt la correction de noeud "noeud_z"
    """
    gamma_z = fct_aux.voisins(liste_aretes_E_C_i, noeud_z)
    c_z = C_z( noeud_z, ens_C_i)
    liste_aretes_ens_C_i = aretes_ens_C_i(ens_C_i)
    l_cliq_contr = cliques_contractables(noeud_z, c_z, ens_C_i, liste_aretes_E_C_i, liste_aretes_ens_C_i)
    s_z = S_z(noeud_z, gamma_z, ens_C_i, liste_aretes_ens_C_i)
        
    l_augm_z = augmentation_z(noeud_z, gamma_z, c_z, s_z, ens_C_i, liste_aretes_E_C_i, 2)
    l_C1_S1 = set_C1_S1(noeud_z, l_cliq_contr, s_z, ens_C_i)
    l_pi1_pi2 = pi1_pi2(noeud_z, l_augm_z, l_C1_S1, number_items_pi1_pi2)
    l_4uplet_pi1_pi2_pis_cliqsASupp = pis(noeud_z, gamma_z, l_pi1_pi2, ens_C_i)
        
#    #print("gamma_z : ", gamma_z," c_z : ", c_z," s_z : ", s_z," l_cliq_contr : ",l_cliq_contr," list_augm_z: ",l_augm_z)

    min_tuple_solution = None; min_cout = 0;
    min_cout, min_tuple_solution, min_liste_arete_Aajouter, min_liste_aretes_aSuppr, min_cliques_aSupprDe_ens_C = \
    cout_min(noeud_z, l_4uplet_pi1_pi2_pis_cliqsASupp, liste_aretes_E_C_i, \
             liste_aretes_E0, dico_proba_cases, coef_fct_cout, critere_selection_pi1_pi2)
    
    dico_correction_noeud_z = dict()
    dico_correction_noeud_z["cout"] = min_cout;
    dico_correction_noeud_z["tuple_cliques_pi1_pi2"] = min_tuple_solution;
    dico_correction_noeud_z["aretes_a_ajouter"] = min_liste_arete_Aajouter;
    dico_correction_noeud_z["aretes_a_supprimer"] = min_liste_aretes_aSuppr;
    dico_correction_noeud_z["cliques_a_supprimer_ens_C"] = min_cliques_aSupprDe_ens_C; 
    
    return dico_correction_noeud_z

def choix_noeud_a_corriger(dico_noeuds_cout_degre, critere_selection_noeuds_1):
    """
    but: choisir un noeud a corriger selon le critere de selection 
         le critere de selection est soit:
             * le cout minimum  (de degre min)
             * le degre minimum (de cout min)
             * les 2: son cout minimum et son degree
         on choisit le noeud ayant le cout le plus petit et 
         lorsque 2 noeuds ont les memes couts, on prend celui ayant le plus grand degree
         NB: 
             si cout_y < cout_z ==> on choisit z
             si cout_y == cout_z et degree_y < degree_z ==> on choisit z
             si cout_y == cout_z et degree_y > degree_z ==> on choisit y
             si cout_y == cout_z et degree_y == degree_z ==> on choisit y
    caracteristiques = [cout_noeud, degree, dico_correction_noeud_z] avec 
        dico_correction_noeud_z = {cout:, tuple_cliques_pi1_pi2:, aretes_a_ajouter:, aretes_a_supprimer:, cliques_a_supprimer_ens_C:}

    return noeud_z
    """
    noeud_choisi = None; noeud_cout = 1000; noeud_degree = 1000;
    #print("dico_noeuds_cout_degre = ")
    if critere_selection_noeuds_1 == "cout_minimum":
        for noeud, caracteristiques in dico_noeuds_cout_degre.items():
            #print("noeud= ",noeud," cout:", caracteristiques[0], " degre:", caracteristiques[1] )
            if noeud_cout > caracteristiques[0]: # caracteristiques[0] = caracteristiques[cout_noeud]
                noeud_choisi = noeud;
                noeud_cout = caracteristiques[0];
                noeud_degree = caracteristiques[1]; # caracteristiques[1] = caracteristiques[degree]
            elif noeud_cout == caracteristiques[0] and noeud_degree < caracteristiques[1]:
                noeud_choisi = noeud;
                noeud_cout = caracteristiques[0];
                noeud_degree = caracteristiques[1];
    elif critere_selection_noeuds_1 == "degre_minimum":
        for noeud, caracteristiques in dico_noeuds_cout_degre.items():
            if noeud_degree > caracteristiques[1]: # caracteristiques[1] = caracteristiques[degree]
                noeud_choisi = noeud;
                noeud_cout = caracteristiques[0];
                noeud_degree = caracteristiques[1]; # caracteristiques[0] = caracteristiques[cout]
            elif noeud_degree == caracteristiques[1] and noeud_cout < caracteristiques[0]:
                noeud_choisi = noeud;
                noeud_cout = caracteristiques[0];
                noeud_degree = caracteristiques[1];
    elif critere_selection_noeuds_1 == "aleatoire":
        # son cout minimum et son degree (identique  cout_min)
        for noeud, caracteristiques in dico_noeuds_cout_degre.items():
            #print("noeud= ",noeud," cout:", caracteristiques[0], " degre:", caracteristiques[1] )
            if noeud_cout > caracteristiques[0]: # caracteristiques[0] = caracteristiques[cout_noeud]
                noeud_choisi = noeud;
                noeud_cout = caracteristiques[0];
                noeud_degree = caracteristiques[1]; # caracteristiques[1] = caracteristiques[degree]
            elif noeud_cout == caracteristiques[0] and noeud_degree < caracteristiques[1]:
                noeud_choisi = noeud;
                noeud_cout = caracteristiques[0];
                noeud_degree = caracteristiques[1];
    else:
        # son cout minimum et son degree (identique  cout_min)
        for noeud, caracteristiques in dico_noeuds_cout_degre.items():
            #print("noeud= ",noeud," cout:", caracteristiques[0], " degre:", caracteristiques[1] )
            if noeud_cout > caracteristiques[0]: # caracteristiques[0] = caracteristiques[cout_noeud]
                noeud_choisi = noeud;
                noeud_cout = caracteristiques[0];
                noeud_degree = caracteristiques[1]; # caracteristiques[1] = caracteristiques[degree]
            elif noeud_cout == caracteristiques[0] and noeud_degree < caracteristiques[1]:
                noeud_choisi = noeud;
                noeud_cout = caracteristiques[0];
                noeud_degree = caracteristiques[1];
    
    #print("noeud_choisi= ", noeud_choisi," noeud_cout= ",noeud_cout," noeud_degree= ", noeud_degree)
    return noeud_choisi        

def noeuds_a_corriger(dico_noeuds_cout_degre):
    """
    but: retourne une liste de noeuds ayant les memes couts et le meme degre 
    caracteristiques = [cout_noeud, degree, dico_correction_noeud_z] avec 
    dico_correction_noeud_z = {cout:, tuple_cliques_pi1_pi2:, aretes_a_ajouter:,\
                               aretes_a_supprimer:, cliques_a_supprimer_ens_C:}
    """
    noeuds_choisis = []; noeud_cout = 1000; noeud_degree = 1000;
    #print("noeudsACorriger dico_noeuds_cout_degre = ")
    for noeud, caracteristiques in dico_noeuds_cout_degre.items():
        #print("noeud= ",noeud," cout:", caracteristiques[0], " degre:", caracteristiques[1] )
        if noeud_cout > caracteristiques[0]:
            noeud_cout = caracteristiques[0];
            noeud_degree = caracteristiques[1];
    for noeud, caracteristiques in dico_noeuds_cout_degre.items():
        if noeud_cout == caracteristiques[0] and noeud_degree == caracteristiques[1]:
            noeuds_choisis.append(noeud)
    
    return noeuds_choisis
    pass    

def appliquer_correction( noeud_z, dico_noeuds_cout_degre, dico_graphe):
    """
    but : appliquer les modifs trouves au noeud "noeud_z" sur tout le graphe.
    NB: 
    dico_correction_noeud_z = {cout:, tuple_cliques_pi1_pi2:, aretes_a_ajouter:, aretes_a_supprimer:, cliques_a_supprimer_ens_C:}
    dico_graphe = {"ens_C_i": , "aretes_E_C_i": aretes_E_C_i, "aretes_E0": aretes_E_0}
    """
    dico_correction_noeud_z = dico_noeuds_cout_degre[noeud_z][2]
    if dico_correction_noeud_z["tuple_cliques_pi1_pi2"] == None:
        return dico_graphe["ens_C_i"], [], dico_graphe["aretes_E_C_i"], dico_graphe["dico_cliq"]
    
    # supprimer aretes dans liste_arcs_E_C_i
    for arete in dico_correction_noeud_z["aretes_a_supprimer"]:
        if (arete[0],arete[1]) in dico_graphe["aretes_E_C_i"]:
            dico_graphe["aretes_E_C_i"].remove( (arete[0],arete[1]) )
            for Clik in dico_graphe["ens_C_i"]:
                if arete[0] in Clik and arete[1] in Clik:
                    dico_graphe["ens_C_i"].remove(Clik)
        if (arete[1],arete[0]) in dico_graphe["aretes_E_C_i"]:
            dico_graphe["aretes_E_C_i"].remove( (arete[1],arete[0]) )
            for Clik in dico_graphe["ens_C_i"]:
                if arete[0] in Clik and arete[1] in Clik:
                    dico_graphe["ens_C_i"].remove(Clik)
    
    dico_graphe["aretes_E_C_i"].extend(dico_correction_noeud_z["aretes_a_ajouter"])
    dico_graphe["dico_cliq"][noeud_z] = 1

    # supprimer les cliques de min_cliques_aSupprDe_ens_C dans ens_C ===> Ne pas oublier
    for clik_a_supprimer in dico_correction_noeud_z["cliques_a_supprimer_ens_C"]:
        if  clik_a_supprimer in dico_graphe["ens_C_i"]:
            dico_graphe["ens_C_i"].remove(clik_a_supprimer)
            #print("3 cliques_a_Supprimer = ", clik_a_supprimer)
            
    #ajouter nouvelle clique dans ens_c
    if len(dico_correction_noeud_z["tuple_cliques_pi1_pi2"][0]) != 0:
        dico_graphe["ens_C_i"] = supprimer_cliq_couvert_min_solution( dico_correction_noeud_z["tuple_cliques_pi1_pi2"][0], dico_graphe["ens_C_i"])
        dico_graphe["ens_C_i"].append(dico_correction_noeud_z["tuple_cliques_pi1_pi2"][0])
    if len( dico_correction_noeud_z["tuple_cliques_pi1_pi2"][1]) != 0:
        dico_graphe["ens_C_i"] = supprimer_cliq_couvert_min_solution( dico_correction_noeud_z["tuple_cliques_pi1_pi2"][1], dico_graphe["ens_C_i"])
        dico_graphe["ens_C_i"].append( dico_correction_noeud_z["tuple_cliques_pi1_pi2"][1])
    
    #print("Ajout ens_C_i: ", dico_graphe["ens_C_i"])
    return dico_graphe["ens_C_i"], dico_correction_noeud_z ["cliques_a_supprimer_ens_C"], dico_graphe["aretes_E_C_i"], dico_graphe["dico_cliq"]


def appliquer_correction_new( noeud_z, dico_noeuds_cout_degre, dico_graphe, matE):
    """
    but : appliquer les modifs trouves au noeud "noeud_z" sur tout le graphe.
    NB: 
    dico_correction_noeud_z = {cout:, tuple_cliques_pi1_pi2:, aretes_a_ajouter:, aretes_a_supprimer:, cliques_a_supprimer_ens_C:}
    dico_graphe = {"ens_C_i": , "aretes_E_C_i": aretes_E_C_i, "aretes_E0": aretes_E_0}
    """
    dico_correction_noeud_z = dico_noeuds_cout_degre[noeud_z][2]
    #print("*** 1 appliquer_correction_new")
    if dico_correction_noeud_z["tuple_cliques_pi1_pi2"] == None:
        return dico_graphe["ens_C_i"], [], dico_graphe["aretes_E_C_i"], dico_graphe["dico_cliq"], matE;
    
    #print("*** 2 appliquer_correction_new")
    # supprimer aretes dans matE et dans liste_arcs_E_C_i
    for arete in dico_correction_noeud_z["aretes_a_supprimer"]:
        if (arete[0],arete[1]) in dico_graphe["aretes_E_C_i"]:
            dico_graphe["aretes_E_C_i"].remove( (arete[0],arete[1]) )
            for Clik in dico_graphe["ens_C_i"]:
                if arete[0] in Clik and arete[1] in Clik:
                    dico_graphe["ens_C_i"].remove(Clik)
        if (arete[1],arete[0]) in dico_graphe["aretes_E_C_i"]:
            dico_graphe["aretes_E_C_i"].remove( (arete[1],arete[0]) )
            for Clik in dico_graphe["ens_C_i"]:
                if arete[0] in Clik and arete[1] in Clik:
                    dico_graphe["ens_C_i"].remove(Clik)
        #print("----arete(s) a supprimer= ", arete)
        #print("----(avant)= mat[",arete[0],",",arete[1],"]=", matE.loc[arete[0]][arete[1]])
        matE.loc[arete[0]][arete[1]] = 0; matE.loc[arete[1]][arete[0]] = 0;
        #print("----(avant)= mat[",arete[0],",",arete[1],"]=", matE.loc[arete[0]][arete[1]])
        
    # ajouter aretes dans matE et dans liste_arcs_E_C_i
    dico_graphe["aretes_E_C_i"].extend(dico_correction_noeud_z["aretes_a_ajouter"])
    for arete in dico_correction_noeud_z["aretes_a_ajouter"]:
        #print("----arete(s) ajoute= ", arete)
        matE.loc[arete[0]][arete[1]] = 1; matE.loc[arete[1]][arete[0]] = 1;
    
    # supprimer les cliques de min_cliques_aSupprDe_ens_C dans ens_C ===> Ne pas oublier
    for clik_a_supprimer in dico_correction_noeud_z["cliques_a_supprimer_ens_C"]:
        if  clik_a_supprimer in dico_graphe["ens_C_i"]:
            dico_graphe["ens_C_i"].remove(clik_a_supprimer)
            #print("3 cliques_a_Supprimer = ", clik_a_supprimer)
            
    cliques_a_ajouter = list()
    if len(dico_correction_noeud_z["tuple_cliques_pi1_pi2"][0]) != 0:
        dico_graphe["ens_C_i"] = supprimer_cliq_couvert_min_solution( dico_correction_noeud_z["tuple_cliques_pi1_pi2"][0], dico_graphe["ens_C_i"])
        cliques_a_ajouter.append(dico_correction_noeud_z["tuple_cliques_pi1_pi2"][0])
    if len( dico_correction_noeud_z["tuple_cliques_pi1_pi2"][1]) != 0:
        dico_graphe["ens_C_i"] = supprimer_cliq_couvert_min_solution( dico_correction_noeud_z["tuple_cliques_pi1_pi2"][1], dico_graphe["ens_C_i"])
        cliques_a_ajouter.append( dico_correction_noeud_z["tuple_cliques_pi1_pi2"][1])
    
    som_cout = 0;
    som_cout += dico_correction_noeud_z["cout"] if dico_correction_noeud_z["cout"] != None else 0;
    return matE, cliques_a_ajouter, som_cout;    
    
import simulation_50_graphes_PARALLELE as simu_50;
def corriger_noeud(liste_noeuds_1, ens_C_i, dico_cliq, liste_aretes_E0, \
                   matE, number_items_pi1_pi2, critere_selection_noeuds_1, \
                   dico_proba_cases, coef_fct_cout, critere_selection_pi1_pi2):
    """
    selectionne un noeud puis applique les corrections autour de ce noeud
    """
    liste_noeuds_non_traites = liste_noeuds_1.copy();
    dico_noeuds_cout_degre = dict();
    liste_aretes_E_C_i = liste_aretes_E0.copy();
    #print("1");
    dico_correction_noeud_z = dict()
    for noeud_z in liste_noeuds_non_traites:
        dico_correction_noeud_z = calcul_cout_correction(noeud_z, ens_C_i, \
                                                         liste_aretes_E_C_i, \
                                                         liste_aretes_E0, dico_cliq, \
                                                         number_items_pi1_pi2, \
                                                         dico_proba_cases, coef_fct_cout, \
                                                         critere_selection_pi1_pi2)
        degre_z = fct_aux.degre_noeud(liste_aretes_E_C_i, noeud_z)
        #print("noeud_z = ",noeud_z, "degre_z= ", degre_z," dico_correction_noeud_z",dico_correction_noeud_z)
        dico_noeuds_cout_degre[noeud_z] = [dico_correction_noeud_z["cout"], degre_z, dico_correction_noeud_z]
    
        dico_correction_noeud_z = dico_noeuds_cout_degre[noeud_z][2]
    #print("*** dico_correction_noeud_z = ", dico_correction_noeud_z)        
    noeud_z = choix_noeud_a_corriger(dico_noeuds_cout_degre, critere_selection_noeuds_1)
    #print("*** noeud_choisi ", noeud_z)    
    dico_graphe = {"ens_C_i": ens_C_i, "aretes_E_C_i": liste_aretes_E_C_i, \
                        "aretes_E0": liste_aretes_E0, "dico_cliq":dico_cliq}
    som_cout = 0;
    matE_old = matE.copy();
    matE, cliques_a_ajouter, som_cout = appliquer_correction_new( noeud_z, dico_noeuds_cout_degre, \
                                                       dico_graphe, matE.copy())
    nb_dh, l_dh = simu_50.distance_hamming(fct_aux.liste_arcs(matE), fct_aux.liste_arcs(matE_old))
#    print("noeud_z:",noeud_z," nb_dh:",nb_dh," l_dh:",l_dh)
    return noeud_z, matE, cliques_a_ajouter, som_cout ;
############### glouton  Fin #############################################################################


############### debut selection noeuds -1 par cout minimum #############################################################################
def calculer_cout_noeud_z(noeud_z, ens_C_i, liste_aretes_E_C_i, liste_aretes_E0, \
                          number_items_pi1_pi2, dico_proba_cases, coef_fct_cout, \
                          correl_seuil, critere_selection_pi1_pi2):
    """
    but: calcul le cout de la correction d'un noeud
    return  som_cout = min_cout: cout correction du noeud, 
            min_tuple_solution: les 2 cliques qui correspondent a la correction du noeud "noeud_z", contiennent aussi le noeud "noeud_z", 
            min_liste_arete_Aajouter: liste des aretes a ajouter pdt la correction de noeud "noeud_z", 
            min_liste_aretes_aSuppr: liste des aretes a supprimer pdt la correction de noeud "noeud_z", 
            min_cliques_aSupprDe_ens_C: liste des cliques a supprimer de l'ensemble de cliques "ens_C_i" pdt la correction de noeud "noeud_z"
    dico_cout_noeud_z = {"noeud_z":,"cout":,"degre":,"tuple_cliques_pi1_pi":,"aretes_a_ajouter":,\
                         "aretes_a_supprimer":,"cliques_a_supprimer_ens_C":,"":,}
    """
    gamma_z = fct_aux.voisins(liste_aretes_E_C_i, noeud_z)
    c_z = C_z( noeud_z, ens_C_i)
    liste_aretes_ens_C_i = aretes_ens_C_i(ens_C_i)
    l_cliq_contr = cliques_contractables(noeud_z, c_z, ens_C_i, liste_aretes_E_C_i, liste_aretes_ens_C_i)
    s_z = S_z(noeud_z, gamma_z, ens_C_i, liste_aretes_ens_C_i)
        
    l_augm_z = augmentation_z(noeud_z, gamma_z, c_z, s_z, ens_C_i, liste_aretes_E_C_i, 2)
    l_C1_S1 = set_C1_S1(noeud_z, l_cliq_contr, s_z, ens_C_i)
    l_pi1_pi2 = pi1_pi2(noeud_z, l_augm_z, l_C1_S1, number_items_pi1_pi2)
    l_4uplet_pi1_pi2_pis_cliqsASupp = pis(noeud_z, gamma_z, l_pi1_pi2, ens_C_i)
        
#    print("gamma_z : ", gamma_z," c_z : ", c_z," s_z : ", s_z," l_cliq_contr : ",l_cliq_contr," list_augm_z: ",l_augm_z)

    min_tuple_solution = None; min_cout = 0;
    min_cout, min_tuple_solution, min_liste_arete_Aajouter, min_liste_aretes_aSuppr, min_cliques_aSupprDe_ens_C = \
    cout_min(noeud_z, l_4uplet_pi1_pi2_pis_cliqsASupp, liste_aretes_E_C_i, \
             liste_aretes_E0, dico_proba_cases, coef_fct_cout, correl_seuil,\
             critere_selection_pi1_pi2)
    
    dico_cout_noeud_z = dict()
    dico_cout_noeud_z["noeud_z"] = noeud_z;
    dico_cout_noeud_z["cout"] = min_cout;
    dico_cout_noeud_z["degre"] = len(gamma_z)
    dico_cout_noeud_z["tuple_cliques_pi1_pi2"] = min_tuple_solution;
    dico_cout_noeud_z["aretes_a_ajouter"] = min_liste_arete_Aajouter;
    dico_cout_noeud_z["aretes_a_supprimer"] = min_liste_aretes_aSuppr;
    dico_cout_noeud_z["cliques_a_supprimer_ens_C"] = min_cliques_aSupprDe_ens_C; 
    
    return dico_cout_noeud_z
    
def calculer_cout_corrections(noeuds_1, ens_C_i, l_aretes_E_C_i, l_aretes_E0, \
                              number_items_pi1_pi2, dico_proba_cases, \
                              coef_fct_cout, correl_seuil, critere_selection_pi1_pi2):
    """
    calculer les couts de tous les noeuds dans noeuds_1 et les mettre ds un dico dico_cout
    dico_cout = {cout: [dico_cout_noeud_u]} avec 
    dico_cout_noeud_z = {"noeud_z":,"cout":,"degre":,"":,"":,"":,"":,}
    """
    dico_cout = dict()
    for noeud_z in noeuds_1:
        dico_cout_noeud_z = dict();
        dico_cout_noeud_z = calculer_cout_noeud_z(noeud_z, ens_C_i, l_aretes_E_C_i, \
                                                  l_aretes_E0, number_items_pi1_pi2, \
                                                  dico_proba_cases, coef_fct_cout, \
                                                  correl_seuil, critere_selection_pi1_pi2)
        if len(dico_cout) == 0 or dico_cout_noeud_z["cout"] not in dico_cout_noeud_z.keys():
            dico_cout[dico_cout_noeud_z["cout"]] = [dico_cout_noeud_z];
        elif dico_cout_noeud_z["cout"] in dico_cout_noeud_z.keys():
            dico_cout[dico_cout_noeud_z["cout"]].append(dico_cout_noeud_z);
    return dico_cout;

def selectionner_noeud_z(dico_cout_noeuds):
    """
    selctionner le noeud_z selon 
        * soit le cout minimum
        * soit le cout et le degree minimum
    le critere est le mode de selection des noeuds_z
        * critere = 1 ==> cout minimum
        * critere = 2 ==> degre minimum
        * critere = 3 ==> cout et le degree minimum
        
    dico_cout_noeud_z = {"noeud_z":,"cout":,"degre":,"tuple_cliques_pi1_pi":,"aretes_a_ajouter":,\
                         "aretes_a_supprimer":,"cliques_a_supprimer_ens_C":,"":,}
    dico_cout_noeuds[cout] = [dico_cout_noeud_z,...]                     
    """
    critere = 1;
    
    if critere == 1:
        min_liste_dico_couts = dico_cout_noeuds[min(dico_cout_noeuds.keys())]
        return min_liste_dico_couts[0]
    elif critere == 3:
        min_liste_dico_couts = dico_cout_noeuds[min(dico_cout_noeuds.keys())]
        min_degre = pow(10,9); min_dico_cout_degre = dict()
        for dico_cout in min_liste_dico_couts:
            if min_degre > dico_cout["degre"]:
                min_degre = dico_cout["degre"]
                min_dico_cout_degre = dico_cout;
        return min_dico_cout_degre
        
    elif critere == 2:
        pass
    
def appliquer_correction_noeud_z(dico_cout_noeud_z, dico_graphe, dico_cliq ):
    """
    but : appliquer les modifs trouves au noeud "noeud_z" sur tout le graphe.
    NB: 
    dico_cout_noeud_z = {cout:, tuple_cliques_pi1_pi2:, aretes_a_ajouter:, aretes_a_supprimer:, cliques_a_supprimer_ens_C:}
    dico_graphe = {"ens_C_i": ens_C_i, "aretes_E_C_i": aretes_E_C_i, "aretes_E0": aretes_E_0}
    
    return cout, ens_C_i, l_aretes_E_C_i, dico_cliq
    """
    if dico_cout_noeud_z["tuple_cliques_pi1_pi2"] == None:
        return dico_cout_noeud_z["cout"],dico_graphe["ens_C_i"], \
                dico_graphe["aretes_E_C_i"],dico_cliq;
    
    # supprimer aretes dans matE et dans liste_arcs_E_C_i
    for arete in dico_cout_noeud_z["aretes_a_supprimer"]:
        if (arete[0],arete[1]) in dico_graphe["aretes_E_C_i"]:
            dico_graphe["aretes_E_C_i"].remove( (arete[0],arete[1]) )
            for Clik in dico_graphe["ens_C_i"]:
                if arete[0] in Clik and arete[1] in Clik:
                    dico_graphe["ens_C_i"].remove(Clik)
        if (arete[1],arete[0]) in dico_graphe["aretes_E_C_i"]:
            dico_graphe["aretes_E_C_i"].remove( (arete[1],arete[0]) )
            for Clik in dico_graphe["ens_C_i"]:
                if arete[0] in Clik and arete[1] in Clik:
                    dico_graphe["ens_C_i"].remove(Clik)
                    
        if arete in dico_graphe["aretes_E_C_i"]:
            dico_graphe["aretes_E_C_i"].remove(arete)
            print("----noeud_z = ", dico_cout_noeud_z["noeud_z"]," arete(s) ", arete," supprimes")
        
        
    # ajouter aretes dans matE et dans liste_arcs_E_C_i
    dico_graphe["aretes_E_C_i"].extend(dico_cout_noeud_z["aretes_a_ajouter"])
    
    # supprimer les cliques de min_cliques_aSupprDe_ens_C dans ens_C ===> Ne pas oublier
    for clik_a_supprimer in dico_cout_noeud_z["cliques_a_supprimer_ens_C"]:
        if  clik_a_supprimer in dico_graphe["ens_C_i"]:
            dico_graphe["ens_C_i"].remove(clik_a_supprimer)
            #print("3 cliques_a_Supprimer = ", clik_a_supprimer)
            
    if len(dico_cout_noeud_z["tuple_cliques_pi1_pi2"][0]) != 0:
        dico_graphe["ens_C_i"] = supprimer_cliq_couvert_min_solution( dico_cout_noeud_z["tuple_cliques_pi1_pi2"][0], dico_graphe["ens_C_i"])
        dico_graphe["ens_C_i"].append(dico_cout_noeud_z["tuple_cliques_pi1_pi2"][0])
    if len( dico_cout_noeud_z["tuple_cliques_pi1_pi2"][1]) != 0:
        dico_graphe["ens_C_i"] = supprimer_cliq_couvert_min_solution( dico_cout_noeud_z["tuple_cliques_pi1_pi2"][1], dico_graphe["ens_C_i"])
        dico_graphe["ens_C_i"].append(dico_cout_noeud_z["tuple_cliques_pi1_pi2"][1])
    
    return dico_cout_noeud_z["cout"], dico_graphe["ens_C_i"], dico_graphe["aretes_E_C_i"], dico_cliq; 
    
def corriger_noeuds_1(noeuds_1, ens_C, l_aretes_E0, dico_cliq, \
                      ordre_noeuds_traites, number_items_pi1_pi2, \
                      dico_proba_cases, coef_fct_cout, correl_seuil, \
                      critere_selection_pi1_pi2):
    """
    1. calculer les couts des noeuds a -1
    2. selection du noeud noeud_u dont le cout est minimal
    3. supprimer noeud_u de noeuds_1
    refaire 1,2,3 tant que noeuds_1 != 0
    """
    ens_C_i = ens_C.copy(); l_aretes_E_C_i = l_aretes_E0.copy(); 
    l_noeuds_traites = []; som_cout_min = 0; 
    min_cliques_aSupprDe_ens_C = []
    while noeuds_1:
        #TODO use dico_cliq to select labelled nodes to -1
        dico_cout_noeuds = dict();
        dico_cout_noeuds = calculer_cout_corrections(noeuds_1, ens_C_i, l_aretes_E_C_i, \
                                                     l_aretes_E0.copy(),number_items_pi1_pi2,\
                                                     dico_proba_cases, coef_fct_cout, \
                                                     correl_seuil, critere_selection_pi1_pi2)
        dico_cout_noeud_z = selectionner_noeud_z(dico_cout_noeuds)
        
        dico_graphe = {"ens_C_i": ens_C_i, "aretes_E_C_i": l_aretes_E_C_i, "aretes_E0": l_aretes_E0}
        cout, ens_C_i, l_aretes_E_C_i, dico_cliq = \
        appliquer_correction_noeud_z(dico_cout_noeud_z, dico_graphe, dico_cliq)
        l_noeuds_traites.append(dico_cout_noeud_z["noeud_z"])
        som_cout_min += cout
        noeuds_1 = [ x for x in noeuds_1 if x != dico_cout_noeud_z["noeud_z"]]
#        print("ici noeud_z = ",  dico_cout_noeud_z["noeud_z"], "cout = ", cout," som_cout_min = ", som_cout_min)
        
    noeuds_traites = tuple(l_noeuds_traites)
    return ens_C_i, dico_cliq, l_aretes_E_C_i, min_cliques_aSupprDe_ens_C, \
            noeuds_traites, ordre_noeuds_traites, som_cout_min, ens_C;  
    pass
############### fin selection noeuds -1 par cout minimum #############################################################################


if __name__ == '__main__':
    start= time.time()
    chemin_datasets = "data/datasets/"
    nbre_ts = 10
    effet_joule = 0.1
    epsilon = 1#0.9# 0.75
    seuil_U = 0.5
    nbre_lien = (2,5)
    dimMat = 5
    test = "EN TEST"
    ascendant_1 = True;
    number_items_pi1_pi2 = 1;
     
#    matE = graph_part.matriceE_particuliere()
#    liste_aretes = graph_part.G0_k(matE, 0)
#    dico_cliq = dict()
#    for noeud in matE.columns.tolist():
#        dico_cliq[noeud] = -1
#        
#    C =[]; som_min = 0;
#    liste_noeuds_1 = decouvClique.liste_noeuds_1(liste_aretes, dico_cliq, ascendant_1 = 1)
#    liste_noeuds_1 = ['B', 'I', 'C', 'E', 'A', 'H', 'G', 'K', 'J', 'L', 'D', 'F']
#    liste_noeuds_1 = ['A', 'E', 'B', 'H', 'K', 'L', 'G', 'J', 'D', 'I', 'C', 'F'] # som = 7
#    #liste_noeuds_1.reverse() # som = 7
#    liste_noeuds_1 = ['A', 'E', 'B', 'H', 'L', 'G', 'J', 'I', 'D', 'C', 'F', 'K']
#    noeud_z, matE_, cliques_a_ajouter, som_cout = corriger_noeud(liste_noeuds_1, C, \
#                                                       dico_cliq, liste_aretes,\
#                                                       matE, number_items_pi1_pi2) 
     
#    #print("liste_aretes= ", liste_aretes)
#    #print("liste_arcs= ",liste_arcs)
#    #print("aretes retires= ", [x for x in liste_arcs if x not in liste_aretes])
#     
#    #print("ens_C = ",ens_C)
#    #print("som_min = ",som_min)
#    print (time.time() - start)
    noeuds_1 = ['A', 'E', 'B', 'H', 'L', 'G', 'J', 'I', 'D', 'C', 'F', 'K']
    ens_C =[]; som_min = 0;
    matE = graph_part.matriceE_particuliere()
    l_aretes_E0 = graph_part.G0_k(matE, 0);
    ordre_noeuds_traites = []; number_items_pi1_pi2 = 1;
    dico_cliq = dict()
    for noeud in matE.columns.tolist():
        dico_cliq[noeud] = -1
    
#    ens_C_i, dico_cliq, l_aretes_E_C_i, min_cliques_aSupprDe_ens_C, \
#    noeuds_traites, ordre_noeuds_traites, som_cout_min, ens_C =\
#    corriger_noeuds_1(noeuds_1, ens_C, l_aretes_E0, dico_cliq, ordre_noeuds_traites, number_items_pi1_pi2)
#    print(" cout = ", som_cout_min)
#    print(" noeuds_traites= ", noeuds_traites)