#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 15:23:36 2018

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

################### fonctions de bases pour la correction ==> debut ###########
def mise_a_jour_cliques(C_new, sommets_a_corriger, dico_sommets_par_cliqs):
    """ mettre a jour les sommets par cliques puis 
        verifier les sommets couverts par plus de deux cliques.
        
    """
    
    dico_sommets_par_cliqs_new = fct_aux.couverture_par_sommets(C_new);
    dico_sommets_corriges = dict(); dico_sommets_non_corriges = dict();
    for id_sommet, sommet_a_corriger in enumerate(sommets_a_corriger):
        cliques_sommet_a_corr = dico_sommets_par_cliqs_new[sommet_a_corriger];
        if len(cliques_sommet_a_corr) == 0:
            dico_sommets_non_corriges[id_sommet] = sommet_a_corriger;
        elif len(cliques_sommet_a_corr) == 1:
            dico_sommets_non_corriges[id_sommet] = sommet_a_corriger;
        elif len(cliques_sommet_a_corr) == 2:
            dico_sommets_corriges[id_sommet] = sommet_a_corriger;
        elif len(cliques_sommet_a_corr) > 2:
            dico_sommets_non_corriges[id_sommet] = sommet_a_corriger;
    return dico_sommets_corriges, dico_sommets_par_cliqs_new;
    
def aretes_differente(aretes_Ec, aretes_cible):
    """ retourner le nombre d'aretes differente entre aretes_Ec, aretes_cible. """
    res = set()
    for arete in aretes_cible:
        if (arete[0], arete[1]) not in aretes_Ec or \
            (arete[1], arete[0]) not in aretes_Ec:
            res.add((arete[0], arete[1]))
    return res;
    
def aretes_dans_cliques(C):
    """ retourne les aretes de tous les cliques. """
    aretes_cliques = list()
    for c in C:
        aretes_cliques.extend(it.combinations(c,2));
    return aretes_cliques;

def cliques_sommet(sommet_z, dico_sommets_par_cliqs):
    """ retourne les cliques contenant le sommet sommet_z. 
        Elle est note C(z) dans le manuscrit de these.
        
    """
    if not dico_sommets_par_cliqs[sommet_z]:
        return [];
    else:
        cliques = list();
        for cliq in dico_sommets_par_cliqs[sommet_z]:
            if len(cliq) > 2:
                cliques.append(cliq);
        return cliques;

def S_sommet(sommet_z, gamma_z, aretes_Ec, C, aretes_cliques):
    """ voisins v de sommet_z tels que 
        * {v, sommet_z} est une clique de C
        * {v, sommet_z} de aretes_Ec n'est couverte par aucune clique de C.
        
    """
    S_z = list();
    for voisin_z in gamma_z :
        if {voisin_z, sommet_z} in C :
            S_z.append(voisin_z);
        elif ((voisin_z, sommet_z) in aretes_Ec or \
            (sommet_z, voisin_z) in aretes_Ec) and \
            ((voisin_z, sommet_z) not in aretes_cliques and \
            (sommet_z, voisin_z) not in aretes_cliques):
            S_z.append(voisin_z);
    return S_z;
    
def is_contractable(clique1, clique2, aretes_Ec, aretes_cliques):
    """ determine si deux cliques sont contractables. 
    
    if true : cliques 1 et 2 sont contractables
    if false : sinon.
    """
    boolean_contractable = True;
    for noeud1, noeud2 in it.product(clique1, clique2):
        if ((noeud1, noeud2) in aretes_Ec or \
            (noeud2, noeud1) in aretes_Ec) and \
            ((noeud1, noeud2) in aretes_cliques or \
            (noeud2, noeud1) in aretes_cliques):
            boolean_contractable = False;
            break;
    return boolean_contractable;
        
def cliques_contractables(sommet_z, aretes_Ec, aretes_cliques, cliques_sommet_z):
    """ retourne la liste des cliques contractables autour du sommet_z. """

    cliq_contractables = [];
    for c1, c2 in it.combinations(cliques_sommet_z.append({})):
        if not c1 or not c2 :
            cliq_contractables.append((c1, c2));
        else:
            if is_contractable(c1, c2, aretes_Ec, aretes_cliques):
                cliq_contractables.append((c1, c2))
    return cliq_contractables;            
    
def voisine_sommet(sommet_z, C, cliques_sommet_z, s_z):
    """ cliques voisines du sommet_z. """
#    C = [frozenset(c) for c in C];
    cliques_sommet_z = [frozenset(c) for c in cliques_sommet_z]
    cliques_not_sommet_z = set(C).union(set(cliques_sommet_z)) - \
                            set(C).intersection(set(cliques_sommet_z))
                            
    cliques_voisines = [c for c in cliques_not_sommet_z 
                        if len(c.intersection(set(s_z))) >=1]
    return cliques_voisines;

def dependance_sommet(sommet_z, gamma_z, cliques_sommet_z, clique_voisine):
    """ retourner les cliques dependantes d une clique clique_voisine. """
    return [cliq for cliq in cliques_sommet_z \
            if len(cliq.intersection(clique_voisine.intersection(gamma_z))) == 0]

def augmentation(sommet_z, gamma_z, cliques_sommet_z, s_z, args):
    """ retourne des augmentations possibles autour du sommet sommet_z. """
    
    # cliques_sommet_z = cliques_sommet(sommet_z, args["dico_sommets_par_cliqs"]);
    cpt = 0;
    dico_cliques_augmentante = dict();
    for cliq in args["C"]:
        cliques_voisine = voisine_sommet(sommet_z, cliq, cliques_sommet_z, s_z);
        for clique_voisine in cliques_voisine:
            cliques_dependante = list();
            cliques_dependante = dependance_sommet(sommet_z, gamma_z, \
                                                     cliques_sommet_z, \
                                                     clique_voisine)
            if not cliques_dependante:
#                dico_cliques_augmentante[cpt] = {"cliq":cliq, 
#                                                "voisine":clique_voisine,
#                                                "dependante":frozenset(),
#                                                "sommet_z":sommet_z}
                dico_cliques_augmentante[(cpt, clique_voisine,\
                                                 frozenset())] = {
                                                "cliq":cliq, 
                                                "voisine":clique_voisine,
                                                "dependante":frozenset(),
                                                "sommet_z":sommet_z}                                  
                cpt += 1;
            else:
                for clique_dependante in cliques_dependante:
                    if is_contractable(clique_voisine, 
                                       clique_dependante, 
                                       args["aretes_Ec"], 
                                       args["aretes_cliques"]):
#                        dico_cliques_augmentante[cpt] = {"cliq":cliq, 
#                                                "voisine":clique_voisine,
#                                                "dependante":clique_dependante,
#                                                "sommet_z":sommet_z}
                        dico_cliques_augmentante[(cpt, clique_voisine,\
                                                  clique_dependante)] = {
                                                "cliq":cliq, 
                                                "voisine":clique_voisine,
                                                "dependante":clique_dependante,
                                                "sommet_z":sommet_z}                        
                        cpt += 1;
    
    return dico_cliques_augmentante;              

def compression_sommet(id_sommet_z, sommet_z, sommets_a_corriger, 
                       cliques_sommet_z, args):
    """ retourne la compression d'un sommet sommet_z. 
    
    la compression est le triplet (pi1, pi2, ps) dans lequel 
        * pi1, pi2 sont des cliques qui fusionnent 
            - des cliques augmentantes C1, C2 ou 
            - des cliques contractables C1, C2 ou 
            - un ensemble S1 tel que S1 n'est contenu par aucune clique C1 ou C2
        * pi1, pi2 sont des augmentations
        * ps est un ensemble de sommets u tel que (z,u) doit etre supprime de aretes_Ec
        
    """
    
    s_z = S_sommet(sommet_z, 
                   args["dico_gamma_sommets"][sommet_z][1], 
                   args["aretes_Ec"], 
                   args["C"], 
                   args["aretes_cliques"]);
    
    # determination de C1 = (C_1,C_2) avec C_1, C_2 contratables
    dico_C1_C2_S1 = dict(); cpt = 0;
    for C1, C2 in cliques_contractables(sommet_z, 
                                       args["aretes_Ec"], 
                                       args["aretes_cliques"], 
                                       cliques_sommet_z):
        S1 = C1.union(C2) - C1.union(C2).intersection(s_z);
        bool_sommet_a_exclu = True; S1_new = frozenset();
        for sommet_S1 in S1:
            for s1_, c1_c2 in it.product(frozenset(sommet_S1), C1.union(C2)):
                if frozenset(s1_, c1_c2) not in args["C"]:
                    bool_sommet_a_exclu = False;
                    break;
            if bool_sommet_a_exclu :
                S1_new.union(set(s1_))
        dico_C1_C2_S1[(cpt, C1, C2, S1_new)] = {
                      "cliques_contratables":(C1,C2),
                      "S1":S1,
                      "clique_possible": 
                          C1.union(C2.union(S1_new.union(frozenset(sommet_z))))
                                            }
        cpt += 1;
    
    # determination de pi1_pi2_ps
    dico_cliques_augmentante = dict();
    dico_cliques_augmentante = augmentation(
                                    sommet_z,
                                    args["dico_gamma_sommets"][sommet_z][1], 
                                    cliques_sommet_z, 
                                    s_z, 
                                    args);
    nb_prod_cartesien = pow(len(dico_C1_C2_S1), len(dico_cliques_augmentante))
    nbre_elts_pi1_pi2 = math.ceil( nb_prod_cartesien * 
                                  args["number_items_pi1_pi2"])
    cpt_prod_cartesien = 0;
    dico_p1_p2_ps = dict();
    #TODO NOK: A refaire car pas de melange de solution """
    ##################33 test combinaision de dico
#    """
#    dico_C1_C2_S1[(cpt, C1, C2, S1_new)] = {
#                      "cliques_contratables":,"S1":,"clique_possible": }
#    dico_cliques_augmentante[(cpt, clique_voisine,\
#                             clique_dependante)] = {
#                              "cliq":, "voisine":,
#                              "dependante":,"sommet_z":}  
#    """                                        
#    for dico_c1c2s1_augm in it.islice(map(dict, it.product(dico_C1_C2_S1.items(), 
#                                         dico_cliques_augmentante.items())),
#                              nbre_elts_pi1_pi2):  
    ##################33 test combinaision de dico
    for k_c1_c2_s1, val_cpt_c1_c2_s1 in dico_C1_C2_S1.items():
        for k_cpt_vois_depend, val_cpt_vois_depend in dico_cliques_augmentante.items():                                        
            cpt_prod_cartesien += 1;
            
            inter_p1_p2 = val_cpt_c1_c2_s1["clique_possible"].intersection(
                            k_cpt_vois_depend[1].union(k_cpt_vois_depend[2])
                            )
            if len(inter_p1_p2) <= 1 and inter_p1_p2 == frozenset(sommet_z):
                
                p1 = val_cpt_c1_c2_s1["clique_possible"];
                p2 = val_cpt_vois_depend["voisine"].union(
                            val_cpt_vois_depend["dependante"].union(
                            frozenset(sommet_z)))
                gamma_z = args["dico_gamma_sommets"][1];
                ps = gamma_z - p1.intersection(gamma_z).union(
                                p2.intersection(gamma_z));
                aretes_ps = frozenset( frozenset((sommet_z, sommet_ps)) 
                                        for sommet_ps in ps)
                # TODO OK: calcul nombre aretes a ajouter pour p1, p2, nombre sommets a corriger
                aretes_p1 = aretes_dans_cliques(p1);
                aretes_ajoute_p1 = aretes_differente(args["aretes_Ec"], 
                                                     aretes_p1);
                
                aretes_p2 = aretes_dans_cliques(p2);
                aretes_ajoute_p2 = aretes_differente(args["aretes_Ec"], 
                                                     aretes_p2);
                
                aretes_Ec_new = set(args["aretes_Ec"]) + aretes_ajoute_p1 + \
                                aretes_ajoute_p2;
                C_new = args["C"].copy();
                #TODO NOK: a verifier C_new les cliques retirees
                C_new = set(C_new) - val_cpt_c1_c2_s1["cliques_contratables"][0] - \
                val_cpt_c1_c2_s1["cliques_contratables"][1] - \
                val_cpt_vois_depend["voisine"] - \
                val_cpt_vois_depend["dependante"]
                
                C_new.add( p1 );
                C_new.add( p2 );
                
                dico_sommets_corriges = dict()
                dico_sommets_non_corriges = dict(), \
                dico_sommets_par_cliqs_new = dict();
                dico_sommets_corriges, dico_sommets_non_corriges, \
                dico_sommets_par_cliqs_new = \
                    mise_a_jour_cliques(C_new,
                                        sommets_a_corriger, 
                                        args["dico_sommets_par_cliqs"])
                
                dico_p1_p2_ps[cpt_prod_cartesien] = {
                    "id_sommet_1": id_sommet_z,
                    "sommet_1": sommet_z,
                    "p1": val_cpt_c1_c2_s1["clique_possible"],
                    "p2": val_cpt_vois_depend["voisine"].union(
                            val_cpt_vois_depend["dependante"].union(
                            frozenset(sommet_z))),
                    "ps": ps,
                    "voisine": val_cpt_vois_depend["voisine"],
                    "dependante": val_cpt_vois_depend["dependante"],
                    "contractable1": val_cpt_c1_c2_s1["cliques_contratables"][0],
                    "contractable2": val_cpt_c1_c2_s1["cliques_contratables"][1],
                    "S1": val_cpt_c1_c2_s1["S1"],
                    "S_z": s_z,
                    "aretes_ajoute_p1": aretes_ajoute_p1,
                    "aretes_ajoute_p2": aretes_ajoute_p2,
                    "aretes_supprimes_ps": aretes_ps,
                    "aretes_Ec_new": aretes_Ec_new,
                    "C_new": C_new,
                    "sommets_corriges": dico_sommets_corriges,
                    "sommets_non_corriges": dico_sommets_non_corriges,
                    "dico_sommets_par_cliqs_new": dico_sommets_par_cliqs_new
                    }
            if cpt_prod_cartesien <= nbre_elts_pi1_pi2:
                break;
        if cpt_prod_cartesien <= nbre_elts_pi1_pi2:
            break;
            
    return dico_p1_p2_ps;
         
################### fonctions de bases pour la correction ==> fin   ###########

################### critere selection  compression ==> debut ##################
def critere_C2_C1(dico_compression, args) :
    """ selectionner le dico selon C2 puis C1
    
    C2 : la maximum de sommets corriges
    C1 : le minimum d'aretes corriges
    
    dico_compression : dictionnaire contenant les compression (p1,p2,ps) du 
                        sommet sommet_z
    """
    
    max_c2 = 0;
    min_c1 = 0;
    dico_c1_c2 = dict();
    
    # definition de C2
    if args["critere_selection_compression"] == "voisins_corriges":             # C2
        for id_sommet_sommet_z, dico_p1_p2_ps in dico_compression.items():
            if len(dico_p1_p2_ps["sommets_corriges"]) >= max_c2 :
                max_c2 = len(dico_p1_p2_ps["sommets_corriges"]);
                min_c1 = max_c2;
                if min_c1 in dico_c1_c2:
                    dico_c1_c2[min_c1] = [dico_p1_p2_ps];
                else:
                    dico_c1_c2[min_c1].append(dico_p1_p2_ps);
    # definition de C1
    elif args["critere_selection_compression"] == "nombre_aretes_corriges":     # C1
        for id_sommet_sommet_z, dico_p1_p2_ps in dico_compression.items():
            nbre_aretes_corriges = len(dico_p1_p2_ps["aretes_ajoute_p1"]) + \
                                    len(dico_p1_p2_ps["aretes_ajoute_p2"]) + \
                                    len(dico_p1_p2_ps["aretes_supprimes_ps"]);
            min_c1 = nbre_aretes_corriges if min_c1 > nbre_aretes_corriges \
                                          else min_c1;
            if min_c1 in dico_c1_c2:
                dico_c1_c2[min_c1] = [dico_p1_p2_ps];
            else:
                dico_c1_c2[min_c1].append(dico_p1_p2_ps);
        pass
    # definition de C2 puis de C1
    elif args["critere_selection_compression"] == "voisins_nombre_aretes_corriges": # C2_C1
        for id_sommet_sommet_z, dico_p1_p2_ps in dico_compression.items():
            if len(dico_p1_p2_ps["sommets_corriges"]) >= max_c2 :
                max_c2 = len(dico_p1_p2_ps["sommets_corriges"]);
                nbre_aretes_corriges = len(dico_p1_p2_ps["aretes_ajoute_p1"]) + \
                                    len(dico_p1_p2_ps["aretes_ajoute_p2"]) + \
                                    len(dico_p1_p2_ps["aretes_supprimes_ps"]);
                min_c1 = nbre_aretes_corriges if min_c1 > nbre_aretes_corriges \
                                                else min_c1;
                if min_c1 in dico_c1_c2:
                    dico_c1_c2[min_c1] = [dico_p1_p2_ps];
                else:
                    dico_c1_c2[min_c1].append(dico_p1_p2_ps);
                    
    if not dico_c1_c2:
        return min_c1, max_c2, dico_c1_c2;
    else:
        return min_c1, max_c2, dico_c1_c2[min_c1];
################### critere selection  compression ==> fin ##################

################### application de la compression ==> debut ##################
def appliquer_correction(dico_sol_C2_C1, sommets_a_corriger, args):
    """ appliquer la compression choisie dans le graphe.
    """
    C = list();
    C = dico_sol_C2_C1["C_new"];
    aretes_Ec = list();
    aretes_Ec = dico_sol_C2_C1["aretes_Ec_new"];
    
    id_sommets_1 = dico_sol_C2_C1["sommets_corriges"].keys();
    id_sommets_1.append(dico_sol_C2_C1["id_sommet_1"]);
    sommets_corriges = dico_sol_C2_C1["sommets_corriges"].values();
    sommets_a_corriger = np.delete(sommets_a_corriger, id_sommets_1).tolist();
    if set(sommets_a_corriger).intersection(set(sommets_corriges)) :
        print("---ERROR : sommets {} suppression : NOK -----".
              format(sommets_corriges))

    return C, aretes_Ec, sommets_a_corriger;
################### application de la compression ==> fin ####################
    

def correction_graphe_correlation(args):
    """ corrige un graphe de correlation en ajoutant ou supprimant des aretes
    
    """
    dico_sommets_corriges = dict();
    dico_compression = dict()
    sommets_a_corriger = list();
    sommets_a_corriger = [sommet for sommet, etat in args["dico_cliq"].items() 
                            if etat == -1]
    if args["critere_selection_compression"] == "voisins_corriges":
        # correction sans remise avec le critere "nombre de voisins corriges"
        cpt_noeud = 0;
        while(sommets_a_corriger):
            for id_sommet_1, sommet_1 in enumerate(sommets_a_corriger):
                cliques_sommet_1 = cliques_sommet(sommet_1, 
                                                  args["dico_sommets_par_cliqs"]);
                dico_compression[(id_sommet_1,sommet_1)] = compression_sommet(
                                                             id_sommet_1,
                                                             sommet_1,
                                                             sommets_a_corriger,
                                                             cliques_sommet_1,
                                                            args);
                                 
            # dico_sol_C2_C1 = {
            #        "id_sommet_1":,"sommet_1":,"p1":,"p2":,"ps":,
            #        "voisine":,"dependante":,
            #        "contractable1":,"contractable2":,
            #        "S1":,"S_z":,
            #        "aretes_ajoute_p1":,"aretes_ajoute_p2":,
            #        "aretes_supprimes_ps":,"aretes_Ec_new":,"C_new":,
            #        "sommets_corriges":,"sommets_non_corriges":,
            #        "dico_sommets_par_cliqs_new": 
            #        }
            dico_sol_C2_C1 = dict();
            min_c1 = 0; max_c2 = 0;
            min_c1, max_c2, dico_sol_C2_C1 = critere_C2_C1(dico_compression,
                                                           args)                # C2 : nombre maximum de voisins corriges par un sommet, C1 : nombre minimum d'aretes a corriger au voisinage d'un sommet  
            C, aretes_Ec, sommets_a_corriger = appliquer_correction(
                                                dico_sol_C2_C1,
                                                sommets_a_corriger,
                                                args)
            args["C"] = C;
            args["aretes_cliques"] = aretes_dans_cliques(C);
            args["aretes_Ec"] = aretes_Ec;
            cout_T = {"aretes_ajoutes_p1":dico_sol_C2_C1["aretes_ajoute_p1"],
                      "aretes_ajoutes_p2":dico_sol_C2_C1["aretes_ajoute_p2"],
                      "aretes_supprimes":dico_sol_C2_C1["aretes_supprimes_ps"],
                      "min_c1":min_c1,"max_c2":max_c2};
            cpt_noeud += 1;
            dico_sommets_corriges[(cpt_noeud, dico_sol_C2_C1["sommet_1"])] = {
                        "compression_p1":dico_sol_C2_C1["p1"],
                        "compression_p2":dico_sol_C2_C1["p2"],
                        "compression_ps":dico_sol_C2_C1["ps"],
                        "sommets_corriges":dico_sol_C2_C1["sommets_corriges"], # voisins_corriges = {"id_voisin_ds_sommets_a_corriger":voisin}
                        "cout_T": cout_T
                        }
            
            # mettre a jour les cliques couvrants les sommets.
            args["dico_sommets_par_cliqs"] = dico_sol_C2_C1[
                                                "dico_sommets_par_cliqs_new"
                                                ]
        return args, dico_sommets_corriges;
        
if __name__ == '__main__':
    ti = time.time();