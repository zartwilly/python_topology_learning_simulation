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
def mise_a_jour_aretes_cliques(C_new, aretes_Ec, aretes_ps, 
                               sommets_a_corriger, dico_sommets_par_cliqs):
    """ mettre a jour les sommets par cliques puis 
        verifier les sommets couverts par plus de deux cliques.
        
    """
#    print("PS: aretes_ps={}, cpt_aretes_Ec={}".format(aretes_ps, len(aretes_Ec)))
    # suppression des aretes_ps dans aretes_Ec
    aretes_ps = set(map(tuple, aretes_ps));
    aretes_ps_new = set();
    for arete in aretes_ps:
        aretes_ps_new.add(arete)
        aretes_ps_new.add( (arete[1],arete[0]) )
    aretes_Ec.difference_update(aretes_ps_new);
#    print("PS new cpt_aretes_Ec={}".format(len(aretes_Ec)))
    
    # suppression cliques dont on a supprime des aretes_ps
    C_nouvelle = C_new.copy();
    for c in C_new:
        for arete_ps in aretes_ps:
            if set(arete_ps).issubset(c):
                C_nouvelle.difference_update({c});
                temp_aretes = [arete_ps, (arete_ps[1],arete_ps[0])]
                aretes_Ec.difference_update(set(temp_aretes));             
                break;
#    print("PS cpt_C_new={},cpt_C_nouvelle={}".format(len(C_new), len(C_nouvelle)))
    
    dico_sommets_par_cliqs_new = fct_aux.couverture_par_sommets(C_nouvelle);
    dico_sommets_corriges = dict(); dico_sommets_non_corriges = dict();
    
    for id_sommet, sommet_a_corriger in enumerate(sommets_a_corriger):
        cliques_sommet_a_corr = dico_sommets_par_cliqs_new[sommet_a_corriger];
        gamma_sommet_a_corriger = set(fct_aux.gamma(aretes_Ec,sommet_a_corriger));
        if len(cliques_sommet_a_corr) == 0:
            dico_sommets_non_corriges[id_sommet] = sommet_a_corriger;
        elif len(cliques_sommet_a_corr) == 1 and \
            cliques_sommet_a_corr[0] == gamma_sommet_a_corriger:
            dico_sommets_corriges[id_sommet] = sommet_a_corriger;
#            print("ICI1")
        elif len(cliques_sommet_a_corr) == 1 and \
            cliques_sommet_a_corr[0] != gamma_sommet_a_corriger:
            dico_sommets_non_corriges[id_sommet] = sommet_a_corriger;
#            print("ICI2")
        elif len(cliques_sommet_a_corr) == 2:
            dico_sommets_corriges[id_sommet] = sommet_a_corriger;
        elif len(cliques_sommet_a_corr) > 2:
            dico_sommets_non_corriges[id_sommet] = sommet_a_corriger;
#    print("??? dico_sommets_corriges={},dico_sommets_non_corriges={}, dico_sommets_par_cliqs_new={}".\
#          format(dico_sommets_corriges, dico_sommets_non_corriges, dico_sommets_par_cliqs_new))
    return C_nouvelle, \
            aretes_Ec,\
            dico_sommets_corriges, \
            dico_sommets_non_corriges, \
            dico_sommets_par_cliqs_new;
    
def aretes_differente(aretes_Ec, aretes_cible):
    """ retourner le nombre d'aretes differente entre aretes_Ec, aretes_cible. """
    res = set()
    for arete in aretes_cible:
        if (arete[0], arete[1]) not in aretes_Ec and \
            (arete[1], arete[0]) not in aretes_Ec:
            res.add((arete[0], arete[1]))
#    res = aretes_Ec.union(aretes_cible) - aretes_Ec.intersection(aretes_cible)         
    return res;

def cliques_sommet(sommet_z, dico_sommets_par_cliqs):
    """ retourne les cliques contenant le sommet sommet_z. 
        Elle est note C(z) dans le manuscrit de these.
        
    """
    if not dico_sommets_par_cliqs[sommet_z]:
        return [];
    else:
        cliques = list();
        for cliq in dico_sommets_par_cliqs[sommet_z]:
            if len(cliq) >= 2:
                cliques.append(cliq);
        return cliques;

def S_sommet(sommet_z, gamma_z, aretes_Ec, C, aretes_cliques):
    """ voisins v de sommet_z tels que 
        * {v, sommet_z} est une clique de C
        * {v, sommet_z} de aretes_Ec n'est couverte par aucune clique de C.
        
    """
    logger = logging.getLogger('S_sommet');
    S_z = list();
    for voisin_z in gamma_z :
        if {voisin_z, sommet_z} in C :
            S_z.append(voisin_z);
        elif ((voisin_z, sommet_z) in aretes_Ec or \
            (sommet_z, voisin_z) in aretes_Ec) and \
            ((voisin_z, sommet_z) not in aretes_cliques and \
            (sommet_z, voisin_z) not in aretes_cliques):
            S_z.append(voisin_z);
    logger.debug(" * S_z: {}".format(S_z))
    return S_z;
    
def is_contractable_old(clique1, clique2, aretes_Ec, aretes_cliques):
    """ determine si deux cliques sont contractables. 
    
    if true : cliques 1 et 2 sont contractables
    if false : sinon.
    """
    boolean_contractable = True;
    sommet_inter = clique1.intersection(clique2)
    for noeud1, noeud2 in it.product(clique1-sommet_inter, 
                                     clique2-sommet_inter):
        if noeud1 != noeud2 and \
            ((noeud1, noeud2) in aretes_Ec or \
             (noeud2, noeud1) in aretes_Ec) and \
            ((noeud1, noeud2) in aretes_cliques or \
             (noeud2, noeud1) in aretes_cliques) and \
             ((noeud1 in clique1 and noeud2 in clique2) or \
              (noeud2 in clique1 and noeud2 in clique2)):
            boolean_contractable = False;
            break;
    return boolean_contractable;
    
def is_contractable(clique1, clique2, aretes_Ec, aretes_cliques, C):
    """ determine si deux cliques sont contractables. 
    
    if true : cliques 1 et 2 sont contractables
    if false : sinon.
    """
    boolean_contractable = True; 
    cliques_suppl_contractables = set()
    sommet_inter = clique1.intersection(clique2)
    for noeud1, noeud2 in it.product(clique1-sommet_inter, 
                                     clique2-sommet_inter):    
        if noeud1 != noeud2 and \
           ((noeud1,noeud2) in aretes_Ec or (noeud2,noeud2) in aretes_Ec) and \
           ( frozenset((noeud1,noeud2)) not in C) :
               boolean_contractable = False; 
               cliques_suppl_contractables = set()
               break;
        if frozenset((noeud1,noeud2)) in C:
           cliques_suppl_contractables.add( frozenset((noeud1,noeud2)) ) 
    return boolean_contractable, cliques_suppl_contractables;
    
def cliques_contractables(sommet_z, aretes_Ec, 
                          aretes_cliques, cliques_sommet_z, C):
    """ retourne la liste des cliques contractables autour du sommet_z. """
    
    logger = logging.getLogger('cliques_contractables');
    cliq_contractables = [];
    """
    TODO a poser la question a un forum.
    how to do itertools with Nonetype inside a list.
    for c1, c2 in it.combinations(cliques_sommet_z.append({}),2):
        if not c1 or not c2 :
            cliq_contractables.append((c1, c2));
        else:
            if is_contractable(c1, c2, aretes_Ec, aretes_cliques):
                cliq_contractables.append((c1, c2))
    return cliq_contractables;
    """
    for c1, c2 in it.combinations(cliques_sommet_z,2):
        boolean_contractable, cliques_suppl_contractables = is_contractable(
                                                            c1, 
                                                            c2, 
                                                            aretes_Ec, 
                                                            aretes_cliques, 
                                                            C)
        if boolean_contractable:
            cliq_contractables.append((c1, c2, cliques_suppl_contractables))
    for c in cliques_sommet_z:
        cliq_contractables.append((c,frozenset(), set()))
        
    logger.debug(" * cliq_contractables: {}".format( len(cliq_contractables)) )
    return cliq_contractables;            
    
def voisine_sommet(sommet_z, cliques_sommet_z, cliques_not_sommet_z, s_z) :
    """ cliques voisines du sommet_z. """
#    cliques_sommet_z = [frozenset(c) for c in cliques_sommet_z]
    logger = logging.getLogger('voisine_sommet');
    cliques_voisines = [c for c in cliques_not_sommet_z 
                        if len(c.intersection(set(s_z))) >=1]
    logger.debug(" * voisine_sommet: {}".format(len(cliques_voisines)))
    return cliques_voisines;

def dependance_sommet(sommet_z, gamma_z, cliques_sommet_z, clique_voisine):
    """ retourner les cliques dependantes d une clique clique_voisine. """
    return [cliq for cliq in cliques_sommet_z \
            if len(cliq.intersection(clique_voisine.intersection(gamma_z))) != 0]

def augmentation(sommet_z, gamma_z, cliques_sommet_z, s_z, args):
    """ retourne les augmentations possibles autour du sommet sommet_z. """
    
    logger = logging.getLogger('augmentation');
    # cliques_sommet_z = cliques_sommet(sommet_z, args["dico_sommets_par_cliqs"]);
    cpt = 0;
    dico_cliques_augmentante = dict();
    cliques_not_sommet_z = list();
    cliques_not_sommet_z = [ c for c in args["C"] 
                             if not c.intersection({sommet_z})]
    
    cliques_voisine = voisine_sommet(sommet_z, cliques_sommet_z, 
                                     cliques_not_sommet_z, s_z);
    for clique_voisine in cliques_voisine:
        cliques_dependante = list();
        cliques_dependante = dependance_sommet(sommet_z, gamma_z, \
                                                 cliques_sommet_z, \
                                                 clique_voisine)
#        print("!!!clique_voisine={},\n cliques_dependante={}".\
#              format(clique_voisine,cliques_dependante))
        if not cliques_dependante:
#                dico_cliques_augmentante[cpt] = {"cliq":cliq, 
#                                                "voisine":clique_voisine,
#                                                "dependante":frozenset(),
#                                                "sommet_z":sommet_z}
            dico_cliques_augmentante[(cpt, clique_voisine,\
                                     frozenset(), frozenset())] = {
                                      "voisine":clique_voisine,
                                      "dependante":frozenset(),
                                      "cliques_suppl_contractables": frozenset(),
                                      "sommet_z":sommet_z}                                  
            cpt += 1;
        else:
            for clique_dependante in cliques_dependante:
                boolean_contractable, cliques_suppl_contractables = \
                        is_contractable(clique_voisine, 
                                   clique_dependante, 
                                   args["aretes_Ec"], 
                                   args["aretes_cliques"],
                                   args["C"])
                if boolean_contractable:
                    dico_cliques_augmentante[(cpt, clique_voisine,\
                        clique_dependante, \
                        frozenset(cliques_suppl_contractables))] = {
                        "voisine":clique_voisine,
                        "dependante":clique_dependante,
                        "cliques_suppl_contractables": \
                            frozenset(cliques_suppl_contractables),
                        "sommet_z":sommet_z}                        
                    cpt += 1;
#    print("??? dico_cliques_augmentante={}".format(dico_cliques_augmentante));
    logger.debug(" * augmentation: {}".format(len(dico_cliques_augmentante)))
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
    logger = logging.getLogger('compression_sommet');
#    print("X cliques_sommet_z={}".format(cliques_sommet_z))
    s_z = S_sommet(sommet_z, 
                   args["dico_gamma_sommets"][sommet_z][1], 
                   args["aretes_Ec"], 
                   args["C"], 
                   args["aretes_cliques"]);
    
    # determination de C1 = (C_1,C_2) avec C_1, C_2 contratables
    dico_C1_C2_S1 = dict(); cpt = 0;
    for C1, C2, cliques_suppl_contractables in cliques_contractables(sommet_z, 
                                       args["aretes_Ec"], 
                                       args["aretes_cliques"], 
                                       cliques_sommet_z.copy(), 
                                       args["C"]):
        S1 = C1.union(C2) - C1.union(C2).intersection(s_z);
        bool_sommet_a_exclu = True; S1_new = frozenset();
        for sommet_S1 in S1:
            for s1_, c1_c2 in it.product(frozenset({sommet_S1}), C1.union(C2)):
                if frozenset({s1_, c1_c2}) not in args["C"]:
                    bool_sommet_a_exclu = False;
                    break;
            if bool_sommet_a_exclu :
                S1_new.union({s1_})
        dico_C1_C2_S1[(cpt, C1, C2, S1_new)] = {
                      "cliques_contratables":(C1,C2),
                      "cliques_suppl_contractables":cliques_suppl_contractables,
                      "S1":S1,
                      "clique_possible": 
                          C1.union(C2.union(S1_new.union(frozenset({sommet_z}))))
                                            }
        cpt += 1;
    
    # determination de pi1_pi2_ps
    dico_cliques_augmentante = dict();
    dico_cliques_augmentante = augmentation(
                                    sommet_z,
                                    args["dico_gamma_sommets"][sommet_z][1], 
                                    cliques_sommet_z.copy(), 
                                    s_z, 
                                    args);
    nb_prod_cartesien = pow(len(dico_C1_C2_S1), len(dico_cliques_augmentante)) \
                        if len(dico_C1_C2_S1) >= len(dico_cliques_augmentante) \
                        else pow(len(dico_cliques_augmentante),len(dico_C1_C2_S1))
    nbre_elts_pi1_pi2 = math.ceil( nb_prod_cartesien * 
                                  args["number_items_pi1_pi2"])
    cpt_prod_cartesien = 0;
    dico_p1_p2_ps = dict();
    print(" sommet_z ={}, ".format(sommet_z)+\
          " nbre_elts_pi1_pi2:{}, ".format(nbre_elts_pi1_pi2) + \
          " dico_C1_C2_S1:{}, ".format(len(dico_C1_C2_S1)) + \
          " dico_cliques_augmentante:{}, ".format(len(dico_cliques_augmentante)) + \
          " nbre_elts_pi1_pi2:{}".format(nbre_elts_pi1_pi2))
    logger.debug(" * compression_sommet : "+\
                 " sommet_z ={}, ".format(sommet_z)+\
                 " nbre_elts_pi1_pi2:{}, ".format(nbre_elts_pi1_pi2)+ \
                 " dico_C1_C2_S1:{}, ".format(len(dico_C1_C2_S1)) + \
                 " dico_cliques_augmentante:{}, ".format(len(dico_cliques_augmentante)) + \
                 " nbre_elts_pi1_pi2:{}".format(nbre_elts_pi1_pi2))
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
    if not dico_C1_C2_S1 and not dico_cliques_augmentante :
        logger.debug(" * compression_sommet : *** dico_C1_C2_S1, "+\
                         "dico_cliques_augmentante VIDE ")
        dico_sommets_non_corriges = dict();
        dico_sommets_corriges = dict();
        for id_sommet_1, sommet_1 in enumerate(sommets_a_corriger):
            dico_sommets_non_corriges[id_sommet_1] = sommet_1;
                
        dico_p1_p2_ps[cpt_prod_cartesien] = {
                    "id_sommet_1": id_sommet_z,
                    "sommet_1": sommet_z,
                    "p1": frozenset(),
                    "p2": frozenset(),
                    "ps": frozenset(),
                    "voisine": frozenset(),
                    "dependante": frozenset(),
                    "contractable1": frozenset(),
                    "contractable2": frozenset(),
                    "S1": frozenset(),
                    "S_z": s_z,
                    "aretes_ajoutees_p1": frozenset(),
                    "aretes_ajoutees_p2": frozenset(),
                    "aretes_supprimees_ps": frozenset(),
                    "aretes_Ec_new": args["aretes_Ec"],
                    "C_new": args["C"],
                    "sommets_corriges": dico_sommets_corriges,
                    "sommets_non_corriges": dico_sommets_non_corriges,
                    "dico_sommets_par_cliqs_new": args["dico_sommets_par_cliqs"]
                        }
    elif not dico_C1_C2_S1 and dico_cliques_augmentante :
        logger.debug(" * compression_sommet : *** dico_C1_C2_S1 VIDE, "+\
                         "dico_cliques_augmentante NON VIDE ")
        for k_cpt_vois_depend, val_cpt_vois_depend in dico_cliques_augmentante.items():
            cpt_prod_cartesien += 1;
            p1 = val_cpt_vois_depend["voisine"].union(
                                val_cpt_vois_depend["dependante"].union(
                                frozenset({sommet_z})));
            p2 = frozenset();
            gamma_z = args["dico_gamma_sommets"][sommet_z][1];
            ps = gamma_z - p1.intersection(gamma_z);
            aretes_ps = set( frozenset((sommet_z, sommet_ps)) 
                             for sommet_ps in ps)
            aretes_p1 = fct_aux.aretes_dans_cliques(p1);
            aretes_ajoutees_p1 = aretes_differente(args["aretes_Ec"], 
                                                   aretes_p1);
            aretes_Ec_new = set(args["aretes_Ec"]).union(aretes_ajoutees_p1);
            
            C_new = set(args["C"].copy());
            ens_cliq_a_supprimer = set();
            ens_cliq_a_supprimer.add(val_cpt_vois_depend["voisine"]);
            ens_cliq_a_supprimer.add(val_cpt_vois_depend["dependante"]);
            for cliq in val_cpt_vois_depend["cliques_suppl_contractables"] : 
                ens_cliq_a_supprimer.add(cliq);
            C_new.difference_update(ens_cliq_a_supprimer);
            
            C_new.add( p1 );
            dico_sommets_corriges = dict();
            dico_sommets_non_corriges = dict();
            dico_sommets_par_cliqs_new = dict();
            C_new, aretes_Ec,\
            dico_sommets_corriges, \
            dico_sommets_non_corriges, \
            dico_sommets_par_cliqs_new = \
                mise_a_jour_aretes_cliques(C_new.copy(), aretes_Ec_new.copy(), 
                                    aretes_ps,\
                                    sommets_a_corriger.copy(), \
                                    args["dico_sommets_par_cliqs"].copy())
            dico_p1_p2_ps[cpt_prod_cartesien] = {
                        "id_sommet_1": id_sommet_z,
                        "sommet_1": sommet_z,
                        "p1": p1,
                        "p2": p2,
                        "ps": ps,
                        "voisine": val_cpt_vois_depend["voisine"],
                        "dependante": val_cpt_vois_depend["dependante"],
                        "contractable1": frozenset(),
                        "contractable2": frozenset(),
                        "S1": frozenset(),
                        "S_z": s_z,
                        "aretes_ajoutees_p1": aretes_ajoutees_p1,
                        "aretes_ajoutees_p2": list(),
                        "aretes_supprimees_ps": aretes_ps,
                        "aretes_Ec_new": aretes_Ec_new,
                        "C_new": C_new,
                        "sommets_corriges": dico_sommets_corriges,
                        "sommets_non_corriges": dico_sommets_non_corriges,
                        "dico_sommets_par_cliqs_new": dico_sommets_par_cliqs_new
                        }
    elif dico_C1_C2_S1 and not dico_cliques_augmentante :
        logger.debug(" * compression_sommet : *** dico_C1_C2_S1 NON VIDE, "+\
                         "dico_cliques_augmentante VIDE ")
        for k_c1_c2_s1, val_cpt_c1_c2_s1 in dico_C1_C2_S1.items():
            cpt_prod_cartesien += 1;
            p1 = val_cpt_c1_c2_s1["clique_possible"];
            p2 = frozenset();
            gamma_z = args["dico_gamma_sommets"][sommet_z][1];
            ps = gamma_z - p1.intersection(gamma_z);
            aretes_ps = set( frozenset((sommet_z, sommet_ps)) 
                                    for sommet_ps in ps)
            aretes_p1 = fct_aux.aretes_dans_cliques(p1);
            aretes_ajoutees_p1 = aretes_differente(args["aretes_Ec"], 
                                                 aretes_p1);
            
            aretes_p2 = fct_aux.aretes_dans_cliques(p2);
            aretes_ajoutees_p2 = aretes_differente(args["aretes_Ec"], 
                                                 aretes_p2);
                                                
            aretes_Ec_new = set(args["aretes_Ec"]).union(
                            aretes_ajoutees_p1.union(aretes_ajoutees_p2));
            
            C_new = set(args["C"].copy());
            #TODO OK: a verifier C_new les cliques retirees
            ens_cliq_a_supprimer = set();
            ens_cliq_a_supprimer.add(val_cpt_c1_c2_s1["cliques_contratables"][0]);
            ens_cliq_a_supprimer.add(val_cpt_c1_c2_s1["cliques_contratables"][1]);
            for cliq in val_cpt_c1_c2_s1["cliques_suppl_contractables"] :
                ens_cliq_a_supprimer.add(cliq);
            C_new.difference_update(ens_cliq_a_supprimer);
            
            C_new.add( p1 );
            dico_sommets_corriges = dict();
            dico_sommets_non_corriges = dict();
            dico_sommets_par_cliqs_new = dict();
            C_new, aretes_Ec,\
            dico_sommets_corriges, \
            dico_sommets_non_corriges, \
            dico_sommets_par_cliqs_new = \
                mise_a_jour_aretes_cliques(C_new.copy(), aretes_Ec_new.copy(), 
                                    aretes_ps,\
                                    sommets_a_corriger.copy(), \
                                    args["dico_sommets_par_cliqs"].copy())
            dico_p1_p2_ps[cpt_prod_cartesien] = {
                        "id_sommet_1": id_sommet_z,
                        "sommet_1": sommet_z,
                        "p1": val_cpt_c1_c2_s1["clique_possible"],
                        "p2": frozenset(),
                        "ps": ps,
                        "voisine": frozenset(),
                        "dependante": frozenset(),
                        "contractable1": val_cpt_c1_c2_s1["cliques_contratables"][0],
                        "contractable2": val_cpt_c1_c2_s1["cliques_contratables"][1],
                        "S1": val_cpt_c1_c2_s1["S1"],
                        "S_z": s_z,
                        "aretes_ajoutees_p1": aretes_ajoutees_p1,
                        "aretes_ajoutees_p2": aretes_ajoutees_p2,
                        "aretes_supprimees_ps": aretes_ps,
                        "aretes_Ec_new": aretes_Ec_new,
                        "C_new": C_new,
                        "sommets_corriges": dico_sommets_corriges,
                        "sommets_non_corriges": dico_sommets_non_corriges,
                        "dico_sommets_par_cliqs_new": dico_sommets_par_cliqs_new
                        }
    else:
        logger.debug(" * compression_sommet : *** dico_C1_C2_S1, "+\
                         "dico_cliques_augmentante NON VIDE ")
        for k_c1_c2_s1, val_cpt_c1_c2_s1 in dico_C1_C2_S1.items():
            for k_cpt_vois_depend, val_cpt_vois_depend in dico_cliques_augmentante.items():                                        
                cpt_prod_cartesien += 1;
                
                inter_p1_p2 = val_cpt_c1_c2_s1["clique_possible"].intersection(
                                k_cpt_vois_depend[1].union(k_cpt_vois_depend[2])
                                )
                logger.debug(" * compression_sommet : *** "+\
                             "cpt_prod_cart={},".format(cpt_prod_cartesien)+\
                             "inter_p1_p2={},".format(len(inter_p1_p2)))
                print(" ***{} inter_p1_p2={},".format(cpt_prod_cartesien,inter_p1_p2)+\
                      "cliq_possible={}".format(val_cpt_c1_c2_s1["clique_possible"])+\
                      "vois_dep={}".format(k_cpt_vois_depend[1].union(k_cpt_vois_depend[2])
                      ))
                if len(inter_p1_p2) <= 1 and inter_p1_p2 == frozenset({sommet_z}):
                    
                    p1 = val_cpt_c1_c2_s1["clique_possible"];
                    p2 = val_cpt_vois_depend["voisine"].union(
                                val_cpt_vois_depend["dependante"].union(
                                frozenset({sommet_z})))
                    gamma_z = args["dico_gamma_sommets"][sommet_z][1];
                    ps = gamma_z - p1.intersection(gamma_z).union(
                                    p2.intersection(gamma_z));
                    aretes_ps = set( frozenset((sommet_z, sommet_ps)) 
                                            for sommet_ps in ps)
                    # TODO OK: calcul nombre aretes a ajouter pour p1, p2, nombre sommets a corriger
                    aretes_p1 = fct_aux.aretes_dans_cliques(p1);
                    aretes_ajoutees_p1 = aretes_differente(args["aretes_Ec"], 
                                                         aretes_p1);
                    
                    aretes_p2 = fct_aux.aretes_dans_cliques(p2);
                    aretes_ajoutees_p2 = aretes_differente(args["aretes_Ec"], 
                                                         aretes_p2);
                                                        
                    aretes_Ec_new = set(args["aretes_Ec"]).union(
                                    aretes_ajoutees_p1.union(aretes_ajoutees_p2));
                    
                    C_new = set(args["C"].copy());
                    #TODO OK: a verifier C_new les cliques retirees
                    ens_cliq_a_supprimer = set();
                    ens_cliq_a_supprimer.add(val_cpt_c1_c2_s1["cliques_contratables"][0]);
                    ens_cliq_a_supprimer.add(val_cpt_c1_c2_s1["cliques_contratables"][1]);
                    ens_cliq_a_supprimer.add(val_cpt_vois_depend["voisine"]);
                    ens_cliq_a_supprimer.add(val_cpt_vois_depend["dependante"]);
                    for cliq in [item for subitem in [val_cpt_c1_c2_s1["cliques_suppl_contractables"], 
                                 val_cpt_vois_depend["cliques_suppl_contractables"]] \
                                 for item in subitem]:
                        ens_cliq_a_supprimer.add(cliq);
                    C_new.difference_update(ens_cliq_a_supprimer);
                    
                    C_new.add( p1 );
                    C_new.add( p2 );
                    print("C_new={}".format(C_new))
                    dico_sommets_corriges = dict();
                    dico_sommets_non_corriges = dict();
                    dico_sommets_par_cliqs_new = dict();
                    C_new, aretes_Ec,\
                    dico_sommets_corriges, \
                    dico_sommets_non_corriges, \
                    dico_sommets_par_cliqs_new = \
                        mise_a_jour_aretes_cliques(C_new.copy(), aretes_Ec_new.copy(), 
                                            aretes_ps,\
                                            sommets_a_corriger.copy(), \
                                            args["dico_sommets_par_cliqs"].copy())
                    ##
    #                print("!!!!C_new={},\n".format(C_new)+ 
    #                      "dico_sommets_corriges={}, \n".format(dico_sommets_corriges)+
    #                      "dico_sommets_non_corriges={}, \n".format(dico_sommets_non_corriges)+
    #                      "dico_sommets_par_cliqs_new={},\n".format(dico_sommets_par_cliqs_new)+
    #                      " aretes_ps={}".format(aretes_ps))
                    ##
                    dico_p1_p2_ps[cpt_prod_cartesien] = {
                        "id_sommet_1": id_sommet_z,
                        "sommet_1": sommet_z,
                        "p1": val_cpt_c1_c2_s1["clique_possible"],
                        "p2": val_cpt_vois_depend["voisine"].union(
                                val_cpt_vois_depend["dependante"].union(
                                frozenset({sommet_z}))),
                        "ps": ps,
                        "voisine": val_cpt_vois_depend["voisine"],
                        "dependante": val_cpt_vois_depend["dependante"],
                        "contractable1": val_cpt_c1_c2_s1["cliques_contratables"][0],
                        "contractable2": val_cpt_c1_c2_s1["cliques_contratables"][1],
                        "S1": val_cpt_c1_c2_s1["S1"],
                        "S_z": s_z,
                        "aretes_ajoutees_p1": aretes_ajoutees_p1,
                        "aretes_ajoutees_p2": aretes_ajoutees_p2,
                        "aretes_supprimees_ps": aretes_ps,
                        "aretes_Ec_new": aretes_Ec_new,
                        "C_new": C_new,
                        "sommets_corriges": dico_sommets_corriges,
                        "sommets_non_corriges": dico_sommets_non_corriges,
                        "dico_sommets_par_cliqs_new": dico_sommets_par_cliqs_new
                        }
                if cpt_prod_cartesien >= nbre_elts_pi1_pi2:
                    break;
            if cpt_prod_cartesien >= nbre_elts_pi1_pi2:
                break;
#    print("@@@cpt_prod_cartesien={}, dico_C1_C2_S1={}, dico_cliques_augmentante={}".\
#          format(cpt_prod_cartesien, len(dico_C1_C2_S1), len(dico_cliques_augmentante)))   
    logger.debug(" * compression_sommet : *** fin compression_sommet, "+\
                 " dico_p1_p2_ps:{}, ".format(len(dico_p1_p2_ps)))
    return dico_p1_p2_ps;
         
################### fonctions de bases pour la correction ==> fin   ###########

################### critere selection  compression ==> debut ##################
def rechercher_min_max(liste_tuples, critere):
    """ retourne la tuple (min, max)
    """
    min_c1 = np.inf;
    #max_c2 = 0;
    if critere == "C1":
        return min(liste_tuples)
    elif critere == "C2":
        return max(liste_tuples)
    elif critere == "C2_C1":
        liste_intermediaires = [];
        min_c1 = min(liste_tuples)[0]
        for tuple_ in liste_tuples:
           if tuple_[0] == min_c1:
               liste_intermediaires.append(tuple_)
        return max(liste_intermediaires)
        
    
def critere_C2_C1_local(dico_compression, args) :
    """ selectionner le dico selon C2 puis C1 parmi les compressions possibles 
        sommet a corriger (sommet a -1)
    
    C2 : le maximum de sommets corriges
        * choisir le sommet a -1 qui corrige le max de sommets a -1 possibles
    C1 : le minimum d'aretes corriges
    
    dico_compression : dictionnaire contenant les compression (p1,p2,ps) du 
                        sommet sommet_z
    """
    
    max_c2 = 0;
    min_c1 = np.inf;
    dico_c1_c2 = dict();
    
    # definition de C2
#    if args["critere_selection_compression"] == "voisins_corriges":             # C2
#        for id_sommet_sommet_z, solutions in dico_compression.items():
#            for cpt_prod_cartesien, dico_p1_p2_ps in solutions.items():
#                if len(dico_p1_p2_ps["sommets_corriges"]) >= max_c2 :
#                    max_c2 = len(dico_p1_p2_ps["sommets_corriges"]);
#                    min_c1 = max_c2;
#                    if min_c1 not in dico_c1_c2:
#                        dico_c1_c2[min_c1] = [dico_p1_p2_ps];
#                    else:
#                        dico_c1_c2[min_c1].append(dico_p1_p2_ps);
# faux => suppose que chaque sommet a une seule solution PAS CORRECTE A EFFACER
    # definition de C2
    if args["critere_selection_compression"] == "voisins_corriges":             # C2
        for cpt_prod_cartesien, dico_p1_p2_ps in dico_compression.items():
            print("dico_p1_p2_ps={}".format(len(dico_p1_p2_ps)))
            if len(dico_p1_p2_ps["sommets_corriges"]) >= max_c2 :
                max_c2 = len(dico_p1_p2_ps["sommets_corriges"]);
                min_c1 = max_c2;
                if min_c1 not in dico_c1_c2:
                    dico_c1_c2[min_c1] = [dico_p1_p2_ps];
                else:
                    dico_c1_c2[min_c1].append(dico_p1_p2_ps);
    # definition de C1
    elif args["critere_selection_compression"] == "nombre_aretes_corrigees":     # C1
        for cpt_prod_cartesien, dico_p1_p2_ps in dico_compression.items():
            nbre_aretes_corriges = len(dico_p1_p2_ps["aretes_ajoutees_p1"]) + \
                                    len(dico_p1_p2_ps["aretes_ajoutees_p2"]) + \
                                    len(dico_p1_p2_ps["aretes_supprimees_ps"]);
            print("@@ min={}, nbre_aretes_corriges={}".format(min_c1, nbre_aretes_corriges))
            min_c1 = nbre_aretes_corriges if min_c1 >= nbre_aretes_corriges \
                                          else min_c1;
            if min_c1 not in dico_c1_c2:
                dico_c1_c2[min_c1] = [dico_p1_p2_ps];
            else:
                dico_c1_c2[min_c1].append(dico_p1_p2_ps);
        pass
    # definition de C2 puis de C1
    elif args["critere_selection_compression"] == "voisins_nombre_aretes_corrigees": # C2_C1
        for cpt_prod_cartesien, dico_p1_p2_ps in dico_compression.items():
            if len(dico_p1_p2_ps["sommets_corriges"]) >= max_c2 :
                max_c2 = len(dico_p1_p2_ps["sommets_corriges"]);
                nbre_aretes_corriges = len(dico_p1_p2_ps["aretes_ajoutees_p1"]) + \
                                    len(dico_p1_p2_ps["aretes_ajoutees_p2"]) + \
                                    len(dico_p1_p2_ps["aretes_supprimees_ps"]);
                min_c1 = nbre_aretes_corriges if min_c1 >= nbre_aretes_corriges \
                                                else min_c1;
                if min_c1 not in dico_c1_c2:
                    dico_c1_c2[min_c1] = [dico_p1_p2_ps];
                else:
                    dico_c1_c2[min_c1].append(dico_p1_p2_ps);
#    print("dico_c1_c2={}".format(dico_c1_c2))                
    if not dico_c1_c2:
        return min_c1, max_c2, dico_c1_c2;
    else:
        return min_c1, max_c2, dico_c1_c2[min_c1];
        
def critere_C2_C1_global(dico_compression, args) :
    """ recherche la compression optimale parmi tous les sommets a corriger 
        selon les criteres C1 et C2.
        
    C2 : le maximum de sommets corriges
        * choisir le sommet a -1 qui corrige le max de sommets a -1 possibles
    C1 : le minimum d'aretes corriges    
    """
    
    max_c2_global = 0;
    min_c1_global = np.inf;
    dico_c1_c2 = dict();
    
    # selection C2
    # je cherche le min local de c1 pour tous les sommets a corriger
    # parmi les min locaux, je cherche le max global de c2
    # une fois la liste des (min_global,max_global), je prends le 1er element.
    if args["critere_selection_compression"] == "voisins_corriges":             # C2
        for id_sommet_z, dicos_p1_p2_ps in dico_compression.items():
            # selection de dico selon C1
            min_c1_local = np.inf;
            max_c2_local = dicos_p1_p2_ps[1]
            for dico_p1_p2_ps in dicos_p1_p2_ps[2] :
                nbre_aretes_corriges = len(dico_p1_p2_ps["aretes_ajoutees_p1"]) + \
                                        len(dico_p1_p2_ps["aretes_ajoutees_p2"]) + \
                                        len(dico_p1_p2_ps["aretes_supprimees_ps"]);
                print("@@ min={}, nbre_aretes_corriges={}".format(min_c1_local, 
                      nbre_aretes_corriges))
                min_c1_local = nbre_aretes_corriges \
                                if min_c1_local >= nbre_aretes_corriges \
                                else min_c1_local;
                if (min_c1_local,max_c2_local) not in dico_c1_c2:
                    dico_c1_c2[(min_c1_local, max_c2_local)] = [dico_p1_p2_ps];
                else:
                    dico_c1_c2[(min_c1_local, max_c2_local)].append(dico_p1_p2_ps);
        # selection selon C2
        cle_min_max_c2 = rechercher_min_max(dico_c1_c2.keys(), "C2");
        min_c1_global = cle_min_max_c2[0];
        max_c2_global = cle_min_max_c2[1]
        print("")
        return min_c1_global, max_c2_global, dico_c1_c2[cle_min_max_c2][0];
    
    elif args["critere_selection_compression"] == "nombre_aretes_corrigees":     # C1
        for id_sommet_z, dicos_p1_p2_ps in dico_compression.items():
            # selection de dico selon C2
            max_c2_local = 0;
            min_c1_local = dicos_p1_p2_ps[0]
            for dico_p1_p2_ps in dicos_p1_p2_ps[2] :
                nbre_sommets_corriges = len(dico_p1_p2_ps["sommets_corriges"]);
                max_c2_local = nbre_sommets_corriges \
                                if nbre_sommets_corriges > max_c2_local \
                                else max_c2_local;
                if (min_c1_local,max_c2_local) not in dico_c1_c2:
                    dico_c1_c2[(min_c1_local, max_c2_local)] = [dico_p1_p2_ps];
                else:
                    dico_c1_c2[(min_c1_local, max_c2_local)].append(dico_p1_p2_ps);
        # selection selon C1
        cle_min_max_c2 = rechercher_min_max(dico_c1_c2.keys(), "C1");
        min_c1_global = cle_min_max_c2[0];
        max_c2_global = cle_min_max_c2[1]
        return min_c1_global, max_c2_global, dico_c1_c2[cle_min_max_c2][0];
        
    elif args["critere_selection_compression"] == "voisins_nombre_aretes_corrigees": # C2_C1
        for id_sommet_z, dicos_p1_p2_ps in dico_compression.items():
            min_c1_local = np.inf
            max_c2_local = 0
            for dico_p1_p2_ps in dicos_p1_p2_ps[2] :
                nbre_sommets_corriges = len(dico_p1_p2_ps["sommets_corriges"]);
                max_c2_local = nbre_sommets_corriges \
                                if nbre_sommets_corriges > max_c2_local \
                                else max_c2_local;
                nbre_aretes_corriges = len(dico_p1_p2_ps["aretes_ajoutees_p1"]) + \
                                        len(dico_p1_p2_ps["aretes_ajoutees_p2"]) + \
                                        len(dico_p1_p2_ps["aretes_supprimees_ps"]);
                min_c1_local = nbre_aretes_corriges \
                                if min_c1_local >= nbre_aretes_corriges \
                                else min_c1_local;
                if (min_c1_local,max_c2_local) not in dico_c1_c2:
                    dico_c1_c2[(min_c1_local, max_c2_local)] = [dico_p1_p2_ps];
                else:
                    dico_c1_c2[(min_c1_local, max_c2_local)].append(dico_p1_p2_ps);
                
        # selection selon C2_C1
        cle_min_max_c2 = rechercher_min_max(dico_c1_c2.keys(), "C2_C1");
        min_c1_global = cle_min_max_c2[0];
        max_c2_global = cle_min_max_c2[1]
        return min_c1_global, max_c2_global, dico_c1_c2[cle_min_max_c2][0];
################### critere selection  compression ==> fin ##################

################### application de la compression ==> debut ##################
def appliquer_correction(dico_sol_C2_C1, sommets_a_corriger, args):
    """ appliquer la compression choisie dans le graphe.
    """
    C = list();
    C = dico_sol_C2_C1["C_new"];
    aretes_Ec = list();
    aretes_Ec = dico_sol_C2_C1["aretes_Ec_new"];
    
    id_sommets_1 = list(dico_sol_C2_C1["sommets_corriges"].keys());
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
    logger = logging.getLogger('correction_graphe_correlation');
    dico_sommets_corriges = dict();
    dico_compression = dict();
    sommets_a_corriger = list();
    sommets_a_corriger = [sommet for sommet, etat in args["dico_cliq"].items() 
                            if etat == -1]
    if args["critere_selection_compression"] == "voisins_corriges" and \
        args["mode_correction"] == "aleatoire_sans_remise":
        # correction sans remise avec le critere "nombre de voisins corriges"
        logger.debug(" * critere_selection_compression : {}".\
                     format(args["critere_selection_compression"]))
        logger.debug(" * mode_correction : {}".\
                     format(args["mode_correction"]))
        cpt_noeud = 0;
        while(sommets_a_corriger):
            for id_sommet_1, sommet_1 in enumerate(sommets_a_corriger):
                cliques_sommet_1 = cliques_sommet(sommet_1, 
                                            args["dico_sommets_par_cliqs"]);
                logger.debug(" * sommet_1:{}".format(sommet_1) +\
                        ", cliques_sommets_1:{}".format(len(cliques_sommet_1)))
                #print("XX cliques_sommet_1={}".format(cliques_sommet_1))
#                dico_compression[(id_sommet_1,sommet_1)] = compression_sommet(
#                                                             id_sommet_1,
#                                                             sommet_1,
#                                                             sommets_a_corriger,
#                                                             cliques_sommet_1.copy(),
#                                                            args);
                dico_p1_p2_ps = dict();
                dico_p1_p2_ps = compression_sommet(id_sommet_1,
                                                     sommet_1,
                                                     sommets_a_corriger,
                                                     cliques_sommet_1.copy(),
                                                    args);
                
                dico_compression[(id_sommet_1,sommet_1)] = critere_C2_C1_local(
                                                            dico_p1_p2_ps,
                                                            args)
            # dico_sol_C2_C1 = {
            #        "id_sommet_1":,"sommet_1":,"p1":,"p2":,"ps":,
            #        "voisine":,"dependante":,
            #        "contractable1":,"contractable2":,
            #        "S1":,"S_z":,
            #        "aretes_ajoutees_p1":,"aretes_ajoutees_p2":,
            #        "aretes_supprimees_ps":,"aretes_Ec_new":,"C_new":,
            #        "sommets_corriges":,"sommets_non_corriges":,
            #        "dico_sommets_par_cliqs_new": 
            #        }
            
#            print("\n dico_compr={}".format(dico_compression))
            dico_sol_C2_C1 = dict();
            min_c1 = 0; max_c2 = 0;
            min_c1, max_c2, dico_sol_C2_C1 = critere_C2_C1_global(
                                            dico_compression,
                                            args)                               # C2 : nombre maximum de voisins corriges par un sommet, C1 : nombre minimum d'aretes a corriger au voisinage d'un sommet  
            
            logger.debug(" * correction_graphe : sommet_1:{}".format(sommet_1))
            logger.debug(" *** min_c1:{}".format(min_c1))
            logger.debug(" *** max_c2:{}".format(max_c2))
            logger.debug(" *** dico_sol_C2_C1:{}".format( len(dico_sol_C2_C1)))
            logger.debug(" *** C_old:{}, aretes_Ec_old:{}".format( len(args["C"]),
                         len(args["aretes_Ec"])))
            
            C, aretes_Ec, sommets_a_corriger = appliquer_correction(
                                                dico_sol_C2_C1,
                                                sommets_a_corriger,
                                                args)
            
            logger.debug(" *** C:{}, aretes_Ec:{}".format( len(C), 
                         len(aretes_Ec)))
            logger.debug(" *** p1:{}, p2:{}, ps:{}".format(dico_sol_C2_C1["p1"],
                         dico_sol_C2_C1["p2"], 
                         dico_sol_C2_C1["aretes_supprimees_ps"]))
            logger.debug(" *** sommet_corriges:{}, sommets_a_corriger={}".format(
                         dico_sol_C2_C1["sommets_corriges"],
                         sommets_a_corriger))
            for sommet, cliques in dico_sol_C2_C1[
                                    "dico_sommets_par_cliqs_new"].items():
                logger.debug(" *** sommet_par_cliques {}={}".format(
                             sommet, cliques));
            
            args["C"] = C;
            args["aretes_cliques"] = fct_aux.aretes_dans_cliques(C);
            args["aretes_Ec"] = aretes_Ec;
            cout_T = {"aretes_ajoutees_p1":dico_sol_C2_C1["aretes_ajoutees_p1"],
                      "aretes_ajoutees_p2":dico_sol_C2_C1["aretes_ajoutees_p2"],
                      "aretes_supprimees":dico_sol_C2_C1["aretes_supprimees_ps"],
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
        
    elif args["critere_selection_compression"] == "voisins_corriges" and \
        args["mode_correction"] == "aleatoire_avec_remise":
         pass 
    
    elif args["critere_selection_compression"] == "nombre_aretes_corrigees" and \
        args["mode_correction"] == "aleatoire_sans_remise":
        pass
    elif args["critere_selection_compression"] == "nombre_aretes_corrigees" and \
        args["mode_correction"] == "aleatoire_avec_remise":
        pass
    elif args["critere_selection_compression"] == "voisins_nombre_aretes_corrigees" and \
        args["mode_correction"] == "aleatoire_sans_remise":
        pass
    elif args["critere_selection_compression"] == "voisins_nombre_aretes_corrigees" and \
        args["mode_correction"] == "aleatoire_avec_remise":
        pass
    
    # degre min {avec, sans} remise 
    # TODO ecrire cette condition
    elif args["critere_selection_compression"] == "voisins_corriges" and \
        args["mode_correction"] == "degre_min_sans_remise":
         pass 
    elif args["critere_selection_compression"] == "voisins_corriges" and \
        args["mode_correction"] == "degre_min_avec_remise":
         pass 
    elif args["critere_selection_compression"] == "nombre_aretes_corrigees" and \
        args["mode_correction"] == "degre_min_sans_remise":
        pass
    elif args["critere_selection_compression"] == "nombre_aretes_corrigees" and \
        args["mode_correction"] == "degre_min_avec_remise":
        pass
    elif args["critere_selection_compression"] == "voisins_nombre_aretes_corrigees" and \
        args["mode_correction"] == "degre_min_sans_remise":
        pass
    elif args["critere_selection_compression"] == "voisins_nombre_aretes_corrigees" and \
        args["mode_correction"] == "degre_min_avec_remise":
        pass   
    
    # cout min {avec, sans} remise
    # TODO ecrire cette condition
    elif args["critere_selection_compression"] == "voisins_corriges" and \
        args["mode_correction"] == "cout_min_sans_remise":
         pass 
    elif args["critere_selection_compression"] == "voisins_corriges" and \
        args["mode_correction"] == "cout_min_avec_remise":
         pass 
    elif args["critere_selection_compression"] == "nombre_aretes_corrigees" and \
        args["mode_correction"] == "cout_min_sans_remise":
        pass
    elif args["critere_selection_compression"] == "nombre_aretes_corrigees" and \
        args["mode_correction"] == "cout_min_avec_remise":
        pass
    elif args["critere_selection_compression"] == "voisins_nombre_aretes_corrigees" and \
        args["mode_correction"] == "cout_min_sans_remise":
        pass
    elif args["critere_selection_compression"] == "voisins_nombre_aretes_corrigees" and \
        args["mode_correction"] == "cout_min_avec_remise":
        pass  
    
    else:
        print("ERROR critere_selection_compression={},mode_correction={} dont exist".\
              format(args["critere_selection_compression"], args["mode_correction"]))
    
if __name__ == '__main__':
    ti = time.time();
    
    C1 = {"1","2","11"}; C2 = {"z","2","3"}; C3 = {"z","4","5"}; 
    C4 = {"z","6","7","8"}; C5 = {"8","9","10"}; C6 = {"z","10"}; C7 = {"z","1"}
    C = [frozenset(C1), frozenset(C2), frozenset(C3), frozenset(C4),
         frozenset(C5), frozenset(C6), frozenset(C7)]
    dico_cliq = {"z":-1, "1":1, "2":1, "3":1, "4":1, "5":1, "6":1,
                 "7":1, "8":1, "9":1,"10":1, "11":1}
    
    dico_matE = {"z":["1","2","3","4","5","6","7","8","10"], "1":["11","2","z"], # graphe du manuscrit
                 "2":["1","z","11","3"], "3":["z","2"], "4":["z","5"], 
                 "5":["z","4"], "6":["z","7","8"], "7":["z","8","6"], 
                 "8":["z","10","9","7","6"], "9":["10","8"], "10":["z","9","8"],
                 "11":["1","2"]}
    matE_k_alpha = pd.DataFrame(index = dico_matE.keys(), 
                                columns = dico_matE.keys());
    for k, vals in dico_matE.items():
        for v in vals:
            matE_k_alpha.loc[k, v] = 1;
            matE_k_alpha.loc[v, k] = 1;
    matE_k_alpha.fillna(value=0, inplace=True);
    
    aretes_Ec = fct_aux.liste_arcs(matE_k_alpha);
    dico_sommets_par_cliqs = fct_aux.couverture_par_sommets(C);
    dico_gamma_sommets = fct_aux.gamma_noeud(matE_k_alpha, aretes_Ec)
    
    sommet_z = "z";
    sommets_a_corriger = [sommet_z];
    
    cliques_sommet_z = cliques_sommet(sommet_z, 
                                      dico_sommets_par_cliqs)
    number_items_pi1_pi2 = 1;
    args = dict();
    args["C"] = C;
    args["dico_sommets_par_cliqs"] = dico_sommets_par_cliqs;
    args["dico_cliq"] = dico_cliq;
    args["aretes_Ec"] = fct_aux.liste_arcs(matE_k_alpha);
    args["aretes_cliques"] = fct_aux.aretes_dans_cliques(C);
    args["dico_gamma_sommets"] = dico_gamma_sommets;
    args["number_items_pi1_pi2"] = number_items_pi1_pi2;
    
    
    
    test_cliques_sommet = False;
    test_S_sommet = False;
    test_is_contractable = False#True#False;
    test_cliques_contractables = False;
    test_voisine_sommet = False;
    test_dependance_sommet = False#True;
    test_augmentation = False#True;
    test_compression_sommet = False#True #False;
    test_critere_C2_C1 = False#True#False;
    test_mise_a_jour_cliques = False;                                           # pas defini
    test_appliquer_correction = False #True#False;                                    
    
    
    if test_cliques_sommet :
        cliques = cliques_sommet(sommet_z, args["dico_sommets_par_cliqs"])
        if C3 in args["C"] and C2 in args["C"] and C4 in args["C"] :
            print("test_cliques_sommet : *** OK ***")
        else:
            for c in [C2, C3, C4]:
                if c not in args["C"]:
                   print("test_cliques_sommet : ***ERROR: {} not in C ***".\
                         format(C)) 
    
    if test_S_sommet:
        s_z = S_sommet(sommet_z, args["dico_gamma_sommets"][sommet_z][1], 
                       aretes_Ec, C, args["aretes_cliques"]);
        if not s_z :
            print("test_S_sommet : ***ERROR: s_z = {}  ***");
        elif "1" in s_z and "10" in s_z:
            print("test_S_sommet : ***OK: s_z={} ***".\
                         format(s_z))
        else:
            for sommet in s_z:
                if sommet != "1" and sommet != "10":
                    print("test_S_sommet : ***ERROR: sommet={} not in  {1,10}  ***".\
                          format(sommet))
                    
    if test_is_contractable:
        for c1, c2 in [(C2,C3), (C3,C4), (C4,C6), (C1,C2), (C4,C5)]:
            boolean_contractable, cliques_suppl_contractables = \
                is_contractable(c1, c2, args["aretes_Ec"], \
                               args["aretes_cliques"], args["C"])
            if boolean_contractable:
                print("test_is_contractable:*** {},{},{} contractables".\
                      format(c1,c2, cliques_suppl_contractables))
            else:
                print("test_is_contractable:*** {},{} non contractables".\
                      format(c1,c2))
            
    if test_cliques_contractables :
        cliq_contractables = cliques_contractables(sommet_z, args["aretes_Ec"], 
                                                   args["aretes_cliques"], 
                                                   cliques_sommet_z, 
                                                   args["C"])
        print("test_cliques_contractables: ***cliq_contractables {}".\
                      format(cliq_contractables))
        cliq_contr_connu = [(C2,C3),(C3,C4)]
        cliq_non_contr_connu = [(C4,C6),(C1,C2)]
        if not cliq_contractables:
            print("test_cliques_contractables: ***ERROR: ens. vide")
        for c1_c2 in cliq_contr_connu:
            if (c1_c2[0],c1_c2[1]) in cliq_contractables or \
                (c1_c2[1],c1_c2[0]) in cliq_contractables:
                print("test_cliques_contractables: ***contraction {}: OK".\
                      format(c1_c2))
            else:
                print("test_cliques_contractables: ***contraction {}: NOK".\
                      format(c1_c2))
                
        for c1_c2 in cliq_non_contr_connu:
            if (c1_c2[0],c1_c2[1]) not in cliq_contractables and \
                (c1_c2[1],c1_c2[0]) not in cliq_contractables:
                print("test_cliques_contractables: ***non_contraction {}: OK".\
                      format(c1_c2))
            else:
                print("test_cliques_contractables: ***non_contraction {}: NOK".\
                      format(c1_c2))
                
                
    if test_voisine_sommet:
        s_z = S_sommet(sommet_z, args["dico_gamma_sommets"][sommet_z][1], 
                       aretes_Ec, args["C"], args["aretes_cliques"]);
        cliques_not_sommet_z = [ c for c in args["C"] 
                                 if not c.intersection({sommet_z})]
        cliques_voisines= voisine_sommet(sommet_z, \
                                         cliques_sommet_z, \
                                         cliques_not_sommet_z, s_z)
        cliq_vois_connus = [C1, C5];
        print("test_voisine_sommet: *** cliques_voisines={}".\
              format(cliques_voisines))
        for c in cliq_vois_connus:
            if c not in cliques_voisines:
                print("test_voisine_sommet: ***ERROR : {} not in cliques_voisines".\
                      format(c))
        for c in cliques_voisines:
            if c not in cliq_vois_connus:
                print("test_voisine_sommet: ***ERROR : ATTENTION {}  inconnu theoriquement".\
                      format(c))
                
                
    if test_dependance_sommet :
        s_z = S_sommet(sommet_z, args["dico_gamma_sommets"][sommet_z][1], 
                       aretes_Ec, args["C"], args["aretes_cliques"]);
        cliques_not_sommet_z = [ c for c in args["C"] 
                                 if not c.intersection({sommet_z})]
        cliques_voisines= voisine_sommet(sommet_z, 
                                         cliques_sommet_z, 
                                         cliques_not_sommet_z, s_z)
        for clique_voisine in cliques_voisines:
            depends_sommet = dependance_sommet(sommet_z, 
                                args["dico_gamma_sommets"][sommet_z][1], 
                                cliques_sommet_z, clique_voisine)
            if not depends_sommet:
                print("test_dependance_sommet: ***ERROR : ens. vide")
            else:
                print("test_dependance_sommet: *** cliq_vois ={}, dependances={}".\
                      format(clique_voisine, depends_sommet))
                
                
    if test_augmentation :
        s_z = S_sommet(sommet_z, args["dico_gamma_sommets"][sommet_z][1], 
                       aretes_Ec, args["C"], args["aretes_cliques"]);
        dico_cliques_augmentante = augmentation(sommet_z, 
                                    args["dico_gamma_sommets"][sommet_z][1], 
                                    cliques_sommet_z, s_z, 
                                    args)
        if not dico_cliques_augmentante :
            print("test_augmentation: ***ERROR : dico. vide")
        else:
            for _ in range(10):
                print("test_augmentation: ***")
                for key, val in dico_cliques_augmentante.items():
                    print("test_augmentation: *** cpt:{}, sommet:{}, voisine:{}, dependante:{}, cliques_suppl_contractables:{}".\
                          format(key[0], val["sommet_z"], val["voisine"], 
                                 val["dependante"], 
                                 val["cliques_suppl_contractables"]))
      
    if test_compression_sommet :
        id_sommet_z = 0;
        dico_p1_p2_ps = compression_sommet(id_sommet_z, sommet_z, 
                                sommets_a_corriger, 
                                cliques_sommet_z, args)
        if not dico_p1_p2_ps :
            print("test_compression_sommet: ***ERROR : dico. vide")
        else:
            for _ in range(3):
                print("test_compression_sommet: ***")
                for key, val in dico_p1_p2_ps.items():
                    print("test_compression: *** p1:{},\n".format(val["p1"])+
                          "p2:{},\n".format(val["p2"])+ 
                          "ps:{},\n".format(val["ps"])+ 
                          "voisine:{},\n".format(val["voisine"])+  
                    "dependante:{},\n".format(val["dependante"])+ 
                    "contractable1:{},\n".format(val["contractable1"])+ 
                    "contractable2:{},\n".format(val["contractable2"])+
                    "S1:{},\n".format(val["S1"])+
                    "S_z:{},\n".format(val["S_z"])+
                    "aretes_ajoutees_p1:{},\n".format(val["aretes_ajoutees_p1"])+
                    "aretes_ajoutees_p2:{},\n".format(val["aretes_ajoutees_p2"])+
                    "aretes_supprimees_ps:{},\n".format(val["aretes_supprimees_ps"])+
                    "sommets_corriges:{},\n".format(val["sommets_corriges"])+
                    "sommets_non_corriges:{},\n".format(val["sommets_non_corriges"])+
                    "C_new:{},\n".format(val["C_new"])+
                    "dico_sommets_par_cliqs_new:{}\n".format(val["dico_sommets_par_cliqs_new"])
                    )
                          
    if test_critere_C2_C1 :
        id_sommet_z = 0;
        critere_selection_compression = "voisins_corriges"; # voisins_corriges => TESTER
                                                        #"nombre_aretes_corrigees" => TESTER
                                                        #"voisins_nombre_aretes_corrigees" => TESTER
        args["critere_selection_compression"] = critere_selection_compression;
        
        dico_compression = compression_sommet(id_sommet_z, sommet_z, 
                                sommets_a_corriger, 
                                cliques_sommet_z, args)
        min_c1, max_c2, dico_sol_C2_C1 = critere_C2_C1(dico_compression,args) 
        if not dico_sol_C2_C1:
            print("test_critere_C2_C1: ***ERROR dico. vide");
        else:
            print("test_critere_C2_C1: *** min_c1={}".format(min_c1))
            print("test_critere_C2_C1: *** max_c2={}".format(max_c2))
            print("test_critere_C2_C1: *** type={}".format(type(dico_sol_C2_C1)))
            for dico_sol_c2_c1 in dico_sol_C2_C1 :
                print("test_critere_c1_c2: *** p1:{},\n".format(dico_sol_c2_c1["p1"]) +\
                          "p2:{},\n".format(dico_sol_c2_c1["p2"])+ 
                          "ps:{},\n".format(dico_sol_c2_c1["ps"])+ 
                          "voisine:{},\n".format(dico_sol_c2_c1["voisine"])+  
                    "dependante:{},\n".format(dico_sol_c2_c1["dependante"])+ 
                    "contractable1:{},\n".format(dico_sol_c2_c1["contractable1"])+ 
                    "contractable2:{},\n".format(dico_sol_c2_c1["contractable2"])+
                    "S1:{},\n".format(dico_sol_c2_c1["S1"])+
                    "S_z:{},\n".format(dico_sol_c2_c1["S_z"])+
                    "aretes_ajoutees_p1:{},\n".format(dico_sol_c2_c1["aretes_ajoutees_p1"])+
                    "aretes_ajoutees_p2:{},\n".format(dico_sol_c2_c1["aretes_ajoutees_p2"])+
                    "aretes_supprimees_ps:{},\n".format(dico_sol_c2_c1["aretes_supprimees_ps"])+
                    "sommets_corriges:{},\n".format(dico_sol_c2_c1["sommets_corriges"])+
                    "sommets_non_corriges:{},\n".format(dico_sol_c2_c1["sommets_non_corriges"])+
                    "C_new:{},\n".format(dico_sol_c2_c1["C_new"])+
                    "dico_sommets_par_cliqs_new:{}\n".format(dico_sol_c2_c1["dico_sommets_par_cliqs_new"])
                    )
                    
    if test_appliquer_correction:
        id_sommet_z = 0;
        critere_selection_compression = "voisins_corriges"; 
        args["critere_selection_compression"] = critere_selection_compression;
        
        dico_compression = compression_sommet(id_sommet_z, sommet_z, 
                                sommets_a_corriger, 
                                cliques_sommet_z, args)
        min_c1, max_c2, dico_sol_C2_C1 = critere_C2_C1(dico_compression, args) # C2 : nombre maximum de voisins corriges par un sommet, C1 : nombre minimum d'aretes a corriger au voisinage d'un sommet  
        print("test_appliquer_correction:***Avant \n"+
              "C={}".format(args["C"])+
              "\n aretes_Ec={}".format(len(aretes_Ec))+
              "\n sommets_a_corriger={}".format(sommets_a_corriger))
        C, aretes_Ec, sommets_a_corriger = appliquer_correction(
                                            dico_sol_C2_C1[0],
                                            sommets_a_corriger,
                                            args)
        print("test_appliquer_correction:***Apres "+
              "\n C={}".format(C)+
              "\n aretes_Ec={}".format(len(aretes_Ec))+
              "\n sommets_a_corriger={}".format(sommets_a_corriger))