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
        elif (voisin_z, sommet_z) not in aretes_cliques and \
            (sommet_z, voisin_z) not in aretes_cliques:
            S_z.append(voisin_z);
    return S_z;
    
def cliques_contractables(sommet_z, aretes_Ec, aretes_cliques, cliques_sommet_z):
    """ retourne la liste des cliques contractables autour du sommet_z. """

    cliq_contractables = [];
    for c1, c2 in it.combinations(cliques_sommet_z.append({})):
        if not c1 or not c2 :
            cliq_contractables.append((c1, c2));
        else:
            boolean_constractable = True;
            for noeud1, noeud2 in it.product(c1, c2):
                if ((noeud1, noeud2) in aretes_Ec or \
                    (noeud2, noeud1) in aretes_Ec) and \
                    {noeud1, noeud2} not in aretes_cliques :
                        boolean_constractable = False;
                        break;
            if boolean_constractable : 
               cliq_contractables.append((c1, c2)); 
            
    return cliq_contractables;
    
def voisine_sommet():
    pass

def dependance_sommet():
    pass
        
################### fonctions de bases pour la correction ==> fin   ###########


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
            cpt_noeud += 1;
            for id_sommet_1, sommet_1 in enumerate(sommets_a_corriger):
                dico_compression[(id_sommet_1,sommet_1)] = compression_sommet(
                                                            args, 
                                                            sommets_a_corriger);
            
            # dico_sol_C2_C1 = {id_sommet_1:, sommet_1:, compression:(pi1,p2,ps), voisins_corriges:{"id_voisin_ds_sommets_a_corriger":voisin}, cout_T:{"":[],"":[]}} 
            dico_sol_C2_C1 = dict();
            dico_sol_C2_C1 = critere_C2_C1(dico_compression)                   # C2 : nombre maximum de voisins corriges par un sommet, C1 : nombre minimum d'aretes a corriger au voisinage d'un sommet  
            C, aretes_Ec = appliquer_correction(dico_sol_C2_C1["sommet_1"], 
                                                dico_sol_C2_C1["compression"], 
                                                args)
            args["C"] = C;
            args["aretes_matE0"] = aretes_Ec;
            dico_sommets_corriges[(cpt_noeud, dico_sol_C2_C1["sommet_1"])] = {
                        "compression_p1":dico_sol_C2_C1["compression"][0],
                        "compression_p2":dico_sol_C2_C1["compression"][1],
                        "compression_ps":dico_sol_C2_C1["compression"][2],
                        "voisins_corriges":dico_sol_C2_C1["voisins_corriges"], # voisins_corriges = {"id_voisin_ds_sommets_a_corriger":voisin}
                        "cout_T":dico_sol_C2_C1["cout_T"],                     # cout_T={"aretes_ajoutes":[],"aretes_supprimees":[]}
                                                }
            sommets_a_corriger.pop(dico_sol_C2_C1["id_sommet_1"])
        # mettre a jour les cliques couvrants les sommets.
        args["dico_sommets_par_cliqs"] = fct_aux.couverture_par_sommets(
                                            args["C"]);
        return args, dico_sommets_corriges;
        
if __name__ == '__main__':
    ti = time.time();