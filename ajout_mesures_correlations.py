#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 09:53:31 2017

@author: willy
completer/ajouter les mesures et correlations sur les arcs n'ayant pas de mesures 
"""
import math;
import os, re;
import numpy as np;
import pandas as pd;
import fonctions_auxiliaires as fct_aux;
import itertools as it;

def comparateur_correlation(matE_gr, matE_gr_cp, cpt_modif):
    """
    retourner le nombre de correlations modifies
    """
    if set(matE_gr.columns.tolist()) != set(matE_gr_cp.columns.tolist()):
        print("matE_gr={} et matE_gr_cp={} sont de dimension differentes".format(\
              len(set(matE_gr.columns.tolist())), len(set(matE_gr_cp.columns.tolist()))))
        return np.nan, []
    cpt_diff = 0; correls_diff = []
    cpt_ident = 0; correls_ident = []
    cpt_case =0
    for row, col in fct_aux.range_2d(matE_gr.columns):
        if matE_gr.loc[row,col] != matE_gr_cp.loc[row,col]:
            cpt_diff += 1; correls_diff.append((row,col)); #print("ICI21 cpt_diff={}".format(cpt_diff))
        elif matE_gr.loc[row,col] == matE_gr_cp.loc[row,col]:
            cpt_ident += 1; correls_ident.append((row,col)); #print("ICI22 cpt_ident={}".format(cpt_ident))
        cpt_case += 1;
    print("%correl_modif={}, %correl_ident={}, modif_faites={}".format(round(cpt_diff/cpt_case,3), \
          round(cpt_ident/cpt_case,3), round(cpt_modif/cpt_case,3)))
    dico_diff = dict()
    for corr in correls_diff:
        dico_diff[corr] = [ matE_gr_cp.loc[corr[0], corr[1]],round(matE_gr.loc[corr[0], corr[1]],4) ]
#    print("%aretes_modif={}".format(correls_diff))
    print("aretes_modif={}".format(dico_diff))
        
def gamma(reseauG,bool_entr_sort, df_gr, matE_gr):
    """
    recherche les aretes sortant et entrants de tous les sommets du graphe/sous graphe reel
    bool_entr_sort = True si entrant gamma- = entrant
                   = False si sortant gamma+ = sortant
                        cols
                     | y | t |          x---> y
                   -----------            \__t
              rows x | 1   1             z__/
                   z |     1
    """
    ###
#    print("len reseauG={}, len matE_gr={}".format( len(reseauG.columns.tolist()), len(matE_gr.columns.tolist())))
#    cols_r = reseauG.columns.tolist(); cols_gr = matE_gr.columns.tolist(); 
#    print("col_r = {}, col_gr = {}".format(cols_r[:3], cols_gr[:3]))
#    arcs_diff = set(cols_r).union(set(cols_gr)) - set(cols_r).intersection(set(cols_gr))  
#    print("arc_diff = {}, union={}, inter={}".format(arcs_diff, len(set(cols_r).union(set(cols_gr))), \
#          len(set(cols_r).intersection(set(cols_gr)))))
#    if set(reseauG.columns.tolist()).issubset( set(matE_gr.columns.tolist()) ):
#        print("ICI1 = len reseauG={}, len matE_gr={}".format( len(reseauG.columns.tolist()), len(matE_gr.columns.tolist())))
#    elif set(matE_gr.columns.tolist()).issubset( set(reseauG.columns.tolist()) ):
#        print("ICI2")
#    else:
#        print("ICI3")
    ###
    
    dico, dico_bar = dict(), dict();
    if bool_entr_sort == False:
        # m++ ==> arcs sortants
        for row in reseauG.index.tolist():
            aretes = list(); aretes_bar = list(); m_plus2 = 0
            for col in reseauG.columns.tolist():
                if reseauG.loc[row, col] == 1:
                    aretes.append( row+"->"+col )
                    if row+"->"+col in df_gr.columns and df_gr.loc[:,row+"->"+col].isnull().all():
                        aretes_bar.append(row+"->"+col);
                    elif col+"->"+row in df_gr.columns and df_gr.loc[:,col+"->"+row].isnull().all():
                        aretes_bar.append(row+"->"+col);
                    
            if len(aretes) == 0:
                dico[row] = [0, aretes]; 
                dico_bar[row] = [0, aretes_bar];
#                print("row={}".format(row))
            else:
                for arete1, arete2 in fct_aux.range_2d(aretes):
                    if arete1 in matE_gr.columns.tolist() and arete2 in matE_gr.columns.tolist():
                        m_plus2 += matE_gr.loc[arete1, arete2];                
                dico[row] = [m_plus2/len(aretes), aretes]
                dico_bar[row] = [0, aretes_bar];
    else:
        # m-- ==> arcs entrants
        for col in reseauG.columns.tolist():
            aretes = list(); aretes_bar = list(); m_minus2 = 0;
            for row in reseauG.index.tolist():
                if reseauG.loc[row, col] == 1:
                    aretes.append( row+"->"+col );
                    if row+"->"+col in df_gr.columns and df_gr.loc[:,row+"->"+col].isnull().all():
                        aretes_bar.append(row+"->"+col);
                    elif col+"->"+row in df_gr.columns and df_gr.loc[:,col+"->"+row].isnull().all():
                        aretes_bar.append(row+"->"+col);
#                    if df_gr.loc[:,row+"->"+col].isnull().all() or df_gr.loc[:,col+"->"+row].isnull().all():
#                        aretes_bar.append(row+"->"+col);
            if len(aretes) == 0:
                dico[col] = [0, aretes]; 
                dico_bar[col] = [0, aretes_bar];
#                print("col={}".format(col))
            else:
                for arete1, arete2 in fct_aux.range_2d(aretes):
                    if arete1 in matE_gr.columns.tolist() and arete2 in matE_gr.columns.tolist():
                        m_minus2 += matE_gr.loc[arete1, arete2];                
                dico[col] = [m_minus2/len(aretes), aretes]
                dico_bar[col] = [0, aretes_bar];
    return dico, dico_bar; 

def sommer_correl_sortantes_or_entrantes(aretes_plus, matE_gr):
    """
    faire la somme des correlations des aretes entrants/sortantes
    """
    m_plus = 0;cpt = 0
    for arete1_, arete2_ in fct_aux.range_2d(aretes_plus):
        arete1, arete2 = verifier_row_col(arete1_, arete2_, matE_gr);
        if arete1 != None and arete2 != None:
            cpt+=1;
            m_plus += matE_gr.loc[arete1, arete2];
        else:
            print("m_plus/minus: {} et {} dont belong to matE_gr".format(arete1_,arete2_) )
    if cpt== 0:
        return 0;
    return m_plus/cpt;

def verifier_row_col(arete1, arete2, matE_gr):
    """
    retourne les bonnes aretes  appartenant a matE_gr
    """
    if arete1.split("->")[0]+"->"+arete1.split("->")[1] in matE_gr.columns or \
         arete2.split("->")[0]+"->"+arete2.split("->")[1] in matE_gr.columns:
        return arete1.split("->")[0]+"->"+arete1.split("->")[1], \
                arete2.split("->")[0]+"->"+arete2.split("->")[1];
    elif arete1.split("->")[0]+"->"+arete1.split("->")[1] in matE_gr.columns or \
            arete2.split("->")[1]+"->"+arete2.split("->")[0] in matE_gr.columns:
        return arete1.split("->")[0]+"->"+arete1.split("->")[1], \
                arete2.split("->")[1]+"->"+arete2.split("->")[0];
    elif arete1.split("->")[1]+"->"+arete1.split("->")[0] in matE_gr.columns or \
            arete2.split("->")[0]+"->"+arete2.split("->")[1] in matE_gr.columns:
        return arete1.split("->")[1]+"->"+arete1.split("->")[0], \
                arete2.split("->")[0]+"->"+arete2.split("->")[1];
    elif arete1.split("->")[1]+"->"+arete1.split("->")[0] in matE_gr.columns or \
            arete2.split("->")[1]+"->"+arete2.split("->")[0] in matE_gr.columns:
        return arete1.split("->")[1]+"->"+arete1.split("->")[0], \
                arete2.split("->")[1]+"->"+arete2.split("->")[0];
    return None, None;
    pass

def sommer_correl_sortantes_entrantes(aretes_plus, aretes_minus, matE_gr):
    """
    sommer les correlations entre aretes sortantes et aretes entrantes pour un sommet defini
    """
    m_plus_minus = 0; cpt = 0;
    for tuple_ in it.product(aretes_plus, aretes_minus):
        arete1, arete2 = verifier_row_col(tuple_[0], tuple_[1], matE_gr);
        if arete1 != None and arete2 != None:
            m_plus_minus += matE_gr.loc[arete1, arete2];
            cpt+=1;
        else:
            print("m_plus_minus: {} et {} dont belong to matE_gr".format(tuple_[0],tuple_[1]) )
    if cpt == 0:
        return 0;
    return m_plus_minus/cpt;
    
def gamma_sommet_bis(sommet, reseauG, df_gr, matE_gr):
    """
    donner pour un noeud son degre entrant et aretes entrantes et aussi 
    son degre sortant et aretes sortantes
    """
    if sommet not in reseauG.columns:
        return [0,[]]
    if sommet not in reseauG.index:
        return [0,[]]
    
    #degre sortant d+,m+ ==> il correspond au parcours des index (row)
    m_plus2 = 0;aretes_plus=[]; aretes_plus_bar=[]; 
    for row in reseauG.index:
        if reseauG.loc[sommet, row] == 1:
            aretes_plus.append((sommet+"->"+row));
            if sommet+"->"+row in df_gr.columns and df_gr.loc[:,sommet+"->"+row].isnull().all():
                aretes_plus_bar.append(sommet+"->"+row);
            elif row+"->"+sommet in df_gr.columns and df_gr.loc[:,row+"->"+sommet].isnull().all():
                aretes_plus_bar.append(sommet+"->"+row);
    if len(aretes_plus) < 2 and len(aretes_plus_bar) < 1:
        m_plus2 = 0; m_plus_minus = 0;
    else:
        m_plus2 = sommer_correl_sortantes_or_entrantes(aretes_plus, matE_gr)
                
    #degre entrant d-,m- ==> il correspond au parcours des colonnes (col)
    m_minus2 = 0; aretes_minus=[]; aretes_minus_bar=[]; 
    for col in reseauG.columns:
        if reseauG.loc[sommet, col] == 1:
            aretes_minus.append((sommet+"->"+col));
            if sommet+"->"+col in df_gr.columns and df_gr.loc[:,sommet+"->"+col].isnull().all():
                aretes_minus_bar.append(sommet+"->"+col);
            elif col+"->"+sommet in df_gr.columns and df_gr.loc[:,col+"->"+sommet].isnull().all():
                aretes_minus_bar.append(sommet+"->"+col);
    if len(aretes_minus) < 2 and len(aretes_minus_bar) < 1:
        m_minus2 = 0; m_plus_minus = 0;
    else:
        m_minus2 = sommer_correl_sortantes_or_entrantes(aretes_minus, matE_gr);
        
    # calcul m_plus_minus:
    m_plus_minus = 0; 
    m_plus_minus = sommer_correl_sortantes_entrantes(aretes_plus, aretes_minus, matE_gr)
    
    return {"m_plus2":m_plus2, "aretes_plus":aretes_plus, "m_minus2":m_minus2,\
            "aretes_minus":aretes_minus, "m_plus_minus":m_plus_minus, \
            "aretes_plus_bar":aretes_plus_bar, "aretes_minus_bar":aretes_minus_bar}
            
def arete_existante(sommet,col,matE_gr):
    """
    retourne le nom de l arete contenu dans matE_gr
    """
    if sommet+"->"+col in matE_gr.columns:
        return sommet+"->"+col;
    elif col+"->"+sommet in matE_gr.columns:
        return col+"->"+sommet;
    else:
        return None;
    
def gamma_sommet(sommet, reseauG, df_gr, matE_gr):
    """
    donner pour un noeud son degre entrant et aretes entrantes et aussi 
    son degre sortant et aretes sortantes
    """
    if sommet not in reseauG.columns:
        return [0,[]]
    if sommet not in reseauG.index:
        return [0,[]]
    
    #degre sortant d+,m+ ==> il correspond au parcours des cols (cols)
    m_plus2 = 0;aretes_plus=[]; aretes_plus_bar=[]; 
    for col in reseauG.columns:
        if reseauG.loc[sommet, col] == 1:
            arete = arete_existante(sommet,col,matE_gr)
            aretes_plus.append(arete) if arete != None else None;
            if arete != None and arete in df_gr.columns and df_gr.loc[:,arete].isnull().all():
                aretes_plus_bar.append(arete);
    if len(aretes_plus) < 2 and len(aretes_plus_bar) < 1:
        m_plus2 = 0; m_plus_minus = 0;
    else:
        m_plus2 = sommer_correl_sortantes_or_entrantes(aretes_plus, matE_gr)
                
    #degre entrant d-,m- ==> il correspond au parcours des index (row)
    m_minus2 = 0; aretes_minus=[]; aretes_minus_bar=[]; 
    for row in reseauG.index:
        if reseauG.loc[row,sommet] == 1:
            arete = arete_existante(sommet,row,matE_gr)
            aretes_minus.append(arete) if arete != None else None;
            if arete != None and arete in df_gr.columns and df_gr.loc[:, arete].isnull().all():
                aretes_minus_bar.append(arete);
    if len(aretes_minus) < 2 and len(aretes_minus_bar) < 1:
        m_minus2 = 0; m_plus_minus = 0;
    else:
        m_minus2 = sommer_correl_sortantes_or_entrantes(aretes_minus, matE_gr);
        
    # calcul m_plus_minus:
    m_plus_minus = 0; 
    m_plus_minus = sommer_correl_sortantes_entrantes(aretes_plus, aretes_minus, matE_gr)
    
    return {"m_plus2":m_plus2, "aretes_plus":aretes_plus, "m_minus2":m_minus2,\
            "aretes_minus":aretes_minus, "m_plus_minus":m_plus_minus, \
            "aretes_plus_bar":aretes_plus_bar, "aretes_minus_bar":aretes_minus_bar}            
            
def aretes_orientes(matE_reel, oriente=True):
    """
    retourner tous les arcs du graphe/sous graphe reel
    """
    arcs = list()
    for row, col in fct_aux.range_2d(matE_reel.columns.tolist()):
        if matE_reel.loc[row,col] == 1 and oriente:
            arcs.append((row,col))
        elif matE_reel.loc[row,col] == 1 and not oriente:
            arcs.append((row,col)); arcs.append((col,row));
        else:
            print("row={},col={} do not have a corresponding case in matE_reel".format(row,col))
    return arcs;
    
def intersection_arcs(aretes_minus, aretes_minus_bar):
    """
    aretes appartenant aux deux listes
    """
    return aretes_minus.intersection(aretes_minus_bar)
    
####### ajout corr moyenne
def ajout_corr_moyenne(matE_gr, aretes_minus, aretes_minus_bar, \
                       aretes_plus, m_minus2, m_plus_minus, cpt_modif):
    """
    """
    if len(set(aretes_minus).intersection(set(aretes_minus_bar))) > 0:
        for arete_minus_inter in set(aretes_minus).intersection(set(aretes_minus_bar)):
            aretes_ = set(aretes_minus) - set(aretes_minus).intersection(set([arete_minus_inter]))
            for arete_minus in aretes_:
                arete_minus_inter_, arete_minus = verifier_row_col(arete_minus_inter, arete_minus, matE_gr)
                if arete_minus_inter_ != None and arete_minus != None:
#                    print("10 ({},{})={}".format(arete_minus_inter_, arete_minus, m_minus2))
                    matE_gr.loc[arete_minus_inter_, arete_minus] = m_minus2;
                    matE_gr.loc[arete_minus, arete_minus_inter_] = m_minus2;
                    cpt_modif += 1;
                else:
#                    print("10 aretes {},{} dont belong to matE_gr".format(arete_minus_inter_, arete_minus));
                    pass
                          
            for arete_plus in aretes_plus:
                arete_minus_inter_, arete_plus = verifier_row_col(arete_minus_inter, arete_plus, matE_gr);
                if arete_minus_inter_ != None and arete_plus != None:
#                    print("11 ({},{})={}".format(arete_minus_inter_, arete_plus, m_plus_minus))
                    matE_gr.loc[arete_minus_inter_, arete_plus] = m_plus_minus;
                    matE_gr.loc[arete_plus, arete_minus_inter_] = m_plus_minus;
                    cpt_modif += 1;
                else:
#                    print("11 aretes {},{} dont belong to matE_gr".format(arete_minus_inter_, arete_plus));
                    pass
                          
    return matE_gr, cpt_modif;
####### ajout corr moyenne

####### ajout corr moyenne bar
def modifier_matE_same_degre(matE,aretes,M, cpt_modif):
    for arete1, arete2 in fct_aux.range_2d(aretes):
        arete1_, arete2_ = verifier_row_col(arete1, arete2, matE);
        if arete1_ != None and arete2_ != None and arete1_ != arete2_:
            print("12 ({},{})={}".format(arete1_, arete2_, M))
            matE.loc[arete1_, arete2_] = M;
            matE.loc[arete2_, arete1_] = M;
            cpt_modif += 1;
        elif arete1_ == arete2_:
            print("12 {}, {} identiques".format(arete1_,arete2_))
            pass
        else:
            print("12 aretes {},{} dont belong to matE_gr".format(arete1, arete2));
            pass
    return matE, cpt_modif;
def modifier_matE_inter_degre(matE,aretes_plus,aretes_minus,M, cpt_modif):
    for tuple_ in it.product(aretes_plus, aretes_minus):
        arete1, arete2 = verifier_row_col(tuple_[0], tuple_[1], matE);
        if arete1 != None and arete2 != None and arete1 != arete2:
            print("13 ({},{})={}".format(arete1, arete2, M))
            matE.loc[arete1, arete2] = M;
            matE.loc[arete2, arete1] = M;
            cpt_modif+= 1
        elif arete1 == arete2:
            print("13 {}, {} identiques".format(tuple_[0],tuple_[1]))
            pass
        else:
            print("13 aretes {} et {} dont belong to matE_gr".format(tuple_[0],tuple_[1]) ) 
            pass
    return matE, cpt_modif;
def ajout_corr_moyenne_bar(matE_gr, aretes_minus, aretes_minus_bar, \
                           aretes_plus, aretes_plus_bar, M_minus2, M_plus2,\
                           M_plus_minus, m_plus2, m_minus2, cpt_modif):
    """
    """
    if len(aretes_minus_bar) == len(aretes_minus) and len(aretes_plus_bar) == len(aretes_plus):
        print("12-13 1")
        matE_gr, cpt_modif = modifier_matE_same_degre(matE_gr.copy(), aretes_minus, M_minus2, cpt_modif)
        matE_gr, cpt_modif = modifier_matE_same_degre(matE_gr.copy(), aretes_plus, M_plus2, cpt_modif)
        matE_gr, cpt_modif = modifier_matE_inter_degre(matE_gr.copy(), aretes_plus_bar, aretes_minus_bar, M_plus_minus, cpt_modif)
    
    elif len(aretes_minus_bar) == 0 and len(aretes_plus_bar) != 0:
        print("12-13 2")
        matE_gr, cpt_modif = modifier_matE_same_degre(matE_gr.copy(), aretes_plus_bar, M_plus2, cpt_modif);
        matE_gr, cpt_modif = modifier_matE_inter_degre(matE_gr.copy(), aretes_plus_bar, aretes_minus, M_plus_minus, cpt_modif)
        matE_gr, cpt_modif = modifier_matE_inter_degre(matE_gr.copy(), aretes_plus_bar, aretes_plus, m_plus2, cpt_modif)
        
    elif len(aretes_plus_bar) == 0 and len(aretes_minus_bar) != 0:
        print("12-13 3")
        matE_gr, cpt_modif = modifier_matE_same_degre(matE_gr.copy(), aretes_minus_bar, M_minus2, cpt_modif);
        matE_gr, cpt_modif = modifier_matE_inter_degre(matE_gr.copy(), aretes_minus_bar, aretes_plus, M_plus_minus, cpt_modif)
        matE_gr, cpt_modif = modifier_matE_inter_degre(matE_gr.copy(), aretes_minus_bar, aretes_minus, m_minus2, cpt_modif)
        
    elif len(aretes_plus_bar) != 0 and len(aretes_minus_bar) != 0:
        print("12-13 4")
        matE_gr, cpt_modif = modifier_matE_same_degre(matE_gr.copy(), aretes_plus_bar, M_plus2, cpt_modif);
        matE_gr, cpt_modif = modifier_matE_same_degre(matE_gr.copy(), aretes_minus_bar, M_minus2, cpt_modif);
        matE_gr, cpt_modif = modifier_matE_inter_degre(matE_gr.copy(), aretes_minus_bar, aretes_plus_bar, M_plus_minus, cpt_modif)
        
    return matE_gr, cpt_modif;
####### ajout corr moyenne bar

##### calcul M_plus2, M_minus2, M_plus_minus
def calcul_M_plus_minus2(sommet,info_sommets,M_plus2s,M_minus2s,M_plus_minuss, nbre_sommets, diviseAllItems):
    """
    calcul/mise a jour de M_plus2, M_minus2, M_plus_minus
    info_sommets = {m_plus2":, "aretes_plus":, "m_minus2":,\
                            "aretes_minus":, "m_plus_minus":,"aretes_plus_bar":, "aretes_minus_bar":}
                            
    diviseAllItems: diviser par tous les elements de laliste sans tenir compte des items 0 de cette liste
    """
    M_plus2s[sommet] = info_sommets["m_plus2"];
    M_minus2s[sommet] = info_sommets["m_minus2"];
    M_plus_minuss[sommet] = info_sommets["m_plus_minus"];
    
    nbre_plus2_0 = 0; nbre_minus2_0 = 0; nbre_plus_minus_0 = 0;
    if diviseAllItems:
        nbre_plus2_0 = M_plus2s.values().count(0); nbre_minus2_0 = M_minus2s.values().count(0);
        nbre_plus_minus_0 = M_plus_minuss.values().count(0);
    
    M_plus2 = sum(M_plus2s.values())/(nbre_sommets-nbre_plus2_0) if nbre_sommets != nbre_plus2_0 \
                                                                else 0;
    M_minus2 = sum(M_minus2s.values())/(nbre_sommets-nbre_minus2_0) if nbre_sommets != nbre_minus2_0 \
                                                                    else 0;
    M_plus_minus = sum(M_plus_minuss.values())/(nbre_sommets-nbre_plus_minus_0) \
                                                        if nbre_sommets != nbre_plus_minus_0 \
                                                        else 0;
    
    return M_plus2, M_minus2, M_plus_minus;
##### calcul M_plus2, M_minus2, M_plus_minus


def ajouter_mesures_correlations(args):
    """
    pour chaque grandeur, ajouter les correlations pour des paires d'aretes manquantes
    """
    matE_gr = None;
    if "matE_gp" not in args.keys():
        matE_gr = pd.read_csv(args["chemin_matrices"]+"matrice_adjacence_proba_"+args["grandeur"]+".csv",index_col = "Unnamed: 0");
    else:
        matE_gr = args["matE_gp"];
    
    reseauG = pd.read_csv(args["chemin_matrices"]+"reseauG_reelles"+".csv",index_col = "Unnamed: 0");
    df_gr = pd.read_csv(args["chemin_datasets"]+"dataset_"+args["grandeur"]+".csv");
    if "timestamp" in df_gr.columns.tolist():
        df_gr = df_gr.set_index("timestamp");
    else:
        df_gr = df_gr.set_index("Unnamed: 0");
        
    if args["dbg_ajoutCorrel"]:
        dico={('x->y'):[0,1,0,0,0],('x->t'):[1,0,1,1,1],('z->t'):[0,1,0,1,1],('t->w'):[0,1,1,0,1],('t->v'):[0,1,1,1,0]}
        matE_gr = pd.DataFrame(dico, columns=[('x->y'),('x->t'),('z->t'),('t->w'),('t->v')], index=[('x->y'),('x->t'),('z->t'),('t->w'),('t->v')]).T
        dico={"x":[0,1,1,0,0,0],'y':[0,0,0,0,0,0],'t':[0,0,0,0,1,1],'z':[0,0,1,0,0,0],'v':[0,0,0,0,0,0],'w':[0,0,0,0,0,0]}
        reseauG = pd.DataFrame(dico, columns=['x','y','t','z','v','w'], index=['x','y','t','z','v','w']).T
        dico={'timestamp':[1,2,3,4,5],'x->y':[2]*5,'x->t':[3]*5,'t->v':[4]*5,\
              'z->t':[np.nan]*5,'t->w':[np.nan]*5}
        df_gr = pd.DataFrame(dico).set_index('timestamp');
    
    matE_gr_cp = matE_gr.copy();
    print("matE_gr={} et matE_gr_cp={} ".format(\
              len(set(matE_gr.columns.tolist())), len(set(matE_gr_cp.columns.tolist()))))
    
#    M_plus2 = 0; M_minus2 = 0; M_plus_minus = 0; dico_sommets = dict();
    M_plus2s = dict(); M_minus2s = dict(); M_plus_minuss = dict(); dico_sommets = dict();
    diviseAllItems = False;
    for sommet in reseauG.index:
        info_sommets = gamma_sommet(sommet, reseauG, df_gr, matE_gr)
        dico_sommets[sommet] = info_sommets;
        print("sommet={},aretes--={}, aretes++={}, aretes_bar--={}, aretes_bar++={}".\
                  format(sommet,info_sommets["aretes_minus"],info_sommets["aretes_plus"],\
                         info_sommets["aretes_minus_bar"],info_sommets["aretes_plus_bar"]))
        M_plus2s[sommet] = info_sommets["m_plus2"];
        M_minus2s[sommet] = info_sommets["m_minus2"];
        M_plus_minuss[sommet] = info_sommets["m_plus_minus"];
    cpt_modif=0;
    for sommet in reseauG.columns:
        if sommet == "R486":
            print("aretes--={}, aretes++={}, aretes_bar--={}, aretes_bar++={}".\
                  format(dico_sommets[sommet]["aretes_minus"],dico_sommets[sommet]["aretes_plus"],\
                         dico_sommets[sommet]["aretes_minus_bar"],dico_sommets[sommet]["aretes_plus_bar"]))
        #### calcul M_plus2, M_minus2, M_plus_minus 
        M_plus2, M_minus2, M_plus_minus = calcul_M_plus_minus2(sommet,dico_sommets[sommet],\
                                                             M_plus2s,M_minus2s,M_plus_minuss,\
                                                             len(reseauG.index), diviseAllItems);
        print("sommet={}, M++={}, M--={}, M+-={}, m++={}, m--={}".format(sommet,round(M_plus2,3),round(M_minus2,3),round(M_plus_minus,3),round(dico_sommets[sommet]["m_plus2"],3), round(dico_sommets[sommet]["m_minus2"],3)))
        
        matE_gr, cpt_modif = ajout_corr_moyenne(matE_gr.copy(), dico_sommets[sommet]["aretes_minus"],\
                                   dico_sommets[sommet]["aretes_minus_bar"], \
                                   dico_sommets[sommet]["aretes_plus"], \
                                   dico_sommets[sommet]["m_minus2"],\
                                   dico_sommets[sommet]["m_plus_minus"], cpt_modif)
        matE_gr, cpt_modif = ajout_corr_moyenne(matE_gr.copy(), dico_sommets[sommet]["aretes_plus"],\
                                   dico_sommets[sommet]["aretes_plus_bar"], \
                                   dico_sommets[sommet]["aretes_minus"], \
                                   dico_sommets[sommet]["m_plus2"],\
                                   dico_sommets[sommet]["m_plus_minus"], cpt_modif)
        matE_gr, cpt_modif = ajout_corr_moyenne_bar(matE_gr.copy(), dico_sommets[sommet]["aretes_minus"],\
                                   dico_sommets[sommet]["aretes_minus_bar"], \
                                   dico_sommets[sommet]["aretes_plus"], \
                                   dico_sommets[sommet]["aretes_plus_bar"], \
                                   M_minus2,M_plus2,M_plus_minus, \
                                   dico_sommets[sommet]["m_plus2"],dico_sommets[sommet]["m_minus2"],\
                                   cpt_modif)
        dico_sommets[sommet] = gamma_sommet(sommet, reseauG, df_gr, matE_gr);
    comparateur_correlation(matE_gr, matE_gr_cp, cpt_modif);
#    matE_gr.to_csv(args["chemin_matrices"]+"test_comprehension.csv")
    return matE_gr;

######### debug grave ############33   
def find_correl_arc_entrant_sortant_sommet(sommet,matE_gr):
    arcs_plus2 = []; arcs_minus2 = [];
    for arete in matE_gr.columns:
        if arete.split("->")[0] == sommet:
            arcs_plus2.append(arete)
        elif arete.split("->")[1] == sommet:
            arcs_minus2.append(arete)
    print("sommet={}, arcs_minus2={}, arcs_plus2={}".format(sommet,arcs_minus2,arcs_plus2))
    m_minus2 = 0; m_plus2=0; dico_minus=dict(); dico_plus=dict();
    for arete1, arete2 in fct_aux.range_2d(arcs_minus2):
        m_minus2 += matE_gr.loc[arete1,arete2]
        dico_minus[(arete1,arete2)] = matE_gr.loc[arete1,arete2]
    for arete1, arete2 in fct_aux.range_2d(arcs_plus2):
        m_plus2 += matE_gr.loc[arete1,arete2]
        dico_plus[(arete1,arete2)] = matE_gr.loc[arete1,arete2]
    print("sommet={}, m_minus2={}, m_plus2={}".format(sommet,m_minus2,m_plus2))
    print("plus={}, minus={}".format(dico_plus,dico_minus))
    pass

######### debug grave ############33    
#import distribution_correlation as distrib_correls;

if __name__ == '__main__':    
    sub = True;
    reseau = "champlan";
    rep = "champlan_newDataset"; # or champlan or velizy
    root = "/home/willy/topologyLearning/datas/"+rep+"/";
    root = root+"sous_graphe/" if sub else root;
    chemin_datasets = root+"datasets/";
    chemin_matrices = root+"matrices/";
    chemin_equipements = root+"matrices_equipements/";
    args_chemins={chemin_datasets:chemin_datasets,
                  chemin_matrices:chemin_matrices,chemin_equipements:chemin_equipements};
    
    methode_correl = "correlation_par_grandeur" #"correlation_par_metrique_fusion";
    metrique_distance = "lb_keogh" #"fdtw_problem"#"lb_keogh" #"sax"#"pearson";
    mode_correlation = "correlationParMorceaux" #"correlationGlissante"#"correlationParMorceaux" # ou liste_mode = ["correlationParMorceaux","correlationGlissante"]
    indice_max_correlation = 0; derivate = False;
    interval_derivate = 10; epsilon_sax = 0.8;
    selected_grandeurs = list(); fenetre = 50;
    if reseau == "velizy":
        selected_grandeurs = ['AvgI1', 'AvgI']
    elif reseau == "champlan":
        selected_grandeurs = ["P","I","U","S","E","U1","U2","U3","I1","I2","I3"]#["P","I"]
        selected_grandeurs = ["P"];
    else:
        print("reseau ={} inexistant".format(reseau))
        
    matE_gr_res = None
    dbg = False#True#False#True;
    args = {"chemin_datasets":chemin_datasets,"chemin_matrices":chemin_matrices,\
            "chemin_equipements":chemin_equipements,\
            "methode_correl": methode_correl,"mode_correlation": mode_correlation,\
            "indice_max_correlation":indice_max_correlation,\
            "selected_grandeurs": selected_grandeurs,\
            "fenetre": fenetre,"derivate": derivate, "epsilon_sax": epsilon_sax,\
            "interval_deriv":interval_derivate,\
            "metrique_distance":metrique_distance, \
            "dbg_ajoutCorrel":dbg}
    metrique_distances = ["lb_keogh","pearson","lcs"]
    metrique_distances = ["metrique_wil"]
    for metrique_distance in metrique_distances:
        for grandeur in selected_grandeurs:
            args["grandeur"] = grandeur; args["metrique_distance"] = metrique_distance;
            matE_gr_res = ajouter_mesures_correlations(args);
            
            
    # debug  fonctions 
    matE_gr = pd.read_csv(args["chemin_matrices"]+"matrice_adjacence_proba_"+args["grandeur"]+".csv",index_col = "Unnamed: 0");
    reseauG = pd.read_csv(args["chemin_matrices"]+"reseauG_reelles"+".csv",index_col = "Unnamed: 0");
    bool_entr_sort = False;
    df_gr = pd.read_csv(args["chemin_datasets"]+"dataset_"+args["grandeur"]+".csv");
    if "timestamp" in df_gr.columns.tolist():
        df_gr = df_gr.set_index("timestamp");
    else:
        df_gr = df_gr.set_index("Unnamed: 0");
        
#    dico_plus, dico_bar_plus = gamma(reseauG,False, df_gr, matE_gr);
#    dico_minus, dico_bar_minus = gamma(reseauG,True, df_gr, matE_gr)
#    sommets = reseauG.index.tolist()
#    M_plus2, M_minus2, M_plus_minus = calculer_moyenne_correlations(dico_plus, dico_bar_plus,\
#                                                      dico_minus, dico_bar_minus, sommets, matE_gr)
#    find_correl_arc_entrant_sortant_sommet("TGBT2",matE_gr);
    
#    # plot distribution valeurs correl selon correl vrai positif et correl vrai negatif  ==> A EFFACER
#    import distribution_correlation as distrib_correls;
#    args["matE_proba"] = matE_gr_res
##    args["matE_proba"] = matE_gr
#    args["matE_reel"] = pd.read_csv(chemin_matrices+"matE_reelles.csv",index_col = "Unnamed: 0");
#    args["path_save"] = chemin_matrices;
#    args["dbg_0_1"] = False;
#    distrib_correls.distribution_case_0_1_graphe_reelle(args)