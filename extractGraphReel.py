#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 18:04:58 2017

@author: willy

NB  executer sur 
    * VM  ssh wil@verif-data.com
    * cd /home/dev/soft/titan-0.5.3-hadoop2
    * bin/gremlin.sh conf/{champlan, velizy} 
    * . /home/wil/extractGraph.groovy 
    * extractGraphBis(g)
    * scp wil@verif-data.com:/home/wil/graphDico_undirected.json /home/willy/Documents/erdfIdLogisticScripts/fromGremlinTopython/
    
"""
import pandas as pd;
import numpy as np;
import time;
import fonctions_auxiliaires as fct_aux;
import matplotlib.pyplot as plt;

def range_2d_str(cols):
    """Produce a stream of two-D coordinates."""
    treated = list()
    for col in cols:
        for row in cols:
            if row != col:
                yield row, col;
#            if row != col and (row,col) not in treated and (col,row) not in treated:
#                treated.append( (row,col) )
#                yield row, col;

##### extraction du graphe d'adjacence des aretes ===> debut
def matrice_adjacence_aretes_from_groovyFile(args):
    """
    creation de la matrice d'adjacence des arcs a partir du reseauG_reelles.json
    """
    dico_json = dict();
    if "dico_json" not in args.keys():
        dico_json = fct_aux.export_from_json_to_dico(args["path_reseauG"]);
    else:
        dico_json = args["dico_json"]
    dico_adj = dict(); aretes_G_tmp = set()
    for ext_source, ext_dests in dico_json.items():
        for ext_dest in ext_dests:
            if (ext_dest,ext_source) not in aretes_G_tmp:
                aretes_G_tmp.add((ext_source,ext_dest))
    aretes_G = set()
    for src_dest1, src_dest2 in fct_aux.range_2d(aretes_G_tmp):
        aretes_G.add(src_dest1[0]+"->"+src_dest1[1]);
        aretes_G.add(src_dest2[0]+"->"+src_dest2[1])
        if src_dest1[0] == src_dest2[0] or src_dest1[0] == src_dest2[1] or \
           src_dest1[1] == src_dest2[0] or src_dest1[1] == src_dest2[1] :
            dico_adj[(src_dest1[0]+"->"+src_dest1[1], src_dest2[0]+"->"+src_dest2[1])] = 1;
    
    l_aretes_G = list(aretes_G)
    matE_reel = pd.DataFrame(columns = l_aretes_G, index = l_aretes_G);
    for arete, val in dico_adj.items():
        matE_reel.loc[arete[0],arete[1]] = 1;
        matE_reel.loc[arete[1],arete[0]] = 1;
        
    matE_reel.fillna(0, inplace = True);
    matE_reel.to_csv(args["chemin_matrices"]+"matE_reelles.csv")
    return matE_reel;
##### extraction du graphe d'adjacence des aretes ===> fin

##### extraction du graphe d'adjacence des sommets ===> debut
def matrice_adjacence_sommets_from_groovyFile(args,directed=True):
    """
    creation de la matrice d'adjacence des sommets du graphe/sous graphe oriente/non-Oriente
    a partir du reseauG_reelles.json
    """
    # test si file on specific path exists 
    if not fct_aux.is_file_exist(args["path_reseauG"]):
        print("ERROR: file = {} dont exist".format(args["path_reseauG"]))
        return None;
        
    dico_json = dict();
    dico_json = fct_aux.export_from_json_to_dico(args["path_reseauG"]);
    nodes = set()
    nodes = set([val for k,vals in dico_json.items() for val in vals])
    nodes = nodes.union(set(dico_json.keys()))
    print("matrice_adjacence_from_groovyFile nodes = ", len(nodes))
        
    reseauG = pd.DataFrame(index = nodes, columns = nodes);
    for node in nodes:
        for key, vals in dico_json.items():
            for val in vals:
                if val == node:
                    if directed == False:
                        reseauG.loc[key,node] = 1;
                        reseauG.loc[node,key] = 1;
                    else:
                        reseauG.loc[key,node] = 1;
    reseauG.fillna(0, inplace = True);
    reseauG.to_csv(args["chemin_matrices"]+"reseauG_reelles.csv");
    args["dico_json"] = dico_json;
    print("creation reseauG_reelles.csv ==> TERMINE, nodes={},index={}".format(len(reseauG.columns),len(reseauG.index)))
    matE_reel = matrice_adjacence_aretes_from_groovyFile(args);
    print("creation matE_reelles.csv ==> TERMINE, aretes ={}".format(len(matE_reel.columns)))
    return matE_reel;
##### extraction du graphe d'adjacence des sommets ===> FIN

#### reduire datasets ==> debut
def modifier_index(df_gr):
    """
    supprimer les 3 chiffres de la fin de chaque timestamp
    """
    dico = dict()
    for ind in df_gr.index:
        dico[ind] = int(str(ind)[:10])
    df_gr.rename(index = dico, inplace = True)
    return df_gr;
def reduire_dataset(chemin_datasets,date_debut, date_fin, dbg):
    """
    ne conserver que les dates comprises entre date_debut et date_fin
    debut = 
    """
    if dbg:
        date_debut = 1358104756; date_fin = 1359658746;
    grandeurs = fct_aux.liste_grandeurs(chemin_datasets)
    for grandeur in grandeurs:
        df_gr = pd.read_csv(chemin_datasets+"dataset_"+grandeur+".csv");
#        if "timestamp" in df_gr.columns.tolist():
#            df_gr = df_gr.set_index("timestamp");
#        else:
#            df_gr = df_gr.set_index("Unnamed: 0")
        df_gr = df_gr.set_index("timestamp") if "timestamp" in df_gr.columns \
                                             else df_gr.set_index("Unnamed: 0")
        df_gr = modifier_index(df_gr)
        df_gr = df_gr.loc[date_debut:date_fin]
        df_gr.to_csv(args["datasets"]+"dataset_"+grandeur+".csv", index=True);
    pass
#### reduire datasets ==> FIN


####### distribution selon le seuil des fausses positives et fausses negatives correlations
def distribution_faux_negatives_positives_correlations(sub_graph):
    """
    pour chaque correl_seuil, on definit une matE contenant des 0 et 1 puis
    on classe les valeurs 0 et 1 selon les correlations {vrai, faux} negatives et positives et enfin
    plot les distogrammes de chaque type de correlation "0->0,0->1,1->0,1->1"
    """
    chemin_matE = "/home/willy/topologyLearning/matrices/Champlan/";
    chemin_matE_correl = "/home/willy/topologyLearning/matrices/Champlan/matrices_equipements/";
    path_save = "/home/willy/Documents/courbePython/"
    
    correl_seuils = [0.5, 0.6, 0.7, 0.8];
#    correl_seuils = [0.5]
    matE_correl = pd.read_csv(chemin_matE_correl+"matE.csv", index_col = "Unnamed: 0")
    matE_correl = matE_correl.apply(lambda x: abs(x))
    matE_titan = pd.read_csv(chemin_matE+sub_graph+"matE_reelles.csv", index_col = "Unnamed: 0")
    for correl_seuil in correl_seuils:
        dico_0_0 = dict(); dico_0_1 = dict() # faux positives
        dico_1_1 = dict(); dico_1_0 = dict() # faux negatives
        dico_2 = dict() # nodes n'appartenant a aucune linegraph
        matE = matE_correl.copy();
        matE[matE < correl_seuil] = 0; matE[matE >= correl_seuil] = 1;
        for node1, node2 in range_2d_str(matE.columns.tolist()):
            if matE_titan.loc[node1, node2] == matE.loc[node1,node2] and \
               matE.loc[node1,node2] == 0 :
                dico_0_0[(node1,node2)] = matE_correl.loc[node1,node2];
            elif matE_titan.loc[node1, node2] == matE.loc[node1,node2] and \
               matE.loc[node1,node2] == 1:
                dico_1_1[(node1,node2)] = matE_correl.loc[node1,node2];
            elif matE_titan.loc[node1, node2] == 0 and matE.loc[node1,node2] == 1:
                dico_0_1[(node1,node2)] = matE_correl.loc[node1,node2];
            elif matE_titan.loc[node1, node2] == 1 and matE.loc[node1,node2] == 0:
                dico_1_0[(node1,node2)] = matE_correl.loc[node1,node2];
            else :
                dico_2[(node1,node2)] = -2
        
        # distribution selon le seuil
        fig = plt.figure(1);
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
        fig.set_figheight(10); fig.set_figwidth(10);
        
        bins = np.arange(0,1,0.1)
#        sns.distplot(list(dico_0_0.values()),ax = ax1, bins = bins, kde = False);
        ax1.hist(list(dico_0_0.values()), bins )
        ax1.plot([correl_seuil, correl_seuil], (0, 30), 'r--' )
        ax1.set_title('distributions vrai negatives')
        
#        sns.distplot(list(dico_0_1.values()),ax = ax2, bins = bins, kde = False);
        ax2.hist(list(dico_0_1.values()), bins)
        ax2.plot([correl_seuil, correl_seuil], (0, 30), 'r--' )
        ax2.set_title('distributions faux positives')
        
#        sns.distplot(list(dico_1_0.values()),ax = ax3, bins = bins, kde = False);
        ax3.hist(list(dico_1_0.values()),bins )
        ax3.plot([correl_seuil, correl_seuil], (0, 30), 'r--' )
        ax3.set_title('distributions faux negatives')
        
#        sns.distplot(list(dico_1_1.values()),ax = ax4, bins = bins, kde = False);
        ax4.hist(list(dico_1_1.values()), bins)
        ax4.plot([correl_seuil, correl_seuil], (0, 30), 'r--' )
        ax4.set_title('distributions vrai positives')
        
        plt.savefig(path_save+"distributions_{vrai,faux}_{negatives_positives}_correlations_seuil_"\
                    +str(correl_seuil)+".jpeg", dpi= 190)
        plt.clf()
    pass
####### distribution selon le seuil des fausses positives et fausses negatives correlations

if __name__ == '__main__':
#    d = matrice_adjacence_subgraph_from_groovyFile()  # commenter le 04/11/2017
#    d = matrice_adjacence_from_groovyFile()
    sub_graph = "" # "" or "sub_"
#    distribution_faux_negatives_positives_correlations(sub_graph)
    
    start= time.time(); sub = True;
    reseau = "champlan"; fichier = "GraphDico_directed_"+reseau+".json";
    rep = "champlan_newDataset"; # or champlan
    root = "/home/willy/topologyLearning/datas/"+rep+"/";
    root = root+"sous_graphe/" if sub else root;
    chemin_datasets = root+"datasets/";
    chemin_matrices = root+"matrices/";
    chemin_equips = root+"matrices_equipements/";
    path_reseauG= chemin_equips+fichier;
    args={"path_reseauG":path_reseauG, "chemin_matrices":chemin_matrices};
#    matrice_adjacence_aretes_from_groovyFile(args)
    directed = True
    matrice_adjacence_sommets_from_groovyFile(args,directed)
    
    print("fin => ", time.time() - start)