#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 08:25:14 2017

@author: willy
creation dataset par grandeur
1- transforme le json en un dico
2- creer un dataFrame tel que columns = nom_arc et ligne = timestamp
3- determiner les grandeurs du dataset
4- construire le dataset par grandeur
"""
import json, ujson;
import time;
import re
import pandas as pd;
import numpy as np;


def export_from_json_to_dico(path_file):
    """
    exporter le json en dico;
    """
    return json.load(open(path_file,"r"));

def trouver_grandeurs(arcs):
    """
    """
    grandeurs = set([ arc.split("_")[1] for arc in arcs if arc != "timestamp"])
    return list(grandeurs);
    
def  old_new_cols(data, grandeur):
    """ renommer les noms des colonnes du dataframe data
    """
    dico_old_new_col = dict()
    for col in data.columns.tolist():
        dico_old_new_col[col] = col.split("_")[0];
    data.rename(columns = dico_old_new_col, inplace = True)
    return data;
    
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
def reduire_dataset(df_gr, date_debut, date_fin, dbg):
    """
    ne conserver que les dates comprises entre date_debut et date_fin
    debut = 
    """
    if dbg:
        date_debut = 1358104756; date_fin = 1359658746;
    else:
        return df_gr;
    
    df_gr = modifier_index(df_gr)
    return df_gr.loc[date_debut:date_fin]
    
#### reduire datasets ==> FIN
def transform_dico_to_dataFrame(args):
    """
    transformer le dico en un dataFrame dont 
        * les columns sont les noms des arcs
        * les lignes sont les valeurs par timestamp
    """
    dico_json = export_from_json_to_dico(args["path_file"]);
    df = pd.DataFrame(dico_json, columns = dico_json.keys())
    df.index.names = ["timestamp"];
    grandeurs = trouver_grandeurs(df.columns.tolist())
    for grandeur in grandeurs:
        print("creation dataset grandeur={} DEBUT".format(grandeur))
        pattern = "_"+grandeur;
        selected_cols = [ arc for arc in df.columns.tolist() if re.search(pattern, arc)];
        df_gr = df[selected_cols].copy();
        df_gr = old_new_cols(df_gr, grandeur);
#       TODO supprimer toutes les lignes n'ayant pas de valeurs.
        df_gr = reduire_dataset(df_gr, 0, 0, args["dbg"]);
        df_gr.to_csv(args["datasets"]+"dataset_"+grandeur+".csv", index=True);
        print("creation dataset grandeur={} TERMINE".format(grandeur))
        
##### test sur compter le nombre de timestamp qui ont plus de 65 items = np.nan  
def compteur_65(df):
    """
    compter tous les timestamps qui ont plus de 65 items = np.nan
    """
#    TODO
    pass
##### test sur compter le nombre de timestamp qui ont plus de 65 items = np.nan     
if __name__ == '__main__':
    
    start= time.time(); sub = True;
    reseau = "champlan"; fichier = "datasets_"+reseau+".json";
    rep = "champlan_newDataset"; # or champlan
    root = "/home/willy/topologyLearning/datas/"+rep+"/";
    root = root+"sous_graphe/" if sub else root;
    datasets = root+"datasets/"
    path_file = datasets+fichier;
    args={"path_file":path_file, "datasets":datasets, "dbg":True};
    transform_dico_to_dataFrame(args);
    
    print("fin => ", time.time() - start)