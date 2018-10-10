# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 15:20:24 2016

@author: Max
"""

import pandas as pd
import numpy as np
import sax as sax_
import fonctions_auxiliaires as fct_aux
import matrice_compatible as mat_compa
import time;
import logging;


def Sax_transform(time_series=None,epsilon = 0.01,sub_sequence_size=None,nb_letters = 4,words_size=8):
    '''
    epsilon etant la hauteur du mot
    '''    
        
    if time_series == None :
        print('ERROR : No data in parameter')
        
    if sub_sequence_size == None :
        sub_sequence_size = 4*words_size
    
    #epsilon = 0.01
    print ('epsilon, = ', epsilon)
    sax = sax_.SAX(words_size,nb_letters, epsilon)

    vocabulaire = []
    sentences = []
    for ts in time_series :
        sentence = []
        sub_ts = []
        for i in range(1,len(ts)+1) :
            if i%sub_sequence_size == 0 :
                word = sax.to_letter_rep(sub_ts)[0]
                sentence.append(word)
                vocabulaire.append(word)
                sub_ts = []
            else :
                sub_ts.append(ts[i-1])
        sentences.append(sentence) 
    
    vocabulaire = list(set(vocabulaire))
    
    return sentences, vocabulaire
    
def Sax_to_Vec(sentences,vocabulaire):
    vectors = []
    for sentence in sentences :
        ts = []
        for w in vocabulaire :
            ts.append(float(sentence.count(w)))
        vectors.append(ts)
    return vectors
    
def distance(x,y):
    if type(x) == list :
        x = np.asarray(x)
    if type(y) == list :
        y = np.asarray(y)
    return np.sqrt(np.sum(np.square(x-y)))
    
def create_dict_similitude(names_list,vectors):
    dict_similitude = dict()
    for i in range(len(names_list)):
        reste_vectors = vectors[:]
        reste_names = names_list[:]
        vector = reste_vectors.pop(i)
        name = reste_names.pop(i)
        
        dict_vector = dict()
        for j in range(len(reste_names)) :
            dict_vector[reste_names[j]] = distance(vector,reste_vectors[j])
        
        dict_similitude[name] = dict_vector
        
    return dict_similitude
    
#### ajouter par wilfried
def quicksort(arr):
    if len(arr) <= 1:
        return arr
    pivot = arr[len(arr) // 2]
    left = [x for x in arr if x < pivot]
    middle = [x for x in arr if x == pivot]
    right = [x for x in arr if x > pivot]
    return quicksort(left) + middle + quicksort(right)
    
def creation_matriceAdjacence_symetrique(dico, grandeur, chemin):
    column = dico.keys()
    ind = column
    df = pd.DataFrame(index = ind, columns = column)
    for cle, dico_valeur in dico.items():
        for cle_dico_valeur, valeur in dico_valeur.items():
            df[cle][cle_dico_valeur] = valeur
            df[cle_dico_valeur][cle] = valeur
    
    df = df.apply(lambda x: (100-x)/100)    
    df[df == 100] = 1
    df[df < 0.95] = 0
    
    df.fillna(0, inplace = True)
    df.to_csv(chemin+'matrice_adjacence_proba_'+grandeur+'.csv')

    df[df >= 0.95] = 1    
    df = mat_compa.creer_matrice_compatible( df.astype(int) )
    df.to_csv(chemin+'matrice_adjacence_'+grandeur+'.csv')
    
    
def create_dico_similarite(data, epsilon):
    data_sum = data.copy().sum()
    data_sum.sort_values(ascending=False,inplace=True)
    cles = data_sum.index.tolist()
    time_series = []
    for c in cles:
        time_series.append(data[c].tolist())
        
    #Testing code :
    sentences, vocabulaire = Sax_transform(time_series,epsilon,sub_sequence_size=None,nb_letters = 4,words_size=8)
    new_timeseries = Sax_to_Vec(sentences, vocabulaire)
    dict_similitude = create_dict_similitude(cles,new_timeseries)
    
    return dict_similitude

def supprimer_les_grands_distances(dico, n):
    '''
    n etant le nombre de cle ayant les plus petites distances. 
    #en d'autres mots, n est le nombre de voisins max d'un noeud ==> faux
    methode:
        trouver les n plus petites valeurs et les mettre dans une liste dans un ordre croissant( liste_valeurs_plus_petites )
        chercher l'ecart entre les 2 valeurs consecutifs de 'liste_valeurs_plus_petites' et s'il est superieure a 1 SEUIL alors effacer tous ceux apres lui
        si la valeur d'une cle n'est pas dans 'liste_valeurs_plus_petites' alors mettre la valeur a 0
    '''
    for cle, dico_valeur in dico.items():
        liste_valeurs_dico = list( set(dico_valeur.values()) )
        liste_valeurs_plus_petites = quicksort( liste_valeurs_dico )[:n]
        my_liste_valeurs_plus_petites = np.asarray(liste_valeurs_plus_petites)
        my_liste_valeurs_plus_petites = my_liste_valeurs_plus_petites[ my_liste_valeurs_plus_petites < 10 ]
        for k, v in dico_valeur.items():
            if dico_valeur[k] not in my_liste_valeurs_plus_petites:
                dico_valeur[k] = 1000
    return dico
    
def create_matriceAdjacence_toutes_grandeurs(epsilon, chemin_data = "data/datasets/", chemin_mat = "data/matrices/"):
    liste_grandeur = fct_aux.liste_grandeurs(chemin_data)
    for grandeur in liste_grandeur:
        t1 = time.time()
        location_dataset = chemin_data+'dataset_'+grandeur+'.csv'
        data = pd.read_csv(location_dataset)
        dict_similitude = create_dico_similarite(data, epsilon)
        supp_dict_similitude = supprimer_les_grands_distances(dict_similitude, 10)
        creation_matriceAdjacence_symetrique(supp_dict_similitude, grandeur, chemin_mat)
        t2 = time.time() -t1
        print( grandeur, ' => matrice adjacence  cree en ', t2)

##### new version ==============>
def creation_matriceAdjacence_symetrique_new(dico, dico_aretes, grandeur, chemin):
    columns = dico.keys()
    
    ind = column = list(set([col.split("_"+grandeur)[0].upper() for col in columns]))
    df = pd.DataFrame(index = ind, columns = column)
    for cle, dico_valeur in dico.items():
        for cle_dico_valeur, valeur in dico_valeur.items():
            arete0 = cle.split("_")[0].upper();
            arete1 = cle_dico_valeur.split("_")[0].upper();
            df.loc[arete0][arete1] = valeur
            df.loc[arete1][arete0] = valeur
            dico_aretes = fct_aux.ajouter_correlation_dico(dico_aretes, \
                                                           (arete0,arete1),\
                                                           valeur)
            
    df = df.apply(lambda x: (100-x)/100)    
    df[df == 100] = 1
    
    df.fillna(0, inplace = True)
    df.to_csv(chemin+'matrice_adjacence_proba_'+grandeur+'.csv')

    df[df >= 0.95] = 1;
    df[df < 0.95] = 0;
    df = mat_compa.creer_matrice_compatible( df.astype(int) )
    df.to_csv(chemin+'matrice_adjacence_'+grandeur+'.csv')
    return dico_aretes;
    
def create_matriceAdjacence_toutes_grandeurs_new(epsilon, chemin_data = "data/datasets/",\
                                                 chemin_mat = "data/matrices/"):
    """
    """
    dico_aretes = dict(); # exple dico_aretes = {(a,b):[...],(c,d):[...],...};
    liste_grandeur = fct_aux.liste_grandeurs(chemin_data)
    for grandeur in liste_grandeur:
        t1 = time.time()
        location_dataset = chemin_data+'dataset_'+grandeur+'.csv'
        data = pd.read_csv(location_dataset)
        dict_similitude = create_dico_similarite(data, epsilon)
        supp_dict_similitude = supprimer_les_grands_distances(dict_similitude, 10)
        dico_aretes = creation_matriceAdjacence_symetrique_new(supp_dict_similitude, \
                                                               dico_aretes, grandeur, \
                                                               chemin_mat)
        t2 = time.time() -t1
        print( grandeur, ' => matrice adjacence  cree en ', t2)
    return dico_aretes;
##### new version ==============>
        
### creation matE pour une grandeur ===> debut
def create_matE_grandeur(grandeur, epsilon = 0.8, chemin_data = "data/datasets/",\
                         chemin_mat = "data/matrices/"):
    """
    return dico des arcs avec leur valeur de correlation
    """
    t1 = time.time();
    data = pd.read_csv(chemin_data+'dataset_'+grandeur+'.csv')
    if "timestamp" in data.columns.tolist():
        data = data.set_index("timestamp");
    else:
        data = data.set_index("Unnamed: 0")
        
    logging.debug("correlation_sax = %s columns = %s",grandeur, data.columns.tolist())
        
    dict_similitude = create_dico_similarite(data, epsilon)
    supp_dict_similitude = supprimer_les_grands_distances(dict_similitude, 10);
    dico_gp = dict(); cols = data.columns.tolist(); 
    df = pd.DataFrame(index = cols, columns = cols);
    for cle, dico_valeur in supp_dict_similitude.items():
        for cle_dico_valeur, valeur in dico_valeur.items():
            arete0 = cle.split("_")[0].upper();
            arete1 = cle_dico_valeur.split("_")[0].upper();
            df.loc[arete0, arete1] = valeur;
            df.loc[arete1, arete0] = valeur;
            
    df = df.apply(lambda x: (100-x)/100)    
    df[df == 100] = 1
    for row, col in fct_aux.range_2d(cols):
        dico_gp[(row, col)] = df.loc[row, col]
    
    df.fillna(0, inplace = True)
    df.to_csv(chemin_mat+'matrice_adjacence_proba_'+grandeur+'.csv')
    logging.debug("correlation_sax: grandeur = %s termine en = %s",grandeur, time.time()-t1)
    return dico_gp;
    
### creation matE pour une grandeur ===> fin

        
### ajouter par wilfried    
    
if __name__ == "__main__":
    
    start= time.time()
    chemin = "data/"
    create_matriceAdjacence_toutes_grandeurs()
    
    print (time.time() - start)  
    
    
