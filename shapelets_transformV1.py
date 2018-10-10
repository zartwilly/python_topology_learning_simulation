#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 17:17:00 2017

@author: willy

shapelets transform v1
"""
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt        # to plot any graph
import fonctions_auxiliaires as fct_aux;
import time;
from multiprocessing import cpu_count, Process, Queue
from psutil import virtual_memory;

#from pylab import rcParams
#rcParams['figure.figsize'] = (6, 6)      # setting default size of plots

#def calculate_entropy(probabilities):
#    return sum([-prob * np.log(prob)/np.log(2) if prob != 0 else 0 for prob in probabilities])
#    
#def manhattan_distance(a, b, min_dist=float('inf')):
#    dist = 0
#    for x, y in zip(a, b):
#        dist += np.abs(float(x)-float(y))
#        if dist >= min_dist: return None
#    return dist;
#
#def calculate_dict_entropy(data):
#    counts = {}
#    for entry in data:
#        if entry[1] in counts: counts[entry[1]] += 1
#        else: counts[entry[1]] = 1
#    return calculate_entropy(np.divide(list(counts.values()), float(sum(list(counts.values())))))
#
#def find_best_split_point(histogram):
#    histogram_values = list(itertools.chain.from_iterable(list(histogram.values())))
#    prior_entropy = calculate_dict_entropy(histogram_values)
#    best_distance, max_ig = 0, 0
#    best_left, best_right = None, None
#    for distance in histogram:
#        data_left = []
#        data_right = []
#        for distance2 in histogram:
#            if distance2 <= distance: data_left.extend(histogram[distance2])
#            else: data_right.extend(histogram[distance2])
#        ig = prior_entropy - (float(len(data_left))/float(len(histogram_values))*calculate_dict_entropy(data_left) + \
#             float(len(data_right))/float(len(histogram_values)) * calculate_dict_entropy(data_right))
#        if ig > max_ig: best_distance, max_ig, best_left, best_right = distance, ig, data_left, data_right
#    return max_ig, best_distance, best_left, best_right
#    
#def subsequence_dist(time_serie, sub_serie):
#    if len(sub_serie) < len(time_serie):
#        min_dist, min_idx = float("inf"), 0
#        for i in range(len(time_serie)-len(sub_serie)+1):
#            dist = manhattan_distance(sub_serie, time_serie[i:i+len(sub_serie)], min_dist)
#            if dist is not None and dist < min_dist: min_dist, min_idx = dist, i
#        return min_dist, min_idx
#    else:
#        return None, None  
#        
#def check_candidate(data, shapelet):
#    histogram = {} 
#    for entry in data:
#        # TODO: entropy pre-pruning in each iteration
#        time_serie, label = entry[0], entry[1]
#        d, idx = subsequence_dist(time_serie, shapelet)
#        if d is not None and not np.isinf(d):
#            histogram[d] = [(time_serie, label)] \
#            if d not in histogram \
#            else histogram[d].append((time_serie, label))
#    return find_best_split_point(histogram);
#    
#def generate_candidates(data, max_len=5, min_len=2):
#    candidates, l = [], max_len
#    while l >= min_len:
#        for i in range(len(data)):
#            time_serie, label = data[i][0], data[i][1]
#            for k in range(len(time_serie)-l+1): 
#                candidates.append((time_serie[k:k+l], label))
#        l -= 1
#    return candidates
#    
#def find_shapelets_bf(data, max_len=100, min_len=1, plot=True, verbose=True):
#    candidates = generate_candidates(data, max_len, min_len)
#    print("len data={}".format(len(data)))
#    bsf_gain, bsf_shapelet = 0, None
#    if verbose: candidates_length = len(candidates)
#    for idx, candidate in enumerate(candidates):
#        gain, dist, data_left, data_right = check_candidate(data, candidate[0])
#        if verbose: print(idx, '/', candidates_length, ":", gain, dist) if idx%candidates_length == 0 \
#                                                                        else None;
#        if gain > bsf_gain:
#            bsf_gain, bsf_shapelet = gain, candidate[0]
#            if verbose:
#                print('Found new best shapelet with gain & dist:', \
#                      bsf_gain, dist, [x[1] for x in data_left], \
#                      [x[1] for x in data_right])
#            if plot:
#                plt.plot(bsf_shapelet)
#                plt.show()
#            plt.show()
#    return bsf_shapelet
#    
#def extract_shapelets(data, min_len=10, max_len=15, verbose=1):
#    _classes = np.unique([x[1] for x in data])
#    shapelet_dict = {}
#    for _class in _classes:
#        print('Extracting shapelets for', _class)
#        transformed_data = []
#        for entry in data:
#            time_serie, label = entry[0], entry[1]
#            print("label={}, time_serie={}".format(label,time_serie))
#            if label == _class: transformed_data.append((time_serie, 1))
#            else: transformed_data.append((time_serie, 0))
#        shapelet_dict[_class] = find_shapelets_bf(transformed_data, max_len=max_len, \
#                                    min_len=min_len, plot=0, verbose=1)
#    return shapelet_dict

########################################333 ddebug ====> DEBUT
def normalize(df):
    result = df.copy()
#    print("shape result = {}".format(result.shape) )
    for feature_name in df.columns.tolist():
        mean = df[feature_name].mean(skipna=True);
        std = df[feature_name].std(skipna=True);
        if df[feature_name].isnull().all() != True and std !=0 and np.isnan(std) and np.isnan(mean):
            result[feature_name] = (df[feature_name] - mean) / (std)
#    print("result %s ", result.isnull().all().any())
    return result;
def subsequence_dist(time_serie, sub_serie):
    if len(sub_serie) < len(time_serie):
        min_dist, min_idx = float("inf"), 0
        for i in range(len(time_serie)-len(sub_serie)+1):
            dist = manhattan_distance(sub_serie, time_serie[i:i+len(sub_serie)], min_dist)
            if dist is not None and dist < min_dist: min_dist, min_idx = dist, i
        return min_dist, min_idx
    else:
#        return None, None;
        return np.nan, np.nan;
def generate_candidates(data, max_len=5, min_len=2):
    candidates, l = [], max_len
    while l >= min_len:
        for i in range(len(data)):
            time_serie, label = data[i][0], data[i][1]
            for k in range(len(time_serie)-l+1): candidates.append((time_serie[k:k+l], label))
        l -= 1
    return candidates
def check_candidate(data, shapelet):
    histogram = dict() 
    for entry in data:
        # TODO: entropy pre-pruning in each iteration
        time_serie, label = entry[0], entry[1]
        d, idx = subsequence_dist(time_serie, shapelet)
        if d is not None and not np.isinf(d):
            histogram[d] = [(time_serie, label)] if d not in histogram.keys() \
                                                 else histogram[d].append((time_serie, label))
    return find_best_split_point(histogram)
def calculate_dict_entropy(data):
    counts = {}
    for entry in data:
        if entry[1] in counts: counts[entry[1]] += 1
        else: counts[entry[1]] = 1
    return calculate_entropy(np.divide(list(counts.values()), float(sum(list(counts.values())))))
def find_best_split_point(histogram):
    histogram_values = list(itertools.chain.from_iterable(list(histogram.values())))
    prior_entropy = calculate_dict_entropy(histogram_values)
    best_distance, max_ig = 0, 0
#    best_left, best_right = None, None
    best_left, best_right = np.nan, np.nan
    for distance in histogram:
        data_left = []
        data_right = []
        for distance2 in histogram:
            if distance2 <= distance: data_left.extend(histogram[distance2])
            else: data_right.extend(histogram[distance2])
        ig = prior_entropy - (float(len(data_left))/float(len(histogram_values))*calculate_dict_entropy(data_left) + \
             float(len(data_right))/float(len(histogram_values)) * calculate_dict_entropy(data_right))
        if ig > max_ig: best_distance, max_ig, best_left, best_right = distance, ig, data_left, data_right
    return max_ig, best_distance, best_left, best_right

def manhattan_distance(a, b, min_dist=float('inf')):
    dist = 0
    for x, y in zip(a, b):
        dist += np.abs(float(x)-float(y))
        if dist >= min_dist: return np.nan #None
    return dist
def calculate_entropy(probabilities):
    return sum([-prob * np.log(prob)/np.log(2) if prob != 0 else 0 for prob in probabilities])
def find_shapelets_bf(data, max_len=100, min_len=1, plot=True, verbose=True):
    candidates = generate_candidates(data, max_len, min_len)
#    bsf_gain, bsf_shapelet = 0, None
    bsf_gain, bsf_shapelet = 0, np.nan
    if verbose: candidates_length = len(candidates)
    for idx, candidate in enumerate(candidates):
        gain, dist, data_left, data_right = check_candidate(data, candidate[0])
        if verbose: print(idx, '/', candidates_length, ":", gain, dist)
        if gain > bsf_gain:
            bsf_gain, bsf_shapelet = gain, candidate[0]
            if verbose:
                print('Found new best shapelet with gain & dist:', \
                      bsf_gain, dist, [x[1] for x in data_left], \
                      [x[1] for x in data_right])
            if plot:
                plt.plot(bsf_shapelet)
                plt.show()
            plt.show()
    return bsf_shapelet
def extract_shapelets(data, min_len=10, max_len=15, verbose=1):
    _classes = np.unique([x[1] for x in data])
    shapelet_dict = {}
    for _class in _classes:
        print('Extracting shapelets for', _class)
        transformed_data = []
        for entry in data:
            time_serie, label = entry[0], entry[1]
            if label == _class: transformed_data.append((time_serie, 1))
            else: transformed_data.append((time_serie, 0))
        shapelet_dict[_class] = find_shapelets_bf(transformed_data, max_len=max_len, min_len=min_len, plot=0, verbose=1)
    return shapelet_dict
########################################333 ddebug  ====> FIN

def data_thread(entries, _class, q):
    for entry in entries:
        time_serie, label = entry[0], entry[1]
        if label == _class: q.put((time_serie, 1))
        else: q.put((time_serie, 0))  
    
def extract_shapelets_v1(data, min_len=10, max_len=15, verbose=1):
#    print("type_d={}, data={}".format(type(data),len(data)))
    _classes = np.unique([x[1] for x in data])
    shapelet_dict = {};
    for _class in _classes:
        print('Extracting shapelets for', _class)
        transformed_data = []
        indices  = [int(i) for i in np.linspace(0,len(data),cpu_count())]
        entries = []
        processes = []
        for i in range(len(indices)-1):
             d = data[indices[i]:indices[i+1]]
#             print("len d={}".format(len(d)))
             q = Queue()
             p = Process(target=data_thread, args=(d, _class, q))
             p.start()
             entries.append(q); processes.append(p);
        for process in processes :
             process.join()  
        
        n_entries = [entry.get(0) for entry in entries]
        for entry in n_entries:
            transformed_data.append(entry);
#            print("entry={}".format(entry))
        
        shapelet_dict[_class] = find_shapelets_bf(transformed_data, max_len=max_len, \
                                    min_len=min_len, plot=0, verbose=1)
    return shapelet_dict
    
###### implementation de shapelet transform ####### debut
def label_arcs(arcs):
    dico = dict();
    for arc in arcs:
        dico[arc]= arc; #cpt+=1
    return dico;

###### implementation de shapelet transform ####### fin

if __name__ == '__main__':
    
    start= time.time();
    bool_creation_datasets = False;
    sub = True;
    reseau = "champlan";
    rep = "champlan_newDataset"; # or champlan or velizy
    root = "/home/willy/topologyLearning/datas/"+rep+"/";
    root = root+"sous_graphe/" if sub else root;
    chemin_datasets = root+"datasets/";
    chemin_matrices = root+"matrices/";
    chemin_equipements = root+"matrices_equipements/";
    
    selected_grandeurs = fct_aux.liste_grandeurs(chemin_datasets);
    selected_grandeurs = ["P"];
    
    df_data_bis = pd.DataFrame()    
    for grandeur in selected_grandeurs:
        df_gr = pd.read_csv(chemin_datasets+"dataset_"+grandeur+".csv");
        df_gr = df_gr.set_index("timestamp") if "timestamp" in df_gr.columns \
                                             else df_gr.set_index("Unnamed: 0");
        cols_to_delete = list(df_gr.columns[df_gr.isnull().all()])
        df_gr.drop(cols_to_delete, axis = 1, inplace = True);
        df_gr_norm = (df_gr - df_gr.mean(skipna=True))/df_gr.std(skipna=True);
        df_gr_norm.fillna(method='pad',inplace = True);
        df_gr_norm = df_gr_norm.rolling(window=20).mean();
        df_gr_norm = df_gr_norm.loc[1359390657:1359477057]
        df_gr_norm.drop(df_gr_norm.columns[df_gr_norm.isnull().all()], axis=1,inplace=True)
        
        label_aretes = label_arcs(df_gr_norm.columns.tolist());
        min_len = 99; max_len = 100;
        min_len = 140; max_len = 141;
        min_len = 48; max_len = 50;
        data = []
        for arete, label in label_aretes.items():
            data.append((df_gr_norm[arete].values, label))
#        shapelet_dict = extract_shapelets_v1(data, min_len, max_len)
        shapelet_dict = extract_shapelets(data, min_len, max_len)

        df_data_bis = pd.DataFrame(shapelet_dict)
    print (time.time() - start)
    