#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 10:22:11 2016

@author: willy
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 14:24:12 2016

@author: willy
"""


# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 21:20:37 2016

@author: willy
"""

import numpy as np
import pandas as pd
import time


#########################################3

def gamma(matE, xi):
    '''
    ensemble de noeuds adjacents a xi
    '''
    ens = {}
    ens =  set([xj for xj in matE.columns.tolist() if matE.loc[xi][xj] == 1 ]) 
    return ens
    
def maximize_cand_inter_gamma(matE, cand):
    maxi = 0
    for u in cand:
        gama = gamma(matE, u)
        if len(gama) >= maxi:
            maxi = len(gama)
            arc = u
#            print ('u =',u,' gama=',gama,' len(gama)=',len(gama),' maxi=',maxi)
#    print ('noeud ', noeud)
    return arc

def difference_cand_gamma(matE, cand, arc):
    
    ens_res = cand.difference( gamma(matE, arc) )
    return ens_res
        

def clique(matE, l=[]):
    '''
     code bon => tester
     A : ensemble des arcs
    '''
    print ('#####################3')
    Q = set()
    A = set( matE.columns.tolist() )

    Q = expand(matE, A, A, Q, l )
    return l

def find_clique(matE, A, l=[]):
    Q = set()
    Q = expand(matE, A, A, Q, l)
    if len(l) != 0:
        #print("cliques trouves")
        pass
    return l
    
def expand(matE, subg, cand, Q, l ):
    '''
     code bon => tester
    '''
    if len( subg ) == 0:
        #print (' clique ')
        if Q not in l:
            l.append(Q)
        return Q
    else:
        u = maximize_cand_inter_gamma( matE, cand )
        ext_u = difference_cand_gamma( matE, cand, u)
        #print('u = ',u ,' ext_u = ', ext_u)
        while( ext_u ):
            q = ext_u.pop()
            Q.add(q)
            #print (' q =', q, 'Q = ',Q)
            cand_q = cand.intersection( gamma(matE, q) )
            subg_q = subg.intersection( gamma(matE, q) )
            #print ('cand_q =', cand_q, ' subg_q =', subg_q)
                                    
            Q = expand( matE, subg_q, cand_q, Q, l)
                        
            cand_q = cand_q.difference( {q} ) # cand_q = cand_q - set([q])
            #print('1=> q = ',q, ' avant diff Q= ',Q )
            Q = Q.difference( {q} )
            #print('2=> cand_q = ', cand_q, ' FINI = ',FINI, ' Q =',Q)
            #print ( 'back,' )
        return Q

def clique_max(l):
    first_max = len(l[0])
    ens_max = set()
    for elt in l:
        if first_max < len(elt):
            first_max = len(elt)
            ens_max = elt
    return ens_max
    
def is_clique(listeAretes, set_noeuds):
    """
    determine if set_noeuds est une clique
    
    return True s'il est une clique, False sinon
    
    test:
        l = [(1,2),(1,3),(1,4),(2,3),(2,4),(3,4),(2,5),(4,5)]
        set_l = {1,2,3,4}
        is_clique(l,set_l)
    """ 
    l_noeuds = list(set_noeuds)
    for i in range(len(l_noeuds)):
        for j in range(i+1, len(l_noeuds)):
            if (l_noeuds[i],l_noeuds[j]) not in listeAretes and \
                (l_noeuds[j],l_noeuds[i]) not in listeAretes:
                    return False
    return True
############################################

if __name__ == "__main__" :
    
    start= time.time()
    
    location_file = 'data/matrices/matrice_adjacence_grapheDeBase_bon.csv'
    matE = pd.read_csv(location_file, index_col = "arcs")     
    print ('matE = ', matE.columns.tolist() )
    liste_arcs = matE.columns.tolist()
    e = gamma(matE, liste_arcs[0] )
    print ('e = ', e)
    cand = set( liste_arcs ) 
    resu = maximize_cand_inter_gamma(matE, cand)
    ens_res = difference_cand_gamma( matE, cand, resu)
    print ('resu = ', resu)
    print ('ens_res = ', ens_res)
    l = []
    print ( 'len(l) =', clique(matE, l) )
    print ('clique max = ', clique_max(l))
    
    li = []
    A = {"x","y","v"}
    li = find_clique(matE, A, [])
    print ("li = ", li)
    print (time.time() - start) 
    
    path = "testClique/graphe_dual.csv" 
    mat = pd.read_csv(path, index_col = "arcs")
    A = {"A","B","C","D"}
    li = find_clique(mat, A, [])
    print ("li = ", li)