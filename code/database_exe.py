# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 20:48:14 2023

@author: linus
"""
from iterate_over_DB_multi import iterate_over_database3
if __name__ == '__main__':
    safe_object=[]
    
    iterate_over_database3(path="C:/Users/linus/python/BioModels-files/",use_comma=True,exelname="biomodels_pulp_multi22",end_index=1000)
    #iterate_over_database3(path="C:/Users/linus/python/BioModels-files/",use_comma=True,exelname="Kaleta_1_0_0_new",end_index=184)

    #iterate_over_database3(path="C:/Users/linus/python/BioModels-files/",use_comma=True,exelname="test_new3",end_index=200)
    print(safe_object)
