# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 16:40:35 2023

@author: linus
"""

from analyse_class import Analyse
from Hasse import analyze_DO,check_ORG_LOR,check_ORG_LOR_expanded
from reaction_ERC import Reaction
import time 
from LP_compartments import get_compartments, get_min_compartments
from iterate_over_DB_multi import get_species_closure
time1=time.time()
#example1=Analyse("C:/Users/linus/Downloads/ntop.sbml")
#x=Reaction("r9",["d1","b1"],["b1","e2"],[1,1],[2,1])
#example1.change_network(add_reaction=x)
#example1=Analyse("C:/Users/linus/python/kaleta_old_files/BIOMD0000000072.xml",consider_reverse=True, consider_constant=True, consider_init_ammount=True)
#example1=Analyse("C:/Users/linus/python/Biomodels-files/BIOMD0000000595.xml")
#example1=Analyse()
example1=Analyse("C:/Users/linus/python/Biomodels-files/BIOMD0000000862.xml",consider_reverse=True,alternative_reverse=False)
#extracting name from sbml-file
#print(example1.dicto)
example1.print_RN()
#for i in example1.reaction_network:
 #   i.reversible=False
"""example1=Analyse("C:/Users/linus/python/Biomodels-files/BIOMD0000000009.xml",alternative_reverse=False)
"""

for i in example1.reaction_network:
    if i.reversible==True:
        print(i.defined_name)
example1.print_ERC()
#example1.change_network(remove_species_from_reactions="EmptySet")
example1.get_species()
print("wolf")
print(example1.species)
#example1.change_network(change_inflow=True)
#calculating DOS and SARs respectively
xax=example1.print_ERC_len()
print()
print(sum(xax))
print(xax)
print("shit")
#example1.change_network(change_inflow=True)
#example1.change_network(change_inflow=True)
#example1.SARs()
#example1.SARs(solve_parameter="trans",second_optimize=True)
#example1.DOs()
example1.SARs(second_optimize=True,use_gurobi=False)
print("nod")
#print(example1.SAR_solution)
#print(len(example1.SAR_solution))
sache=[]
sache_len=[]
#for i in example1.SAR_solution:
#    
#    bol,boos=check_ORG_LOR_expanded(i[0],example1.reaction_network)
#    if not bol:
#        sache_len.append(len(i))
##        solution_gc=get_compartments(example1.reaction_network,i[0])
 #       sache.append(len(solution_gc[1]))
        #compartments=get_min_compartments(*solution_gc)
 #   print(bol)
#print(example1.SAR_solution)   
#print(example1.SAR_solution)
print("hier")
print()
#print(example1.DO_solution)
#example1.draw_hasse(solution="DOs",use_naive=False,second_value=True,show_compartments=False)
#example1.draw_hasse(solution="DOs",use_naive=True,show_compartments=False, second_value=False)
example1.draw_hasse(solution="SARs",second_value=False, use_naive=False,show_compartments=True, show_new=False,shortform=True)
print("öff")
#print(len(example1.SAR_solution))
olaf=get_species_closure(example1.reaction_network,example1.SAR_solution)
#print(example1.SAR_solution)
#analyze_dict= analyze_DO(example1.SAR_solution, example1.reaction_network)
#print((olaf))
#print(len((olaf)))
print(time.time()-time1)
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
# Erstellen einer Liste mit zufälligen Werten
data = sache
print("what?")
#print(sache)
print()
current_SBML=example1
compartment_3=0
compartment_4=0
compartment_5=0
compartment_6=0
for i in range(len(current_SBML.SAR_solution)):
    if not check_ORG_LOR(current_SBML.SAR_solution[i][0],current_SBML.reaction_network):
        if current_SBML.SAR_solution[i][0]:
            compartments=[]
            solution_gc=get_compartments(current_SBML.reaction_network,current_SBML.SAR_solution[i][0])
            #maxnumbercomp=max(maxnumbercomp,len(solution_gc[1]))
            #print("shiyy")
            #print(len(solution_gc[1]))
            #print(solution_gc[1])
            #print("cheese")
            #print(i)
            #print(current_SBML.SAR_solution[i])
            #print("MRCS:")
            #print(solution_gc[1])

            compartments=get_min_compartments(*solution_gc)
            #print("COMPS")
            #print(compartments)
        else: compartments=[]
        #if len(compartments)>2:
         #   compartment_counter.append([i,len(compartments)])
        if len(compartments)==3:
            compartment_3+=1
        if len(compartments)==4:
            compartment_4+=1
        if len(compartments)==5:
            compartment_5+=1
        if len(compartments)>5:
            compartment_6+=1
print(compartment_3)
print(compartment_4)
print(compartment_5)
print(compartment_6)
print("yello")
# Berechnen der kumulativen Häufigkeit
cumulative_freq = np.cumsum(np.unique(data, return_counts=True)[1]) / len(data)
cumulative_freq = 1 - cumulative_freq
# Erstellen des Seaborn-Plots
sns.set()
sns.lineplot(x=np.unique(data), y=cumulative_freq)
#plt.yscale("log")
plt.xscale("log")
print()
# Anzeigen des Plots
plt.show()
