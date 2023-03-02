# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 20:19:25 2023

@author: linus
"""

import os
import time
import csv
from Hasse import analyze_DO,check_ORG_LOR
from analyse_class import Analyse
from LP_compartments import get_compartments, get_min_compartments
from reaction_ERC import create_closures
import multiprocessing
import pickle


def iterate_over_database3(path="C:/Users/linus/python/nice_BM2/", get_time=False, get_constr=False, start_index=0,end_index=3000, use_comma=False, exelname=False):
    
   # header = ['ID',"name",'number_species', 'number_reactions','numberInflow', 'termination', 'number_Os', "number_unique_DOs",\
    #          "type1(O)","type2", "type3(DO)","type1_new(O)","type2_new","type3_new(DO)","Os","unique_DOs","avg_ERC_lenght",\
     #             "largest_SAR_is","compartment_lengths","Timer_ERC","Timer_LP","timer_MRC","number_constr"]
    header = ['ID',"name",'number_species', 'number_reactions',"reaction_orders",'numberInflow', 'numberReversible','termination', 'number_SAR_O', "number_SAR_DO_only",\
              "new_reactions","new_reactions_O","s(DO)=s(O)", "s(DO)=s(DO)","DO_is_largest","Os","unique_DOs","avg_ERC_lenght",\
              "largest_SAR_is","max_MRC","comp3","comp4","comp5","comp6","Timer_ERC","Timer_LP","timer_MRC","number_constr","kaleta_O","timer_all"]
    
    
    inflow_memory=[]
#     cores=mp.cpu_count()
    cores=8
    job_queueueue=[]
    global safe_object
    with os.scandir(path) as sbml_folder:
        #pen('test-output/output_table.csv', 'w', newline='') as exel:
        
        current_index=-1
        memory_runtime_error1=[]
        memory_runtime_error2=[]
        for entry in sbml_folder:
            #skips files before reaching start_index
            if start_index-1 >current_index:
                current_index+=1
                continue
            if entry.name.endswith("560.xml"):
                continue
            
            if entry.name.endswith("595.xml"):
                continue
            
            if entry.name.endswith("909.xml"):
                continue
            #builds up to large number of MRC candidates
            if entry.name.endswith("469.xml"):
                continue
            if entry.name.endswith("470.xml"):
                continue
            if entry.name.endswith("471.xml"):
                continue
            if entry.name.endswith("472.xml"):
                continue
            if entry.name.endswith("473.xml"):
                continue
            ### these models have to long variable names and also have just a O
            if entry.name.endswith("908.xml"):
                continue
            if entry.name.endswith("909.xml"):
                continue
            if end_index <=current_index:
                break
#             if j%25==0:
#                 print(j)
            current_index+=1
            #check for file integrity
            if  ".orgml" in entry.name:
                continue
            if (entry.name.endswith(".sbml") or entry.name.endswith(".xml")) and entry.is_file():
                    job_queueueue.append([entry.path,current_index,use_comma])
    if exelname:
        exel_insert=str(exelname)
    else:
        exel_insert="output_table_200"
    if not os.path.exists("test-output/"):
        os.makedirs("test-output/")
    with open('test-output/' + exel_insert + '.csv',     'w', newline='',encoding="utf-8") as exel:            
        
        print(job_queueueue)
        with multiprocessing.Pool(cores) as p:
            safe_object = p.map(process_file, job_queueueue)
        safe_object=[]    
        writer_csv = csv.writer(exel, delimiter=';')
        writer_csv.writerow(header)
        writer_csv2 = csv.DictWriter(exel, header,delimiter=';') 
        for ele in job_queueueue:
            safe_object=process_file(ele)
            if type(safe_object)==list:
                for i in range(len(safe_object)):
                    writer_csv2.writerow(safe_object[i])
            else:
                writer_csv2.writerow(safe_object)
            
    ##    with mp.Pool(cores) as pool:
    ##   # apply the process_input function to each input in the list
    ##        safe_object = pool.map(process_file, job_queueueue)
            
            """ for ele in range(len(safe_object)-1,-1,-1):
            if type(safe_object[ele])==list:
    
                betweee=safe_object[ele].copy()
                safe_object[ele]=safe_object[ele][1]
                safe_object.insert(ele, betweee[0])
                print("change")
       
                writer_csv = csv.writer(exel, delimiter=';')
                writer_csv.writerow(header)
                writer_csv2 = csv.DictWriter(exel, header,delimiter=';')  
                for i in range(len(safe_object)):
            #        if type(safe_object[i])==list:
           #             print("ciao")
            #            print(safe_object[i])
                    if safe_object[i] is None:
                        pass
                    else:
                      #  print(i)
                      #  print(safe_object[i])
                        print(safe_object[i])
                        try:
                            writer_csv2.writerow(safe_object[i])
                        except:
                            pass"""
        
    
    print("error1")
    print(memory_runtime_error1)
    print("error2")
    print(memory_runtime_error2)
    print(inflow_memory)
def timer(start_time):
    pass
def process_file(input_list):   
    
#outputdict saves output for different rows, the headers work as keys
    entry,current_index,use_comma=input_list
    print(current_index)
  # Pickling the object
    time_new=time.time()
    #counter for each sbml file
    inflow_counter=0
    numberReversible=0
    if current_index%50==0:
        print("index="+str(current_index))
    current_SBML=Analyse(str(entry), consider_reverse=True)
    if current_SBML.reaction_network is None:
        print(entry)
        return()

    #current_SBML.print_RN()
#    if "EmptySet" in current_SBML.species:
 #       current_SBML.change_network(remove_species_from_reactions="EmptySet")
    current_SBML.get_species()
    current_SBML.get_name()
    for i in current_SBML.reaction_network:
        if i.always==True:
            inflow_counter+=1
        if "_reverse" in i.defined_name:
            numberReversible+=1
    solution=[]
    solution.append(calculate_RN(current_SBML,entry,inflow_counter=inflow_counter,reverse_counter=numberReversible,use_comma=use_comma))
    
    if inflow_counter>0:
        current_SBML.change_network(change_inflow=True)
        solution.append(calculate_RN(current_SBML,entry,inflow_counter=None, reverse_counter=numberReversible, use_comma=use_comma))
    if numberReversible>0:
        current_SBML=Analyse(str(entry), alternative_reverse=True)
        current_SBML.get_species()
        current_SBML.get_name()
        solution.append(calculate_RN(current_SBML,entry,inflow_counter=inflow_counter, reverse_counter=numberReversible,use_alt_rev=True, use_comma=use_comma))
        
        if inflow_counter>0:
            current_SBML=Analyse(str(entry),alternative_reverse=True)
            current_SBML.change_network(change_inflow=True)
            current_SBML.get_species()
            current_SBML.get_name()
            solution.append(calculate_RN(current_SBML,entry,inflow_counter=None, reverse_counter=numberReversible,use_alt_rev=True, use_comma=use_comma))
    
    if len(solution)==1:
        solution=solution[0]
    return(solution)
def calculate_RN(current_SBML,entry,inflow_counter=None,reverse_counter=None, use_comma=False,use_alt_rev=False):
#call of analysis 
    timer_all=time.time()
    outputdict={}
    outputdict['ID']=os.path.basename(entry)
    current_SBML.SARs(get_time=True,see_constraints=True,solve_parameter="trans",second_optimize=False)
    timeout_counter=0
    outputdict['name']=current_SBML.name
    if inflow_counter==None:
        outputdict['name']=current_SBML.name+str("inflow")
        inflow_counter=0
    if use_alt_rev:
        outputdict['name']=outputdict['name']+str("alt_reverse")
        inflow_counter=0
    if hasattr(current_SBML,"abort_run"):
        outputdict["termination"]="ERC timeout"
        return(outputdict)
    spec_output=set()
    

    ERC_lengths=current_SBML.print_ERC_len()
    if len(current_SBML.reaction_network)==0:
        return outputdict
    outputdict['avg_ERC_lenght']=(sum(ERC_lengths)/len(ERC_lengths))/len(current_SBML.reaction_network)
    outputdict['number_species']=len(current_SBML.species)
    outputdict['number_reactions']=len(current_SBML.reaction_network)
    numbers_array=[0,0,0,0,0,0,0,0,0]
    for ele in current_SBML.reaction_network:
        numbers=len(ele.listOfReactants)
        try:
            numbers_array[numbers]+=1
        except IndexError:
            numbers_array=numbers
            break
    outputdict["reaction_orders"]=numbers_array
    outputdict['numberInflow']=inflow_counter
    outputdict['numberReversible']=reverse_counter
    
    #erc_timer=current_SBML.ERC_dict[current_SBML.reaction_network[0].defined_name].timer
    if use_comma==True:
        outputdict["Timer_ERC"]=str(current_SBML.time_ERC).replace('.',',')
        outputdict["Timer_LP"]=str(current_SBML.time_solve).replace('.',',')
    else:
        outputdict["Timer_ERC"]=current_SBML.time_ERC
        outputdict["Timer_LP"]=current_SBML.time_solve
    outputdict["number_constr"]=current_SBML.number_constr
#                 print(current_SBML.SAR_solution)
    
    if hasattr(current_SBML,"SAR_solution_wrong"):
        outputdict['termination']=current_SBML.SAR_solution_wrong
        timeout_counter=1
        if current_SBML.SAR_solution_wrong=="Solution_timeout":
            outputdict['termination']
#                 if "error" in outputdict:
#                     timeout_counter=1
#                     outputdict['termination']=outputdict["error"]
    if len(current_SBML.SAR_solution)==0 or current_SBML.SAR_solution[0]==([], []):
        outputdict['termination']="empty_solution"
        timeout_counter=1
        print("empty")
    
    if timeout_counter==0:
        
        #analyze_dict= analyze_DO(current_SBML.SAR_solution, current_SBML.reaction_network)
        try:
            analyze_dict= analyze_DO(current_SBML.SAR_solution, current_SBML.reaction_network)
        except KeyError:
            analyze_dict={"unique_DOs":"error","Os":"","unique_DOs":"","unique_DOs":"",}
            pass
        outputdict.update(analyze_dict)
        
        
        
        
        
        o_kaleta=get_species_closure(current_SBML.reaction_network, current_SBML.SAR_solution)
        o_kaleta_new=[]
        
        """for ele in o_kaleta:
            appender=[]
            
            for spec in ele:
                if current_SBML.dicto[spec]!= "":
                    appender.append(current_SBML.dicto[spec])
                else:
                    appender.append(spec)
            o_kaleta_new.append(appender)"""
        
        outputdict["kaleta_O"]=o_kaleta_new
        compartment_counter=[]
        start_MRC=time.time()
        cancel_mrc=False
        get_comp=True
        compartment_3=0
        compartment_4=0
        compartment_5=0
        compartment_6=0
        maxnumbercomp=0
        if not use_alt_rev:
            if get_comp==True:
                for i in range(len(current_SBML.SAR_solution)):
                    if not check_ORG_LOR(current_SBML.SAR_solution[i][0],current_SBML.reaction_network):
                        if current_SBML.SAR_solution[i][0]:
                            compartments=[]
                            solution_gc=get_compartments(current_SBML.reaction_network,current_SBML.SAR_solution[i][0])
                            maxnumbercomp=max(maxnumbercomp,len(solution_gc[1]))
                            #print("shiyy")
                            #print(len(solution_gc[1]))
                            #print(solution_gc[1])
                            compartments=get_min_compartments(*solution_gc)
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
                    
                        if int(time.time()-start_MRC)>180:
                            cancel_mrc=True
                            print("cANCELA")
                            break
                if cancel_mrc:   
                    outputdict["timer_MRC"]="timeout"
                else:
                    if use_comma==True:
                        outputdict["timer_MRC"]=str(time.time()-start_MRC).replace('.',',')
                    else:
                        outputdict["timer_MRC"]=(time.time()-start_MRC)
        
        
            outputdict["max_MRC"]=maxnumbercomp
            outputdict["comp3"]=compartment_3
            outputdict["comp4"]=compartment_4
            outputdict["comp5"]=compartment_5
            outputdict["comp6"]=compartment_6
        else:
        #outputdict["compartment_lengths"]=compartment_counter
            outputdict["max_MRC"]=""
            outputdict["comp3"]=""
            outputdict["comp4"]=""
            outputdict["comp5"]=""
            outputdict["comp6"]=""
        if len(str(outputdict["unique_DOs"]))>32755:
            outputdict["unique_DOs"]=str(outputdict["unique_DOs"])[:32755]
    #         outputdict["unique_DOs"]=str(outputdict["unique_DOs"])[:32760]
        if len(str(outputdict["Os"]))>32755:
            outputdict["Os"]=str(outputdict["Os"])[:32755]
        #print(outputdict)
    else:
        try:
            if hasattr(current_SBML, "info_so_far"):
                outputdict["number_SAR_O"]=current_SBML.info_so_far[0]
                outputdict["number_SAR_DO_only"]=current_SBML.info_so_far[1]
        except KeyError:
            pass
    outputdict["timer_all"]=str(time.time()-timer_all).replace('.',',')
    return outputdict
    
def get_species_closure(loR,list_of_SARs):
    o_set=set()
    
    for SAR in list_of_SARs:
        if check_ORG_LOR(SAR[0],loR):
            species_set=set()
            for reaction in loR:
                if reaction.defined_name in SAR[0]:
                    species_set.update(reaction.listOfReactants)
                    species_set.update(reaction.listOfProducts)
            species_set_frozen=frozenset(species_set)
            o_set.add(species_set_frozen)
    return(o_set)

def calculate_number_DO(listofSpecies,SARs,LOR):
    species_set_SAR=set()
    for rea in SARs[0][0]:
        for reaction in LOR:
            if reaction.defined_name==rea:
                species_set_SAR.add(set(reaction.listOfProducts)+set(reaction.listOfReactants))
    free_spec=set()
    for spec in listofSpecies:
        if spec not in species_set_SAR:
            free_spec.add(spec)
    ERC_dict1= create_closures(LOR, species=free_spec, Timer_in_sec=120)

    for ele in ERC_dict1:
        if len(ele!=0):
            for solution in SARs:
                pass


