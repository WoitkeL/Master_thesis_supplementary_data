# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 17:07:39 2023

@author: linus
"""
import pulp
import time
import regex as re
from Hasse import check_ORG_LOR
from reaction_ERC import create_closures, ERC_meets_transitivity
#import gurobipy as gp
from gurobipy import GRB
import gurobipy
import os
import sys
from reaction_ERC import Reaction
def Sorting(lst):
    lst2 = sorted(lst, key=len)
    return lst2
#def setup_LP_Var(lOR_solve, job="SAR", solve_parameter="",second_optimize="Default" ,only_largest_solution=False, see_constraints=False, get_time=False, use_new=False, usenew_length=0, time_restart=False): 
def setup_LP_Var(lOR_solve_input, job="SAR", solve_parameter="",second_optimize="Default" ,only_largest_solution=False, \
                 see_constraints=False, get_time=False, new_loop=False,try_loop=False, use_gurobi=False):
    lOR_solve=lOR_solve_input.copy()
    noHasse=True
    ERC_dict={}
    ERC_dict1={}
    speciesdict={}
    setOfSpecies=set()
    outputdict={}
    reverse_map={}
    
    alternative_reverse=False
    if get_time==True:
        time_ERC_start=time.process_time()
    if not new_loop:
        ERC_dict= create_closures(lOR_solve, Timer_in_sec=120)
        #[print(key,':',[ele.defined_name for ele in ERC_dict[key].reactions ]) for key in ERC_dict.keys()]
        ERC_dict=ERC_meets_transitivity(ERC_dict)
        #[print(key,':',[ele.defined_name for ele in ERC_dict[key].reactions ]) for key in ERC_dict.keys()]
        if type(ERC_dict)==str:
            print("timeout ERC")
            outputdict["error"]="timeout_ERC"
            return(outputdict)
    else:
        ERC_dict= new_loop[2]
#generates species set
    for ele in lOR_solve:
        setOfSpecies.update(ele.listOfReactants)  
        setOfSpecies.update(ele.listOfProducts)
        
        if ele.reversible==True:
                alternative_reverse=True
                lOR_solve.append(Reaction(ele.defined_name+str("_reverse"),ele.listOfProducts,ele.listOfReactants,ele.prod_stoich,ele.reac_stoich))
 ##########         #######               #lOR_solve[-1].reverse=True
                reverse_map[ele.defined_name]=ele.defined_name+str("_reverse")
                #alternative_reverse
                #solve_dist_org+=Rbool[ele]+Rbool[ele]<=1
    if get_time==True:
        time_ERC=time.process_time()-time_ERC_start
##########################################################################################################################
#declaration of the required dictionaries, which link names of elements of the reaction environment to the pulp objects
#also defines the pulp variables and their range of numbers

#generates the reaction names of the list of reactions
    listOfReactionsname = [i.defined_name for i in lOR_solve]
#the exact stoichiometric parameters for reactions and the resulting sspecies quantity

    reac=pulp.LpVariable.dicts("r",listOfReactionsname,lowBound=0) #
    spec=pulp.LpVariable.dicts("s",setOfSpecies,lowBound=0)
#the boolean value for reactions
    Rbool=pulp.LpVariable.dicts("Rbool",listOfReactionsname,lowBound=0, upBound=1,cat="Integer")
    
    
#     if not bool(re.search('SAR', job)): new
    if job=='DO':
        species_exist=pulp.LpVariable.dicts("species_exist",setOfSpecies,lowBound=0, upBound=1,cat="Integer")
#these values are the boolean value for the overproduction of a species, only used when looking for SARs
    if job=='SAR' and bool(second_optimize)==True:
        overprod=pulp.LpVariable.dicts("overprod",setOfSpecies,lowBound=0, upBound=1,cat="Integer")
        spec=pulp.LpVariable.dicts("s",setOfSpecies,lowBound=0)
    
    
    
    
    
    #these values represent the boolean value of existence of a species in a DO and is therefore only used when looking for DOs 
#     species_exist=pulp.LpVariable.dicts("species_exist",setOfSpecies,lowBound=0, upBound=1,cat="Integer")
    
#defines optimization goal for LP problem
#always maximization problem except when looking for minimal SAR to sustain a given reaction
    if job=="minSAR":
        solve_dist_org=pulp.LpProblem("Solver Dist Org",pulp.LpMinimize)
    else:
        solve_dist_org=pulp.LpProblem("Solver Dist Org",pulp.LpMaximize)
   
     #optimization in regard to command given
    
    #looping for all solution:
    #order giving parameter has larger influence, then criteria to find specific DO/SAR
    if job=="DO":
        objective_function=pulp.lpSum(species_exist.values())+pulp.lpSum(Rbool.values())/1000
    else:
        if bool(second_optimize)==True:
            objective_function=pulp.lpSum(Rbool.values())+pulp.lpSum(overprod.values())/1000
        else:
            objective_function=pulp.lpSum(Rbool.values())
#     if job=="loop_DO" or job=="DO":
#         objective_function=pulp.lpSum(species_exist.values())+pulp.lpSum(Rbool.values())/1000
#     if job=="loop_SAR" or job=="max_SAR":
#         objective_function=pulp.lpSum(Rbool.values())+pulp.lpSum(overprod.values())/1000


                         
    
    #minSAR shows smallest solution for the given reaction
    #the reaction is set active 
    if job=="minSAR":
        try:
            solve_dist_org+=Rbool[solve_parameter]==1
            objective_function=pulp.lpSum(Rbool.values())-pulp.lpSum(overprod.values())/1000
        except KeyError: 
            print("given reaction does not exist")
            return
    
    #the objective function is then implented in pulp
    solve_dist_org+=objective_function    
    
##########################################################################################################################        
#managing contraints regarding species:
#links species to their reactions in regard to being product or support of it
    for species in setOfSpecies:
            speciesdict[species+'_r']=[]
            speciesdict[species+'_p']=[]  

    for reaction in lOR_solve:
        for i in range(len(reaction.listOfReactants)):
            speciesdict[reaction.listOfReactants[i]+'_r'].append((reaction.defined_name,reaction.reac_stoich[i]))
        for i in range(len(reaction.listOfProducts)):
            speciesdict[reaction.listOfProducts[i]+'_p'].append((reaction.defined_name,reaction.prod_stoich[i]))    
#the stoichiometric matrix is implemented as equations for each species as the sum of consumption and 
#production over all reactions
    
    for ele in setOfSpecies:
        educt=[]
        product=[]
        for combi in speciesdict[ele+'_r']:
            educt.append(int(combi[1])*reac[combi[0]])
        for combi in speciesdict[ele+'_p']:
            product.append(int(combi[1])*reac[combi[0]])
         #lpSum groups list of reactions to one expression
        educt_sum=pulp.lpSum(educt)
        product_sum=pulp.lpSum(product)
    #the equations for the changes of species are implemented
        solve_dist_org+= 0<=product_sum-educt_sum 
    
    #if bool(second_optimize)==True:
    wrong=False
    if wrong:
        for key, element in spec.items():
            educt=[]
            product=[]
            
            for combi in speciesdict[key+'_r']:
                educt.append(int(combi[1])*reac[combi[0]])
            for combi in speciesdict[key+'_p']:
                product.append(int(combi[1])*reac[combi[0]])
        
           
        #lpSum groups list of reactions to one expression
            educt_sum=pulp.lpSum(educt)
            product_sum=pulp.lpSum(product)
        #the equations for the changes of species are implemented
            solve_dist_org+= element==product_sum-educt_sum 
        #overproduction is inserted, if we look for SARs
        #         if job=="minSAR" or job=="loop_SAR" or job=="get_SAR":
            if job=='SAR' and bool(second_optimize)==True:
                solve_dist_org+= overprod[key]<=element                      
                solve_dist_org+= element<=overprod[key]*1000  
        
#########################################################################################################################
#managing contraints regarding reactions:
   # if filter_reversible:
    #    extra_var=[]
     #   for reaction1 in lOR_solve:
      #      if reaction1.reversible==True:
       #         extra_var.append(reaction1.defined_name)
                
        
    for ele in lOR_solve:
        if ele.reversible==True:
            if "_reverse" not in ele.defined_name:
                solve_dist_org+=Rbool[ele.defined_name]+Rbool[reverse_map[ele.defined_name]]<=1
    
        
    #direction=pulp.LpVariable.dicts("direction",extra_var,  -1, 1,cat="Binary")
    #the relationship between the boolean reaction value and the real value are set up    
    for key, element in reac.items():
        solve_dist_org+=Rbool[key]<=element
        solve_dist_org+=element<=Rbool[key]*1000
    
   
    if type(ERC_dict)==str:
        print("timeout ERC")
        outputdict["error"]="timeout_ERC"
        return(outputdict)
        
    
    for reaction1 in lOR_solve:
    #ERCs are implemented as pulp constraints. the dependance of all reactions to the first entry of the ERC is 
    #managed by a simple loop over all elements besides the first element. since the dependance is one sided we use
    #element<= first element
    #closure of reactions
  #      if "_reverse" in reaction1.defined_name:
   #         
    #        cut_len=len("_reverse")
     #       name_without_rev=reaction1.defined_name[:-cut_len]
      #      shit_fuck.append(name_without_rev)
       #     for reaction2 in lOR_solve:
        #        if reaction2.defined_name==name_without_rev:
         #           #solve_dist_org+= Rbool[reaction1.defined_name]+Rbool[reaction2.defined_name]<=1
          #          #geht nicht
           #         #solve_dist_org+=reac[reaction1.defined_name]-reac[reaction2.defined_name]!=0
            #        a=pulp.LpVariable(str(name_without_rev)+"decider",lowBound=0, upBound=1,cat="Integer")
             #       m=100
              #      solve_dist_org+=reac[reaction1.defined_name]-reac[reaction2.defined_name]+m*a>=1
               #     solve_dist_org+=reac[reaction1.defined_name]-reac[reaction2.defined_name]+m*a<=m-1
                #    
                 #   print("no_rev")
                  #  print(reaction1.defined_name)
                   # print(reaction2.defined_name)
        #if "_reverse" not in ele.defined_name and alternative_reverse:
        if alternative_reverse:
            if "_reverse" not in reaction1.defined_name:
                if reaction1.reversible==True:
                    term1= Rbool[reaction1.defined_name]+Rbool[reverse_map[reaction1.defined_name]]
                else:
                    term1= Rbool[reaction1.defined_name]
                for x in range (len(ERC_dict[reaction1.defined_name].reactions)-1):
                    
                        if ERC_dict[reaction1.defined_name].reactions[x+1].reversible==True:
                            term2= Rbool[ERC_dict[reaction1.defined_name].reactions[x+1].defined_name]+Rbool[reverse_map[ERC_dict[reaction1.defined_name].reactions[x+1].defined_name]]
                        else:
                            term2=Rbool[ERC_dict[reaction1.defined_name].reactions[x+1].defined_name]
                            #solve_dist_org+= Rbool[reaction1.defined_name]<= Rbool[ERC_dict[reaction1.defined_name].reactions[x+1].defined_name]+Rbool[reverse_map[ERC_dict[reaction1.defined_name].reactions[x+1].defined_name]]
                            
                        
                        #else:
                            #solve_dist_org+= Rbool[reaction1.defined_name]<= Rbool[ERC_dict[reaction1.defined_name].reactions[x+1].defined_name]
                        solve_dist_org+=term1<=term2
          #  else:
           #     for x in range (len(ERC_dict[reaction1.defined_name].reactions)-1):
            #        solve_dist_org+= Rbool[reaction1.defined_name]<= Rbool[ERC_dict[reaction1.defined_name].reactions[x+1].defined_name]      
        else:
            for x in range (len(ERC_dict[reaction1.defined_name].reactions)-1):
                solve_dist_org+= Rbool[reaction1.defined_name]<= Rbool[ERC_dict[reaction1.defined_name].reactions[x+1].defined_name]      
                       # Rbool[ERC_dict[reaction1.defined_name].reactions[x+1].defined_name]
                       # Rbool[ele]+Rbool[reverse_map[ele]]
                   
            #test of reversibility
#            for x in range (len(ERC_dict[reaction1.defined_name].reactions)-1):
 #               
  #              if ERC_dict[reaction1.defined_name].reactions[x+1].defined_name!=name_without_rev:
   #                 solve_dist_org+= Rbool[reaction1.defined_name]<= Rbool[ERC_dict[reaction1.defined_name].reactions[x+1].defined_name]
    #    else:
#
 #           for x in range (len(ERC_dict[reaction1.defined_name].reactions)-1):
  #              
   #             if ERC_dict[reaction1.defined_name].reactions[x+1].defined_name
    #            solve_dist_org+= Rbool[reaction1.defined_name]<= Rbool[ERC_dict[reaction1.defined_name].reactions[x+1].defined_name]
           
            
        #print(solve_dist_org)  
            
            #solve_dist_org+= reac[reaction1.defined_name]/100<= reac[ERC_dict[reaction1.defined_name].reactions[x+1].defined_name]
#########################################################################################################################
    #closure from each species of the DO 
        if job=="DO":   
    #new ERC initialization with the alternative input of solospecies
    #same loop as used in first ERCs, but now the set of reaction is not dependent to first reaction of the ERC
    #but the solospecies used in the class initialization
            ERC_dict1= create_closures(lOR_solve, species=setOfSpecies, Timer_in_sec=120)
            if type(ERC_dict1)==str:
                print("timeout ERC")
                outputdict["error"]="timeout_ERC"
                return(outputdict)
            for species in setOfSpecies:
                for x in range (len(ERC_dict1[species].reactions)):
                    solve_dist_org+= species_exist[species]<= Rbool[ERC_dict1[species].reactions[x].defined_name]
                    
    # also sets dependance of species_exist to the occurence in reactions. Since if a reaction is active, the
    #involved species are part of the DO
            for reaction2 in lOR_solve:
                for species in reaction2.listOfReactants:
                    solve_dist_org+= species_exist[species]>= Rbool[reaction2.defined_name]
                for species in reaction2.listOfProducts:
                    solve_dist_org+= species_exist[species]>= Rbool[reaction2.defined_name]
                    
    #input reaction active for all LP except when solving for minimal circle, to better showcase their excluded behavior
        if job!="minSAR":
            if reaction1.always:
                solve_dist_org+= Rbool[reaction1.defined_name]==1
                
                
                
    #a reaction which is not closed in regard to the species set has to be inactive
        if not reaction1.closed:
            solve_dist_org+= Rbool[reaction1.defined_name]==0
            
    #command to see constraints of the LP 
    if see_constraints==True:
 #       print(solve_dist_org)
         pass
    if new_loop:
        solve_dist_org+=pulp.lpSum(Rbool.values())<=new_loop[0]
        timeout_start=new_loop[1]
    else:
        timeout_start = time.time()
 #######################################################################################################################
    #specific commands used for getting desired output
    #LOOP
    #implemented timer to kill Lp after specific time
    timeout = 180   # [seconds]
    output_list=[]
    #check for all solutions
    i1=0
    new_var=True
    outputlist=[]
    directory = "gurobi_solutions"

    max_element=[]
    current_solution=[]
    i1=i1+1



    if use_gurobi==True:
        path_to_gurobi = "C:/gurobi1000/win64/bin/gurobi_cl.exe"
        sol_pool_var=3000
        time1=time.process_time()
        solve_dist_org.writeLP("problem.lp")
        gurobi_model = gurobipy.read("problem.lp")
        gurobi_model.params.TimeLimit = 300
        gurobi_model.params.PoolSearchMode = 2
        gurobi_model.params.Threads = 8
        gurobi_model.params.LogToConsole  = 0
        gurobi_model.params.LogFile=""
        gurobi_model.setParam(GRB.Param.PoolSearchMode, 2)
        gurobi_model.setParam(GRB.Param.OutputFlag, 0)
        gurobi_model.params.PoolSolutions = int(sol_pool_var)
        result= gurobi_model.optimize()
        # Get the time taken to solve the LP problem
        solve_time = gurobi_model.Runtime
   
        print("Time taken to solve LP:", solve_time, "seconds")
        if gurobi_model.status == gurobipy.GRB.TIME_LIMIT:
            print("Time limit reached")

        if gurobi_model.status == gurobipy.GRB.OPTIMAL:
            print("The model has an optimal solution with objective value", gurobi_model.objVal)
        print(gurobi_model.status) 
        #sys.stdout.flush()
        pool_solutions = gurobi_model.SolCount
        #print(pool_solutions)
        #pool_solutions = gurobi_model.getAttr("PoolSolutions")
        if job=="DO":
            sort_para="species_exist"
        else:
            sort_para="Rbool"
        pattern1 = re.compile(sort_para)
        for i in range(pool_solutions):
            current_solution=[]
            gurobi_model.setParam(gurobipy.GRB.Param.SolutionNumber, i)
            gurobi_model.setParam(GRB.Param.SolutionNumber, i)
            for v in gurobi_model.getVars():
                if pattern1.match(v.VarName) and int(v.Xn)==1:
            
                        
                    current_solution.append(v.VarName[(len(sort_para)+1):])
            output_list.append([current_solution,[]])

    
        
    else:
        obj_val_dict={}
        while True:
            
            max_element=[]
            current_solution=[]
            i=i+1
    #check to break if timer exceeded
            if time.time() > timeout + timeout_start:
            #if time.process_time()+time_LP > timeout + timeout_start:
                print("timeout")
                outputdict["error"]="Solution_timeout"
                outputdict["number_SAR_O"]=int(0)
                outputdict["number_SAR_DO_only"]=int(0)
                for knot_reactions in output_list: 
                                
                                
                    if check_ORG_LOR(knot_reactions[0],lOR_solve):
                        outputdict["number_SAR_O"]+=1
                    else:
                        outputdict["number_SAR_DO_only"]+=1
                break
    
    #solve LP
            solve_dist_org.solve()

    
    #check if solution is valid
            
            if (pulp.LpStatus[solve_dist_org.status]=="Infeasible") or (pulp.LpStatus[solve_dist_org.status]=="Undefined"):
                if pulp.LpStatus[solve_dist_org.status]=="Undefined":
                    print("satus is undefined.")
                    print("an alternate solution may exist")
                    outputdict["error"]="Solution_Undefined"
                break   
            
    #both Loop checks differ only in objective and output variables
    #depending on the case, these are inserted here
            if job=="DO":
                sort_para="species_exist"
                max_para="Rbool"
    
            else:
                sort_para="Rbool"
                max_para="overprod"
            pattern1 = re.compile(sort_para)
            pattern2 = re.compile(max_para)
    #all boolean variables, which are 1 and match the pattern are extracted  
            for variable in solve_dist_org.variables():
                if pattern1.match(variable.name) and variable.varValue==1:
                    current_solution.append(variable.name[(len(sort_para)+1):])
                if pattern2.match(variable.name) and variable.varValue==1:
                    max_element.append(variable.name[(len(max_para)+1):])
    #to see progress in loop   
            if i%20==0:
                print(i)
            if try_loop==True:
                if (i>20 and len(output_list[-1][0])>len(current_solution)):
                    new_loop=[len(current_solution),timeout_start,ERC_dict]
                    intermediate_output=setup_LP_Var(lOR_solve, job=job, solve_parameter=solve_parameter,second_optimize=second_optimize, \
                        see_constraints=see_constraints,get_time=get_time, new_loop=new_loop, try_loop=try_loop)
                    output_list.extend(intermediate_output["output_list"])
                    try:
                        if intermediate_output["error"]=="Solution_timeout":
                            outputdict["error"]="Solution_timeout"
                    except:
                        pass
                    return()
    #append solution to output object
            output_list.append((current_solution, max_element))
            locs = locals()
    
            if only_largest_solution==True:
                break
    #solution is excluded by integer cut
            solve_dist_org+=pulp.lpSum([(eval(eval("sort_para",locs),locs)[i]) for i in current_solution]) <=len(current_solution)-1


    #print("Total end")
    #output_list.sort(key=lambda x: len(x[0]))
    #output_list.reverse()
    outputdict["output_list"]=output_list
    outputdict["ERC_dict"]=ERC_dict
    if len(output_list)==0:
        outputdict["error"]="no_solution"
#                             outputdict['termination']="empty solution"

    if see_constraints==True:
        outputdict["constraints"]=len(solve_dist_org.constraints)
    if get_time==True:
        outputdict["time_ERC"]=time_ERC
    #print(len(lOR_solve))
    if use_gurobi==False:
        solve_time=time.time()-timeout_start
    outputdict["time_solve"]=str(solve_time)
    
    return(outputdict)
######################################################################################################################       
