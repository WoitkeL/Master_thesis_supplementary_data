# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 17:14:02 2023

@author: linus
"""

"""HASSE"""
"""HASSE"""
"""HASSE"""

from LP_compartments import get_min_compartments,get_compartments
import regex as re
from graphviz import Digraph, Source
import copy
#function to check a subset of reactions for closure and therefore the attribute of being DO lvl 1 (being a O)
def check_ORG_LOR(knot_reactions,LOR):
    LOR_dict={}
    alternate_reverse=False
    for reac in LOR:
        LOR_dict[reac.defined_name]=set(reac.listOfReactants+ reac.listOfProducts)   
        if reac.reversible==True:
            alternate_reverse=True
    speciesset=set()
    
    #generate species set
    knot_reactions_remover=[]
    for reaction in knot_reactions:
        if "_reverse" in reaction and alternate_reverse:
            knot_reactions.append(reaction[:-8])
            speciesset.update(LOR_dict[reaction[:-8]])
            knot_reactions_remover.append(knot_reactions)
        else:
            speciesset.update(LOR_dict[reaction])
    if len(knot_reactions_remover)!=0:
        knot_reactions = [knot for knot in knot_reactions if knot not in knot_reactions_remover]    
    #check for support of inactive reactions
    for reaction in LOR:
        if reaction.defined_name not in knot_reactions and set(reaction.listOfReactants).issubset(speciesset):
            return(False)
        if reaction.reversible==True:
            if reaction.defined_name not in knot_reactions and set(reaction.listOfProducts).issubset(speciesset):
                return(False)
    return(True)

def check_ORG_LOR_expanded(knot_reactions,LOR):
    LOR_dict={}
    alternate_reverse=False
    for reac in LOR:
        LOR_dict[reac.defined_name]=set(reac.listOfReactants+ reac.listOfProducts) 
        if reac.reversible==True:
            alternate_reverse=True
    speciesset=set()
    
    #generate species set
    knot_reactions_remover=[]
    for reaction in knot_reactions:
        if "_reverse" in reaction and alternate_reverse:
            knot_reactions.append(reaction[:-8])
            speciesset.update(LOR_dict[reaction[:-8]])
            knot_reactions_remover.append(knot_reactions)
        else:
            speciesset.update(LOR_dict[reaction])
    if len(knot_reactions_remover)!=0:
        knot_reactions = [knot for knot in knot_reactions if knot not in knot_reactions_remover]
    #check for support of inactive reactions
    for reaction in LOR:
        if reaction.defined_name not in knot_reactions and set(reaction.listOfReactants).issubset(speciesset):
            return(False,speciesset)
            if reaction.reversible==True:
                if reaction.defined_name not in knot_reactions and set(reaction.listOfProducts).issubset(speciesset):
                    return(False,speciesset)
    return(True,speciesset)

def analyze_DO(x,y):
    niceDOs=[]
    niceOs=[]
    state_dict={}
    check_if_Org_dict={}
    outputdict={}
    #print(x)
    #get dictionaries for new active reactions and new active overproduction
    print()
    alternate_reverse=False
    for reaction in y:
        if reaction.reversible==True:
            alternate_reverse=True
    #dict1 conatains all reactions, which did not appear in subsets of the SAR
    if alternate_reverse:
        dict1,dict2= create_Hasse(x,y,analyze_only=True)
    else:
        dict1,dict2= create_Hasse3(x,y,analyze_only=True)
    
    
#    for i in range(len(x)):
 #       add_new_counter=0
  #      if bool(dict1[i])==True:
   #         add_new_counter=3
    #    #check if DO is also ORG
     #   #bool_result=check_ORG_LOR_expanded(x[i][0],y)
      #  bool_result,speciesout=check_ORG_LOR_expanded(x[i][0],y)
       # speciesout=frozenset(speciesout)
        #if bool_result==True:
         #   niceOs.append(x[i][0])
          #  state_dict[i]=1+add_new_counter
           # check_if_Org_dict[speciesout]=i
#        else:
 #           #safe DOs which are not Os
  #          niceDOs.append(x[i][0])
   #         if speciesout in check_if_Org_dict:
    #            state_dict[i]=2+add_new_counter
     #           
      #      else:
       #         state_dict[i]=3+add_new_counter
#                 for key in check_if_Org_dict:
#                     if speciesout.issubset(key):
#                         state_dict[i]=3+add_new_counter
#                         break
#                 if state_dict[i]=
    new_reactions=0
    new_reactions_O=0
    equal_DO_DO=[0,0]
    equal_DO=[0,0]
    smaller_DO=[0,0]
    DO_large=0
    new_reactions_O=0
    check_if_spec_new={}
    check_if_Dist_Org_dict={}
    for i in range(len(x)):
        add_new_counter=0
        new_reac_bool=0
        
        if bool(dict1[i])==True:
            new_reac_bool=1
            
            
        #check if DO is also ORG
        #bool_result=check_ORG_LOR_expanded(x[i][0],y)
        bool_result,speciesout=check_ORG_LOR_expanded(x[i][0],y)
        speciesout=frozenset(speciesout)
        #for every org
        if i==0 and bool_result==False:
            DO_large=1
        if bool_result==True:
            
            niceOs.append(x[i][0])
            if new_reac_bool:
                new_reactions_O+=1
            #saves species of sar belonging to org as key, and 
            #value determines if sar of org has new rea
            if bool(dict1[i])==True:
                check_if_Org_dict[speciesout]=1
            else:
                check_if_Org_dict[speciesout]=0
        else:
            if new_reac_bool:
                new_reactions+=1
            #safe DOs which are not Os
            niceDOs.append(x[i][0])
            if speciesout in check_if_Dist_Org_dict:
                if check_if_Dist_Org_dict[speciesout]:
                    equal_DO_DO[1]+=1
                else:
                    equal_DO_DO[0]+=1
            #saves species of sar belonging to org as key and 
            #value determines if sar of org has new rea
            if bool(dict1[i])==True:
                check_if_Dist_Org_dict[speciesout]=1
            else:
                check_if_Dist_Org_dict[speciesout]=0
            
            #check for sar of DO, if the species corresponds to SAR for DO or O respectively
            if speciesout in check_if_Org_dict.keys():
                if check_if_Org_dict[speciesout]:
                    equal_DO[1]+=1
                else:
                    equal_DO[0]+=1
            else:
                for species in check_if_Org_dict.keys():
                    if species in speciesout: 
                        if check_if_Org_dict[speciesout]:
                            smaller_DO[1]+=1
                        else:
                            smaller_DO[0]+=1
        
        
                    
        if bool(dict1[i])==True:
            if bool_result==True:
                check_if_spec_new[speciesout]="O"
            else:
                check_if_spec_new[speciesout]="DO"
    #    'number_SAR_O', "number_SAR_DO_only", "new_reactions","s(DO)=s(O)", "s(DO)=s(DO)","DO_is_largest","Os","unique_DOs",
        
#    new_reactions=0
 #   equal_DO_DO=0
  #  equal_DO=0
   # smaller_DO=0
    #DO_large=0
    #check_if_Dist_Org_dict={}
    #for i in range(len(x)):
    #    add_new_counter=0
    #    new_reac_bool=0
    #    if bool(dict1[i])==True:
    #        new_reac_bool=1
    #        new_reactions+=1
    #        
    ##    #check if DO is also ORG
     #   #bool_result=check_ORG_LOR_expanded(x[i][0],y)
     #   bool_result,speciesout=check_ORG_LOR_expanded(x[i][0],y)
     #   speciesout=frozenset(speciesout)
     #   if bool_result==True:
     #       niceOs.append(x[i][0])
     #       state_dict[i]=1+add_new_counter
     #       check_if_Org_dict[speciesout]=i
     #   else:
     #       #safe DOs which are not Os
     #       niceDOs.append(x[i][0])
     ##       check_if_Dist_Org_dict[speciesout]=i
      ##      if speciesout in check_if_Org_dict:
      #          equal_DO+=1
      #          
      #      else:
      ###          for species in check_if_Org_dict.keys():
      #              if species in speciesout: 
     #                   smaller_DO+=1
     #   
     #       if speciesout in check_if_Dist_Org_dict:
     #           equal_DO_DO+=1
      #      for solution in check_if_Dist_Org_dict.keys():
       #         if speciesout in solution:
        #            DO_large+=1    
    
 #   if state_dict[0]==3:
  #      outputdict["largest_SAR_is"]= "DO"
   # if state_dict[0]==6:
    #    outputdict["largest_SAR_is"]= "DO new rea"
 #   bla={1:0, 2:0 ,3:0 ,4:0 ,5:0 ,6:0}
  #  for i in range(len(x)):
   #     bla[state_dict[i]]+=1
    outputdict["number_SAR_O"]=len(niceOs)
    outputdict["number_SAR_DO_only"]=len(niceDOs)
    outputdict["new_reactions"]=new_reactions
    outputdict["new_reactions_O"]=new_reactions_O
    outputdict["s(DO)=s(O)"]=str(equal_DO[0])+"/"+str(equal_DO[1])
    outputdict["s(DO)=s(DO)"]=str(equal_DO_DO[0])+"/"+str(equal_DO_DO[1])
    #outputdict["s(DO)<s(O)"]=str(smaller_DO[0])+"/"+str(smaller_DO[1])
    outputdict["DO_is_largest"]=str(DO_large)
    outputdict["unique_DOs"]=niceDOs
    outputdict["Os"]=niceOs
    
    
    #outputdict["type1(O)"]=bla[1]
    #outputdict["type2"]=bla[2]
    #outputdict["type3(DO)"]=bla[3]
    #outputdict["type1_new(O)"]=bla[4]
    #outputdict["type2_new"]=bla[5]
    #outputdict["type3_new(DO)"]=bla[6]
    #outputdict["unique_DOs"]=niceDOs
    #outputdict["Os"]=niceOs
    
    return(outputdict)
#simple function to check if first list is subset of secondlist, when bot are sorted
def sublist(list1, list2):
    x=0
    #for each element, we check if it matches an element of the other list
    for ele in list1:
        
        while True:
            try:
                if ele==list2[x]:
    #increase index in list2 if either match or no match
                    x+=1
    #if match break loop and therefore check for next element
                    break
                else:
                    x+=1
    #indexerror shows if list2 ends and therefor list1 isnt subset
            except IndexError:
                return(False)
    #if correct loop through all elements, return true
    return(True)

#function to generate highlighted tet by giving whole knot name and new/to be highlighted terms
def get_highlighted_text(name1, name2):
    name2 = set(name2)-set(name1)
    #goes over all cases of name1 or name 2 being empty
    #uses html as language in digraph file
    
    if name1: 
        name1=list(name1)
        name1.sort()
        name1="<font color='green3'>"+'%s'%', '.join(map(str, name1))+"</font>"
    else: 
        name1=""
    if name1 and name2:
        name1=name1+","
    if name2:
        name2=list(name2)
        name2.sort()
        name2='%s'%', '.join(map(str, name2))
    else:
        name2=""
    name= name1+str(name2)
    return(name)

def draw_hasse(label_highlight_reaction_dict,label_highlight_overprod_dict,dot,list_of_solutions,LOR, DOs=False, shortform=False, show_species=False, show_new=False,show_compartments="length", second_value=False, pdf=True, display_bool=True):
    comp=[]
    namedict={}
    if shortform:
        short_form_knot={}    
        for i in range(len(LOR)):
            short_form_knot[LOR[i].defined_name]="R"+str(i+1)
    for i in range(len(list_of_solutions)):
        
#         name_merge=[]
        addinfo=""
        name="" 
        
        prod_name=""
        start=""
        #change names if shortform is active
        if shortform:  
            name1= [short_form_knot[element] for element in label_highlight_reaction_dict[i]]
            name2= [short_form_knot[element] for element in list_of_solutions[i][0]]
        else:
            name1= label_highlight_reaction_dict[i]
            name2= list_of_solutions[i][0]
        name1_overprod= label_highlight_overprod_dict[i]
        name2_overprod= list_of_solutions[i][1]
        
        if show_new:
            name= ','.join(map(str,name1))
            prod_name= ','.join(map(str,name1_overprod))
        else:
            start="<"
            name= get_highlighted_text(name1,name2)
            if name1==[] and name2==[]:
                name="none"
            name=start+name    
            prod_name= get_highlighted_text(name1_overprod,name2_overprod)
        if prod_name and second_value==True:
            if DOs==False:
                name=name+" OP: "+prod_name
            else:
                name=name+" R: "+prod_name
        if start: name=name+"<br/>" 
        else: name=name+"\n"
        if show_species:
            LOR_spec_dict={}
            for reac in LOR:
                LOR_spec_dict[reac.defined_name]=set(reac.listOfReactants+ reac.listOfProducts)
            speciesset=set()
            for reaction in list_of_solutions[i][0]:
                if "_reverse" in reaction:
                    speciesset.update(set(LOR_spec_dict[reaction[:-8]]))
                else:
                    speciesset.update(set(LOR_spec_dict[reaction]))
            speciesset_ordered=list(speciesset)
            speciesset_ordered.sort()
            addinfo=addinfo+"species: "+'%s'%', '.join(map(str, speciesset_ordered))
            if start: addinfo=addinfo+"<br/>" 
            else: addinfo=addinfo+"\n"

####################################################################
#compartments
        if DOs==True and show_compartments and list_of_solutions[i][0]:
            comp=get_min_compartments(*get_compartments(LOR,list_of_solutions[i][1], extraspecies=list_of_solutions[i][0]))
        if show_compartments and list_of_solutions[i][0] and DOs==False:
            comp=get_min_compartments(*get_compartments(LOR,list_of_solutions[i][0]))
        if not list_of_solutions[i][0]:
            comp={}
        if show_compartments==True:
            comp=[list(comp_ele) for comp_ele in comp ]
            for comp_ele in comp:
                comp_ele.sort()
            comp1= ['%s'%', '.join(map(str, comp_ele))for comp_ele in comp]
            addinfo=addinfo+"{"+'%s'%' | '.join(map(str, comp1))+"}"
            #if start: addinfo=addinfo+"<br>" else: addinfo=addinfo+"\n"
        elif show_compartments=="length":
            addinfo=addinfo+"#"+ str(len(comp))
      
        name=name+addinfo
        if re.search("\<",name):
            name= name+ ">"
        namedict['node'+str(i)]=name
    
    for i in range(len(list_of_solutions)-1,-1,-1):
        print(list_of_solutions[i][0])

        is_o=True
        if DOs==True:
            true_rea_DO=[]
            knot_reac=list_of_solutions[i][1]
            for r in LOR:
                        if set(r.listOfReactants).issubset(set(list_of_solutions[i][0])):
                            
                            true_rea_DO.append(r)
            print(true_rea_DO)
            print(list_of_solutions[i][1])
            for reac in true_rea_DO:
                    if reac.defined_name not in knot_reac:
                        is_o=False
        else:
            knot_reac=list_of_solutions[i][0]
        
        if check_ORG_LOR(knot_reac,LOR) and is_o:
            dot.node('node'+str(i), namedict['node'+str(i)], shape='box')
            
        else:
            dot.node('node'+str(i), namedict['node'+str(i)])
   
    line=dot.source
    line=[e+"]" for e in line.split("]")]
    line[-1]= re.sub("]","", line[-1] )
    for i in range(len(line)):
        if re.search("\<",line[i]):
            
            line[i]= re.sub("\"","", line[i] )
          
    if shortform:
        shortform_label_string="reaction: shortform \n"
        for ele in short_form_knot:
            shortform_label_string=shortform_label_string+ele+": "+short_form_knot[ele]+"\n"
        #dot.node('node'+str(number_of_nodes+3), shortform_label_string, pos = '1000,1000!',shape ="box")       
    addstring_merged = '\n'.join(map(str, line))
    
    dot = Source(addstring_merged)
#     print(dot.source)
    if pdf==True:
        dot.render('test-output/Hasse_Diagram.gv', view=True)  
        #only works in jupyter notebook
#    if display_bool==True:
#        display(dot)

    
    
def create_Hasse3(list_of_solutions,LOR, DOs=False, shortform=False, show_species=False, show_new=False,show_compartments="length", second_value=False, analyze_only=False, pdf=True):
    #inverse solution
    #print(list_of_solutions)
    max_knot=set(list_of_solutions[0][0])
    max_second=set(list_of_solutions[0][1])
    inverse_list=[(max_knot.difference(set(list_of_solutions[i][0])),max_second.difference(set(list_of_solutions[i][1]))) for i in range(len(list_of_solutions))]
    knot_to_index={}
    
    #extralists for special marking of new reactions/production
    label_highlight_overprod_dict={}
    label_highlight_reaction_dict={}
    
    #create elements for edge algorithm
    facedict={}
    borderlist=[]
    candidatelist=[]
    Os_list=[]
    #create objects for hasse
    dot = Digraph(strict=True,comment='Hasse_Diagram')
    nodedict={}
    Os_dict={}
    check_if_Org_dict2={}
    check_if_Org_dict={}   
    for i in range(len(list_of_solutions)):
        #print("letsgo")
        #print(list_of_solutions[i][0])
        #setup of objects, which are changed before looping over its knot
        knot_to_index[frozenset(inverse_list[i][0])]=i
        label_highlight_reaction_dict[i]=list_of_solutions[i][0][:]
        label_highlight_overprod_dict[i]=list_of_solutions[i][1][:]
        setconverter=frozenset(inverse_list[i][0])
        """lw"""
#         nodedict[frozenset(inverse_list[i][0])]=str('node'+str(i))
        nodedict[frozenset(inverse_list[i][0])]=i
        if analyze_only:
            if check_ORG_LOR(list_of_solutions[i][0], LOR):
                Os_dict[frozenset(list_of_solutions[i][0])]=True
            else: Os_dict[frozenset(list_of_solutions[i][0])]=False
        if DOs==False:
            bool_result,speciesout=check_ORG_LOR_expanded(list_of_solutions[i][0],LOR)
            check_if_Org_dict2[i]=speciesout
#         if bool_result:
#             check_if_Org_dict[speciesout]=i
            
#             print(speciesout)
    #initialization of first knot
    borderlist.append(inverse_list[0][0])
    facedict[frozenset(inverse_list[0][0])]=[]

    for i in range(1,len(list_of_solutions)):
        
        #create cacndidates by intersection of border and current knot
        candidatelist=[]
        for borderelement in borderlist:
            candidatelist.append(set(inverse_list[i][0]) & (set(borderelement)))
            
    
            
        facedict[frozenset(inverse_list[i][0])]=[]    
        
        #check if candidates have larger elements in face (elements smaller than knot and larger than candidate)
        for candidate in candidatelist:
            candidate_froz= frozenset(candidate)

            for faceobject in facedict[candidate_froz]:
                if faceobject.issubset(inverse_list[i][0]):
                    break
            else: 
                #create edge
#                 dot.edge(nodedict[frozenset(tuple(inverse_list[i][0]),)], nodedict[candidate_froz], arrowhead = "none")
                if DOs==False:
                    if check_if_Org_dict2[nodedict[candidate_froz]]==check_if_Org_dict2[i]:
                        dot.edge( str('node'+str(nodedict[candidate_froz])),str('node'+str(nodedict[frozenset(tuple(inverse_list[i][0]),)])), arrowhead = "none",color="red")
                    else:
                        dot.edge( str('node'+str(nodedict[candidate_froz])),str('node'+str(nodedict[frozenset(tuple(inverse_list[i][0]),)])), arrowhead = "none")
                else:
                    dot.edge( str('node'+str(nodedict[candidate_froz])),str('node'+str(nodedict[frozenset(tuple(inverse_list[i][0]),)])), arrowhead = "none")
                
                j=knot_to_index[candidate_froz]
                
                #save new reactions
                if analyze_only:
                    if Os_dict[frozenset(list_of_solutions[i][0])]:
                        liste=list_of_solutions[i][0]
                    else:
                        liste=[e for e in list_of_solutions[i][0] if e not in label_highlight_reaction_dict[i]]
                    label_highlight_reaction_dict[j]= [e for e in label_highlight_reaction_dict[j] if e not in list_of_solutions[i][0]]
                else:
                    
                    label_highlight_reaction_dict[j]= [e for e in label_highlight_reaction_dict[j] if e not in list_of_solutions[i][0]]
                    label_highlight_overprod_dict[j]= [e for e in label_highlight_overprod_dict[j] if e not in list_of_solutions[i][1]]
                #add knot to face of candidate
                add_element_to_facedict=inverse_list[i][0].copy()
                extract=facedict[candidate_froz]
                extract.append(frozenset(tuple(add_element_to_facedict)))
                facedict[candidate_froz]=extract
                #remove candidate from border
                #these candidates dont have to be in the border if there they are subsets of bigger knots in the border
                try:
                    borderlist.remove(candidate)
                except ValueError:
                    pass


        borderlist.append(inverse_list[i][0])
    if analyze_only:
        return(label_highlight_reaction_dict,label_highlight_overprod_dict)
    else:
        draw_hasse(label_highlight_reaction_dict,label_highlight_overprod_dict,dot,list_of_solutions,LOR,DOs, shortform, show_species, show_new, show_compartments, second_value, pdf)
        
def create_Hasse(list_of_solutions,LOR, DOs=False, shortform=False, show_species=False, show_new=False,show_compartments="length", second_value=False,analyze_only=False, pdf=True):
    dot = Digraph(strict=True,comment='Hasse_Diagram')
    #dot.attr( size='8,5')
    number_of_nodes=0
    nodedict={}
    comp=[]
    knot_reac=[]
    cover={}
    list_of_solutions.reverse()
    edge_memory=[]
    label_highlight_overprod_dict={}
    label_highlight_reaction_dict={}
    index_reactions_dict={}
    namedict={}
    ################################################################################################################
    #naming
    list_of_solutions_empty=copy.deepcopy(list_of_solutions)
    if list_of_solutions[0][0]==[]:
        [s[0].append("empty_reaction") for s in list_of_solutions_empty]
    else:
        list_of_solutions_empty=list_of_solutions
    if shortform:
        short_form_knot={}    
        list_of_reactions=list_of_solutions[-1][0]
        for i in range(len(list_of_reactions)):
            short_form_knot[list_of_reactions[i]]="R"+str(i)
    ##########################################################################################################
    #setup first elements
    for i in range(len(list_of_solutions)):
        cover[i]=[]
        #umwandlung in frozenset um für dictionary vorzubereiten
        index_reactions_dict[i]=list_of_solutions_empty[i][0]
    cover[0]=[0]
    label_highlight_reaction_dict[0]=list_of_solutions[0][0]
    label_highlight_overprod_dict[0]=list_of_solutions[0][1]
##################################################################################################
#loop over all knots
    for i in range (1, len(list_of_solutions)):
        reactions_to_check=list_of_solutions[i][0][:]
        production_to_check=list_of_solutions[i][1][:]
        checklist=[index for index in range(i)]
        checkdic={}
        for j in range(i):
            checkdic[j]=""
        startcheck=i
        while True:
            for k in range(startcheck,-2,-1):
                if k in checkdic:
                    check_index=k
                    break
            
            startcheck=check_index
            if sublist(list_of_solutions_empty[check_index][0],list_of_solutions_empty[i][0]):
                dot.edge("node"+str(i),"node"+str(check_index),  arrowhead = "none")
                reactions_to_check = [e for e in reactions_to_check if e not in list_of_solutions[check_index][0]]
                production_to_check = [e for e in production_to_check if e not in list_of_solutions[check_index][1]]
                
                #add all elements of linked knots to cover of current knot
                for ele in cover[check_index]:
                    if ele not in cover[i]:
                        cover[i]=cover[i]+[ele]
    
                for del_element in sorted(cover[i], reverse=True):
                    checkdic.pop(del_element, None)
            else:
                checkdic.pop(check_index, None)
            if len(checkdic)==0:
                cover[i]=cover[i]+[i]
                break
        label_highlight_reaction_dict[i]=reactions_to_check
        label_highlight_overprod_dict[i]=production_to_check
    if analyze_only:
        return(label_highlight_reaction_dict,label_highlight_overprod_dict)
    else:
   # draw_hasse(label_highlight_reaction_dict,label_highlight_overprod_dict,dot,list_of_solutions,LOR,DOs, shortform, show_species, show_new, show_compartments, overproduction)
        draw_hasse(label_highlight_reaction_dict,label_highlight_overprod_dict,dot,list_of_solutions,LOR,DOs, shortform, show_species, show_new, show_compartments, second_value, pdf)
