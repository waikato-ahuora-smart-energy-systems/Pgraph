import subprocess
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
from networkx.drawing.nx_pydot import pydot_layout
from lxml import etree
import networkx as nx
import platform
import tempfile


class Pgraph():
    def __init__(self, problem_network, mutual_exclusion=[], solver="INSIDEOUT",max_sol=100, input_file=None):
        ''' 
        Pgraph(problem_network, mutual_exclusion=[], solver="INSIDEOUT",max_sol=100)
                
        Description
        This function initializes the Pgraph object and prepares the solver.
        
        Arguments
        problem_network: (DiGraph() object) Directed graph object specified using networkx
        mutual_exclusion: (list of list) List of lists containing mutually excluded elements. Symbols of nodes should be used. e.g. "M1"
        solver: (str): Solver type that is used. Possibilities include "MSG", "SSG", "SSGLP" (for SSG+LP), "INSIDEOUT" (for ABB)
        max_sol: (int): Maximum number of solutions required for the solver.
        
        '''
    
        #In case names is not specified, revert to symbol of graph
        for n in problem_network:
            if 'names' not in list(problem_network.nodes()[n].keys()):
                problem_network.nodes()[n]['names']=n
            
        self.G=problem_network
        self.ME=mutual_exclusion
        self.solver=solver
        self.max_sol=max_sol
        self.library_path=os.path.dirname(os.path.realpath(__file__))+r"/solver/"
        self.gmatlist=[]
        self.goplist=[]
        self.goolist=[]
        self.wine_installed=False #For Linux Only
        self.input_file=input_file
        
    def plot_problem(self,figsize=(5,10),padding=0.25,titlepos=0.95,rescale=2,box=True,node_size=3000):
        '''
        plot_problem(self,figsize=(5,10),padding=0,titlepos=0.95,rescale=2,box=True)
                
        Description
        This function plots the network figure of the problem. 
        
        Arguments
        figsize: (tuple) Size of figure
        padding: (float) Padding of the sides of the image
        titlepos: (float) Position of title in the figure
        rescale: (float) Rescaling the axis of the figures. Makes nodes further apart and appear smaller.
        box: (boolean) Whether the figure should appear square.
        
        Return:
        ax: (matplotlib.axes) Axes of the figure. Can be manipulated further before plotting.
        '''
        
        G=self.G.copy()
        for n in G.nodes():
            if n[0]=="O":
                G.nodes[n]['s']=mpl.markers.MarkerStyle(marker='s', fillstyle='top')
            else:
                G.nodes[n]['s']='o'
        plt.rc('figure',figsize=figsize)
        label_options = {"ec": "k", "fc": "white", "alpha": 0.8}
        edges=G.edges()
        weights = [G[u][v]['weight'] for u,v in edges]
        labels = nx.get_edge_attributes(G,'weight')
        labels={k:round(v,2) for k,v in labels.items()}
        node_labels={n:G.nodes()[n]['names']  for n in G.nodes()}

        nodeShapes = set((aShape[1]["s"] for aShape in G.nodes(data = True)))
        pos=pydot_layout(G,prog='dot')
        pos2=pydot_layout(G,prog='dot')
        for key, (v1,v2) in pos2.items():
            if key[0]=="O":
                pos2[key]=(v1,v2-3)
            
        pos=nx.rescale_layout_dict(pos,scale=rescale)
        pos2=nx.rescale_layout_dict(pos2,scale=rescale)
        nx.draw_networkx(G, pos=pos, node_color='white',alpha=0.9,node_shape="o", edge_color='black',labels=node_labels, with_labels = True,node_size=node_size,bbox=label_options,width=weights,font_size=12)
        for aShape in nodeShapes:
            if aShape=="o":
                node_size1=node_size
            else:
                node_size1=node_size*2
            nx.draw_networkx_nodes(G,pos=pos2,node_color='black',node_shape = aShape, nodelist = [sNode[0] for sNode in filter(lambda x: x[1]["s"]==aShape,G.nodes(data = True))],node_size=node_size1)
        
        
        for key, values in nx.get_node_attributes(G,'type').items():
            if values=="raw_material":
                nx.draw_networkx_nodes(G,pos=pos,node_color='white',node_shape = 'v', nodelist = [key],node_size=node_size/3*1.6)
                
            elif values=="product":
                nx.draw_networkx_nodes(G,pos=pos,node_color='white',node_shape = 'o', nodelist = [key],node_size=node_size/3*2)
                nx.draw_networkx_nodes(G,pos=pos,node_color='black',node_shape = 'o', nodelist = [key],node_size=node_size/3*1.25)
                nx.draw_networkx_nodes(G,pos=pos,node_color='white',node_shape = 'o', nodelist = [key],node_size=node_size/3*0.75)
                
                
        nx.draw_networkx_edge_labels(G, pos=pos,edge_labels=labels)
        ax= plt.gca()
        plt.axis('off')
        ax.autoscale()
        ax.set_xlim([ax.get_xlim()[0]+padding*ax.get_xlim()[0],ax.get_xlim()[1]+padding*ax.get_xlim()[1]])
        ax.set_ylim([ax.get_ylim()[0]+padding*ax.get_ylim()[0],ax.get_ylim()[1]+padding*ax.get_ylim()[1]])
        if box:
            ax.set_aspect('equal', adjustable='box')
        ax.set_title("Original Problem ",y=titlepos)
        return ax
    
    def create_solver_input(self, input_file=None):
        '''
        create_solver_input()
        
        Description
        This function creates the solver input from the networkx DiGraph() object specified.
        
        Arguments
        input_file: (string) Path of the input file to be generated. 

        '''

        G=self.G
        ME=self.ME
        path=self.library_path
        ### MAKE INPUT FILE #############
        
        # Declare default values
        material_flow_rate_lower_bound = 0
        material_flow_rate_upper_bound =1000000000
        material_price = 0
        operating_unit_capacity_lower_bound = 0
        operating_unit_capacity_upper_bound = 1000000000
        operating_unit_fix_cost = 0
        operating_unit_proportional_cost = 0
        material_type = 'raw_material'
        mass_unit = 't'
        time_unit = 'y'
        money_unit = 'NZD'
        
        prelines=[
        'file_type=PNS_problem_v1','\n',
        'file_name=Graph_1','\n',
        '\n',
        'measurement_units:','\n',
        'mass_unit=',mass_unit,'\n',
        'time_unit=',time_unit,'\n',
        'money_unit=',money_unit,'\n',
        '\n',
        'defaults:','\n',
        'material_type=',str(material_type),'\n',
        'material_flow_rate_lower_bound=',str(material_flow_rate_lower_bound),'\n',
        'material_flow_rate_upper_bound=',str(material_flow_rate_upper_bound),'\n',
        'material_price=',str(material_price),'\n',
        'operating_unit_capacity_lower_bound=',str(operating_unit_capacity_lower_bound),'\n',
        'operating_unit_capacity_upper_bound=',str(operating_unit_capacity_upper_bound),'\n',
        'operating_unit_fix_cost=',str(operating_unit_fix_cost),'\n',
        'operating_unit_proportional_cost=',str(operating_unit_proportional_cost),'\n',        '\n',
        'materials:','\n'
        ]

        for n in G.nodes():
            if n[0]=="M":
                add_list=[n+": "]
                for k,v in G.nodes[n].items():
                    if k=='type':
                        add_list[0]=add_list[0]+v
                    elif k not in ["names", "units"]:
                        add_list.append(str(k)+"="+str(v))
                n_=", ".join(add_list)
                prelines.append(n_)
                prelines.append('\n')

        prelines.append('\n')
        prelines.append('operating_units:')
        prelines.append('\n')
        
        self.problem_materials={'ID':[],'names':[],'type':[],'flow_rate_lower_bound':[],'flow_rate_upper_bound':[],'price':[],'time_unit':[],'money_unit':[],'mass_unit':[]}
        for n in G.nodes():
            if n[0]=="M":
                #print(n)
                self.problem_materials['ID'].append(n)
                for k,v in G.nodes[n].items():
                    #print('items',G.nodes[n])
                    if k not in ['units']:
                        self.problem_materials[k].append(v)
                        # print(k,':-',v)
                    elif k in ['units']:
                        #print(v)
                        for i,j in enumerate(v):
                            self.problem_materials[j].append(v[j])
                            # print(i,'yy',j,'jj',v[j])
                    
                    # else:

                        #self.problem_materials['units'][]
                if 'flow_rate_lower_bound' not in G.nodes[n]:
                    self.problem_materials['flow_rate_lower_bound'].append(material_flow_rate_lower_bound)
                if 'flow_rate_upper_bound' not in G.nodes[n]:
                    self.problem_materials['flow_rate_upper_bound'].append(material_flow_rate_upper_bound)
                if 'price' not in G.nodes[n]:
                    self.problem_materials['price'].append(material_price)
                if 'type' not in G.nodes[n]:
                    self.problem_materials['type'].append(material_type)
                if 'units' not in G.nodes[n] or not all(key in G.nodes[n]['units'] for key in ['time_unit', 'money_unit', 'mass_unit']):
                    self.problem_materials['money_unit'].append(money_unit)
                    self.problem_materials['time_unit'].append(time_unit)
                    self.problem_materials['mass_unit'].append(mass_unit)                    
                    


        self.namelist_mat={"ID":[],"names":[],'type':[]}
        for n in G.nodes():

            if n[0] =="M":
                for k, v in G.nodes[n].items():
                    #print(G.nodes[n].items())
                    #print(v)
                    if k =='names':
                        self.namelist_mat['ID'].append(n)
                        self.namelist_mat['names'].append(v)
                    elif k=='type':
                        self.namelist_mat['type'].append(v)


        #print(self.namelist_mat)

        self.namelist_op = {"ID": [], "names": []}
        for n in G.nodes():
            if n[0] == "O":
                for k, v in G.nodes[n].items():
                    # print(G.nodes[n].items())
                    if k == 'names':
                        self.namelist_op['ID'].append(n)
                        self.namelist_op['names'].append(v)

        #print(self.namelist_op)
        self.problem_operating={'ID':[],'names':[],'fix_cost':[],'proportional_cost':[],'capacity_lower_bound':[],'capacity_upper_bound':[]}
        for n in G.nodes():
            # print('try',n)
            if n[0]=="O":
                self.problem_operating['ID'].append(n)
                for k,v in G.nodes[n].items():
                    # print(k,'gg',v)
                    # if k in ['capacity_lower_bound'] or k in ['capacity_upper_bound']:
                    #     pass

                    
                    if k not in ['units']:
                        self.problem_operating[k].append(v)
                        #print(self.problem_operating[k])
                if 'capacity_lower_bound' not in G.nodes[n]:
                    self.problem_operating['capacity_lower_bound'].append(operating_unit_capacity_lower_bound)
                if 'capacity_upper_bound' not in G.nodes[n]:
                    self.problem_operating['capacity_upper_bound'].append(operating_unit_capacity_upper_bound)
                if 'fix_cost' not in G.nodes[n]:
                    self.problem_operating['fix_cost'].append(operating_unit_fix_cost)
                if 'proportional_cost' not in G.nodes[n]:
                    self.problem_operating['proportional_cost'].append(operating_unit_proportional_cost)

        for n in G.nodes():
            if n[0]=="O":        
                add_list=[n+": "]
                first=True
                for k,v in G.nodes[n].items():
                    if k=='type':
                        add_list[0]=add_list[0]+v
                    elif len(add_list)==1 and first and k!="names":
                        add_list[0]=add_list[0]+str(k)+"="+str(v)
                        first=False
                    elif k not in ["names", "units"]:
                        add_list.append(str(k)+"="+str(v))
                n_=", ".join(add_list)
                prelines.append(n_)
                prelines.append('\n')

        prelines.append('\n')
        prelines.append('material_to_operating_unit_flow_rates:')
        prelines.append('\n')

        for n in G.nodes():
            if n[0]=="O": 
                add_str=n+": "
                for x1, x2 in G.in_edges(n):
                    add_str=add_str+str(G[x1][x2]["weight"])+" "+x1
                    add_str=add_str+" + "
                if G.in_edges(n):
                    add_str=add_str[:-3]
                add_str=add_str+" => "
                for x1, x2 in G.out_edges(n):
                    add_str=add_str+str(G[x1][x2]["weight"])+" "+x2
                    add_str=add_str+" + "
                add_str=add_str[:-3]
                prelines.append(add_str)
                prelines.append("\n")
        if ME!=[]:
            prelines.append("\n")

        if ME !=[]:
            prelines.append("mutually_exlcusive_sets_of_operating_units:")
            prelines.append("\n")
        for i in range(len(ME)):
            MEstring="ME"+str(i)+": "
            for x in ME[i]:
                MEstring=MEstring+x+", "
            MEstring=MEstring[:-2]  
            prelines.append(MEstring)
            prelines.append("\n")
                
        prelines.append("\n")
        

        system = platform.system()  # Detect the current operating system
        
        if system=="Windows": #support for windows
            with open(path+'input.in', 'w') as f:
                for line in prelines:
                    f.write(line)    
        elif system=="Linux":      
            with open(input_file, 'w') as f:
                for line in prelines:
                    f.write(line)

    def solve(self,system=None,skip_wine=False, solver_name='pgraph_solver.exe', input_file: str = "", output_file: str = "", path=None):
        '''
        solve(system=None,skip_wine=False)
        
        Description
        Runs the solver.
        
        Arguments:
        system: (string) (optional) Operating system. Options of "Windows", "Linux". MacOS is not supported yet. Specifying this makes function slightly faster.
        skip_wine: (boolean) Only relevent for Linux. Skip the dependency "wine" if it is already installed. 
        solver_name: (string) For advanced users only. Choose your customized solver. 'pgraph_solver.exe' or 'pgraph_solver_new.exe'
        path: (string) path to the custom solver. If None, then the default library installation path will be used.
        '''
        path=self.library_path
        max_sol=self.max_sol
        solver=self.solver
        # solver_dict={0:"MSG",1:"SSG",2:"SSGLP",3:"INSIDEOUT"}
        
        if system==None:
            system=platform.system()
            
        if system=="Windows": #support for windows
            if type(self.input_file)==str:
                rc=subprocess.run([path+solver_name,solver, self.input_file, path+"test_out.out", str(max_sol)])
            else:
                rc=subprocess.run([path+solver_name,solver, path+"input.in", path+"test_out.out", str(max_sol)])                
        elif system=="Linux":
            # detect architecture, ARM or x86
            if platform.machine() in ["arm64", "aarch64"]:
                # use the ARM-compiled pns solver
                solver_name = "pns_solver_2018_ARM"
            elif platform.machine() in ["x86_64"]:
                # use the x86-compiled pns solver
                solver_name = "pns_solver_2018_x86"
            else:
                # unknown architecture, raise an error
                raise Exception("Could not run solver, unknown architecture: " + platform.machine())
            
            subprocess.run([
                path + solver_name,
                solver,
                input_file,
                output_file,
                str(max_sol)
            ])

        ################
    
    def read_solutions(self, output_file=None):
        '''
        read_solutions()
        
        Description
        Reads the solution from the solver.     
        
        Arguments
        output_file: (string) Path to generate the input file in. If None, then the default library installation path will be used.       
        '''
    
        path=self.library_path
        gmatlist=[]
        goplist=[]
        goolist=[]
        
        system = platform.system()  # Detect the current operating system
        
        #clean strings
        if system=="Windows": #support for windows
            with open(path+"test_out.out","r") as f:
                lines = f.readlines()
        elif system=="Linux":    
            with open(output_file, 'r') as f:    
                lines = f.readlines()
        for i in range(len(lines)-1,1,-1):
            if lines[i-1].strip() == "Operating units(1):":
                continue
            if lines[i][0]==" " or (len(lines[i].strip())>0 and ":" not in lines[i] and "," not in lines[i] and "= " not in lines[i] and "End." not in lines[i]): ## Attention for possible future changes
                if lines[i-1][-3:]!="\n":
                    lines[i-1]=lines[i-1].rstrip()+" "+lines[i].strip()
                    lines[i]=""
                else:
                    lines[i-1]=lines[i-1][:-3].rstrip()+" "+lines[i].strip()
                    lines[i]=""
        lines=list(map(lambda x:x.strip(),lines))
        lines=list(filter(None, lines))
        #Lines are now clean with next lines combined as elements of list.
        
        ###### Read for the case of SSGLP and INSIDEOUT (ABB) ######
        if self.solver in ["SSGLP","INSIDEOUT",2,3]:
            mat_list=lines[1]
            op_list=lines[3] 
            used_mat_list=lines[6] 

            #Find solutions via Feasible Structure tag
            sol_start_index=[]
            for i in range(len(lines)):
                if lines[i][:18]=="Feasible structure":
                    sol_start_index.append(i)

            sol_start_index.append(len(lines)-1)

            sol_list=[]
            for i in range(1,len(sol_start_index)):
                sol_list.append(lines[sol_start_index[i-1]:sol_start_index[i]])

            comp=["Materials:","Operating units:","Total annual cost="]
            for i in range(len(sol_list)): #loop through solution number
                comp_ind=-1
                s=False
                tmatlist=[]
                toplist=[]
                for j in range(len(sol_list[i])):
                    if sol_list[i][j][:len(comp[0])]==comp[0]: #Materials
                        comp_ind=0
                        s=True
                    elif sol_list[i][j][:len(comp[1])]==comp[1]:    #Operating units
                        comp_ind=1
                        s=True
                    elif sol_list[i][j][:len(comp[2])]==comp[2]:   # Total annual cost
                        comp_ind=2
                        s=True
                    if s==False:
                        if comp_ind==0: #Materials
                           
                            tlist=sol_list[i][j].replace('(',' ')
                            tlist=tlist.replace(')','')
                            tlist=tlist.split()
                            if tlist[0][-1]==":":
                                tlist[0]=tlist[0][:-1] #correct for semicolon  
                            if tlist[1]=="balanced": #correct for balanced
                                tlist=[tlist[0],0,0,0]
                            tmatlist.append(tlist)
                        elif comp_ind==1: #Operating units
                            glist=sol_list[i][j].split(')')
                            glist=glist[0].replace('*',' ')
                            glist=glist.replace('(', ' ')
                            glist=glist.split()
                            toplist.append(glist)
                    if comp_ind==2:  #Total annual cost
                        goolist.append(sol_list[i][j].split()[3])
                    s=False

                goplist.append(toplist)
                gmatlist.append(tmatlist)
                
            for solution in gmatlist:
                for material in solution:
                    # print(material)
                    material_name = material[0]
                    node = self.G.nodes[material_name]
                    #print(len(material))
                    if len(material) == 5:
                        if not node.get('units'):
                            material[2] = "NZD/t"
                        else:
                            units = node.get("units")
                            time_unit = units.get("time_unit")
                            money_unit = units.get("money_unit")
                            mass_unit = units.get("mass_unit")
                            material[2] = money_unit + "/" + time_unit
                            material[4] = mass_unit + "/" + time_unit
                            # print(material)

            for solution in goplist:
                for operator in solution:
                    operator_name = operator[1]
                    node = self.G.nodes[operator_name]
                    #if node['units']:
                    if not node.get('units'):
                        operator[3] = "NZD/t"
                    else:
                        units = node.get("units")
                        time_unit = units.get("time_unit")
                        money_unit = units.get("money_unit")
                        operator[3] = money_unit + "/" + time_unit
            self.goplist=goplist
            self.gmatlist=gmatlist
            self.goolist=goolist
            if len(goolist)==0: 
                print("No Feasible Solution Found!")
        
        ###### Read for the case MSG ######
        if self.solver in ["MSG",0,"SSG",1]:
            for i in range(len(lines)):
                #### maximal structure ####
                if self.solver in ["MSG", 0] and lines[i]=="Maximal Structure:":
                    if lines[i+1][9:]!="(0):":
                        gmatlist.append(lines[i+2].split(", "))  # materials
                        goplist.append(lines[i+4].split(", "))   # operating units
                    else:
                        gmatlist.append([])
                        goplist.append([])
                    goolist.append("0")
                
                #### solution structure ####
                if lines[i][:19]=="Solution structure ":
                    goolist.append(lines[i].split("#")[1][:-1]) #SSG number
                    if lines[i+1][9:]!="(0):":
                        gmatlist.append(lines[i+2].split(", "))  # materials
                        goplist.append(lines[i+4].split(", "))   # operating units
                    else:
                        gmatlist.append([])
                        goplist.append([])
            self.goplist=goplist
            self.gmatlist=gmatlist
            self.goolist=goolist
    
    def get_solution_as_network(self, sol_num=0, return_excluded_nodes=False):
        '''
        get_solution_as_network(sol_num=0)
                
        Description
        This function returns the DiGraph() object of the solution
        
        Arguments
        sol_num: (int) Index of the solution to be returned
        return_excluded_nodes: (boolean) If True, the returnedd network will also contain the excluded nodes (ignored by algorithms only generating structures, e.g., SSG)
        
        Return:
        H: (networkx DiGraph() object) Directed Graph object of the solution.
        '''
        sol_num=sol_num
        H=self.G.copy()
        gmatlist=self.gmatlist
        goplist=self.goplist
        goolist=self.goolist
        if self.solver in ["SSGLP","INSIDEOUT",2,3] :
            if not return_excluded_nodes:
                total_node=[node[0] for node in gmatlist[sol_num]]+[node[1] for node in goplist[sol_num]]
                node_list=list(H.nodes()).copy()
                for n in node_list:
                    if n not in total_node:
                        H.remove_node(n)
            for n in H.nodes():
                if n[0]=="O":
                    H.nodes[n]['s']=mpl.markers.MarkerStyle(marker='s', fillstyle='top')
                else:
                    H.nodes[n]['s']='o'
            attr_op={x[1]:{"Capacity":x[0],"Cost":x[2]} for x in goplist[sol_num]}
            attr_mat={x[0]:{"Flow":x[3],"Cost":x[1]} for x in gmatlist[sol_num]}
            nx.set_node_attributes(H,attr_op)
            nx.set_node_attributes(H,attr_mat)
        else:
            total_node=gmatlist[sol_num]+goplist[sol_num]
            node_list=list(H.nodes()).copy()
            for n in node_list:
                if n not in total_node:
                    H.remove_node(n)
       
        return H
    
    def plot_solution(self,sol_num=0,figsize=(5,10),padding=0.25,titlepos=0.95,rescale=2,box=True,node_size=3000):
        '''
        plot_solution(sol_num=0,figsize=(5,10),padding=0,titlepos=0.95,rescale=2,box=True)
                
        Description
        This function plots the network figure of the solution via its index. 
        
        Arguments
        sol_num: (int) Index of the solution to be plotted
        figsize: (tuple) Size of figure
        padding: (float) Padding of the sides of the image
        titlepos: (float) Position of title in the figure
        rescale: (float) Rescaling the axis of the figures. Makes nodes further apart and appear smaller.
        box: (boolean) Whether the figure should appear square.
        
        Return:
        ax: (matplotlib.axes) Axes of the figure. Can be manipulated further before plotting.
        '''
        
        sol_num=sol_num
        H=self.G.copy()
        gmatlist=self.gmatlist
        goplist=self.goplist
        goolist=self.goolist
        

        if self.solver in ["SSGLP","INSIDEOUT",2,3] and len(goolist)!=0:
            for n in H.nodes():
                if n[0]=="O":
                    H.nodes[n]['s']=mpl.markers.MarkerStyle(marker='s', fillstyle='top')
                else:
                    H.nodes[n]['s']='o'
            attr_op={x[1]:{"Capacity":x[0],"Cost":x[2]} for x in goplist[sol_num]}
            attr_mat={x[0]:{"Flow":x[3],"Cost":x[1]} for x in gmatlist[sol_num]}
            nx.set_node_attributes(H,attr_op)
            nx.set_node_attributes(H,attr_mat)
            plt.rc('figure',figsize=figsize)
            label_options = {"ec": "k", "fc": "white", "alpha": 0.8}
            edges=H.edges()
            weights = [H[u][v]['weight'] for u,v in edges]
            labels = nx.get_edge_attributes(H,'weight')
            labels={k:round(v,2) for k,v in labels.items()}
            labels_flow=nx.get_node_attributes(H,'Flow')
            labels_cap=nx.get_node_attributes(H,'Capacity')
            labels_cost=nx.get_node_attributes(H,'Cost')
            all_node=list(set(list(labels_flow.keys())+list(labels_cap.keys())+list(labels_cost.keys())))
            labels1={}
            for x in all_node:
                string=H.nodes()[x]['names']
                if labels_flow.get(x) is not None:
                    string=string+"\nFlow="+str(abs(float(labels_flow.get(x))))
                
                if labels_cap.get(x) is not None:
                    string=string+"\nCap.="+str(labels_cap.get(x))
                    
                if labels_cost.get(x) is not None:    
                    string=string+"\nCost="+str(labels_cost.get(x))
                
                labels1.update({x:string})

            for e in H.edges():
                if e[0] in all_node and e[1] in all_node:
                    H[e[0]][e[1]]['color']='black'
                else:
                    H[e[0]][e[1]]['color']='lightgrey'
                    
            for n in H.nodes():
                if n in all_node:
                    H.nodes()[n]['color']='black'
                else:
                    H.nodes()[n]['color']='lightgrey'
                if n[0]=='M':
                    H.nodes()[n]['shape']='o'
                elif n[0]=='O':
                    H.nodes()[n]['shape']='s'

            shape_list=[H.nodes()[n]['shape'] for n in H.nodes()]  
            node_color_list=[H.nodes()[n]['color'] for n in H.nodes()]        
            edge_color_list=[H[e[0]][e[1]]['color'] for e in H.edges()]    
            
            nodeShapes = set((aShape[1]["s"] for aShape in H.nodes(data = True)))

            pos=pydot_layout(H,prog='dot')
            pos2=pydot_layout(H,prog='dot')
            

            for key, (v1,v2) in pos2.items():
                if key[0]=="O":
                    pos2[key]=(v1,v2-3)
            pos=nx.rescale_layout_dict(pos,scale=rescale)
            pos2=nx.rescale_layout_dict(pos2,scale=rescale)
            nx.draw_networkx(H, pos=pos,labels=labels1, node_color='white',alpha=0.9,node_shape='o', edge_color=edge_color_list, with_labels = True,node_size=node_size,bbox=label_options,width=weights,font_size=10)
            for aShape in nodeShapes:
                node_list=[sNode[0] for sNode in filter(lambda x: x[1]["s"]==aShape,H.nodes(data = True))]
                if aShape=="o":
                    node_size1=node_size
                else:
                    node_size1=node_size*2
                for node in node_list:
                    nx.draw_networkx_nodes(H,pos=pos2,node_color=H.nodes()[node]['color'],node_shape = aShape, nodelist =[node] ,node_size=node_size1) 
            
            for key, values in nx.get_node_attributes(H,'type').items():
                if values=="raw_material":
                    nx.draw_networkx_nodes(H,pos=pos,node_color='white',node_shape = 'v', nodelist = [key],node_size=node_size/3*1.6)
                elif values=="product":
                    nx.draw_networkx_nodes(H,pos=pos,node_color='white',node_shape = 'o', nodelist = [key],node_size=node_size/3*2)
                    nx.draw_networkx_nodes(H,pos=pos,node_color=H.nodes()[key]['color'],node_shape = 'o', nodelist = [key],node_size=node_size/3*1.25)
                    nx.draw_networkx_nodes(H,pos=pos,node_color='white',node_shape = 'o', nodelist = [key],node_size=node_size/3*0.75)
            
            nx.draw_networkx_edge_labels(H, pos=pos,edge_labels=labels)
            ax= plt.gca()
            plt.axis('off')
            ax.autoscale()

            ax.set_xlim([ax.get_xlim()[0]+padding*ax.get_xlim()[0],ax.get_xlim()[1]+padding*ax.get_xlim()[1]])
            ax.set_ylim([ax.get_ylim()[0]+padding*ax.get_ylim()[0],ax.get_ylim()[1]+padding*ax.get_ylim()[1]])
            if box:
                ax.set_aspect('equal', adjustable='box')
            ax.set_title("Solution #"+str(sol_num+1)+" Total Costs="+str(goolist[sol_num]),y=titlepos)
            
        if self.solver in ["SSG","MSG",0,1]:
            for n in H.nodes():
                if n[0]=="O":
                    H.nodes[n]['s']=mpl.markers.MarkerStyle(marker='s', fillstyle='top')
                else:
                    H.nodes[n]['s']='o'
            plt.rc('figure',figsize=figsize)
            label_options = {"ec": "k", "fc": "white", "alpha": 0.8}
            edges=H.edges()
            weights = [H[u][v]['weight'] for u,v in edges]
            labels = nx.get_edge_attributes(H,'weight')
            labels={k:round(v,2) for k,v in labels.items()}

            all_node=self.goplist[sol_num]+self.gmatlist[sol_num]
                        
            labels1={}
            for x in all_node:
                string=H.nodes()[x]['names']
                labels1.update({x:string})
            
            for e in H.edges():
                if e[0] in all_node and e[1] in all_node:
                    H[e[0]][e[1]]['color']='black'
                else:
                    H[e[0]][e[1]]['color']='lightgrey'
                    
            for n in H.nodes():
                if n in all_node:
                    H.nodes()[n]['color']='black'
                else:
                    H.nodes()[n]['color']='lightgrey'
                if n[0]=='M':
                    H.nodes()[n]['shape']='o'
                elif n[0]=='O':
                    H.nodes()[n]['shape']='s'

            shape_list=[H.nodes()[n]['shape'] for n in H.nodes()]  
            node_color_list=[H.nodes()[n]['color'] for n in H.nodes()]        
            edge_color_list=[H[e[0]][e[1]]['color'] for e in H.edges()]    
            
            nodeShapes = set((aShape[1]["s"] for aShape in H.nodes(data = True)))

            pos=pydot_layout(H,prog='dot')
            pos2=pydot_layout(H,prog='dot')
            

            for key, (v1,v2) in pos2.items():
                if key[0]=="O":
                    pos2[key]=(v1,v2-3)
            pos=nx.rescale_layout_dict(pos,scale=rescale)
            pos2=nx.rescale_layout_dict(pos2,scale=rescale)
            nx.draw_networkx(H, pos=pos,labels=labels1, node_color='white',alpha=0.9,node_shape='o', edge_color=edge_color_list, with_labels = True,node_size=node_size,bbox=label_options,width=weights,font_size=10)
            for aShape in nodeShapes:
                node_list=[sNode[0] for sNode in filter(lambda x: x[1]["s"]==aShape,H.nodes(data = True))]
                if aShape=="o":
                    node_size1=node_size
                else:
                    node_size1=node_size*2
                for node in node_list:
                    nx.draw_networkx_nodes(H,pos=pos2,node_color=H.nodes()[node]['color'],node_shape = aShape, nodelist =[node] ,node_size=node_size1) 
            
            for key, values in nx.get_node_attributes(H,'type').items():
                if values=="raw_material":
                    nx.draw_networkx_nodes(H,pos=pos,node_color='white',node_shape = 'v', nodelist = [key],node_size=node_size/3*1.6)
                elif values=="product":
                    nx.draw_networkx_nodes(H,pos=pos,node_color='white',node_shape = 'o', nodelist = [key],node_size=node_size/3*2)
                    nx.draw_networkx_nodes(H,pos=pos,node_color=H.nodes()[key]['color'],node_shape = 'o', nodelist = [key],node_size=node_size/3*1.25)
                    nx.draw_networkx_nodes(H,pos=pos,node_color='white',node_shape = 'o', nodelist = [key],node_size=node_size/3*0.75)
            
            nx.draw_networkx_edge_labels(H, pos=pos,edge_labels=labels)
            ax= plt.gca()
            plt.axis('off')
            ax.autoscale()

            ax.set_xlim([ax.get_xlim()[0]+padding*ax.get_xlim()[0],ax.get_xlim()[1]+padding*ax.get_xlim()[1]])
            ax.set_ylim([ax.get_ylim()[0]+padding*ax.get_ylim()[0],ax.get_ylim()[1]+padding*ax.get_ylim()[1]])
            if box:
                ax.set_aspect('equal', adjustable='box')

            sol_id = self.goolist[sol_num]
            if sol_id=="0":
                ax.set_title("Maximal Structure",y=titlepos)  
            else:
                ax.set_title("Solution Structure #"+sol_id,y=titlepos)
   
            plt.tight_layout()

        
        return ax
    
    def to_studio(self, path=None,file_name="studio_file.pgsx",verbose=False):
        '''
        to_studio(path=None,file_name="studio_file.pgsx",verbose=False)
        
        Description
        Convert the current status of the problem (with potential solutions) to .pgsx (P-graph Studio File). 
        
        Arguments
        path: (string)(optional) Path of the expected output file. Default path is the directory of the file.
        file_name: (string)(optional) Name of the file. Default is "studio_file.pgsx"
        verbose: (boolean) Whether to print the content of the file.       
        
        '''
        
        G=self.G
        gmatlist=self.gmatlist
        goplist=self.goplist
        goolist=self.goolist
        ME= self.ME
        if path==None:
            path= os.path.dirname(__file__)+"/"
        header=r'<?xml version="1.0" encoding="utf-16"?>'+"\n"
        xsi="http://www.w3.org/2001/XMLSchema-instance"
        xsd="http://www.w3.org/2001/XMLSchema"
        Type="PNS"
        Visible="true"
        PeriodExtension="true"
        NSMAP={"xsi":xsi,"xsd":xsd}
        attrib={"Type":Type,"Visible":Visible,"PeriodExtension":PeriodExtension}
        PG="PGraph"

        PGraph=etree.Element(PG,nsmap=NSMAP,attrib=attrib)

        #Default
        Default=etree.SubElement(PGraph,"Default")

        ##Default Material
        Def_Material=etree.SubElement(Default,"Material")
        DM_FRLB=etree.SubElement(Def_Material,"FlowRateLowerBound")
        DM_FRLB.text="0"
        DM_FRUP=etree.SubElement(Def_Material,"FlowRateUpperBound")
        DM_FRUP.text="1000000000"
        DM_Price=etree.SubElement(Def_Material,"Price")
        DM_Price.text="0"
        DM_Type=etree.SubElement(Def_Material,"Type")
        DM_Type.text="1"
        DM_Deadline=etree.SubElement(Def_Material,"Deadline")
        DM_Deadline.text="31536000"
        DM_EA=etree.SubElement(Def_Material,"EarliestAvability")
        DM_EA.text="0"
        DM_Storage=etree.SubElement(Def_Material,"StorageStrategy")
        DM_Storage.text="default"


        ## Operating Unit
        Def_Op=etree.SubElement(Default,"OperatingUnit")
        DO_OFC=etree.SubElement(Def_Op,"OperatingFixCost")
        DO_OFC.text="0"
        DO_IFC=etree.SubElement(Def_Op,"InvestmentFixCost")
        DO_IFC.text="0"
        DO_OpFC=etree.SubElement(Def_Op,"OpUnitFixCost")
        DO_OpFC.text="0"
        DO_OPC=etree.SubElement(Def_Op,"OperatingPropCost")
        DO_OPC.text="0"
        DO_IPC=etree.SubElement(Def_Op,"InvestmentPropCost")
        DO_IPC.text="0"
        DO_OpUPC=etree.SubElement(Def_Op,"OpUnitPropCost")
        DO_OpUPC.text="0"
        DO_CLB=etree.SubElement(Def_Op,"CapacityLowerBound")
        DO_CLB.text="0"
        DO_CUB=etree.SubElement(Def_Op,"CapacityUpperBound")
        DO_CUB.text="1000000000"
        DO_PP=etree.SubElement(Def_Op,"PaybackPeriod")
        DO_PP.text="10"
        DO_WHPY=etree.SubElement(Def_Op,"WorkingHoursPerYear")
        DO_WHPY.text="8000"
        DO_FOT=etree.SubElement(Def_Op,"FixOperTime")
        DO_FOT.text="0"
        DO_POT=etree.SubElement(Def_Op,"PropOperTime")
        DO_POT.text="0"
        DO_EA=etree.SubElement(Def_Op,"EarliestAvability")
        DO_EA.text="0"
        DO_LA=etree.SubElement(Def_Op,"LatestAvability")
        DO_LA.text="31536000"
        DO_RM=etree.SubElement(Def_Op,"RelaxMode")
        DO_RM.text="strong"

        ## Edge
        Def_E=etree.SubElement(Default,"Edge")
        DE_FR=etree.SubElement(Def_E,"FlowRate")
        DE_FR.text="1"

        ## Quantity
        Def_Q=etree.SubElement(Default,"Quantity")
        DQ_DM=etree.SubElement(Def_Q,"default_mes")
        DQ_DM.text="gram (g)"
        DQ_TM=etree.SubElement(Def_Q,"time_mu")
        DQ_TM.text="y"
        DQ_QT=etree.SubElement(Def_Q,"quant_type")
        DQ_QT.text="Mass"
        DQ_MM=etree.SubElement(Def_Q,"money_mu")
        DQ_MM.text="EUR"

        ## SolverParameter
        DEF_SP=etree.SubElement(Default,"SolverParameter")
        DS_MC=etree.SubElement(DEF_SP,"MakespanCoefficient")
        DS_MC.text="0"
        DS_CC=etree.SubElement(DEF_SP,"CostCoefficient")
        DS_CC.text="1"
        DS_TC=etree.SubElement(DEF_SP,"TimeCoefficient")
        DS_TC.text="0"
        DS_TCLB=etree.SubElement(DEF_SP,"TotalCostLowerBound")
        DS_TCLB.text="-1000000000"
        DS_TCUB=etree.SubElement(DEF_SP,"TotalCostUpperBound")
        DS_TCUB.text="1000000000"

        # Materials
        Materials=etree.SubElement(PGraph,"Materials")

        ## Material
        type_converter={"raw_material":0,"intermediate":1,"product":2}

        Mats_list=[]
        MPar_list=[]
        for n in G.nodes():
            if n[0]=="M":
                Mats_list.append(n)
                if 'type' in G.nodes()[n]:
                    type_ind=type_converter[G.nodes()[n]['type']]
                    attr={"ID":str(list(G.nodes()).index(n)+1),"Name":G.nodes()[n]["names"],"Type":str(type_ind)}
                    Mats_list.append(etree.SubElement(Materials,'Material',attrib=attr))
                    MPar_list.append(etree.SubElement(Mats_list[-1],'ParameterList'))
                    if 'price' in G.nodes()[n]:
                        price=G.nodes()[n]['price']
                    else:
                        price=-1
                    if 'flow_rate_lower_bound' in G.nodes()[n]:
                        flow_rate_lower_bound=G.nodes()[n]['flow_rate_lower_bound']
                    else:
                        flow_rate_lower_bound=-1            
                    
                    if 'flow_rate_upper_bound' in G.nodes()[n]:
                        flow_rate_upper_bound=G.nodes()[n]['flow_rate_upper_bound']
                    else:
                        flow_rate_upper_bound=-1             
                    
                    etree.SubElement(MPar_list[-1],'Parameter',attrib={"Name":"price", "Prefix":"Price: ", "Value":str(price), "Visible":"false"})
                    etree.SubElement(MPar_list[-1],'Parameter',attrib={"Name":"reqflow", "Prefix":"Required flow: ", "Value":str(flow_rate_lower_bound), "Visible":"false"})
                    etree.SubElement(MPar_list[-1],'Parameter',attrib={"Name":"maxflow", "Prefix":"Maximum flow: ", "Value":str(flow_rate_upper_bound), "Visible":"false"})
                    etree.SubElement(MPar_list[-1],'Parameter',attrib={"Name":"quantitytype", "Prefix":"Quantity type: ", "Value":"Mass", "Visible":"false"})
                    etree.SubElement(MPar_list[-1],'Parameter',attrib={"Name":"measurementunit", "Prefix":"Measurement unit: ", "Value":"gram (g)", "Visible":"false"})
                    
                    etree.SubElement(Mats_list[-1], 'Label',attrib={"Text":G.nodes()[n]["names"]})
        # Edges
        Edges=etree.SubElement(PGraph,"Edges")

        ## Edge
        global_edge_count=len(list(G.nodes()))+1
        edge_list=[]
        nodes_list=[]
        label_list=[]
        offset_list=[]
        X_list=[]
        Y_list=[]
        for n in G.nodes():
            if n[0]=="O":
                ee=G.in_edges(n)
                for e in ee:
                    ratio=G[e[0]][e[1]]['weight']
                    attr={"ID":str(global_edge_count),"BeginID":str(list(G.nodes()).index(e[0])+1),"EndID":str(list(G.nodes()).index(e[1])+1),"Rate":str(ratio),"Title":str(ratio), "ArrowOnCenter":"true","ArrowPosition":"50"}
                    edge_list.append(etree.SubElement(Edges,'Edge', attrib=attr))
                    nodes_list.append(etree.SubElement(edge_list[-1],'Nodes')) #Hanging
                    label_list.append(etree.SubElement(edge_list[-1],'Label',attrib={"Text":str(ratio)}))
                    offset_list.append(etree.SubElement(label_list[-1],'Offset'))
                    X_list.append(etree.SubElement(offset_list[-1],"X"))
                    Y_list.append(etree.SubElement(offset_list[-1],"Y"))
                    X_list[-1].text="0"
                    Y_list[-1].text="0"
                    etree.SubElement(label_list[-1],"FontSize").text="-1"
                    etree.SubElement(label_list[-1],"Color").text="-16777216"
                    etree.SubElement(edge_list[-1],'Color').text="-16777216"
                    etree.SubElement(edge_list[-1],'LongFormat').text="false"
                    etree.SubElement(edge_list[-1],'Comment') #hanging
                    etree.SubElement(edge_list[-1],'CommentVisible').text="false"
                    
                    global_edge_count+=1
                ff=G.out_edges(n)
                for f in ff:
                    ratio=-G[f[0]][f[1]]['weight']
                    attr={"ID":str(global_edge_count),"BeginID":str(list(G.nodes()).index(f[0])+1),"EndID":str(list(G.nodes()).index(f[1])+1),"Rate":str(ratio),"Title":str(ratio), "ArrowOnCenter":"true","ArrowPosition":"50"}
                    edge_list.append(etree.SubElement(Edges,'Edge', attrib=attr))
                    nodes_list.append(etree.SubElement(edge_list[-1],'Nodes')) #Hanging
                    label_list.append(etree.SubElement(edge_list[-1],'Label',attrib={"Text":str(ratio)}))
                    offset_list.append(etree.SubElement(label_list[-1],'Offset'))
                    X_list.append(etree.SubElement(offset_list[-1],"X"))
                    Y_list.append(etree.SubElement(offset_list[-1],"Y"))
                    X_list[-1].text="0"
                    Y_list[-1].text="0"
                    etree.SubElement(label_list[-1],"FontSize").text="-1"
                    etree.SubElement(label_list[-1],"Color").text="-16777216"
                    etree.SubElement(edge_list[-1],'Color').text="-16777216"
                    etree.SubElement(edge_list[-1],'LongFormat').text="false"
                    etree.SubElement(edge_list[-1],'Comment') #hanging
                    etree.SubElement(edge_list[-1],'CommentVisible').text="false"
                    global_edge_count+=1

        # Operating Units
        OperatingUnits=etree.SubElement(PGraph,"OperatingUnits")
         
        ## Operating Unit
        opu_list=[]
        OPar_list=[]
        for n in G.nodes():
            if n[0]=="O":
                attr={"ID":str(list(G.nodes()).index(n)+1),"Name":G.nodes()[n]["names"],"Title":""}
                opu_list.append(etree.SubElement(OperatingUnits,'OperatingUnit', attrib=attr))
                OPar_list.append(etree.SubElement(opu_list[-1],'ParameterList'))
                if 'capacity_lower_bound' in G.nodes()[n]:
                    capacity_lower_bound=G.nodes()[n]['capacity_lower_bound']
                else:
                    capacity_lower_bound=-1
                if 'capacity_upper_bound' in G.nodes()[n]:
                    capacity_upper_bound=G.nodes()[n]['capacity_upper_bound']
                else:
                    capacity_upper_bound=-1                
                if 'fix_cost' in G.nodes()[n]:
                    fix_cost=G.nodes()[n]['fix_cost']
                else:
                    fix_cost=-1 
                if 'proportional_cost' in G.nodes()[n]:
                    proportional_cost=G.nodes()[n]['proportional_cost']
                else:
                    proportional_cost=-1 
                
                etree.SubElement(OPar_list[-1],'Parameter',attrib={"Name":"caplower", "Prefix":"Capacity, lower bound: ", "Value":str(capacity_lower_bound), "Visible":"false"})
                etree.SubElement(OPar_list[-1],'Parameter',attrib={"Name":"capupper", "Prefix":"Capacity, upper bound: ", "Value":str(capacity_upper_bound), "Visible":"false"})
                etree.SubElement(OPar_list[-1],'Parameter',attrib={"Name":"investcostfix", "Prefix":"Investment cost, fix: ", "Value":str(fix_cost), "Visible":"false"})
                etree.SubElement(OPar_list[-1],'Parameter',attrib={"Name":"investcostprop", "Prefix":"Investment cost, proportional: ", "Value":str(proportional_cost), "Visible":"false"})
                etree.SubElement(OPar_list[-1],'Parameter',attrib={"Name":"opercostfix", "Prefix":"Operating cost, fix: ", "Value":"-1", "Visible":"false"})
                etree.SubElement(OPar_list[-1],'Parameter',attrib={"Name":"opercostprop", "Prefix":"Operating cost, proportional: ", "Value":"-1", "Visible":"false"})
                etree.SubElement(OPar_list[-1],'Parameter',attrib={"Name":"workinghour", "Prefix":"Working hours per year: ", "Value":"-1", "Visible":"false"})
                etree.SubElement(OPar_list[-1],'Parameter',attrib={"Name":"payoutperiod", "Prefix":"Payout Period: ", "Value":"-1", "Visible":"false"})
                
                etree.SubElement(opu_list[-1], 'Label',attrib={"Text":G.nodes()[n]["names"]})
        # MutualExclusions
        MutualExclusions=etree.SubElement(PGraph,"MutualExclusions")

        ## MutualExclusion
        ME_list=[]
        MOP_list=[]
        
        for M in ME:
            #New version P-graph studio change
            nameM=[G.nodes()[mm]['names'].replace(" ", "_") for mm in M]
            Name="_".join(nameM)+'_'+"_".join(M)
            
            
            #Name="_".join(M)
            
            attr={"ID":str(global_edge_count),"Name":Name}
            ME_list.append(etree.SubElement(MutualExclusions,"MutualExclusion",attrib=attr))
            MOP_list.append(etree.SubElement(ME_list[-1],"OperatingUnits"))
            #for x in M:
            for x in nameM:
                etree.SubElement(MOP_list[-1],"OperatingUnit").text=x
            global_edge_count+=1
        if self.solver in ["SSGLP","INSIDEOUT",2,3]:
            # Solutions
            Solutions=etree.SubElement(PGraph,"Solutions")

            ## Solution 
            ssol_list=[]
            smats_list=[]
            sops_list=[]
            sop_list=[]

            for i in range(len(goolist)):
                snum=i+1
                attr={"Index":str(i), "Title":"Solution #"+str(snum), "OptimalValue":str(goolist[i]),"TotalTime":"0", "TotalMakespan":"0", "ObjectiveValue":str(goolist[i]), "AlgorithmUsed":"INSIDEOUT"}
                ssol_list.append(etree.SubElement(Solutions,"Solutions",attrib=attr))
                smats_list.append(etree.SubElement(ssol_list[-1],"Materials"))
                for x in gmatlist[i]:
                    attr={"Name":x[0],"Flow":str(x[3]),"Cost":str(x[1]), "MU":""}
                    etree.SubElement(smats_list[-1],"Material",attrib=attr)
                
                sops_list.append(etree.SubElement(ssol_list[-1],"OperatingUnits"))
                for x in goplist[i]:
                    attr={"Name":x[1],"Size":str(x[0]),"Cost":str(x[2]),"MU":""}
                    sop_list.append(etree.SubElement(sops_list[-1],"OperatingUnit",attrib=attr))
                    etree.SubElement(sop_list[-1],"Input")
                    etree.SubElement(sop_list[-1],"Output")
                
        xml=etree.tostring(PGraph,pretty_print=True,encoding='unicode', method='xml')

        with open(path+file_name,"w") as f:
            f.write(header)
            f.write(xml)
        if verbose:
            print(header+xml)
            print("Generated P-graph Studio File at ", path)
        return header+xml    
        
    def run(self, skip_wine=False, solver_name='pgraph_solver.exe', path=None, input_file: str | None = None, output_file: str | None = None):
        '''
        run(system=None,skip_wine=False)
        
        Description
        Create input, solve problem and read solution.
        
        Arguments
        system: (string) (optional) Operating system. Options of "Windows", "Linux". MacOS is not supported yet. Specifying this makes function slightly faster.
        skip_wine: (boolean) Only relevent for Linux. Skip the dependency "wine" if it is already installed. 
        solver_name= (string) For advanced users only. Choose your customized solver. 'pgraph_solver.exe' or 'pgraph_solver_new.exe'
        solver_path = (string) path to the custom solver. If None, then the default library installation path will be used.
        '''
        system = platform.system()
        
        if system=="Windows":
            self.create_solver_input(input_file=input_file)
            self.solve(system=system,skip_wine=skip_wine,solver_name=solver_name,path=path)
            self.read_solutions()
            # print('chong',path)
        elif system=="Linux":
            if input_file is None:
                # set input file to a temp file
                input_file = tempfile.NamedTemporaryFile(delete=False).name
            if output_file is None:
                # set output file to a temp file
                output_file = tempfile.NamedTemporaryFile(delete=False).name
            self.create_solver_input(input_file=input_file)
            self.solve(
                system=system,
                skip_wine=skip_wine,
                solver_name=solver_name,
                input_file=input_file,
                output_file=output_file,
            )
            self.read_solutions(output_file=output_file)
        
    def get_sol_num(self):
        '''
        get_sol_num()
        
        Description
        Get number of solutions generated by solver
        
        Return
        num_sol: (int) Number of solutions
        '''
        num_sol=len(self.goolist)
        return num_sol
if __name__=="__main__":
    pass
    ##TEST1########################################
    ### Prepare Network Structure #############
    # G = nx.DiGraph()
    # G.add_node("M1",names="Product A",type='product',flow_rate_lower_bound=100)
    # G.add_node("M2",names="Chemical A",type='raw_material',price=200,flow_rate_upper_bound=50)
    # G.add_node("M3",names="Chemical B", type='raw_material',price=100)
    # G.add_node("O1",names="Reactor A",fix_cost=2000, proportional_cost=400)
    # G.add_node("O2", names="Reactor B",fix_cost=1000, proportional_cost=400)
    # G.add_edge("O1","M1", weight = 3) 
    # G.add_edge("O2","M1", weight = 1) 
    # G.add_edge("M2","O1", weight = 2)
    # G.add_edge("M3","O2", weight = 4)
    # ME=[["O1","O2"]]
    # P=Pgraph(problem_network=G, mutual_exclusion=ME, solver="SSG",max_sol=10)
    # ax1=P.plot_problem()
    # plt.show()
    # P.run(solver_name='pgraph_solver.exe')
    # total_sol_num=P.get_sol_num()
    # print(total_sol_num)
    # for i in range(total_sol_num): #show all solutions in plot
    #     ax=P.plot_solution(sol_num=i)
        
    #     plt.show()
    
    # P.to_studio(verbose=True)
    #####################################
    
    '''
    ##### STEP 1 : Problem Specification ######
    G = nx.DiGraph()
    G.add_node("M1",names="Product D",type='product',flow_rate_lower_bound=100, flow_rate_upper_bound=100)
    G.add_node("M2",names="Chemical A",type='raw_material',price=200,flow_rate_lower_bound=0)
    G.add_node("M3",names="Chemical B", type='raw_material',price=100,flow_rate_lower_bound=0)
    G.add_node("M4",names="Chemical C", type='raw_material',price=10,flow_rate_lower_bound=0)
    G.add_node("O1",names="Reactor 1",fix_cost=2000, proportional_cost=400)
    G.add_node("O2", names="Reactor 2",fix_cost=1000, proportional_cost=400)
    G.add_edge("M2","O1", weight = 1)
    G.add_edge("M3","O2", weight = 1)
    G.add_edge("M4","O2", weight = 2)
    G.add_edge("O1","M1", weight = 0.7) 
    G.add_edge("O2","M1", weight = 0.9) 
    ME=[["O1","O2"]]  #Reactor 1 and Reactor 2 are mutually excluded. Only one can be chosen as solution.
    P=Pgraph(problem_network=G, mutual_exclusion=ME, solver="MSG",max_sol=100)
    ax1=P.plot_problem()
    P.run(solver_name='pgraph_solver_new.exe')
    plt.show()
    total_sol_num=P.get_sol_num() 
    for i in range(total_sol_num): # Here we loop through all the solutions to plot everything
        ax=P.plot_solution(sol_num=i) #Plot Solution Function
        plt.show()
    #############################
    '''
    '''
    ## TEST 2 ###################
    from sklearn.datasets import load_diabetes
    from chemsy.predict import *
    from sklearn.preprocessing import MinMaxScaler
    diabetes = load_diabetes()
    X, y = diabetes.data, diabetes.target
    xscaler=MinMaxScaler()
    yscaler=MinMaxScaler()

    X=xscaler.fit_transform(X)
    y=yscaler.fit_transform(y.reshape(-1, 1))

    model=PartialLeastSquaresCV()
    model.fit(X,y)
    B=model.model.coef_ #Get B vector of PLSCV

    ##### STEP 1 : Problem Specification ######
    G = nx.DiGraph()
    G.add_node("M1",names="combine",type='intermediate',flow_rate_lower_bound=0, flow_rate_upper_bound=0,price=0)
    global_count=2
    for i in range(B.shape[0]):
      if B[i][0]>=0:
        G.add_node("M"+str(global_count),names=diabetes.feature_names[i],type='raw_material',flow_rate_lower_bound=0, flow_rate_upper_bound=1,price=0)
        G.add_node("O"+str(global_count),names=diabetes.feature_names[i])
        G.add_edge("M"+str(global_count),"O"+str(global_count), weight = 1, prices=0)
        G.add_edge("O"+str(global_count),"M1", weight = B[i][0])
        global_count=global_count+1
      else:
        G.add_node("M"+str(global_count),names=diabetes.feature_names[i],type='product',flow_rate_lower_bound=0, flow_rate_upper_bound=1,price=0)
        G.add_node("O"+str(global_count),names=diabetes.feature_names[i])
        G.add_edge("O"+str(global_count),"M"+str(global_count), weight = 1, price=0)
        G.add_edge("M1","O"+str(global_count), weight = -B[i][0])
        global_count=global_count+1

    G.add_node("O"+str(global_count),names="sum")  
    G.add_edge("M1","O"+str(global_count), weight = 1)
    G.add_edge("O"+str(global_count),"M"+str(global_count+1), weight = 1)
    global_count=global_count+1
    G.add_node("M"+str(global_count),names="Prediction",type='product',price=100,flow_rate_lower_bound=0, flow_rate_upper_bound=10000)

    #### Step 2:  Setup Solver ####
    P=Pgraph(problem_network=G,  solver="SSG",max_sol=100)
    ###############################
    ax1=P.plot_problem(figsize=(20,20))

    plt.show()
    #### Step 3: Run ####
    P.run()
    ####################
    #### Step 3.1: Plot Solution########
    total_sol_num=P.get_sol_num() 
    for i in range(1): # Here we plot top 3 solutions
        ax=P.plot_solution(sol_num=i,figsize=(20,20)) #Plot Solution Function
        #ax.set_xlim(-100,800)
        #plt.show()
    #####################################  
    
    '''
    '''
    P.to_studio()
    
    
    ##NEW FUNCTIONS #####
    for i in range(total_sol_num):
        HH=P.get_solution_as_network(sol_num=i)
        #print(HH.nodes())
        
    a,b,c=P.get_info()
    
    
    #get_solution_as_network(1)
    '''