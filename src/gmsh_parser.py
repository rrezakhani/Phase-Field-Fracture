##############################################################################
# Author: Roozbeh Rezakhani
# Email:  rrezakhani@gmail.com
#
# 2-dimensional mesh file generate by gmsh file (with version 2)
#
##############################################################################

class gmsh_parser:
    
    def __init__(self, path, dim):
    
        self.path = path
        self.dim = dim
        self.num_phys_names = 0
        self.phys_array = []
        self.num_nodes = 0
        self.nodes = []
        self.num_bndry_elems = 0
        self.bndry_elems = []
        self.num_elems = 0
        self.elems = []
                            
        mesh_file = open(self.path, 'r')
        line = mesh_file.readline()
        while (line != ''):  # breaks when EOF is reached
            line = mesh_file.readline()

            if("$PhysicalNames" in  line):
                self.num_phys_names = mesh_file.readline()
                for i in range(int(self.num_phys_names)):
                    line = mesh_file.readline()
                    self.phys_array.append([line.split(' ')[0],
                                            line.split(' ')[1],
                                            line.split(' ')[2].split('\n')[0][1:-1]]) 
                
            if("$Nodes" in  line):
                self.num_nodes = mesh_file.readline()
                for i in range(int(self.num_nodes)):
                    line = mesh_file.readline()
                    if(self.dim == 2):
                        self.nodes.append([line.split(' ')[1],
                                           line.split(' ')[2],
                                           line.split(' ')[3]]) 
                    elif(self.dim == 3):
                        self.nodes.append([line.split(' ')[1],
                                           line.split(' ')[2],
                                           line.split(' ')[3]]) 
                
            if("$Elements" in  line):
                num_elems_total = mesh_file.readline()
                for i in range(int(num_elems_total)):
                    line = mesh_file.readline()
                    if(self.dim == 2):
                        if (line.split(' ')[1] == '1'):
                            self.bndry_elems.append([line.split(' ')[1], # gmsh element type
                                                     line.split(' ')[3], # phys. tag 
                                                     line.split(' ')[-2],
                                                     line.split(' ')[-1].split('\n')[0]])
                        elif (line.split(' ')[1] == '3'):
                            self.elems.append([line.split(' ')[1], # gmsh element type
                                               line.split(' ')[3], # phys. tag
                                               line.split(' ')[-4],
                                               line.split(' ')[-3],
                                               line.split(' ')[-2],
                                               line.split(' ')[-1]])
        self.num_bndry_elems = len(self.bndry_elems) 
        self.num_elems = len(self.elems)        
               
                
    def get_dim(self):
        return self.dim
    
    def get_nodes(self):
        return self.nodes
    
    def get_num_nodes(self):
        return int(self.num_nodes)
    
    def get_elems(self):
        return self.elems
    
    def get_num_elems(self):
        return int(self.num_elems)
    
    def get_bndry_elems(self):
        return self.bndry_elems
        
    def get_phys_array(self):
        return self.phys_array
    