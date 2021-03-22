##############################################################################
# Author: Roozbeh Rezakhani
# Email:  rrezakhani@gmail.com
#
# This is the phase_field_model class.
#
##############################################################################

import numpy as np
from src.element import element

class phase_field_model:
    
    def __init__(self, mat, mesh, dim):        
        print('Solid Mechanics Model object is constructed!')
        self.mat  = mat
        self.mesh = mesh
        self.dim  = dim
        
        # construct elements object list
        elem_list = np.array(self.mesh.get_elems()).astype(np.int)
        node_list = np.array(self.mesh.get_nodes()).astype(np.float)
        self.elem_obj_list = []
        for e in range(len(elem_list)): 
            nodes = node_list[elem_list[e][2:]-1][:,:self.dim]
            neN = len(nodes)
            element_gmsh_type = elem_list[e,0]
    
            Le = np.array([])
            LePF = np.array([])
            for i in range(elem_list.shape[1]-2): 
                Le = np.concatenate((Le, [2*elem_list[e][i+2]-1, 2*elem_list[e][i+2]]))
                LePF = np.concatenate((LePF, [elem_list[e][i+2]]))
            Le = Le.astype(np.int)
            LePF = LePF.astype(np.int)
    
            elem_obj = element(element_gmsh_type, neN, self.dim, nodes, Le, LePF)
            self.elem_obj_list.append(elem_obj)
            
    def compute_disp_stiffness_matrix(self, phi):        
        num_nodes_total = self.mesh.get_num_nodes()
        C = self.mat.get_C()
        K = np.zeros((num_nodes_total*self.dim, num_nodes_total*self.dim))
        
        for elem_obj in self.elem_obj_list:                                    
            LePF = elem_obj.get_elem_PF_connectivity()
            phi_elem = phi[LePF-1]
            K_elem = elem_obj.compute_element_stiffness(C, phi_elem)
            Le = elem_obj.get_elem_connectivity()
                
            for i in range(len(Le)):
                for j in range(len(Le)): 
                    K[Le[i]-1][Le[j]-1] = K[Le[i]-1][Le[j]-1] + K_elem[i][j]
        return K
    
    def compute_disp_internal_forces(self, U, dU, phi):
        num_nodes_total = self.mesh.get_num_nodes()
        C = self.mat.get_C()
        F_int = np.zeros(num_nodes_total*self.dim)
        
        for elem_obj in self.elem_obj_list:             
            Le = elem_obj.get_elem_connectivity()
            u_elem  = U[Le-1]
            du_elem = dU[Le-1]
            LePF = elem_obj.get_elem_PF_connectivity()
            phi_elem = phi[LePF-1]
            F_int_elem = elem_obj.compute_element_internal_forces(C, u_elem, du_elem, phi_elem)
                
            for i in range(len(Le)):
                F_int[Le[i]-1] = F_int[Le[i]-1] + F_int_elem[i]

        return F_int
    
    def compute_PF_stiffness_matrix(self):        
        num_nodes_total = self.mesh.get_num_nodes()
        K_phi = np.zeros((num_nodes_total, num_nodes_total))
        
        for elem_obj in self.elem_obj_list:                                    
            LePF = elem_obj.get_elem_PF_connectivity()
            K_phi_elem = elem_obj.compute_element_PF_stiffness(self.mat)
                
            for i in range(len(LePF)):
                for j in range(len(LePF)): 
                    K_phi[LePF[i]-1][LePF[j]-1] = K_phi[LePF[i]-1][LePF[j]-1] + K_phi_elem[i][j]
    
        return K_phi
    
    def compute_PF_residual(self, phi):
        num_nodes_total = self.mesh.get_num_nodes()
        res_phi = np.zeros(num_nodes_total)
        
        for elem_obj in self.elem_obj_list:                                    
            LePF = elem_obj.get_elem_PF_connectivity()
            phi_elem = phi[LePF-1]
            res_phi_elem = elem_obj.compute_element_PF_residual(phi_elem, self.mat)
            
            for i in range(len(LePF)):
                res_phi[LePF[i]-1] = res_phi[LePF[i]-1] + res_phi_elem[0,i]        
        
        return res_phi
    

