##############################################################################
# Author: Roozbeh Rezakhani
# Email:  rrezakhani@gmail.com
#
# This is the Element class.
#
##############################################################################

import numpy as np
from numpy.linalg import inv
from numpy.linalg import det

class element:
    
        def __init__(self, element_gmsh_type, neN, dim, nodes, Le, LePF):
            if (element_gmsh_type == 1):
                self.elem_type = '2node-line'
            if (element_gmsh_type == 2):
                self.elem_type = '3node-triangular'
                self.nqp = 1
            if (element_gmsh_type == 3):
                self.elem_type = '4node-quadrilateral'
                self.nqp = 4
            self.neN = neN
            self.dim = dim
            self.nodes = nodes
            self.Le = Le
            self.K_elem = np.zeros((self.neN*self.dim, self.neN*self.dim))
            self.F_int_elem = np.zeros(self.neN*self.dim)
            self.elas_str_ene = np.zeros(self.nqp)
            if (self.dim == 2):
                self.stress_total = np.zeros(self.nqp*3)
                self.strain_total = np.zeros(self.nqp*3)
                self.stress_incr  = np.zeros(self.nqp*3)
                self.strain_incr  = np.zeros(self.nqp*3)                
            else:
                self.stress_total = np.zeros(self.nqp*6)
                self.strain_total = np.zeros(self.nqp*6)
                self.stress_incr  = np.zeros(self.npp*6)
                self.strain_incr  = np.zeros(self.nqp*6)
            self.LePF = LePF    
            self.K_phi_elem = np.zeros((self.neN, self.neN))
            self.res_phi_elem = np.zeros(self.neN)
        
        def compute_element_stiffness(self, C, phi_elem):   
            self.K_elem = np.zeros((self.neN*self.dim, self.neN*self.dim))
            w_qp, qp = self.get_w_qp()
                                
            for p in range(self.nqp):             
                xi  = qp[p][0]
                eta = qp[p][1]
                
                J = self.compute_jacobian(xi, eta)         
                B = self.compute_B_matrix(xi, eta)
                N = self.get_shape_function(xi, eta)
                
                # update material property matrix
                C_degraded = (1-np.dot(np.transpose(N),phi_elem))**2 * C
                
                K_elem_qp = np.dot(np.dot(np.transpose(B), C_degraded), B) * w_qp[p] * det(J)
                self.K_elem = self.K_elem + K_elem_qp  
                
            return self.K_elem
        
        def compute_element_internal_forces(self, C, u_elem, du_elem, phi_elem):
            self.F_int_elem = np.zeros(self.neN*self.dim)
            w_qp, qp = self.get_w_qp()
                                
            for p in range(self.nqp):            
                xi  = qp[p][0]
                eta = qp[p][1]
                
                J = self.compute_jacobian(xi, eta)         
                B = self.compute_B_matrix(xi, eta)
                N = self.get_shape_function(xi, eta)
                
                # update material property matrix
                C_degraded = (1-np.dot(np.transpose(N),phi_elem))**2 * C
                
                # compute element total and incremental strain
                self.strain_incr[p*3:(p+1)*3]  = np.dot(B, du_elem)
                self.strain_total[p*3:(p+1)*3] = np.dot(B, u_elem)
                
                # compute element total and incremental stress
                self.stress_incr[p*3:(p+1)*3]  = np.dot(C_degraded, self.strain_incr[p*3:(p+1)*3])
                self.stress_total[p*3:(p+1)*3] = np.dot(C_degraded, self.strain_total[p*3:(p+1)*3])
                
                # compute elastic strain energy at the current qp
                self.elas_str_ene[p] = 0.5 * np.dot(np.dot(self.strain_total[p*3:(p+1)*3], C), 
                                                    self.strain_total[p*3:(p+1)*3])
                
                # compute element internal forces at current qp
                F_int_elem_qp = np.dot(np.transpose(B), self.stress_total[p*3:(p+1)*3]) \
                                * w_qp[p] * np.linalg.det(J)
                self.F_int_elem = self.F_int_elem + F_int_elem_qp    
                
            return self.F_int_elem
            
        def compute_jacobian(self, xi, eta):   
            dNdxi = self.get_iso_shape_function_derivative(xi, eta) 
            J = np.dot(np.transpose(dNdxi), self.nodes)

            return J
        
        def compute_shape_function_derivate(self, xi, eta):
            dNdxi = self.get_iso_shape_function_derivative(xi, eta)   
            J = self.compute_jacobian(xi, eta)
            invJ = inv(J)
            dNdx = np.dot(invJ, np.transpose(dNdxi))
            
            return dNdx
        
        def compute_B_matrix(self, xi, eta): 
            dNdx = self.compute_shape_function_derivate(xi, eta)
            B = np.zeros((3, 2*self.neN))
            B[0, 0:2*self.neN+1:2] = dNdx[0,:]
            B[1, 1:2*self.neN+1:2] = dNdx[1,:]
            B[2, 0:2*self.neN+1:2] = dNdx[1,:]
            B[2, 1:2*self.neN+1:2] = dNdx[0,:]    
            
            return B
        
        def compute_B_PF_matrix(self, xi, eta):   
            dNdx = self.compute_shape_function_derivate(xi, eta)
            B_PF = dNdx  
            
            return B_PF
        
        def get_shape_function(self, xi, eta):
            # Isoparametric quadrilateral element shape functions
            N = 1/4 * np.array([(1-xi)*(1-eta), 
                                (1+xi)*(1-eta), 
                                (1+xi)*(1+eta), 
                                (1-xi)*(1+eta)])
            
            return N
        
        def get_iso_shape_function_derivative(self, xi, eta):
            # Derivate of isoparametric quadrilateral element shape functions                   
            dNdxi = 1/4 * np.array([[-(1-eta), -(1-xi)],
                                    [ (1-eta), -(1+xi)],
                                    [ (1+eta),  (1+xi)],
                                    [-(1+eta),  (1-xi)]])
            return dNdxi
        
        def get_w_qp(self):
            w_qp = np.array([1, 1, 1, 1])
            qp = np.array([[ 0.5774,  0.5774],
                            [ 0.5774, -0.5774],
                            [-0.5774,  0.5774],
                            [-0.5774, -0.5774]])
            return w_qp, qp
        
        def compute_qp_stress(self, u_elem, du_elem):
            #RR link to material class goes here!
            return u_elem    
        
        def get_elem_connectivity(self):
            return self.Le
        
        def get_elem_PF_connectivity(self):
            return self.LePF
            
        def compute_element_PF_residual(self, phi_elem, mat):
            self.res_phi_elem = np.zeros(self.neN)
            w_qp, qp = self.get_w_qp()
            
            Gc = mat.Gc
            el = mat.el
                                
            for p in range(self.nqp):            
                xi  = qp[p][0]
                eta = qp[p][1]
                
                J = self.compute_jacobian(xi, eta)
                B_PF = self.compute_B_PF_matrix(xi, eta)
                N = self.get_shape_function(xi, eta)
                N = N[np.newaxis]
                
                # elastic strain energy at the current qp
                Eels = self.elas_str_ene[p]

                # compute element internal forces at current qp
                res_phi_elem_qp = (Gc/el * np.dot(np.dot(np.transpose(N), N), phi_elem) + \
                                   Gc*el * np.dot(np.dot(np.transpose(B_PF), B_PF), phi_elem) + \
                                   2*Eels * np.dot(np.dot(np.transpose(N), N), phi_elem) - \
                                   2*Eels * N) \
                                   * w_qp[p] * np.linalg.det(J)
                                   
                self.res_phi_elem = self.res_phi_elem + res_phi_elem_qp    
                
            return self.res_phi_elem
        
        def compute_element_PF_stiffness(self, mat):  
            self.K_phi_elem = np.zeros((self.neN, self.neN))
            w_qp, qp = self.get_w_qp()
            
            Gc = mat.Gc
            el = mat.el

            for p in range(self.nqp):             
                xi  = qp[p][0]
                eta = qp[p][1]
                
                J = self.compute_jacobian(xi, eta)
                B_PF = self.compute_B_PF_matrix(xi, eta)
                N = self.get_shape_function(xi, eta)
                N = N[np.newaxis]
                
                # elastic strain energy at the current qp
                Eels = self.elas_str_ene[p]
                
                K_phi_elem_qp = (Gc/el * np.dot(np.transpose(N), N) + \
                                 Gc*el * np.dot(np.transpose(B_PF), B_PF) + \
                                 2*Eels * np.dot(np.transpose(N), N)) \
                                 * w_qp[p] * np.linalg.det(J)
                self.K_phi_elem = self.K_phi_elem + K_phi_elem_qp  
                
            return self.K_phi_elem