##############################################################################
# Author: Roozbeh Rezakhani
# Email:  rrezakhani@gmail.com
#
# This is the main file, which serves as the driver of the finite element 
# simulation. 
#
##############################################################################
import time
import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve

from src.phase_field_model import phase_field_model
from src.gmsh_parser import gmsh_parser
from src.material import material
from src.vtk_writer import vtk_writer
from termcolor import colored

##############################################################################
# Read the input file 
input_file = open("./input-PhaseField.in", 'r')
line = input_file.readline()
disp_BC = []
trac_BC = []
while (line != ''):  # breaks when EOF is reached
    key = line.split(' ')[0]
    if (key == 'dim'):
        dim = int(line.split(' ')[-1])
    if (key == 'two_dimensional_problem_type'):
        two_dimensional_problem_type = line.split(' ')[-1].split('\n')[0]
    if (key == 'E'):
        E = float(line.split(' ')[-1])
    if (key == 'nu'):
        nu = float(line.split(' ')[-1])    
    if (key == 'Gc'):
        Gc = float(line.split(' ')[-1])
    if (key == 'el'):
        el = float(line.split(' ')[-1])         
    if (key == 'mesh_path'):
        mesh_path = line.split(' ')[-1].split('\n')[0]
    if (key == 'Dirichlet_BC'):
        disp_BC.append([line.split(' ')[-3][1:], 
                        line.split(' ')[-2], 
                        line.split(' ')[-1].split('\n')[0][1:-2]])
    if (key == 'analysis_type'):
        analysis_type = line.split(' ')[-1].split('\n')[0]
    if (key == 'solver'):
        solver = line.split(' ')[-1].split('\n')[0]
    line = input_file.readline()    
    
out_dir = "./examples/"+mesh_path.split('/')[-1][:-4]

##############################################################################
# Create the mesh object
mesh = gmsh_parser(mesh_path, dim)
num_nodes = mesh.get_num_nodes()

##############################################################################
# Instantiate material class and initialize material properties
mat = material(E, nu, Gc, el, dim, two_dimensional_problem_type)
C = mat.get_C()

##############################################################################
# Instantiate model class to build the general framework
phase_field_model = phase_field_model(mat, mesh, dim)

##############################################################################
# Initialize global displacement vector
U = np.zeros(num_nodes*dim)
dU = np.zeros(num_nodes*dim)

##############################################################################
# Boundary condition arrays
blocked = np.zeros(num_nodes*dim)
u_bar = np.zeros(num_nodes*dim)

# Apply prescribed boundary conditions
phys_array = mesh.get_phys_array()
bndry_elems = np.array(mesh.get_bndry_elems()).astype(np.int)   
 
for d_BC in disp_BC:
    val = float(d_BC[0])
    comp = int(d_BC[1])
    phys_tag = d_BC[2]
    for i in range(len(phys_array)):
        if(phys_array[i][2] == phys_tag):
            phys_index = int(phys_array[i][1])
    phys_elems = bndry_elems[bndry_elems[:,1]==phys_index][:,2:]
    blocked[2*(np.unique(phys_elems)-1)+comp-1] = 1
    u_bar[2*(np.unique(phys_elems)-1)+comp-1] = val

##############################################################################   
# External force vector
F_ext = np.zeros(num_nodes*dim)

# Reaction force vector
R_ext = np.zeros(num_nodes*dim)

# Internal force vector
F_int = np.zeros(num_nodes*dim)

##############################################################################   
# Construct PF (phase field) related variables
phi = np.zeros(num_nodes)
K_phi = np.zeros((num_nodes, num_nodes))
res_phi = np.zeros(num_nodes)

##############################################################################   
# Iterative solve of the equilibruim equation
# Newton Raphson Method
num_load_steps = 30
itr_max = 10
itr_tol = 1E-5
F_ext_ = np.zeros(num_nodes*dim)
u_bar_ = np.zeros(num_nodes*dim)
res = np.zeros(num_nodes*dim)

vtk_writer(out_dir, '0', mesh, U, phi, F_int)

# starting time
start = time.time()

for l in range(num_load_steps):
    
    print("Load step {}:".format(l+1))
    
    #=====================================================================
    # update essential and neumann boundary conditions
    u_bar_ += 1/num_load_steps * u_bar
    F_ext_ += 1/num_load_steps * F_ext
    
    #=====================================================================
    # update stiffness matrix of the displacement governing equation
    K = phase_field_model.compute_disp_stiffness_matrix(phi)
    
    #=====================================================================
    # compute residual of the displacement governing equation
    res = F_ext_ + R_ext - F_int
    
    #=====================================================================
    # Impose essential boundary conditions
    K_temp = np.copy(K)
    for m in range(num_nodes*dim):
        if(blocked[m]==1):
            for n in range(num_nodes*dim):
                if(blocked[n]==0):
                    res[n] -= K[n,m] * u_bar_[m]
            K_temp[m,:] = 0.0
            K_temp[:,m] = 0.0
            K_temp[m,m] = 1.0
            res[m] = u_bar_[m]
    
    #=====================================================================
    # Loop on iterations for the DISPLACEMENT governing equation
    max_itr_reached = False
    for k in range(itr_max):
        
        # calculate displacement vector increment
        #U = np.dot(inv(K_temp), res) 
        
        # calculate displacement vector increment using scipy csc matrix
        K_temp_csc = csc_matrix(K_temp)
        U = spsolve(K_temp_csc, res)
        
        # update the reaction forces
        for m in range(num_nodes*dim):
            if(blocked[m]==1):
                R_ext[m] = np.dot(K[m,:], U)       
               
        # compute internal forces vector
        F_int = phase_field_model.compute_disp_internal_forces(U, dU, phi)
        
        # update residual and check convergence
        res = F_ext_ + R_ext - F_int        
        tol = np.linalg.norm(res)/np.linalg.norm(F_ext_ + R_ext)
        
        print("Displacement solve: Iteration {} - tolerance = {}".format(k+1, tol))
        if (tol < itr_tol):
            print(colored("Displacement solve converged!",'green'))
            break 
        if (k == itr_max-1):
            print(colored("Maximum number of iteration in displacement solve is reached!",'red'))
            print(colored("Displacement solve did NOT converge!",'red'))
            max_itr_reached = True
            break 
    if(max_itr_reached):
        break 
      
    #=====================================================================
    # compute residual of the PF (phase field) governing equation
    res_phi = phase_field_model.compute_PF_residual(phi)

    #=====================================================================
    # update stiffness matrix of the PF (phase field) governing equation
    K_phi = phase_field_model.compute_PF_stiffness_matrix()  
    
    #=====================================================================
    # Loop on iterations for the PF (phase field) governing equation
    for k in range(itr_max):
        
        # calculate phase field vector increment
        #dphi = np.dot(inv(K_phi), -res_phi) 
        
        # calculate phase field vector increment using scipy csc matrix
        K_phi_csc = csc_matrix(K_phi)
        dphi = spsolve(K_phi_csc, -res_phi)
        
        # update the phase field vector
        phi = phi + dphi
               
        # update residual and check convergence
        res_phi= phase_field_model.compute_PF_residual(phi)        
        tol = np.linalg.norm(res_phi)
        
        print("Phase field solve: Iteration {} - tolerance = {}".format(k+1, tol))
        if (tol < itr_tol):
            print(colored("Phase field solve converged!",'green'))
            break 
        if (k == itr_max-1):
            print(colored("Maximum number of iteration in displacement solve is reached!",'red'))
            print(colored("Phase field solve did NOT converge!",'red'))
            max_itr_reached = True
            break 
    if(max_itr_reached):
        break 
    
    #=====================================================================
    # write vtk file   
    vtk_writer(out_dir, l+1, mesh, U, phi, F_int)
 
# end time
end = time.time()

# total time taken
print(f"Runtime of the program is {end - start}")
    
    
    
    
    
    