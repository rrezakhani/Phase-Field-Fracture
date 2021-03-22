##############################################################################
# Author: Roozbeh Rezakhani
# Email:  rrezakhani@gmail.com
#
# 2-dimensional mesh file generate by gmsh file (with version 2)
#
##############################################################################

import numpy as np

def vtk_writer(path, l, mesh, U, phi, F_int):   

    dim = mesh.get_dim()    
    num_nodes = mesh.get_num_nodes()    
    num_elems = mesh.get_num_elems()   
    
    # get the node list coordinates 
    node_list = np.array(mesh.get_nodes()).astype(np.float)    
    # 1 is subtracted from elem_list because node indices start from 0
    elem_list = np.array(mesh.get_elems()).astype(np.int) - 1 
    
    f = open(path+"/out-"+str(l)+".vtu", "w")   
    
    ######## unstructured grids syntax ########
    f.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">\n')
    f.write('\t<UnstructuredGrid>\n')
    f.write('\t\t<Piece NumberOfPoints="'+str(num_nodes)+'"'+' NumberOfCells="'+str(num_elems)+'">\n')      
    
    ######## points ########
    f.write('\t\t\t<Points>\n')
    f.write('\t\t\t\t<DataArray type="Float32" NumberOfComponents="3" Format="ascii">\n')   
    for row in node_list:
        f.write('\t\t\t\t\t'+' '.join(map(lambda x: "{:.3f}\t".format(x), row))+'\n')
    f.write('\t\t\t\t</DataArray>\n')
    f.write('\t\t\t</Points>\n')
    
    ######## cells ########
    f.write('\t\t\t<Cells>\n')
    f.write('\t\t\t\t<DataArray type="Int32" Name="connectivity" Format="ascii">\n')
    for row in elem_list:
        f.write('\t\t\t\t\t'+' '.join(map(lambda x: "{:d}\t".format(x), row[2:]))+'\n')   
    f.write('\t\t\t\t</DataArray>\n')
    f.write('\t\t\t\t<DataArray type="Int32" Name="offsets" Format="ascii">\n')
    cnt = 0
    for row in elem_list:
        cnt += len(row[2:])
        f.write('\t\t\t\t\t'+str(cnt)+'\n')
    f.write('\t\t\t\t</DataArray>\n')
    f.write('\t\t\t\t<DataArray type="Int32" Name="types" Format="ascii">\n')
    for row in elem_list:
        if (row[0]==2): # index of triangular element in gmsh
            f.write('\t\t\t\t\t'+str(5)+'\n')
        if (row[0]==3): # index of quadrilateral element in gmsh
            f.write('\t\t\t\t\t'+str(9)+'\n')
    f.write('\t\t\t\t</DataArray>\n')   
    f.write('\t\t\t</Cells>\n')     
     
    ######## point data ########
    f.write('\t\t\t<PointData Scalars="scalars">\n');
    f.write('\t\t\t\t<DataArray type="Float32" Name="{}" NumberOfComponents="3" Format="ascii">\n'.format('disp'))    
    for i in range(num_nodes):
        if (dim==2):
            f.write('\t\t\t\t\t'+str(U[2*i])+'\t'+str(U[2*i+1])+'\t'+'0.000\n')
        if (dim==3):
            f.write('\t\t\t\t\t'+str(U[3*i])+'\t'+str(U[3*i+1])+'\t'+str(U[3*i+2])+'n')
    f.write('\t\t\t\t</DataArray>\n')  
    f.write('\t\t\t\t<DataArray type="Float32" Name="{}" Format="ascii">\n'.format('d'))    
    for i in range(num_nodes):
        f.write('\t\t\t\t\t'+str(phi[i])+'\n')
    f.write('\t\t\t\t</DataArray>\n') 
    f.write('\t\t\t\t<DataArray type="Float32" Name="{}" NumberOfComponents="3" Format="ascii">\n'.format('Fint'))    
    for i in range(num_nodes):
        if (dim==2):
            f.write('\t\t\t\t\t'+str(F_int[2*i])+'\t'+str(F_int[2*i+1])+'\t'+'0.000\n')
        if (dim==3):
            f.write('\t\t\t\t\t'+str(F_int[3*i])+'\t'+str(F_int[3*i+1])+'\t'+str(F_int[3*i+2])+'n')
    f.write('\t\t\t\t</DataArray>\n')  
    f.write('\t\t\t</PointData>\n')      
    
    ######## close the file ########
    f.write('\t\t</Piece>\n')
    f.write('\t</UnstructuredGrid>\n')
    f.write('</VTKFile>')
    f.close()
