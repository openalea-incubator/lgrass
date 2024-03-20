'''
    VTK
    ***

    Writes VTK files from LightVegeManager geometry and lighting data. Used for visualisation
    We recommend the software Paraview for visualisation
'''
import itertools
import numpy

def VTKtriangles(trimesh, var=[], varname=[], filename=""):
    """Writes VTK files from a triangulation mesh. Possibility to associate physical values to the triangles

    :param trimesh: triangles mesh aggregated by indice elements
        
        .. code-block:: 

            { id : [triangle1, triangle2, ...]}

    :type trimesh: dict of list
    :param var: list of physical values associated to each triangle. Each element of var is a physical value for each triangle.
        Then, for n triangles and m values, var looks like:

        .. code-block:: python

            var = [[var1_1, ..., var1_n], ..., [varm_1, ..., varm_n]]

    :type var: list of list
    :param varname: list of variable names
    :type varname: lits of string
    :param filename: name and path for the output file
    :type filename: string
    """    
    # correct and remove spaces from varname
    varname = [n.replace(" ", "_") for n in varname]

    if isinstance(trimesh, dict) :
        list_triangles = list(itertools.chain(*[v for v in trimesh.values()]))
        ndefaultfiedl = 2
    elif isinstance(trimesh, list) : 
        list_triangles = trimesh
        ndefaultfiedl = 1

    nbtr = len(list_triangles)
    
    f = open(filename, 'w')
    f.write('# vtk DataFile Version 3.0\n')
    f.write('vtk output\n')
    f.write('ASCII\n')
    f.write('DATASET UNSTRUCTURED_GRID\n')
    f.write('POINTS '+str(nbtr * 3)+' float\n')

    for tr in list_triangles:
        for i in range(3):
            f.write(str(tr[i][0])+' '+str(tr[i][1])+' '+str(tr[i][2])+'\n')
    f.write('\n')

    f.write('CELLS '+str(nbtr)+' '+str(nbtr*4)+'\n')
    for i in range(nbtr):
        f.write('3 '+str(3*i)+' '+str(1+3*i)+' '+str(2+3*i)+'\n')
    f.write('\n')

    f.write('CELL_TYPES '+str(nbtr)+'\n')
    for i in range(nbtr):
        f.write('5\n')
    f.write('\n')

    f.write('CELL_DATA '+str(nbtr)+'\n')
    f.write('FIELD FieldData '+str(len(varname) + ndefaultfiedl)+'\n')

    # input fields
    for i, name in enumerate(varname):
        f.write(name+' 1 '+str(nbtr)+' float\n')
        for j in range(nbtr):
            f.write(str(var[i][j])+'\n')
        f.write('\n')
    f.write('\n')

    # ID field : 0 to len(triangles)
    f.write("ID 1 "+str(nbtr)+' float\n')
    for i in range(nbtr):
        f.write(str(i)+'\n')
    f.write('\n')

    # Element ID
    if ndefaultfiedl == 2 :
        f.write("Organ 1 "+str(nbtr)+' float\n')
        for id, triangles in trimesh.items():
            for t in triangles :
                f.write(str(id)+'\n')
        f.write('\n')
    
    f.close()

    var=[]
    varname=[]

def VTKline(start, end, filename):
    """Writes a VTK file representing a line

    :param start: starting point ``[x, y, z]`` of the line
    :type start: list
    :param end: ending point ``[x, y, z]`` of the line
    :type end: list
    :param filename: name and path for the output file
    :type filename: string
    """    
    f=open(filename, 'w')
    f.write('# vtk DataFile Version 3.0\n')
    f.write('vtk Polygon\n')
    f.write('ASCII\n')
    f.write('DATASET POLYDATA\n')
    f.write('POINTS 2 float\n')
    f.write(str(start[0])+" "+str(start[1])+" "+str(start[2])+"\n")
    f.write(str(end[0])+" "+str(end[1])+" "+str(end[2])+"\n")
    f.write("LINES 1 3\n")
    f.write("2 0 1")

    f.close()

def ratp_prepareVTK(ratpgrid, filename, columns=[], df_values=None) :
    """Sets and prepare data to write a voxel grid in VTK format from a RATP grid
    
    Possibility to sets physical values to each voxel. Otherwise, it will assign leaf area density
    for each voxel.

    Then, it runs :func:VTKvoxels to write a VTK file with the grid of voxels.
    
    :param ratpgrid: RAPT grid of voxels
    :type ratpgrid: pyratp.grid
    :param filename: name and path for the output file
    :type filename: string
    :param columns: list of columns name of ``df_values`` to associate to each voxel, defaults to []
    :type columns: list, optional
    :param df_values: dataframe with a ``"Voxel"`` and a ``"VegetationType"`` columns, matching ratpgrid. It stores physical values associated to each voxel that you want to print, defaults to None
    :type df_values: pandas.Dataframe, optional
    """    
    # default prints LAD
    kxyz, entities, lad = [], [], []
    for i in range(ratpgrid.nveg)  :
        for j in range(ratpgrid.nje[i]):
            lad.append(ratpgrid.leafareadensity[j, i])
            entities.append(ratpgrid.nume[j,i])
            kxyz.append(int(i)+1)
    datafields = [numpy.array(kxyz), numpy.array(entities), numpy.array(lad)]

    # other variables
    if columns and df_values is not None:
        for name in columns :
            v = []
            for i in range(ratpgrid.nveg)  :
                for j in range(ratpgrid.nje[i]):
                    filter = (df_values.Voxel == i+1) & (df_values.VegetationType == j+1)
                    if not df_values[filter].empty :
                        v.append(df_values[filter][name].values[0])
                    else:
                        v.append(0.)
                        
            datafields.append(numpy.array(v))
    
    columns.insert(0, "LAD")
    # correct and remove spaces from columns
    columns = [n.replace(" ", "_") for n in columns]

    VTKvoxels(ratpgrid, datafields, columns, filename)

def VTKvoxels(grid, datafields, varnames, nomfich):
    """Writes a VTK file with a grid of voxels and associated physical values.
    
    Display Voxels colored by variable with Paraview
    
    RATP Grid is written in VTK Format as a structured grid

    :param grid: the RATP grid
    :type grid: pyratp.grid
    :param datafields: a list of 3 arrays composed of the a RATP variable to be plotted, corresponding entities, and Voxel ID

        * [0] : kxyz (numpy.array)
        * [1] : entities (numpy.array)
        * [2, ..., n] : variable_0, ... variable_n-2 (numpy.array)

    :type datafields: list of list
    :param varnames:  name of the variable to be plotted
    :type varnames: list of string
    :param nomfich: the VTK filename and path
    :type nomfich: string
    """    
    f=open(nomfich,'w')

    f.write('# vtk DataFile Version 3.0\n')
    f.write('vtk output\n')
    f.write('ASCII\n')
    f.write('DATASET RECTILINEAR_GRID\n')
    f.write('DIMENSIONS '+str(grid.njx+1)+' '+str(grid.njy+1)+' '+str(grid.njz+2)+'\n')

    f.write('Z_COORDINATES '+str(grid.njz+2)+' float\n')
    for i in range(grid.njz):
       z = grid.dz[i:grid.njz].sum()-grid.zorig
       f.write(str(z)+' ')
    f.write(str(-grid.zorig)+' ')
    f.write(str(-10*grid.dz[0]-grid.zorig)+' ')
    f.write('\n')

    f.write('Y_COORDINATES '+str(grid.njy+1)+' float\n')
    for i in range(grid.njy, -1, -1):
       y = grid.dy*i+grid.yorig
       f.write(str(y)+' ')
    f.write('\n')

    f.write('X_COORDINATES '+str(grid.njx+1)+' float\n')
    for i in  range(grid.njx+1):
       x = grid.dx*i+grid.xorig
       f.write(str(x)+' ')
    f.write('\n')

    numVoxels = (grid.njx)*(grid.njy)*(grid.njz+1)
    f.write('CELL_DATA '+str(numVoxels)+'\n')

    # Set the number of entities to write - NbScalars
    numberofentities = len(set(datafields[1]))
    numberofdatafields = len(datafields) - 2

    for ent in range(numberofentities) :
        for nvar in range(numberofdatafields) : 
            f.write('SCALARS ' + varnames[nvar] + '_entity_' + \
                                            str(int(ent)) + ' float  1 \n')
            f.write('LOOKUP_TABLE default\n')
            for ik in range(grid.njz+1):
                for ij in range(grid.njy):
                    for ii in range(grid.njx): 
                        k =grid.kxyz[ii,ij,ik]
                        if (k>0): 
                            kidDummy = numpy.where(numpy.array(datafields[0])==k)
                            kid = kidDummy[0]
                            entity = numpy.array(datafields[1])[kid]
                            ent_in_vox =numpy.where(entity == ent+1)
                            # pas dans le voxel
                            if numpy.alen(ent_in_vox[0])<1:
                                f.write(str(-9999.0)+'\n')
                            else:
                                # couche du sol
                                if ik == grid.njz:
                                    f.write(str(-9999.0)+'\n')
                                else:
                                    where = (kid[numpy.where(entity==ent+1)])
                                    value = datafields[2+nvar][where[0]]
                                    f.write(str(value)+'\n')
                        # voxel vide
                        else:
                            f.write(str(-9999.0)+'\n')
            f.write('\n')

        # indice des voxels
        f.write('SCALARS Voxel_ID float  1 \n')
        f.write('LOOKUP_TABLE default\n')
        for ik in range(grid.njz+1):
            for ij in range(grid.njy):
                for ii in range(grid.njx):
                    k =grid.kxyz[ii,ij,ik]
                    f.write(str(k)+'\n')
        f.write('\n')
    f.close()


def PlantGL_to_VTK(scenes, path, i=0, in_unit="m", out_unit="m"):
    """Directly converts a list plantGL scenes to a single VTK file
    The routine aggregates all the scenes in one global triangulation and then calls :func:VTKtriangles
    Possbility to rescale measure unit scenes

    :param scenes: list of scenes to be written
    :type scenes: lits of plantgl.Scene
    :param path: path of the file. File name is set by default with ``"triangles_plantgl_"+str(i)+".vtk"``
    :type path: string
    :param i: id of the file. Used if you want to batch a large number of scenes, defaults to 0
    :type i: int, optional
    :param in_unit: input measure unit, defaults to "m"
    :type in_unit: str, optional
    :param out_unit: output measure unit, defaults to "m"
    :type out_unit: str, optional
    :raises ValueError: in_unit or out_unit not a valid entry upon the listed measure units
    """  
    from alinea.caribu import plantgl_adaptor
    
    # rescale the scenes 
    units = {'mm': 0.001, 
                'cm': 0.01, 
                'dm': 0.1, 
                'm': 1, 
                'dam': 10, 
                'hm': 100,
                'km': 1000}
    if in_unit not in units or out_unit not in units:
        raise ValueError("Unknown unit: select one in this \
                    list ['mm','cm','dm','m','dam','hm','km']")
    
    rescale=False
    if (in_unit != out_unit) : rescale=True
    
    trimesh = {}
    # aggregates the scenes
    if type(scenes) == list :
        for s in scenes :   
            trimesh = trimesh + plantgl_adaptor.scene_to_cscene(s)
    else :
       trimesh = plantgl_adaptor.scene_to_cscene(scenes)

    if rescale :
        h = units[in_unit]/units[out_unit]
        for id in trimesh.keys() : trimesh[id] = rescale(trimesh[id], h)

    VTKtriangles(trimesh, [], [], path+"triangles_plantgl_"+str(i)+".vtk")