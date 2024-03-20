'''

 A generic 3D soil grid wrapper from IMPULSE project
 ******************************
 Authors: G. Louarn, F. Boudon




'''


import numpy as np
import openalea.plantgl.all as pgl
from math import floor

default_properties = ['Qwater', 'Volume_vox', 'Qnorg' , 'Qcorg', 'QNO3', 'QNH4', 'DA', 'Resist' ]

class Soil3D_wrapper(object):
    def __init__(self, origin = (0,0,0), size = (100,100,100), dxyz = (1,1,1), properties = {}, toricities = (True,True,False)):
        """ 
            origin: position of the first voxel (0,0,0) 
            size: number of voxels along z:x:y axes
            dxyz: voxel dimensions along x,y,z
            toricities: define for each dimension if the grid is toric
            properties: dictionnary of 3D numpy grids of the different quantitative properties
        (filled with ones)

            By default, dimensions are expected to be expressed in cm.
            """

        self.maxdimension = 3
        self.toricities = toricities
        self.origin = np.array(origin)
        self.size = np.array(size)
        self.dxyz = np.array(dxyz)
        self.ijk = np.array([1,1,-1])
        #self.gridindexing = pgl.Grid3Indexing( dxyz, origin, self.upper())
        self.m = {}
        for p, v in properties:
            self.add_property(p,v)

    def iindex(self, i, coord):
        val = int(floor((coord - self.origin[i]) / (self.ijk[i] * self.dxyz[i])))
        if self.toricities[i] : val = val % self.size[i]
        return val

    def indexFromPoint(self, pos, unit='m'):
        res = [self.iindex(i, posi) for i, posi in enumerate(pos)]
        return tuple(res)

    #def indexFromPoint(self, pos):
    #    return self.gridindexing.indexFromPoint(pos)

    def upper(self):
        return [self.origin[i] + self.size[i]*self.dxyz[i]*self.ijk[i] if self.ijk[i] > 0 else self.origin[i] for i in range(self.maxdimension)]

    def lower(self):
        return [self.origin[i] + self.size[i]*self.dxyz[i]*self.ijk[i] if self.ijk[i] < 0 else self.origin[i] for i in range(self.maxdimension)]

    def getVoxelCenter(self, index):
        return self.origin + self.ijk * (self.dxyz*index + self.dxyz*0.5)

    def properties(self):
        return self.m

    def add_property(self, name, default_value = 1, type=np.float):
        try :
            assert default_value.shape == self.size
            self.m[name] = default_value
        except:
            self.m[name]= np.ones(self.size, dtype=type) * default_value

    def property(self, name):
        return self.m[name]

    def property_names(self):
        return list(self.m.keys())

    def pgl_representation_property(self, cm='jet', property_name='QWater', sizeratio = 0.1, transparency = 0, minvalue = None, maxvalue = None, scalefunc = None, cmview = False, scaling = 1):
        """ return Plantgl scene """
        if property_name not in self.property_names(): return

        mproperty = self.property(property_name)
        if (not scalefunc is None) and ((not minvalue is None) or (not maxvalue is None)):
            vscalefunc = np.vectorize(scalefunc)
            mproperty = vscalefunc(mproperty)
        minvalue = mproperty.min() if minvalue is None else minvalue
        maxvalue = mproperty.max() if maxvalue is None else maxvalue
        colormap = pgl.PglMaterialMap(minvalue, maxvalue, cm)
        sc = pgl.Scene()
        vsize = np.array(self.dxyz)*sizeratio/2.
        it = np.nditer(self.property(property_name), flags=['multi_index'])
        while not it.finished:                  
            idx, value = it.multi_index, it[0]
            if not scalefunc is None : value = scalefunc(value)
            if minvalue <= value <= maxvalue:
                mat = colormap(value)
                mat.transparency = transparency
                sc += pgl.Shape(pgl.Translated(self.getVoxelCenter(idx)*scaling,
                                               pgl.Box(vsize*scaling)),
                                mat)
            it.iternext()

        if cmview:
            sc += colormap.pglrepr()
        return sc

    def getValueAt(self, property_name, pos):
        cid = self.indexFromPoint(pos)
        return self.property(property_name)[cid]

    def setValueAt(self, property_name, pos, value):
        self.property(property_name)[self.indexFromPoint(pos)] = value

    def incValueAt(self, property_name, pos, value):
        #print pos, self.indexFromPoint(pos)
        self.property(property_name)[self.indexFromPoint(pos)] += value

    def setLayerValue(self, property_name, layerdimension, layervalue, value):
        """
        An operator to assign directly a layer of the soil
        Example :
            # To assign a Z layer with a given value of Water
            soil.setLayerValue('QWater', 2, range(0,3), 5)
        """
        assert 0 <= layerdimension < self.maxdimension
        indices = [slice(None) for i in range(self.maxdimension)]
        indices[layerdimension] = layervalue
        indices = tuple(indices)

        expectedshape = list(self.size)
        try:
            expectedshape[layerdimension] = len(layervalue)
        except:
            expectedshape[layerdimension] = 1

        self.property(property_name)[indices] = value

    def setSliceValue(self, property_name, xslice = slice(None), yslice = slice(None), zslice = slice(None), value = None):
        self.property(property_name)[xslice,yslice,zslice] = value

    #def set_3ds_properties(self, Smodel_obj, ls_properties):
    #    """ set a list of properties for 3ds soil model from  Smodel_obj"""
    #    for p in ls_properties:
    #        vals = getattr(Smodel_obj, p)#comme Smodel_obj.property_name, mais acces par nom
    #        reshaped_vals = np.reshape(vals, self.size)#reshape x,y,z
    #        self.add_property(p, reshaped_vals)

    def set_3ds_properties(self, Smodel_obj, ls_properties):
        """ set a list of properties for 3ds soil model from  Smodel_obj"""
        for p in ls_properties:
            vals = getattr(Smodel_obj, p)  # comme Smodel_obj.property_name, mais acces par nom
            # reshaped_vals = np.reshape(vals, self.size)  # reshape x,y,z
            reshaped_vals = np.zeros(self.size)
            for z in range(self.size[2]):
                reshaped_vals[:, :, z] = vals[z, :, :]

            self.add_property(p, reshaped_vals)


    def __getattr__(self, name):
        try:
            self.__getattribute__(name)
        except:
            try:
                return self.m[name]
            except:
                raise AttributeError(name)

    def __contains__(self,pos):
        return pgl.BoundingBox(self.getLowerCorner(),self.getUpperCorner()-(1e-5,1e-5,1e-5)).contains(pos)


def soil3Dw2s3DSprop(struct1, struct2, propname):
    return np.transpose(struct1.m[propname], (2,1,0) )

def s3DS2soil3Dw(struct1, struct2, propname):
    return struct2.add_property(propname, np.transpose(getattr(struct1,propname), (2,1,0) ))


def pgl_representation(S, m_property, cm='jet', sizeratio=0.1, transparency=0, minvalue=None, maxvalue=None, scalefunc=None, cmview=False, mask=None, dxyz = (1, 1, 1), scaling=1):
    """ return Plantgl scene for 3ds soil grid

    :param S:
    :param m_property:
    :param cm:
    :param sizeratio:
    :param transparency:
    :param minvalue:
    :param maxvalue:
    :param scalefunc:
    :param cmview:
    :param mask:
    :param dxyz:
    :param scaling:
    :return:
    """

    if mask is None:#mask==None:
        mask = S.m_1

    mproperty = m_property
    if (not scalefunc is None) and ((not minvalue is None) or (not maxvalue is None)):
        vscalefunc = np.vectorize(scalefunc)
        mproperty = vscalefunc(mproperty)

    minvalue = mproperty.min() if minvalue is None else minvalue
    maxvalue = mproperty.max() if maxvalue is None else maxvalue
    colormap = pgl.PglMaterialMap(minvalue, maxvalue, cm)
    sc = pgl.Scene()
    vsize = np.array(dxyz) * sizeratio / 2.
    nz, nx, ny = S.m_1.shape[0], S.m_1.shape[1], S.m_1.shape[2]

    for z in range(nz):
        for x in range(nx):
            for y in range(ny):
                value = mproperty[z,x,y]
                idpos = S.get_vox_coordinates(z, x, y)

                if not scalefunc is None: value = scalefunc(value)
                if minvalue <= value <= maxvalue and mask[z,x,y]>=1:
                    mat = colormap(value)
                    mat.transparency = transparency
                    sc += pgl.Shape(pgl.Translated(idpos * scaling, pgl.Box(vsize * scaling)), mat)


    # it = np.nditer(self.property(property_name), flags=['multi_index'])
    # while not it.finished:
    #     idx, value = it.multi_index, it[0]
    #     if not scalefunc is None: value = scalefunc(value)
    #     if minvalue <= value <= maxvalue:
    #         mat = colormap(value)
    #         mat.transparency = transparency
    #         sc += pgl.Shape(pgl.Translated(self.getVoxelCenter(idx) * scaling, pgl.Box(vsize * scaling)), mat)
    #
    #     it.iternext()

    if cmview:
        sc += colormap.pglrepr()

    return sc
