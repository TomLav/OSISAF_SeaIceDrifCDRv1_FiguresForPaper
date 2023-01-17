import scipy as sp           # scientific Python
import numpy as np
import numpy.ma as ma                             # handling of masked arrays
from numpy import random

def polygon_area(p):
    '''computes the area of a non-intersecting polygon'''
    '''p is a list containing the vertices of the polygon'''
    p.append(p[0]) #close polygon
    area=0.
    for i in sp.arange(len(p)-1):
        area=area+p[i][1]*p[i+1][0]-p[i+1][1]*p[i][0]#Surveyor's Formula
    p.pop()
    return area/2

def grid_area(x,y):
    area=ma.masked_all_like(x)
    for xi in sp.arange(x.shape[1]-1):
        for yi in sp.arange(x.shape[0]-1):
            if not x.mask[yi,xi] and not x.mask[yi,xi+1] and not x.mask[yi+1,xi] and not x.mask[yi+1,xi+1]:
                p=[(y[yi+1,xi],x[yi+1,xi]),(y[yi+1,xi+1],x[yi+1,xi+1]),(y[yi,xi+1],x[yi,xi+1]),(y[yi,xi],x[yi,xi])]
                area[yi,xi]=-polygon_area(p)#minus because points are in clockwise order here
    return area

def calc_divergence(x0,y0,x1,y1):
    '''Compute divergence of a drift field'''
    yi,xi=1,1
    p0=[(y0[yi+1,xi],x0[yi+1,xi]),(y0[yi+1,xi+1],x0[yi+1,xi+1]),(y0[yi,xi+1],x0[yi,xi+1]),(y0[yi,xi],x0[yi,xi])]
    area0=abs(polygon_area(p0))
    area1=abs(grid_area(x1,y1))
    div=ma.masked_all_like(x1)
    div=(area1-area0)/area0
    return div

def calc_diffs(f,grid_space=62.5):

    # Try something smart without loops
    #empty_column = ma.masked_all_like(f[:,0])
    #f_l = ma.column_stack((f[:,1:],empty_column))
    #f_r = ma.column_stack((empty_column,f))[:,1:]
    #empty_row    = ma.masked_all_like(f[0,:])
    #f_b = ma.row_stack((f[1:,:],empty_row))
    #f_u = ma.row_stack((empty_row,f))[1:,:]
    #dfdx = (f_r - f_l)/grid_space
    #dfdy = (f_u - f_b)/grid_space

    # but fall back to loops for now
    dfdx = ma.masked_all_like(f[1:,1:])
    dfdy = ma.masked_all_like(dfdx)
    f_sizey = f.shape[0]
    f_sizex = f.shape[1]
    for row in range(f_sizey-1):
        for col in range(f_sizex-1):
            ul = (row,col)
            ur = (row,col+1)
            bl = (row+1,col)
            br = (row+1,col+1)
            if f.mask[ul] or f.mask[ur] or f.mask[bl] or f.mask[br]:
                continue # already masked
            dfdx.mask[row,col] = False
            dfdy.mask[row,col] = False
            dfdx[row,col] = (((f[ur]-f[ul])+(f[br]-f[bl])))/grid_space
            dfdy[row,col] = (((f[ul]-f[bl])+(f[ur]-f[br])))/grid_space

    return dfdx, dfdy

def calc_deform(dX,dY,grid_space=62.5):
    ''' Compute deformations (divergence, curl, and shear) of a drift field '''

    dudx, dudy = calc_diffs(dX,grid_space)
    dvdx, dvdy = calc_diffs(dY,grid_space)

    div  = dudx+dvdy
    vort = dvdx-dudy
    NDR  = dudx-dvdx
    SDR  = dvdx+dudy
    shr  = (NDR**2 + SDR**2)**0.5

    return div,vort,shr

def calc_deform_random(dX,dY,grid_space=62.5,nb_samples=1000,sdev=2):

    div, vort, shr = calc_deform(dX,dY,grid_space)

    dX_random_shape = list(dX.shape)
    dX_random_shape.extend([nb_samples-1])
    dX_random_shape = tuple(dX_random_shape)
    print(dX_random_shape)

    dX_random = random.normal(0,sdev,dX_random_shape)
    dY_random = random.normal(0,sdev,dX_random_shape)

    dX_local  = dX.copy()
    dY_local  = dY.copy()
    for s in range(nb_samples-1):
        dX_local.data[:]  = dX.data[:] + dX_random[:,:,s]
        dY_local.data[:]  = dY.data[:] + dY_random[:,:,s]
        divR, vortR, shrR = calc_deform(dX_local, dY_local, grid_space)
        div  += divR
        vort += vortR
        shr  += shrR
        if (divmod(s,100)[1] == 0):
            print("Sample %d" % s)

    div  /= nb_samples
    vort /= nb_samples
    shr  /= nb_samples

    return div, vort, shr
