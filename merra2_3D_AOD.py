#!/usr/bin/env python 
'''
aguyonnet@fas.harvard.edu
read, parse and extract parameters above given sites from MERRA-2 files NetCDF4 format
'''

from mpl_toolkits.basemap import Basemap
import os, sys, re
import numpy as np
import matplotlib.pyplot as pl
import h5py
import pandas as pd
import time as ti
import argparse
from scipy import interpolate
from matplotlib import cm
import astropy.time
from mpl_toolkits.mplot3d import Axes3D

g       = 9.81 # pesanteur


def DumpTuple(names, list, file):
    out  = open(file, 'w')
    for i in names :
        out.write('# '+i + ' :' +'\n')
    out.write('#end'+ '\n')
    list = zip(*list)
    for i in list:
        out.write(' '.join(map("{}".format, i))+'\n')
    out.close()
    return




def Beta(rh, type):
    unit  = 1000. 
    scale = [0.0, 0.8, 0.95]
    if type == 'DU001':
        beta_table = [2.02, 2.02, 2.02]
    if type == 'DU002':
        beta_table = [0.64, 0.64, 0.64]
    if type == 'DU003':
        beta_table = [0.33, 0.33, 0.33]
    if type == 'DU004':
        beta_table = [0.17, 0.17, 0.17]
    if type == 'DU005':
        beta_table = [0.08, 0.08, 0.08]
    if type == 'SS001':
        beta_table = [0.73, 4.54, 25.98]
    if type == 'SS002':
        beta_table = [3.48, 10.01, 24.02]
    if type == 'SS003':
        beta_table = [0.74, 2.04, 4.85]
    if type == 'SS004':
        beta_table = [0.30, 0.86, 2.02]
    if type == 'SS005':
        beta_table = [0.10, 0.30, 0.72]
    if type == 'BCPHOBIC':
        beta_table = [9.28, 9.28, 9.28]
    if type == 'BCPHILIC':
        beta_table = [9.28, 11.27, 15.77]
    if type == 'OCPHOBIC':
        beta_table = [2.67, 2.67, 2.67]
    if type == 'OCPHILIC':
        beta_table = [2.67, 7.01, 16.04]
    if type == 'SO4':
        beta_table = [3.15, 14.29, 22.53]
    return np.interp(rh, scale, beta_table) * unit     


# http://www.ifa.hawaii.edu/88inch/manuals/user.pdf
def geography(site):
    area = []
    if (site == 'HAWAII'):
        hPaBin = 53
        altitude = 4214. 
        area.append(15)   # lat
        area.append(23)   # lat
        area.append(-162)   # lon
        area.append(-150)   # lon
        city = 'Honolulu'
        #Site
        lon_site = -155.472 
        lat_site = 19.826 # Exact location UH teleso
        #City 
        lon_city = -157.8 
        lat_city = 21.35
        offset   = 0.3 # (size / 40)
    if (site == 'CTIO'):
        altitude = 2200.
        hPaBin = 59
        area.append(-40)   # lat
        area.append(-20)   # lat
        area.append(-77)   # lon
        area.append(-60)   # lon
        city = 'Santiago'
        #CTIO Site
        lon_site = -70.815 
        lat_site = -30.16528
        #Santiago 
        lon_city = -70.674 
        lat_city = -33.447
        offset   = 0.5 # (size / 40)
    if (site == 'ANTOFAGASTA'):
        altitude =0.
        hPaBin = 72
        area.append(-30)   # lat
        area.append(-10)   # lat
        area.append(-77)   # lon
        area.append(-60)   # lon
        city = 'Antofagasta'
        #CTIO Site
        lon_site = -70.441 
        lat_site = -23.450
        #Santiago 
        lon_city = -70.441
        lat_city = -23.450
        offset   = 0.5 # (size / 40)
    if (site == 'OHP'):
        altitude = 600.
        hPaBin = 70
        # this is stupid, I should use input data and enlarge them
        area.append(40)   # lat
        area.append(49)   # lat
        area.append(1)   # lon
        area.append(10)   # lon
        city = 'Marseille'
        #OHP Site
        lon_site = 5.718683
        lat_site  = 43.923179
        #Marseille
        lon_city = 5.408641
        lat_city = 43.342758
        offset   = 0.25
    print 'site, lat_site, lon_site, hPaBin ',\
        site, lat_site, lon_site, hPaBin
    return area, offset, lat_site, lon_site, lat_city, lon_city, city, hPaBin, altitude


def Date(str):
    dates = str.split('.')[-3]
    year = dates[0:4]
    month = dates[4:6]
    day = dates[6:8]
    return year, month, day


def closest(X,Y, x_site, y_site, Z):
    ref = 2.
    for enum_i, i in enumerate(X):
        for enum_j, j in enumerate(Y):
            dist = np.sqrt((i-x_site)**2 + (j-y_site)**2)
            if ((dist<ref) and (str(Z[enum_j, enum_i])!='nan')):
                print enum_j, enum_i, i,j, Z[enum_j, enum_i]
                ref   = dist
                closest_x = enum_i
                closest_y = enum_j
                lon = i
                lat = j
    return ref, lon, lat, closest_x, closest_y


'''
If nan surrounds the site, value is not interpolated. 
Non Nearest nan neighbor is returned instead.
'''
def InterpOrNearest(lon_site, lat_site, X, Y, Z):
    val = interpolate.griddata((X.flatten(), Y.flatten()), Z.flatten(),
                               (lon_site, lat_site),
                               method='linear')

    if(str(val)=='nan'):        
        ref, lon, lat, closest_x, closest_y = closest(X[0],Y[:,0], lon_site, lat_site, Z)
        #print "the Z ", Z[closest_y, closest_x]
        method = 'nearest'       
    else:
        lat = lat_site
        lon = lon_site
        method = 'linear'
    print 'method is ', method, ' at lat, lon = ', lat, lon
    return lat, lon, method



'''
Integrating value down to level in between below and above site bin
'''
def integValue(data, vert_bins, coefs):
    lines   = len(data[0,0,:])
    columns = len(data[0,:,0])
    integ0  = np.zeros(shape=(columns, lines), dtype=float)
    integ1  = np.zeros(shape=(columns, lines), dtype=float)
    for line in range(lines):
        for column in range(columns):
            if vert_bins[column, line]>0:
                pwv0 = 0.
                for i in range(int(vert_bins[column, line])+1):
                    pwv0 += data[i, column, line] 
                    if i == int(vert_bins[column, line]-1): #last bin above, to interp at z
                        pwv1 = pwv0
                integ0[column, line] = pwv0
                integ1[column, line] = pwv1
            else: #if the column does not go down to the required altitude.
                integ0[column, line] = np.nan
                integ1[column, line] = np.nan
    #print 'The integrated values (below z) to interpolate through are : ',  integ0
    #print
    #print 'The high integrated values (above z) to interpolate through are : ',  integ1
    delta = integ1 - integ0    # interpolating between two adjacent z bins
    integ = integ0 + coefs*delta  # check ok
    return integ


'''
Integrating value down to level in between below and above site bin
allRh is used to determined the beta value in each bin
'''
def integIsoAOD(data, vert_bins, allBp, allRh, coefs, variable):
    lines   = len(data[0,0,:])
    columns = len(data[0,:,0])
    integ0  = np.zeros(shape=(columns, lines), dtype=float)
    integ1  = np.zeros(shape=(columns, lines), dtype=float)    
    for line in range(lines):
        for column in range(columns):
            if vert_bins[column, line]>0:
                value = 0.
                for i in range(int(vert_bins[column, line])+1):
                    beta =  Beta(allRh[i, column, line], variable)
                    value += data[i, column, line] *  allBp[i, column, line] * beta / g 
                    if i == int(vert_bins[column, line]-1): #last bin above, to interp at z
                        value_above = value
                integ0[column, line] = value
                integ1[column, line] = value_above
            else: #if the column does not go down to the required altitude.
                integ0[column, line] = np.nan
                integ1[column, line] = np.nan
    #print 'The integrated values (below z) to interpolate through are : ',  integ0
    #print
    #print 'The high integrated values (above z) to interpolate through are : ',  integ1
    delta = integ1 - integ0    # interpolating between two adjacent z bins
    integ = integ0 + coefs*delta  # check ok
    return integ




'''
return the ndarray of the data at the same altitude
'''
def getIsoAOD(data, vert_bins, allBp, allRh, variable):
    lines   = len(data[0,0,:])
    columns = len(data[0,:,0])
    Z       = np.zeros(shape=(columns, lines), dtype=float)
    for line in range(lines):
        for column in range(columns):
            if vert_bins[column, line]>0:
                beta =  Beta(allRh[int(vert_bins[column, line]), column, line], variable)
                Z[column, line] = data[int(vert_bins[column, line]), column, line] *  allBp[int(vert_bins[column, line]), column, line] * beta / g                 
            else: #if the column does not go down to the required altitude.
                Z[column, line] = np.nan
    print 'getIsoAOD: The values to interpolate through are : ',  Z
    return Z


'''
return the ndarray of the data at the same altitude
'''
def getIsoAltitude(data, vert_bins):
    lines   = len(data[0,0,:])
    columns = len(data[0,:,0])
    Z       = np.zeros(shape=(columns, lines), dtype=float)
    for line in range(lines):
        for column in range(columns):
            if vert_bins[column, line]>0:
                #print line, column,  vert_bins[column, line]
                Z[column, line] = data[int(vert_bins[column, line]), column, line]
            else: #if the column does not go down to the required altitude.
                Z[column, line] = np.nan
    print 'getIsoAltitude: The values to interpolate through are : ',  Z
    return Z




'''
I just need to determine the bins for all (lat,lon) coordinates once and for all
--> orography(height)
return array of bin numbers.
'''
def atSite(lon_site, lat_site, X, Y, data, steps, plot,
           otime, xpt_site, ypt_site, site, offset,
           vert_bins, single, coefs, allBp, allRh, variable,
           allaod, interpolated, species):
    if single != 0:
        steps = [single]
    else:
        steps = np.arange(steps)        

 
    "determination of the nearest non nan point, both for data and gradient(data)"
    #Z   = integIsoAOD(data[0], vert_bins, allBp[0], allRh[0], coefs, variable)
    Z   = integValue(data[0], vert_bins, coefs)

    print 'nearest value'
    lat, lon, method = InterpOrNearest(lon_site, lat_site, X, Y, Z)
    g_Z  = np.gradient(Z)
    print 'nearest SN gradient'
    gSN_lat, gSN_lon, gSN_method = InterpOrNearest(lon_site, lat_site, X, Y, g_Z[0])
    print 'nearest WE gradient'
    gWE_lat, gWE_lon, gWE_method = InterpOrNearest(lon_site, lat_site, X, Y, g_Z[1])


    " looping on all timestamps "
    out = []
    for tt in steps:
        if allaod is True :        
            Z = integIsoAOD(data[tt], vert_bins,
                                 allBp[tt], allRh[tt],
                                 coefs, variable)
        else:
            print 'check that PL and integ delp give same result (I think they don t and its depl that is wrong -> does this indicate a problem for QV and AOD?)'          
            if species is True: # for aod species at a given altitude
                print variable, " aod species at a given altitude"
                Z0    = getIsoAOD(data[tt], vert_bins, allBp[tt], allRh[tt], variable)
                Z1    = getIsoAOD(data[tt], vert_bins-1, allBp[tt], allRh[tt], variable)
                delta = Z1 - Z0    # interpolating between two adjacent z bins
                Z     = Z0 + coefs*delta  # check ok
            else:
                if interpolated :  # for PL
                    Z0    = getIsoAltitude(data[tt], vert_bins)
                    Z1    = getIsoAltitude(data[tt], vert_bins-1)
                    delta = Z1 - Z0    # interpolating between two adjacent z bins
                    Z     = Z0 + coefs*delta  # check ok
                else: 
                    Z = integValue(data[tt], vert_bins, coefs) 


        if (np.shape(data[tt][0]) != np.shape(Z)):
            print "SHAPES : ", np.shape(data[tt][0]),' must be equal ', np.shape(Z)
            sys.exit()

        time = str(otime[tt])
        d    = vert_bins[len(X)/2, len(Y)/2] # returning central bin number
        print 'time = ', time, ' JD = ', (astropy.time.Time(time)).jd
        time = (astropy.time.Time(time)).jd

        
        val  = interpolate.griddata((X.flatten(), Y.flatten()),
                                    Z.flatten(),
                                    (lon, lat),
                                    method = method)
        print 'val at lon lat ', lon, lat, ' is ', val
        if plot :
            fig  = pl.figure(2)
            ax   = fig.gca(projection='3d')    
            surf = ax.scatter(X, Y, Z, color='blue')
            ax.scatter(lon, lat, 0., marker='x',color='r') 
        
        '''The west-east North_south gradient'''
        gradient    = np.gradient(Z)
        WE_gradient = interpolate.griddata((X.flatten(), Y.flatten()),
                                           gradient[1].flatten(),
                                           (gWE_lon, gWE_lat),
                                           method = gWE_method)

        
      
        WE_gradient = WE_gradient * 10. / 69.4
        print 'West-East gradient on a 10km scale' ,WE_gradient
        pl.show()
        
        SN_gradient = interpolate.griddata((X.flatten(), Y.flatten()),
                                           gradient[0].flatten(),
                                           (gSN_lon, gSN_lat),
                                           method = gSN_method)
        SN_gradient = SN_gradient * 10 / 55.5
        print 'South-North gradient on a 10km scale' ,SN_gradient
        out.append([d, time, val, SN_gradient, WE_gradient])
    out = np.array(out)
    return out



'''
The bins corresponding to the altitude immediately below the site are returned
as well as the coefficient to interpolate with the bin immediately above
The problem may be that the coefs are determined from the scene of the first map
but the coefs  maybe slightly different for the subsequent images -> Mm, 2nd order

When table has no altitude entry, use a grid form file
'''

def loadarray(cat):
    array = []
    file  = open(cat, "r")
    lines = file.readlines()
    for line in lines :
        array.append([float(x) for x in line.split()])
    return array

          


'''
The aer files have no H key. The bins and coefs are found from mean Barometric pressure at site
'''
def getBinsfromBp(file, pressure, table=0, **kwargs):
    nc4f    = h5py.File(file, mode='r')
    bp      = nc4f['DELP'][table,:,:,:]      # time, altitude, lat, lon
    lines   = len(bp[0,0,:])
    columns = len(bp[0,:,0])
    bins    = np.zeros(shape=(columns, lines), dtype=float)
    bins.fill(np.nan)
    coefs   = np.zeros(shape=(columns, lines), dtype=float)
    for line in range(lines):
        for column in range(columns):
            press = 0
            for item, alt in enumerate(bp[:,column,line]):
                press += alt
                if press > pressure:
                    print 'line, column, pressure ', line, column, press
                    bins[column, line] = item # or press to know the pressure, which is needed to interpolate between the two adjacent bins in pressure
                    sumabove = 0
                    for i in range(item):                        
                        sumabove += bp[i,column,line]
                    coef = (pressure - press) /(sumabove - press)
                    coefs[column, line] = coef # to interpolate                
                    break
    print 'The bins corresponding to the same layers are : ', bins
    print 'The coef to interpolate at the pressure of the site : ', coefs
    return bins, coefs





def grabargs():
    usage = "usage: [%prog] [options]\n"
    usage += "fit a continuum and extract EW"
   
    parser = argparse.ArgumentParser(description='Usage',
                                     epilog="extract data from nc4 files")
    parser.add_argument('-p',"--plot", 
		        help = "show control plots", 
		        action='store_true')
    parser.add_argument('-g',"--single", type =int, 
		        help = "only first time stamp", 
		        default= 0)
    parser.add_argument('-m',"--movie", 
		        help = "run movie", 
		        action='store_true')
    parser.add_argument('-v',"--value", type=str, 
	                help = "name of observable", 
	                default='TQV')
    parser.add_argument('-s',"--site", type=str, 
	                help = "name of site", 
	                default='CTIO')
    parser.add_argument('-a',"--altitude", type=int, 
	                help = "altitude from 0 to 71", 
	                default=71)
    parser.add_argument('-i',"--files", nargs='+', type=str, 
	                help = "input files", 
	                default=None)
    parser.add_argument("--interpolated", type =bool, 
		        help = "for a value that is interpolated at site rather than integrated down to site", 
		        default = False)
    parser.add_argument("--allaod", action='store_true', 
		        help = "set to false to measure something else than allAOD")
    parser.add_argument("--species", action='store_true',  
		        help = "to measure a given species at a given altitude")
    parser.add_argument("--hack_pressure", type=float, 
	                help = "input pressure from 0. to 100000.", 
	                default=-1)
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    #if 3 hours slices, if 1hour slices then -> 'H'
    slices = '3H'
    args   = grabargs()
    # Set up formatting for the movie files
    files  = args.files
    plot   = args.plot
    single = args.single
    value  = args.value
    site   = args.site
    movie  = args.movie
    altitude      = args.altitude
    species       = args.species
    allaod        = args.allaod
    interpolated  = args.interpolated
    hack_pressure = args.hack_pressure

    
    animeData = []
    allBp     = [] # barometric pressure
    allRh     = [] # relative humidity
    LAT       = []
    LON       = []
    STEPS     = []
    TIME      = []


    area, offset, lat_site, lon_site,\
        lat_city, lon_city, city, hPaBin, altitude = geography(site)
    frame = Basemap( llcrnrlat=area[0], urcrnrlat=area[1],
                     llcrnrlon=area[2], urcrnrlon=area[3], resolution='i')
    xpt_site, ypt_site = frame(lon_site, lat_site)
    xpt_city, ypt_city = frame(lon_city, lat_city)
  
    '''
    the bin corresponding to the altitude of the site is determined from a table
    '''
    if site.lower() == 'hawaii':
        pressure = 61600. # needed in pascal 
    elif site.lower() == 'ctio':
        pressure = 77600. # needed in pascal
    elif site.lower() == 'ohp':
        pressure = 94200. # needed in pascal




    if hack_pressure >= 0.:
        pressure = hack_pressure
        print "hacked pressure = ", pressure

    
        
    vert_bins, coefs = getBinsfromBp(files[0], pressure)
    
     
    for item, file in enumerate(files):
        year, month, day = Date(file)
        print year, month, day

        
        nc4f = h5py.File(file, mode='r')
        unit = nc4f[value].attrs['units']
        print nc4f.items()

        
        # Get the time data.
        time       = nc4f['/time'][:]
        time_units = nc4f['/time'].attrs['units']
        start_time = re.findall(r'\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}', time_units)
        time_rng   = pd.date_range(start_time[0], periods=time.shape[0], freq=slices)
        start_time = ti.time()

        
        if single !=0 :
            steps      = single
        else:
            steps      = len(time)

        # timestamps
        for i in range(steps):
            print i, time_rng[i]       


        # Get the geolocation data.
        latitude  = nc4f['/lat'][:]
        longitude = nc4f['/lon'][:]
        lat_start = min(latitude)
        lon_start = min(longitude)
        lat_end   = max(latitude)
        lon_end   = max(longitude)
        print "longitude = ", longitude
        print "latitude  = ", latitude
            
        '''
        retrieving the variable
        '''
        name       = value
        data       = nc4f[name][:,:,:,:]   # time, altitude, lat, lon
        _FillValue = nc4f[name].attrs['_FillValue']
        data[data == _FillValue] = np.nan
        data       = np.ma.masked_where(np.isnan(data), data)


        '''
        needs pressure thickness and relative humidity
        '''

        for p in ('DELP','delp'):
            print p
            try:
                deltaP   = nc4f[p][:,:,:,:] #  lat, lon are the same.
            except:
                print p,' not found'
        allBp.extend(deltaP)


        allRh.extend(nc4f['RH'][:,:,:,:])
        
        '''
        Recording several days
        '''
        if item ==0:
            LAT = latitude
            LON = longitude
        

        STEPS.append(steps)
        animeData.extend(data)
        TIME.append(time_rng)
     

    X,Y = np.meshgrid(LON,LAT)
    steps = sum(STEPS)
    otime = []
    for t in TIME:
        for i in t:
            otime.append(i)


    s_time = ti.time()
    outname = year+month+day+'_'+value
   

    out = atSite(lon_site, lat_site, X, Y, animeData, steps, plot, otime,
                 xpt_site, ypt_site, site, offset,
                 vert_bins, single, coefs,
                 allBp, allRh, value, allaod, interpolated, species)

        
    names   = ['entry', 'jd', 'value', 'SN_gradient', 'WE_gradient']
    tb.DumpTuple(names,
                 zip(*out),
                 outname+'.list')
    print 'writing : ', str(outname+'.list')

    
