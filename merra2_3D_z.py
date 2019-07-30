#!/usr/bin/env python 
'''
aguyonnet@fas.harvard.edu
read, parse and 2-D plot of MERRA-2 files NetCDF4 format
'''



import os, sys, re
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as pl
import h5py
import pandas as pd
import time as ti
from matplotlib.animation import Animation
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import argparse
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import astropy.time




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



# http://www.ifa.hawaii.edu/88inch/manuals/user.pdf
def geography(site, sealevel):
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

    if sealevel >=0:
        altitude = sealevel
    print 'site, lat_site, lon_site, hPaBin , altitude',\
        site, lat_site, lon_site, hPaBin, altitude
    return area, offset, lat_site, lon_site, lat_city, lon_city, city, hPaBin, altitude


def Date(str):
    dates = str.split('.')[-3]
    year = dates[0:4]
    month = dates[4:6]
    day = dates[6:8]
    return year, month, day


'''
Integrating value down to level in between below and above site bin
this is to determine PWV from QV.
'''
def integIsoALtitude(data, vert_bins, allBp, coefs):
    g       = 9.81 # pesanteur
    lines   = len(data[0,0,:])
    columns = len(data[0,:,0])
    integ0  = np.zeros(shape=(columns, lines), dtype=float)
    integ1  = np.zeros(shape=(columns, lines), dtype=float)
    for line in range(lines):
        for column in range(columns):
            if vert_bins[column, line]>0:
                pwv0 = 0.
                for i in range(int(vert_bins[column, line])+1):
                    pwv0 += data[i, column, line] *  allBp[i, column, line] / g 
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
    print 'method is ', method
    print 'which is at ', lat, lon
    return lat, lon, method


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
                print line, column,  vert_bins[column, line]
                Z[column, line] = data[int(vert_bins[column, line]), column, line]
            else: #if the column does not go down to the required altitude.
                Z[column, line] = np.nan
    print 'The values to interpolate through are : ',  Z
    return Z




'''
I just need to determine the bins for all (lat,lon) coordinates once and for all
--> orography(height)
return array of bin numbers.
'''
def atSite(lon_site, lat_site, X, Y, data, steps, plot,
           otime, xpt_site, ypt_site, site, offset,
           single, allBp, files, orographies, altitude, interpolated):
    '''
    there are three periodicity at which the bin corresponding to the altitude of the site is determined for each (lon,lat):
    from the first timestamp of the first file,
    from the first timestamp of each file,
    for each timespamp.
    
    '''
    if orographies == True:
        modulo = steps/len(files)
        print steps, len(files)
        item= 0
        new_file = np.arange(0,steps,modulo)
        print new_file
        single_orographie = False
    else:
        vert_bins, coefs = getBins(files[0], altitude, 0)
        print 'orography from ', files[0], ' and time ' , 0
        



 
    "determination of the nearest non nan point, both for data and gradient(data)"
    vert_bins, coefs = getBins(files[0], altitude, 0)
    Z   = integValue(data[0], vert_bins, coefs)
    lat, lon, method = InterpOrNearest(lon_site, lat_site, X, Y, Z)
    print ' lat, lon, method ', lat, lon, method

    g_Z  = np.gradient(Z)
    print 'nearest SN gradient'
    gSN_lat, gSN_lon, gSN_method = InterpOrNearest(lon_site, lat_site, X, Y, g_Z[0])
    print 'nearest WE gradient'
    gWE_lat, gWE_lon, gWE_method = InterpOrNearest(lon_site, lat_site, X, Y, g_Z[1])

    
        
    if single != 0:
        steps = [single]
    else:
        steps = np.arange(steps)        
    out    = []

    print " looping on all timestamps "
    print steps
    for tt in steps:
        if orographies :
            nb = tt / modulo
            if single_orographie == True:
                '''
                One orography per file is determined
                '''
                if tt == new_file[item]:
                    if item < (len(files)-1):
                        item+=1
                    table = 0
                    vert_bins, coefs = getBins(files[nb], altitude, table)
                    print 'orography from ', files[nb], ' and time ' , table
            else:
                table = tt % modulo
                vert_bins, coefs = getBins(files[nb], altitude, table)            
                print 'orography from ', files[nb], ' and time ' , table
            
        "Integrating down to altitude"
        if allBp:
            print 'integrating QV'
            Z = integIsoALtitude(data[tt], vert_bins, allBp[tt], coefs)
        else:
            print 'check that PL and integ delp give same result (I think they don t and its depl that is wrong -> does this indicate a problem for QV and AOD?)'
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
        print 'time jd ', time, (astropy.time.Time(time)).jd
        time = (astropy.time.Time(time)).jd
        print len(X.flatten()), len(Y.flatten()), len(Z.flatten())
        val  = interpolate.griddata((X.flatten(), Y.flatten()), Z.flatten(),
                                    (lon, lat),
                                    method = method)

        '''The west-east North_south gradient'''
        gradient = np.gradient(Z)

        WE_gradient = interpolate.griddata((X.flatten(), Y.flatten()),
                                           gradient[1].flatten(),
                                           (gWE_lon, gWE_lat),
                                           method = gWE_method)
        
        WE_gradient = WE_gradient * 10. / 69.4
        print 'West-East gradient on a 10km scale' ,WE_gradient 
        SN_gradient = interpolate.griddata((X.flatten(), Y.flatten()),
                                           gradient[0].flatten(),
                                           (gSN_lon, gSN_lat),
                                           method = gSN_method)
        
        SN_gradient = SN_gradient * 10 / 55.5
        print 'South-North gradient on a 10km scale' ,SN_gradient

        
        if plot :
            fig  = pl.figure(2)
            ax   = fig.gca(projection='3d')
            surf = ax.scatter(X, Y, Z, color='blue')

            figg = pl.figure(3)
            ax2  = figg.gca(projection='3d')
            surf = ax2.scatter(X,Y, gradient[0], color='blue', label='SN')
            surf = ax2.scatter(X,Y, gradient[1], color='red', label='WE')
            pl.legend()
            #gradient x
            fig3,gx = pl.subplots()
            surf = gx.contourf(X, Y, gradient[1], 20)
            cbar = pl.colorbar(surf)
            pl.plot([xpt_site],
                    [ypt_site],
                    label = "W-E grad. over CTIO/10km = "+str(WE_gradient),
                    marker='x',color='r',markersize=12)  # plot a red dot there
            pl.legend()
            pl.title("West-East gradient ")
            #gradient y
            fig4,gy = pl.subplots()
            surf = gy.contourf(X, Y, gradient[0], 100)
            cbar = pl.colorbar(surf)
            pl.title("South-North gradient ")
            pl.plot([xpt_site],
                    [ypt_site],
                    label = "S-N grad. over CTIO/10km = "+str(SN_gradient),
                    marker='x',color='r',markersize=12)  # plot a red dot there
            pl.legend()
            pl.show()
        
        out.append([d, time, val, SN_gradient, WE_gradient])

    out = np.array(out)

    return out



'''
The bins corresponding to the altitude immediately below the site are returned
as well as the coefficient to interpolate with the bin immediately above
The problem may be that the coefs are determined from the scene of the first map
but the coefs  maybe slightly different for the subsequent images -> Mm, 2nd order
'''
def getBins(file, altitude, table):
    nc4f    = h5py.File(file, mode='r')

    height  = nc4f['H'][table,:,:,:]      # time, altitude, lat, lon
    print 'table' , table,  len(nc4f['H'][:,0,0,0])

    #for i in range(72):
    #    print 'height bin ', i,  height[i,:,:]
    lines   = len(height[0,0,:])
    columns = len(height[0,:,0])
    bins    = np.zeros(shape=(columns, lines), dtype=float)
    bins.fill(np.nan)
    coefs   = np.zeros(shape=(columns, lines), dtype=float)
    for line in range(lines):
        for column in range(columns):
            for item, alt in enumerate(height[:,column,line]):
                #print 'bin, line, column, alt, altitude ', item, line, column, alt, altitude
                if alt<altitude:
                    print 'FOUND  ', line, column, alt, altitude
                    bins[column, line] = item # or alt to know the alt, which is needed to interpolate between the two adjacent bins in altitude
                    coef = (altitude - alt) /(height[item-1,column,line] - alt)
                    #print 'coef ', coef
                    coefs[column, line] = coef # to interpolate
                    break
    print 'The bins corresponding to the same layers are : ', bins
    print 'The coef to interpolate at the altitude of the site : ', coefs
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
    parser.add_argument("--orographies", type =bool, 
		        help = "weather or not the orography is redetermined for each file", 
		        default = False)
    parser.add_argument("--interpolated", type =bool, 
		        help = "for a value that is interpolated at site rather than integrated down to site", 
		        default = False)
    parser.add_argument('-m',"--movie", 
		        help = "run movie", 
		        action='store_true')
    parser.add_argument('-v',"--value", type=str, 
	                help = "name of observable", 
	                default='TQV')
    parser.add_argument('-se',"--sealevel", type=float, 
	                help = "Integrate down to sea level", 
	                default=-1.)
    parser.add_argument('-s',"--site", type=str, 
	                help = "name of site", 
	                default='CTIO')
    parser.add_argument('-a',"--hack_altitude", type=float, 
	                help = "input altitude from 0. to 70000.", 
	                default=-1)
    parser.add_argument('-i',"--files", nargs='+', type=str, 
	                help = "input files", 
	                default=None)
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    #if 3 hours slices, if 1hour slices then -> 'H'
    slices = '3H'
    args   = grabargs()
    files         = args.files
    plot          = args.plot
    single        = args.single
    value         = args.value
    sealevel      = args.sealevel
    site          = args.site
    movie         = args.movie
    hack_altitude = args.hack_altitude
    orographies   = args.orographies
    interpolated  = args.interpolated
    
    animeData = []
    allBp     = []
    LAT       = []
    LON       = []
    STEPS     = []
    TIME      = []


    area, offset, lat_site, lon_site,\
        lat_city, lon_city, city, hPaBin, altitude = geography(site, sealevel)
    

    if hack_altitude >= 0.:
        altitude = hack_altitude
        print altitude

    
    fig,ax = pl.subplots()
    frame = Basemap( llcrnrlat=area[0], urcrnrlat=area[1],
                     llcrnrlon=area[2], urcrnrlon=area[3], resolution='i')
    frame.drawcoastlines(color="blue")
    frame.drawcountries(color="black")
    frame.drawparallels(np.arange(area[0], area[1], 2.),
                        labels=[True,False,False,False],
                        linewidth=0.5)
    frame.drawmeridians(np.arange(area[2], area[3], 2.),
                        labels=[True,False,False,True],
                        linewidth=0.5)
    
    xpt_site, ypt_site = frame(lon_site, lat_site)
    xpt_city, ypt_city = frame(lon_city, lat_city)
    frame.plot(xpt_site,
               ypt_site,
               marker='x',color='r',markersize=4)  # plot a red dot there
    
    pl.text(xpt_site + offset,
            ypt_site,
            '%s' % (site),
            #'%s (%5.1fW,%3.1fN)' % (site, lon_site, lat_site),
            color='red', fontsize=12)

    frame.plot(xpt_city,
               ypt_city,
               marker='o',color='r',markersize=3)  # plot a red dot there
    pl.text(xpt_city + offset,
            ypt_city - offset,
            str(city), color='red', fontsize=14)


  
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
        print "longitude,latitude ", longitude,latitude
            
        '''
        retrieving the variable
        '''
        name       = value
        data       = nc4f[name][:,:,:,:]   # time, altitude, lat, lon
        _FillValue = nc4f[name].attrs['_FillValue']
        data[data == _FillValue] = np.nan
        data       = np.ma.masked_where(np.isnan(data), data)


        '''
        When integrating QV, needs also pressure thickness
        '''
        if ((name == 'QV') and (interpolated == False)):            
            for p in ('DELP','delp'):
                print p
                try:
                    deltaP   = nc4f[p][:,:,:,:] #  lat, lon are the same.
                except:
                    print p,' not found'
            allBp.extend(deltaP)


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
    print("after reading data --- %s seconds ---" % ( s_time))
    outname = year+month+day+'_'+value
   

    
    out = atSite(lon_site, lat_site, X, Y, animeData, steps, plot, otime,
                 xpt_site, ypt_site, site, offset,
                 single, allBp, files,
                 orographies, altitude, interpolated)

 
    #for i in out:
    #    print i
    names   = ['entry', 'jd', 'value', 'SN_gradient', 'WE_gradient']
    DumpTuple(names,
                 zip(*out),
                 outname+'.list')
    print 'writing : ', str(outname+'.list')
    print " Note that X,Y from lat,lon are determined only once and for all on the first file, so never solve files of various lat lon extent !!!"
    
