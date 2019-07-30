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



def geography(site):
    area = []
    if (site == 'HAWAII'):
        area.append(15)   # lat
        area.append(23)   # lat
        area.append(-162)   # lon
        area.append(-150)   # lon
        city = 'Honolulu'
        #Site
        lon_site = -155.5 
        lat_site = 19.9
        #City 
        lon_city = -157.8 
        lat_city = 21.35
        offset   = 0.3 # (size / 40)
    if ((site == 'CTIO') or (site == 'LSST')):
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
    if (site == 'OHP'):
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
    return area, offset, lat_site, lon_site, lat_city, lon_city, city


def Date(str):
    dates = str.split('.')[-3]
    year = dates[0:4]
    month = dates[4:6]
    day = dates[6:8]
    return year, month, day


def init():
    foo = ax.contourf(X, Y, start, 100)
    return foo,


def update(num, data, ln, occurences, X, Y, value, s_time, site, lon_site, lat_site):
    data = data[num]
    val  = interpolate.griddata((X.flatten(), Y.flatten()), data.flatten(),(lon_site, lat_site))
    print 'Value at ', site ,' = ', val
    title.set_text(value +' ' +str(occurences[num]))
    minmax.set_text('min ' +str(int(np.min(data)))+ ' max ' +str(int(np.max(data))))
    locVal.set_text('%s %3.2f' % (name, val))
    print 'time : ', occurences[num],\
        'min ' +str((np.min(data)))+ ' max ' +str((np.max(data)))
    print("after 1 --- %s seconds ---" % (ti.time() - s_time))
    ln = ax.contourf(X, Y, data, 100)
    print("after 2--- %s seconds ---" % (ti.time() - s_time))
    return ln,


def atSite(lon_site, lat_site, X, Y, data, steps, plot, otime, xpt_site, ypt_site, site, offset):
    out    = []
    frames = np.arange(steps)
    print frames

    for d in frames:
        Z    = data[d]
        time = str(otime[d])
        print time, (astropy.time.Time(time)).jd
        time = (astropy.time.Time(time)).jd
        val  = interpolate.griddata((X.flatten(), Y.flatten()), Z.flatten(),(lon_site, lat_site))

        '''The west-east North_south gradient'''
        gradient = np.gradient(Z)

        WE_gradient = interpolate.griddata((X.flatten(), Y.flatten()), gradient[1].flatten(),(lon_site, lat_site))
        WE_gradient = WE_gradient * 10. / 69.4
        print 'West-East gradient on a 10km scale' ,WE_gradient 
        SN_gradient = interpolate.griddata((X.flatten(), Y.flatten()), gradient[0].flatten(),(lon_site, lat_site))
        SN_gradient = SN_gradient * 10 / 55.5
        print 'South-North gradient on a 10km scale' ,SN_gradient
        
        if plot:
            fig  = pl.figure(2)
            ax   = fig.gca(projection='3d')
            surf = ax.scatter(X, Y, Z, color='blue')


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


def grabargs():
    usage = "usage: [%prog] [options]\n"
    usage += "fit a continuum and extract EW"
   
    parser = argparse.ArgumentParser(description='Usage',
                                     epilog="extract data from nc4 files")
    parser.add_argument('-p',"--plot", 
		        help = "show control plots", 
		        action='store_true')
    parser.add_argument('-m',"--movie", 
		        help = "run movie", 
		        action='store_true')
    parser.add_argument('-v',"--value", type=str, 
	                help = "name of observable", 
	                default='TO3')
    parser.add_argument('-s',"--site", type=str, 
	                help = "name of site", 
	                default='LSST')
    parser.add_argument('-i',"--files", nargs='+', type=str, 
	                help = "input files", 
	                default=None)
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args   = grabargs()
    # Set up formatting for the movie files
    Writer = animation.writers['ffmpeg']
    files  = args.files
    plot   = args.plot
    value  = args.value
    site  = args.site
    movie  = args.movie
 
    animeData = []
    LAT       = []
    LON       = []
    STEPS     = []
    TIME      = []


    area, offset, lat_site, lon_site, lat_city, lon_city, city = geography(site)
    
 
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
        time = nc4f['/time'][:]
        time_units = nc4f['/time'].attrs['units']

        sampling =  str(24/int(len(time)))+"H"
        print 'sampling is ', sampling
        
        if(value=='TAUTOT' or value=="TOTANGSTR"):
            start_time = re.findall("^minutes since[ ]([0-9.].+[0-9.].+[0-9.].+)[ ]00:30:00$",time_units)
        else:
            start_time = re.findall("^minutes since[ ]([0-9.].+[0-9.].+[0-9.].+)[ ]00:00:00$",time_units)
        
        time_rng   = pd.date_range(start_time[0], periods=time.shape[0], freq=sampling)
        print start_time
        start_time = ti.time()
        steps      = len(time)
        # Hourly map
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
            
 
        # retrieve a variable to plot
        name = value
        data = nc4f[name][:,:,:]
     
        _FillValue = nc4f[name].attrs['_FillValue']
        data[data == _FillValue] = np.nan
        data = np.ma.masked_where(np.isnan(data), data)

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
    
    if movie:
        start = animeData[0]
        img  = frame.contourf(X,Y, start,100)
        title = ax.text(0.2,1.05, "", bbox={'facecolor':'w', 'alpha':0.5, 'pad':5},
                transform=ax.transAxes, ha="center")
        minmax = ax.text(1.,1.05, "", bbox={'facecolor':'w', 'alpha':0.5, 'pad':5},
                transform=ax.transAxes, ha="center")
        locVal = ax.text(xpt_site + 3*offset,
                ypt_site - 2*offset,
                'val %s' % (0.),
                color='red', fontsize=12)


        ani = FuncAnimation(fig, update,
                            frames=np.arange(steps),
                            fargs=(animeData, img, otime,
                                   X, Y, value,
                                   s_time, site, lon_site, lat_site),
                            init_func=init,
                            interval=400,
                            save_count =1,
                            repeat = False)


        cbar = frame.colorbar(img,"right", size="5%", pad="2%")
        cbar.set_label(name + ' ('+unit+')')


        if not plot:           
            print 'writing movie : ', outname
            writer = Writer(fps=5, metadata=dict(artist='Me'), bitrate=1800)
            ani.save(outname+'.mp4', writer=writer)
        else:
            pl.show()

    pl.close()
    out = atSite(lon_site, lat_site, X, Y, animeData, steps, plot, otime,
                 xpt_site, ypt_site, site, offset)
    #for i in out:
    #    print i
    names   = ['entry', 'jd', 'value', 'SN_gradient', 'WE_gradient']
    DumpTuple(names,
              zip(*out),
              outname+'.list')
    print 'writing : ', str(outname+'.list')
     

    divide =  24.
    if value=='AODANA':
        divide = 8.
    pl.plot(out[:,0]/divide, out[:,2])
    pl.xlabel('time (days)')
    pl.ylabel(value)
    pl.show()
    pl.xlabel('time (days)')
    pl.ylabel(value+' '+ site+' gradient on 10km scale')
    pl.plot(out[:,0]/divide, out[:,3], label = 'SN_gradient')
    pl.plot(out[:,0]/divide, out[:,4], label = 'WE_gradient')
    pl.legend()
    pl.show()

    
    
