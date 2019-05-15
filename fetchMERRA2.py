#!/usr/bin/env python 
'''
aguyonnet@fas.harvard.edu
wget from NASA MERRA-2 online tables 
inst3_3d_asm_Nv (M2I3NVASM): Assimilated Meteorological Fields 
given parameter [pwv,] [date] [site=(lsst,...)] 
list of tables :
https://gmao.gsfc.nasa.gov/pubs/docs/Bosilovich785.pdf
data can be access from :
https://disc.sci.gsfc.nasa.gov/datasets?page=1&keywords=MERRA-2
'''

import os, sys, re
import numpy as np
import argparse
import datetime


'''
latest available data are from previous month

'''

'''
Template for requests : 
http://goldsmr4.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2FMERRA2%2FM2I1NXASM.5.12.4%2F2018%2F01%2FMERRA2_400.inst1_2d_asm_Nx.20180101.nc4&FORMAT=bmM0Lw&BBOX=-32.587%2C-73.213%2C-28.017%2C-68.643&LABEL=MERRA2_400.inst1_2d_asm_Nx.20180101.SUB.nc&SHORTNAME=M2I1NXASM&SERVICE=SUBSET_MERRA2&VERSION=1.02&DATASET_VERSION=5.12.4

https://goldsmr4.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2FMERRA2%2FM2T1NXAER.5.12.4%2F2012%2F02%2FMERRA2_400.tavg1_2d_aer_Nx.20120201.nc4&FORMAT=bmM0Lw&BBOX=-33.554%2C-72.686%2C-28.808%2C-67.939&LABEL=MERRA2_400.tavg1_2d_aer_Nx.20120201.SUB.nc&SHORTNAME=M2T1NXAER&SERVICE=SUBSET_MERRA2&VERSION=1.02&DATASET_VERSION=5.12.4

http://goldsmr4.gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2FMERRA2%2FM2I1NXASM.5.12.4%2F2018%2F01%2FMERRA2_400.inst1_2d_asm_Nx.20180101.nc4&FORMAT=bmM0Lw&BBOX=-32.587%2C-73.213%2C-28.017%2C-68.643&LABEL=MERRA2_400.inst1_2d_asm_Nx.20180101.SUB.nc&SHORTNAME=M2I1NXASM&SERVICE=SUBSET_MERRA2&VERSION=1.02&DATASET_VERSION=5.12.4
'''


'''
create, name and edit the file submitting the request to the merra-2 server
for LSST site, ozone table, at the day == date
'''


def buildRequest(request_directory, site, date, parameter):
    if parameter == 'pwv':
        parameter = 'I3NVASM'
    if parameter == 'ozone':
        parameter = 'I1NXASM'
    if parameter == 'aod':
        parameter = 'I3NVGAS'
    if parameter == 'mmr':
        parameter = 'I3NVAER'
    if parameter == 'angstrom':
        parameter = 'T1NXAER'

    filename = "subset_M2%s_V5.12.4_%s.txt"\
               % (parameter, datetime.date.today().strftime('%Y%m%d'))
    request  = os.path.join(request_directory, filename)
    writeFile(request, date, site, parameter)
    return request, filename


def writeFile(request, date, site="hawaii", parameter="pwv", **kwargs):
    if (site == 'hawaii'):
        box = "16.402%2C-159.214%2C22.862%2C-152.139"
    if (site == 'ctio' or site == 'lsst'):
        box= "-32.499%2C-74.399%2C-26.698%2C-67.983"
    time    = '%04d%02d%02d'%(date.year, date.month, date.day)
    if int(time) >= 20110101:
        merra_xxx = "400"
    else:
        merra_xxx = "300"
        
    if parameter == 'I3NVASM':
        fname = 'inst3_3d_asm_Nv'
        version = "5"
    if parameter == 'I3NVGAS':
        fname = 'inst3_3d_gas_Nv'
        version = "5"
    if parameter == 'I3NVAER':
        fname = 'inst3_3d_aer_Nv'
        version = "5"
    if parameter == 'I1NXASM':
        version = "4"
        fname ='inst1_2d_asm_Nx'
    if parameter == 'T1NXAER':
        fname = 'tavg1_2d_aer_Nx'
        version = "4"
    f = open(request,'w')

 
    command = "http://goldsmr%s"%version
    command += ".gesdisc.eosdis.nasa.gov/daac-bin/OTF/HTTP_services.cgi?FILENAME=%2Fdata%2FMERRA2%2FM2"
    command += "%s" % parameter
    command += ".5.12.4%2F"
    command += "%s"   %(date.year)
    command += "%2F"
    command += "%02d" %(date.month)
    command += "%2F" 
    command += "MERRA2_%s." % merra_xxx
    command += "%s"         % fname
    command += ".%s"        % time
    command += ".nc4&FORMAT=bmM0Lw&BBOX=%s&LABEL="%(box)
    command += "MERRA2_%s." % merra_xxx
    command += "%s"         % fname
    command += ".%s"        % time
    command += ".SUB.nc&SHORTNAME=M2"
    command += "%s" % parameter
    command += "&SERVICE=SUBSET_MERRA2&VERSION=1.02&DATASET_VERSION=5.12.4https://goldsmr%s"%version
    if not parameter == 'T1NXAER':
        command += ".gesdisc.eosdis.nasa.gov/data/MERRA2/M2"
        command += "%s" % parameter
        command += ".5.12.4/doc/MERRA2.README.pdf"
    f.write(command + '\n')
    f.close()
    return
    

def grabargs():
    usage = "usage: [%prog] [options]\n"
    usage += "wget from NASA MERRA-2 online tables given parameter [ozone,] [date] [site=(lsst,...)] "
   
    parser = argparse.ArgumentParser(description='Usage',
                                     epilog="wget merra-2 data")
    parser.add_argument('-p',"--parameters", type=str,
		        help = "indicate a parameter : pwv, aod, mmr", 
		        default= ['ozone', 'aod', 'pwv'], nargs='+')
    parser.add_argument('-d', '--dates',  # either of this switches
                        help = "indicate date yyyy-mm-dd or yesterday", 
                        type = str,
                        default ='latest',
                        dest = 'dates',     
                        nargs='+')
    parser.add_argument('-s',"--site", type=str, 
	                help = "indicate a site", 
	                default='lsst')
    parser.add_argument('-r',"--dir", type=str, 
	                help = "downloading in directory", 
	                default='./')
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args              = grabargs()
    parameters        = args.parameters
    dates             = args.dates
    site              = args.site
    dir               = args.dir

    '''
    latest available date is last day of last month
    '''
    if dates == 'latest': 
        today = datetime.date.today()
        last_day_of_last_month = datetime.datetime(today.year, today.month, 1) - datetime.timedelta(days=1)
        dates = [last_day_of_last_month.strftime("%Y-%m-%d")]
    
    print dates
        
    for day in dates:
        print 'requested date : ', day
        for parameter in parameters:
            print "requested parameter : ", parameter
            date = datetime.datetime.strptime(day, '%Y-%m-%d')
            request, filename  = buildRequest(dir, site, date, parameter)
            wget    = "wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --auth-no-challenge=on --keep-session-cookies --content-disposition -i %s  --directory-prefix=%s" % (request, dir)
            print wget
            os.system('%s' %(wget))
            rm = "rm -f %s"%filename
            print rm
            os.system('%s' %(rm))
