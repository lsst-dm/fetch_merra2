#!/usr/bin/env python 
'''
aguyonnet@fas.harvard.edu
run the extraction of the AOD components
'''

import os, sys, re, glob



def fromRep(rep, date, **kwargs):
    string = "*inst3_3d_aer_Nv."+date+"*"
    files  = glob.glob(os.path.join(rep,string))
    aod_files  = ' '.join(files)
    string = "*inst1_2d_asm_Nx."+date+"*"
    files  = glob.glob(os.path.join(rep,string))
    ozone_files  = ' '.join(files)
    string = "*inst3_3d_asm_Nv."+date+"*"
    files  = glob.glob(os.path.join(rep,string))
    pwv_files  = ' '.join(files)
    return aod_files, ozone_files, pwv_files 
    

if __name__ == "__main__":
    parameters =[
        "DU001",
        "DU002",
        "DU003",
        "DU004",
        "DU005",
        "SS001",
        "SS002",
        "SS003",
        "SS004",
        "SS005",
        "BCPHILIC",
        "BCPHOBIC",
        "OCPHILIC",
        "OCPHOBIC",
        "SO4"]

    site = sys.argv[1]
    rep  = sys.argv[2]
    date = sys.argv[3]
   
    aod_files, ozone_files, pwv_files   = fromRep(rep, date)

    if site =='hawaii':
        site = 'HAWAII'
    if ((site =='ctio') and (site =='lsst')):
        site = 'CTIO'



    cmd = "merra2-chain.py --site %s --value TO3 --files %s"%(site, ozone_files)
    print cmd
    os.system('%s' %cmd)


    cmd = "merra2_3D_z.py --site %s --value QV --files %s"%(site, pwv_files)
    print cmd
    os.system('%s' %cmd)

        
    for p in parameters:
        cmd = "merra2_3D_AOD.py --allaod --site %s --value %s --files %s"%(site, p, aod_files)
        print cmd
        os.system('%s' %cmd)


