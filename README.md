# fetch_merra2


Accessing and processing MERRA-2 data


      -   1 .   -
      
NASA provides an online map of worldwide atmospheric parameters :
https://worldview.earthdata.nasa.gov/?p=geographic&l=VIIRS_SNPP_CorrectedReflectance_TrueColor(hidden),MODIS_Aqua_CorrectedReflectance_TrueColor(hidden),MODIS_Terra_CorrectedReflectance_TrueColor,MODIS_Terra_Aerosol,Reference_Labels,Reference_Features,Coastlines&t=2017-05-02&z=3&v=-125.4375,-48.8671875,81.84375,55.3359375

Tables can be access from this website :

https://disc.sci.gsfc.nasa.gov/datasets?page=1&keywords=MERRA-2

Open an account is required to access the data :

https://urs.earthdata.nasa.gov/users/new?client_id=C_kKX7TXHiCUqzt352ZwTQ&redirect_uri=https%3A%2F%2Fdisc.gsfc.nasa.gov%2Flogin%2Fcallback&response_type=code

The list of files and their variables are described in this publication :

https://gmao.gsfc.nasa.gov/pubs/docs/Bosilovich785.pdf

     -     2.     -
     
Once an account has been set, the script fetchMERRA2.py in this repository (git clone https://github.com/lsst-dm/fetch_merra2) can be used to upload all the atmospheric parameters above the LSST site for a list of date.
For exemple :

fetchMERRA2.py --site lsst --parameters aod pwv ozone --dir . --date 2018-3-27 2018-3-28
By default 

fetchMERRA2.py

fetchMERRA2.py --site lsst --parameters aod pwv ozone --dir . --date [last day of last month (which is the the latest available date)]

     -    3 .    -
     
     
The extration from the HDF5 files of the parameters' value above the site is performed :    

     runAll.py lsst . 20180504

which returns the lat-lon interpolation and z-integration of Ozone, PWV and the 15 AOD species (17 files)
Ozone, PWV, AOD :



