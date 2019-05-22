# fetch_merra2


Accessing and processing MERRA-2 data

NASA provides an online map of worldwide atmospheric parameters :
https://worldview.earthdata.nasa.gov/?p=geographic&l=VIIRS_SNPP_CorrectedReflectance_TrueColor(hidden),MODIS_Aqua_CorrectedReflectance_TrueColor(hidden),MODIS_Terra_CorrectedReflectance_TrueColor,MODIS_Terra_Aerosol,Reference_Labels,Reference_Features,Coastlines&t=2017-05-02&z=3&v=-125.4375,-48.8671875,81.84375,55.3359375

Tables can be access from this website :
https://disc.sci.gsfc.nasa.gov/datasets?page=1&keywords=MERRA-2

Open an account is required to access the data :
https://urs.earthdata.nasa.gov/users/new?client_id=C_kKX7TXHiCUqzt352ZwTQ&redirect_uri=https%3A%2F%2Fdisc.gsfc.nasa.gov%2Flogin%2Fcallback&response_type=code

The list of files and their variables are described in this publication :

https://gmao.gsfc.nasa.gov/pubs/docs/Bosilovich785.pdf

git clone https://github.com/lsst-dm/fetch_merra2

fetchMERRA2.py --date 2018-05-04 --site ctio
fetchMERRA2.py --site ctio --parameters aod --dir . --date 2018-4-27
fetchMERRA2.py --site ctio --date 2018-4-27


Ozone, PWV, AOD :
merra2-chain.py --site CTIO --value TO3 --files MERRA2_400.inst1_2d_asm_Nx.20180504.SUB.nc
equivalent to :
merra2-chain.py --files MERRA2_400.inst1_2d_asm_Nx.20180504.SUB.nc
merra2_3D_z.py --site CTIO --value QV --files MERRA2_400.inst3_3d_asm_Nv.20180504.SUB.nc

-->
runAll.py ctio . 20180504


