# fetch_merra2

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


