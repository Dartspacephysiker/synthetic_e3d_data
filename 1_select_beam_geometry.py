"""
In this script one defines the radar site locations, and designates each site as a transmitter, a receiver, or both.

By default, the sites are Skibotn (Tx/Rx), Kaiseniemi (Rx), and Kaaresuvanto (Rx)

This script will output TX_ECEF.txt, RX_ECEF.txt, and scatterpoints_ECEF.txt, which will then be used in 2_get_noise_estimates.R

S. M. Hatch
2023/05/02
"""

import numpy as np

# MPL stuff
import matplotlib as mpl
mplBkgrnd = 'QtAgg'
mpl.use(mplBkgrnd)
from mpl_toolkits.mplot3d import Axes3D

import matplotlib.pyplot as plt
plt.ion()

from radar import *

outdir = './output/'

txoutfile = 'TX_ECEF.txt'
rxoutfile = 'RX_ECEF.txt'
scatterptoutfile = 'scatterpoints_ECEF.txt'

saverxout = True
savetxout = True
savescatterout = True

##
# Geomagnetic reference radius:
RE = 6371.2 # km

# Skibotn
SKI = (69.58333, 20.46667)
SKI = (69.3401035, 20.315087)  # Anders's suggestion from Google Maps

KAI = (68.29,19.45)        # Kaiseniemi
KRS = (68.44933, 22.48325) # Karesuvanto

# Use Skibotn as transmitter
TX = SKI
RXs = [SKI,KAI,KRS]

## Select beam azimuths and elevations
config = 1

if config == 1:
    outputdir = outdir+''

    # use points given in Y. Ogawa-san's common mode presentation
    # https://www.space.irfu.se/workshops/EISCAT-3D_User2021-3/20211130_EISCAT_3D_User_meeting_yogawa_v4s.pdf

    el_arr1=np.array([64, 61, 60, 58, 57, 55, 54, 54, 57, 59, 61, 61])
    az_arr1=np.array([0, 35, 69, 101, 130, 156, 180, 204, 231, 258,288, 323])
    el_arr2=np.array([30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30])
    az_arr2=np.array([0, 30, 60, 90, 120, 150, 180, 210, 240, 270,300, 330])
    el_arr3=np.array([66, 77.8, 90])
    az_arr3=np.array([180, 180, 180])

elif config == 2:
    # Here we tighten up the elevations
    outputdir = outdir+'config2_tighter_beam/'

    el_arr1=np.array([64, 61, 60, 58, 57, 55, 54, 54, 57, 59, 61, 61])+12
    az_arr1=np.array([0, 35, 69, 101, 130, 156, 180, 204, 231, 258,288, 323])
    el_arr2=np.array([30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30])+24
    az_arr2=np.array([0, 30, 60, 90, 120, 150, 180, 210, 240, 270,300, 330])
    el_arr3=np.array([66+10, 77.8+5, 90])
    az_arr3=np.array([180, 180, 180])
    
elif config == 3:
    # Here we tighten up the elevations AND shift the EISCAT_3D sites 3 degrees northward
    outputdir = outdir+'config3_tighter_beam_shift_loc/'

    el_arr1=np.array([64, 61, 60, 58, 57, 55, 54, 54, 57, 59, 61, 61])+12
    az_arr1=np.array([0, 35, 69, 101, 130, 156, 180, 204, 231, 258,288, 323])
    el_arr2=np.array([30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30])+24
    az_arr2=np.array([0, 30, 60, 90, 120, 150, 180, 210, 240, 270,300, 330])
    el_arr3=np.array([66+10, 77.8+5, 90])
    az_arr3=np.array([180, 180, 180])
    
    # Shift everything up three degrees
    SKI = (SKI[0]+3,SKI[1])
    KAI = (KAI[0]+3,KAI[1])
    KRS = (KRS[0]+3,KRS[1])

    TX = SKI
    RXs = [SKI,KAI,KRS]


print(f"config: {config}")


az_arr = np.concatenate((az_arr1,az_arr2,az_arr3))
el_arr = np.concatenate((el_arr1,el_arr2,el_arr3))

#Get lots of heights
h_km = [80,90,100,110,120,200,250,300,350,400]
heightcol = ['blue','blue','blue','blue','blue',
             'orange','orange','orange','orange','orange']
heightm = ['o','o','o','o','o',
           '^','^','^','^','^']
    
# h_km = [80,90,100,110,120,130,140,150,160,170,180,190,
#         200,225,250,275,300,325,350,375,400]
# heightcol = ['blue']*len(h_km)
# heightm = ['o']*len(h_km)

r_los = []
cols = []
ms = []
azs = []
els = []
hs = []
for ih,h in enumerate(h_km):
    tmprlos = get_range_line(*TX, az_arr, el_arr, h,returnbonus=False)
    nelem = len(tmprlos)
    r_los.append(tmprlos)
    azs.append(az_arr)
    els.append(el_arr)
    hs.append([h]*nelem)

    cols.append([heightcol[ih]]*nelem)
    ms.append([heightm[ih]]*nelem)

r_los = np.vstack(r_los)
cols = np.hstack(cols)
ms = np.hstack(ms)
azs = np.hstack(azs)
hs = np.hstack(hs)
els = np.hstack(els)

xlos, ylos, zlos = r_los.T

## Write out receiver locations

R_Rs = []
for RX in RXs:
    gclatRec, gclonRec = RX
    eR, nR, uR = get_enu_vectors_cartesian(gclatRec,gclonRec,degrees=True)
    R_R = RE * uR
    R_Rs.append(R_R)
R_Rs = np.vstack(R_Rs)

Rxnames = np.array(['SKI','KAI','KRS'])
Rxarr = np.zeros(Rxnames.size, dtype=[('name', 'U6'), ('xECEF', float), ('yECEF', float), ('zECEF', float)])
Rxarr['name'] = Rxnames
Rxarr['xECEF'] = R_Rs[:,0]
Rxarr['yECEF'] = R_Rs[:,1]
Rxarr['zECEF'] = R_Rs[:,2]
if saverxout:
    print(f"Saving to {rxoutfile}")
    np.savetxt(outputdir+rxoutfile, Rxarr,
               fmt="%s,%8.3f,%8.3f,%8.3f",
               header="STATION,xECEF,yECEF,zECEF")

## Write out transmitter location

TXs = [SKI]

R_Ts = []
for TX in TXs:
    gclatRec, gclonRec = TX
    eR, nR, uR = get_enu_vectors_cartesian(gclatRec,gclonRec,degrees=True)
    R_T = RE * uR
    R_Ts.append(R_T)
R_Ts = np.vstack(R_Ts)

Txnames = np.array(['SKI'])
Txarr = np.zeros(Txnames.size, dtype=[('name', 'U6'), ('xECEF', float), ('yECEF', float), ('zECEF', float)])
Txarr['name'] = Txnames
Txarr['xECEF'] = R_Ts[:,0]
Txarr['yECEF'] = R_Ts[:,1]
Txarr['zECEF'] = R_Ts[:,2]
if savetxout:
    print(f"Saving to {txoutfile}")
    np.savetxt(outputdir+txoutfile, Txarr,
               fmt="%s,%8.3f,%8.3f,%8.3f",
               header="STATION,xECEF,yECEF,zECEF")


## Write out scatter locations

scatterarr = np.zeros(xlos.size, dtype=[('point', int), ('xECEF', float), ('yECEF', float), ('zECEF', float),
                                ('az', float), ('el', float),('h', float)])
scatterarr['point'] = np.arange(len(xlos),dtype=int)
scatterarr['xECEF'] = xlos
scatterarr['yECEF'] = ylos
scatterarr['zECEF'] = zlos
scatterarr['az']    = azs
scatterarr['el']    = els
scatterarr['h']     = hs
if savescatterout:
    print(f"Saving to {scatterptoutfile}")
    np.savetxt(outputdir+scatterptoutfile, scatterarr,
               fmt="%i,%8.3f,%8.3f,%8.3f,%7.2f,%7.2f,%7.2f",
               header="POINT,xECEF,yECEF,zECEF,AZ,EL,H_KM")

## Plot sampling locations
# So that we can plot a section of a sphere in our fancy 3D plot
def sphpoints(gclat,gclon,dlat=10,dlon=10):
    dlathaslen = hasattr(dlat,'__len__')
    dlonhaslen = hasattr(dlon,'__len__')

    if not dlathaslen:
        minlat,maxlat = gclat-dlat/2,gclat+dlat/2
    else:
        minlat,maxlat = gclat-dlat[0],gclat+dlat[1]
    minclat,maxclat = 90-maxlat,90-minlat

    if not dlonhaslen:
        minlon,maxlon = gclon-dlon/2,gclon+dlon/2
    else:
        minlon,maxlon = gclon-dlon[0],gclon+dlon[1]
        
    sphclat = np.deg2rad(np.linspace(minclat,maxclat,100))
    sphlon = np.deg2rad(np.linspace(minlon,maxlon,100))

    x = RE * np.outer(np.cos(sphlon), np.sin(sphclat))
    y = RE * np.outer(np.sin(sphlon), np.sin(sphclat))
    z = RE * np.outer(np.ones(np.size(sphlon)), np.cos(sphclat))

    return x,y,z

x,y,z = sphpoints(TX[0],TX[1],dlat=5,dlon=10)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

plotindices = slice(200,None,None)
plotindices = hs <= 120
plotindices = slice(None,None,None)
ax.scatter(xlos[plotindices], ylos[plotindices], zlos[plotindices],color=cols[plotindices])


# Plot the sphere surface
ax.plot_surface(x, y, z,color='orange')

ax.set_xlabel('$x_{\mathrm{ECEF}}$ [km]')
ax.set_ylabel('$y_{\mathrm{ECEF}}$ [km]')
ax.set_zlabel('$z_{\mathrm{ECEF}}$ [km]')

plt.show()

