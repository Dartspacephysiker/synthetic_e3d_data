"""

"""

#OPPLEGG
# 1. Sample alle punktene hvor vi vil ha data ved å bruke punktene vi ble enige om inn i GEMINI, sjekk hvor vi har data
#    a. Get the points we need
#    b. Now get model quantities at relevant locations
# 2. Gir det mening å snakke om et gjennomsnitt inne i hver sampling-volum? Hvor stor er ellipsoide ift GEMINI-oppløsning?
# 3. Sammenligne tetthetene med tettheten som Ilkka har brukt, høre med ham
# 4. Gjøre et forsøk med rekonstruksjon, kanskje strøm?

import numpy as np

import pandas as pd
import xarray as xr

from utils.utils import ECEF2geodetic,spherical_to_dipole,dipole_to_spherical
from utils.uncertainty import func_gammaT_div_gamma0, vError, isoError

import gemini_utils.read as read
from gemini_utils.gridmodeldata import model2pointsgeogcoords
from gemini_utils.convert import unitvecs_geographic,geog2geomag
from gemini_utils.read_jsons import read_cfg,read_xg

subdir = ''
inputdir = './'+subdir
plotdir = './'+subdir

isotropic_rvalue_file = inputdir+"output/isotropicnoise.txt"
velocity_rvalue_file = inputdir+"output/velocitynoise.txt"
scatterptoutfile = inputdir+'output/scatterpoints_ECEF.txt'

Re = 6370e3                     # from gemini_utils/convert.py

gemini_data_dir = inputdir+"gemini_data/"


# Stuff for uncertainty calculation
errtablefile = inputdir+'ErrorTable_E3D_4params.txt'

integ_T_min = 5
integ_T_sec_iso = 60*integ_T_min
integ_T_sec_vel = 60*integ_T_min

fwhmRange = 3
dutyCycle = 0.25

# Output file
outnetcdf = inputdir+f'output/E3D__GEMINI_aurora_EISCAT3D_simple_wide_samples__{integ_T_min}min_integ.nc'

elem_charge = 1.602176634e-19

###############################################################################
# load GEMINI config and grid
# cfg = read.config(gemini_data_dir)
cfg = read_cfg()
# xg = read.grid(gemini_data_dir)
xg = read_xg()

# for key in xg2.keys():
#     print(key,np.all(np.isclose(xg[key],xg2[key])))

getems = ['ne','Ti','Te','v1','v2','v3','J1','J2','J3',
          'vperpgalt','vperpglat','vperpglon',
          'vperpr','vperptheta','vperpphi',
          've2','ve3',
          'Opratio']

[egalt,eglon,eglat]=unitvecs_geographic(xg)    
#^ returns a set of geographic unit vectors on xg; these are in ECEF geomag comps
#    like all other unit vectors in xg

###############################################################################
# 1. Sample all points where we want data by using points vi ble enige om inn i GEMINI, sjekk hvor vi har data

##############################
# a. Get the points we need
Txname = 'SKI'
Rxname = 'KAI' #Or 'KAI' or 'KRS"

pt = pd.read_csv(scatterptoutfile)
npt = pt.shape[0]
coordvars = ['# POINT', 'xECEF', 'yECEF', 'zECEF', 'AZ', 'EL', 'H_KM','beam','gdlat','glon']  # Use these to define coordinates in xarray dataset

coordvars = ['# POINT', 'xECEF', 'yECEF', 'zECEF', 'AZ', 'EL', 'H_KM','beam','gdlat','glon',
             'mlon','mlat','phi','theta','r','x1','x2','x3','q','p',
             'gemini_x1_i','gemini_x2_i','gemini_x3_i']  # Use these to define coordinates in xarray dataset


nh = pt['H_KM'].unique().size

beamarr = np.zeros(npt,dtype=int)
nbeams = npt//np.unique(pt['H_KM']).size
for i in range(nbeams):
    beamarr[(pt['AZ'].values == pt['AZ'].values[i]) & (pt['EL'].values == pt['EL'].values[i])] = i
pt['beam'] = beamarr

h_km,gdlat,glon = ECEF2geodetic(pt['xECEF'].values, pt['yECEF'].values, pt['zECEF'].values, degrees = True)

pt['gdlat'] = gdlat
pt['glon'] = glon

pt['phi'], pt['theta'] = geog2geomag(pt['glon'].values, pt['gdlat'].values)
pt['r'] = np.sqrt((pt[['xECEF','yECEF','zECEF']]**2).sum(axis=1)).values*1e3
pt['mlon'], pt['mlat'] = np.degrees(pt['phi']), 90.-np.degrees(pt['theta'])  # 'cause it's a dipole field

# 'x1' and 'x2' correspond to coordinates 'q' and 'p' (first two dipole coordinates) in GEMINI document
# dipole coordinates
pt['x1'] = (Re/pt['r'].values)**2 * np.cos(pt['theta'])  # eq 81, Zettergren (2019)
pt['x2'] = (pt['r']/Re)/ (np.sin(pt['theta']))**2        # eq 82, Zettergren (2019)
pt['x3'] = pt['phi']

pt['q'] = (Re/pt['r'].values)**2 * np.cos(pt['theta'])# eq 81, Zettergren (2019)
pt['p'] = (pt['r']/Re)/ (np.sin(pt['theta']))**2      # eq 82, Zettergren (2019)

# Check conversion
rback, thetaback = dipole_to_spherical(pt['x1'].values,pt['x2'].values)
assert np.all(np.isclose(pt['r'].values/Re,rback))
assert np.all(np.isclose(pt['theta'].values,thetaback))

##############################
# b. Now get model quantities at relevant locations

ntottimes = len(cfg['time'])

ds = []

# for i_time in range(ntimes):
for i_time in range(30):

    tim = cfg["time"][i_time]
    print(f"{i_time:02d} {tim}")

    try:
        dat = read.frame(gemini_data_dir, tim, cfg=cfg, xg=xg)
    except:
        print(f"Couldn't get data for {tim.strftime('%Y%m%d %H:%M:%S')}, skipping!")

    if 'ns' not in list(dat.keys()):
        print("Don't have ion densities, skipping this time step!")
        continue

    # Velocity, B-field, E-field in ENU coords
    v1=dat["v1"]; v2=dat["v2"]; v3=dat["v3"];

    # #B vectors
    # Bgalt= np.sum(xg["e1"]*egalt*xg['Bmag'][...,np.newaxis],3)
    # Bglat= np.sum(xg["e1"]*eglat*xg['Bmag'][...,np.newaxis],3)
    # Bglon= np.sum(xg["e1"]*eglon*xg['Bmag'][...,np.newaxis],3)
    
    # each of the components in models basis projected onto geographic unit vectors
    # vgalt=( np.sum(xg["e1"]*egalt,3)*v1 + np.sum(xg["e2"]*egalt,3)*v2 + 
    #         np.sum(xg["e3"]*egalt,3)*v3 )
    # vglat=( np.sum(xg["e1"]*eglat,3)*v1 + np.sum(xg["e2"]*eglat,3)*v2 +
    #         np.sum(xg["e3"]*eglat,3)*v3 )
    # vglon=( np.sum(xg["e1"]*eglon,3)*v1 + np.sum(xg["e2"]*eglon,3)*v2 + 
    #         np.sum(xg["e3"]*eglon,3)*v3 )

    vperpgalt=( np.sum(xg["e2"]*egalt,3)*v2 + 
            np.sum(xg["e3"]*egalt,3)*v3 )
    vperpglat=( np.sum(xg["e2"]*eglat,3)*v2 +
            np.sum(xg["e3"]*eglat,3)*v3 )
    vperpglon=( np.sum(xg["e2"]*eglon,3)*v2 + 
            np.sum(xg["e3"]*eglon,3)*v3 )
    
    vperpr = ( np.sum(xg["e2"]*xg['er'],3)*v2 + 
            np.sum(xg["e3"]*xg['er'],3)*v3 )
    vperptheta = ( np.sum(xg["e2"]*xg['etheta'],3)*v2 + 
            np.sum(xg["e3"]*xg['etheta'],3)*v3 )
    vperpphi = ( np.sum(xg["e2"]*xg['ephi'],3)*v2 + 
            np.sum(xg["e3"]*xg['ephi'],3)*v3 )
    
    # dat = dat.assign(vgalt=dat['v1']*0+vgalt)
    # dat = dat.assign(vglat=dat['v1']*0+vglat)
    # dat = dat.assign(vglon=dat['v1']*0+vglon)
    # dat = dat.assign(oplus=dat['v1']*0+dat['ns'][0,:,:,:])

    dat = dat.assign(vperpgalt=dat['v1']*0+vperpgalt)
    dat = dat.assign(vperpglat=dat['v1']*0+vperpglat)
    dat = dat.assign(vperpglon=dat['v1']*0+vperpglon)
    dat = dat.assign(vperpr=dat['v1']*0+vperpr)
    dat = dat.assign(vperptheta=dat['v1']*0+vperptheta)
    dat = dat.assign(vperpphi=dat['v1']*0+vperpphi)
    
    dat = dat.assign(ve2=dat['v1']*0+dat['v2']-dat['J2']/(elem_charge*dat['ne']))
    dat = dat.assign(ve3=dat['v1']*0+dat['v3']-dat['J3']/(elem_charge*dat['ne']))
    
    dat = dat.assign(Opratio=dat['ne']*0+dat['ns'].values[0]/dat['ne'])

    # # Make B vectors and v vectors
    # Bvec_ENU = np.stack([Bglon,Bglat,Bgalt],axis=-1)
    # vvec_ENU = np.stack([vglon,vglat,vgalt],axis=-1)

    # # Get E-field 
    # Evec_ENU = -np.cross(vvec_ENU,Bvec_ENU,axis=-1)*1e3  # mV/m

    # model2pointsgeogcoords wants altitude in METERS
    for getem in getems:
        pt[getem] = model2pointsgeogcoords(xg,dat[getem],h_km*1000,glon,gdlat)

    tmpdict = dict()
    tmpcoords = dict()
    for col in pt.columns:
        if col in coordvars:
            if col == 'H_KM':
                tmpcoords[col] = pt[col].unique()
            elif col == 'beam':
                tmpcoords[col] = np.unique(beamarr)
            else:
                tmpcoords[col] = (['H_KM','beam'],pt[col].values.reshape(nh,nbeams))
        else:
            tmpdict[col] = (['H_KM','beam','time'],pt[col].values.reshape(nh,nbeams,1))

    tmpcoords['time'] = [tim]
    tmpdict['itime'] = ('time',np.array([i_time]))

    oneds = xr.Dataset(
        data_vars=tmpdict,
        coords=tmpcoords,
    )

    ds.append(oneds)

    # if i_time > 1:
    #     break

ds = xr.concat(ds,'time')
ds = ds.rename({"# POINT":"ptindex"})

ds = ds.assign(vmag=np.sqrt(ds['v1']**2+ds['v2']**2+ds['v3']**2))

##############################
# Now calculate uncertainties

# Load noise estimates γ/γ_0 outputted from get_noise_estimates.R
iso_rvalue = pd.read_csv(isotropic_rvalue_file,header=None)
vel_rvalue = pd.read_csv(velocity_rvalue_file,header=None)

# Load error table
print(f"Loading error table file '{errtablefile}'")
errtab = pd.read_csv(errtablefile,delimiter='\s+',header=11)
uniq = {lab:errtab[lab].unique() for lab in errtab.columns}
whereuniq = dict()              # use this dict to play matchie-match
for lab in errtab.columns:
    nuniq = len(uniq[lab])
    tmpdict = dict()
    for i in range(nuniq):
        tmpdict[i] = np.where(errtab[lab] == uniq[lab][i])[0]
    whereuniq[lab] = tmpdict


# Need to find where GEMINI quantities are closest to values in table
neinds = np.argmin(np.abs(ds['ne'].values.ravel()[:,np.newaxis]-uniq['Ne'][np.newaxis,:]),axis=1)
Tiinds = np.argmin(np.abs(ds['Ti'].values.ravel()[:,np.newaxis]-uniq['Ti'][np.newaxis,:]),axis=1)
Teinds = np.argmin(np.abs((ds['Te']/ds['Ti']).values.ravel()[:,np.newaxis]-uniq['Te/Ti'][np.newaxis,:]),axis=1)
Viinds = np.argmin(np.abs(ds['vmag'].values.ravel()[:,np.newaxis]-uniq['Vi'][np.newaxis,:]),axis=1)
Oinds = np.argmin(np.abs(ds['Opratio'].values.ravel()[:,np.newaxis]-uniq['[O+]'][np.newaxis,:]),axis=1)
collinds = np.zeros_like(ds['ne'].values.ravel(),dtype=int)
collinds[:] = 1                 # corresponds to collfreq = 500


NeN0err = np.zeros_like(ds['ne'].values.ravel())
TiT0err = np.zeros_like(ds['ne'].values.ravel())
TeTierr = np.zeros_like(ds['ne'].values.ravel())
ViV0err = np.zeros_like(ds['ne'].values.ravel())
finalcandidates = []
for i in range(len(neinds)):
    neind = neinds[i]
    Tiind = Tiinds[i]
    Teind = Teinds[i]
    Viind = Viinds[i]
    Oind = Oinds[i]
    collind = collinds[i]

    candidates = whereuniq['Ne'][neind]
    candidates = candidates[np.in1d(candidates,whereuniq['Ti'][Tiind])]
    candidates = candidates[np.in1d(candidates,whereuniq['Te/Ti'][Teind])]
    candidates = candidates[np.in1d(candidates,whereuniq['Vi'][Viind])]
    candidates = candidates[np.in1d(candidates,whereuniq['[O+]'][Oind])]
    candidates = candidates[np.in1d(candidates,whereuniq['Coll'][collind])]

    finalcandidates.append(candidates)
    
finalcandidates = np.array(finalcandidates).ravel()
    # assert 2<0    

ntimes = ds['time'].size
intermediate_shape = (nh*nbeams,ntimes)
final_shape = (nh,nbeams,ntimes)

ds = ds.assign(errNeN0=(('H_KM','beam','time'),
                        errtab.iloc[finalcandidates]['Err(Ne/N0)'].values.reshape(final_shape)))
ds = ds.assign(errTiT0=(('H_KM','beam','time'),
                        errtab.iloc[finalcandidates]['Err(Ti/T0)'].values.reshape(final_shape)))
ds = ds.assign(errTeTi=(('H_KM','beam','time'),
                        errtab.iloc[finalcandidates]['Err(Te/Ti)'].values.reshape(final_shape)))
ds = ds.assign(errViV0=(('H_KM','beam','time'),
                        errtab.iloc[finalcandidates]['Err(Vi/Vi0)'].values.reshape(final_shape)))
    
ds = ds.assign(integ_T_sec_vel=integ_T_sec_vel)
ds = ds.assign(integ_T_sec_iso=integ_T_sec_iso)
ds = ds.assign(fwhmRange=fwhmRange)
ds = ds.assign(dutyCycle=dutyCycle)

# put actual uncertainties into dataset BUT NOTE, THESE AREN'T QUITE RIGHT
print("Calculating uncertainties ...")
ds = ds.assign(vel_unc=(['H_KM','beam','time'],
                        vError(func_gammaT_div_gamma0(
                            vel_rvalue.values.ravel()[:,np.newaxis],
                            integ_T_sec_vel,fwhmRange=fwhmRange,dutyCycle=dutyCycle),
                               ds['errViV0'].values.reshape((intermediate_shape)),
                               Vi0=1).reshape(final_shape)))

ds = ds.assign(ne_unc=(['H_KM','beam','time'],
                       isoError(
                           func_gammaT_div_gamma0(vel_rvalue.values.ravel()[:,np.newaxis],
                                                  integ_T_sec_iso,fwhmRange=fwhmRange,dutyCycle=dutyCycle),
                           ds['ne'].values.reshape(intermediate_shape),
                           ds['errNeN0'].values.reshape((intermediate_shape))).reshape((final_shape))))

gammaT_div_gamma0_iso = func_gammaT_div_gamma0(iso_rvalue.values.ravel()[:,np.newaxis],
                                               integ_T_sec_iso,fwhmRange=fwhmRange,dutyCycle=dutyCycle)

ds = ds.assign(ne_unc=(['H_KM','beam','time'],
                       isoError(
                           gammaT_div_gamma0_iso,
                           ds['ne'].values.reshape(intermediate_shape),
                           ds['errNeN0'].values.reshape((intermediate_shape))).reshape((final_shape))
                       ))
ds = ds.assign(Ti_unc=(['H_KM','beam','time'],
                       isoError(gammaT_div_gamma0_iso,
                                ds['Ti'].values.reshape(intermediate_shape),
                                ds['errTiT0'].values.reshape((intermediate_shape))).reshape((final_shape))
                       )
               )
ds = ds.assign(Te_unc=(['H_KM','beam','time'],
                       (ds['Te']/ds['Ti']*ds['Ti_unc'].values+ \
                        ds['Ti'] * \
                        isoError(gammaT_div_gamma0_iso,
                                 (ds['Te']/ds['Ti']).values.reshape(intermediate_shape),
                                 ds['errTeTi'].values.reshape((intermediate_shape))).reshape((final_shape))
                       ).values)
               )


nstride = int(integ_T_min*60/int(np.diff(ds['time'][:2].values.ravel())[0]/1e9))

dsroll = ds.rolling(time=nstride+1, center=False, min_periods=nstride+1).mean()

dsroll = dsroll.isel(time=slice(None,None,nstride))
print(f"Saving it all to {outnetcdf}")
dsroll.to_netcdf(outnetcdf)


