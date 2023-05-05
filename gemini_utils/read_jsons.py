import numpy as np
import sys
import json
from datetime import datetime,timedelta
from utils.utils import spherical_to_dipole,dipole_to_spherical,dipole_basis_vectors
from gemini_utils.convert import geomag2geog

Re = 6370e3                     # from convert.py

def read_cfg(cfgjson="./gemini_data/cfg.json"):
    with open(cfgjson, "r") as cfgfile:
        cfg = json.load(cfgfile)

    for key in cfg.keys():
        if isinstance(cfg[key],str):
            if cfg[key].startswith("deltat"):
                print("Converting "+key)
                cfg[key] = timedelta(int(cfg[key][6:]))
        elif key == 'time':
            cfg[key] = [datetime.strptime(cfg[key][i],"%Y%m%d %H:%M:%S") for i in range(len(cfg[key]))]
    
    return cfg


def read_xg(xgjson="./gemini_data/grid.json"):

    with open(xgjson, "r") as xgfile:
        xgin = json.load(xgfile)

    xg = {}
    for key in xgin:
        if key.endswith('shape'):
            continue
        print(f"Reading in {key}")
        xg[key] = np.array(xgin[key]).reshape(xgin[key+'_shape'])


    xg['r'], xg['theta'] = dipole_to_spherical(xg['x1'][2:-2,np.newaxis],xg['x2'][np.newaxis,2:-2])
    xg['r'] = xg['r'] * Re

    # xg['r'] = np.broadcast_to(xg['r'],(*xg['r'].shape,len(xg['x3'][2:-2])))
    # xg['theta'] = np.broadcast_to(xg['theta'],(*xg['theta'].shape,len(xg['x3'][2:-2])))
    # xg['phi'] = np.broadcast_to(xg['x3'][2:-2],(*xg['theta'].shape,len(xg['x3'][2:-2])))

    xg['r'] = np.repeat(xg['r'][:,:,np.newaxis],len(xg['x3'][2:-2]),axis=2)
    xg['theta'] = np.repeat(xg['theta'][:,:,np.newaxis],len(xg['x3'][2:-2]),axis=2)
    xg['phi'] = np.repeat(xg['x3'][2:-2][np.newaxis,:],len(xg['x1'][2:-2]),axis=0)
    xg['phi'] = np.repeat(xg['phi'][:,np.newaxis,:],len(xg['x2'][2:-2]),axis=1)

    xg['h1'] = np.repeat(xg['h1'],xg['x3'].size,axis=1).reshape((*xgin['h1_shape'],len(xg['x3'])))

    xg['glon'], xg['glat'] = geomag2geog(xg['phi'], xg['theta'])

    t, p = xg['theta'], xg['phi']
    xg['x'] = xg['r']*np.sin(t)*np.cos(p)
    xg['y'] = xg['r']*np.sin(t)*np.sin(p)
    xg['z'] = xg['r']*np.cos(t)

    # spherical basis vectors
    xg['er'    ] = np.array([ np.sin(t)*np.cos(p), np.sin(t)*np.sin(p),  np.cos(t)])
    xg['etheta'] = np.array([ np.cos(t)*np.cos(p), np.cos(t)*np.sin(p), -np.sin(t)])
    xg['ephi'  ] = np.array([-np.sin(p)          ,           np.cos(p),  np.sin(t)*0.])

    xg['er'] = np.transpose(xg['er'],axes=[1,2,3,0])
    xg['etheta'] = np.transpose(xg['etheta'],axes=[1,2,3,0])
    xg['ephi'] = np.transpose(xg['ephi'],axes=[1,2,3,0])

    # dipole basis vectors
    xg['e1'], xg['e2'] = dipole_basis_vectors(90.-np.degrees(xg['theta']),
                                              lon=np.degrees(xg['phi']),
                                              cartesian=True)
    xg['e3'] = xg['ephi']

    xg['alt'] = xg['r']-Re

    return xg

    
if __name__ == '__main__':
    xg = read_xg()
