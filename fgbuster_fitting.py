

import numpy as np

import pylab as pl
pl.rcParams['figure.figsize'] = 12, 16

import healpy as hp


import fgbuster.separation_recipes as sr
from  fgbuster.observation_helpers import get_instrument
from fgbuster.visualization import corner_norm
# Imports needed for component separation
from fgbuster import (CMB,  Synchrotron, FreeFree,  # sky-fitting model
                      MixingMatrix)  # separation routine
import warnings
warnings.filterwarnings("ignore")
import argparse
import os
from astropy.io import fits

def main(args) :
    try :
        os.makedirs(args.output_dir )
    except  FileExistsError:
        print (f"Warning: Overwriting files in {args.output_dir}")
    import pandas as pd


    hdul = fits.open(args.input_data)
    df = pd.DataFrame()
    df['frequency'] = np.float_(hdul[1].columns.names)/1e3
    data = hdul[1].data
    cols = hdul[1].columns.names
    nfreq=len(cols)
    df['depth_p'] =  np.ones(nfreq)
    df['depth_i'] = np.ones(nfreq)
    instrument = df.dropna(axis=1, how='all')
    hdul.close()

    #freq_maps = np.vstack([  data[c] for c in cols]  )
    freq_maps=hp.read_map(args.input_data, field=None )
    for i  in range(len(cols)):
        mask =np.ma.masked_invalid(freq_maps[i]).mask
        freq_maps[i][mask]= hp.UNSEEN

    freq_maps= freq_maps.reshape( freq_maps.shape[0],1,  freq_maps.shape[1])

    nside =args.nside
    options={'disp':False   , 'gtol': 1e-18, 'eps': 1e-18,
                'maxiter': 1000, 'ftol': 1e-18 }
    tol = 1e-18
    method='TNC'
    #freefree=  hp.read_map(filename=f"./lwa_data/COM_CompMap_freefree-commander_0256_R2.00.fits" ,field= ['EM_ML', 'TEMP_ML'] )
    components =[  Synchrotron(nu0=.040 , running=None  ,  units='K_RJ'  ) ]
    nsidepatches = [nside,4, nside  ]
    nsidegains = (nfreq-1)*[8]
    #nsidepatches += nsidegains

    import time
    s= time.perf_counter()

    results   = sr.multi_res_comp_sep(components, instrument, freq_maps, nsides= nsidepatches ,
                                    method=method,
                                    tol = tol, options=options)
    e= time.perf_counter()
    print(f"compsep took {e-s} sec. ")

    import pdb
    #pdb.set_trace()

    np.savez(f'{args.output_dir}/fgbuster_params_{args.label}.npz',
                    **{n: a for n, a in zip(results.params, results.x)})



if __name__=="__main__":
    parser = argparse.ArgumentParser( description="prepare training and testing dataset from a healpix map " )
    parser.add_argument("--output-dir" ,    help='path for outputs', default='./')
    parser.add_argument("--label" ,    help='path for outputs', default='res')
    parser.add_argument("--input-data" ,    help='path of frequency maps ', required=True)
    parser.add_argument("--nside", help="nside of output maps" ,required=True ,  type=np.int_)

    """
    parser.add_argument("--Bs-clusters" ,    help='path of Bs cluster patches', required=True)
    parser.add_argument("--Bd-clusters" ,   help='path of Bd cluster patches', required=True)
    parser.add_argument("--Td-clusters" ,    help='path of Td cluster patches', required=True)
    parser.add_argument("--polarization", help='compsep on polarization data ',
                                        action='store_true')
    parser.add_argument("--include-temperature",
                            help='add temperature data to compsep', action='store_true')
    parser.add_argument("--galmask", help = 'path to the  healpix galactic mask ', )
    """

    args = parser.parse_args()


    main( args)
