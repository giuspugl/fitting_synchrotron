

import numpy as np

import pylab as pl
pl.rcParams['figure.figsize'] = 12, 16
import pandas as pd

import healpy as hp

import pdb
from astropy.table import Table

import fgbuster.separation_recipes as sr
# Imports needed for component separation
from fgbuster import (   Synchrotron, FreeFree,AnalyticComponent,  Component , # sky-fitting model
                      MixingMatrix)  # separation routine
import warnings
warnings.filterwarnings("ignore")
import argparse
import os
from astropy.io import fits
from  matplotlib import pyplot as plt
import time




class  LogFreeFree(AnalyticComponent):

    """ Log  Power law

    Parameters
    ----------
    nu0: float
        Reference frequency
    beta_pl: float
        Spectral index
    units:
        Output units (K_CMB and K_RJ available)
    """
    _REF_BETA = -2.14
    
    
    #logsynch=AnalyticComponent(analytic_expr, nu0 =0.408  )
    def __init__(self, nu0, beta_pl=None, 
                 units='K_CMB'):
        

        # Prepare the analytic expression
        analytic_expr = ('log(nu / nu0)*  beta_pl ')
        if 'K_CMB' in units:
            analytic_expr += ' * ' + K_RJ2K_CMB_NU0
        elif 'K_RJ' in units:
            pass
        else:
            raise ValueError("Unsupported units: %s"%units)

        kwargs = {'nu0': nu0,  
                  'beta_pl': beta_pl }

        super(LogFreeFree, self).__init__(analytic_expr, **kwargs)

        self._set_default_of_free_symbols(
            beta_pl=self._REF_BETA )


class  LogSynchrotron (AnalyticComponent):

    """ Log  Power law

    Parameters
    ----------
    nu0: float
        Reference frequency
    beta_pl: float
        Spectral index
    nu_pivot: float
        Pivot frequency for the running
    running: float
        Curvature of the power law
    units:
        Output units (K_CMB and K_RJ available)
    """
    _REF_BETA = -3
    _REF_RUN = 0.
    _REF_NU_PIVOT = 70.
    
    #logsynch=AnalyticComponent(analytic_expr, nu0 =0.408  )
    def __init__(self, nu0, beta_pl=None, nu_pivot=None, running=0.,
                 units='K_CMB'):
        if nu_pivot == running == None:
            print('Warning: are you sure you want both nu_pivot and the running'
                  'to be free parameters?')

        # Prepare the analytic expression
        analytic_expr = ('log(nu / nu0)* ( (beta_pl)  + running * log (nu/nu_pivot) ) ')
        if 'K_CMB' in units:
            analytic_expr += ' * ' + K_RJ2K_CMB_NU0
        elif 'K_RJ' in units:
            pass
        else:
            raise ValueError("Unsupported units: %s"%units)

        kwargs = {'nu0': nu0, 'nu_pivot': nu_pivot,
                  'beta_pl': beta_pl, 'running': running}

        super(LogSynchrotron, self).__init__(analytic_expr, **kwargs)

        self._set_default_of_free_symbols(
            beta_pl=self._REF_BETA, running=self._REF_RUN, nu_pivot=self._REF_NU_PIVOT)

def main(args) :
 

    try :
        os.makedirs(args.output_dir )
    except  FileExistsError:
        print (f"Warning: Overwriting files in {args.output_dir}")


    hdul = fits.open(args.input_data)
    df = pd.DataFrame()
    df['frequency'] = np.float_(hdul[1].columns.names )/1e3
    data = hdul[1].data
    cols = hdul[1].columns.names
    nfreq=len(cols)
    mapserr_table = np.load(args.input_errmaps)

    df['depth_i']= [v for v in mapserr_table.values() ]
    instrument = df.dropna(axis=1, how='all')
    hdul.close()

    freq_maps=hp.read_map(args.input_data, field=cols  )

    print(instrument)
    nu0 = 0.408
    idmap =  np.argmin(np.fabs(df['frequency'] - nu0 ))
    freq_maps = np.log(freq_maps/freq_maps[idmap]) 
    j=0
    freq_maps_ud=[]
    for i  in range(nfreq ):
        hp.mollview(   freq_maps[i ]     ,  title=f'{cols[i] } MHz'  , sub=(5,4,1+j) , notext=True)
        freq_maps_ud.append( hp.ud_grade(freq_maps[i] , nside_out=64))

        j+=1
    plt.show()

    freq_maps=np.vstack(freq_maps_ud)
    freq_maps= freq_maps.reshape( freq_maps.shape[0],1,  freq_maps.shape[1])

    nside =args.nside
    options={'disp':False   , 'gtol': 1e-18, 'eps': 1e-18,
                'maxiter': 1000, 'ftol': 1e-18 }
    tol = 1e-18
    method='TNC'
    components =[  LogSynchrotron(nu0=.408 , running=None   ,  units='K_RJ'  ) ,
                  LogFreeFree(nu0=2. ,  units='K_RJ'  ) ]
    pdb.set_trace() 

    g0 = list(np.ones(nfreq-1))
    sp0 = [ -3,60,0. ]
    x0 = sp0 + g0

    bounds_g = [(None, None) for i in g0]
    bounds_sp = [(-4, 0.), (10., 200.), (-.5, .5 )]
    bounds = bounds_sp  + bounds_g


    nsidegains = list(np.zeros(len(g0)))

    nsidepatches = [nside,4  , nside  ]

    nsidepatches = nsidepatches + nsidegains

    s= time.perf_counter()
    #pdb.set_trace()
    gain= False 
    if gain :
        g0 = list(np.ones(nfreq-1))
        sp0 = [ -3,60,0. ]
        x0 = sp0 + g0

        bounds_g = [(None, None) for i in g0]
        bounds_sp = [(-4.5, 0.), (10., 200.), (-.5, .5 )]
        bounds = bounds_sp  + bounds_g


        nsidegains = list(np.zeros(len(g0)))

        nsidepatches = [nside,4  , nside  ]

        nsidepatches = nsidepatches + nsidegains
        results   = sr.multi_res_comp_sep_gain(components, instrument, freq_maps, nsides= nsidepatches ,
                                        x0=x0,  
                                        method=method,tol = tol,  options=options,
                                        #bounds=bounds,
                                        )
    else : 
        x0 = [ -3,60,0. ]
        x0 = sp0  

        bounds = [(None , 0.), (None , None ), (-.5, .5 )]
        nsidepatches = [nside,4  , nside ,4  ]

        results   = sr.multi_res_comp_sep(components, instrument, freq_maps, nsides= nsidepatches ,
                                    method=method,tol = tol, options=options,
                                    #bounds=bounds,
                                    )
    e= time.perf_counter()
    print(f"compsep took {e-s} sec. ")


    np.savez(f'{args.output_dir}/fgbuster_params_{args.label}_nside{args.nside}.npz',
                    **{n: a for n, a in zip(results.params, results.x)})



if __name__=="__main__":
    parser = argparse.ArgumentParser( description="prepare training and testing dataset from a healpix map " )
    parser.add_argument("--output-dir" ,    help='path for outputs', default='./')
    parser.add_argument("--label" ,    help='path for outputs', default='res')
    parser.add_argument("--input-data" ,    help='path of frequency maps ', required=True)
    parser.add_argument("--input-errmaps" ,    help='path of frequency maps errors ', required=True)
    parser.add_argument("--nside", help="nside of output maps" ,required=True ,  type=np.int_)

    args = parser.parse_args()


    main( args)
