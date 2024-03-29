# Copyright (C) 2015 Associated Universities, Inc. Washington DC, USA.
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope. that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
# 
# Correspondence concerning GBT software should be addressed as follows:
#       GBT Operations
#       National Radio Astronomy Observatory
#       P. O. Box 2
#       Green Bank, WV 24944-0002 USA

import math
import numpy as np
import sys, os
from scipy import interpolate
import scipy
from ctypes import *

# Loading functions in C++ library
array_1d_float  = np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,flags="CONTIGUOUS")
array_1d_double = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,flags="CONTIGUOUS")

libdir = os.path.dirname(os.path.realpath(__file__))
try:
    libobj = np.ctypeslib.load_library("grid_otflib.so",libdir)
    C_available = True
    libobj.grid_otf_C.restype = POINTER(c_float)
    libobj.grid_otf_C.argtypes = [array_1d_float,c_size_t,c_size_t,
                                  array_1d_double,array_1d_double,array_1d_double,
                                  c_size_t,c_size_t,c_char_p,c_double,c_double,array_1d_double,
                                  c_int,c_double,c_double,c_double,c_int]
except:
    print ("Cannot find C library, C functions will not be available.")
    C_available = False


def grid_otf(data, xsky, ysky, wcsObj, nchan, xsize, ysize, pix_scale, weight=None, beam_fwhm=None,
             kern="gaussbessel", gauss_fwhm=None, verbose=4, useC=True, threads=1):
    """
    Grid individual spectra onto a specified regular grid following the 
    recommendations of Mangum et al. (2007).  This is writtent to be of general
    use in gridding OTF data.

    Adapted from IDL code provided by Adam Leroy (aleroy@nrao.edu).

    Inputs:
       data - an nspec by nchan shaped 2-d array containing the spectra to be gridded
       xsky - an nspec length vector of X positions on the sky for each spectra in data (deg)
       ysky - an nspec length vector of Y positions on the sky for each spectra in data (deg)
       wcsObj - a wcs object appriopriate for gridding xsky and ysky on to
       weight - (optional) an nspec length vector of weights for each spectra in data.
                if not supplied, equal weights are assumed
       beam_fwhm - the fwhm in decimal degrees of the telescope beam.  If not
                   supplied it must be found in the target_hdr.
       gauss_fwhm - the fwhm in decimal degrees of the gaussian used in the
                    convolution kernel.  Used only when kern="gauss".
                    Defaults to beam_fwhm/3 if not set.
       kern - specify the gridding kernel to use from "gaussbessel", "gauss", or
              "nearest".
              defaults to "gaussbessel"

    Returns: (cube, weight, beam_fwhm) where cube is the cube array after gridding and
       weight is the related weight array and beam_fwhm is effective fwhm of the beam
       after the convolution.  Returns (None, None, None) on failure.

    Note: apodize and median averaging are not implemented here
          weight is current per spectra only

    Original IDL credits to CLASS's xy_map, E. Rosolowsky
    """

    result = (None,None,None)

    # argument checking
    # data must be 2D, nspec and nchan must be non-zero, nspec must
    # agree with length of xsky and ysky
    # it's assumed these are np arrays
    if len(data.shape) != 2 or len(xsky.shape) != 1 or len(ysky.shape) != 1:
        if verbose > 1:
            print("data, sky coordinates have unexpected shapes")
            print("data : ", data.shape)
            print("xsky : ", xsky.shape)
            print("ysky : ", ysky.shape)
        return result

    nspec, nchan_data = data.shape
    if nspec == 0 or nchan_data == 0:
        if verbose > 1:
            print("no data given")
        return result

    if nspec != len(xsky) or nspec != len(ysky):
        if verbose > 1:
            print("Number of sky position values does not match number of spectra in data")
        return result

    if beam_fwhm is None:
        if 'BMAJ' not in target_hdr:
            if verbose > 1:
                print("BMAJ must be in target_hdr or supplied as an argument")
            return result
        beam_fwhm = target_hdr['BMAJ']

    if kern not in ["gaussbessel","gauss","nearest"]:
        if verbose > 1:
            print("kern must be one of gaussbessel or gauss")
        return result

    # use the python shape - fastest changing axis is last
    cubeShape = (nchan,ysize,xsize)

    # frequency axis should agree with data
    if cubeShape[0] != nchan_data:
        if verbose > 1:
            print("Frequency axis in target header and spectra length do not match")
        return result

    # does not yet support per channel weights
    if weight is None:
        # equal weight = all 1s
        weight = np.ones(nspec)

    # need to watch for NaNs here - weight those to zero
    # current managed upstream when the data are read ...

    # determine the convolution function
    if kern == "gauss":
        if gauss_fwhm is None:
            gauss_fwhm = beam_fwhm / 3.0
        # original calculation
        # r_support = 3.0 * gauss_fwhm
        # this is what we do in AIPS - 5 pixels
        r_support = 5.0 * pix_scale
        # Adam has max_conv_fn = 0.5 here, that seems wrong as this function peaks at 1.0
        max_conv_fn = 1.0
        # this seems extreme - from Adam, especially since has the max at 0.5
        # cutoff_conv_fn = 0.25*max_conv_fn
        # try this instead - I think this is closer to what AIPS does
        cutoff_conv_fn = 0.0
        scale_fwhm = math.sqrt(beam_fwhm**2+gauss_fwhm**2)/beam_fwhm
    elif kern == "gaussbessel":
        # gaussian*bessel from Mangum, Emerson, Greisen (2007)
        # this still needs to be checked against the default AIPS values
        a = 1.55 * (beam_fwhm) / 3.0
        b = 2.52 * (beam_fwhm) / 3.0
        # original calculation
        # r_support = beam_fwhm
        # aips - 3 pixels
        r_support = 3.0 * pix_scale
        # this is true for this function
        max_conv_fn = 0.5
        # again, this seems extreme
        # cutoff_conv_fn = 0.25 * max_conv_fn
        cutoff_conv_fn = 0.0
        # in image.py we don't do this scaling for this kernel
        # which is clearly wrong as the beam is wider as a result of the
        # convolution.  This scale factor comes from Adam.
        # We should verify the appropriateness of this number by 
        # simulating a point source and a given beam size and then
        # measuring the source in the resulting image using this
        # kernel to determine the scaling.
        # It would be good to also know how this scales with the
        # parameters if we decide to expose those parameters to
        # the user.
        scale_fwhm = 1.09        # approximate
    else:
        # just makes the code easier to read, not actually used
        r_support = 1.0
        scale_fwhm = 1.0
        cutoff_conv_fn = 0.0

    # make 2 empty cubes: cube and weight
    # IDL routine set unset values to NaN, which are overwritten later
    # leaving truely unused pixels at NaN.  Need to sort out what's
    # appropriate here
    # keep these as 32-bit floats here, or as 64-bit floats that
    # get reduced to 32-bit at the end?
    data_cube = np.zeros(cubeShape,dtype=np.float32)
    weight_cube = np.zeros(cubeShape,dtype=np.float32)

    # some precalculations

    # calculate the pixel coordinates for each spectrum
    # the choice of velocity and stokes pixel shouldn't matter
    # use "0" as origin pixel here to match indexing convention in python
    zeros = np.zeros(len(xsky))
    x_pix, y_pix, v_pix = wcsObj.wcs_world2pix(xsky, ysky, zeros, 0)

    # the support radius squared
    r_support_sqrd = r_support**2

    # the support radius in pixels - at the center
    # assumes the pixel size is the same in both dimensions
    r_support_pix = r_support / pix_scale
    r_support_pix_sqrd = r_support_pix**2

    # the convolution function as a function of distance squared
    if kern != "nearest":
        n_pre = 10000
        pre_delta_sqrd = (1./n_pre)*r_support_sqrd
        pre_delta_pix_sqrd = pre_delta_sqrd / (pix_scale**2)
        pre_dist_sqrd = np.arange(n_pre*1.01)*pre_delta_sqrd
        pre_dist_sqrd_pix = pre_dist_sqrd / (pix_scale**2)
        pre_dist = np.sqrt(pre_dist_sqrd)
    
    if kern=="gauss":
        pre_conv_fn = max_conv_fn * np.exp(-0.5 * (pre_dist/(gauss_fwhm/2.354))**2)
    elif kern == "gaussbessel":
        # protect against pre_dist = 0.0 in the division
        x = scipy.pi * pre_dist/a
        zeroElem = np.where(x == 0.0)
        # at most, it can happen once
        if len(zeroElem) == 1:
            x[zeroElem[0]] += 0.00001
        pre_conv_fn = (scipy.special.j1(x)/x) * np.exp(-1.0 * (pre_dist/b)**2)
        if len(zeroElem) == 1:
            pre_conv_fn[zeroElem[0]] = max_conv_fn

    # the interpolator function where pre_conv_fn is used
    if kern != "nearest":
        x = np.arange(len(pre_conv_fn))
        interp_pre_conv_fn = interpolate.interp1d(x,pre_conv_fn)

        # DISTANCE INSIDE WHICH WE CAP THE CONVOLUTION FUNCTION
        cap_dist_sqrd = pre_dist_sqrd[1]
        cap_dist_sqrd_pix = pre_dist_sqrd_pix[1]

    # LOOP OVER THE SPATIAL AXES OF THE CUBE AND GRID

    # NOTE PRECOMPUTE ON THE R.A. TERM AND CONVOLUTION FUNCTION
    # *** ASSUMES THAT SPHERICAL TERMS ARE NEGLIGIBLE OVER THE KERNEL ***        
    
    nx = cubeShape[2]
    ny = cubeShape[1]

    unitSpectrum = np.ones((1,nchan))
    wtSpectrum = np.zeros(nchan)

    ##############################################################################
    if useC and C_available:
        if kern=="nearest":
            pre_delta_pix_sqrd = 0.
            pre_conv_fn = np.arange(2).astype(np.float64)
            cap_dist_sqrd_pix = 0.
            max_conv_fn = 0.
        
        ret = libobj.grid_otf_C(np.ravel(data).astype(np.float32),data.shape[0],nchan,
                                weight,np.asarray(x_pix,order='C'),np.asarray(y_pix,order='C'),
                                nx,ny,kern.encode('utf-8'),
                                r_support_pix_sqrd,pre_delta_pix_sqrd,pre_conv_fn, 
                                len(pre_conv_fn),cap_dist_sqrd_pix,max_conv_fn,
                                cutoff_conv_fn,threads)

        data_cube = np.ctypeslib.as_array(ret, shape=(nchan,ny,nx))
    #######################################################################
    else:
        counterMax = "%d" % nx
        counterFormat = "Row %%%dd out of %s" % (len(counterMax),counterMax)
        for i in range(nx):
            if verbose > 3:
                # update the counter string
                counterStr = counterFormat % (i+1)
                sys.stdout.write("\r%s" % counterStr)
                sys.stdout.flush()
            for j in range(ny):
          
                if kern == "nearest":
                    xdist = x_pix-i
                    ydist = y_pix-j
                    keep = (np.where((xdist>=-0.5) & (xdist<0.5) & (ydist>=-0.5) & (ydist<0.5)))[0]
                else:
                    pix_dist_sqrd = ((x_pix - i)*(x_pix - i) + (y_pix - j)*(y_pix - j))

                    #       EXTRACT PIXELS WITHIN THE RADIUS OF INTEREST OF THIS PIXEL 
                    keep = (np.where(pix_dist_sqrd <= r_support_pix_sqrd))[0]

                keep_ct = len(keep)

                if keep_ct > 1:
                    if kern == "nearest":
                        conv_fn = np.ones(keep_ct)
                    else:
                        #          INTERPOLATE TO GET THE CONV. FN. FOR EACH DATA POINT
                        pre_grid_x = pix_dist_sqrd[keep] / pre_delta_pix_sqrd

                        conv_fn = interp_pre_conv_fn(pre_grid_x)

                        #          DATA RIGHT ON TOP OF THE PIXEL TO THE MAX OF THE CONV. FN.
                        ind = (np.where(pix_dist_sqrd[keep] < cap_dist_sqrd_pix))[0]
                        if len(ind) > 0: 
                            conv_fn[ind] = max_conv_fn

                    #          PLACE A MINIMUM THRESHOLD NEEDED TO CONSIDER A GRID POINT
                    coverage = conv_fn.sum()
                    if coverage > cutoff_conv_fn:

                        #             GET A VECTOR OF NORMALIZED WEIGHTS
                        combined_weight = conv_fn * weight[keep]

                        # only go on if there's something there to grid
                        if combined_weight.sum() > 0.0:
                            
                            # reshape to 2D
                            combined_weight.shape = (len(combined_weight),1)
                            # dot product with unitSpectrum to get 2D of same shape as
                            # data to be gridded 
                            wtMatrix = combined_weight.dot(unitSpectrum)
                            # the chunk of data to be gridded here
                            spectrumMatrix = data[keep,:]
                            # the location of any NaNs
                            nanMask = np.isnan(spectrumMatrix)
                            # and set those locations to 0.0 in both the data and weights
                            wtMatrix[nanMask] = 0.0

                            # and again, only go on if there's something there to grid
                            if wtMatrix.sum() > 0.0:
                                spectrumMatrix[nanMask] = 0.0
                                # multiply by data by weights
                                spectrumMatrix *= wtMatrix
                                # sum and normalize spectrum to be inserted into data
                                wtSpectrum = wtMatrix.sum(axis=0)
                                spectrum = spectrumMatrix.sum(axis=0)

                                # the normalization should happen at the end if this is a re-entrant routine
                                # adding on to an existing image and weights.
                                spectrum /= wtSpectrum

                                # insert into the matrix
                                data_cube[:,j,i] = spectrum
                                weight_cube[:,j,i] = wtSpectrum

                #       HANDLE THE CASE OF ONLY ONE DATA POINT
                # this does the same as the above code, but with fewer steps since it's known
                # that there is only one spectrum contributing to this grid point
                #       ignore data with a weight <= 0.0
                if keep_ct == 1 and weight[keep[0]] > 0.0:
                    if kern == "nearest":
                        conv_fn = 1.0
                    else:
                        if pix_dist_sqrd[keep[0]] < cap_dist_sqrd_pix:
                            conv_fn = max_conv_fn 
                        else:
                            conv_fn = interp_pre_conv_fn(pix_dist_sqrd[keep[0]] / pre_delta_pix_sqrd)

                    if (conv_fn > cutoff_conv_fn):
                        thisSpectrum = data[keep[0],:]
                        wtSpectrum[:] = weight[keep[0]]
                        nanMask = np.isnan(thisSpectrum)
                        wtSpectrum[nanMask] = 0.0
                        if wtSpectrum.sum() > 0.0:
                            thisSpectrum[nanMask] = 0.0
                            data_cube[:,j,i] = thisSpectrum
                            weight_cube[:,j,i] = conv_fn*wtSpectrum

    # set 0.0 to NaN
    data_cube[data_cube==0.0] = float('nan')
    
    if verbose > 3:
        # finish the counter so output is back to it's usual place
        print('')

    # scale the beam size to reflect the convolution
    beam_fwhm = scale_fwhm*beam_fwhm

    return (data_cube, weight_cube, beam_fwhm)
