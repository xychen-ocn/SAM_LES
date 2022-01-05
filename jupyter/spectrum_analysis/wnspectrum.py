#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module contains functions to compute cross-spectrum (co, quadrature) for
spatial (x,y) data in LES. Reference from Maria's spacetime package.

List of function:
    cross_spectrum = cospectrum + i quadrature spectrum

Created on Tue Jan  4 14:33:25 2022

@author: X.Y.Chen

"""


import numpy as np
import scipy.interpolate as intp
import matplotlib.pyplot as plt

def cross_spectrum(XX, YY, dx, dy, opt= False):
    """
    Compute the FFT to get the power and cross-spectra for one time segment.
    :param XX: Input array (time, lat, lon)
    :param YY: Input array (time, lat, lon)
    :param opt: Optional parameter, not currently used. Set to False.
    :return STC: Spectra array of shape (8, nfreq, nwave). Last 4 entries are blank and need to be computed by calling
    mjo_cross_coh2pha. The first 4 entries contain power spectra for XX, power spectra for YY, co-spectra between XX
    and YY, quadrature spectra between XX and YY.
    """

    
    NT, NM, NL = XX.shape
    # compute fourier decomposition in x(NL) and y(NM) direciton 
    Xfft = np.fft.fft2(XX, axes=(1, 2))
    Yfft = np.fft.fft2(YY, axes=(1, 2))
    # normalize by # time samples
    #Xfft = Xfft / (NM * NL)
    #Yfft = Yfft / (NM * NL)
    # shift 0 frequency and 0 wavenumber to the center
    Xfft = np.fft.fftshift(Xfft, axes=(1, 2))
    Yfft = np.fft.fftshift(Yfft, axes=(1, 2))
    
    
    ############################

    # average the power spectra across all latitudes (I don't need averaging)
    # so PX and PY has three dimensions. [time, y, x]
    PX = np.square(np.abs(Xfft))  
    PY = np.square(np.abs(Yfft))

    # compute co- and quadrature spectrum
    #PXY = np.average(np.conj(Yfft) * Xfft, axis=1)  # latitudinal averaged again.
    PXY = np.conj(Yfft) * Xfft
    CXY = np.real(PXY)
    QXY = np.imag(PXY)

    # is reversing the longitudinal dim. a requirement?   (it is reversing negative and positive frequency order on the axis) 
    PX = PX[:, ::-1]
    PY = PY[:, ::-1]
    CXY = CXY[:, ::-1]
    QXY = QXY[:, ::-1]

    # test if latitude[NM] and longitude[NL] are odd or even, fft algorithm
    # returns the Nyquist frequency once for even NT or NL and twice
    # if they are odd; (XYC modified the shifting a little bit, so that the Power spectrum
    # is symmetric relative to the origin.
    if NM % 2 == 1:    # odd
        nky = NM
        if NL % 2 == 1:
            nkx = NL
            STC = np.zeros([8, NT, nky, nkx], dtype='double')
            STC[0, :NT, :NM, :NL] = PX
            STC[1, :NT, :NM, :NL] = PY
            STC[2, :NT, :NM, :NL] = CXY
            STC[3, :NT, :NM, :NL] = QXY
        
        else:         # even
            nkx = NL + 1
            STC = np.zeros([8, NT, nky, nkx], dtype='double')
            STC[0, :NT, :NM, :NL ] = PX
            STC[1, :NT, :NM, :NL ] = PY
            STC[2, :NT, :NM, :NL ] = CXY
            STC[3, :NT, :NM, :NL ] = QXY
            STC[:, :, :, NL] = STC[:, :, :, 0]
    else:
        nky = NM + 1
        if NL % 2 == 1:
            nkx = NL
            STC = np.zeros([8, NT, nky, nkx], dtype='double')
            STC[0, :NT, 1:NM+1, :NL] = PX
            STC[1, :NT, 1:NM+1, :NL] = PY
            STC[2, :NT, 1:NM+1, :NL] = CXY
            STC[3, :NT, 1:NM+1, :NL] = QXY
            STC[:, :, 0, :] = STC[:, :, NM, :]
        else:
            nkx = NL + 1
            STC = np.zeros([8, NT, nky, nkx], dtype='double')
            STC[0, :NT, 1:NM+1, :NL ] = PX
            STC[1, :NT, 1:NM+1, :NL ] = PY
            STC[2, :NT, 1:NM+1, :NL ] = CXY
            STC[3, :NT, 1:NM+1, :NL ] = QXY
            STC[:, :, 0, :] = STC[:, :, NM, :]   # last element = first entry 
            STC[:, :, :, NL] = STC[:, :, :, 0]
            
    # define wavenumber space:
    ikx = np.arange(-int(nkx/ 2), int(nkx / 2)+ 1, 1.)
    iky = np.arange(-int(nky/ 2), int(nky / 2)+1, 1.)

    kx = 2*np.pi/(NL*dx) * ikx   # the last kx is the nyquist frequency (k at the edge)
    ky = 2*np.pi/(NM*dy) * iky 
        
        
    
    # compute the coherence and phase related parameters:
    cross_coh2pha(STC)

    return {'STC': STC, 'kx': kx, 'ky': ky}


def cross_coh2pha(STC):
    """
    Compute coherence squared and phase spectrum from averaged power and
    cross-spectral estimates.
    :param STC: Spectra array.
    :return STC: Spectra array of the same size with entries 4-7 (coherence squared, phase angle, phase component 1,
    phase component 2) recomputed based on the power and cross-spectra in entries 0-3.
    
    author: M. Gehne
    """

    nvar, ntime, nkx, nky = STC.shape

    PX = STC[0, :, :, :]
    PY = STC[1, :, :, :]
    CXY = STC[2, :, :, :]
    QXY = STC[3, :, :, :]

    PY[PY == 0] = np.nan
    COH2 = (np.square(CXY) + np.square(QXY)) / (PX * PY)   # PXY/(PX * PY)
    PHAS = np.arctan2(QXY, CXY)

    # why V1 is negative QXY?
    V1 = -QXY / np.sqrt(np.square(QXY) + np.square(CXY))   # sine
    V2 = CXY / np.sqrt(np.square(QXY) + np.square(CXY))    # cosine

    STC[4, :, :] = COH2
    STC[5, :, :] = PHAS
    STC[6, :, :] = V1
    STC[7, :, :] = V2

    return (STC)


def omnidirectional_wavenumber_spectrum(Sin, kx, ky, dx, dy, checkflag=True):
    """
    obtain the omni-directional wavenumber spectrum from wavenumber-wavenumber spectrum
    Steps:  kx-ky spectrum --> K-theta spectrum (wavenumber-direction spectrum)
            --> integrate in directional space 
    Inputs: S --> a kx-ky spectrum 
    author: X.Y. Chen
    """
    
    NT, NY, NX = np.shape(Sin)
    KX, KY = np.meshgrid(kx,ky)
    # construct k-theta space:
    K = np.sqrt(KX**2 + KY**2)
    THETA = np.arctan2(KY,KX)*180/np.pi
    
    # bin kx-ky spectrum according to wawvenumber and direction:
    dx_nyquist = 2* (np.sqrt(dx**2 + dy**2))

    # ###### Work Zone #######
    k_nyq = 2*np.pi/dx_nyquist
    Lmax = np.sqrt((NX*dx)**2 + (NY*dy)**2)
    k0 =  2*np.pi/Lmax           # minimum wave number (longest wavelength)
    
    #exp_ind = (np.log10(k0)):0.05:(np.log10(k_nyq))
    exp_ind = np.arange(np.log10(k0), np.log10(k_nyq)+0.05, 0.05)
    wn_edges = 10**exp_ind
    
    Indx_kbined = np.digitize(K, wn_edges)

    DirEdge = np.linspace(-180,180,num=37)
    #Indx_thbined = np.digitize(THETA, DirEdge)
    
    wn_bins = 0.5*(wn_edges[:-1] + wn_edges[1:])
    theta_bins = 0.5*(DirEdge[:-1] + DirEdge[1:])
    
    nk = len(wn_bins)
    nth = len(theta_bins)
    S_kth = np.zeros([nk, nth])
    PK = np.zeros([NT, nk])
    
    # initialize some record holding np arrays:
    #rec_k = np.zeros(nk)
    #rec_th = np.zeros([nk,nth])
    for it in range(NT):
        S = Sin[it,:,:]
        
        for ik in np.linspace(1,nk,num=nk, dtype=np.int32):
            ids = np.where(Indx_kbined==ik)[0]
            
            if np.size(ids)>0:     # not empty
                #rec_k[ik] = len(ids)
                Sk_tmp = S[ids]    # S has two dimensions [NY, NX] same as K
                theta_tmp = THETA[ids]
                #K_tmp = K[ids]
                
                
                Ysub = np.digitize(theta_tmp, DirEdge)
                    
                for ith in np.linspace( 1, nth, num=nth, dtype=np.int32):
                    idth = np.where(Ysub==ith)[0]
                        
                    if np.size(idth)>0:
                        #rec_th[ik,ith] = len(idth)
                        # bin averaged in direcitonal bins
                        S_kth[ik-1, ith-1] = np.nanmean(Sk_tmp[idth])
                    else:
                            
                        S_kth[ik-1, ith-1] = np.nan
        
        # % ---- build a regular grid:
        HH, KK = np.meshgrid(theta_bins, wn_bins)   # in this order, HH & KK will have dim: [nk, nth]
        
        mask = ~np.isnan(S_kth)    # eqiv. to ~isnan in matlab
        theta_bins_valid = HH[mask]
        wnbins_valid = KK[mask]
        S_kth_valid = S_kth[mask]

        gridin = np.zeros([len(theta_bins_valid),2])
        gridin[:,0] = theta_bins_valid
        gridin[:,1] = wnbins_valid
        
        pnts = np.zeros([np.size(HH),2])
        pnts[:,0] = np.reshape(HH,[np.size(HH)])
        pnts[:,1] = np.reshape(KK,[np.size(KK)])
        S_kth_filled = intp.griddata(gridin, S_kth_valid, pnts, method='cubic')
        #print(np.shape(S_kth_filled))
        
        S_kth_filled = np.reshape(S_kth_filled, [nk, nth])
       
        # % 3. integrate the 2D spectrum in the directional space:
        # %    using trapezoidal method:
        Sk= np.trapz(S_kth_filled, theta_bins, axis =1)
        
        # % do another interpolation:
        mask0 = (Sk>np.nanmax(Sk)*0.001)
        Sk_valid = Sk[mask0]
        wnbins_valid = wn_bins[mask0]
        F = intp.interp1d(wnbins_valid, Sk_valid, kind='linear', fill_value='extrapolate')
        Sk_filled = F(wn_bins)
        
        
        # check results:
        if checkflag:
            if it%10==0:
                wvlen  = 2*np.pi/wn_bins
                fig, ax = plt.subplots(1,2, figsize=(12, 8))
                pm =ax[0].pcolormesh(theta_bins, wn_bins, S_kth_filled,shading ='auto')
                plt.colorbar(pm, ax=ax[0])
                plt.yscale('log')
                
                ax[1].plot(wvlen, Sk_filled, '-r')
                ax[1].plot(wvlen, Sk, '--k')
                plt.gca().invert_xaxis()
                plt.xscale('log')
                plt.yscale('log')
                fig.tight_layout()
        
        
    PK[it,:] =  Sk_filled
    
    return {'PK': PK, 'wvnum': wn_bins, 'wvlen': 2*np.pi/wn_bins}


def compute_Ogive_Length(Pk, k):
    
    return L_Ogive
    