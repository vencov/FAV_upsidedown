#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 28 19:10:26 2021

@author: vencov
"""

from scipy.io.wavfile import read

import matplotlib.pyplot as plt
import numpy as np

plt.close('all') # close figures
fs, soae_signal = read('s015_soae_prave.wav')
#fs, soae_signal = read('s015_soae_prave.wav')
#fs, soae_signal = read('s013_soae_leve.wav')

def rfft(x):
    if x.ndim==1:  # check whether x is 1D
        N = int(np.ceil(len(x)/2)) # half of the spectrum
        xc = np.fft.fft(x)
        X = 2*xc[0:N]/len(x)  # take first half of the spectrum and normalize
    elif x.ndim==2: # for 2D array, do it column-wise
        m,n = x.shape
        N = np.ceil(m/2)
        xc = np.fft.fft(x)
        X = 2*xc[0:N,:]/len(x)  # take 
    else:
        print('Input must be either 1 or 2D array')
        X = None
    
    return X


def irfft(X,nfft=None):
    if X.ndim == 1:
        nf = len(X)
        X[0] = np.real(X[0])
        X[nf-1] = np.real(X[nf-1])
        if nf%2:  # liche
            X = np.concatenate((X, np.conj(X[:0:-1])))  # conjugate and flip except the first sample X[0]
        else:   # sude
            X = np.concatenate((X, [np.conj(X[0])], np.conj(X[:0:-1])))
    elif X.ndim == 2: # assume that spectrum is along the column
        nf,_ = X.shape
        X[0,:] = np.real(X[0,:])
        X[nf-1,:] = np.real(X[nf-1,:])
        X = np.vstack((X, np.conj(X[1:-1,:])))
    else:
        print('Input must be either 1 or 2D array')
        return None
    
    x = np.fft.ifft(X,nfft)
    
    x = np.real(x) + np.imag(x)
    return len(x)/2*x


'''
    else:
        np.fft
  [k,l]=size(X);
  if (k==1)
    nf=l;
    % zero out Imag part in DC and Fmax
    X(:,1)=real(X(:,1));
    X(:,nf)=real(X(:,nf));
    X=[X,conj(X(:,nf-1:-1:2))];
  else
    nf=k;
    % zero out Imag part in DC and Fmax
    X(1,:)=real(X(1,:));
    X(nf,:)=real(X(nf,:));
    X=[X;conj(X(nf-1:-1:2,:))];
  end


  if (nargin == 2)
    x = ifft(X,nfft);
  else
    x = ifft(X);
  end
  x = real(x) + imag(x);

  [m,n] = size(x);
  if m == 1
    m = n;
  end
  x = (m/2) * x;

'''

N = len(soae_signal)

Nwin = 8192  # time frame for spectrum calculation

# average magnitude spectra in Nwin frames
for i in range(N//Nwin):
    if i==0:
        #AvgSp = np.abs(np.fft.fft(soae_signal[i*Nwin:(i+1)*Nwin]))/(Nwin*N//Nwin)
        AvgSp = np.abs(rfft(soae_signal[i*Nwin:(i+1)*Nwin]))/(N//Nwin)
    else:
        #AvgSp = AvgSp + np.abs(np.fft.fft(soae_signal[i*Nwin:(i+1)*Nwin]))/(Nwin*N//Nwin)
        AvgSp = AvgSp + np.abs(rfft(soae_signal[i*Nwin:(i+1)*Nwin]))/(N//Nwin)



    
fxw = np.arange(0,Nwin)*fs/Nwin # frequency axis
fxw =fxw.flatten()


fx = np.arange(0,N)*fs/N
fx =fx.flatten()

SpAll = rfft(soae_signal)

#fig,ax = plt.subplots()
#ax.plot(fx[:len(SpAll)],20*np.log10(np.abs(SpAll)/(np.sqrt(2)*2e-5)))

soae_inv = irfft(SpAll)
#fig,ax = plt.subplots()
#ax.plot(soae_signal)c
#ax.plot(irfft(SpAll))


# extract one of the peaks in the data and perform analysis
# for extraction, use filtration in the frequency domain


def roex_in_freq(fx,fc,N=10):
    g0 = 1
    for i in range(2,N+1):
        g0 = np.log(g0+1)
    T = np.sqrt(g0)*fx/fc
    G = np.exp(T**2)
    for i in range(2,N+1):
        G = np.exp(G-1)
    G = 1/G
    idx = np.where(fx>=1.2*fc)
    G = np.concatenate((G[idx[0][0]-1:0:-1],G[:idx[0][0]]))
    
    
    return G

# create a window



#fig,ax = plt.subplots()
#ax.plot(roex_in_freq(fx,1012,10))

CF = 1740
idx = np.where(fx>=CF)
BW = 40
G = roex_in_freq(fx,BW/2,10)
# add zeros between 0 and begining of window

nulypred = np.zeros(idx[0][0]-int(len(G)/2))
nulyvzad = np.zeros(len(fx)-len(nulypred)-len(G))
okno = np.concatenate((nulypred,G,nulyvzad))


#fig,ax = plt.subplots()

#ax.plot(fx,okno)
#ax.plot(fx[:len(SpAll)],np.abs(SpAll)*1e5)



fig, ax = plt.subplots()

ax.plot(fxw[:len(AvgSp)],20*np.log10(AvgSp/(np.sqrt(2)*2e-5)))
ax.plot(fx,20*np.log10(20*okno))
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel('Amplitude of sound pressure (dB re 20 $\mu$Pa)')
plt.xlim((100, 16000))
plt.ylim((-30, 30))

# nasob spektrum okenkem

vysek = SpAll*okno[:len(SpAll)]

#fig,ax = plt.subplots()
#ax.plot(irfft(vysek))

tdsignal = irfft(vysek)
# obalka


from scipy.signal import hilbert

hv = hilbert(tdsignal)

tx = np.arange(0,len(tdsignal)/fs,1/fs)
               
fig,ax = plt.subplots()
ax.plot(tx,tdsignal,'b',lw=0.5,alpha=0.6)
ax.plot(tx,abs(hv),'k',tx,-abs(hv),'k',alpha=1)
ax.set_xlabel('time (seconds)')
ax.set_ylabel('amplitude (Pa)')



bins = np.linspace(-5e-4,5e-4,1000)
h2d = np.histogram2d(np.real(hv),np.imag(hv),bins)



fig, ax = plt.subplots()
ax.imshow(h2d[0])
#X, Y = np.meshgrid(h2d[1], h2d[2])
#fix,ax = plt.subplots(subplot_kw={"projection": "3d"})
#ax.plot_surface(X,Y,h2d[0])
#surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
 #                      linewidth=0, antialiased=False)



