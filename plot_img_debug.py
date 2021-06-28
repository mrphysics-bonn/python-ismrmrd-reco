"""
Created on Mon Sep 28 16:46:04 2020

@author: ehsesp
"""

import matplotlib.pyplot as plt
import numpy as np
import os


os.chdir('debug/')

raw = np.load('raw.npy')
img = np.load('img.npy')
sensmaps = np.load('sensmaps.npy')
acs = np.load('acs.npy')
acs_img = np.load('acs_img.npy')
trj = np.load('trj.npy')


if img.ndim>2:
    img = img[:,:,img.shape[-1]//2]

if sensmaps.ndim == 4:
    columns = 4
    rows = 4
    maxsig = max(abs(sensmaps.flatten()))
    fig=plt.figure(figsize=(12, 6))
    for i in range(columns*rows):
        fig.add_subplot(rows, columns, i+1)
        plt.imshow(abs(sensmaps[:,:,0,int(i)]))
        plt.axis('off')
        plt.tight_layout(pad=0.2)
    
    fig=plt.figure(figsize=(12, 6))
    for i in range(columns*rows):
        fig.add_subplot(rows, columns, i+1)
        plt.imshow(np.angle(sensmaps[:,:,0,int(i)]))
        plt.axis('off')
        plt.tight_layout()
else:
    fig=plt.figure(figsize=(12, 6))
    plt.imshow(abs(sensmaps))

plt.figure(figsize=(10,10))
plt.imshow(img, 'gray', interpolation='none')
plt.axis('off')

acs_fft = acs.copy()
for dim in range(0,3):
    acs_fft = np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(acs_fft, axes=[dim]),axis=dim),axes=[dim])
acs_fft = np.sqrt(np.sum(np.abs(acs_fft)**2, -1))
acs_fft = acs_fft[0]

plt.figure(figsize=(12,6))
plt.imshow(acs_fft, 'gray', interpolation='none')
plt.axis('off')

if acs_img.ndim == 4:
    acs_img = np.sqrt(np.sum(np.abs(acs_img)**2, -1))
plt.figure(figsize=(12,6))
plt.imshow(abs(acs_img), 'gray', interpolation='none')
plt.axis('off')