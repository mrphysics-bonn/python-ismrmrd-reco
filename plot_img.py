#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot images from FIRE Reco
"""

import matplotlib.pyplot as plt
import numpy as np
import h5py
import os

show_set = 0 # number of image set in reco
show_slice = 0 # number of slice shown in single view

recodir = os.getcwd() + "/recon/"
file = "out.h5"
data = h5py.File(recodir+file, 'r+')

imgs = []
for key in data['images']:
    if key=='xml':
        continue
    imgs.append(data['images'][key]['data'][:])
    
img_set = imgs[show_set][:]
if np.iscomplexobj(img_set):
    img_set = abs(img_set)

slices = img_set.shape[0]
if slices > 1:
    columns = int(np.sqrt(slices))
    rows = int(slices/columns+1)
    maxsig = max(abs(img_set.flatten()))
    fig=plt.figure(figsize=(12, 12))
    for i in range(slices):
        fig.add_subplot(rows, columns, i+1)
        plt.imshow(img_set[i,0,0],cmap='gray', vmin=0, vmax=maxsig, interpolation='none')
        plt.axis('off')
        plt.tight_layout(pad=0.2)


plt.figure(figsize=(10,10))
plt.imshow(img_set[show_slice,0,0], 'gray', interpolation='none')
plt.axis('off')

data.close()
