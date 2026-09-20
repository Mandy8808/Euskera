"""Spectral plotting implementation."""

import numpy as np
import matplotlib.pyplot as plt

def plotImag(r_dis2, Sigma0, Ueff):
   fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 4),
                       sharex=False, sharey=False,
                       gridspec_kw=dict(hspace=0.0, wspace=.28))

   ax[0].plot(r_dis2, np.diagonal(Sigma0), 'o', markersize=2, mfc='white')
   ax[1].plot(r_dis2, np.diagonal(Ueff), 'o', markersize=2, mfc='white')
   #ax[0].set_xlim(0, 30)
   plt.show()
   return
