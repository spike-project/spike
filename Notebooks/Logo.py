# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # How to make a Logo
#
# ### first import good tools

# %%
import numpy as np
import matplotlib.pyplot as plt
import spike
import spike.NMR

# %% [markdown]
# ### import (or *simulate*) a nice time-signal
# for instance a 1 sec FID sampled at 1kHz - with a multiplet structure

# %%
# set-up scene
sw = 1000.0    # spectral width in Hz
w = 314.15926  # line position in Hz
J = 1.3        # coupling in Hz
tau = 0.9      # damping in second

# generate main frequency
t = np.arange(1000)/sw
fid = (np.cos(2*np.pi*w*t)+1j*np.sin(2*np.pi*w*t))*np.exp(-t/tau)

# add a pattern - here a heptuplet.
for i in range(6):
    fid *= np.cos(2*np.pi*J*t)

# load into spike
data = spike.NMR.NMRData(buffer=fid) 
data.axis1.specwidth = sw
data.set_unit('sec').display()

# %% [markdown]
# # And process
# Here, comparing a simple FFT to a apodisation with a shifted sine-bell

# %%
# compute
d1 = data.copy().zf(16).fft()
d2 = data.copy().apod_sin().zf(16).fft()

# %% [markdown]
# ### and draw the picture with nice settings

# %%
# define zoom and draw
z = (w+30, w-30)
darkblue = '#0A2B3D'
with plt.xkcd(length=150):   # add some effect
    d1.set_unit('Hz').display(zoom=z, linewidth=3, color='#3995C8')
    d2.set_unit('Hz').display(zoom=z, linewidth=3, color='orange', new_fig=False) #figure=d1.mplfigure) # other possibility
# add some text
fd = {'fontfamily': 'open sans',
      'color':  darkblue,
      'weight': 'semibold',
      'size': 110}
plt.text(w+28,15,'SP', fontdict=fd)
plt.text(w-5,15,'KE',  fontdict=fd)
# clean picture and change axes
plt.yticks([])
plt.ylabel('')
plt.xlabel('')
gca = plt.gca()
for side in ['top','left','right']:
    gca.spines[side].set_color('white')
gca.spines['bottom'].set_color(darkblue)
gca.xaxis.set_tick_params(color=darkblue, labelcolor=darkblue)
# and save
plt.savefig('Logo.png')


# %%
