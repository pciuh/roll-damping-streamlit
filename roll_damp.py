#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''

@date:   Fri, 03 Mar 2023, 17:32

@author: pciuh

This script estimates roll damping coefficient acc. to IKEDA
'''

import os
import sys
sys.path.append(os.getcwd()+'/modules/')
from ikeda import *

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import streamlit as st
#iDir = 'input/'
oDir = 'output/'
cMap = 'GnBu_r'


st.title('Roll Damping Coefficient by Ikeda')

fnam = 'ship'
LBL = ['LPP','B','T','CB','CM','KG','Tfi','lBK','hBK']
PLH = ['Length Between Perpendiculars (m)',
       'Breadth (m)',
       'Draught (m)',
       'Block Coefficient',
       'Midship Coefficient',
       'Vertical Centre of Gravity (m)',
       'Natural Roll Period (s)',
       'Bilge Keels Length (m)',
       'Bilge Keels Height (m)']

vec = []
col1,col2 = st.columns(2)
nl = int(len(LBL)/2)+1
with col1:
    for lbl,plh in zip(LBL[:nl],PLH[:nl]):
        vec.append(st.number_input(lbl,value=None,placeholder=plh,help=plh))

with col2:
    for lbl,plh in zip(LBL[nl:],PLH[nl:]):
        vec.append(st.number_input(lbl,value=None,placeholder=plh,help=plh))
##### Define Computations Domain
VS = np.arange(0,25,2)  #### Vector of speeds, knots
phi = np.arange(1,20,2) #### Vector of roll angles, degrees

##### Default Enviromental Conditions
rho,G = 1025.9,9.80665

if True in (x==None for x in vec):
    df = pd.read_csv('particulars/ship.csv',index_col=0)
    df = df.fillna(0)
    print(df)
    lpp,brth,t,cb,cm,KG,tw,lbk,bbk = df['Value'].values
else:
    lpp,brth,t,cb,cm,KG,tw,lbk,bbk = vec

ogd = 1-KG/t
b44,b44d = [],[]
for vs in VS:

    u = vs*1852/3600

    bbkd = bilke(cb,cm,ogd,brth,bbk,lbk,lpp,t,phi,tw)
    bfd  = fric(lpp,brth,t,cb,ogd,phi,tw)
    bwd  = wave(cb,cm,ogd,brth,t,phi,tw,u)
    bed  = eddy(cb,cm,ogd,brth,t,phi,tw)
    bld  = lift(lpp,brth,t,cb,cm,ogd,u)

    b44t = bfd + bwd + bed + bld + bbkd
    b44d = np.append(b44d,b44t)

    b44  = np.append(b44,b44t*rho*cb*brth**3*t*lpp*np.sqrt(2*G/brth))


b44d = b44d.reshape(-1,len(phi))
b44 = b44.reshape(-1,len(phi))

cnam = ['%4.1f knt'%x for x in VS]
inam = ['%4.1f deg'%x for x in phi]
dfa = pd.DataFrame( b44.T,columns=cnam,index=inam).round(1)
#dfa.to_csv(oDir + 'csv/' + fnam+'-b44.csv',sep=';')
dfd = pd.DataFrame(b44d.T,columns=cnam,index=inam).round(8)
#dfd.to_csv(oDir + 'csv/' + fnam+'-b44d.csv',sep=';')

phig,vg = np.meshgrid(phi,VS)

fig,ax = plt.subplots(figsize=(7,5))
cnt = ax.contourf(phig,vg,b44d,cmap=cMap,levels=11)
cbar = fig.colorbar(cnt,shrink=.75,orientation='horizontal',label=r'$\nu_{\varphi\varphi}$')
cbar.formatter.set_powerlimits((0, 0))
cbar.formatter.set_useMathText(True)
ax.set_ylabel(r'$V_S$ (knots)')
ax.set_xlabel(r'$\varphi$ (deg)')
fig.tight_layout()
fig.savefig(oDir + 'png/' + fnam+'-plt.png',dpi=300)
st.pyplot(fig)
st.subheader(r"Dimesionless Damping Coefficient $\nu_{\varphi\varphi}$", divider="gray")
st.dataframe(dfd.style.highlight_max(axis=0))

st.subheader(r"Damping Coefficient $B_{\varphi\varphi}$", divider="gray")
st.dataframe(dfa.style.highlight_max(axis=0))

#st.image(oDir + 'png/'+fnam+'-plt.png')
