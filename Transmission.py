import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns

# Set Parameters
sigh = 19
sigv = 1
betahv = .7
betavh = .7
infecth = 1/5.5
infectv = 1/8.2
recoverh = 1/6
birthv = 1/14
deathv = 1/14
Nh = 15276566
Nv = 200000000

# Define Differential Equations

def Sh(exposeh,Sh):
    val = -exposeh * Sh
    return val

def Eh(exposeh,Sh,infecth,Eh):
    val = exposeh * Sh - infecth * Eh
    return val

def Ih(infecth,Eh,recoverh,Ih):
    val = infecth*Eh - recoverh*Ih
    return val

def Rh(recoverh,Ih):
    val = recoverh * Ih
    return val

def Sv(birthv,Nv,exposev,Sv,deathv):
    val = birthv * Nv - exposev * Sv - deathv * Sv
    return val

def Ev(exposev,Sv,infectv,Ev,deathv):
    val = exposev * Sv - infectv * Ev - deathv * Ev
    return val

def Iv(infectv,Ev,deathv,Iv):
    val = infectv * Ev - deathv * Iv
    return val

def exposeh(sigh,sigv,Nv,Nh,betahv,Iv):
    val = ((sigh*sigv*Nv)/(sigv*Nv + sigh*Nh))*betahv*(Iv/Nv)
    return val

def exposev(sigh,sigv,Nv,Nh,betavh,Iv):
    val = ((sigh*sigv*Nv)/(sigv*Nv + sigh*Nh))*betavh*(Ih/Nh)
