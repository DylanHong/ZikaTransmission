import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns

# Set Parameters
sigh = 19
sigv = .15
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
def dSh(exposeh,Sh):
    val = -exposeh * Sh
    return val

def dEh(exposeh,Sh,infecth,Eh):
    val = exposeh * Sh - infecth * Eh
    return val

def dIh(infecth,Eh,recoverh,Ih):
    val = infecth*Eh - recoverh*Ih
    return val

def dRh(recoverh,Ih):
    val = recoverh * Ih
    return val

def dSv(birthv,Nv,exposev,Sv,deathv):
    val = birthv * Nv - exposev * Sv - deathv * Sv
    return val

def dEv(exposev,Sv,infectv,Ev,deathv):
    val = exposev * Sv - infectv * Ev - deathv * Ev
    return val

def dIv(infectv,Ev,deathv,Iv):
    val = infectv * Ev - deathv * Iv
    return val

def dexposeh(sigh,sigv,Nv,Nh,betahv,Iv):
    val = ((sigh*sigv*Nv)/(sigv*Nv + sigh*Nh))*betahv*(Iv/Nv)
    return val

def dexposev(sigh,sigv,Nv,Nh,betavh,Ih):
    val = ((sigh*sigv*Nv)/(sigv*Nv + sigh*Nh))*betavh*(Ih/Nh)
    return val


# big diff eq representing the whole equation
def model(y, t):
    exposeh = dexposeh(sigh, sigv, Nv, Nh, betahv, y[6])
    exposev = dexposev(sigh, sigv, Nv, Nh, betahv, y[2])
    sh = dSh(exposeh, y[0])
    eh = dEh(exposeh, y[0], infecth, y[1])
    ih = dIh(infecth, y[1], recoverh, y[2])
    rh = dRh(recoverh, y[2])
    sv = dSv(birthv, Nv, exposev, y[4], deathv)
    ev = dEv(exposev, y[4], infectv, y[5], deathv)
    iv = dIv(infectv, y[5], deathv, y[6])

    return [sh,eh,ih,rh,sv,ev,iv,exposeh,exposev]


# initial conditions
Sh0 = 1000000
Eh0 = 0
Ih0 = 0
Rh0 = 0
Sv0 = 10000000
Ev0 = 10000
Iv0 = 0
exposeh0 = dexposeh(sigh, sigv, Nv, Nh, betahv, Iv0)
exposev0 = dexposev(sigh, sigv, Nv, Nh, betahv, Ih0)

y0 = [Sh0, Eh0, Ih0, Rh0, Sv0, Ev0, Iv0, exposeh0, exposev0]


# time points
t = np.linspace(0,200,num=100000)

# solve ODE
y = odeint(model, y0, t)

infectedH = y[:, 2]
infectedV = y[:, 6]

sns.set()

# plot results
fig, ax = plt.subplots(2)
ax[0].plot(t, infectedH)
ax[1].plot(t, infectedV)
plt.show()
