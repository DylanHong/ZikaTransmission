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


# Define Differential Equations
def dSh(exposeh,Sh):
    val = -exposeh * Sh
    return val

def dEh(exposeh,Sh,infecth,Eh):
    val = (exposeh * Sh) - (infecth * Eh)
    return val

def dIh(infecth,Eh,recoverh,Ih):
    val = (infecth * Eh) - (recoverh * Ih)
    return val

def dRh(recoverh,Ih):
    val = (recoverh * Ih)
    return val

def dSv(birthv,Nv,exposev,Sv,deathv):
    val = (birthv * Nv) - (exposev * Sv) - (deathv * Sv)
    return val

def dEv(exposev,Sv,infectv,Ev,deathv):
    val = (exposev * Sv) - (infectv * Ev) - (deathv * Ev)
    return val

def dIv(infectv,Ev,deathv,Iv):
    val = (infectv * Ev) - (deathv * Iv)
    return val

def dexposeh(sigh,sigv,Nv,Nh,betahv,Iv):
    val = ((sigh*sigv*Nv)/(sigv*Nv + sigh*Nh))*betahv*(Iv/Nv)
    return val

def dexposev(sigh,sigv,Nv,Nh,betavh,Ih):
    val = ((sigh*sigv*Nv)/(sigv*Nv + sigh*Nh))*betavh*(Ih/Nh)
    return val


# big diff eq representing the whole equation
def model(y, t):
    Nh = y[0]+y[1]+y[2]+y[3]
    Nv = y[4]+y[5]+y[6]
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
Sh0 = 15276566
Eh0 = 0
Ih0 = 0
Rh0 = 0
Sv0 = 200000000
Ev0 = 1000
Iv0 = 0

Nh0 = Sh0
Nv0 = Sv0 + Ev0

exposeh0 = dexposeh(sigh, sigv, Nv0, Nh0, betahv, Iv0)
exposev0 = dexposev(sigh, sigv, Nv0, Nh0, betahv, Ih0)

y0 = [Sh0, Eh0, Ih0, Rh0, Sv0, Ev0, Iv0, exposeh0, exposev0]


# time points
t = np.linspace(0,1000,num=100000)

# solve ODE
y = odeint(model, y0, t)

infectedH = y[:, 2]
infectedV = y[:, 6]
exposedH = y[:, 1]
susceptibleH = y[:, 0]
recoveredH = y[:, 3]

sns.set()

maxes = []
x = np.linspace(0,.5,200)
for i in x:
    # solve ODE
    sigv = i
    y = odeint(model, y0, t)
    maxes.append(max(y[:, 2]))

plt.figure()
plt.plot(x,maxes)
plt.title('Effect of Bite Rate on Severity of Outbreak')
plt.xlabel('Bite Rate')
plt.ylabel('Maximum Number of Infected Humans')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.savefig('maxes.png')

maxes1 = []
x = np.linspace(10000000,250000000,200)
for i in x:
    # solve ODE
    Sv0 = i
    y0 = [Sh0, Eh0, Ih0, Rh0, Sv0, Ev0, Iv0, exposeh0, exposev0]
    y = odeint(model, y0, t)
    maxes1.append(max(y[:, 2]))

plt.figure()
plt.plot(x,maxes1)
plt.title('Effect of Bite Rate on Severity of Outbreak')
plt.xlabel('Bite Rate')
plt.ylabel('Maximum Number of Infected Humans')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.savefig('maxes1.png')


sums = []
x = np.linspace(0,.5,200)
for i in x:
    # solve ODE
    sigv = i
    y = odeint(model, y0, t)
    sums.append(sum(y[:, 2]))

plt.figure()
plt.plot(x,sums)
plt.title('Effect of Bite Rate on Severity of Outbreak')
plt.xlabel('Bite Rate')
plt.ylabel('Maximum Number of Infected Humans')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.savefig('sums.png')


duration = []
vals = np.linspace(0.05,.5,200)
for i in vals:
    # solve ODE
    sigv = i
    y = odeint(model, y0, t)
    foundStart = False
    start = None
    end = None
    for j, value in enumerate(y[:, 2]):
        if foundStart is False and value >= 1000:
            start = j
            foundStart = True
        elif foundStart is True and value <= 1000:
            end = j
            break

    dur = t[end] - t[start]
    duration.append(dur)

plt.figure()
plt.plot(vals,duration)
plt.title('Effect of Bite Rate on Duration of Outbreak')
plt.ylabel('Duration (Days)')
plt.xlabel('Bite Rate')
plt.savefig('durations.png')

duration1 = []
vals = np.linspace(10000000,250000000,200)

for i in vals:
    # solve ODE
    Sv0 = i
    y0 = [Sh0, Eh0, Ih0, Rh0, Sv0, Ev0, Iv0, exposeh0, exposev0]
    y = odeint(model, y0, t)
    foundStart = False
    start = None
    end = None
    for j, value in enumerate(y[:, 2]):
        if foundStart is False and value >= 1000:
            start = j
            foundStart = True
        elif foundStart is True and value <= 1000:
            end = j
            break

    dur = t[end] - t[start]
    duration1.append(dur)

plt.figure()
plt.plot(vals,duration1)
plt.title('Effect of Bite Rate on Duration of Outbreak')
plt.ylabel('Duration (Days)')
plt.xlabel('Bite Rate')
plt.savefig('durations1.png')

# time points
t = np.linspace(0,100,num=10000)

# solve ODE
y = odeint(model, y0, t)

infectedH = y[:, 2]
infectedV = y[:, 6]
exposedH = y[:, 1]
susceptibleH = y[:, 0]
recoveredH = y[:, 3]

# plot results
plt.figure()
plt.plot(t, infectedH)
plt.plot(t, exposedH)
plt.plot(t, susceptibleH)
plt.plot(t, recoveredH)
plt.legend(['Infected', 'Exposed', 'Susceptible', 'Recovered'])
plt.title('Zika Infection Dynamics over Time')
plt.ylabel('Number (in each category)')
plt.xlabel('Days')
plt.savefig('overall.png')
