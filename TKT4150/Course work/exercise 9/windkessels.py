
# coding: utf-8

# # TKT 4150 Biomechanics: Problem set 9
# 
# ### Introduction
# This is an IPython notebook, for a more general introduction please see this [introduction](https://nbviewer.jupyter.org/github/jupyter/notebook/blob/master/docs/source/examples/Notebook/Notebook%20Basics.ipynb) and [guide](https://nbviewer.jupyter.org/github/jupyter/notebook/blob/master/docs/source/examples/Notebook/Running%20Code.ipynb) where you can download a similar notebooks that are more instructional. This website gives a similar [guide](https://jupyter-notebook-beginner-guide.readthedocs.io/en/latest/). 
# 
# ### Save your work!
# * If you are running this on [tmpnb.org](http://tmpnb.org), remember to save your work by downloading the notebook via clicking `File->Download As` in the toolbar.
# 
# 
# 
# ## Exercise 1 Deriving the Windkessel Model
# ### Part A Derivation from 1st principles of conservation of mass
# Derive the 2 element Windkessel by assuming a compliant tube with inflow Q_1 such the $A(P) = C P$ where the distal (downstream) flow has a fluid resistance $P=RQ_2$ where $Q_2$ is the flow out of the tube.
# 
# ### Part B Derivation of 4-element Windkessel Model via averaging of the 1-D Navier Stokes
# Consider the governing equations for flow through a compliant tube
# $$ \frac{\partial A}{\partial t} + \frac{\partial Q}{\partial x} = 0$$
# $$ \frac{\partial Q}{\partial t} + \frac{\partial}{\partial x}(\frac{Q^2}{A}) + \frac{A}{\rho} \frac{\partial P}{\partial x} + K_R \frac{Q}{A} = 0$$
# completed with the constitutive equation $P=\beta(\sqrt{A}-\sqrt{A_0})$
# 
# from $x=0$ to $x=l$ the average values of $P$, $A$ and $Q$, i.e. $\bar{P} = \frac{1}{l}\int_0^l P dx$ may be related as 
# 
# $$ k_1 l \frac{d \bar{P}}{dt} = Q_1 - Q_2$$ 
# and
# 
# $$ \frac{d\bar{Q}}{dt}\frac{\rho l}{A_0} + \frac{K_{R} l}{{A_0}} \bar{Q} = P(0) - P(l) $$ 
# where $K_R$ is the coefficient of viscous friction for a given velocity profile and $A_0$ is the characteristic cross sectional area.
# 
# Assume the convective term $\frac{\partial}{\partial x} (\frac{Q^2}{A})$ is negligible and that spatial variations in $A$ are much smoother than variations in $P$ and $Q$, i.e. treat $A$ as a constant over space.
# 
# ## Exercise 2 Frequency domain impedance analysis comparison of Windkessel models
# Windkessel models are often used to interpret pressure and flow signals in the frequence domain. In this domain the relationship between pressure and flow is given as $\hat{p} = Z \hat{q}$ where $\hat{p}$ and $\hat{q}$ are the fourier transform of $p$ and $q$.
# Impedances are linear and thus may be added accordning to the formula $Z_{eq} = \sum Z_{i}$ for impedances in series, i.e. same flow through a sequence of elements, and $Z_{eq}^{-1} = \sum Z_{i}^{-1}$ for impedances in parallel, i.e. same pressure across a sequence of elements by split of flow. Capacitance has an impedance $Z_{C} = \frac{1}{j \omega C}$, resistance has impedance $Z_{R} = R$, and inertance has impedance $Z_{L} = j \omega L$ where $\omega= 2\pi f$ is the angular frequency.
# 
# 
# ![Windkessel Schematics](./WindkesselSchematics.jpg)
# 
# ### Derive the impedance for the 2-Element, 3-Element and 4-Element Windkessel models
# Using the above relations find equivalent impedances for each of the schematics. Compare the input impedances over frequencies of inputs. What interpretation can you make about the impedances at low and high frequencies?
# 
# ### Compare these relations
# Below you will find implementation of the functions necessary to fit the impedance of the 2 element Windkessel model to a data set, e.g. './55arteryNetwork_age19.csv' and './55arteryNetwork_age75.csv'. If you prefer to use `Matlab` the attached files have code to get you started. The same python code is also available as a python script.
# 
# Use this to compare compliance, characteristic aortic impediance, and characteristic aortic intertance between the 19 year old and the 75 year old. 

# In[1]:

get_ipython().magic('matplotlib notebook')
import numpy as np
import scipy as sp
import scipy.optimize 
from scipy import io as sio
from matplotlib import pyplot as plt


# In[2]:

unit_Pa_to_mmHg = 1/133.32
unit_m3_to_ml = 1e6


# In[3]:

# This code may be used to import ./55arteryNetwork_age19.csv or ./55arteryNetwork_age75.csv
csv = np.genfromtxt ('./55arteryNetwork_age19.csv', delimiter=";", names=True)
p = csv["P_Pa"]*unit_Pa_to_mmHg
q = csv["Q_m3s"]*unit_m3_to_ml
t = csv["t_s"]
dt = dt = t[1]-t[0]
N = len(t)
T = N*dt
w = np.arange(N//2+1)*2*np.pi/T # Frequencies corresponding to an FFT of the data.


# In[4]:

def ZWK2(R,C,w):
    """ Calculate and return impedance of a 2elt WK 
    Args:
        R Resistance
        C Compliance
        w the frequencies at which the impediance is calculated.
    """
    Z=R/(1+1j*w*R*C)
    return Z


# In[5]:

def WKresidual(C,R,p,q,t):
    """ Calculates error between Windkessel model and Data
        Args:
            C - compliance
            R - Resistance
            p - pressure measurements in time domain
            Qfft - Discrete fourier transform of flows corresponding to p
            w - frequencies of DFT
            t - times of p
        Returns:
            r : the sum of squared errors between pwk and p
    """
    N=len(q)
    dt=t[1]-t[0]
    T = N*dt;
    w = np.arange(N//2+1)*2*np.pi/T
    Qfft = np.fft.rfft(q)
    Pfft = np.fft.rfft(p)
    Zdata=Pfft/Qfft
    Z=ZWK2(R,C,w)
    Pwk_fft=Z*Qfft
    r=np.sum(np.abs(Pfft-Pwk_fft)**2)/N
    return r


# In[7]:

#===============================================================================
# Guess a compliance value
#===============================================================================
R=np.mean(p)/np.mean(q) # Total arterial reistance is given by the ratio of 
                        # average pressure to average flow. This is not always 
                        # the same as the peripheral resistance
C_guess=0.9


#---Compute the pressure pwk0 from impedance based on the guessed compliance value
Qfft=np.fft.rfft(q)
Z = ZWK2(R,C_guess,w)
pwk0 = np.fft.irfft(Qfft*Z,len(q))  #the len(q) has has to be specified when len(q) is odd

def residual(opt_args):
    """Transform WKresidual to a residual function of C, 
    which is the only scalar argument"""
    return WKresidual(opt_args[0], R, p, q,t)

#===============================================================================
#  Estimate the compliance 
#===============================================================================
import scipy.optimize

res = scipy.optimize.minimize(residual, [C_guess,])
#res = scipy.optimize.minimize_scalar(lambda C: WKresidual(C, R, p, q,t))
C=res.x[0]
print ('C=%f'%C)


Z = ZWK2(R,C,w)
pwk = np.fft.irfft(Z*Qfft,len(q))


#===============================================================================
# Plot solutions
#===============================================================================
plt.figure()
plt.title('Predicted pressures based on estimated compliances')
_=plt.plot(t,pwk0, label="C guess")
_=plt.plot(t,pwk, label="C estimated fft")
_=plt.plot(t,p,'k', label="data")
_=plt.legend()
plt.xlabel("Time")
plt.ylabel("Pressure")

