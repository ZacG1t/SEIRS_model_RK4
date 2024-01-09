import numpy as np 
import pandas as pd 
import matplotlib as mpl 
import matplotlib.pyplot as plt 
import seaborn as sns

# Define arrays 
t = []
S = []
E = []
I = [] 
R = []
dS = [] 
dE = [] 
dI = [] 
dR = []

# Define initial variables 
S.append(8771761)    # No of population susceptible to Rubella
E.append(142)        # No of population with symptoms of Rubella 
I.append(37)        # No of population infected with Rubella 
R.append(30)         # No of population recovered from Rubella 

# Parameter values 
h = 1
N = 1
sigma = 0.309 
beta = 0.00000002 
alpha = 0.4 
gamma = 0.1667 
delta = 0.14 
mu = 0.0012 
theta = 0.0187

# Define system of differential equations 
# dS/dt 
def dS(R, E, S):
    return (1 - sigma)*theta*N + delta*R - beta*S*E - mu*S 
# dE/dt
def dE(S, E):
    return beta*S*E - alpha*E - mu*E
# dI/dt 
def dI(E, I):
    return alpha*E - gamma*I - mu*I 
# dR/dt 
def dR(I, R):
    return gamma*I + sigma*theta*N - delta*R - mu*R


for i in range(0,100):
    t.append(i)
    # Runge-Kutta Method Order 4 for dS/dt
    k1 = h * dS(R[i], E[i], S[i])
    k2 = h * dS(R[i]+(h/2), E[i]+(h/2), S[i]+(k1/2))
    k3 = h * dS(R[i]+(h/2), E[i]+(h/2), S[i]+(k2/2))
    k4 = h * dS(R[i]+h, E[i]+h, S[i]+k3)
    S.append(S[i] + (1/6)*(k1 + 2*k2 + 2*k3 + k4))

    # Runge-Kutta Method order 4 for dE/dt 
    l1 = h * dE(S[i], E[i]) 
    l2 = h * dE(S[i]+(h/2), E[i]+(l1/2)) 
    l3 = h * dE(S[i]+(h/2), E[i]+(l2/2)) 
    l4 = h * dE(S[i]+h, E[i]+l3)
    E.append(E[i] + (1/6)*(l1 + 2*l2 + 2*l3 + l4))

    # Runge-Kutta Method order 4 for dI/dt 
    m1 = h * dI(E[i], I[i]) 
    m2 = h * dI(E[i]+(h/2), I[i]+(m1/2))
    m3 = h * dI(E[i]+(h/2), I[i]+(m2/2))
    m4 = h * dI(E[i]+h, I[i]+m3)
    I.append(I[i] + (1/6)*(m1 + 2*m2 + 2*m3 + m4))

    # Runge-Kutta Method order 4 for dR/dt 
    n1 = h * dR(I[i], R[i]) 
    n2 = h * dR(I[i]+(h/2), R[i]+(n1/2)) 
    n3 = h * dR(I[i]+(h/2), R[i]+(n2/2)) 
    n4 = h * dR(I[i]+h, R[i]+n3)
    R.append(R[i] + (1/6)*(n1 + 2*n2 + 2*n3 + n4))

t = list(range(0,101))
d = {'t':t, 'S(t)':S, 'E(t)':E, 'I(t)':I, 'R(t)':R}
df = pd.DataFrame(data=d)
df

# Plot 
figure, axis = plt.subplots(2, 2, figsize=(12, 7))
axis[0, 0].scatter(t, S, s=10)
axis[0, 0].set_title("No of population susceptible to Rubella over time")
axis[0, 1].scatter(t, E, s=10)
axis[0, 1].set_title("No of population with symptoms of Rubella over time")
axis[1, 0].scatter(t, I, s=10)
axis[1, 0].set_title("No of population infected with Rubella over time")
axis[1, 1].scatter(t, R, s=10)
axis[1, 1].set_title("No of population recovered from Rubella over time")
plt.show()

# Export data to csv file
df.to_csv('data.csv')