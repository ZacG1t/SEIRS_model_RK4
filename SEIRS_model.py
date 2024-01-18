import numpy as np 
import pandas as pd 
import matplotlib as mpl 
import matplotlib.pyplot as plt 

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
S0 = 8771761
E0 = 142
I0 = 37
R0 = 30

# Add initial variables to list
S.append(S0)
E.append(E0) 
I.append(I0) 
R.append(R0) 

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
def dS(T, S, E, I, R):
    return (1 - sigma)*theta*N + delta*R - beta*S*E - mu*S 
# dE/dt
def dE(T, S, E, I, R):
    return beta*S*E - alpha*E - mu*E
# dI/dt 
def dI(T, S, E, I, R):
    return alpha*E - gamma*I - mu*I 
# dR/dt 
def dR(T, S, E, I, R):
    return gamma*I + sigma*theta*N - delta*R - mu*R

# Runge-Kutta method calculation
# k for dS/dt
# l for dE/dt 
# m for dI/dt
# n for dR/dt 
for i in range(0,100): 
    
    k1 = h*dS(i, S[i], E[i], I[i], R[i])
    l1 = h*dE(i, S[i], E[i], I[i], R[i])
    m1 = h*dI(i, S[i], E[i], I[i], R[i])
    n1 = h*dR(i, S[i], E[i], I[i], R[i])

    k2 = h*dS(i+0.5*h, S[i]+0.5*(h*k1), E[i]+0.5*(h*l1), I[i]+0.5*(h*m1), R[i]+0.5*(h*n1))
    l2 = h*dE(i+0.5*h, S[i]+0.5*(h*k1), E[i]+0.5*(h*l1), I[i]+0.5*(h*m1), R[i]+0.5*(h*n1))
    m2 = h*dI(i+0.5*h, S[i]+0.5*(h*k1), E[i]+0.5*(h*l1), I[i]+0.5*(h*m1), R[i]+0.5*(h*n1))
    n2 = h*dR(i+0.5*h, S[i]+0.5*(h*k1), E[i]+0.5*(h*l1), I[i]+0.5*(h*m1), R[i]+0.5*(h*n1))

    k3 = h*dS(i+0.5*h, S[i]+0.5*(h*k2), E[i]+0.5*(h*l2), I[i]+0.5*(h*m2), R[i]+0.5*(h*n2))
    l3 = h*dE(i+0.5*h, S[i]+0.5*(h*k2), E[i]+0.5*(h*l2), I[i]+0.5*(h*m2), R[i]+0.5*(h*n2))
    m3 = h*dI(i+0.5*h, S[i]+0.5*(h*k2), E[i]+0.5*(h*l2), I[i]+0.5*(h*m2), R[i]+0.5*(h*n2))
    n3 = h*dR(i+0.5*h, S[i]+0.5*(h*k2), E[i]+0.5*(h*l2), I[i]+0.5*(h*m2), R[i]+0.5*(h*n2))

    k4 = h*dS(i+h, S[i]+h*k3, E[i]+h*l3, I[i]+h*m3, R[i]+h*n3)
    l4 = h*dE(i+h, S[i]+h*k3, E[i]+h*l3, I[i]+h*m3, R[i]+h*n3)
    m4 = h*dI(i+h, S[i]+h*k3, E[i]+h*l3, I[i]+h*m3, R[i]+h*n3)
    n4 = h*dR(i+h, S[i]+h*k3, E[i]+h*l3, I[i]+h*m3, R[i]+h*n3)

    S.append(S[i] + (1/6)*h*(k1 + 2*k2 + 2*k3 + k4))
    E.append(E[i] + (1/6)*h*(l1 + 2*l2 + 2*l3 + l4))
    I.append(I[i] + (1/6)*h*(m1 + 2*m2 + 2*m3 + m4))
    R.append(R[i] + (1/6)*h*(n1 + 2*n2 + 2*n3 + n4))

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
df.to_csv('result.csv')