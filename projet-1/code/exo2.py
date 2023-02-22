from cProfile import label
import numpy as np
import matplotlib.pyplot as plt 

L= [np.log(1+10**(-i)) for i in range(0, 7)]
A= [np.arctan(10**(-i)) for i in range (0,5)]

def calcule_exp(x):
    k,y=0,1
    while k <= 6:
        while x >= L[k]:
            x -= L[k]
            y += y*10**(-k)
        k+=1
    return y+y*x

def calcule_log(x):
    k,y,p=0,0,1
    while k <= 6:
        while x >= p+p*10**(-k):
            y += L[k]
            p += p*10**(-k)
        k+=1
    return y+(x/p-1)

def calcule_arctan(x):
    if x < 0:
        return -calcule_arctan(-x)
    k,y,r = 0, 1, 0
    while (k < 4):  
        while (x < y*10**(-k)): 
            k += 1
            if (k >= 4): # pour s'assurer que k ne dépasse pas la taille de A
                break
        xp = x - y*10**(-k)
        y += x*10**(-k)
        x = xp 
        r += A[k]
    return r+(x/y)

def calcule_tan(x:int):
    k,n,d = 0, 0, 1
    while (k <= 4):
        while (x >= A[k]):
            x-= A[k]
            np = n + d*10**(-k)
            d -= n*10**(-k)
            n = np
        k += 1
    return (n+x*d)/(d-n*x)

# Question 5

def Tests():
    print("--Tests begin--")
    
    assert calcule_exp(0) == 1
    assert calcule_log(1) == 0
    assert calcule_tan(np.pi/4) == 1

    print("--Tests ended successfullly--")

Tests()

fig = plt.figure()
fig,ax = plt.subplots(2,1)

y_exp = np.empty((1,100))
x_exp = np.linspace(0,20,100)
for i in range (0,len(x_exp)):
    y_exp[0][i] = calcule_exp(x_exp[i])

ax[0].plot(x_exp, y_exp[0], label="y = exp_cordic(x)")
ax[0].plot(x_exp, np.exp(x_exp), label="y = exp(x)")

ax[0].legend()
#plt.show()
################################################
erreur_exp = np.empty((1,100))
for i in range (0,len(x_exp)):
    erreur_exp[0][i] = abs((calcule_exp(x_exp[i])-np.exp(x_exp[i])))/abs(np.exp(x_exp[i]))
ax[1].plot(x_exp, erreur_exp[0], label="Erreur relative de la fonction exp avec Cordic")
ax[1].legend()
plt.show()

#################################################
# comparaison avec np.log

fig = plt.figure()
fig,ax = plt.subplots(2,1)

y_log =np.empty((1,1000))
x_log = np.linspace(1/1000,2,1000)
for i in range (0,len(x_log)):
    y_log[0][i] = calcule_log(x_log[i])
ax[0].plot(x_log, y_log[0], label="y = log(x) avec Cordic")
ax[0].plot(x_log, np.log(x_log), label="y = ln(x)")  # Les courbes sont équivalentes !! meme pour [1, 100]
ax[0].legend()

##################################################
erreur_log = np.empty((1,1000))
for i in range (0,len(x_log)-1):
    erreur_log[0][i] = abs((calcule_log(x_log[i])-np.log(x_log[i])))/abs(np.log(x_log[i]))
ax[1].plot(x_log, erreur_log[0], label="Erreur relative de la fonction log avec Cordic")
ax[1].legend()
plt.show()

##################################################
# comparaison avec np.tan
y_tan =np.empty((1,100))
x_tan = np.linspace(-np.pi/2+0.1,np.pi/2-0.1,100)
for i in range (0,len(x_tan)):
    y_tan[0][i] = calcule_tan(x_tan[i])
plt.plot(x_tan, y_tan[0], label="calcul_tan(x)")
plt.plot(x_tan, np.tan(x_tan), label="tan(x)") # Les courbes sont équivalentes !! 
plt.legend()
plt.show()
##################################################
erreur_tan = np.empty((1,100))
for i in range (0,len(x_tan)):
    erreur_tan[0][i] = abs((calcule_tan(x_tan[i])-np.tan(x_tan[i])))/abs(np.tan(x_tan[i]))
plt.plot(x_tan, erreur_tan[0], label="courbe d'erreur de la fonction tan")
plt.legend()
plt.show()

##################################################

fig = plt.figure()
fig,ax = plt.subplots(2,1)

y_arctan = np.empty((1,100))
x_arctan = np.linspace(-10,30,100)
for i in range (0,len(x_arctan)):
    y_arctan[0][i] = calcule_arctan(x_arctan[i])
ax[0].plot(x_arctan, y_arctan[0], label="y = arctan(x) avec Cordic")
ax[0].plot(x_arctan, np.arctan(x_arctan), label="y = arctan(x)") # Les courbes sont équivalentes !! 
ax[0].legend()

##################################################
erreur_arctan = np.empty((1,100))
for i in range (0,len(x_arctan)):
    erreur_arctan[0][i] = abs((calcule_arctan(x_arctan[i])-np.arctan(x_arctan[i])))/abs(np.arctan(x_arctan[i]))
ax[1].plot(x_arctan, erreur_arctan[0], label="Erreur relative de la fonction arctan avec Cordic")
##ax[1].yscale('log')
ax[1].legend()
plt.show()
