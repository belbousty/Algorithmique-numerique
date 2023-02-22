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

# comparaison avec np.log
plt.ylabel("y")
plt.xlabel("x")

y_log =np.empty((1,100))
x_log = np.linspace(1/1000,99/100,100)
for i in range (0,len(x_log)):
    y_log[0][i] = calcule_log(x_log[i])
plt.plot(x_log, y_log[0], label="calcul_ln(x)")
plt.plot(x_log, np.log(x_log), label="ln(x)")  # Les courbes sont équivalentes !! meme pour [1, 100]
plt.legend()
plt.show()
##################################################
erreur_direct_log = np.empty((1,100))
for i in range (0,len(x_log)):
    erreur_direct_log[0][i] = abs((calcule_log(x_log[i])-np.log(x_log[i])))/abs(np.log(x_log[i]))
#plt.plot(x_log, erreur_direct_log[0], label="Courbe de l'erreur relative de la fonction ln")

erreur_inverse = []
for i in range(len(x_log)):
    erreur_inverse.append(abs(np.exp(y_log[0][i]) - x_log[i]) / x_log[i])
    #erreur_inverse.append(np.exp(y_log[0][i]))
plt.plot(x_log, erreur_inverse, label = "Courbe de l'erreur inverse de la fonc ln")
plt.xscale('log')
plt.legend()
plt.show()

#print("x = ", x_log[3], " np.log = ", np.log(x_log[3]), " faux ln = ", y_log[0][3], " x de l'erreur ", np.exp(y_log[0][3]) )


