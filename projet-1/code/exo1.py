import numpy as np
import matplotlib.pyplot as plt

#Premier essaie de la fonction rp

def rp(x,p):
    signe = 1

    if abs(x) < 10**(-10):
        return 0
    if x < 0 :
        x = -x
        signe = -1

    k = 0
    while x >= 1:
        k = k + 1
        x = x / 10.0
   

    while x < 0.1:
        x = x * 10
        k = k - 1
     

    x = x * 10**p
    reste = x%1
    if reste >= 0.5:
        x = x + 1
    x = int(x)

    x = x * 10**(k-p)
        
    return x

#Second essaie de la fonction rp 

def rp1(x,p):
    x = str(x)
    x = list(x)
    taille = len(x)
    nb_carac = 0
    cmp = 0

    while ((cmp < taille) and ((x[cmp] == '0') or (x[cmp] == '.'))):
        cmp+=1
    
    while (cmp < taille):
        if (nb_carac >= p) and (x[cmp] != '.'):
            x[cmp] = '0'

        if x[cmp] != '.':
            nb_carac+=1 

        cmp+=1
   

    
    y = ""
    for k in range(cmp):
        y+=str(x[k])

    y = float(y)


    return y

#Graphes pour regarder le fonctionnement de rp1

def graph_rp():

    L0 = []
    L1 = []
    L2 = [] 
    L3 = []
    L4 = []
    L5 = []
    L6 = [] 
    X = np.linspace(0, 1, 1000)
    for i in X:
            L1.append(rp(i, 1))
            L2.append(rp(i, 2))
            L3.append(rp(i, 3))

    Y1 = np.array(L1)
    Y2 = np.array(L2)
    Y3 = np.array(L3)

    plt.plot(X, Y1, label="Avec p = 1")
    plt.plot(X, Y2, label="Avec p = 2")
    plt.plot(X, Y3, label="Avec p = 3")
    plt.legend(prop = {'size' : 12})
    plt.grid()
    plt.show()
        
#Erreur relative sur la somme
 
def err_rel_sum(x,y,p):
    rpxy = round(x+y,p)

    return (abs(x+y - rpxy)/abs(x+y))

#Erreur relative sur le produit

def err_rel_prod(x,y,p):
    rpxy = round(x*y,p)
    if x == 0:
        return 0
    else:
        return (abs(x*y - rpxy)/abs(x*y))

#Affichage des gaphes pour l'erreur relative sur 
#la somme et le produit autour de x

def graph_erreur():
    x = 0.193267

    y1 = np.linspace(-0.193267-0.193267,-0.193267+0.193267,200)

    s1 = []
    s2 = []
    s3 = []
    s4 = []

    p1 = []
    p2 = []
    p3 = []
    p4 = []

    for i in y1:
        s1.append(err_rel_sum(x,i,1))
        s2.append(err_rel_sum(x,i,2))
        s3.append(err_rel_sum(x,i,3))
        s4.append(err_rel_sum(x,i,4))

        p1.append(err_rel_prod(x,i,1))
        p2.append(err_rel_prod(x,i,2))
        p3.append(err_rel_prod(x,i,3))
        p4.append(err_rel_prod(x,i,4))

    s1 = np.array(s1)
    s2 = np.array(s2)
    s3 = np.array(s3)
    s4 = np.array(s4)

    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = np.array(p3)
    p4 = np.array(p4)


    fig = plt.figure()
    fig,ax = plt.subplots(2,1)
    ax[0].set_yscale('log')
    ax[1].set_yscale('log')
    #ax.plot(y,p)

    ax[0].plot(y1, s1, label='Avec p = 1')
    ax[0].plot(y1, s2, label='Avec p = 2')
    ax[0].plot(y1, s3, label='Avec p = 3')
    ax[0].plot(y1, s4, label='Avec p = 4')
    ax[0].set_title('Erreur relative sur la somme')

    ax[1].plot(y1, p1, label='Avec p = 1')
    ax[1].plot(y1, p2, label='Avec p = 2')
    ax[1].plot(y1, p3, label='Avec p = 3')
    ax[1].plot(y1, p4, label='Avec p = 4')
    ax[1].set_title('Erreur relative sur le produit')

    ax[0].legend(prop = {'size' : 12})
    ax[0].grid()

    ax[1].legend(prop = {'size' : 12})
    ax[1].grid()

    plt.show()

graph_rp()
graph_erreur()