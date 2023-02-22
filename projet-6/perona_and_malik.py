import numpy as np
import matplotlib.pyplot as mp
import scipy.signal as signal



##compute a function matrix
##Input : a matrix
##output : a function matrix
def matrix_to_function(matrix):
    def u(x=0, y=0, get_size=False):
        if get_size:
            return np.shape(matrix)
        return matrix[x,y]
    return u


##do the opposite of matrix_to_function
##Input : a function u(x,y)
##Output : a matrix of size nxm as matrix[k,l] = u(k,l)
def function_to_matrix(u):
    matrix = np.zeros((u(get_size=True)))
    s = np.shape(matrix)
    for k in range(s[0]):
        for l in range(s[1]):
            matrix[k][l] = u(k,l)
    return matrix


##Create a function image
##Input : an array of size nxm with values in [0,1]
##Output : a function u(x,y) that returns the value of the image at (x,y)
def create_image(image):
    return matrix_to_function(image)



###############################################################################################################################################################################
###############################################################################################################################################################################



##return the gradient matrix of the image U
##Input : a matrix such as u[x,y] that returns the color pixel of the image at (x,y)
##Output : a matrix of size nxm with the gradient of u at each point
def gradient(u):
    (N,M) = np.shape(u)
    gradu = np.zeros((N,M,2))
    for k in range(N):
        for l in range(M):
            gradu[k][l] = [0,0]
            if k < N-1:
                gradu[k][l][0] = u[k+1,l] - u[k,l]
            else:
                gradu[k][l][0] = 0
            
            if l < M-1:
                gradu[k][l][1] = u[k,l+1] - u[k,l]
            else:
                gradu[k][l][1] = 0
    return gradu



##it is the conjugate operation of the gradient
##Input : a 2 dimensional vector field
##output : the gradient of the vector field
def divergence(p):
    (M, N, d) = np.shape(p)
    divp = np.zeros((M,N))
    for k in range(M):
        for l in range(N):
            if 0 < k and k < N-1:
                s1 = p[k][l][0]
            elif k == 0:
                s1 = p[k][l][0]
            else:
                s1 = -p[k-1][l][0]
            
            if 0 < l and l < M-1:
                s2 = p[k][l][1] - p[k][l-1][1]
            elif l == 0:
                s2 = p[k][l][1]
            else:
                s2 = -p[k][l-1][1]
            
            divp[k][l] = s1 + s2
    
    return divp


##implementation of the heat equation
##input: u0 the image you want to apply the heat equation to, h the time step and N the number of iterations
##output: the image after the heat equation has been applied
def heat_equation(u0, h, N):
    u = u0
    for i in range(N):
        u = u + h * divergence(gradient(u))
    return u



###############################################################################################################################################################################
###############################################################################################################################################################################

##return a gaussian function from an image
## Input : a matrix
##Output : a matrix reprensting the gaussian function of the image function
def gaussian(image):
    (M,N) = np.shape(image)
    G = np.zeros((M,N))
    sigma2 = image[M//2,0] - 2*np.log(0.1) 
    sigma2 *= sigma2
    C = 2*np.pi*sigma2
    for k in range(M):
        for l in range(N):
            # G[k,l] = np.exp(-((k-M/2)**2 + (l-N/2)**2) / (2*(M/2)**2)) ##maybe it has to be modified
            G[k,l] = np.exp(-(k**2 + l**2)/(2*sigma2))/C
    return G


##return the convolution of an image with a gaussian
##Input : a matrix
##Output : a matrix reprenting the convulation function
def convolution_with_gaussian(image):
    G = gaussian(image)
    F = signal.convolve2d(image, G, mode='same', boundary='symm')
    return F


##Return a function f such as it is close to one near the contours of the image
##input : a matrix
##output : a matrix reprenseting the function f
def function_modulator(image, already_made = False):
    if already_made:
        return np.genfromtxt("func.txt")
    (M,N) = np.shape(image)
    f = np.zeros((M,N))
    F = convolution_with_gaussian(image)
    gradF = gradient(F)
    for k in range(M):
        for l in range(N):
            f[k,l] = np.exp(-np.linalg.norm(gradF[k,l])**2)
    return f



def term_product(f, x):
    (M, N) = np.shape(f)
    result = np.zeros((M,N,2))
    for k in range(M):
        for l in range(N):
            result[k,l] = f[k,l] * x[k,l]
    return result


##implementation of the heat equation
##input: u0 the image you want to apply the heat equation to, h the time step and N the number of iterations
##output: the image after the heat equation has been applied
def perona_and_malik(u0, h = 0.2, N = 200):
    u = u0
    f = function_modulator(u, True)
    for i in range(N):
        u = u + h * divergence(term_product(f, gradient(u)))
    return u

###############################################################################################################################################################################
###############################################################################################################################################################################


##plot the image
##Input : a matrix with coefficient between 0 and 1
def plot_image(image):
    mp.imshow(image, cmap='gray')
    mp.show()


##plot an image modified with the heat equation filter
#Input : a matrix with coefficient between 0 and 1
def plot_heat_map_filter(image):
    plot_image(heat_equation(image, 0.1, 10))


##plot an image modified with the heat equation filter and the orignal image
##input : a matrix with coefficient between 0 and 1
def compare_plot_image_normal_and_heat(image):
    fig, axs = mp.subplots(1, 2)
    axs[0].imshow(image, cmap='gray')
    axs[0].set_title("Original image")

    image_heat = heat_equation(image, 0.1, 20)
    axs[1].imshow(image_heat, cmap='gray')
    axs[1].set_title("Original image with the heat equation filter")

    mp.show()

def plot_perona_malik_filter(image):
    plot_image(perona_and_malik(image, N = 20))


##plot the modulator needed for the perona and malik equations
##Inpput : a matrix reprenting the image
def plot_function_modulator(image):
    f = function_modulator(image)
    mp.imshow(f, cmap='gray')
    mp.show()


def compare_plot_image_normal_and_perona_malik(image):
    fig, axs = mp.subplots(1, 2)
    axs[0].imshow(image, cmap='gray')
    axs[0].set_title("Original image")

    image_pema = perona_and_malik(image, 0.2, 20)
    axs[1].imshow(image_pema, cmap='gray')
    axs[1].set_title("Original image with the Perona and Malik equation filter")

    mp.show()


def compare_plot_heat_and_perona_malik(image):
    fig, axs = mp.subplots(1, 2)
    axs[0].imshow(heat_equation(image, 0.2, 20), cmap='gray')
    axs[0].set_title("Original image with the heat equation filter")

    image_pema = perona_and_malik(image, 0.2, 20)
    axs[1].imshow(image_pema, cmap='gray')
    axs[1].set_title("Original image with the Perona and Malik equation filter")

    mp.show()


def compare_plot_normal_heat_and_perona_malik(image):
    fig, axs = mp.subplots(1, 3)
    axs[0].imshow(image, cmap='gray')
    axs[0].set_title("Original image")

    image_heat = heat_equation(image, 0.2,30)
    axs[1].imshow(image_heat, cmap='gray')
    axs[1].set_title("Original image with the heat equation filter")

    image_pema = perona_and_malik(image, 0.2, 30)
    axs[2].imshow(image_pema, cmap='gray')
    axs[2].set_title("Original image with the Perona and Malik equation filter")

    mp.show()

if __name__ == "__main__":
    image = mp.imread("./p6_peronamalik_original.png")
    # plot_image(image)
    # plot_heat_map_filter(image)
    # compare_plot_image_normal_and_heat(image)
    # plot_function_modulator(image)
    # plot_perona_malik_filter(image)
    # compare_plot_heat_and_perona_malik(image)
    compare_plot_normal_heat_and_perona_malik(image)