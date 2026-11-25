import numpy as np
from matplotlib import pyplot as plt
import h5py



# Output = h5py.File('./HDF_Files/ConvergenceTest_N3200_K1p0.hdf5', 'r')
Output = h5py.File('./HDF_Files/test_atmosphere.hdf5', 'r')
t = Output['Time'][:]
boundaries = Output['Fluid_boundaries'][:,:]
soln = Output['Solution/Solution_Array'][:,:,:]
coordinates = Output['x_coordinates'][:]
dx = np.abs(coordinates[1]-coordinates[0])
N = np.shape(soln)[2]
K = Output['Aux'][0]


def compute_integral_abs(soln,dx):

    norm = np.abs(dx*np.sum(soln[:]))
        
    return norm


def compute_l2norm(soln,dx):

    norm = np.sqrt(dx*np.sum(soln[:]**2))
        
    return norm

def DxF(f,dx):
    """ This function calculates the first derivative
        using a 2nd order forward finite difference stencil and backward stencil at right boundary"""
    h = dx
    n = np.shape(f)[0]
    df = np.zeros(np.shape(f))

    # #2nd Order
    # df[0:n-2] = (-3/2*f[0:n-2]+2*f[1:n-1]-0.5*f[2:n])/h
    # df[0] = (-3/2*f[0]+2*f[1]-0.5*f[2])/h
    # df[n-2] = (3/2*f[n-2]-2*f[n-3]+0.5*f[n-4])/h
    # df[n-1] = (3/2*f[n-1]-2*f[n-2]+0.5*f[n-3])/h
    

    # #4th Order 
    df[0] = (-25/12*f[0]+4*f[1]-3*f[2]+4/3*f[3]-1/4*f[4])/h
    df[1] = (-25/12*f[1]+4*f[2]-3*f[3]+4/3*f[4]-1/4*f[5])/h
    df[2:n-2] = (1/12*f[0:n-4]-2/3*f[1:n-3]+2/3*f[3:n-1]-1/12*f[4:n])/h
    df[n-2] = (25/12*f[n-2]-4*f[n-3]+3*f[n-4]-4/3*f[n-5]+1/4*f[n-6])/h
    df[n-1] = (25/12*f[n-1]-4*f[n-2]+3*f[n-3]-4/3*f[n-4]+1/4*f[n-5])/h

    return df

###################################
# Recover Primitive Variables
###################################
T00 = soln[:,0,:]
T01 = soln[:,1,:]
a = soln[:,2,:]
v = np.zeros_like(T00)
mu = np.zeros_like(T00)
Gamma2 = np.zeros_like(T00)



#-------- Explicit Analytic Formula
for i in range(np.shape(t)[0]):

    mask1 = T00[i,:] >= 0.5*(1.0/K) + 0.5*np.sqrt((1.0 + 4.0*K**2*T01[i,:]**2)/(K**2))

    mu[i,mask1] = (1.0 + (K - 1.0)*T00[i,mask1] + np.sqrt((-1.0 + T00[i,mask1] + K*T00[i,mask1])**2 - 4.0*K*T01[i,mask1]**2))/(2.0*K)

    v[i,mask1] = ((2.0*T01[i,mask1]))/(-1.0 + (K + 1.0)*T00[i,mask1] + np.sqrt((-1.0 + (K+1.0)*T00[i,mask1])**2 - 4.0*K*T01[i,mask1]**2))

    mask2 = T00[i,:] < 0.5*(1.0/K) + 0.5*np.sqrt((1.0 + 4.0*K**2*T01[i,:]**2)/(K**2))
    v[i,mask2] = T01[i,mask2]/T00[i,mask2]
    Gamma2[i,mask2] = 1.0/(1.0 - v[i,mask2]**2)
    mu[i,mask2] = T00[i,mask2]/Gamma2[i,mask2]

        


# plt.rcParams.update({"text.usetex": True,
#     "font.family": "serif",
#     "font.serif": "Computer Modern",
#     "savefig.bbox": "tight",
#     "savefig.format": "pdf"})
# plt.rc('font', size=16)


### Plot conservation
# norm_list = []
# for i in range(0,np.shape(t)[0]):
#     norm = compute_integral_abs(soln[i,0,:]-soln[0,0,:],dx)
#     norm_list.append(norm)


# plt.plot(t,norm_list)
# plt.yscale('log')
# plt.show()

# plt.plot(coordinates,v[0,:])
# plt.plot(t,boundaries[:,1])
# plt.ylim([1.099,1.1001])
# plt.show()

# for i in range(0,np.shape(t)[0],1):

#     mask1 = coordinates>boundaries[i,0]
#     mask2 = coordinates<boundaries[i,1]
#     mask = mask1==mask2

#     plt.plot(t[i],mu[i,mask][-1],'x',color='tab:blue')

# plt.show()    


### Animation of solution
for i in range(0,np.shape(t)[0]):

   

    # plt.plot(coordinates[:],mu[i,:][:],'--')
    mask = coordinates<2
    plt.plot(coordinates[mask],v[i,mask][:],'--')
    # plt.plot(coordinates[:],a[i,:],'--')
    # plt.plot(coordinates[mask][mask_pressure],mu[i,mask][mask_pressure],'x',color='tab:red')


    # plt.xlim([0,1.2])
    plt.legend()

    # if np.any(K*mu[i,mask]-1<0):
        # print(np.where(K*mu[i,mask]-1))
        # print(t[i])
        # breaks



    # if t[i]>t_crit:
    #     break

    # print(soln[i,0,mask] <= 0.5*K + 0.5*np.sqrt((1.0 + 4.0*K**2*soln[i,1,mask]**2)/(K**2)))
    plt.title(t[i])
    plt.draw()
    plt.pause(0.01)
    plt.cla()








