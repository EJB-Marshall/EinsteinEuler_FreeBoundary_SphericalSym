import numpy as np
from matplotlib import pyplot as plt
import h5py

import numpy as np
import h5py
from scipy import integrate
import matplotlib.pyplot as plt
from scipy.optimize import newton
from scipy.interpolate import CubicSpline

def Dx(f,dx):
    """ This function calculates the first derivative
        using a 2nd or 4th order central finite difference stencil"""
    h = dx
    n = np.shape(f)[0]
    df = np.zeros(np.shape(f))

    # # #2nd Order Central Finite Difference
    df[1:n-1] = (-f[0:n-2]+f[2:n])/(2*h)

    #Periodic Boundary Condition
    df[0] = (-f[n-1]+f[1])/(2*h)
    df[n-1] = (-f[n-2]+f[0])/(2*h)

    # #4th Order Central Finite Difference
    # df[2:n-2] = (f[0:n-4]-8*f[1:n-3]+8*f[3:n-1]-f[4:n])/(12*h)

    # #Periodic Boundary Condition
    # df[1] = (f[n-1]-8*f[0]+8*f[2]-f[3])/(12*h)
    # df[0] = (f[n-2]-8*f[n-1]+8*f[1]-f[2])/(12*h)
    # df[n-2] = (f[n-4]-8*f[n-3]+8*f[n-1]-f[0])/(12*h)
    # df[n-1] = (f[n-3]-8*f[n-2]+8*f[0]-f[1])/(12*h)

    return df



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

def Dx3(f,dx):
    """ This function calculates the third derivative
        using a 4th order central finite difference stencil"""
    h = dx
    n = np.shape(f)[0]
    df = np.zeros(np.shape(f))

    
    #4th Order Central Finite Difference
    df[3:n-3] = (1/8*f[0:n-6]-1*f[1:n-5]+13/8*f[2:n-4]\
                 -13/8*f[4:n-2]+1*f[5:n-1]-1/8*f[6:n])/(h**3)

    #Periodic Boundary Condition
    df[2] = (1/8*f[-1]-1*f[0]+13/8*f[1]\
                 -13/8*f[3]+1*f[4]-1/8*f[5])/(h**3)
    df[1] = (1/8*f[-2]-1*f[-1]+13/8*f[0]\
                 -13/8*f[2]+1*f[3]-1/8*f[4])/(h**3)
    df[0] = (1/8*f[-3]-1*f[-2]+13/8*f[-1]\
                 -13/8*f[1]+1*f[2]-1/8*f[3])/(h**3)
    
    df[n-3] = (1/8*f[n-6]-1*f[n-5]+13/8*f[n-4]\
                 -13/8*f[n-2]+1*f[n-1]-1/8*f[0])/(h**3)
    df[n-2] = (1/8*f[n-5]-1*f[n-4]+13/8*f[n-3]\
                 -13/8*f[n-1]+1*f[0]-1/8*f[1])/(h**3)
    df[n-1] = (1/8*f[n-4]-1*f[n-3]+13/8*f[n-2]\
                 -13/8*f[0]+1*f[1]-1/8*f[2])/(h**3)
    
    return df

def Dx4(f,dx):
    """ Computes periodic 4th derivative for artificial viscosity terms
      Inputs: f - function to be differentiated
                r - Grid instance"""
    
    n = np.shape(f)[0]
    df = np.zeros(np.shape(f))

    h = dx
    df[2:n-2] = (1*f[0:n-4]-4*f[1:n-3]+6*f[2:n-2]-4*f[3:n-1]+1*f[4:n])/h**4

    df[1] = (1*f[n-1]-4*f[0]+6*f[1]-4*f[2]+1*f[3])/h**4
    df[0] = (1*f[n-2]-4*f[n-1]+6*f[0]-4*f[1]+1*f[2])/h**4
    df[n-2] = (1*f[n-4]-4*f[n-3]+6*f[n-2]-4*f[n-1]+1*f[0])/h**4
    df[n-1] = (1*f[n-3]-4*f[n-2]+6*f[n-1]-4*f[0]+1*f[1])/h**4

    return df

def DxF2(f,dx):
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
    
    # # #4th Order 
    df[0:n-4] = (-25/12*f[0:n-4]+4*f[1:n-3]-3*f[2:n-2]+4/3*f[3:n-1]-1/4*f[4:n])/h
    df[n-4] = (25/12*f[n-4]-4*f[n-5]+3*f[n-6]-4/3*f[n-7]+1/4*f[n-8])/h
    df[n-3] = (25/12*f[n-3]-4*f[n-4]+3*f[n-5]-4/3*f[n-6]+1/4*f[n-7])/h
    df[n-2] = (25/12*f[n-2]-4*f[n-3]+3*f[n-4]-4/3*f[n-5]+1/4*f[n-6])/h
    df[n-1] = (25/12*f[n-1]-4*f[n-2]+3*f[n-3]-4/3*f[n-4]+1/4*f[n-5])/h

    return df


############################################################################
#Evaluate constraint quantity and calculate L2 norm functions
############################################################################

def computel2norm(f,dx):
    """This function outputs a vector whose 
    ith entry is the L2 norm at timestep i"""
    l2norm = np.zeros_like(f[:,1])
    for i in range(0,np.shape(f)[0]):
        l2norm[i] = np.sqrt(dx*np.sum(f[i,:]**2))
    return l2norm

def computel1norm(f,dx):
    """This function outputs a vector whose 
    ith entry is the L2 norm at timestep i"""
    l1norm = np.zeros_like(f[:,1])
    for i in range(0,np.shape(f)[0]):
        l1norm[i] = dx*np.sum(np.abs(f[i,:]))
    return l1norm

def compute_integral_abs(f,dx):
    """This function outputs a vector whose 
    ith entry is the absolute value of spatial integral at timestep i"""
    int_eval = np.zeros_like(f[:,1])
    for i in range(0,np.shape(f)[0]):
        int_eval[i] = np.abs(dx*np.sum(f[i,:]))
    return int_eval

def computeH3norm(f,dx):
    """This function outputs a vector whose 
    ith entry is the H3 norm at timestep i"""
    l2norm = np.zeros_like(f[:,1])
    for i in range(0,np.shape(f)[0]):
        l2norm[i] = np.sqrt(dx*np.sum(Dx3(f[i,:],dx)**2))
    return l2norm


def computel2norm_single(f,dx):    
    return np.sqrt(dx*np.sum(f[:]**2))

def computel1norm_single(f,dx):    
    return dx*np.sum(np.abs(f[:]))

def computelinfnorm_single(f):    
    return np.max(np.abs(f[:]))

def compute_integral_abs_single(f,dx):
    int_eval = np.abs(dx*np.sum(f[:]))
    return int_eval


def restrict_grid(soln,num_grids):
    """ In a cell-centred scheme, the cell centres will not be aligned
        for different resolutions. Thus, we must average the cell-centred values from
         high resolution runs to align with the low resolution cells. """
    
    n = soln.shape[0] # Number of grid points for high res run
    averaged_runs = []

    averaged_runs.append(soln)

    for k in range(1,num_grids+1):

        averaged_runs.append(np.zeros((int(n/(2**(k))))))

        for parent_index in range(int(n/(2**(k)))):
            child_index = 2 * parent_index

            ### Conservative update
            averaged_runs[k][parent_index] = 0.5*(averaged_runs[k-1][child_index] + averaged_runs[k-1][child_index+1])

    return averaged_runs



def conserved_to_primitive(T00,T01,grid):

    mu = np.zeros_like(T00)
    v = np.zeros_like(T00)

    for i in range(np.shape(t)[0]):
        mask1 = grid>boundaries[i,0]
        mask2 = grid<boundaries[i,1]
        mask = mask1==mask2

        mu[i,mask] = (1.0 + (K - 1.0)*T00[i,mask] + np.sqrt((-1.0 + T00[i,mask] + K*T00[i,mask])**2 - 4.0*K*T01[i,mask]**2))/(2.0*K)

        v[i,mask] = ((2.0*T01[i,mask]))/(-1.0 + (K + 1.0)*T00[i,mask] + np.sqrt((-1.0 + (K+1.0)*T00[i,mask])**2 - 4.0*K*T01[i,mask]**2))

    return [mu,v]


def constraint_eval(soln,boundaries,dx,grid,K):

    tau = soln[:,0,:]
    S = soln[:,1,:]
    a = soln[:,2,:]
    alpha = soln[:,3,:]

    mu,v = conserved_to_primitive(tau,S,grid)

    Gamma2 = 1.0/(1.0 - v**2)

    C_H = np.zeros_like(a)

    for i in range(np.shape(t)[0]):

        star_mask = (grid < boundaries[i,1])
        vacuum_mask = (boundaries[i,1] < grid)

        C_H[i,star_mask] = a[i,star_mask]*(1-a[i,star_mask]**2)/(2*grid[star_mask]) + (grid[star_mask]/2)*a[i,star_mask]**3*(1 - K*mu[i,star_mask] + Gamma2[i,star_mask]*(-1 + (1+K)*mu[i,star_mask])) - DxF(a[i,star_mask],dx)

        C_H[i,vacuum_mask] = a[i,vacuum_mask]*(1-a[i,vacuum_mask]**2)/(2*grid[vacuum_mask]) - DxF(a[i,vacuum_mask],dx)

    return C_H



############################################################################
# Recover Primitive Variables
############################################################################
grid_list = []
convergence_list = []
ID_list = []
t_list = []
dx_list = []
norm_list = []
soln_list = []
N_list = []
boundaries_list = []


# file_names = ['convergence_static_N800_K1p0.hdf5','convergence_static_N1600_K1p0.hdf5','convergence_static_N3200_K1p0.hdf5','convergence_static_N6400_K1p0.hdf5']

# file_names = ['LongEvolution_t5e3_N1600_K1p0_CD1p1_eps1e-3.hdf5','LongEvolution_t5e3_N3200_K1p0_CD1p1_eps1e-3.hdf5','LongEvolution_t5e3_N6400_K1p0_CD1p1_eps1e-3.hdf5']
# file_names = ['test_N800_K1p0_CD1p1_eps1e-3.hdf5','test_N1600_K1p0_CD1p1_eps1e-3.hdf5','test_N3200_K1p0_CD1p1_eps1e-3.hdf5','test_N6400_K1p0_CD1p1_eps1e-3.hdf5']
file_names = ['test_N800_K1p0_CD1p1_Static.hdf5','test_N1600_K1p0_CD1p1_Static.hdf5','test_N3200_K1p0_CD1p1_Static.hdf5','test_N6400_K1p0_CD1p1_Static.hdf5']


for k in range(len(file_names)):
    Output = h5py.File('./HDF_Files/'+file_names[k], 'r')

    time = -1
    Var = 0
    t = Output['Time'][:]
    boundaries = Output['Fluid_boundaries'][:,:]
    soln = Output['Solution/Solution_Array'][:,:,:]
    grid = Output['x_coordinates'][:]
    dx = np.abs(grid[1]-grid[0])
    N = np.shape(soln)[2]
    K = Output['Aux'][0]

    C_H = constraint_eval(soln,boundaries,dx,grid,K)
    l2norm = computel2norm(C_H,dx)

    # Norm for exact TOV initial data test
    l2norm = computel2norm(soln[:,Var,:]-soln[0,Var,:],dx)
    # l2norm = np.abs(boundaries[:,1]-boundaries[0,1]) # Test boundary position against exact solution

    t_list.append(t)
    grid_list.append(grid)
    convergence_list.append(soln[time,Var,:])
    ID_list.append(soln[0,Var,:])
    dx_list.append(dx)
    norm_list.append(l2norm)
    soln_list.append(soln[:,Var,:])
    N_list.append(grid.shape[0])
    boundaries_list.append(boundaries)

    print("The bdy position at t=100 is,", boundaries[-1,1])


### Compute restricted high resolution run
avg_high_res = restrict_grid(convergence_list[3],3)


plt.rcParams.update({"text.usetex": True,
    "font.family": "serif",
    "font.serif": "Computer Modern",
    "savefig.bbox": "tight",
    "savefig.format": "pdf"})
plt.rc('font', size=16)





plt.figure()
plt.plot(t_list[0],np.log2(norm_list[0]),label=r'$800$')
plt.plot(t_list[1],np.log2(norm_list[1]),label=r'$1600$')
plt.plot(t_list[2],np.log2(norm_list[2]),label=r'$3200$')
plt.plot(t_list[3],np.log2(norm_list[3]),label=r'$6400$')
# plt.plot(t_list[3],np.max(np.log(np.abs(E11)),axis=1),label='E11')
# plt.plot(t_list[3],np.max(tau_log,axis=1),label='T^00')
# plt.plot(t_list[4],np.log2(norm_list[4]),label=r'$3200$')
plt.xlabel(r'$t$')
# plt.ylabel(r'$\log_{2}\|C\|_{2}$')
# plt.ylabel(r'$\log_{2}\big(|x_{R}(t) - x_{R}(0)|\big)$')
plt.ylabel(r'$\log_{2}\big(\|\tau(t) - \tau(0)\|_{2}\big)$')
plt.legend()
plt.grid(alpha=0.8)
# plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.tight_layout()
# plt.yscale("log")
# plt.ylim([10**(-15),1])
# print(t_original[time])
plt.show()


# plt.figure()
# plt.plot(grid_list[0],np.log2(np.abs(test_list[0][time,:])),label='400')
# plt.plot(grid_list[1],np.log2(np.abs(convergence_list[1][:])),label='800')
# plt.plot(grid_list[2],np.log2(np.abs(convergence_list[2][:])),label='1600')
# # plt.plot(grid_list[3],np.log2(np.abs(convergence_list[3][:])),label='1600')
# # plt.plot(grid_list[4],np.log2(np.abs(convergence_list[4][:])),label='3200')
# # plt.plot(grid_list[5],convergence_list[5][:],label='3200')
# plt.xlabel(r'$x$')
# plt.ylabel(r'$\log_{2}|\Delta|$')
# plt.legend()
# plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
# plt.tight_layout()
# plt.show()



# print(convergence_list[0][:].shape)
# print(avg_high_res[-1][:].shape)
# print(grid_list[0].shape)

plt.figure(1)
plt.plot(grid_list[0],np.log2(np.abs(convergence_list[0][:]-avg_high_res[-1][:])),label='800')
plt.plot(grid_list[1],np.log2(np.abs(convergence_list[1][:]-avg_high_res[-2][:])),label='1600')
plt.plot(grid_list[2],np.log2(np.abs(convergence_list[2][:]-avg_high_res[-3][:])),label='3200')
# plt.plot(grid_list[3],np.log2(np.abs(convergence_list[3][:]-avg_high_res[-4][:])),label='3200')
# plt.plot(grid_list[4],np.log2(np.abs(convergence_list[4][:]-avg_high_res[-5][:])),label='800')
# plt.plot(grid_list[5],np.log2(np.abs(convergence_list[5][:]-avg_high_res[-6][:])),label='1600')
# plt.plot(grid_list[6],np.log2(np.abs(convergence_list[6][:]-avg_high_res[-7][:])),label='3200')
plt.xlabel(r'$r$')
plt.ylabel(r'$\log_{2}|\Delta|$')
plt.legend()
plt.grid(alpha=0.8)
# plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.tight_layout()
# print(t_original[time])
plt.show()


# plt.figure(1)
# plt.plot(grid_list[0],convergence_list[0][:],'.',markersize=0.75,label='400')
# plt.plot(grid_list[1],convergence_list[1][:],'.',markersize=0.75,label='800')
# plt.plot(grid_list[2],convergence_list[2][:],'.',markersize=0.75,label='1600')
# plt.plot(grid_list[3],convergence_list[3][:],'.',markersize=0.75,label='3200')
# plt.plot(grid_list[4],convergence_list[4][:],'.',markersize=0.75,label='6400')
# plt.xlabel(r'$x$')
# plt.ylabel(r'$\log_{2}|\Delta|$')
# plt.legend()
# plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
# plt.tight_layout()
# # print(t_original[time])
# plt.show()


# norm1_L2 = computel2norm_single(convergence_list[0][:]-avg_high_res[-1][:],dx_list[0])
# norm2_L2 = computel2norm_single(convergence_list[1][:]-avg_high_res[-2][:],dx_list[1])
# norm3_L2 = computel2norm_single(convergence_list[2][:]-avg_high_res[-3][:],dx_list[2])
# norm4_L2 = computel2norm_single(convergence_list[3][:]-avg_high_res[-4][:],dx_list[3])
# L2_errors = [norm1_L2,norm2_L2,norm3_L2,norm4_L2]
# L2_order = np.log2([0.0, L2_errors[0]/L2_errors[1], L2_errors[1]/L2_errors[2], L2_errors[2]/L2_errors[3]])
# L2_order[0] = 0.0


# norm1_L1 = computel1norm_single(convergence_list[0][:]-avg_high_res[-1][:],dx_list[0])
# norm2_L1 = computel1norm_single(convergence_list[1][:]-avg_high_res[-2][:],dx_list[1])
# norm3_L1 = computel1norm_single(convergence_list[2][:]-avg_high_res[-3][:],dx_list[2])
# norm4_L1 = computel1norm_single(convergence_list[3][:]-avg_high_res[-4][:],dx_list[3])
# L1_errors = [norm1_L1,norm2_L1,norm3_L1,norm4_L1]
# L1_order = np.log2([0.0, L1_errors[0]/L1_errors[1], L1_errors[1]/L1_errors[2], L1_errors[2]/L1_errors[3]])
# L1_order[0] = 0.0

# norm1_Linf = computelinfnorm_single(convergence_list[0][:]-avg_high_res[-1][:])
# norm2_Linf = computelinfnorm_single(convergence_list[1][:]-avg_high_res[-2][:])
# norm3_Linf = computelinfnorm_single(convergence_list[2][:]-avg_high_res[-3][:])
# norm4_Linf = computelinfnorm_single(convergence_list[3][:]-avg_high_res[-4][:])
# Linf_errors = [norm1_Linf,norm2_Linf,norm3_Linf,norm4_Linf]
# Linf_order = np.log2([0.0, Linf_errors[0]/Linf_errors[1], Linf_errors[1]/Linf_errors[2], Linf_errors[2]/Linf_errors[3]])
# Linf_order[0] = 0.0

# # Print error table 
# print('===============================================================================')
# print('  N     L_2 Error     Order     L_1 Error     Order       L_Inf Error    Order ')
# print('===============================================================================')
# for k in range(len(file_names)-1):
#     print('{:6d}   {:.4e}  {:.4e}   {:.4e}  {:.4e}   {:.4e}  {:.4e} '.format(N_list[k], L2_errors[k], L2_order[k], L1_errors[k], L1_order[k], Linf_errors[k], Linf_order[k]))
# print('===============================================================================')




#### Convergence for exact TOV Initial Data
norm1_L2 = computel2norm_single(convergence_list[0][:]-ID_list[0][:],dx_list[0])
norm2_L2 = computel2norm_single(convergence_list[1][:]-ID_list[1][:],dx_list[1])
norm3_L2 = computel2norm_single(convergence_list[2][:]-ID_list[2][:],dx_list[2])
norm4_L2 = computel2norm_single(convergence_list[3][:]-ID_list[3][:],dx_list[3])
L2_errors = [norm1_L2,norm2_L2,norm3_L2,norm4_L2]
L2_order = np.log2([0.0, L2_errors[0]/L2_errors[1], L2_errors[1]/L2_errors[2], L2_errors[2]/L2_errors[3]])
L2_order[0] = 0.0


norm1_L1 = computel1norm_single(convergence_list[0][:]-ID_list[0][:],dx_list[0])
norm2_L1 = computel1norm_single(convergence_list[1][:]-ID_list[1][:],dx_list[1])
norm3_L1 = computel1norm_single(convergence_list[2][:]-ID_list[2][:],dx_list[2])
norm4_L1 = computel1norm_single(convergence_list[3][:]-ID_list[3][:],dx_list[3])
L1_errors = [norm1_L1,norm2_L1,norm3_L1,norm4_L1]
L1_order = np.log2([0.0, L1_errors[0]/L1_errors[1], L1_errors[1]/L1_errors[2], L1_errors[2]/L1_errors[3]])
L1_order[0] = 0.0

norm1_Linf = computelinfnorm_single(convergence_list[0][:]-ID_list[0][:])
norm2_Linf = computelinfnorm_single(convergence_list[1][:]-ID_list[1][:])
norm3_Linf = computelinfnorm_single(convergence_list[2][:]-ID_list[2][:])
norm4_Linf = computelinfnorm_single(convergence_list[3][:]-ID_list[3][:])
Linf_errors = [norm1_Linf,norm2_Linf,norm3_Linf,norm4_Linf]
Linf_order = np.log2([0.0, Linf_errors[0]/Linf_errors[1], Linf_errors[1]/Linf_errors[2], Linf_errors[2]/Linf_errors[3]])
Linf_order[0] = 0.0

# Print error table 
print('===============================================================================')
print('  N     L_2 Error     Order     L_1 Error     Order       L_Inf Error    Order ')
print('===============================================================================')
for k in range(len(file_names)):
    print('{:6d}   {:.4e}  {:.4e}   {:.4e}  {:.4e}   {:.4e}  {:.4e} '.format(N_list[k], L2_errors[k], L2_order[k], L1_errors[k], L1_order[k], Linf_errors[k], Linf_order[k]))
print('===============================================================================')

