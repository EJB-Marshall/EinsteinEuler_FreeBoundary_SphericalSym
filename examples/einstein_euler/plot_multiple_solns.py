import numpy as np
from matplotlib import pyplot as plt
import h5py



def restrict_grid(soln,num_grids):
    """ In a cell-centred scheme, the cell centres will not be aligned
        for different resolutions. Thus, we must average the cell-centred values from
         high resolution runs to align with the low resolution cells. """
    
    n = soln.shape[1] # Number of grid points for high res run
    averaged_runs = []

    averaged_runs.append(soln)
    t_shape = np.shape(soln)[0]

    for k in range(1,num_grids+1):

        averaged_runs.append(np.zeros((int(t_shape),int(n/(2**(k))))))


        for parent_index in range(int(n/(2**(k)))):
            child_index = 2 * parent_index

            ### Conservative update
            averaged_runs[k][:,parent_index] = 0.5*(averaged_runs[k-1][:,child_index] + averaged_runs[k-1][:,child_index+1])

    return averaged_runs

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


def compute_integral_abs_single(f,dx):
    int_eval = np.abs(dx*np.sum(f[:]))
    return int_eval


def integral_cutcell(tau, r, xR, dx):
    """Computes ∫_0^{xR} tau r^2 dr with cut cell near surface"""
    iR = np.searchsorted(r, xR) - 1          # last cell centre inside the star
    rl = r[:iR+1] - 0.5*dx
    rr = r[:iR+1] + 0.5*dx
    rr[-1] = xR                              # surface cell ends at x_R
    return np.sum(tau[:iR+1]*(rr**3 - rl**3)/3.0)


# filenames = ['test_N800_K1p0_CD1p1_eps1e-3.hdf5','test_N1600_K1p0_CD1p1_eps1e-3.hdf5','test_N3200_K1p0_CD1p1_eps1e-3.hdf5','test_N6400_K1p0_CD1p1_eps1e-3.hdf5']

filenames = ['PerturbTest_N3200_K1p0_CD1p1_eps1e-3.hdf5','PerturbTest_N3200_K1p0_CD1p1_eps5e-3.hdf5','PerturbTest_N3200_K1p0_CD1p1_eps1e-2.hdf5']
eps_vals = [1e-3,5e-3,1e-2]
colours = ['tab:blue','tab:orange','tab:green']
soln_list = []
dx_list = []
boundary_list = []
coords_list = []
CD_list = []

plt.rcParams.update({"text.usetex": True,
        "font.family": "serif",
        "font.serif": "Computer Modern",
        "savefig.bbox": "tight",
        "savefig.format": "pdf"})
plt.rc('font', size=16)


for k in range(0,np.shape(filenames)[0]):
    file_name = './HDF_Files/'+filenames[k]
    Output = h5py.File(file_name, 'r')
    t = Output['Time'][:]
    boundaries = Output['Fluid_boundaries'][:,:]
    soln = Output['Solution/Solution_Array'][:,:,:]
    coordinates = Output['x_coordinates'][:]
    dx = np.abs(coordinates[1]-coordinates[0])
    N = np.shape(soln)[2]
    K = Output['Aux'][0]


    ###################################
    # Recover Primitive Variables
    ###################################
    T00 = soln[:,0,:]
    T01 = soln[:,1,:]
    a = soln[:,2,:]
    v = np.zeros_like(T00)
    mu = np.zeros_like(T00)
    m = np.zeros_like(T00)
    mass = np.zeros_like(t)
    Gamma2 = np.zeros_like(T00)
    p_old = np.ones_like(mu[0,:])
    mu_CD = np.zeros_like(t)
    Cons_Error = np.zeros_like(t)
    cut_cell_integral = np.zeros_like(t)
    # for i in range(np.shape(t)[0]):

    #     for j in range(np.shape(coordinates)[0]):

    #         if T00[i,j]>0:
    #             Gamma2_star = 1/(1 - (T01[i,j]/(T00[i,j] + p_old[j]))**2)
    #             mu_star = (T00[i,j]+p_old[j]*(1-Gamma2_star))/Gamma2_star

    #             p_fun =  p_old[j] - (K*mu_star - 1)

    #             dp_fun = 1 - K*T01[i,j]/((T00[i,j] + p_old[j])**2)

    #             p_new = p_old[j] - p_fun/dp_fun

    #             v[i,j] = T01[i,j]/(T00[i,j] + p_new)

    #             Gamma2[i,j] = 1/(1-v[i,j]**2)

    #             mu[i,j] = (T00[i,j] +p_new*(1-Gamma2[i,j]))/Gamma2[i,j]

    #             p_old[j] = p_new

    #-------- Explicit Analytic Formula
    for i in range(np.shape(t)[0]):
        mask1 = coordinates>boundaries[i,0]
        mask2 = coordinates<boundaries[i,1]
        mask = mask1 & mask2

        mu[i,mask] = (1.0 + (K - 1.0)*T00[i,mask] + np.sqrt((-1.0 + T00[i,mask] + K*T00[i,mask])**2 - 4.0*K*T01[i,mask]**2))/(2.0*K)

        v[i,mask] = ((2.0*T01[i,mask]))/(-1.0 + (K + 1.0)*T00[i,mask] + np.sqrt((-1.0 + (K+1.0)*T00[i,mask])**2 - 4.0*K*T01[i,mask]**2))

        m = (coordinates/2)*(1 - 1/(a[i,:]**2))

        mass[i] = m[mask][-1]

        # Lagrange Interpolation to obtain boundary values
        l0 = -1.0/(6.0*dx**3)*(0.0-coordinates[1])*(0.0-coordinates[2])*(0.0-coordinates[3])
        l1 = 1.0/(2.0*dx**3)*(0.0-coordinates[0])*(0.0-coordinates[2])*(0.0-coordinates[3])
        l2 = -1.0/(2.0*dx**3)*(0.0-coordinates[0])*(0.0-coordinates[1])*(0.0-coordinates[3])
        l3 = 1.0/(6.0*dx**3)*(0.0-coordinates[0])*(0.0-coordinates[1])*(0.0-coordinates[2])

        mu_CD[i] = mu[i,0]*l0 + mu[i,1]*l1  + mu[i,2]*l2  + mu[i,3]*l3 

        cut_cell_integral[i] = integral_cutcell(T00[i,:], coordinates, boundaries[i,1], dx)
        
    Cons_Error = np.abs(cut_cell_integral-cut_cell_integral[0])

    if k == 0:
            avg_high_res = restrict_grid(mu,2)

    if k != 0:
        soln_list.append(mu)
        dx_list.append(dx)
        coords_list.append(coordinates)

    boundary_list.append(boundaries)
    CD_list.append(mu_CD)

    i = 500
    mask1 = coordinates>boundaries[i,0]
    mask2 = coordinates<boundaries[i,1]
    mask = mask1 & mask2
    plt.plot(coordinates[mask],mu[i,mask],color=colours[k],label=f'$\\epsilon = {eps_vals[k]}$')
    # plt.title(f'$t = {t[i]:.2f}$')
    

    # plt.plot(t,np.log2(np.max(np.abs(mu[:,:]-mu[0,:]),axis=1)),label=f'N={N}')
    # plt.plot(t,np.log2(np.abs(boundaries[:,1]-boundaries[0,1])),label=f'N={N}')
    # plt.plot(t,boundaries[:,1],label=f'N={N}')
    # plt.plot(t,np.log2(Cons_Error[:]),label=f'N={N}')


# plt.plot(t,np.log2(computel2norm(soln_list[1]-avg_high_res[-1],dx_list[1])/computel2norm(soln_list[0]-avg_high_res[-2],dx_list[0])))
# plt.plot(t,np.log2(e_c/e_m))
# plt.plot(t,np.log2(np.abs(boundary_list[-1][:,1]-boundary_list[0][:,1])/np.abs(boundary_list[-2][:,1]-boundary_list[0][:,1])))
# plt.plot(t,np.log2(np.abs(CD_list[-2][:]-CD_list[0][:]))+2.0,label=r'$N=3200$, shifted')
# plt.plot(t,np.log2(np.abs(CD_list[-1][:]-CD_list[0][:])),label=r'$N=1600$')
# plt.plot(t,np.log2(np.abs(boundary_list[-2][:,1]-boundary_list[0][:,1]))+2.0,label=r'$N=3200$, shifted')
# plt.plot(t,np.log2(np.abs(boundary_list[-1][:,1]-boundary_list[0][:,1])),label=r'$N=1600$')
# plt.xlabel(r'$t$')
# # plt.ylabel(r'$\log_{2}\Big(|\rho_{N}(t,0)-\rho_{6400}(t,0)|\Big)$')
# plt.ylabel(r'$\log_{2}\Big(|x^{N}_{R}(t)-x_{R}^{6400}(t)|\Big)$')
# plt.ylabel(r'$|\int_{\mathcal{D}} \tau r^{2} \; dr \big(t\big) -  \int_{\mathcal{D}} \tau r^{2} \; dr \big(0\big) |$')
# plt.axvline(boundaries[i,1],alpha=1.0,color='grey',linestyle='--')
plt.ylabel(r'$\rho$')
plt.xlabel(r'$r$')
plt.legend()
plt.show()