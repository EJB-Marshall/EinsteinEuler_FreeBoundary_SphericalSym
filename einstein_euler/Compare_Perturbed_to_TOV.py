import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp
import h5py



#######################################
# Norms etc.
#######################################


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




################################################
# Compare perturbed solution to static TOV
################################################


#-----------------------------------
# TOV Equations
#-----------------------------------

def TOV_Star(r,y,K):

    m, L, rho = y

    drm = (r**2)/2*rho

    drL = (m + 0.5*(r**3)*(K*rho-1))/(r*(r - 2*m))

    drrho = (1/K)*(1 - (K+1)*rho)*drL

    return [drm, drL, drrho]


### Event Function to detect when the star boundary is reached (i.e. pressure = 0)
def Star_Boundary(r,y,K):

    m, L, rho = y

    return rho - 1/K

Star_Boundary.terminal = True # Integration should finish when the star boundary is reached

#-----------------------------------
# Series Expansion Near Origin
#-----------------------------------

def TOV_Series(r,rho0,c0,K):

    m_series = (r**3*rho0)/6. + (r**5*(-1 + (1 + K)*rho0)*(-3 + (-1 + 3*K)*rho0))/(120.*K)

    L_series = (r**2*(-3 + (1 + 3*K)*rho0))/12. - (r**4*(-3 + (1 + 3*K)*rho0)*(-3*(1 + 5*K) + (3 + K*(-2 + 15*K))*rho0))/(1440.*K) + \
                (r**6*(-3 + (1 + 3*K)*rho0)*(9*(3 + K*(-2 + 35*K)) + rho0*(-39 + 3*K*(16 + (47 - 210*K)*K) + (12 + K*(-45 + K*(65 + 3*K*(-41 + 105*K))))*rho0)))/(181440.*K**2) + c0

    rho_series = rho0 - (r**2*(-1 + (1 + K)*rho0)*(-3 + (1 + 3*K)*rho0))/(12.*K) + \
                (r**4*(-1 + (1 + K)*rho0)*(-3 + (1 + 3*K)*rho0)*(-3*(3 + 5*K) + (4 + 3*K*(3 + 5*K))*rho0))/(720.*K**2) - \
                (r**6*(-1 + (1 + K)*rho0)*(-3 + (1 + 3*K)*rho0)*(9*(31 + 3*K*(32 + 35*K)) - 18*(15 + K*(44 + 3*K*(32 + 35*K)))*rho0 + \
                (61 + K*(200 + 3*K*(206 + 9*K*(32 + 35*K))))*rho0**2))/(181440.*K**3)


    return [m_series, L_series, rho_series]



#-----------------------------------------------
# Solve ODEs using Scipy, Save Initial Data
#----------------------------------------------

def fmt(x):
    """Format a float for use in a filename, e.g. 1.5 -> 1p5"""
    return f"{x:.2f}".replace(".", "p")

def Construct_Initial_Data(x_start,x_end,Nx,Ngz,rho0,K):

    dx = (x_end-x_start)/Nx
    x_start_grid = x_start + (0.5 - Ngz)*dx
    x_end_grid = x_end + (Ngz - 0.5)*dx
    grid = np.linspace(x_start_grid, x_end_grid, Nx + 2*Ngz)

    ### Arrays for Initial Data
    m_ID = np.zeros_like(grid)
    L_ID = np.zeros_like(grid)
    a_ID = np.zeros_like(grid)
    alpha_ID = np.zeros_like(grid)
    rho_ID = np.zeros_like(grid)

    ### Solve ODEs
    r0 = grid[Ngz]
    c0 = 1 # Parameter for L ODE
    yinit = TOV_Series(r0,rho0,c0,K) # Initial data for ODE at r=r0

    soln = solve_ivp(TOV_Star, [r0, x_end], yinit, args=(K,),events=Star_Boundary,dense_output=True, method="DOP853", rtol=1e-12, atol=1e-14)


    ### Set up initial data for saving
    star_bdy_r = soln.t_events[0][0]
    m_star_bdy = soln.y_events[0][0][0]
    L_star_bdy = soln.y_events[0][0][1]
    star_mask = (x_start < grid) & (grid < star_bdy_r)
    vacuum_mask = (star_bdy_r < grid) & (grid< x_end)

    soln_vec = soln.sol(grid[star_mask])

    ###--------- Fill ID arrays -------------
    rho_ID[star_mask] = soln_vec[2,:]
    m_ID[star_mask] = soln_vec[0,:]
    m_ID[vacuum_mask] = m_star_bdy # Match to Schwarzschild solution at bdy and beyond
    a_ID[Ngz:-Ngz] = np.sqrt(1/(1 - 2*m_ID[Ngz:-Ngz] /grid[Ngz:-Ngz]))
    # print(2*m_ID[Ngz]/grid[Ngz])
   

    ### Shift L value to match with Schwarzschild solution
    L_ID[star_mask] = soln_vec[1,:]
    L_ID[vacuum_mask] = L_star_bdy
    L_ID[star_mask] = L_ID[star_mask] + (0.5*np.log(1 - 2*m_star_bdy/star_bdy_r) - L_star_bdy)
    L_ID[vacuum_mask] = 0.5*np.log(1 - 2*m_star_bdy/grid[vacuum_mask]) 
    alpha_ID = np.exp(L_ID)


    return np.array([rho_ID,a_ID]), star_bdy_r



#-------------------- Load HDF5 File -----------------------------

# filenames = ['LongEvolution_t1e4_N3200_K1p0_CD1p1_eps1e-3.hdf5','LongEvolution_t1e4_N6400_K1p0_CD1p1_eps1e-3.hdf5']
# filenames = ['test_N800_K1p0_CD1p1_eps1e-3.hdf5','test_N1600_K1p0_CD1p1_eps1e-3.hdf5','test_N3200_K1p0_CD1p1_eps1e-3.hdf5','test_N6400_K1p0_CD1p1_eps1e-3.hdf5']
filenames = ['test_N800_K1p0_CD1p1_Static.hdf5','test_N1600_K1p0_CD1p1_Static.hdf5','test_N3200_K1p0_CD1p1_Static.hdf5','test_N6400_K1p0_CD1p1_Static.hdf5']
rho_TOV_CD_vals = [1.1,1.1,1.1,1.1]


soln_list = []
dx_list = []
coords_list = []

for k in range(len(filenames)):
    file_name = './HDF_Files/'+filenames[k]
    Output = h5py.File(file_name, 'r')
    t = Output['Time'][:]
    boundaries = Output['Fluid_boundaries'][:,:]
    soln = Output['Solution/Solution_Array'][:,:,:]
    coordinates = Output['x_coordinates'][:]
    dx = np.abs(coordinates[1]-coordinates[0])
    N = np.shape(soln)[2]
    K = Output['Aux'][0]



    #--------- Create background TOV solution -----------------------
    Ngz = 3
    rho_TOV_CD = rho_TOV_CD_vals[k]
    TOV_soln, TOV_bdy = Construct_Initial_Data(0,5,N,Ngz,rho_TOV_CD,K)
    rho_TOV = TOV_soln[0,Ngz:-Ngz]
    a_TOV = TOV_soln[1,Ngz:-Ngz]




    ###################################
    # Recover Primitive Variables
    ###################################
    T00 = soln[:,0,:]
    T01 = soln[:,1,:]
    a = soln[:,2,:]
    v = np.zeros_like(T00)
    mu = np.zeros_like(T00)
    Gamma2 = np.zeros_like(T00)
    p_old = np.ones_like(mu[0,:])
    mu_CD = np.zeros_like(t)
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


        # Lagrange Interpolation to obtain boundary values
        l0 = -1.0/(6.0*dx**3)*(0.0-coordinates[1])*(0.0-coordinates[2])*(0.0-coordinates[3])
        l1 = 1.0/(2.0*dx**3)*(0.0-coordinates[0])*(0.0-coordinates[2])*(0.0-coordinates[3])
        l2 = -1.0/(2.0*dx**3)*(0.0-coordinates[0])*(0.0-coordinates[1])*(0.0-coordinates[3])
        l3 = 1.0/(6.0*dx**3)*(0.0-coordinates[0])*(0.0-coordinates[1])*(0.0-coordinates[2])

        mu_CD[i] = mu[i,0]*l0 + mu[i,1]*l1  + mu[i,2]*l2  + mu[i,3]*l3 

    if k == 3:
            avg_high_res = restrict_grid(mu,3)
    
    if k < 3:
        soln_list.append(mu)
        dx_list.append(dx)
        coords_list.append(coordinates)

    ##############################################
    # Plotting Code
    ##############################################

    plt.rcParams.update({"text.usetex": True,
        "font.family": "serif",
        "font.serif": "Computer Modern",
        "savefig.bbox": "tight",
        "savefig.format": "pdf"})
    plt.rc('font', size=16)

    l = -1
    print(t[l])
    plt.plot(coordinates,np.log2(np.abs(mu[l,:]-rho_TOV)),label=f'$N = {N}$')
    plt.xlabel(r'$r$')
    plt.ylabel(r'$\log_{2}\Big(|\rho_{N} - \rho_{TOV}|\Big)$')
    plt.legend()
    plt.grid(alpha=0.8)

    # plt.figure(1)
    # plt.plot(t,mu_CD,'-',label=f'$\\rho_{{0}} = {rho_TOV_CD}$')
    # plt.axhline(rho_TOV_CD,color='tab:red')#,label=r'$\rho^{TOV}_{0}$')
    # # plt.title(f'$N={N}$, $K={K}$, $\\rho^{{TOV}}_{{0}} = {rho_TOV_CD}$')
    # plt.ylabel(r'$\rho(0,t)$')
    # plt.xlabel(r'$t$')
    # plt.grid(alpha=0.8)
    # plt.legend()


    # plt.figure(2)
    # plt.plot(t,boundaries[:,1],'-',label=f'$N={N}$')
    # plt.axhline(TOV_bdy,color='tab:red')#,label=r'$x_{R}^{TOV}$')
    # # plt.title(f'$N={N}$, $K={K}$, $\\rho^{{TOV}}_{{0}} = {rho_TOV_CD}$')
    # plt.ylabel(r'$x_{R}$')
    # plt.xlabel(r'$t$')
    # plt.grid(alpha=0.8)
    # plt.legend()


    # plt.figure(3)
    # plt.plot(t,computel2norm(v,dx),'-',label=f'$\\rho_{{0}} = {rho_TOV_CD}$')
    # # plt.axhline(TOV_bdy,color='tab:red')#,label=r'$x_{R}^{TOV}$')
    # # plt.title(f'$N={N}$, $K={K}$, $\\rho^{{TOV}}_{{0}} = {rho_TOV_CD}$')
    # plt.ylabel(r'$\| v \|_{L^{2}}$')
    # plt.xlabel(r'$t$')
    # plt.grid(alpha=0.8)
    # plt.legend()


# l = 5
# plt.plot(coords_list[0],np.log2(np.abs(soln_list[0][l,:]-avg_high_res[3][l,:])),label=f'$N = 800$')
# plt.plot(coords_list[1],np.log2(np.abs(soln_list[1][l,:]-avg_high_res[2][l,:])),label=f'$N = 1600$')
# plt.plot(coords_list[2],np.log2(np.abs(soln_list[2][l,:]-avg_high_res[1][l,:])),label=f'$N = 3200$')
# plt.title(f't={t[l]}')
# plt.xlabel(r'$r$')
# plt.ylabel(r'$\log_{2}\Big(|\Delta|\Big)$')
# plt.legend()
# plt.grid(alpha=0.8)
plt.show()



