import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp



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

    m_series = (r**3*rho0)/6. - (r**5*(-1 + (1 + K)*rho0)*(-3 + (1 + 3*K)*rho0))/(120.*K)

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


    # plt.plot(grid[Ngz:-Ngz],a_ID[Ngz:-Ngz])
    # plt.plot(grid[197],rho_ID[197],'x')
    # plt.axvline(star_bdy_r,color='tab:green')
    # plt.show()

    ### Save initial data to file:
    rho_ID.tofile(f"rho_ID.bin")
    # alpha_ID.tofile(f"alpha_ID.bin")
    a_ID.tofile(f"a_ID.bin")
    star_bdy_r.tofile(f"star_bdy_ID.bin")


    return


#------------------------------------------------
# Introduce Velocity Perturbation to TOV Star
#------------------------------------------------

# Initial Fluid Velocity Perturbation
def v_perturbation(r,star_bdy_r,eps):

    offset = 0.1

    # Smooth bump function
    peak = (star_bdy_r-offset + offset)/2
    # v = eps*np.exp(-1/((r-offset)*(star_bdy_r-offset -r)))/(np.exp(-1/((peak-offset)*(star_bdy_r-offset - peak))))
    if r<=offset or r>=star_bdy_r-offset:
        v = 0.0
    else:
        v = eps*np.exp(-1/((r-offset)*(star_bdy_r-offset -r)))/(np.exp(-1/((peak-offset)*(star_bdy_r-offset - peak))))

    return v

def m_ODE_velocity(r,y,TOV_soln,v_perturbation,K,star_bdy_r,eps):

    # Velocity Perturbation
    v = v_perturbation(r,star_bdy_r,eps)
    rho = TOV_soln(r)[2]

    Gamma2 = 1/(1-v**2)

    drm = (r**2)/2*(1 - K*rho + Gamma2*(-1 + (1+K)*rho))

    return [drm]


def Construct_Perturbed_Initial_Data(x_start,x_end,Nx,Ngz,rho0,K,eps):

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
    # plt.plot(grid[Ngz:-Ngz],m_ID[Ngz:-Ngz],label='TOV')
   

    ### Shift L value to match with Schwarzschild solution
    L_ID[star_mask] = soln_vec[1,:]
    L_ID[vacuum_mask] = L_star_bdy
    L_ID[star_mask] = L_ID[star_mask] + (0.5*np.log(1 - 2*m_star_bdy/star_bdy_r) - L_star_bdy)
    L_ID[vacuum_mask] = 0.5*np.log(1 - 2*m_star_bdy/grid[vacuum_mask]) 
    alpha_ID = np.exp(L_ID)

    
    # Velocity Perturbation
    v_ID = np.zeros_like(grid)
    for i in range(np.shape(grid)[0]):
        v_ID[i] = v_perturbation(grid[i],star_bdy_r,eps)

    # Re-solve constraint for m
    soln_perturbed_m = solve_ivp(m_ODE_velocity, [r0, star_bdy_r], [yinit[0]], args=(soln.sol,v_perturbation,K,star_bdy_r,eps),dense_output=True, method="DOP853", rtol=1e-12, atol=1e-14)
    m_ID[star_mask] = soln_perturbed_m.sol(grid[star_mask])[0]
    m_ID[vacuum_mask] = soln_perturbed_m.sol(star_bdy_r)[0]# Match to Schwarzschild solution at bdy and beyond
    a_ID[Ngz:-Ngz] = np.sqrt(1/(1 - 2*m_ID[Ngz:-Ngz] /grid[Ngz:-Ngz]))


    # plt.plot(grid[Ngz:-Ngz],v_ID[Ngz:-Ngz],label='perturbed')
    # # plt.plot(grid[197],rho_ID[197],'x')
    # plt.axvline(star_bdy_r,color='tab:green')
    # plt.show()

    ### Save initial data to file:
    rho_ID.tofile(f"rho_perturbed_ID.bin")
    # alpha_ID.tofile(f"alpha_perturbed_ID.bin")
    a_ID.tofile(f"a_perturbed_ID.bin")
    v_ID.tofile(f"v_perturbed_ID.bin")
    star_bdy_r.tofile(f"star_bdy_perturbed_ID.bin")


    return star_bdy_r


