import numpy as np
import matplotlib.pyplot as plt

#Function defining initial conditions
def initialPosBell(x):
	return np.where(x%1. < 0.5 , np.power( np.sin( 2*x*np.pi), 2), 0)

def initialNegBell(x):
	return np.where(x%1. < 0.5 , -1 * np.power( np.sin( 2*x*np.pi), 2), 0)
	
def burger_FTCS(parameters, total_time, initial_cond):
#Parameters: space, time, initial profile
    nx = parameters['nx']           #No. of space steps
    dx = 1./nx         #space steps (m)
    nt = parameters['nt']          #No. of time steps
    dt = total_time/nt          #Time step in seconds
       
    uadv = initial_cond
    uadv_new = uadv.copy()
    uadv_old = uadv.copy()
    
    ucons = initial_cond
    ucons_new = ucons.copy()
    ucons_old = ucons.copy()   
    #FTCS
    #For advective form     
    for n in range(1,nt):
        for j in range(1,nx):
            uadv_new[j] = uadv[j]*( 1 - 0.5*dt/dx * ( uadv[j+1] - uadv[j-1] ) )
            
            #Apply bound conditions
            #Calculate u[0] by using the penultimate u as u[-1]
            uadv_new[0] = uadv[0]*( 1 - 0.5*dt/dx * ( uadv[1] - uadv[nx-1] ) )
            
            #update u for new timestep
            uadv_old = uadv.copy() 
            uadv = uadv_new.copy()
            
            
      #Conservative form      
    for n in range(1,nt):
        for j in range(1,nx):
            ucons_new[j] = ucons[j]*( 1 - 0.125*dt/dx * ( ucons[j+1]**2 - ucons[j-1]**2 ) )
            
            #Apply bound conditions
            #Calculate u[0] by using the penultimate u as u[-1]
            ucons_new[0] = ucons[0]*( 1 - 0.125*dt/dx * ( ucons[1]**2 - ucons[nx-1]**2 ) )
            
            #update u for new timestep
            ucons_old = ucons.copy() 
            ucons = ucons_new.copy()
        
    return uadv, ucons     

def burger_FTBS(parameters, total_time, initial_cond):
#Parameters: space, time, initial profile
    nx = parameters['nx']           #No. of space steps
    dx = 1./nx         #space steps (m)
    nt = parameters['nt']          #No. of time steps
    dt = total_time/nt          #Time step in seconds
  
    uadv = initial_cond
    uadv_new = uadv.copy()
    uadv_old = uadv.copy()
     
    ucons = initial_cond
    ucons_new = ucons.copy()
    ucons_old = ucons.copy()   
    
    #FTBS
         
    for n in range(1,nt):
        for j in range(1,nx):
            uadv_new[j] = uadv[j]*( 1 - 0.5*dt/dx * ( uadv[j] - uadv[j-1] ) )
            
            #Apply bound conditions
            #Calculate u[0] by using the penultimate u as u[-1]
            uadv_new[0] = uadv[0]*( 1 - 0.5*dt/dx * ( uadv[0] - uadv[nx-1] ) )
            
            #update u for new timestep
            uadv_old = uadv.copy() 
            uadv = uadv_new.copy()
            
    for n in range(1,nt):
        for j in range(1,nx):
            ucons_new[j] = ucons[j]*( 1 - 0.125*dt/dx * ( ucons[j]**2 - ucons[j-1]**2 ) )
            
            #Apply bound conditions
            #Calculate u[0] by using the penultimate u as u[-1]
            ucons_new[0] = ucons[0]*( 1 - 0.125*dt/dx * ( ucons[0]**2 - ucons[nx-1]**2 ) )
            
            #update u for new timestep
            ucons_old = ucons.copy() 
            ucons = ucons_new.copy()
            
    return uadv, ucons


def burger_adv_FT_5points(parameters, total_time):
#Parameters: space, time, initial profile
    nx = parameters['nx']           #No. of space steps
    dx = 1./nx         #space steps (m)
    nt = parameters['nt']          #No. of time steps
    dt = total_time/nt          #Time step in seconds
       
    x = np.linspace(0.0, 1.0, nx+1)
    u = initialBell(x)
    u_new = u.copy()
    u_old = u.copy()
        
    #FTBS
         
    for n in range(1,nt):
        for j in range(2,nx):
            u_new[j] = u[j]*( 1 + dt/dx * ( u[j] - u[j] ) )
            
            #Apply bound conditions
            #Calculate u[0] by using the penultimate u as u[-1]
            u_new[0] = u[0]*( 1 + dt/dx * ( u[0] - u[nx-1] ) )
            
            #update u for new timestep
            u_old = u.copy() 
            u = u_new.copy()
            
    return u


def main():
    #Parameters
    param = {'nx': 70, 'nt': int(1e4)}

    #Execute code
    total_time = 0.40    
    x = np.linspace(0.0, 1.0, param['nx']+1)
    pos_cond = initialPosBell(x)
    neg_cond = initialNegBell(x)
    
    u_posCSadv , u_posCScons=  burger_FTCS(param, total_time, pos_cond) 
    u_posBSadv, u_posBScons= burger_FTBS(param, total_time, pos_cond)
    
    #u_negCSadv , u_negCScons=  burger_FTCS(param, total_time, neg_cond) 
    #u_negBSadv, u_negBScons= burger_FTBS(param, total_time, neg_cond)

    #plot results
    plt.figure(figsize=(12,5))

    #plt.plot(x, u_posCSadv, 'yellow', label='adv form, FTCS')
    plt.plot(x, u_posCScons, 'green', label='cons form, FTCS')
    plt.plot(x, u_posBSadv, 'red', label='adv form, FTBS')
    plt.plot(x, u_posBScons, 'blue', label='cons form, FTBS')
    plt.plot(x, pos_cond, 'grey', linestyle = '--', label='initial condition')
    plt.legend()
    plt.ylabel('u (ms-1)')
    time = 'Total time = {} seconds'.format(total_time)
    plt.text(0.5,1,time)
    plt.title('+ve u initial cond.')
    plt.show()
    
main()
    
    
    
    
