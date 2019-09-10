import numpy as np
import matplotlib.pyplot as plt

#Function defining initial conditions
def initialBell(x):
	return np.where(x%1. < 0.5 , np.power( np.sin( 2*x*np.pi), 2), 0)

def burger_adv_FTCS(parameters, total_time):
#Parameters: space, time, initial profile
    nx = parameters['nx']           #No. of space steps
    dx = 1./nx         #space steps (m)
    nt = parameters['nt']          #No. of time steps
    dt = total_time/nt          #Time step in seconds
       
    x = np.linspace(0.0, 1.0, nx+1)
    uadv = initialBell(x)
    uadv_new = u.copy()
    uadv_old = u.copy()
    
    ucons = initialBell(x)
    ucons_new = u.copy()
    ucons_old = u.copy()   
    #FTCS
    #For advective form     
    for n in range(1,nt):
        for j in range(1,nx):
            uadv_new[j] = uadv[j]*( 1 + 2*dt/dx * ( uadv[j+1] - uadv[j-1] ) )
            
            #Apply bound conditions
            #Calculate u[0] by using the penultimate u as u[-1]
            uadv_new[0] = uadv[0]*( 1 + 2*dt/dx * ( uadv[1] - uadv[nx-1] ) )
            
            #update u for new timestep
            uadv_old = uadv.copy() 
            uadv = uadv_new.copy()
            
            
      #Conservative form      
      for n in range(1,nt):
        for j in range(1,nx):
            ucons_new[j] = ucons[j]*( 1 + 0.25*dt/dx * ( ucons[j+1]**2 - ucons[j-1]**2 ) )
            
            #Apply bound conditions
            #Calculate u[0] by using the penultimate u as u[-1]
            ucons_new[0] = ucons[0]*( 1 + 0.25*dt/dx * ( ucons[1]**2 - ucons[nx-1]**2 ) )
            
            #update u for new timestep
            ucons_old = ucons.copy() 
            ucons = ucons_new.copy()
        
    return x, uadv, ucons     

def burger_cons_FTBS(parameters, total_time):
#Parameters: space, time, initial profile
        nx = parameters['nx']           #No. of space steps
    dx = 1./nx         #space steps (m)
    nt = parameters['nt']          #No. of time steps
    dt = total_time/nt          #Time step in seconds
       
nx = 100            #No. of space steps
    dx = 1./nx         #space steps (m)
    nt = 10000           #No. of time steps
    dt = total_time/nt          #Time step in seconds
       
    x = np.linspace(0.0, 1.0, nx+1)
    u = initialBell(x)
    u_new = u.copy()
    u_old = u.copy()
        
    #FTBS
         
    for n in range(1,nt):
        for j in range(1,nx):
            u_new[j] = u[j]*( 1 + 0.5*dt/dx * ( u[j]**2 - u[j-1]**2 ) )
            
            #Apply bound conditions
            #Calculate u[0] by using the penultimate u as u[-1]
            u_new[0] = u[0]*( 1 + 0.5*dt/dx * ( u[0]**2 - u[nx-1]**2 ) )
            
            #update u for new timestep
            u_old = u.copy() 
            u = u_new.copy()
            
    return x, u


def burger_adv_FTBS(total_time):
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
        for j in range(1,nx):
            u_new[j] = u[j]*( 1 + dt/dx * ( u[j] - u[j-1] ) )
            
            #Apply bound conditions
            #Calculate u[0] by using the penultimate u as u[-1]
            u_new[0] = u[0]*( 1 + dt/dx * ( u[0] - u[nx-1] ) )
            
            #update u for new timestep
            u_old = u.copy() 
            u = u_new.copy()
            
    return x, u


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
            
    return x, u



#Parameters
para = {'nx': 100, 'nt': 1e4}

nx = 100            #No. of space steps
dx = 1./nx         #space steps (m)
nt = 10000           #No. of time steps
dt = total_time/nt          #Time step in seconds
   

#Execute code
total_time = 0.07
x_CSadv, u_CSadv =  burger_adv_FTCS(total_time) 
x_CScons, u_CScons = burger_cons_FTCS(total_time) 
x_BSadv, u_BSadv = burger_adv_FTBS(total_time)
x_BScons, u_BScons = burger_cons_FTBS(total_time)

#plot results
plt.plot(x_CSadv, u_CSadv, 'yellow', label='adv form, FTCS')
plt.plot(x_CScons, u_CScons, 'green', label='cons form, FTCS')
plt.plot(x_BSadv, u_BSadv, 'red', label='adv form, FTBS')
plt.plot(x_BScons, u_BScons, 'blue', label='cons form, FTBS')
plt.legend()
plt.ylabel('u (ms-1)')
time = '{} seconds'.format(total_time)
plt.text(0.5,1,time)
plt.show()
    
    
    
    
