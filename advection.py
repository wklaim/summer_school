import numpy as np
import matplotlib.pyplot as plt

#Function defining initial conditions
def initialPosBell(x):
	return np.where(x%1. < 0.5 , np.power( np.sin( 2*x*np.pi), 2), 0)

def initialSin(x):
	return np.where(x%1. < 1 , np.sin( 2*x*np.pi), 0)

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
    print('Calc. FTCS') 
    
    uadv_arr = np.zeros((100,nx+1)) 
    ucons_arr = np.zeros((100,nx+1)) 
    #Outside loop save the u arrays every nt/100 timesteps for ploting a gif
    for n_arr in range(100):  
        
        #Inside loop calculates each imestep
        for n in range(1, int(nt/100) ):
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
                ucons_new[j] = ucons[j] - 0.25*dt/dx * ( ucons[j+1]**2 - ucons[j-1]**2 ) 
            
                #Apply bound conditions
                #Calculate u[0] by using the penultimate u as u[-1]
                ucons_new[0] = ucons[0] - 0.25*dt/dx * ( ucons[1]**2 - ucons[nx-1]**2 ) 
            
                #update u for new timestep
                ucons_old = ucons.copy() 
                ucons = ucons_new.copy()
        
        uadv_arr[n_arr,:] = uadv
        ucons_arr[n_arr,:] = ucons 
    return uadv_arr, ucons_arr     
        

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
    print('Calc. FTBS')     
    for n in range(1,nt):
        for j in range(1,nx):
            uadv_new[j] = uadv[j]*( 1 - dt/dx * ( uadv[j] - uadv[j-1] ) )
            
            #Apply bound conditions
            #Calculate u[0] by using the penultimate u as u[-1]
            uadv_new[0] = uadv[0]*( 1 - dt/dx * ( uadv[0] - uadv[nx-1] ) )
            
            #update u for new timestep
            uadv_old = uadv.copy() 
            uadv = uadv_new.copy()
            
    for n in range(1,nt):
        for j in range(1,nx):
            ucons_new[j] = ucons[j] - 0.5*dt/dx * ( ucons[j]**2 - ucons[j-1]**2 ) 
            
            #Apply bound conditions
            #Calculate u[0] by using the penultimate u as u[-1]
            ucons_new[0] = ucons[0] - 0.5*dt/dx * ( ucons[0]**2 - ucons[nx-1]**2 ) 
            
            #update u for new timestep
            ucons_old = ucons.copy() 
            ucons = ucons_new.copy()
            
    return uadv, ucons


def burger_FT_5points(parameters, total_time, initial_cond):
#Parameters: space, time, initial profile
    nx = parameters['nx']           #No. of space steps
    dx = 1./nx         #space steps (m)
    nt = parameters['nt']          #No. of time steps
    dt = total_time/nt          #Time step in seconds
       
    uadv = initial_cond.copy()
    uadv_new = uadv.copy()
    uadv_old = uadv.copy()
    
    ucons = initial_cond.copy()
    ucons_new = ucons.copy()
    ucons_old = ucons.copy()    
    
    #FTCS 5 point stencil
    print('Calc. FT 5 points') 
         
    for n in range(1,nt):
        for j in range(0,nx+1):
            uadv_new[(j)%nx] = uadv[(j)%nx] - uadv[(j)%nx]*dt/(12*dx) * \
            ( -uadv[(j+2)%nx] +8*uadv[(j+1)%nx] - 8*uadv[(j-1)%nx] + uadv[(j-2)%nx] ) 

            #update u for new timestep
            uadv_old = uadv.copy() 
            uadv= uadv_new.copy()
            
    for n in range(1,nt):
        for j in range(0,nx+1):
            ucons_new[j%nx] = ucons[j%nx] - dt/(24*dx) * \
            ( -ucons[(j+2)%nx]**2 +8*ucons[(j+1)%nx]**2 - 8*ucons[(j-1)%nx]**2 + ucons[(j-2)%nx]**2 ) 
            
            ucons_old = ucons.copy() 
            ucons = ucons_new.copy()     
    
    print('dt/dx =', dt/dx)
    return uadv, ucons





def main():
    #Parameters
    param = {'nx': 100, 'nt': int(1e4)}

    #Execute code
    total_time = 0.470   
    x = np.linspace(0.0, 1.0, param['nx']+1)
    pos_cond = initialPosBell(x)
    neg_cond = initialNegBell(x)
    
    u_posCSadv , u_posCScons=  burger_FTCS(param, total_time, pos_cond) 
    #u_posBSadv, u_posBScons= burger_FTBS(param, total_time, pos_cond)
    
    #u_negCSadv , u_negCScons=  burger_FTCS(param, total_time, neg_cond) 
    #u_negBSadv, u_negBScons= burger_FTBS(param, total_time, neg_cond)

    #u_pos_5pt_adv, u_pos_5pt_cons =  burger_FT_5points(param, total_time, pos_cond)         
    
    #plot results
    plt.figure(figsize=(12,5))

    plt.plot(x, u_posCSadv[-1], 'yellow', label='adv form, FTCS')
    plt.plot(x, u_posCScons[-1], 'green', label='cons form, FTCS')
    #plt.plot(x, u_posBSadv, 'red', label='adv form, FTBS')
    #plt.plot(x, u_posBScons, 'blue', label='cons form, FTBS')
    #plt.plot(x, u_pos_5pt_adv, 'grey', label='adv form, FT 5points')    
    #plt.plot(x, u_pos_5pt_cons, 'pink', label='cons form, FT 5points')
            
    plt.plot(x, pos_cond, 'grey', linestyle = '--', label='initial condition')
    plt.legend()
    plt.ylabel('u (ms-1)')
    plt.xlabel('x (m)')
    
    time = 'Total time = {} seconds'.format(total_time)
    plt.text(0,1.2,time)
    plt.title('+ve u initial cond.')
    #plt.savefig('burgers.pdf', bbox_syle = 'tight')
    plt.show()
    
main()
    
    
    
    
