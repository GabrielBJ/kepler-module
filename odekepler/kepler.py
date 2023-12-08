#!/usr/bin/env python3

# import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib.animation import FuncAnimation
import matplotlib 
from PIL import Image
import argparse
import requests
from io import BytesIO



# create a class that initialises the system

class TwoBodySimulation():
    '''
    This class initialises the system to be integrated. It takes the eccentricity and period of the orbit as arguments.

    Attributes:
    __________

    e       eccentricity of the orbit
    period  period of the orbit in years

    Methods:
    _______

    ics()   returns the initial conditions of the system

    '''
    def __init__(self, e, T):
        '''
        Initialises the system with the eccentricity and period of the orbit.
        '''
        self.e       = e
        self.period  = T
    

    def ics(self):
        '''
        Method that returns the initial conditions of the system given the eccentricity and period of the orbit.

        Returns:
        ________

        u       vector containing the initial conditions of the system (x0, y0, vx0, vy0) in SI units
        '''
        # initialise the initial conditions vector
        u = np.zeros(4)
        
        # compute the initial conditions
        u[0] = 0 
        u[1] = 1.496e11 * (1 - self.e)
        u[2] = -1 * np.sqrt(8.86893e8 * (1 + self.e)/(1 - self.e))
        u[3] = 0
        
        return u 


def slope(t, state):
    '''
    Function that computes the slope of the system given the state of the system and time.

    Arguments:
    _________

    t       time in seconds
    state   vector containing the state of the system (x, y, vx, vy) in SI units at time t

    Returns:
    _______

    du      vector containing the slope of the system (vx, vy, ax, ay) in SI units at time t
    '''

    # compute the distance from the sun
    r = np.sqrt(state[0]**2 + state[1]**2)

    # compute G * M_sun
    k = 1.32679e20 
    
    # initialise the slope vector
    du = np.zeros(4)
    
    # compute the slope
    for i in range(2):
        du[i]       = state[i + 2]
        du[i + 2]   = -k * state[i] * 1/r**3

    return du


def time_array(T, dt):
    '''
    Function that computes the time array given the period of the orbit and the time step.

    Arguments:
    _________

    T       period of the orbit in years
    dt      time step in days

    Returns:
    _______

    t       time array in seconds
    '''
    
    # transform the period and time step to seconds
    period  = T * 3.156e+7 # approx. the length of a year in seconds
    dt      = dt * 86411  # approx. the length of a day in seconds

    # compute the time array
    t = np.arange(0, period + dt, dt)

    return t


class RKIntegrate:
    '''
    This class contains the methods to integrate the system using the Runge-Kutta methods.

    Methods:
    _______

    RK2(system, t_step)     integrates the system using the RK2 method
    RK3(system, t_step)     integrates the system using the RK3 method
    RK4(system, t_step)     integrates the system using the RK4 method
    '''

    def RK2(system, t_step, save_dir = ''):
        '''
        Method that integrates the system using the RK2 method and saves the orbit in the provided directory. If no directory is provided, it will save the file in the current directory.

        Arguments:
        _________

        system      instance of the TwoBodySimulation class containing the system to be integrated
        t_step      time step in days
        save_dir    directory to save the integrated orbit (default is './RK2_integrated_orbit.txt')

        Returns:
        _______

        orbit       array containing the integrated orbit (x, y, vx, vy) in SI units
        '''
        print("Integrating system...") 
        # compute time array                                                                    
        t       = time_array(system.period, t_step)
        t_step  = t_step * 86411 # approx. the length of a day in seconds

        # compute initial conditions
        u = system.ics()

        # initialize orbit array
        orbit = np.zeros((len(t), 4))

        # set initial conditions
        orbit[0, :] = u.T 
        
        # integrate using RK2
        for i in range(0, len(t) - 1):
            
            # define the slopes
            k1 = slope(t[i], orbit[i,:])
            k2 = slope(t[i] + t_step, orbit[i, :] + t_step * k1)
    
            # compute the next state
            orbit[i + 1, :] = orbit[i, :] + t_step * (k1 + k2)/2
        
        
        # save the orbit
        print("Saving orbit...")

        data = np.column_stack((t, orbit))
        header = 't [s], x [m], y [m], vx [m/s], vy [m/s]'

        if save_dir != '':
            np.savetxt(save_dir, data, header = header, delimiter = ',')
        elif save_dir == '':
            np.savetxt(f'RK2_integrated_orbit.txt', data, header = header, delimiter = ',')

        return orbit
    

    
    def RK3(system, t_step, save_dir = ''):
        '''
        Method that integrates the system using the RK3 method, and saves the orbit in the provided directory. If no directory is provided, it will save the file in the current directory.  

        Arguments:
        _________

        system      instance of the TwoBodySimulation class containing the system to be integrated
        t_step      time step in days
        save_dir    directory to save the integrated orbit (default is './RK3_integrated_orbit.txt')

        Returns:
        _______

        orbit       array containing the integrated orbit (x, y, vx, vy) in SI units
        '''
        print("Integrating system...")

        # compute time array
        t       = time_array(system.period, t_step)
        t_step  = t_step * 86411 # approx. the length of a day in seconds

        # compute initial conditions
        u = system.ics()

        # initialize orbit array
        orbit = np.zeros((len(t), 4))

        # set initial conditions
        orbit[0, :] = u.T 
        
        # integrate using RK2
        for i in range(0, len(t) - 1):
            
            # define the slopes
            k1 = slope(t[i], orbit[i,:])
            k2 = slope(t[i] + t_step/2, orbit[i, :] + t_step * k1/2)
            k3 = slope(t[i] + t_step, orbit[i, :] - t_step * k1 + 2 * t_step * k2)
            
            # compute the next state
            orbit[i + 1, :] = orbit[i, :] + t_step * (k1 + 4 * k2 + k3)/6
        
        # save the orbit
        print("Saving orbit...")

        data = np.column_stack((t, orbit))
        header = 't [s], x [m], y [m], vx [m/s], vy [m/s]'

        if save_dir != '':
            np.savetxt(save_dir, data, header = header, delimiter = ',')
        elif save_dir == '':
            np.savetxt(f'RK3_integrated_orbit.txt', data, header = header, delimiter = ',')

        return orbit


    
    def RK4(system, t_step, save_dir = ''):
        '''        
        Method that integrates the system using the RK4 method, and saves the orbit in the provided directory. If no directory is provided, it will save the file in the current directory.

        Arguments:
        _________

        system      instance of the TwoBodySimulation class containing the system to be integrated
        t_step      time step in days
        save_dir    directory to save the integrated orbit (default is './RK4_integrated_orbit.txt')

        Returns:
        _______

        orbit       array containing the integrated orbit (x, y, vx, vy) in SI units
        '''
        print("Integrating system...")

        # compute time array
        t       = time_array(system.period, t_step)
        t_step  = t_step * 86411 # approx. the length of a day in seconds

        # compute initial conditions
        u = system.ics()

        # initialize orbit array
        orbit = np.zeros((len(t), 4))

        # set initial conditions
        orbit[0, :] = u.T 
        
        # integrate using RK2
        for i in range(0, len(t) - 1):
            
            # define the slopes
            k1 = slope(t[i], orbit[i,:])
            k2 = slope(t[i] + t_step/2, orbit[i, :] + t_step * k1/2)
            k3 = slope(t[i] + t_step/2, orbit[i, :] + t_step * k2/2)
            k4 = slope(t[i] + t_step, orbit[i, :] + t_step * k3)
            
            # compute the next state
            orbit[i + 1, :] = orbit[i, :] + t_step * (k1 + 2 * k2 + 2 * k3 + k4)/6
        
                # save the orbit
        print("Saving orbit...")

        data = np.column_stack((t, orbit))
        header = 't [s], x [m], y [m], vx [m/s], vy [m/s]'

        if save_dir != '':
            np.savetxt(save_dir, data, header = header, delimiter = ',')
        elif save_dir == '':
            np.savetxt(f'RK4_integrated_orbit.txt', data, header = header, delimiter = ',')

        return orbit



def plot_initial_system(system, save_dir = ''):
    '''
    Function that plots the initial system state.
    
    Arguments:
    _________

    system      instance of the TwoBodySimulation class containing the system to be integrated
    save_dir    directory to save the initial system plot (default is '')
    
    '''
    
    
    fig, ax = plt.subplots(figsize = (8, 8))

    # get initial state
    u0 = system.ics()
    
    # sun and earth images
    earth_img_url = 'https://raw.githubusercontent.com/GabrielBJ/practice-modules/main/earth.png'
    sun_img_url = 'https://raw.githubusercontent.com/GabrielBJ/practice-modules/main/sun.png'

    def fetch_image(path):
        response = requests.get(path)
        img = plt.imread(BytesIO(response.content))
        return img
    
    #earth_img = plt.imread('earth.png')
    #sun_img = plt.imread('sun.png')
    earth_img = fetch_image(earth_img_url)
    sun_img = fetch_image(sun_img_url)

    # include the images

    sun_marker = AnnotationBbox(OffsetImage(sun_img, zoom=0.07), \
            (0, 0), frameon=False, boxcoords = 'data', pad=0)
    ax.add_artist(sun_marker)

    earth_marker = AnnotationBbox(OffsetImage(earth_img, zoom=0.03), \
            (u0[0], u0[1]), frameon=False, boxcoords="data", pad=0)
    ax.add_artist(earth_marker)
    

    
    vel_box = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    fig.text(0.67, 0.82,f'$x$ = {u0[0]:.2e}, $x$ = {u0[1]:.2e}\n $v_x = ${u0[2]:.2e}, $v_y = $ {u0[3]:.2e}', fontsize=10,
            verticalalignment='bottom', horizontalalignment='right', bbox=vel_box)

    # plotting initial system 
    plt.plot(u0[0], u0[1], linestyle = '-', linewidth = 0.3)
    
    # decorate the plot
    ax.set_title(f"Initial system e ={system.e:.2f}")
    ax.set_xlabel(r'$x$ [m]')
    ax.set_ylabel(r'$y$ [m]')
    ax.grid(True)
    
    ax.set_xlim(-2e11, 2e11)
    ax.set_ylim(-2e11, 2.e11)
   
    if save_dir == '' :
       
        plt.savefig('./initial_system.png')

    else:
        plt.savefig(save_dir)

def animation_orbit(orbit, e, save_dir = ''):
    '''
    Function that animates the integrated orbit and saves the gif in the provided directory. If no directory is provided, it will save the file in the current directory.

    Arguments:
    _________

    orbit       array containing the integrated orbit (x, y, vx, vy) in SI units
    e           eccentricity of the orbit
    save_dir    directory to save the gif of the orbit (default is "./orbit.gif")

    '''
    
    x   = orbit[:,0]
    y   = orbit[:,1]
    vx  = orbit[:,2]
    vy  = orbit[:,3]

    fig, ax = plt.subplots(figsize = (8, 8))

    # sun and earth images
    earth_img_url = 'https://raw.githubusercontent.com/GabrielBJ/practice-modules/main/earth.png'
    sun_img_url = 'https://raw.githubusercontent.com/GabrielBJ/practice-modules/main/sun.png'

    def fetch_image(path):
        response = requests.get(path)
        img = plt.imread(BytesIO(response.content))
        return img
    
    #earth_img = plt.imread('earth.png')
    #sun_img = plt.imread('sun.png')
    earth_img = fetch_image(earth_img_url)
    sun_img = fetch_image(sun_img_url)


    earth_marker = AnnotationBbox(OffsetImage(earth_img, zoom=0.03), \
            (x[0], y[0]), frameon=False, boxcoords="data", pad=0)
    ax.add_artist(earth_marker)
    
    sun_marker = AnnotationBbox(OffsetImage(sun_img, zoom=0.07), \
            (0, 0), frameon=False, boxcoords = 'data', pad=0)
    ax.add_artist(sun_marker)

    line2, = ax.plot(x[0], y[0], linewidth = 0.4, c = 'b', label = 'orbit')

    vel_box = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    text_annotation = ax.text(0.33, 0.93, f'$x$ = {x[0]:.2e}, $y$ = {x[1]:.2e}\n$v_x$ = {vx[0]:.2e}, $v_y$ = {vy[0]:.2e}', fontsize=10, transform=ax.transAxes, bbox=vel_box)



    ax.set_title(f"Earth's orbit e = {e:.2f}")
    ax.set_xlabel(r'$x$ [m]')
    ax.set_ylabel(r'$y$ [m]')
    ax.grid(True)
    
    ax.set_xlim(-2e11, 2e11)
    
    if e < 0.3:
        ax.set_ylim(-2e11, 2.e11)
    elif e >= 0.3 and e < 0.5:
        ax.set_ylim(-3e11, 2.e11)
    elif e >= 0.5 and e < 0.7:
        ax.set_ylim(-3e11, 1.5e11)
    else:
        ax.set_ylim(-3e11, 1e11)


    def update(frame):
    
        line2.set_xdata(x[:frame])
        line2.set_ydata(y[:frame])
        
        new_position        = (x[frame], y[frame])
        earth_marker.xybox  = new_position
        earth_marker.xy     = new_position

        text_annotation.set_text(f'$x$ = {x[frame]:.2e}, $y$ = {y[frame]:.2e}\n$v_x$ = {vx[frame]:.2e}, $v_y$ = {vy[frame]:.2e}')


        return line2, earth_marker, text_annotation
    
    if save_dir == '':
        ani = FuncAnimation(fig = fig, func = update, frames = len(x), interval = 10, blit = True)
        ani.save('./orbit.gif', writer='pillow', fps=30)
        plt.show()
    elif save_dir != '':
        ani = FuncAnimation(fig = fig, func = update, frames = len(x), interval = 10, blit = True)
        ani.save(save_dir, writer='pillow', fps=30)
        plt.show()


def initialise_system(e, T, save_dir = ''):
    '''
    Function that initialises the system and plots it.

    Arguments:
    _________

    e           eccentricity of the orbit
    T           period of the orbit in years
    save_dir    directory to save the initial system plot (default is '')

    '''
    system = TwoBodySimulation(e = e, T = T)
    
    if save_dir != '':
        plot_initial_system(system, save_dir = save_dir)
    else:
        pass

    return system

def orbit_error(*orbits):
    '''
    Function that computes the radius of an orbit after a period T= 1 year.

    Arguments:
    _________

    *orbits      array containing the different orbits to be compared

    Returns:
    _______

    r            array containing the radius of the orbits after a period T = 1 year
    '''
    # compute the radius of the orbit
    r = [] 

    for orbit in orbits:
        x, y, vx, vy = np.split(orbit[-1:], 4, axis = 1)

        radius = np.sqrt(x**2 + y**2)
        r.append(radius)

    r = np.array(r)
    r = r.flatten()

    return r
    




def main():
    '''
    This method initialises the system, integrates it and plots the orbit when the script is run from the command line.

    Arguments:
    _________

    -e, --eccentricity      eccentricity of the orbit       float
    -T, --period            period of the orbit in years    float
    -m, --method            integration method              string
    -s, --savemap           directory to save the initial system plot (default is '') string
    -sg, --savegif          directory to save the gif of the orbit (default is '') string

    '''
    matplotlib.use('TkAgg')

    parser = argparse.ArgumentParser(description='Two-body simulation')
    parser.add_argument('-e', '--eccentricity', type=float, default=0.0167, help='eccentricity of the orbit')
    parser.add_argument('-T', '--period', type=float, default=1, help='period of the orbit in years')
    parser.add_argument('-m', '--method', type=str, default='RK2', help='integration method')
    parser.add_argument('-s', '--savemap', type=str, default='', help='directory to save the initial system')
    #parser.add_argument('-so', '--saveorbit', type=str, default='', help='directory to save the integrated orbit')
    parser.add_argument('-sg', '--savegif', type=str, default='', help='directory to save the gif of the orbit')

    args = parser.parse_args()

    print("Initialising system...")
    system = initialise_system(args.eccentricity, args.period, args.savemap)

    print("Integrating system...")
    if args.method == 'RK2':
        orbit = RKIntegrate.RK2(system, 1)
    elif args.method == 'RK3':
        orbit = RKIntegrate.RK3(system, 1)
    else:
        args.method == 'RK4'
        orbit = RKIntegrate.RK4(system, 1)

    # save the orbit
    print("Saving orbit...")
    t = time_array(args.period, 1)
    data = np.column_stack((t, orbit))

    header = 't [s], x [m], y [m], vx [m/s], vy [m/s]'
    np.savetxt(f'{args.method}_integrated_orbit.txt', data, header = header, delimiter = ',')
    
    print("Plotting orbit...")
    animation_orbit(orbit,args.eccentricity ,args.savegif)



if __name__ == '__main__':
    main()









