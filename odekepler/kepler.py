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
import scienceplots



# create a class that initialises the system
class TwoBodySimulation():
    '''
    This class initialises the system to be integrated. It takes the eccentricity and period of the orbit as arguments.

    Attributes:
            -e (float) : eccentricity of the orbit
            -T (float) : period of the orbit

    Methods:
            -ics() (function) : returns the initial conditions of the system

    '''
    def __init__(self, e, T):
        '''
        Initialises the system with the eccentricity and period of the orbit.
        '''
        self.e       = e
        self.period  = T
    

    def ics(self):
        '''
        Computes the initial conditions of the system given the eccentricity and period of the orbit.

        Returns:
                -u (np.array) : vector containing the initial conditions of the system (x0, y0, vx0, vy0) in astronomic units
        '''
        # initialise the initial conditions vector
        u = np.zeros(4)
        
        # compute the initial conditions
        u[0] = 0 
        u[1] = 1.0 * (1 - self.e)
        u[2] = -1 * np.sqrt(4 * np.pi**2 * (1 + self.e)/(1 - self.e))
        u[3] = 0
        
        return u 


def slope(t, state):
    '''
    Function that computes the slope of the system given the state of the system at a certain time.

    Arguments:
            -state (np.array) : vector containing the state of the system (x, y, vx, vy) in astronomic units.

    Returns:
            -du (np.array) : vector containing the slope of the system (vx, vy, ax, ay) in astronomic units per year
    '''
    # compute the distance from the sun
    r = np.sqrt(state[0]**2 + state[1]**2)

    # compute G * M_sun
    k = 4 * np.pi**2
    
    # initialise the slope vector
    du = np.zeros(4)
    
    # compute the slope
    for i in range(2):
        du[i]       = state[i + 2]
        du[i + 2]   = -k * state[i] * 1/r**3

    return du

def time_array(T , dt):
    '''
    Function that computes the time array given the period of the orbit and the time step.

    Arguments:
            -T  (float) : period of the orbit in years
            -dt (float) : time step in years (values between 0 and 1)

    Returns:
            -t (np.array) : time array in years
    '''
    period  = T  # in years
    dt      = dt # time step in years

    # compute the time array
    t = np.arange(0, period + dt, dt)

    return t


class RKIntegrate:
    '''
    This class contains the methods to integrate the system using the Runge-Kutta methods.

    Methods:
            -RK2(system, t_step) (functions) : integrates the system using the RK2 method
            -RK3(system, t_step) (functions) : integrates the system using the RK3 method
            -RK4(system, t_step) (functions) : integrates the system using the RK4 method
    '''

    def RK2(system, t_step, save_dir = ''):
        '''
        Method that integrates the system using the RK2 method and saves the orbit in 
        the provided directory. If no directory is provided, it will save the file in 
        the current directory (./RK2_integrated_orbit.txt).

        Arguments:
                -system   (class instance) : instance of the TwoBodySimulation class containing the system to be integrated
                -t_step   (float)          : time step in years (values between 0 and 1)
                -save_dir (str)            : directory to save the integrated orbit (default is './RK2_integrated_orbit.txt')

        Returns:
                -orbit (np.array) : array containing the integrated orbit (x, y, vx, vy) in SI units
        '''
        print("Integrating system...") 

        # compute time array                                                                    
        t       = time_array(system.period, t_step)
        t_step  = t_step 

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
        print("Saving integrated orbit...")

        data   = np.column_stack((t, orbit))
        header = 't [yrs], x [au], y [au], vx [au/yrs], vy [au/yrs]'

        if save_dir != '':
            np.savetxt(save_dir, data, header = header, delimiter = ',')
        elif save_dir == '':
            np.savetxt(f'RK2_integrated_orbit.txt', data, header = header, delimiter = ',')

        return orbit
    

    
    def RK3(system, t_step, save_dir = ''):
        '''
        Method that integrates the system using the RK2 method and saves the orbit in 
        the provided directory. If no directory is provided, it will save the file in 
        the current directory (./RK3_integrated_orbit.txt).

        Arguments:
                -system   (class instance) : instance of the TwoBodySimulation class containing the system to be integrated
                -t_step   (float)          : time step in years (values between 0 and 1)
                -save_dir (str)            : directory to save the integrated orbit (default is './RK3_integrated_orbit.txt')

        Returns:
                -orbit (np.array) : array containing the integrated orbit (x, y, vx, vy) in SI units
        '''
        print("Integrating system...")

        # compute time array
        t       = time_array(system.period, t_step)
        t_step  = t_step 

        # compute initial conditions
        u = system.ics()

        # initialize orbit array
        orbit = np.zeros((len(t), 4))

        # set initial conditions
        orbit[0, :] = u.T 
        
        # integrate using RK3
        for i in range(0, len(t) - 1):
            
            # define the slopes
            k1 = slope(t[i], orbit[i,:])
            k2 = slope(t[i] + t_step/2, orbit[i, :] + t_step * k1/2)
            k3 = slope(t[i] + t_step, orbit[i, :] - t_step * k1 + 2 * t_step * k2)
            
            # compute the next state
            orbit[i + 1, :] = orbit[i, :] + t_step * (k1 + 4 * k2 + k3)/6
        
        # save the orbit
        print("Saving integrated orbit...")

        data   = np.column_stack((t, orbit))
        header = 't [yrs], x [au], y [au], vx [au/yrs], vy [au/yrs]'

        if save_dir != '':
            np.savetxt(save_dir, data, header = header, delimiter = ',')

        elif save_dir == '':
            np.savetxt(f'RK3_integrated_orbit.txt', data, header = header, delimiter = ',')

        return orbit


    
    def RK4(system, t_step, save_dir = ''):
        '''
        Method that integrates the system using the RK2 method and saves the orbit in 
        the provided directory. If no directory is provided, it will save the file in 
        the current directory (./RK4_integrated_orbit.txt).

        Arguments:
                -system   (class instance) : instance of the TwoBodySimulation class containing the system to be integrated
                -t_step   (float)          : time step in years (values between 0 and 1)
                -save_dir (str)            : directory to save the integrated orbit (default is './RK4_integrated_orbit.txt')

        Returns:
                -orbit (np.array) : array containing the integrated orbit (x, y, vx, vy) in SI units
        '''
        print("Integrating system...")

        # compute time array
        t       = time_array(system.period, t_step)
        t_step  = t_step 

        # compute initial conditions
        u = system.ics()

        # initialize orbit array
        orbit = np.zeros((len(t), 4))

        # set initial conditions
        orbit[0, :] = u.T 
        
        # integrate using RK4
        for i in range(0, len(t) - 1):
            
            # define the slopes
            k1 = slope(t[i], orbit[i,:])
            k2 = slope(t[i] + t_step/2, orbit[i, :] + t_step * k1/2)
            k3 = slope(t[i] + t_step/2, orbit[i, :] + t_step * k2/2)
            k4 = slope(t[i] + t_step, orbit[i, :] + t_step * k3)
            
            # compute the next state
            orbit[i + 1, :] = orbit[i, :] + t_step * (k1 + 2 * k2 + 2 * k3 + k4)/6
        
                # save the orbit
        print("Saving integrated orbit...")

        data   = np.column_stack((t, orbit))
        header = 't [yrs], x [au], y [au], vx [au/yrs], vy [au/yrs]'

        if save_dir != '':
            np.savetxt(save_dir, data, header = header, delimiter = ',')

        elif save_dir == '':
            np.savetxt(f'RK4_integrated_orbit.txt', data, header = header, delimiter = ',')

        return orbit


def plot_initial_system(system, save_dir = ''):
    '''
    Function that plots the initial system state.
    
    Arguments:
            -system   (class instance) : instance of the TwoBodySimulation class containing the system to be integrated
            -save_dir (str)            : directory to save the initial system plot (default is '')
    
    '''
    with plt.style.context(['notebook', 'no-latex', 'grid']):

        fig, ax = plt.subplots(figsize = (8, 8))

        # get initial state
        u0 = system.ics()
        
        # sun and earth images
        earth_img_url = 'https://raw.githubusercontent.com/GabrielBJ/practice-modules/main/earth.png'
        sun_img_url   = 'https://raw.githubusercontent.com/GabrielBJ/practice-modules/main/sun.png'

        def fetch_image(path):
            response = requests.get(path)
            img      = plt.imread(BytesIO(response.content))

            return img

        earth_img = fetch_image(earth_img_url)
        sun_img   = fetch_image(sun_img_url)

        # include sun and earth images
        sun_marker = AnnotationBbox(OffsetImage(sun_img, zoom=0.07), \
                (0, 0), frameon=False, boxcoords = 'data', pad=0)
        ax.add_artist(sun_marker)

        earth_marker = AnnotationBbox(OffsetImage(earth_img, zoom=0.03), \
                (u0[0], u0[1]), frameon=False, boxcoords="data", pad=0)
        ax.add_artist(earth_marker)
        
        # Add text box with initial conditions
        vel_box = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        fig.text(0.67, 0.82,f'$x$ = {u0[0]:.2e}, $y$ = {u0[1]:.2e}\n$v_x = ${u0[2]:.2e}, $v_y = $ {u0[3]:.2e}', fontsize=10,
                verticalalignment='bottom', horizontalalignment='right', bbox=vel_box)

        # plotting initial system 
        plt.plot(u0[0], u0[1], linestyle = '-', linewidth = 0.3, color = "black")
        
        # decorate the plot
        ax.set_title(f"Init system: e ={system.e:.2f}")
        ax.set_xlabel(r'$x$ [au]')
        ax.set_ylabel(r'$y$ [au]')
        ax.set_xlim(-1.5, 1.5)
        ax.set_ylim(-1.5, 1.5)

        if save_dir != '' :
            plt.savefig(save_dir)
        else:
            plt.savefig('./initial_system.png')
        plt.show()

def animation_orbit(orbit, e, save_dir = ''):
    '''
    Function that animates the integrated orbit and saves the gif 
    in the provided directory. If no directory is provided, it will 
    just display the animation without saving it.

    Arguments:
            -orbit    (np.array) : array containing the integrated orbit (x, y, vx, vy) in SI units
            -e        (float)    : eccentricity of the orbit
            -save_dir (str)      : directory to save the gif of the orbit.

    '''
    x   = orbit[:,0]
    y   = orbit[:,1]
    vx  = orbit[:,2]
    vy  = orbit[:,3]

    with plt.style.context(['notebook', 'no-latex', 'grid']):
        fig, ax = plt.subplots(figsize = (8, 8))

        # sun and earth images
        earth_img_url = 'https://raw.githubusercontent.com/GabrielBJ/practice-modules/main/earth.png'
        sun_img_url   = 'https://raw.githubusercontent.com/GabrielBJ/practice-modules/main/sun.png'

        def fetch_image(path):
            response = requests.get(path)
            img      = plt.imread(BytesIO(response.content))

            return img
        
        earth_img = fetch_image(earth_img_url)
        sun_img   = fetch_image(sun_img_url)

        earth_marker = AnnotationBbox(OffsetImage(earth_img, zoom=0.03), \
                (x[0], y[0]), frameon=False, boxcoords="data", pad=0)
        ax.add_artist(earth_marker)
    
        sun_marker = AnnotationBbox(OffsetImage(sun_img, zoom=0.07), \
                (0, 0), frameon=False, boxcoords = 'data', pad=0)
        ax.add_artist(sun_marker)

        line2, = ax.plot(x[0], y[0], linewidth = 0.4, c = 'black', label = 'orbit')

        vel_box = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        text_annotation = ax.text(0.33, 0.93, f'$x$ = {x[0]:.2e}, $y$ = {x[1]:.2e}\n$v_x$ = {vx[0]:.2e}, $v_y$ = {vy[0]:.2e}', fontsize=10, transform=ax.transAxes, bbox=vel_box)

        ax.set_title(f"Orbit: e = {e:.2f}")
        ax.set_xlabel(r'$x$ [au]')
        ax.set_ylabel(r'$y$ [au]')
        ax.set_xlim(-1.5, 1.5)
        
        if e < 0.3:
            ax.set_ylim(-1.5, 1.5)
        elif e >= 0.3 and e < 0.5:
            ax.set_ylim(-2, 1)
        elif e >= 0.5 and e < 0.7:
            ax.set_ylim(-2, 1)
        else:
            ax.set_ylim(-2, 1)

        def update(frame):
        
            line2.set_xdata(x[:frame])
            line2.set_ydata(y[:frame])
            
            new_position        = (x[frame], y[frame])
            earth_marker.xybox  = new_position
            earth_marker.xy     = new_position

            text_annotation.set_text(f'$x$ = {x[frame]:.2e}, $y$ = {y[frame]:.2e}\n$v_x$ = {vx[frame]:.2e}, $v_y$ = {vy[frame]:.2e}')

            return line2, earth_marker, text_annotation
        
        if save_dir == '':
            ani = FuncAnimation(fig = fig, func = update, frames = len(x), interval = 40, blit = True)
            plt.show()
            #ani.save('./orbit.gif', writer='pillow', fps=30)
        elif save_dir != '':
            ani = FuncAnimation(fig = fig, func = update, frames = len(x), interval = 40, blit = True)
            plt.show()
            ani.save(save_dir, writer='pillow', fps=30)

def initialise_system(e, T, save_dir = ''):
    '''
    Function that initialises the system and plots it.

    Arguments:
            -e        (float) : eccentricity of the orbit
            -T        (float) : period of the orbit in years
            -save_dir (srt)   : directory to save the initial system plot (default is '')

    '''
    system = TwoBodySimulation(e = e, T = T)
    plot_initial_system(system, save_dir = save_dir)
    return system

def orbit_error(*orbits):
    '''
    Function that computes the radius of an orbit after a period T= 1 year.

    Arguments:
            -orbits (np.array) : array containing the different orbits to be compared

    Returns:
            -r (np.array) : array containing the radius of the orbits after a period T = 1 year
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
            -e  (float) --eccentricity  : eccentricity of the orbit      
            -T  (float) --period        : period of the orbit in years   
            -m  (float) --method        : integration method              
            -s  (float) --savemap       : directory to save the initial system plot (default is '') 
            -sg (float) --savegif       : directory to save the gif of the orbit (default is '') 

    '''
    matplotlib.use('TkAgg') # to avoid the error: TclError: no display name and no $DISPLAY environment variable

    parser = argparse.ArgumentParser(description = 'Two-body simulation')
    parser.add_argument('-e', '--eccentricity', type = float, default = 0.0167, help = 'eccentricity of the orbit')
    parser.add_argument('-T', '--period', type = float, default = 1, help = 'period of the orbit in years')
    parser.add_argument('-dt', '--timestep', type = float, default = 0.01, help = 'time step in yrs')
    parser.add_argument('-m', '--method', type = str, default = 'RK2', help = 'integration method')
    parser.add_argument('-s', '--savemap', type = str, default = '', help = 'directory to save the initial system')
    #parser.add_argument('-so', '--saveorbit', type=str, default='', help='directory to save the integrated orbit')
    parser.add_argument('-sg', '--savegif', type = str, default = '', help = 'directory to save the gif of the orbit')

    args = parser.parse_args()

    print("Initialising system...")
    system = initialise_system(args.eccentricity, args.period, args.savemap)

    #print("Integrating system...")
    if args.method == 'RK2':
        orbit = RKIntegrate.RK2(system, args.timestep)
    elif args.method == 'RK3':
        orbit = RKIntegrate.RK3(system, args.timestep)
    else:
        args.method == 'RK4'
        orbit = RKIntegrate.RK4(system, args.timestep)

    # save the orbit
    t    = time_array(args.period, args.timestep)
    data = np.column_stack((t, orbit))

    header = 't [yrs], x [au], y [au], vx [au/yrs], vy [au/yrs]'
    np.savetxt(f'{args.method}_integrated_orbit.txt', data, header = header, delimiter = ',')
    
    print("Plotting orbit...")
    animation_orbit(orbit,args.eccentricity ,args.savegif)

if __name__ == '__main__':
    main()