# odekepler.kepler module
This module generates a simulation of the two-body problem by solving the ODE system of equations of motions of the Earth:

$$\frac{d\vec{r}}{dt}=\vec{v}$$

$$m\frac{d\vec{v}}{dt}=-\frac{G\,m\,M}{r^3}\vec{r}$$

using the initial conditions:

$$x_0 = 0$$

$$y_0 = a\,(1-e)$$

$$v_{x0} = -\sqrt{\frac{G\,M}{a}\frac{1+e}{1-e}}$$

$$v_{y0} = 0$$

## General structure of the module 

```
midtermpart2.tar

    odekepler
    ├── odekepler
        ├── kepler.py
        ├── __init__.py
    ├── analysis.ipynb
    ├── outputfolder
    ├── helpfolder
    └── README.md
   
```

## Using the module as a script 

Download the `midtermpart2.tar` file. Uncompress it and copy the `odekepler` folder to your working directory.

Open in Terminal an run:
```
cd ~/path_to_module/odekepler/odekepler/
```
this will place you in the directory where the `kepler.py` file is placed.

Now, you are ready to use the `kepler.py` module as a script: 
```
python kepler.py -e  -T  -m  -s  -sg
```
here,
```
-e -> this flag takes as input the eccentricity of the orbit.
-T -> this flag takes as input the period in years.
-m -> this flag takes as input the method to perform the integration (RK2, RK3, RK4)
-s -> this flag takes as input the directory to save the initial map
-sg -> this flag takes as input the directory to save the gif animation of the orbit
```
### Example:
```
python kepler.py -e 0.25 -T 1 -m RK4 -s "./intial_map.png" -sg "./animation.gif"
```
This example will generate the simulation for the Pluto eccentricity using the RK4 method. Then, it will save the initial map in the `"./intial_map.png"` directory,  and the animation in the `"./animation.gif"` directory.

It will also generate a file called `RK4_integrated_orbit.txt`, which is saved in the current directory. 

## Importing the module
You can also import the module in a python notebook for example. 

### Installation
For this, download the `midtermpart2.tar` file. Open the file, and go to the `odekepler` folder, which contains the `setup.py` file. 
Then, in Terminal run:
```
python setup.py install --user
```
Now, you are ready to use the `odekepler` module 

### Example:
Import the module: 
```python
import odekpler.kepler as kp
```
initialising a system for e = 0.01671 (Earth's eccentricity) for a period of T=1 year
```python
system = kp.initialise_system(0.01671, 1, './eart_intial_system.png')
```
this will initialise a system with the initial conditions. Now, you can pass this system to one of the three integrators to get the Earth's orbit. : 

```python
orbit = earth_orbit = kp.RKIntegrate.RK2(earth_system, 1., "./earth_orbit.txt")
```
this will return the integrated orbit of the Earth. Now you can use this orbit to generate an animation using the following function:

```python
animation = kp.animation_orbit(orbit, 0.01671, "./earth_orbit.gif")
```
This will generate a gif file containing the simulation.

The module also include a function to get the time array, given the period T and the time step

```python
time = kp.time_array(T, time_step)
```
it will be useful if you want to generate plots vs. time. 

If you want to analyse the convergence of the RK methods, you can use the following function:

```python
rk3_r = kp.orbit_error(*rk3_orbits)
```
this function receives multiple orbits integrated with the same method `RK2, RK3, RK4` with period 1, but using different time step. It will return the radius of the orbit after a T=1, so you can see how the precision of the methods varies as the time step changes. 

## Usage Examples 
You can find examples and the documentation of the methods used here in the `helpfolder`

### Good luck!
