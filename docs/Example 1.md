# Example 1: Scattering off of Three Spheres
In this first tutorial, we will walk through the process of using UFBS to run a scattering simulation featuring three spheres.

## Step 1: Main.py
The set up for a UFBS simulation begins in the Main.py file. Our first order of business will be to specify a name and unit of length for our simulation. All other lengths and constants in the simulation are measured with respect to this chosen unit of length. For example, later on we will have to provide a number for the wavelength of the incident light for this simulation, and that number will be with respect to this unit of length. As another example, we will also have to provide the positions of the objects we wish to scatter off of, and those positions will be in terms of this unit of length. For this tutorial, let's name our simulation "Three Spheres Tutorial" and let's set our unit of length to "um" (micrometers) so the top of the Main.py file looks like:

```python
from Core import ProgramTimer, ScatteringObjectCollection, Simulation, Sensor

# Meta variables:
simulation_name = "Three Spheres Tutorial"
unit = "um"  # Unit which sets length scale of the simulation. Must be one of: 'm', 'dm', 'cm', 'mm', 'um', 'nm'.
```

Next you will see that there are three variables called signal_boost, using_GMRES and max_iters. The signal_boost variable allows us to multiply the signal output of the simulation by a constant to make it larger in value. The using_GMRES variable tells the program if we would like to use the GMRES iterative solver for this simulation. The max_iters variable tells the program how many iterations to use if we are using GMRES. For now, let's not boost the signal and not use GMRES by setting these variables as follows:

```python
from Core import ProgramTimer, ScatteringObjectCollection, Simulation, Sensor

# Meta variables:
simulation_name = "Three Spheres Tutorial"
unit = "um"  # Unit which sets length scale of the simulation. Must be one of: 'm', 'dm', 'cm', 'mm', 'um', 'nm'.
signal_boost = 1  # Factor to multiply the output by to boost the signal.
using_GMRES = False
max_iters = None  # Maximum number of iterations to run the simulation for.
```
Next, let's skip a few lines and travel down to where it says sim.plot_scattering_objs(). This line instructs the program to display a preview of our scattering objects (which we will define later) using gmsh. If you chose not to download gmsh during the installation of this UFBS, please comment-out, or delete this line. 

Now let's skip down to the line where it says sensor1 = Sensor(...). This is where we define the sensors on which we would like to detect the scattered light. In UFBS, a sensor is a flat plane at some position which acts like a camera sensor and detects the intensity of light hitting it. The final output of the simulation will be a plot showing the intensity of scattered light on this sensor. For simplicity, let's set our sensor to be 1cm (10000um) away from the origin, and let's make its width and height be 1cm and 1cm. The width_samples and height_samples variables tell the program what the resolution of our sensor should be; let's set them both to 300 for now. For an example of how to use the theta and gamma parameters, please see Example 3, but for now, we will set them to zero. After all our set up, the Main.py file should look like this:

```python
from Core import ProgramTimer, ScatteringObjectCollection, Simulation, Sensor

# Meta variables:
simulation_name = "Three Spheres Tutorial"
unit = "um"  # Unit which sets length scale of the simulation. Must be one of: 'm', 'dm', 'cm', 'mm', 'um', 'nm'.
signal_boost = 1  # Factor to multiply the output by to boost the signal.
using_GMRES = False  # Whether or not the simulation should use GMRES as the solver.
max_iters = None  # Maximum number of iterations to run the simulation for if using GMRES.


# Beginning a timer to see how long the simulation takes:
timer = ProgramTimer()
timer.start()

# Loading scattering objects from the Scattering Objects folder:
scattering_objs = ScatteringObjectCollection()

# Running simulation:
sim = Simulation(simulation_name=simulation_name, unit=unit, scattering_objs=scattering_objs, GMRES=using_GMRES, max_iters=max_iters)
sim.plot_scattering_objs()  # Please comment out or delete this line if you have not downloaded Gmsh.
sim.solve()

# Plotting scattered field intensity on sensor:
# Each of the length values here are in units of the unit variable.
sensor1 = Sensor(distance=10000, theta=0, gamma=0, width=10000, height=10000, width_samples=300, height_samples=300)
# Remember that below, sensors must be a list of Sensor types, so even if you just want one sensor, you must write sensors = [sensor] not sensors = sensor
sim.plot_sensor_intensity(sensors=[sensor1], signal_boost=signal_boost)

# End timer:
timer.stop()
```

## Step 2: Sphere Informtion
Next, we must tell UFBS what objects we want to scatter off of, and what properties they have. To start out, we will use one of the "Standard Objects" which are already built into the UFBS. There are four standard objects: cubes, cuboids, screens, and spheres. For this tutorial we will instruct UFBS to use three spheres. To do this, we can open the Spheres.txt file by following this path: Scattering Objects $\rightarrow$ Standard Objects $\rightarrow$ Spheres.txt. In this file you will see that the very first line gives instructions on what information you should provide to specify a sphere:

\# x y z (relative permittivity) (relative permeability) radius (max element size) \#

According to these instructions, to specify a sphere we must provide the x, y, and z coordinates of its center in um, its realative permitivity, its relative permeability, its radius in um, and the maximum size of one of it's boundary mesh elements in um. To specify multiple spheres, we must write multiple lines into Spheres.txt, each providing the above information about a single sphere. Since we want to include three spheres, this is what our Spheres.txt will look like for this simulation:

\# x y z (relative permittivity) (relative permeability) radius (max element size) \#
<br>1 0 0.2 3.0 1.2 0.6 0.2 
<br>-1 0 -0.2 1.5.0 1.1 0.5 0.2
<br>0 3 0.5 2.0 1.0 0.7 0.2

Notice that each number is separated by a space. This is essential for UFBS to interpret each line correctly. Multiple spaces is okay too. To ensure that no other standard objects are included in the simulation, the user should look through the other text files for the standard objects and ensure that they only include the first line which gives instructions on what information to provide. Thus each of Cubes.txt, Cuboids.txt, and Screens.txt should only be one line long for this tutorial.

## Step 3: Incident Wave Informtion
The next step is to tell UFBS what the incident light wave hitting the scattering objects looks like. This can be done in the Incident_Wave_Params.py file. For this tutorial, let's use a wave with wavelength 0.5um, directed along the z axis, with polarization along the x axis. To specify this information, the Incident_Wave_Params.py file should look like this:

```python
"""This file contains parameters for the incident wave used in scattering problems.
Unfortunately we cannot just define these in Main.py because of how the underlying Bempp library works"""

import numpy as np

INC_WL = 0.5  # Wavelength of the incident plane wave in units of UNIT as defined in Main.py
INC_DIR = np.array([0, 0, 1.0])  # Direction of the incident plane wave. Must be a unit vector. Has no units.
INC_POL = np.array([1.0, 0, 0])  # Polarization of the incident plane wave. Must be a unit vector. Has no units.
```

## Step 4: Running the Simulation
Now that we have provided all the information necessary to UFBS, we are ready to run the simulation. To achive this, simply run the Main.py file. 

Not done yet...