from Core import ProgramTimer, ScatteringObjectCollection, Simulation, Sensor

# Meta variables:
simulation_name = "Solution Name"
unit = "um"  # Unit which sets length scale of the simulation. Must be one of: 'm', 'dm', 'cm', 'mm', 'um', 'nm'.
signal_boost = 1  # Factor to multiply the output by to boost the signal.
using_GMRES = False  # Whether or not the simulation should use GMRES as the solver.
max_iters = None  # Maximum number of iterations to run the simulation if using GMRES.


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
# Remember that below, sensors must be a LIST of Sensor instances, so even if you just want one sensor, you must write sensors = [sensor] not sensors = sensor
sim.plot_sensor_intensity(sensors=[sensor1], signal_boost=signal_boost)

# End timer:
timer.stop()
