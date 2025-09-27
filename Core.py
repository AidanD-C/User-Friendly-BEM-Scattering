import bempp_cl.api
from bempp_cl.api.operators.boundary.maxwell import multitrace_operator
from bempp_cl.api import ZeroBoundaryOperator, BlockedOperator
import numpy as np
from numpy.typing import NDArray
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.spatial.transform import Rotation
import time
import logging
import os
from pathlib import Path
from Incident_Wave_Params import INC_WL, INC_DIR, INC_POL

K = 2 * np.pi / INC_WL  # wave number in medium

C = 299792458.0  # m/s
MU_0 = 1.256637061e-6  # SI units (MKSA)

root_directory = Path(__file__).parent.resolve()
scat_objs_directory = Path(__file__).parent.resolve() / "Scattering Objects"
if os.name == "nt":  # Windows
    gmsh_name = "gmsh.exe"
else:  # Linux, macOS
    gmsh_name = "gmsh"
gmsh_path = list(root_directory.rglob(gmsh_name))[0]
bempp_cl.api.GMSH_PATH = str(gmsh_path)

logging.basicConfig(level=logging.INFO)
bempp_cl.api.enable_console_logging()


class Sensor:
    """
    A float observation plane on which the field intensity can be measured and plotted.
    The sensor is a rectangular plane at a distance from the origin, with a specified width and height.
    The sensor can be placed at any distance and can be rotated relative to the origin to achieve
    different viewing angles.

    Attributes:
    - distance: Distance from the origin to the middle point of the sensor plane.
    - theta: Rotation angle about the x axis in degrees.
    - gamma: Rotation angle about the y axis in degrees.
    - width: Width of the sensor plane.
    - height: Height of the sensor plane.
    - width_samples: Number of samples along the width of the sensor planes.
    - height_samples: Number of samples along the height of the sensor plane.
    - points: The 2D array of points in 3D space which are located on the sensor plane.
    """

    distance: float
    theta: float
    gamma: float
    width: float
    height: float
    width_samples: int
    height_samples: int
    points: NDArray[np.float64]
    unit_normal: NDArray[np.float64]

    def __init__(self, distance: float, theta: float, gamma: float, width: float, height: float, width_samples: int, height_samples: int):

        self.distance = distance
        self.theta = theta
        self.gamma = gamma
        self.width = width
        self.height = height
        self.width_samples = width_samples
        self.height_samples = height_samples

        num_points = self.width_samples * self.height_samples

        width_array = np.linspace(start=-width / 2, stop=width / 2, num=width_samples)
        height_array = np.linspace(start=-height / 2, stop=height / 2, num=height_samples)

        # making vertices
        self.points = np.zeros(shape=(num_points, 3))

        point_counter = 0
        for x in width_array:
            for y in height_array:

                self.points[point_counter] = np.array([x, y, self.distance])
                point_counter += 1

        x_rot_vec = np.array([theta, 0, 0])
        y_rot_vec = np.array([0, gamma, 0])

        x_rot_matrix = Rotation.from_rotvec(x_rot_vec, degrees=True)
        y_rot_matrix = Rotation.from_rotvec(y_rot_vec, degrees=True)

        full_rot_matrix = y_rot_matrix * x_rot_matrix

        self.points = full_rot_matrix.apply(self.points)

        self.unit_normal = full_rot_matrix.apply(np.array([0, 0, 1]))


class ScatteringObject:
    """
    This the base class for a general object which will be scattered off of.

    Attributes:
    - position: A translation vector which can move the object to a desired position in space.
    - eps_r: Relative permittivity of the object.
    - mu_r: Relative permeability of the object.
    - file_path: Path to a file which contains the object geometry. If None, then the object must also be an instance of a standard shape class (Cube, Sphere, etc.)
    - k_int: Wave number inside the object.
    """

    position: NDArray[np.float64]
    eps_r: float
    mu_r: float
    file_path: Path | None
    k_int: float

    def __init__(self, position: NDArray[np.float64], eps_r: float, mu_r: float, file_path: Path | None = None) -> None:
        """
        Initializing the scattering object attributes.
        """
        self.position = position
        self.eps_r = eps_r
        self.mu_r = mu_r
        self.file_path = file_path
        self.k_int = K * np.sqrt(self.eps_r * self.mu_r)

    def generate_grid(self) -> bempp_cl.api.grid.Grid:
        """
        Based on the geometry info contained in self.file_path, this will generate and return a bempp grid for the object.
        This method will be overridden in the child classes for standard shapes however.
        These grids are required for the bempp simulation.
        """
        return bempp_cl.api.import_grid(str(self.file_path))


class Cube(ScatteringObject):
    """
    A cube scattering object.

    Attributes:
    - position: A translation vector which can move the cube to a desired position in space.
    - eps_r: Relative permittivity of the cube.
    - mu_r: Relative permeability of the cube.
    - side_length: Length of each side of the cube.
    - ele_size: The maximum element size for the grid which will be generated for the cube.
    - k_int: Wave number inside the cube.
    - file_path: See ScatteringObject class for description. For this class, it will always be None.
    """

    position: NDArray[np.float64]
    eps_r: float
    mu_r: float
    side_length: float
    ele_size: float
    k_int: float

    def __init__(self, position: NDArray[np.float64], eps_r: float, mu_r: float, side_length: float, element_size: float) -> None:
        """
        Initializing the cube attributes.
        """
        super().__init__(position, eps_r, mu_r)
        self.side_length = side_length
        self.ele_size = element_size

    def generate_grid(self) -> bempp_cl.api.grid.Grid:
        """
        Returns a bempp grid for this cube with the specified attributes.
        """

        return bempp_cl.api.shapes.cube(length=self.side_length, origin=self.position, h=self.ele_size)


class Cuboid(ScatteringObject):
    """
    A cuboid scattering object.

    Attributes:
    - position: A translation vector which can move the cube to a desired position in space.
    - eps_r: Relative permittivity of the cube.
    - mu_r: Relative permeability of the cube.
    - side_lengths: Lengths of each side of the cube.
    - ele_size: The maximum element size for the grid which will be generated for the cube.
    - k_int: Wave number inside the cube.
    - file_path: See ScatteringObject class for description. For this class, it will always be None.
    """

    position: NDArray[np.float64]
    eps_r: float
    mu_r: float
    side_lengths: tuple[float, float, float]
    ele_size: float
    k_int: float

    def __init__(self, position: NDArray[np.float64], eps_r: float, mu_r: float, side_lengths: tuple[float, float, float], element_size: float) -> None:
        """ """
        super().__init__(position, eps_r, mu_r)
        self.side_lengths = side_lengths
        self.ele_size = element_size

    def generate_grid(self) -> bempp_cl.api.grid.Grid:
        """
        Returns a bempp grid for this cuboid with the specified attributes.
        """

        return bempp_cl.api.shapes.cuboid(length=(self.side_lengths[0], self.side_lengths[1], self.side_lengths[2]), origin=self.position, h=self.ele_size)


class Screen(ScatteringObject):
    """
    A screen scattering object.

    Attributes:
    - position: A translation vector which can move the screen to a desired position in space. This will always be zero for this class. Only corners are used to specify position.
    - corners: A 2D array of 3D corner positions which define the screen vertices.
    - eps_r: Relative permittivity of the screen.
    - mu_r: Relative permeability of the screen.
    - ele_size: The maximum element size for the grid which will be generated for the screen
    - k_int: Wave number on the screen.
    - file_path: See ScatteringObject class for description. For this class, it will always be None.
    """

    position: NDArray[np.float64] = np.array([0, 0, 0])
    eps_r: float
    mu_r: float
    corners: NDArray[np.float64]
    ele_size: float
    k_int: float

    def __init__(self, corners: NDArray[np.float64], eps_r: float, mu_r: float, element_size: float) -> None:
        """
        Initializing the screen attributes.
        """
        super().__init__(np.array([0, 0, 0]), eps_r, mu_r)
        self.corners = corners
        self.ele_size = element_size

    def generate_grid(self) -> bempp_cl.api.grid.Grid:
        """
        Returns a bempp grid for this screen with the specified attributes.
        """

        return bempp_cl.api.shapes.screen(corners=self.corners, h=self.ele_size)


class Sphere(ScatteringObject):
    """
    A shpere scattering object.

    Attributes:
    - position: A translation vector which can move the sphere to a desired position in space.
    - eps_r: Relative permittivity of the sphere.
    - mu_r: Relative permeability of the sphere.
    - radius: Radius of the sphere.
    - ele_size: The maximum element size for the grid which will be generated for the sphere
    - k_int: Wave number inside the sphere.
    - file_path: See ScatteringObject class for description. For this class, it will always be None.
    """

    position: NDArray[np.float64]
    eps_r: float
    mu_r: float
    radius: float
    ele_size: float
    k_int: float

    def __init__(self, position: NDArray[np.float64], eps_r: float, mu_r: float, radius: float, element_size: float) -> None:
        """
        Initializing the sphere attributes.
        """
        super().__init__(position, eps_r, mu_r)
        self.radius = radius
        self.ele_size = element_size

    def generate_grid(self) -> bempp_cl.api.grid.Grid:
        """
        Returns a bempp grid for this sphere with the specified attributes.
        """

        return bempp_cl.api.shapes.sphere(r=self.radius, origin=self.position, h=self.ele_size)


class ScatteringObjectCollection:
    """
    A collection of scattering objects (really just a wrapper for a list of ScatteringObject instances).
    This class will handle the loading of scattering objects from files and the extraction of their attributes.
    It will also handle the overall position of the objects in space, which can be used to translate all objects at once.

    Attributes:
    - file_path: the path to the folder which contains the scattering objects files (the Scattering Objects folder).
    - overall_position: a translation vector which can move all objects to a desired position in space.
    - scat_objs: a list of ScatteringObjects which are extracted from the files in the file_path.
    """

    file_path: Path
    overall_position: NDArray[np.float64]
    scat_objs: list[ScatteringObject]

    def __init__(self, file_path: Path = scat_objs_directory, overall_position: NDArray[np.float64] = np.array([0, 0, 0])) -> None:
        """
        Initializing the scattering object collection attributes and extracting the scattering objects from the files in the file_path.
        """

        print("Initializing Scattering Objects...")

        self.file_path = file_path
        self.overall_position = overall_position
        self.scat_objs = []
        self._extract_cubes()
        self._extract_cuboids()
        self._extract_screens()
        self._extract_spheres()
        self._extract_user_objs()

    def _extract_cubes(self) -> None:
        """
        Extends self.scat_objs with a list of Cubes based on the information in the Cubes.txt file.
        """

        cubes_path = self.file_path / "Standard Objects" / "Cubes.txt"

        cubes = []

        with open(cubes_path, "r") as f:
            lines = f.read().splitlines()[1:]
            for line in lines:
                parts = line.split()
                position = np.array([float(parts[0]), float(parts[1]), float(parts[2])])
                eps_r = float(parts[3])
                mu_r = float(parts[4])
                side_length = float(parts[5])
                ele_size = float(parts[6])
                cube = Cube(position + self.overall_position, eps_r, mu_r, side_length, ele_size)
                cubes.append(cube)

        self.scat_objs.extend(cubes)

    def _extract_cuboids(self) -> None:
        """
        Extends self.scat_objs with a list of Cuboids based on the information in the Cuboids.txt file.
        """

        cuboids_path = self.file_path / "Standard Objects" / "Cuboids.txt"

        cuboids = []

        with open(cuboids_path, "r") as f:
            lines = f.read().splitlines()[1:]
            for line in lines:
                parts = line.split()
                position = np.array([float(parts[0]), float(parts[1]), float(parts[2])])
                eps_r = float(parts[3])
                mu_r = float(parts[4])
                length = float(parts[5])
                width = float(parts[6])
                height = float(parts[7])
                ele_size = float(parts[8])
                cuboid = Cuboid(position + self.overall_position, eps_r, mu_r, (length, width, height), ele_size)
                cuboids.append(cuboid)

        self.scat_objs.extend(cuboids)

    def _extract_screens(self) -> None:
        """
        Extends self.scat_objs with a list of Screens based on the information in the Screens.txt file.
        """

        screens_path = self.file_path / "Standard Objects" / "Screens.txt"

        screens = []

        with open(screens_path, "r") as f:
            lines = f.read().splitlines()[1:]
            for line in lines:
                parts = line.split()
                corner_positions = np.array([[float(parts[0]), float(parts[1]), float(parts[2])], [float(parts[3]), float(parts[4]), float(parts[5])], [float(parts[6]), float(parts[7]), float(parts[8])], [float(parts[9]), float(parts[10]), float(parts[11])]])
                eps_r = float(parts[12])
                mu_r = float(parts[13])
                ele_size = float(parts[14])
                screen = Screen(corner_positions + self.overall_position, eps_r, mu_r, ele_size)
                screens.append(screen)

        self.scat_objs.extend(screens)

    def _extract_spheres(self) -> None:
        """
        Extends self.scat_objs with a list of Spheres based on the information in the Spheres.txt file.
        """

        spheres_path = self.file_path / "Standard Objects" / "Spheres.txt"

        spheres = []

        with open(spheres_path, "r") as f:
            lines = f.read().splitlines()[1:]
            for line in lines:
                parts = line.split()
                position = np.array([float(parts[0]), float(parts[1]), float(parts[2])])
                eps_r = float(parts[3])
                mu_r = float(parts[4])
                radius = float(parts[5])
                ele_size = float(parts[6])
                sphere = Sphere(position + self.overall_position, eps_r, mu_r, radius, ele_size)
                spheres.append(sphere)

        self.scat_objs.extend(spheres)

    def _extract_user_objs(self) -> None:
        """
        Extends self.scat_objs with a list of ScatteringObject based on the objects defined in the User Objects folder within the Scattering Objects folder.
        """

        user_objects_path = self.file_path / "User Objects"

        user_objects = []

        for txt_file in user_objects_path.glob("*.txt"):
            with open(txt_file, "r") as f:

                lines = f.read().splitlines()[1:]
                file_name = txt_file.stem

                for path in user_objects_path.glob(f"{file_name}.*"):
                    if path != txt_file:
                        obj_geom_file = path
                        break
                else:
                    raise FileNotFoundError(f"Unable to find the geometry file corresponding to {txt_file.name}.")

                for line in lines:
                    parts = line.split()

                    position = np.array([float(parts[0]), float(parts[1]), float(parts[2])])
                    eps_r = float(parts[3])
                    mu_r = float(parts[4])

                    scat_object = ScatteringObject(position + self.overall_position, eps_r, mu_r, obj_geom_file)
                    user_objects.append(scat_object)

        self.scat_objs.extend(user_objects)


def plane_wave(self, point: NDArray[np.float64]) -> NDArray[np.complex64]:
    """
    Returns the value of the incident plane wave at a given point in space. The incident wave is not the physical wave,
    but rather the frequency domain representation of the wave, meaning that it will have complex components.
    """
    return INC_POL * np.exp(1j * K * np.dot(point, INC_DIR))


# The @bempp_cl.api.complex_callable decorator needs to be used for both tangential and Neumann traces in bempp.
@bempp_cl.api.complex_callable
def tangential_trace(point, n, domain_index, result):
    """
    Returns the tangential trace of the incident plane wave at a given point on the surface of one of the scattering objects.
    n is a unit normal vector to the surface at the given point.
    """
    value = INC_POL * np.exp(1j * K * np.dot(point, INC_DIR))
    result[:] = np.cross(value, n)


@bempp_cl.api.complex_callable
def neumann_trace(point, n, domain_index, result):
    """
    Returns the Neumann trace of the incident plane wave at a given point on the surface of one of the scattering objects.
    n is a unit normal vector to the surface at the given point.
    The formula below might look a bit strange, but it is the correct formula for the Neumann trace of a plane wave
    when you compute it by hand.
    """
    value = np.cross(INC_DIR, INC_POL) * 1j * K * np.exp(1j * K * np.dot(point, INC_DIR))
    result[:] = 1.0 / (1j * K) * np.cross(value, n)


class Simulation:
    """
    Simulation class which contains all the information and functions needed for the conducting the bempp simulation.

    Attributes:
    - unit: The unit of length used in the simulation. Must be one of: 'm', 'dm', 'cm', 'mm', 'um', 'nm'.
    - scat_objs: A ScatteringObjectCollection which contains all the scattering objects in the simulation.
    - solution: A list of bempp_cl.api.GridFunction instances which will contain the simulation solution grid functions.
    - sim_name: The name of the simulation, used for saving outputs.
    - GMRES: If True, solver will use the GMRES iterative solver. If false, solver will use a scipy direct solver (scipy.linalg.solve or scipy.linalg.lu_solve).
    - max_iters: The maximum number of iterations for the solver if using the GMRES iterative solver. None if using the LU direct solver.
    - _scat_objs_grids: A list of bempp_cl.api.grid.Grid instances which are the bempp grids for the scattering objects in the simulation.
    - _unit_value: The numerical value of unit in meters.
    """

    unit: str
    scat_objs: ScatteringObjectCollection
    solution: list[bempp_cl.api.GridFunction]
    sim_name: str
    GMRES: bool
    max_iters: int | None
    _scat_objs_grids: list[bempp_cl.api.grid.Grid]
    _unit_value: float

    def __init__(self, simulation_name: str, unit: str, scattering_objs: ScatteringObjectCollection, GMRES: bool, max_iters: int | None) -> None:
        """
        Initializing the simulation attributes. Also finding the bempp grids for the scattering objects in the simulation.
        """
        self.sim_name = simulation_name
        self.unit = unit
        self.scat_objs = scattering_objs
        self.GMRES = GMRES
        self.max_iters = max_iters
        self._scat_objs_grids = []
        for obj in scattering_objs.scat_objs:
            self._scat_objs_grids.append(obj.generate_grid())

        unit_conversion = {"m": 1.0, "dm": 0.1, "cm": 0.01, "mm": 0.001, "um": 1e-6, "nm": 1e-9}
        self._unit_value = unit_conversion[unit]

    def solve(self) -> None:
        """
        This is the main workhorse method of the Simulation class. The majority of this method is dedicated to defining
        the operators and grid functions involved in the main equation which need to be solved to find the scattered field.
        The very end of this method tion will actually solve the system of equations and store the solution in self.solution.
        """

        if len(self.scat_objs.scat_objs) == 1:
            self._solve_single_object()
            return

        grids = self._scat_objs_grids

        print("Defining operators and grid functions...")
        # defining multitrace opterators and the full Generalized Blocked Operator (All are boundary integral operators)
        Ai_int = []
        Ai_ext = []

        i = 0
        for obj in self.scat_objs.scat_objs:
            Ai_int.append(multitrace_operator(grid=grids[i], wavenumber=obj.k_int, epsilon_r=obj.eps_r, mu_r=obj.mu_r, space_type="all_rwg"))
            Ai_ext.append(multitrace_operator(grid=grids[i], wavenumber=K, space_type="all_rwg"))

            i += 1

        Aij = []

        i = 0
        for obj in self.scat_objs.scat_objs:
            j = 0
            Aij.append([])
            for obj in self.scat_objs.scat_objs:
                if i == j:
                    Aij[i].append(Ai_int[i] + Ai_ext[i])
                else:
                    Aij[i].append(multitrace_operator(grids[j], K, target=grids[i], space_type="all_rwg"))  # should check that the ordering of j and i is correct here
                j += 1
            i += 1

        A = bempp_cl.api.GeneralizedBlockedOperator(Aij)

        # defining the right hand side of the systme of linear equations that needs solving
        incident_field_traces = []
        num_trace_entries = 2 * len(grids)

        for i in range(num_trace_entries):
            if i % 2 == 0:
                incident_field_traces.append(bempp_cl.api.GridFunction(space=A.range_spaces[i], dual_space=A.dual_to_range_spaces[i], fun=tangential_trace))
            else:
                incident_field_traces.append(bempp_cl.api.GridFunction(space=A.range_spaces[i], dual_space=A.dual_to_range_spaces[i], fun=neumann_trace))

        # first term in the right hand side of the system of equations
        rhs1 = [0.5 * f for f in incident_field_traces]

        # Defining the diagonal operator Ai on the right hand side of the system of equations
        Ai_diag = []
        num_rows_or_cols = len(grids)

        for i in range(num_rows_or_cols):
            Ai_diag.append([])
            for j in range(num_rows_or_cols):
                if i == j:
                    Ai_diag[i].append(Ai_int[i])
                else:
                    # off diagonal blocks are all zero. but we need to be careful to set the correct spaces. Hence the mess of indices in the ZeroBoundaryOperator parameters.
                    zero_block = BlockedOperator(2, 2)

                    domain_space = Ai_int[j].domain_spaces[0]
                    range_space = Ai_int[i].range_spaces[0]
                    dual_to_range_space = Ai_int[i].dual_to_range_spaces[0]

                    zero_block[0, 0] = ZeroBoundaryOperator(domain_space, range_space, dual_to_range_space)
                    zero_block[0, 1] = ZeroBoundaryOperator(domain_space, range_space, dual_to_range_space)
                    zero_block[1, 0] = ZeroBoundaryOperator(domain_space, range_space, dual_to_range_space)
                    zero_block[1, 1] = ZeroBoundaryOperator(domain_space, range_space, dual_to_range_space)

                    Ai_diag[i].append(zero_block)

        Ai_diag_operator = bempp_cl.api.GeneralizedBlockedOperator(Ai_diag)

        # second term in the right hand side of the system of equations
        rhs2 = Ai_diag_operator * incident_field_traces

        # full right hand side of the system of equations
        rhs = [r1 - r2 for r1, r2 in zip(rhs1, rhs2)]

        print("Solving system of equations...")

        # solving the system of linear equations
        if self.GMRES:
            self.solution = bempp_cl.api.linalg.gmres(A, rhs, return_residuals=True, return_iteration_count=True, maxiter=self.max_iters, use_strong_form=False)[0]
        else:
            self.solution = bempp_cl.api.linalg.lu(A, rhs)

        print("Solution found.")

    def _solve_single_object(self) -> None:
        """
        This is a helper method for the solve() method. This method is a simplified version of the solve() method which is
        used when there is only one scattering object in the simulation. Unfortunately, the solve() method will not work for a single
        object on its own. This method will store the solution in self.solution just like the solve() method does.
        """

        print("Defining object grids...")
        # finding object grid
        obj = self.scat_objs.scat_objs[0]
        grid = obj.generate_grid()

        print("Defining operators and grid functions...")
        # defining multitrace opterators and the full Generalized Blocked Operator (All are boundary integral operators)
        A0_int = multitrace_operator(grid=grid, wavenumber=obj.k_int, epsilon_r=obj.eps_r, mu_r=obj.mu_r, space_type="all_rwg")
        A0_ext = multitrace_operator(grid=grid, wavenumber=K, space_type="all_rwg")

        A = A0_int + A0_ext

        # defining the right hand side of the systme of linear equations that needs solving
        incident_field_traces = [bempp_cl.api.GridFunction(space=A.range_spaces[0], dual_space=A.dual_to_range_spaces[0], fun=tangential_trace), bempp_cl.api.GridFunction(space=A.range_spaces[1], dual_space=A.dual_to_range_spaces[1], fun=neumann_trace)]

        # first term in the right hand side of the system of equations
        rhs1 = [0.5 * f for f in incident_field_traces]

        # second term in the right hand side of the system of equations
        rhs2 = A0_int * incident_field_traces

        # full right hand side of the system of equations
        rhs = [r1 - r2 for r1, r2 in zip(rhs1, rhs2)]

        print("Solving system of equations...")

        # solving the system of linear equations
        if self.GMRES:
            self.solution = bempp_cl.api.linalg.gmres(A, rhs, return_residuals=True, return_iteration_count=True, maxiter=self.max_iters, use_strong_form=False)[0]
        else:
            self.solution = bempp_cl.api.linalg.lu(A, rhs)

        print("Solution found.")

    def plot_sensor_intensity(self, sensors: list[Sensor], signal_boost: float) -> None:
        """
        This method will use the solution stored in self.solution and the given sensors to compute the scattered field intensity on the plane of each sensor.
        The method will use the magnetic and electric field potential operators to evaluate the scattered field at all points of each sensor and then compute the intensity at those points
        from there.

        It should be noted that only the scattered electric field is computed. Technically, both the electric and magnetic fields are needed to find the intensity,
        but we can use the Silver-Müller radiation condition to find a very good approximation for the intensity using only the electric field. This approximation gets
        better the further away the sensor is from the scattering objects.
        """

        print("Extracting E and B fields and plotting...")

        images_folder = root_directory / "Outputs"
        subfolder_name = self.sim_name
        sim_folder = images_folder / subfolder_name
        sim_folder.mkdir(exist_ok=True)
        base_batch_folder_name = "Batch"
        batch_folder = sim_folder / f"{base_batch_folder_name}_1"
        counter = 1
        while batch_folder.exists():
            counter += 1
            batch_folder = sim_folder / f"{base_batch_folder_name}_{counter}"

        batch_folder.mkdir(exist_ok=True)

        for sensor in sensors:
            # finding the evaluation points
            evaluation_points = sensor.points
            # .T This is needed specifically for bempp. Without, will encounter errors later.
            evaluation_points_transpose = sensor.points.T

            # finding the potential operators which provide the solution in the bulk
            potential_operators = []
            num_solution_entries = len(self.solution)

            for i in range(num_solution_entries):
                if i % 2 == 0:
                    potential_operators.append(bempp_cl.api.operators.potential.maxwell.magnetic_field(self.solution[i].space, evaluation_points_transpose, K))
                else:
                    potential_operators.append(bempp_cl.api.operators.potential.maxwell.electric_field(self.solution[i].space, evaluation_points_transpose, K))

            # computing the complex E field
            complex_E_field = -potential_operators[0] * self.solution[0]

            for i in range(1, num_solution_entries):
                complex_E_field -= potential_operators[i] * self.solution[i]

            E_magnitudes_squared = np.sum(np.abs(complex_E_field.T) ** 2, axis=1)

            evaluation_points_mags = np.linalg.norm(evaluation_points, axis=1)

            r_hats = evaluation_points / evaluation_points_mags[:, np.newaxis]

            dot_prod_scaling = np.dot(r_hats, sensor.unit_normal)

            scaling = 1 / (2 * C * MU_0)

            intensity = signal_boost * scaling * E_magnitudes_squared * dot_prod_scaling

            # Show the resulting images
            scattered_image = intensity.reshape(sensor.width_samples, sensor.height_samples).T
            fig, axes = plt.subplots()

            f0 = axes.imshow(scattered_image, origin="lower", cmap="magma", extent=(-sensor.width / 2, sensor.width / 2, -sensor.height / 2, sensor.height / 2), vmin=0, vmax=max(intensity))

            axes.set_title(f"Scattered Field Intensity (Watts/m^2)\nSensor distance: {sensor.distance}{self.unit}")
            axes.set_xlabel(f"Sensor X ({self.unit})")
            axes.set_ylabel(f"Sensor Y ({self.unit})")
            divider = make_axes_locatable(axes)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            fig.colorbar(f0, cax=cax)

            # Saving the plot
            base_filename = "sensor_intensity_plot"
            ext = ".png"
            file_path = batch_folder / f"{base_filename}_1{ext}"
            counter = 1
            while file_path.exists():
                counter += 1
                file_path = batch_folder / f"{base_filename}_{counter}{ext}"

            plt.savefig(file_path)

        info_file = batch_folder / "Batch parameter info.txt"

        with open(info_file, "w") as f:
            f.write(f"Incident wavelength: {INC_WL}{self.unit}\n")
            f.write(f"Incident polarization: {INC_POL}\n")
            f.write(f"Incident direction: {INC_DIR}\n")
            f.write(f"Wave number in medium: {K}{self.unit}^-1\n")
            f.write(f"Unsing GMRES: {self.GMRES}\n")

    def plot_scattering_objs(self) -> None:
        """
        This method plots all the scattering objects involved in the simulation using bempp's built-in plotting functionality.
        The plot is done in gmsh.
        """

        full_grid = bempp_cl.api.grid.union(self._scat_objs_grids)
        full_grid.plot()


class ProgramTimer:
    """
    For timing the main program.
    """

    _start_time: float

    def start(self) -> None:
        self._start_time = time.perf_counter()

    def stop(self) -> None:
        end = time.perf_counter()
        elapsed = end - self._start_time  # seconds
        minutes = elapsed // 60
        seconds = round(elapsed - 60 * minutes)
        print(f"Execution time: {minutes}m {seconds}s")
