# User Friendly BEM Scattering (UFBS)
## About
User Friendly BEM Scattering (UFBS for short) is a program designed to make it easy to simulate scattering of Electro-Magnetic waves off of any number of arbitrarily shaped objects with different indices of refraction. The program is built upon Bempp: https://bempp.com/, which is an "open-source computational boundary element platform to solve electrostatic, acoustic and electromagnetic problems" as described on their website home page. UFBS is really a wrapper around some of Bempp's functionalities which make it easier for users to get running scattering simulations quickly.

### Constraints on USBF simulations:
- Each scattering object must have a uniform relative permittivity and permeability (and therefore index of refraction as well) within it. But each object can have a unique relative permittivity and permeability.
- Incident light can only be a plane-wave as of right now.

### Workflow:
For an in-depth explanation about how to use UFBS, please see the "How to use the program" section in the Instruction Manual file in the docs, folder or the Example 1 and 2 files also in the docs folder. As a start however, here is a general overview of the workflow involved in running a simulation using UFBS:

1. Information about the scattering objects can be placed in the Scattering Objects folder. This folder has two sub-folders, Standard Objects where the user can tell the program which pre-made scattering objects they would like to use, and User Objects where the user can give the program new scattering objects they've made.
2. Next, inside the Incident_Wave_Params file the user can specify information about the incident plane wave involved in the simulation.
3. The Main.py file serves as the main entry point of the program. In Main.py, the user should provide the final pieces of information that the program needs to run the simulation.
4. Once steps 1 - 3 are complete, the user can run the simulation by running the Main.py file. The output of the simulation will be contained in the Outputs folder.

## Installation
To begin the installation, we need to clone this repository onto your computer. First, open your terminal and navigate to a folder where you want to put this repository by entering the following in your terminal: 
```bash
cd <path-to-your-folder>
```
Once in the desired folder, enter the following command into your terminal to clone this repository into that folder:
```bash
git clone https://github.com/AidanD-C/User-Friendly-BEM-Scattering
```
Next, you might want to set up a python virtual environment in the User-Friendly-BEM-Scattering folder for dependency control, but this is optional. After that, you will want to download all the necessary dependencies for this repository. These dependencies can be found in the requirements.txt folder. After executing the command above, these next two commands can be used to download the dependencies in requirements.txt using pip:
```bash
cd User-Friendly-BEM-Scattering

# activate your virtual environment here if you chose to use one.

pip install -r requirements.txt
```

In order to visualize scattering objects, you will also need to download gmsh. This is different from using pip to install gmsh since what we're really after is the gmsh.exe file which is not downloaded through "pip install gmsh". Here's a link to the download page: https://gmsh.info/. On that download page, in the "Current stable release" section, there should be a line that says "Download Gmsh for...". Clicking on whichever platform you are using (Windows, Mac, or Linux) will download a zip file. The contents of that zip file should be extracted into the User-Friendly-BEM-Scattering folder so that it can be found by UFBS. For example, at the time of writing this, using the download link for windows it downloads a gmsh-4.14.0-Windows64.zip file which contains a single folder named gmsh-4.14.0-Windows64. This folder should be placed into the User-Friendly-BEM-Scattering folder so that after all is done, the contents of the User-Friendly-BEM-Scattering folder should look like the following (or something similar since this image shows the layout on windows):
<p align="center">
  <img src="Images/Setup Screenshot.png" alt="Image 1" width="27%">
</p>
Though the .venv folder might not be present if you chose not to use a python virtual environment.

As long as the gmsh folder you extracted from the zip file is somwhere in the User-Friendly-BEM-Scattering folder, UFBS will be able to locate the gmsh.exe file it needs.

## Example outputs:
### Scattering off of a long bacteria:
Gmsh preview:

<p align="center">
  <img src="Images\Long Bacteria.png" alt="Image 2" width="40%">
</p>

Output:
<p align="center">
  <img src="Images\Bacteria Sensor.png" alt="Image 2" width="40%">
</p>

### Scattering off of a single sphere:
Gmsh preview:

<p align="center">
  <img src="Images\Sphere.png" alt="Image 3" width="40%">
</p>

Output:
<p align="center">
  <img src="Images\Sphere Sensor.png" alt="Image 4" width="40%">
</p>

## Notes and warnings on the current state of UFBS:
- UFBS has not yet undergone a rigorous set of tests to ensure its results are physically accurate. It has undergone a mild set of tests to verify its results are reasonable, but not enough tests to be considered complete. The outputs of this software should not yet be trusted completely.
- The ability to efficiently use GMRES in UFBS has not yet been implemented. This would require implementing a preconditioning scheme. Thus, for now, its recommended that the using_GMRES variable in the Main.py file is always set to False.
- UFBS can generally only handle up to about four scattering objects in a single simulation due to GMRES not yet being implemented.
- The latest version of python UFBS has been confirmed to work on is python 3.13.2.
