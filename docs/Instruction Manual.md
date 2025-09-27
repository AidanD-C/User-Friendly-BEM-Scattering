User-Friendly BEM Scattering (UFBS) User Manual
==============

**Author:** Aidan Downie-Cheetham

**Contents:**
1. [Introduction](#1-introduction)
2. [How to use the program](#2-how-to-use-the-program)
3. [Troubleshooting: common issues](#3-troubleshooting-common-issues)
4. [Theory Background](#4-theory-background)
5. [How UFBS works on the inside](#5-how-ufbs-works-on-the-inside)
6. [What still needs to be done](#6-what-still-needs-to-be-done)

---

## 1. Introduction
The main goal of this program is to find out how light scatters off any number of arbitrarily shaped objects. The program is based on the Boundary Element Method (BEM), a general method which can be used to solve Maxwell's equations. The computational backbone of the program is Bempp-cl, an open-source software package which can be used to solve a number of different Partial Differential Equation (PDE) problems using the BEM, and which can be specifically applied to solving Maxwell's equations.

Along with solving Maxwell's equations, another goal of this program is for it to be easy to use for non-physicists and non-mathematicians. Using this program involves a straighforward workflow which a user can follow to conduct their scattering simulations.

To be completed...

---

# 2. How to use the program
At a high level, all the user needs to do is give the program the information about the scattering objects (which will be located in the Scattering Objects folder) and give the program the information about the incident light (located in the Incident_Wave_Params.py file) and then click run. If everything goes well, the user should see the output of the simulation in the Outputs folder.

To be completed...

---

# 3. Troubleshooting: common issues
To be completed...

---

# 4. Theory background
Mainly following the theory of this paper: https://www.sciencedirect.com/science/article/abs/pii/S0022407318306320

uses the ... approximation to find the time averaged poynting vector

To be completed...


---

# 5. How UFBS works on the inside
Bempp on its own is a pretty general and mathematics-oriented solver. For example, if you just use Bempp on its own, you must define objects such as "boundary integral operators" and "grid functions" which live on your scattering objects. All this detail is unnecessary for the user who just wants to set up a scattering simulation.

To be completed...

---

# 6. What still needs to be done
- Still need to thoroughly test the software to see how realistic the outputs are. I have run some test simulations and they look pretty good.
- It would be good to add in a model for speckle
- All the equations of the main simulation section and the plotting section should still be triple checked.
- This does not handle absorbing materials or conductors
- Preconditioning
- Capability for non-plane wave incident light
- General rule for size of elements compared to wave length and object size
- Only going to work for evaluation of fields at exterior points.
- Element size is very small could cause system of equations solver to be extremely slow.