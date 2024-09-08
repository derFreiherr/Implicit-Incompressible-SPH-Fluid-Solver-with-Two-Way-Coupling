This is a Implicit Incompressible SPH Fluid Solver with Two Way Coupling.
The visualization is based on the Particles / Instancing tutorial https://www.opengl-tutorial.org/intermediate-tutorials/billboards-particles/particles-instancing/. 
I used the following libraries in the simulator: ImGui, Eigen, OMP, Glew, GLFW3, GLM, and the standard libraries.
The logic of this porgramm is based on the Simulation in Computer Graphics Rigid Bodies course slides from Matthias Teschner. It is aswell based on https://sph-tutorial.physics-simulation.org/pdf/SPH_Tutorial.pdf. 
The user can chose between iisph and sesph using either an uniformgrid neighbour search or an hashmap neighbour search. The user can aswell chose between sph extrapolated boundaries or pressure mirroring boundaries.
There are 11 Scenarios to chose from. technical and visual settings can be setted. The overall stats and the stats of an single particle can be shown. Several stats, or the whole scene can be exported. A scene can be imported aswell.
With W,A,S,D and the arrow keys the camera position can be moved. With the key C the mouse is centered. Wihle pressing the key v or the right mouse button, the view angle can be changed using the mouse. With pressung ESC the programm will close. 
