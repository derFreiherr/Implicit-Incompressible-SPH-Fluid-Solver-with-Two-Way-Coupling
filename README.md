This is an Implicit Incompressible SPH Fluid Solver with Two-Way Coupling. The visualization is based on the Particles/Instancing tutorial from OpenGL Tutorial https://www.opengl-tutorial.org/intermediate-tutorials/billboards-particles/particles-instancing/. I used the following libraries in the simulator: ImGui, Eigen, OMP, GLEW, GLFW3, GLM, and the standard libraries.

The logic of this program is based on the "Simulation in Computer Graphics" Rigid Bodies course slides from Matthias Teschner. It is also based on the SPH Tutorial https://sph-tutorial.physics-simulation.org/pdf/SPH_Tutorial.pdf. The user can choose between IISPH and SESPH, using either a uniform grid neighbor search or a hashmap neighbor search. The user can also choose between SPH extrapolated boundaries or pressure mirroring boundaries.

There are 11 scenarios to choose from. Technical and visual settings can be configured. The overall stats and the stats of a single particle can be displayed. Several stats, or the entire scene, can be exported, and a scene can be imported as well.

The camera position can be moved using W, A, S, D, and the arrow keys. Pressing the C key centers the mouse. While pressing the V key or the right mouse button, the view angle can be adjusted using the mouse. Pressing ESC will close the program.
