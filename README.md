<h1>Thesis Project</h1>

A computational model of stochastic field dynamics. This was submitted as my thesis for BSc Physics at Newcastle University. I received 80%, the highest marked thesis project in my year group, and a £100 award from the School of Mathematics, Statistics and Physics.

The model simulates the potential trajectory of a scalar field by numerically solving a Hamilton-Jacobi PDE using the Euler method. Quantum fluctations in potential are modelled by a stochastic term, creating a Langevin equation. Trajectories are simulated until the field reaches a specified exit potential.

A probability distribution function (PDF) of exit times is obtained by performing 10^7 simulations. Multiple PDFs are then obtained by repeating the process over the desired parameter space - in this case, plateau length (a characteristic of the potential profile) is varied. Statistical methods are used to analyse the data.

In modern physical cosmology, the mass fraction of primordial black holes (PBHs) is a quantity directly tied to the time the universe spent undergoing inflation. The profile of the PDFs (specifically, the tail area) is used to calculate an upper bound of the PBH mass fraction. In my report, I used this calculation along with published observed upper bounds of the PBH mass fraction to constrain the plateau length of the inflaton field potential profile.
