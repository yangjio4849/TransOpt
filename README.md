The package TransOpt makes it possible for VASP users to calculate electrical transport properties (Seebeck coefficients, electrical conductivities, and electronic thermal conductivities) by using (1) the momentum matrix method (Method 1) or (2) derivative method same as adopted in BoltzTrap (Method 2). 
The advantage of Method 1 is the overcome of "band-crossing" problem, especially useful for supercell calculations. However, Method 1 cannot deal with the calculations with spin-orbital coupling.
Method 2 can treat systems with/without spin-orbital coupling.
In 2018, the treatment of relaxation time is implemented in the code (compatible with both Method 1 and 2). 
The interface with Quantum Espresso is done in 2019.
In 2022, TransOpt has been updated to Version 2.0. The ionized impurity scattering has been introduced in additional to the previous deformation potential method for acoustic phonon scattering. Besides, one can type in the deformation potentials for both VBM and CBM simultaneously, and TransOpt will calculate the transport properties considering different deformation potentials at the band edges. The new input file, TransOpt.input, is adopted.
