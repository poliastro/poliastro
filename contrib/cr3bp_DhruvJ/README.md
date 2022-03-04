# Dhruv Jain - Circular Restircted Three Body Problem 

This is a first major update to the proposed contrib/cr3bp_DhruvJ

Goal: To addd the capability to Poliastro investiage the dynamical flows in CR3BP.

**Current Features:**
1. Compute the characterisitic quantities of a system in CR3BP: mu,l*, t* => _cr3bp_char_quant.py_
2. Computing the 5 libration points position and Jacobi Constant => _cr3bp_lib_JC_calc.py_
3. Numerical Integration of CR3BP EOMs (with the option of using events functions => _cr3bp_master.py_
4. Numerical Integration of CR3BP + State Transition Matrix EOMs (with the option of using events functions) => _cr3bp_master.py_
5. Computing the first derivative of pseudo potenial, second-derivative of pseudo potenial and acceleration terms => _cr3bp_master.py_
6. Periodic Orbit Single Shooter Targeter => _cr3bp_PO_targeter.py_
   * Target states, time, JC, pseudo-arc length constraint
   * Can exploit XZ plane symmetry, X-axis symmetry, Periodicity
7. Periodic Orbit Multiple Shooter node setup => _cr3bp_PO_targeter.py_
   * Place nodes after equal time intervals, or equal number of integrated time-steps
8. Periodic Orbit Family Computation => _cr3bp_fam_continuation.py_
   a. Natural Parameter Continuation:  
   	* Use NPC to compute a family of Periodic Orbits
	* line search: To update step size if unable to converge with previous step size
	* Capable of continuing in **'jc'** even if not expliclty defined as a constraint but stated as the _param_conti_
   b. Pseudo-Arc Length Parameter Continuation:  
   	* Use PALC to compute a family of Periodic Orbits
	* line search: To update step size if unable to converge with previous step size
9. Plots family of periodic orbtis => _cr3bp_plot_orbits.py_

The current work requires the widely used numpy, scipy, plotly and matplotlib libraries. 

## Future Work: 
I hope to make a robuts setup for people to play around and undersntand periodic orbits and the general dynamical flows in CR3BP

1. Extend single shooter to be a multiple shooter
2. Stability analysis: Sorting Eigenevalues and Eigenvectors of Monodromy matrix, Stability Index plots, Broucke Stability Diagram
3. Manifold generation 
4. Heteroclinic transfer design using Tau-Alpha method
5. Method to transition and corret in N-body ephemeris model

My long term goal is to build be a capable setup that can compute Quasi-Periodic Orbits(QPOs) in CR3BP, so that users can access the less investaged but more useful QPOs. 

## About Me:
I am Dhruv Jain, dhruvj9922@gmail.com, pursuing Master of Science in Aeroanautics and Astronautics Engienering at Purdue University. I am part of the Multi-Body Dynamics Research Group and I am working on designing Quasi-Periodic Orbits in CR3BP. When I am not working, I like playing Badminton and cooking!
