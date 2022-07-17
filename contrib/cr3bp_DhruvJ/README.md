# Dhruv Jain - Circular Restircted Three Body Problem 

Goal: To add the capability to investiage the dynamical flows in CR3BP to Poliastro.

**Current Features:**
1. Compute the characterisitic quantities of a system in CR3BP: mu,l*, t* => _cr3bp_char_quant.py_
2. Computing the 5 libration points position => _cr3bp_lib_calc.py_
3. Numerical Integration of CR3BP EOMs (with the option of using events functions => _cr3bp_master.py_
4. Numerical Integration of CR3BP + State Transition Matrix EOMs (with the option of using events functions) => _cr3bp_model_master.py_
5. Compute Jacobi Cosntant => _cr3bp_model_master.py_
6. Computing the first derivative of pseudo potenial, second-derivative of pseudo potenial and acceleration terms => _cr3bp_master.py_
7. Periodic Orbit Single Shooter Targeter => _cr3bp_PO_master.py_
   * Target states, time, JC, pseudo-arc length constraint
   * Can exploit XZ plane symmetry, X-axis symmetry, Periodicity
   * Phase constraint for PALC with Periodicity Targeter
   * Compute local and global manifolds of a periodic orbit
9. Periodic Orbit Multiple Shooter node setup => _cr3bp_PO_master.py_
   * Place nodes after equal time intervals, or equal number of integrated time-steps
10. Periodic Orbit Family Computation => _cr3bp_po_fam_continuation.py_
   a. Natural Parameter Continuation:  
   	* Use NPC to compute a family of Periodic Orbits
	* line search: To update step size if targeter unable to converge with previous step size
	* Capable of continuing in **'jc'** even if not expliclty defined as a constraint but stated as the _param_conti_
   b. Pseudo-Arc Length Parameter Continuation:  
   	* Use PALC to compute a family of Periodic Orbits
	* line search: To update step size if targeter unable to converge with previous step size
9. Plots family of periodic orbtis => _cr3bp_po_plot_orbits.py_

There are multiple examples included in the PR. The examples are meant to showcase the various cases that the robust targeter and continuation methods can handle. 

The current work requires the widely used numpy, scipy, plotly and matplotlib libraries. 

## Future Work: 
I hope to make a robuts setup for people to play with and understand periodic orbits and the general dynamical flows in CR3BP.

1. Extend single shooter to be a multiple shooter
2. Incorporate unit tests for regression testing 
3. Stability analysis: Sorting Eigenevalues and Eigenvectors of Monodromy matrix, Stability Index plots, Broucke Stability Diagram
4. Heteroclinic transfer design using Tau-Alpha method

My long term goal is to build be a capable setup that can compute Quasi-Periodic Orbits(QPOs) in CR3BP, so that users can access the less investaged but more useful QPOs. 

## About Me:
I am Dhruv Jain, dhruvj9922@gmail.com, pursuing Master of Science in Aeroanautics and Astronautics Engienering at Purdue University. I am part of the Multi-Body Dynamics Research Group and I am working on designing Quasi-Periodic Orbits in CR3BP. When I am not working, I like playing Badminton and cooking!
