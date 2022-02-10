# Dhruv Jain - Circular Restircted Three Body Problem 

I hope to add the capability to investiage the dynamical flows in CR3BP for Poliastro users. The current files are capable of: 
1. cr3bp_char_quant.py: Computing the characterisitic quantities of a system in CR3BP: mu,l*, t*
2. cr3bp_lib_JC_calc.py: Computing the 5 libration points position and Jacobi Constant
3. cr3bp_master.py: Numerical Integration of CR3BP EOMs (with the option of using events functions)
4. cr3bp_master.py: Numerical Integration of CR3BP + State Transition Matrix EOMs (with the option of using events functions)
5. cr3bp_master.py: Computing the first derivative of pseudo potenial, second-derivative of pseudo potenial and acceleration terms

I will also be working to leverage the preexisisting functions of poliastro to further streamline these features and tweak cr3bp_master.py to be a class. 

The current work requires the widely used scipy, numpy and matplotlib libraries. 

## Future Work: 
I want to make the entire setup more robus. Create a good starting point for users to compute families of periodic orbits and understand the general dynamical flows in CR3BP
1. General purpose Multiple shooter to compute collinear libration point periodic orbit families (Planar and Spatial)
2. Numerical Continuation schemes: Natural parameter continuation, Pseudo-arc length continuation
3. Stability analysis: Sorting Eigenevalues and Eigenvectors of Monodromy matrix, Stability Index plots, Broucke Stability Diagram
4. Manifold generation 
5. Heteroclinic transfer design using Tau-Alpha method
6. Method to transition and corret in N-body ephemeris model

My dream is to build enough features to enable me to include a targeter that I have created to compute Quasi-Periodic Orbits(QPOs) in CR3BP, so that users can access the less investaged but more useful QPOs. 

## About Me:
I am Dhruv Jain, dhruvj9922@gmail.com, pursuing Master of Science in Aeroanautics and Astronautics Engienering at Purdue University. I am part of the Multi-Body Dynamics Research Group and I am working on designing Quasi-Periodic Orbits in CR3BP. When I am not working, I like playing Badminton and cooking!
