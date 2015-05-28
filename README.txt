# Siwir2_Ex1
Run the following commands for specific actions:

To compile evrything: make all
To run for one grid size 1/256 for 10 V-Cycles: make onegrid
To run for multiple grid sizes ranging from {1/8, 1/16,....,1/256} for 20 V-Cycles : make allgrid

It will give output for both Dirichlet and Neumann boundary condition problems and display Residual L2 Norm and Convergence rate for each V-Cycle and a final Error L2 Norm and time spent for each of the 2 problems.

After any execution datafiles and surface plots will be stored in following folders,

./data/Dirichlet/
./data/Neumann/

Please remember all datafiles from previous execution will be deleted.

Please find Siwir2_Ex1_Report.pdf for the detail data analysis. 
