#!/bin/bash 



# Running The simulation of the new Code
../build/sim -c


diff "../Temp_Results/Cavity.10000.vtk" "../Ref_Results/WS_1/cavity100.10000.vtk"
echo "System test finished successfully"