#include <string>
#include <sys/stat.h>
#include <sys/types.h>

/*------------Directory paths for the Folder structure in the project-------------*/
// Paths are relative to the build path (/build/test) of the test executable not to the solver executable(build)
#pragma region 
// Directory of Testing Data, log files for matrices
const std::string Ref_Test_Data_t = "../../Ref_Test_Data/";
// const std::string Ref_Test_Data_t  = "../../Testing_Data/";
const std::string Temp_Test_Data_t = "../../Temp_Test_Data/";
const std::string Cmp_Test_Data_t  = "../../Cmp_Test_Data/";

// Directory of Input Data
const std::string Sim_Data_t = "../../Simulation_Data/";
const std::string Deb_Data_t = "../../Debugging_Data/";
const std::string Test_Data_t = "../../Testing_Data/";

// Directory of Results in vtk format
const std::string Ref_Results_t = "../../Ref_Results/";
const std::string Temp_Results_t = "../../Temp_Results/";
#pragma endregion

/*--------------------- Log file paths-------------------------------------------*/
#pragma region 
const std::string ref_log_U_t = std::string(Ref_Test_Data_t) + "log_U";
const std::string ref_log_V_t = std::string(Ref_Test_Data_t) + "log_V";
const std::string ref_log_F_t = std::string(Ref_Test_Data_t) + "log_F";
const std::string ref_log_G_t = std::string(Ref_Test_Data_t) + "log_G";
const std::string ref_log_P_t = std::string(Ref_Test_Data_t) + "log_P";
const std::string ref_log_RS_t = std::string(Ref_Test_Data_t) + "log_RS";
const std::string ref_log_dt_t = std::string(Ref_Test_Data_t) + "log_dt";
#pragma endregion

/*--------------------- comparision file paths-------------------------------------------*/
#pragma region 
const std::string cmp_log_U_t = std::string(Cmp_Test_Data_t) + "log_U";
const std::string cmp_log_V_t = std::string(Cmp_Test_Data_t) + "log_V";
const std::string cmp_log_F_t = std::string(Cmp_Test_Data_t) + "log_F";
const std::string cmp_log_G_t = std::string(Cmp_Test_Data_t) + "log_G";
const std::string cmp_log_P_t = std::string(Cmp_Test_Data_t) + "log_P";
const std::string cmp_log_T_t = std::string(Cmp_Test_Data_t) + "log_T";
const std::string cmp_log_RS_t = std::string(Cmp_Test_Data_t) + "log_RS";
const std::string cmp_log_dt_t = std::string(Cmp_Test_Data_t) + "log_dt";
#pragma endregion

/*--------------------- comparision file paths-------------------------------------------*/
#pragma region 
const std::string Lid_Driven_Cavity_t = std::string(Sim_Data_t) + "Cavity";
const std::string Plane_Shear_Flow_t = std::string(Sim_Data_t) + "Plane_Shear_Flow";
const std::string Flow_Over_A_Step_t = std::string(Sim_Data_t) + "Flow_Over_A_Step";
const std::string Karman_Vortex_Street_t = std::string(Sim_Data_t) + "Karman_Vortex_Street";
const std::string Natural_Connvection_t = std::string(Sim_Data_t) + "Natural_Convection";
const std::string Fluid_Trap_t = std::string(Sim_Data_t) + "Fluid_Trap";
const std::string RB_Convection_t = std::string(Sim_Data_t) + "RB_Convection";
#pragma endregion

/*----------Output Directories Handling----------------*/
#pragma region 
void clear_output_dir_test_();
#pragma endregion