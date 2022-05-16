//
// Created by Moritz Gnisia on 05.04.20.
//

#ifndef CFDLAB_ENUMS_H
#define CFDLAB_ENUMS_H

enum class velocity_type {
    U,
    V
};

enum class border_position {
    TOP,
    BOTTOM,
    LEFT,
    RIGHT
};

enum class matrix_selection {
    ROW,
    COLUMN
};

enum class read_type {
    INT,
    DOUBLE
};

enum class boundary_type{
 F         = 1,
 NS_BN     = 34,
 NS_BS     = 66,
 NS_BW     = 130,
 NS_BE     = 258,
 NS_BN_BW  = 162,
 NS_BN_BE  = 290,
 NS_BS_BW  = 194,
 NS_BS_BE  = 322,
 FS_BN     = 36,
 FS_BS     = 68,
 FS_BW     = 132,
 FS_BE     = 260,
 FS_BN_BW  = 164,
 FS_BN_BE  = 292,
 FS_BS_BW  = 196,
 FS_BS_BE  = 324,
 OF_BN     = 40,
 OF_BS     = 72,
 OF_BW     = 136,
 OF_BE     = 264,
 OF_BN_BW  = 168,
 OF_BN_BE  = 296,
 OF_BS_BW  = 200,
 OF_BS_BE  = 328,
 IF_BN     = 48,
 IF_BS     = 80,
 IF_BW     = 144,
 IF_BE     = 272,
 IF_BN_BW  = 176,
 IF_BN_BE  = 304,
 IF_BS_BW  = 208,
 IF_BS_BE  = 336
      
};

enum class log_file{

log_dt = 0,
log_U_Input = 1,
log_V_Input = 2,
log_P_Input = 3,
log_T_Input = 4,
log_F = 5,
log_G = 6,
log_RS = 7,
log_U_Output = 8,
log_V_Output = 9,
log_P_Output = 10,
log_T_Output = 11
};

#endif //CFDLAB_ENUMS_H
