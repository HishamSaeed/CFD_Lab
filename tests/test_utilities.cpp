#include "test_utilities.hpp"

/*----------Output Directories Handling----------------*/
#pragma region 
void clear_output_dir_test_()
{
    struct stat st = {0};

    if (stat(Cmp_Test_Data_t.c_str(), &st) == -1)
    {
    	mkdir(Cmp_Test_Data_t.c_str(), 0700);
	}
    else
    {
        system(("exec rm -r " + Cmp_Test_Data_t + "*").c_str());
    }
}
#pragma endregion