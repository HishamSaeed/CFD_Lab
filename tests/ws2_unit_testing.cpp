// Definition of sample test
#define CATCH_CONFIG_RUNNER

#include "catch.hpp"
#include "../init.hpp"
#include "../grid.hpp"
#include "../utilities.hpp"
#include "../boundary_val.hpp"
#include "../uvp.hpp"
#include "../sor.hpp"
#include "../enums.hpp"
#include "test_utilities.hpp"

int main( int argc, char* argv[] ) {
  clear_output_dir_test();  
  int result = Catch::Session().run( argc, argv );

  return result;
}


/*-------------------------------read functions unit Tests----------------------------------------*/
#pragma region 
TEST_CASE("Test grid parameters", "[read_grid_parameters] [!mayfail]"){
    
    SECTION("LID_Driven_Cavity")
    {
      // Dummy Grid Paramters testing1
      #pragma region 
      
      int imax = 10;
      int jmax = 10;

      double UI = 0.0;
      double VI = 0.0;
      double PI = 0.0;
      double TI = 0.0;
      int** geometry;
      #pragma endregion

      // Geometry Initialization
      #pragma region 
      geometry = imatrix(0, imax+1, 0, jmax+1);

      //Setting upper and lower boundary
      for(int i = 0; i <= imax+1; i++)
      {
        geometry[i][0] = 0;
        geometry[i][jmax+1] = 0;
      }
      // Setting righ and left boundary
      for(int j = 0; j <= jmax+1; j++)
      {
        geometry[0][j] = 0 ;
        geometry[imax+1][j] = 0;
      }

      // Setting inner Cell
      for(int i = 1; i <= imax; i++)
      { 
        for(int j = 1; j <= jmax; j++)
        {
          geometry[i][j] = 4;
          
        }
      } 
      #pragma endregion
      
      // Initialize Dummy Grid
      Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry);
      /** Possible Cases For Cell Flag
       * F        f = 1
       * NS/BN    f = 34
       * NS/BS    f = 66
       * NS/BW    f = 130
       * NS/BE    f = 258
       * NS/BN/BW f = 
       * NS/BN/BE f = 
       * NS/BS/BW f =
       * NS/BS/BE f =
       * FS/BN    f = 36
       * FS/BS    f = 68
       * FS/BW    f = 132
       * FS/BE    f = 260
       * FS/BN/BW f = 
       * FS/BN/BE f = 
       * FS/BS/BW f =
       * FS/BS/BE f =
       * OF/BN    f = 40
       * OF/BS    f = 72
       * OF/BW    f = 136
       * OF/BE    f = 264
       * OF/BN/BW f = 
       * OF/BN/BE f = 
       * OF/BS/BW f =
       * OF/BS/BE f =
       * IF/BN    f = 48
       * IF/BS    f = 80
       * IF/BW    f = 144
       * IF/BE    f = 272
       * IF/BN/BW f = 
       * IF/BN/BE f = 
       * IF/BS/BW f =
       * IF/BS/BE f =
       * 
       */
      
      // Test Dummy Grid
      #pragma region
      REQUIRE ( domain.boundary_size() == boundary_size );
      REQUIRE ( domain.imax() == imax );
      REQUIRE ( domain.jmax() == jmax );
      REQUIRE ( domain.imaxb() == (jmax + (2*boundary_size)) );
      REQUIRE ( domain.jmaxb() == (jmax + (2*boundary_size)) );

      for (size_t i = 1; i < domain.imaxb()-1; i++)
      {
          for (size_t j = 1; j < domain.jmaxb()-1; j++)
          {
              REQUIRE ( domain.cell(i, j).velocity(velocity_type::U) == UI );
              REQUIRE ( domain.cell(i, j).velocity(velocity_type::V) == VI );
              REQUIRE ( domain.cell(i, j).pressure() == PI );
              REQUIRE ( domain.cell(i, j).flag() == 1);
          }
          
      }

      for (size_t i = 0; i < domain.imaxb(); i++)
      {
          REQUIRE ( domain.cell(i, 0).border(border_position::BOTTOM) == true );
          REQUIRE ( domain.cell(i, 0).flag() == 34);
          REQUIRE ( domain.cell(i, domain.jmaxb()-1).border(border_position::TOP) == true );
          REQUIRE ( domain.cell(i, domain.jmaxb()-1).flag() == 66);
      }
      
      for (size_t j = 0; j < domain.jmaxb(); j++)
      {
          REQUIRE ( domain.cell(0, j).border(border_position::LEFT) == true );
          REQUIRE ( domain.cell(0, j).flag() == 258);
          REQUIRE ( domain.cell(domain.imaxb()-1, j).border(border_position::RIGHT) == true );
          REQUIRE ( domain.cell(domain.imaxb(), j).flag() == 130);
      }
      #pragma endregion     
    }

    // SECTION("LID_Driven_Cavity")
    // {
    //   // Dummy Grid Paramters
    //   #pragma region 
      
    //   int imax = 10;
    //   int jmax = 10;

    //   double UI = 0.0;
    //   double VI = 0.0;
    //   double PI = 0.0;
    //   int** geometry;
    //   #pragma endregion

    //   // Geometry Initialization
    //   #pragma region 
    //   geometry = imatrix(0, imax+1, 0, jmax+1);

    //   //Setting upper and lower boundary
    //   for(int i = 0; i <= imax+1; i++)
    //   {
    //     geometry[i][0] = 0;
    //     geometry[i][jmax+1] = 0;
    //   }
    //   // Setting righ and left boundary
    //   for(int j = 0; j <= jmax+1; j++)
    //   {
    //     geometry[0][j] = 0 ;
    //     geometry[imax+1][j] = 0;
    //   }

    //   // Setting inner Cell
    //   for(int i = 1; i <= imax; i++)
    //   { 
    //     for(int j = 1; j <= jmax; j++)
    //     {
    //       geometry[i][j] = 4;
          
    //     }
    //   } 
    //   #pragma endregion
      
    //   // Initialize Dummy Grid
    //   Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,geometry);
    //   /** Possible Cases For Cell Flag
    //    * F        f = 1
    //    * NS/BN    f = 34
    //    * NS/BS    f = 66
    //    * NS/BW    f = 130
    //    * NS/BE    f = 258
    //    * NS/BN/BW f = 
    //    * NS/BN/BE f = 
    //    * NS/BS/BW f =
    //    * NS/BS/BE f =
    //    * FS/BN    f = 36
    //    * FS/BS    f = 68
    //    * FS/BW    f = 132
    //    * FS/BE    f = 260
    //    * FS/BN/BW f = 
    //    * FS/BN/BE f = 
    //    * FS/BS/BW f =
    //    * FS/BS/BE f =
    //    * OF/BN    f = 40
    //    * OF/BS    f = 72
    //    * OF/BW    f = 136
    //    * OF/BE    f = 264
    //    * OF/BN/BW f = 
    //    * OF/BN/BE f = 
    //    * OF/BS/BW f =
    //    * OF/BS/BE f =
    //    * IF/BN    f = 48
    //    * IF/BS    f = 80
    //    * IF/BW    f = 144
    //    * IF/BE    f = 272
    //    * IF/BN/BW f = 
    //    * IF/BN/BE f = 
    //    * IF/BS/BW f =
    //    * IF/BS/BE f =
    //    * 
    //    */
      
    //   // Test Dummy Grid
    //   #pragma region
    //   REQUIRE ( domain.boundary_size() == boundary_size );
    //   REQUIRE ( domain.imax() == imax );
    //   REQUIRE ( domain.jmax() == jmax );
    //   REQUIRE ( domain.imaxb() == (jmax + (2*boundary_size)) );
    //   REQUIRE ( domain.jmaxb() == (jmax + (2*boundary_size)) );

    //   for (size_t i = 1; i < domain.imaxb()-1; i++)
    //   {
    //       for (size_t j = 1; j < domain.jmaxb()-1; j++)
    //       {
    //           REQUIRE ( domain.cell(i, j).velocity(velocity_type::U) == UI );
    //           REQUIRE ( domain.cell(i, j).velocity(velocity_type::V) == VI );
    //           REQUIRE ( domain.cell(i, j).pressure() == PI );
    //           REQUIRE ( domain.cell(i, j).flag() == 1);
    //       }
          
    //   }

    //   for (size_t i = 0; i < domain.imaxb(); i++)
    //   {
    //       REQUIRE ( domain.cell(i, 0).border(border_position::BOTTOM) == true );
    //       REQUIRE ( domain.cell(i, 0).flag() == 34);
    //       REQUIRE ( domain.cell(i, domain.jmaxb()-1).border(border_position::TOP) == true );
    //       REQUIRE ( domain.cell(i, domain.jmaxb()-1).flag() == 66);
    //   }
      
    //   for (size_t j = 0; j < domain.jmaxb(); j++)
    //   {
    //       REQUIRE ( domain.cell(0, j).border(border_position::LEFT) == true );
    //       REQUIRE ( domain.cell(0, j).flag() == 258);
    //       REQUIRE ( domain.cell(domain.imaxb()-1, j).border(border_position::RIGHT) == true );
    //       REQUIRE ( domain.cell(domain.imaxb(), j).flag() == 130);
    //   }
    //   #pragma endregion     
    // }
}

TEST_CASE( "Test read parameters", "[read_parameters] [!mayfail]" ) {
    
   // Problem Parameters 
  #pragma region 
  // Solvers Parameters, 
  int itermax;
  double eps;
  double alpha;
  double omg;
  double tau;
  double beta;
  
  // Grid Data
  double xlength;
  double ylength;
  double dx;
  double dy;
  int imax;
  int jmax;

  // Fluid Data, reynolds number, initial conditions for velocity(x and y direction) and pressure
  double Re;
  double PR;
  double UI;
  double VI;
  double PI;
  double TI;

  // No Slip Boundary Values
  double UT;
  // Forces Data, gravitational forces
  double GX;  
  double GY;             
  
  // Time Stepping parameters
  double t_end;                
  double dt;             
  double dt_value;
  #pragma endregion

  

  SECTION( "LID_Driven_Cavity" ) 
  {
    read_parameters(Lid_Driven_Cavity_t + ".dat",&Re,&PR,&UI,&VI,&PI,&TI,&GX,&GY,&t_end,&xlength,&ylength,
                  &dt,&dx,&dy,&imax,&jmax,&alpha,&omg,&tau,&itermax,&eps,&beta,&dt_value,&UT);

    CHECK( xlength == 1 );
    CHECK( ylength == 1 );
    CHECK( imax == 50 );
    CHECK( jmax == 50 );
    CHECK( dt == 0.05 );
    CHECK( t_end == 50.0 );
    CHECK( tau == 0.5 );
    CHECK( dt_value == 0.5 );  
    CHECK( itermax == 100 );
    CHECK( eps == 0.001 );
    CHECK( omg == 1.7 );
    CHECK( alpha == 0.5 );
    CHECK( Re == 100 );
    CHECK( GX == 0.0 );
    CHECK( GY == 0.0 );
    CHECK( PI == 0.0 );
    CHECK( UI == 0.0 );
    CHECK( VI == 0.0 );
    CHECK( UT == 1.0 );
  }

  SECTION( "Plane_Shear_Flow" ) 
  {
    read_parameters(Plane_Shear_Flow_t + ".dat",&Re,&PR,&UI,&VI,&PI,&TI,&GX,&GY,&t_end,&xlength,&ylength,
                  &dt,&dx,&dy,&imax,&jmax,&alpha,&omg,&tau,&itermax,&eps,&beta,&dt_value,&UT);

    CHECK( xlength == 10 );
    CHECK( ylength == 2 );
    CHECK( imax == 100 );
    CHECK( jmax == 20 );
    CHECK( dt == 0.05 );
    CHECK( t_end == 30.0 );
    CHECK( tau == 0.5 );
    CHECK( dt_value == 5.0 );  
    CHECK( itermax == 500 );
    CHECK( eps == 0.001 );
    CHECK( omg == 1.7 );
    CHECK( alpha == 0.9 );
    CHECK( Re == 10 );
    CHECK( GX == 0.0 );
    CHECK( GY == 0.0 );
    CHECK( PI == 0.0 );
    CHECK( UI == 1.0 );
    CHECK( VI == 0.0 );
    CHECK( UT == 0.0 );
  }  

  SECTION( "Flow_Over_A_Step" ) 
  {
    read_parameters(Flow_Over_A_Step_t + ".dat",&Re,&PR,&UI,&VI,&PI,&TI,&GX,&GY,&t_end,&xlength,&ylength,
                  &dt,&dx,&dy,&imax,&jmax,&alpha,&omg,&tau,&itermax,&eps,&beta,&dt_value,&UT);

    CHECK( xlength == 10 );
    CHECK( ylength == 2 );
    CHECK( imax == 100 );
    CHECK( jmax == 20 );
    CHECK( dt == 0.05 );
    CHECK( t_end == 500.0 );
    CHECK( tau == 0.5 );
    CHECK( dt_value == 10.0 );  
    CHECK( itermax == 500 );
    CHECK( eps == 0.001 );
    CHECK( omg == 1.7 );
    CHECK( alpha == 0.9 );
    CHECK( Re == 100 );
    CHECK( GX == 0.0 );
    CHECK( GY == 0.0 );
    CHECK( PI == 0.0 );
    CHECK( UI == 0.0 );
    CHECK( VI == 0.0 );
  }

  SECTION( "Karman_Vortex_Street" ) 
  {
    read_parameters(Karman_Vortex_Street_t + ".dat",&Re,&PR,&UI,&VI,&PI,&TI,&GX,&GY,&t_end,&xlength,&ylength,
                  &dt,&dx,&dy,&imax,&jmax,&alpha,&omg,&tau,&itermax,&eps,&beta,&dt_value,&UT);

    CHECK( xlength == 10 );
    CHECK( ylength == 2 );
    CHECK( imax == 100 );
    CHECK( jmax == 20 );
    CHECK( dt == 0.05 );
    CHECK( t_end == 20.0 );
    CHECK( tau == 0.5 );
    CHECK( dt_value == 2.0 );  
    CHECK( itermax == 500 );
    CHECK( eps == 0.001 );
    CHECK( omg == 1.7 );
    CHECK( alpha == 0.9 );
    CHECK( Re == 10000 );
    CHECK( GX == 0.0 );
    CHECK( GY == 0.0 );
    CHECK( PI == 0.0 );
    CHECK( UI == 1.0 );
    CHECK( VI == 0.0 );
  }

  SECTION( "Natural_Connvection" ) 
  {
    read_parameters(Natural_Connvection_t + ".dat",&Re,&PR,&UI,&VI,&PI,&TI,&GX,&GY,&t_end,&xlength,&ylength,
                  &dt,&dx,&dy,&imax,&jmax,&alpha,&omg,&tau,&itermax,&eps,&beta,&dt_value,&UT);

    CHECK( xlength == 1 );
    CHECK( ylength == 1 );
    CHECK( imax == 50 );
    CHECK( jmax == 50 );
    CHECK( dt == 0.05 );
    CHECK( t_end == 1000.0 );
    CHECK( tau == 0.5 );
    CHECK( dt_value == 10.0 );  
    CHECK( itermax == 100 );
    CHECK( eps == 0.00001 );
    CHECK( omg == 1.7 );
    CHECK( alpha == 0.5 );
    CHECK( beta == 0.00021 );
    CHECK( Re == 1000 );
    CHECK( PR == 7 );
    CHECK( GX == 0.0 );
    CHECK( GY == -1.1 );
    CHECK( PI == 0.0 );
    CHECK( UI == 0.0 );
    CHECK( VI == 0.0 );
    CHECK( TI == 0.0 );
  }    

  SECTION( "Fluid_Trap" ) 
  {
    read_parameters(Fluid_Trap_t + ".dat",&Re,&PR,&UI,&VI,&PI,&TI,&GX,&GY,&t_end,&xlength,&ylength,
                  &dt,&dx,&dy,&imax,&jmax,&alpha,&omg,&tau,&itermax,&eps,&beta,&dt_value,&UT);

    CHECK( xlength == 2 );
    CHECK( ylength == 1 );
    CHECK( imax == 100 );
    CHECK( jmax == 50 );
    CHECK( dt == 0.05 );
    CHECK( t_end == 2000.0 );
    CHECK( tau == 0.5 );
    CHECK( dt_value == 10.0 );  
    CHECK( itermax == 1000 );
    CHECK( eps == 0.00001 );
    CHECK( omg == 1.7 );
    CHECK( alpha == 0.5 );
    CHECK( beta == 0.00063 );
    CHECK( Re == 10000 );
    CHECK( PR == 7 );
    CHECK( GX == 0.0 );
    CHECK( GY == -9.81 );
    CHECK( PI == 0.0 );
    CHECK( UI == 0.0 );
    CHECK( VI == 0.0 );
    CHECK( TI == 0.0 );
  } 

  SECTION( "RB_Convection" ) 
  {
    read_parameters(RB_Convection_t + ".dat",&Re,&PR,&UI,&VI,&PI,&TI,&GX,&GY,&t_end,&xlength,&ylength,
                  &dt,&dx,&dy,&imax,&jmax,&alpha,&omg,&tau,&itermax,&eps,&beta,&dt_value,&UT);

    CHECK( xlength == 8.5 );
    CHECK( ylength == 1 );
    CHECK( imax == 85 );
    CHECK( jmax == 18 );
    CHECK( dt == 0.05 );
    CHECK( t_end == 45000.0 );
    CHECK( tau == 0.5 );
    CHECK( dt_value == 100.0 );  
    CHECK( itermax == 100 );
    CHECK( eps == 0.00001 );
    CHECK( omg == 1.7 );
    CHECK( alpha == 0.5 );
    CHECK( beta == 0.000179 );
    CHECK( Re == 33.73 );
    CHECK( PR == 12500 );
    CHECK( GX == 0.0 );
    CHECK( GY == -0.3924 );
    CHECK( PI == 0.0 );
    CHECK( UI == 0.0 );
    CHECK( VI == 0.0 );
    CHECK( TI == 293.0 );
  }  
}

TEST_CASE( "Test read pgm file", "[read_pgm]")
{
  int imax;
  int jmax;
  int** geometry;

  SECTION("LID_DrivenCavity")
  {
    imax = 50;
    jmax = 50;
    geometry = read_pgm((Lid_Driven_Cavity_t + ".pgm").c_str());

    // Testing Horizontal Boundary Cells
    for(int i = 0; i <= imax+1; i++)
    {
      REQUIRE(geometry[i][0] == 0 );
      REQUIRE(geometry[i][jmax+1] == 0);
    }
    // Testing Vertical Boundary Cells
    for(int j = 0; j <= jmax+1; j++)
    {
      REQUIRE(geometry[0][j] == 0 );
      REQUIRE(geometry[imax+1][j] == 0);
    }

    // Testing Inner Cell
    for(int i = 1; i <= imax; i++)
    { 
      for(int j = 1; j <= jmax; j++)
      {
        REQUIRE(geometry[i][j] == 4);
      }
    }  
  }

  SECTION("Plane_Shear_Flow")
  {
    imax = 100;
    jmax = 20;
    geometry = read_pgm((Plane_Shear_Flow_t + ".pgm").c_str());

    // Testing Horizontal Boundary Cells
    for(int i = 0; i <= imax+1; i++)
    {
      REQUIRE(geometry[i][0] == 0 );
      REQUIRE(geometry[i][jmax+1] == 0);
    }
    // Testing Vertical Boundary Cells
    for(int j = 1; j <= jmax; j++)
    {
      REQUIRE(geometry[0][j] == 3 );
      REQUIRE(geometry[imax+1][j] == 2);
    }

    // Testing Inner Cell
    for(int i = 1; i <= imax; i++)
    { 
      for(int j = 1; j <= jmax; j++)
      {
        REQUIRE(geometry[i][j] == 4);
      }
    }  
  }
  

}

#pragma endregion

#pragma region 
TEST_CASE("Test boundary values", "[boundaryvalues] [!mayfail]")
{
    // temp variables
    double ghost_cell;
    double inner_cell;

    // tolerance
    double eps = 0.00001;

    // Dummy Grid
    double PI = 0.0;
    double UI = 3.0;
    double VI = 4.0;
    double TI = 0.0;
    double UT;

    int imax;
    int jmax;
    int** geometry;

    SECTION("LID_DrivenCavity")
    {
      UT = 1.0;
      imax = 10;
      jmax = 10;
      // Geometry initialization
      #pragma region 
      geometry = imatrix(0, imax+1, 0, jmax+1);
      // Setting upper and lower boundary
      for(int i = 0; i <= imax+1; i++)
      {
        geometry[i][0] = 0;
        geometry[i][jmax+1] = 0;
      }
      // Setting right and left boundary
      for(int j = 0; j <= jmax+1; j++)
      {
        geometry[0][j] = 0 ;
        geometry[imax+1][j] = 0;
      }
      // Setting inner Cell
      for(int i = 1; i <= imax; i++)
      { 
        for(int j = 1; j <= jmax; j++)
        {
          geometry[i][j] = 4;
        }
      } 
      #pragma endregion
      // Grid Initilization
      Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry) ;

      boundaryvalues(imax,jmax,domain,UT);

      #pragma region 
      // Testing upper and lower boundary
      for(int i = 1; i <= imax; i++)
    {
        ghost_cell = domain.cell(i,0).velocity(velocity_type::U);
        inner_cell = domain.cell(i,1).velocity(velocity_type::U);
        REQUIRE( fabs(ghost_cell + inner_cell) < eps);

        ghost_cell = domain.cell(i,0).velocity(velocity_type::V);
        REQUIRE( ghost_cell == 0.0);

        ghost_cell = domain.cell(i,jmax+1).velocity(velocity_type::U);
        inner_cell = domain.cell(i,jmax).velocity(velocity_type::U);
        REQUIRE( fabs(ghost_cell + inner_cell) == 2);

        ghost_cell = domain.cell(i,jmax).velocity(velocity_type::V);
        REQUIRE( ghost_cell == 0.0);
    }

      // Testing left and right boundary
      for(int j = 1; j <= jmax; j++)
      {
          ghost_cell = domain.cell(0,j).velocity(velocity_type::V);
          inner_cell = domain.cell(1,j).velocity(velocity_type::V);
          REQUIRE( fabs(ghost_cell + inner_cell) < eps);

          ghost_cell = domain.cell(0,j).velocity(velocity_type::U);
          REQUIRE( ghost_cell == 0.0);

          ghost_cell = domain.cell(imax+1,j).velocity(velocity_type::V);
          inner_cell = domain.cell(imax,j).velocity(velocity_type::V);
          REQUIRE( fabs(ghost_cell + inner_cell) < eps);

          ghost_cell = domain.cell(imax,j).velocity(velocity_type::U);
          REQUIRE( ghost_cell == 0.0);
      }
      #pragma endregion
    }

    SECTION("Plane_Shear_Flow")
    {
      imax = 10;
      jmax = 5;
      UT = 0.0;
      // Geometry initialization
      geometry = imatrix(0, imax+1, 0, jmax+1);
      // Setting upper and lower boundary
      for(int i = 0; i <= imax+1; i++)
      {
        geometry[i][0] = 0;
        geometry[i][jmax+1] = 0;
      }
      // Setting right and left boundary
      for(int j = 1; j <= jmax; j++)
      {
        geometry[0][j] = 3 ;
        geometry[imax+1][j] = 2;
      }
      // Setting inner Cell
      for(int i = 1; i <= imax; i++)
      { 
        for(int j = 1; j <= jmax; j++)
        {
          geometry[i][j] = 4;
        }
      } 
      Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry) ;

      boundaryvalues(imax,jmax,domain,UT);
      
      // Testing upper and lower boundary
      for(int i = 1; i <= imax; i++)
      { 
          // Check that obstacle cells have zero velocity at the boundary and negative value at ghost values
          ghost_cell = domain.cell(i,0).velocity(velocity_type::U);
          inner_cell = domain.cell(i,1).velocity(velocity_type::U);
          REQUIRE( fabs(ghost_cell + inner_cell) < eps);

          ghost_cell = domain.cell(i,0).velocity(velocity_type::V);
          REQUIRE( ghost_cell == 0.0);
          // Check that obstacle cells have zero velocity at the boundary and negative value at ghost values
          ghost_cell = domain.cell(i,jmax+1).velocity(velocity_type::U);
          inner_cell = domain.cell(i,jmax).velocity(velocity_type::U);
          REQUIRE( fabs(ghost_cell + inner_cell) < eps);

          ghost_cell = domain.cell(i,jmax).velocity(velocity_type::V);
          REQUIRE( ghost_cell == 0.0);
      }

      // Testing left and right boundary
      for(int j = 1; j <= jmax; j++)
      {
          // Check that Inflow cells have the predescribed values for inflow
          ghost_cell = domain.cell(0,j).velocity(velocity_type::V);
          REQUIRE( ghost_cell == 0.0);
  
          ghost_cell = domain.cell(0,j).velocity(velocity_type::U);
          REQUIRE( ghost_cell == 1.0);

          // Check that outflow Cells have the same values as the previous cells for both U and V
          ghost_cell = domain.cell(imax+1,j).velocity(velocity_type::V);
          inner_cell = domain.cell(imax,j).velocity(velocity_type::V);
          REQUIRE( fabs(ghost_cell - inner_cell) < eps);
  
          ghost_cell = domain.cell(imax,j).velocity(velocity_type::U);
          inner_cell = domain.cell(imax,j).velocity(velocity_type::U);
          REQUIRE( fabs(ghost_cell - inner_cell) < eps);
      }

      // Testing Corners
      // Left Bottom Corner
      ghost_cell = domain.cell(0,0).velocity(velocity_type::U);
      inner_cell = domain.cell(0,1).velocity(velocity_type::U);
      REQUIRE( fabs(ghost_cell + inner_cell) < eps);
      ghost_cell = domain.cell(0,0).velocity(velocity_type::V);
      REQUIRE( ghost_cell == 0.0);
      // Right Bottom Corner
      ghost_cell = domain.cell(imax+1,0).velocity(velocity_type::U);
      inner_cell = domain.cell(imax+1,1).velocity(velocity_type::U);
      REQUIRE( fabs(ghost_cell + inner_cell) < eps);
      ghost_cell = domain.cell(imax+1,0).velocity(velocity_type::V);
      REQUIRE( ghost_cell == 0.0);
      // Left Top Corner
      ghost_cell = domain.cell(0,jmax+1).velocity(velocity_type::U);
      inner_cell = domain.cell(0,jmax).velocity(velocity_type::U);
      REQUIRE( fabs(ghost_cell + inner_cell) < eps);
      ghost_cell = domain.cell(0,jmax+1).velocity(velocity_type::V);
      REQUIRE( ghost_cell == 0.0);
      // Right Top Corner
      ghost_cell = domain.cell(imax+1,jmax+1).velocity(velocity_type::U);
      inner_cell = domain.cell(imax+1,jmax).velocity(velocity_type::U);
      REQUIRE( fabs(ghost_cell + inner_cell) < eps);
      ghost_cell = domain.cell(0,0).velocity(velocity_type::V);
      REQUIRE( ghost_cell == 0.0);
    }
    
    // SECTION("Karman_Vortex")
    // {
    //   imax = 8;
    //   jmax = 8;
    //   UT = 0.0;
    //   // Geometry initialization
    //   geometry = read_pgm((std::string(Deb_Data_t) + "Karman_Vortex_Street.pgm").c_str() );
    //   Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry) ;

    //   boundaryvalues(imax,jmax,domain,UT);
    //   boundaryvalues(imax,jmax,domain,UT);
    //   domain.print_velocity(velocity_type::V);
    //   // Checking 8 Obstacle Cells
    //   // Checking NW boundary cells
      
    //   int i;
    //   int j;
    //   for (int counter = 3; counter <= 5; counter++)
    //   {
    //     i = counter;
    //     j = counter + 1;
    //     ghost_cell = domain.cell(i,j).velocity(velocity_type::V);
    //     inner_cell = domain.cell(i-1,j).velocity(velocity_type::U);
    //     REQUIRE( ghost_cell == 0);
    //     REQUIRE( inner_cell == 0);

    //     ghost_cell = domain.cell(i,j).velocity(velocity_type::U);
    //     inner_cell = domain.cell(i,j+1).velocity(velocity_type::U);
    //     REQUIRE( fabs(ghost_cell + inner_cell) < eps);
        
    //     ghost_cell = domain.cell(i,j-1).velocity(velocity_type::V);
    //     inner_cell = domain.cell(i-1,j-1).velocity(velocity_type::V);
    //     REQUIRE( fabs(ghost_cell + inner_cell) < eps);

    //     i = counter + 1; 
    //     j = counter;

    //     // SE
    //     ghost_cell = domain.cell(i,j).velocity(velocity_type::U);
    //     inner_cell = domain.cell(i,j-1).velocity(velocity_type::V);
    //     REQUIRE( ghost_cell == 0);
    //     REQUIRE( inner_cell == 0);
        
    //     ghost_cell = domain.cell(i-1,j).velocity(velocity_type::U);
    //     inner_cell = domain.cell(i-1,j-1).velocity(velocity_type::U);
    //     REQUIRE( fabs(ghost_cell + inner_cell) < eps);

    //     ghost_cell = domain.cell(i,j).velocity(velocity_type::V);
    //     inner_cell = domain.cell(i+1,j).velocity(velocity_type::V);
    //     REQUIRE( fabs(ghost_cell + inner_cell) < eps);
    //   }

    
  
    //   // Checking SW boundary cells
    //   #pragma region 
    //   ghost_cell = domain.cell(2,3).velocity(velocity_type::U);
    //   inner_cell = domain.cell(3,2).velocity(velocity_type::V);
    //   REQUIRE( ghost_cell == 0);
    //   REQUIRE( inner_cell == 0);

    //   ghost_cell = domain.cell(3,3).velocity(velocity_type::U);
    //   inner_cell = domain.cell(3,2).velocity(velocity_type::U);
    //   REQUIRE( fabs(ghost_cell + inner_cell) < eps);

    //   ghost_cell = domain.cell(3,3).velocity(velocity_type::V);
    //   inner_cell = domain.cell(2,3).velocity(velocity_type::V);
    //   REQUIRE( fabs(ghost_cell + inner_cell) < eps);
    //   #pragma endregion
    //   // Checking NE boundary cells
    //   #pragma region 
    //   ghost_cell = domain.cell(6,6).velocity(velocity_type::U);
    //   inner_cell = domain.cell(6,6).velocity(velocity_type::V);
    //   REQUIRE( ghost_cell == 0);
    //   REQUIRE( inner_cell == 0);

    //   ghost_cell = domain.cell(5,6).velocity(velocity_type::U);
    //   inner_cell = domain.cell(5,7).velocity(velocity_type::U);
    //   REQUIRE( fabs(ghost_cell + inner_cell) < eps);

    //   ghost_cell = domain.cell(6,5).velocity(velocity_type::V);
    //   inner_cell = domain.cell(7,5).velocity(velocity_type::V);
    //   REQUIRE( fabs(ghost_cell + inner_cell) < eps);
    //   #pragma endregion
    //   // Testing upper and lower boundary
    //   for(int i = 1; i <= imax; i++)
    //   { 
    //       // Check that obstacle cells have zero velocity at the boundary and negative value at ghost values
    //       ghost_cell = domain.cell(i,0).velocity(velocity_type::U);
    //       inner_cell = domain.cell(i,1).velocity(velocity_type::U);
    //       REQUIRE( fabs(ghost_cell + inner_cell) < eps);

    //       ghost_cell = domain.cell(i,0).velocity(velocity_type::V);
    //       REQUIRE( ghost_cell == 0.0);
    //       // Check that obstacle cells have zero velocity at the boundary and negative value at ghost values
    //       ghost_cell = domain.cell(i,jmax+1).velocity(velocity_type::U);
    //       inner_cell = domain.cell(i,jmax).velocity(velocity_type::U);
    //       REQUIRE( fabs(ghost_cell + inner_cell) < eps);

    //       ghost_cell = domain.cell(i,jmax).velocity(velocity_type::V);
    //       REQUIRE( ghost_cell == 0.0);
    //   }

    //   // Testing left and right boundary
    //   for(int j = 1; j <= jmax; j++)
    //   {
    //       // Check that Inflow cells have the predescribed values for inflow
    //       ghost_cell = domain.cell(0,j).velocity(velocity_type::V);
    //       REQUIRE( ghost_cell == 0.0);
  
    //       ghost_cell = domain.cell(0,j).velocity(velocity_type::U);
    //       REQUIRE( ghost_cell == 1.0);

    //       // Check that outflow Cells have the same values as the previous cells for both U and V
    //       ghost_cell = domain.cell(imax+1,j).velocity(velocity_type::V);
    //       inner_cell = domain.cell(imax,j).velocity(velocity_type::V);
    //       REQUIRE( fabs(ghost_cell - inner_cell) < eps);
  
    //       ghost_cell = domain.cell(imax,j).velocity(velocity_type::U);
    //       inner_cell = domain.cell(imax,j).velocity(velocity_type::U);
    //       REQUIRE( fabs(ghost_cell - inner_cell) < eps);
    //   }

    //   // Testing Corners
    //   // Left Bottom Corner
    //   ghost_cell = domain.cell(0,0).velocity(velocity_type::U);
    //   inner_cell = domain.cell(0,1).velocity(velocity_type::U);
    //   REQUIRE( fabs(ghost_cell + inner_cell) < eps);
    //   ghost_cell = domain.cell(0,0).velocity(velocity_type::V);
    //   REQUIRE( ghost_cell == 0.0);
    //   // Right Bottom Corner
    //   ghost_cell = domain.cell(imax+1,0).velocity(velocity_type::U);
    //   inner_cell = domain.cell(imax+1,1).velocity(velocity_type::U);
    //   REQUIRE( fabs(ghost_cell + inner_cell) < eps);
    //   ghost_cell = domain.cell(imax+1,0).velocity(velocity_type::V);
    //   REQUIRE( ghost_cell == 0.0);
    //   // Left Top Corner
    //   ghost_cell = domain.cell(0,jmax+1).velocity(velocity_type::U);
    //   inner_cell = domain.cell(0,jmax).velocity(velocity_type::U);
    //   REQUIRE( fabs(ghost_cell + inner_cell) < eps);
    //   ghost_cell = domain.cell(0,jmax+1).velocity(velocity_type::V);
    //   REQUIRE( ghost_cell == 0.0);
    //   // Right Top Corner
    //   ghost_cell = domain.cell(imax+1,jmax+1).velocity(velocity_type::U);
    //   inner_cell = domain.cell(imax+1,jmax).velocity(velocity_type::U);
    //   REQUIRE( fabs(ghost_cell + inner_cell) < eps);
    //   ghost_cell = domain.cell(0,0).velocity(velocity_type::V);
    //   REQUIRE( ghost_cell == 0.0);
    // }
}

TEST_CASE("Test calculate dt ", "[calculate_dt]")
{
    // Dummy Grid Parameters
    #pragma region 
    double PI;
    double UI;
    double VI;
    double TI =0.0;

    double Re = 100;

    int imax = 10;
    int jmax = 10;

    double dx = 0.1;
    double dy = 0.1;

    double tau = 0.5;
    double dt = 0.5;

    double PR = 1;
    bool temp_flag = false;
    
    int** geometry;
    #pragma endregion
    

    // Geometry initialization
      geometry = imatrix(0, imax+1, 0, jmax+1);

      // Setting upper and lower boundary
      for(int i = 0; i <= imax+1; i++)
      {
        geometry[i][0] = 0;
        geometry[i][jmax+1] = 0;
      }
      // Setting righ and left boundary
      for(int j = 0; j <= jmax+1; j++)
      {
        geometry[0][j] = 0 ;
        geometry[imax+1][j] = 0;
      }

      // Setting inner Cell
      for(int i = 1; i <= imax; i++)
      { 
        for(int j = 1; j <= jmax; j++)
        {
          geometry[i][j] = 4;
        }
      } 
    // Initialzie Dummy Grid
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry) ;
    // Data Structure Declaration
    #pragma region 
    matrix<double> U;
    matrix<double> V;
    matrix<double> P;
    #pragma endregion

    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();

     
    /** Dummy Scenario 1
     * dx/umax = inf
     * dy/vmax = inf
     * (Re/2)*(1/dx) + (1/dy) = 0.25
     * result tau * 0.25 = 0.125
     * */
    #pragma region 
    PI = 0.0;
    UI = 0.0;
    VI = 0.0;
    U.resize(imaxb,std::vector<double>(jmaxb,UI));
    V.resize(imaxb,std::vector<double>(jmaxb,VI));
    P.resize(imaxb,std::vector<double>(jmaxb,PI));

    domain.set_velocity(U,velocity_type::U);
    domain.set_velocity(V,velocity_type::V);
    domain.set_pressure(P);
    calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,PR,&temp_flag,domain);

    REQUIRE( dt == 0.125);
    
    P.clear();
    U.clear();
    V.clear();
    #pragma endregion

    /** Dummy Scenario 2
     * dx/umax = inf
     * dy/vmax = 0.025
     * (Re/2)*(1/dx) + (1/dy) = 0.25
     * result tau * 0.025 = 0.0125
     * */
    #pragma region 
    PI = 0.0;
    UI = 0.0;
    VI = 4.0;
    P.resize(imaxb,std::vector<double>(jmaxb,PI));
    U.resize(imaxb,std::vector<double>(jmaxb,UI));
    V.resize(imaxb,std::vector<double>(jmaxb,VI));
    
    domain.set_pressure(P);
    domain.set_velocity(U,velocity_type::U);
    domain.set_velocity(V,velocity_type::V);  

    calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,PR,&temp_flag,domain);

    REQUIRE( dt == 0.0125);

    P.clear();
    U.clear();
    V.clear();
    #pragma endregion

    /** Dummy Scenario 3
     * dx/umax = 0.05
     * dy/vmax = inf
     * (Re/2)*(1/dx) + (1/dy) = 0.25
     * result tau * 0.05 = 0.025
     * */
    #pragma region 
    PI = 0.0;
    UI = 2.0;
    VI = 0.0;
    U.resize(imaxb,std::vector<double>(jmaxb,UI));
    V.resize(imaxb,std::vector<double>(jmaxb,VI));
    P.resize(imaxb,std::vector<double>(jmaxb,PI));

    domain.set_velocity(U,velocity_type::U);
    domain.set_velocity(V,velocity_type::V);
    domain.set_pressure(P);
    calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,PR,&temp_flag,domain);

    REQUIRE( dt == 0.025);

    P.clear();
    U.clear();
    V.clear();
    #pragma endregion

    /* Real Scenario*/
    #pragma region 
    double dt_ref;

    // File and directory paths
    std::string logFileName;

    std::ifstream fin(ref_log_dt_t);
    fin >> dt_ref;

    U.resize(imaxb,std::vector<double>(jmaxb,UI));
    V.resize(imaxb,std::vector<double>(jmaxb,VI));
    
    // Reading log files in data structures
    #pragma region 
    logFileName = ref_log_U_t + "1";
    read_matrix(logFileName.c_str(),U,imaxb,jmaxb);
    logFileName = ref_log_V_t + "1";
    read_matrix(logFileName.c_str(),V,imaxb,jmaxb);
    #pragma endregion

    domain.set_velocity(U, velocity_type::U);
    domain.set_velocity(V, velocity_type::V);

    calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,PR,&temp_flag,domain);

    std::ofstream fout(cmp_log_dt_t);
    fout << std::setprecision(4) << dt;
    fout << std::endl;

    const std::string hashdt = get_md5hash(ref_log_dt_t);
    const std::string hash_cmpdt = get_md5hash(cmp_log_dt_t);

    REQUIRE (hashdt == hash_cmpdt);

    U.clear();
    V.clear();
    #pragma endregion
}

TEST_CASE( "Test F and G matrix", "[calculate_FG]" )
{
    SECTION("LID_DrivenCavity")
    {
      std::cout << "  " << std::endl;
    std::cout<< "================================================================="<<std::endl;
    std::cout<< "---------------------Testing calculate_fg()----------------------"<<std::endl;
    std::cout<< "================================================================="<<std::endl;
    std::cout << "  " << std::endl;
    std::string szFileName = "../../cavity100.dat";

    // Dummy Grid Paramters
    #pragma region 
    
    int imax = 10;
    int jmax = 10;

    double UI = 0.0;
    double VI = 0.0;
    double PI = 0.0;
    double TI = 0.0;

    double UT = 1.0;

    int** geometry;
    #pragma endregion
    
    // Geometry initialization
    geometry = imatrix(0, imax+1, 0, jmax+1);
    // Setting upper and lower boundary
    for(int i = 0; i <= imax+1; i++)
      {
        geometry[i][0] = 0;
        geometry[i][jmax+1] = 0;
      }
    // Setting righ and left boundary
    for(int j = 0; j <= jmax+1; j++)
      {
        geometry[0][j] = 0 ;
        geometry[imax+1][j] = 0;
      }
    // Setting inner Cell
    for(int i = 1; i <= imax; i++)
      { 
        for(int j = 1; j <= jmax; j++)
        {
          geometry[i][j] = 4;
        }
      }
    // Initialize Dummy Grid
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry);
    // Dummy Scenario Parameters
    #pragma region 
    double Re = 100;

    double dx = 0.1;
    double dy = 0.1;

    double GX = 0.0;
    double GY = 0.0;

    double alpha = 0.5;
    double dt = 0.0803324;

    #pragma endregion


    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();

    // File and directory paths
    std::string logFileName;

    // Data Structure Declaration
    #pragma region 
    // Data structures for first time step
    matrix<double> F1;
    matrix<double> G1;
    matrix<double> RS1;
    matrix<double> P1;
    matrix<double> U1;
    matrix<double> V1;

    // Data structures for second time step
    matrix<double> F2;
    matrix<double> G2;
    matrix<double> RS2;
    matrix<double> P2;
    matrix<double> U2;
    matrix<double> V2;
    #pragma endregion

    // Data Structure Allocation
    #pragma region 
    // Allocating memory for matrices
    F1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    RS1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    P1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    U1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V1.resize(imaxb,std::vector<double>(jmaxb,0.0));

    // Allocating memory for matrices
    F2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    RS2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    P2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    U2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    #pragma endregion

    // Reading log files in data structures
    #pragma region 
    logFileName = ref_log_F_t + "1";
    read_matrix(logFileName.c_str(),F1,imaxb,jmaxb);
    logFileName = ref_log_G_t + "1";
    read_matrix(logFileName.c_str(),G1,imaxb,jmaxb);
    logFileName = ref_log_RS_t + "1";
    read_matrix(logFileName.c_str(),RS1,imaxb,jmaxb);
    logFileName = ref_log_P_t + "1";
    read_matrix(logFileName.c_str(),P1,imaxb,jmaxb);
    logFileName = ref_log_U_t + "1";
    read_matrix(logFileName.c_str(),U1,imaxb,jmaxb);
    logFileName = ref_log_V_t + "1";
    read_matrix(logFileName.c_str(),V1,imaxb,jmaxb);
    #pragma endregion


    domain.set_velocity(U1, velocity_type::U);
    domain.set_velocity(V1, velocity_type::V);
    domain.set_pressure(P1);

    // std::cout << "-------------" << std::endl;
    logFileName = ref_log_F_t + "2";
    const std::string hashF2 = get_md5hash(logFileName);

    logFileName = ref_log_G_t + "2";
    const std::string hashG2 = get_md5hash(logFileName);
    
    // Dummy Scenario
    #pragma region
    

    boundaryvalues(imax, jmax, domain,UT);
    calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, domain, F2, G2);

    #pragma endregion
    


    write_matrix_less_precision(cmp_log_F_t, F2, imaxb, jmaxb);
    write_matrix_less_precision(cmp_log_G_t, G2, imaxb, jmaxb);

    const std::string hash_cmpF2 = get_md5hash(cmp_log_F_t);
    const std::string hash_cmpG2 = get_md5hash(cmp_log_G_t);

    #pragma region

    REQUIRE (hashF2 == hash_cmpF2);
    std::cout<<hashF2<<" == "<<hash_cmpF2<<std::endl;
    REQUIRE (hashG2 == hash_cmpG2);
    std::cout<<hashG2<<" == "<<hash_cmpG2<<std::endl;

    #pragma endregion
    }
    
    
    
      
}

TEST_CASE( "Test rhs of PPE matrix", "[calculate_rs]" )
{
    std::cout << "  " << std::endl;
    std::cout<< "================================================================="<<std::endl;
    std::cout<< "---------------------Testing calculate_rhs()----------------------"<<std::endl;
    std::cout<< "================================================================="<<std::endl;
    std::cout << "  " << std::endl;
    std::string szFileName = "../../cavity100.dat";


    // Dummy Grid Paramters
    #pragma region 
    
    int imax = 10;
    int jmax = 10;

    double UI = 0.0;
    double VI = 0.0;
    double PI = 0.0;
    double TI = 0.0;

    int** geometry;
    #pragma endregion
    
    // Geometry initialization
      geometry = imatrix(0, imax+1, 0, jmax+1);

      // Setting upper and lower boundary
      for(int i = 0; i <= imax+1; i++)
      {
        geometry[i][0] = 0;
        geometry[i][jmax+1] = 0;
      }
      // Setting righ and left boundary
      for(int j = 0; j <= jmax+1; j++)
      {
        geometry[0][j] = 0 ;
        geometry[imax+1][j] = 0;
      }

      // Setting inner Cell
      for(int i = 1; i <= imax; i++)
      { 
        for(int j = 1; j <= jmax; j++)
        {
          geometry[i][j] = 4;
        }
      }
    // Initialize Dummy Grid
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry);
    // Dummy Scenario Parameters
    #pragma region 
    double Re = 100;

    double dx = 0.1;
    double dy = 0.1;

    double GX = 0.0;
    double GY = 0.0;

    double alpha = 0.5;
    // double dt = 0.0803324;
    double dt = 0.08033235832; 
    #pragma endregion
    
    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();



    // File and directory paths
    std::string logFileName;

    // Data structure Declaration
    #pragma region 
    // Data structures for first time step
    matrix<double> F1;
    matrix<double> G1;
    matrix<double> RS1;
    matrix<double> P1;
    matrix<double> U1;
    matrix<double> V1;

    // Data structures for second time step
    matrix<double> F2;
    matrix<double> G2;
    matrix<double> RS2;
    matrix<double> P2;
    matrix<double> U2;
    matrix<double> V2;
    #pragma endregion

    // Data Structure initialization
    #pragma region 
    // Allocating memory for matrices
    F1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    RS1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    #pragma endregion

    // Reading log files in data structure
    #pragma region 
    logFileName = ref_log_F_t + "1";
    read_matrix(logFileName.c_str(),F1,imaxb,jmaxb);
    logFileName = ref_log_G_t + "1";
    read_matrix(logFileName.c_str(),G1,imaxb,jmaxb);
    logFileName = ref_log_RS_t + "1";
    #pragma endregion    


    logFileName = ref_log_RS_t + "1";
    const std::string hashRS1 = get_md5hash(logFileName);
    
    #pragma region
    calculate_rs(dt, dx, dy, imax, jmax, F1, G1, RS1,domain);

    write_matrix_less_precision(cmp_log_RS_t, RS1, imaxb, jmaxb);

    const std::string hash_cmpRS1 = get_md5hash(cmp_log_RS_t);

    REQUIRE (hashRS1 == hash_cmpRS1);
    std::cout<<hashRS1<<" == "<<hash_cmpRS1<<std::endl;
    #pragma endregion
    
}

TEST_CASE( "Test SOR", "[sor_iteration]" )
{
    std::cout << "  " << std::endl;
    std::cout<< "================================================================="<<std::endl;
    std::cout<< "---------------------------Testing sor()-------------------------"<<std::endl;
    std::cout<< "================================================================="<<std::endl;
    std::cout << "  " << std::endl;
    std::string szFileName = "../../cavity100.dat";


    // Dummy Grid Paramters
    #pragma region 
    
    int imax = 10;
    int jmax = 10;

    double UI = 0.0;
    double VI = 0.0;
    double PI = 0.0;
    double TI = 0.0;
    int** geometry;
    #pragma endregion
    
    // Geometry initialization
      geometry = imatrix(0, imax+1, 0, jmax+1);

      // Setting upper and lower boundary
      for(int i = 0; i <= imax+1; i++)
      {
        geometry[i][0] = 0;
        geometry[i][jmax+1] = 0;
      }
      // Setting righ and left boundary
      for(int j = 0; j <= jmax+1; j++)
      {
        geometry[0][j] = 0 ;
        geometry[imax+1][j] = 0;
      }

      // Setting inner Cell
      for(int i = 1; i <= imax; i++)
      { 
        for(int j = 1; j <= jmax; j++)
        {
          geometry[i][j] = 4;
        }
      }
    // Initialize Dummy Grid
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry);
    // Dummy Scenario Parameters
    #pragma region 
    double Re = 100;

    double dx = 0.1;
    double dy = 0.1;

    double GX = 0.0;
    double GY = 0.0;

    double alpha = 0.5;
    double omg = 1.7;
    double eps = 0.001;
    // double dt = 0.5;
    double dt = 0.08033235832; 

    int itermax = 100;

    #pragma endregion
    
    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();


    // File and directory paths
    std::string logFileName;

    // Data structure Declaration
    #pragma region 
    // Data structures for first time step
    matrix<double> F1;
    matrix<double> G1;
    matrix<double> RS1;
    matrix<double> P1;
    matrix<double> U1;
    matrix<double> V1;

    // Data structures for second time step
    matrix<double> F2;
    matrix<double> G2;
    matrix<double> RS2;
    matrix<double> P2;
    matrix<double> U2;
    matrix<double> V2;
    #pragma endregion

    // Data structure Initialization
    #pragma region 
    // Allocating memory for matrices
    F1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    RS1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    P1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    U1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V1.resize(imaxb,std::vector<double>(jmaxb,0.0));

    // Allocating memory for matrices
    F2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    RS2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    P2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    U2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    #pragma endregion
    

    logFileName = ref_log_RS_t + "1";
    read_matrix(logFileName.c_str(),RS1,imaxb,jmaxb);
    logFileName = ref_log_P_t + "1";
    read_matrix(logFileName.c_str(),P1,imaxb,jmaxb);
    logFileName = ref_log_U_t + "1";
    read_matrix(logFileName.c_str(),U1,imaxb,jmaxb);
    logFileName = ref_log_V_t + "1";
    read_matrix(logFileName.c_str(),V1,imaxb,jmaxb);

    logFileName = ref_log_RS_t + "2";
    read_matrix(logFileName.c_str(),RS2,imaxb,jmaxb);
    logFileName = ref_log_P_t + "2";
    const std::string hashP2 = get_md5hash(logFileName);

    
    domain.set_velocity(U1, velocity_type::U);
    domain.set_velocity(V1, velocity_type::V);
    domain.set_pressure(P1);

    // boundaryvalues(imax, jmax, domain);
    double res = 1.0;
    int it = 0;
    while(it < itermax && res > eps){
      sor(omg,dx,dy,imax,jmax,domain,RS2,&res);
      it++;
    }

    domain.pressure(P2);
    write_matrix_less_precision(cmp_log_P_t, P2, imaxb, jmaxb);

    const std::string hash_cmpP2 = get_md5hash(cmp_log_P_t);

    REQUIRE (hashP2 == hash_cmpP2);
    std::cout<<hashP2<<" == "<<hash_cmpP2<<std::endl;    

}

TEST_CASE("Test velocity (U and v) values", "[calculate_uv]")
{
    std::cout << "  " << std::endl;
    std::cout<< "================================================================="<<std::endl;
    std::cout<< "---------------------Testing calculate_uv()----------------------"<<std::endl;
    std::cout<< "================================================================="<<std::endl;
    std::cout << "  " << std::endl;
    std::string szFileName = "../../cavity100.dat";

    // Dummy Grid Paramters
    #pragma region 
    
    int imax = 10;
    int jmax = 10;

    double UI = 0.0;
    double VI = 0.0;
    double PI = 0.0;
    double TI = 0.0;

    double UT = 1.0;

    int** geometry;
    #pragma endregion

    // Geometry initialization
    geometry = imatrix(0, imax+1, 0, jmax+1);
    // Setting upper and lower boundary
    for(int i = 0; i <= imax+1; i++)
      {
        geometry[i][0] = 0;
        geometry[i][jmax+1] = 0;
      }
    // Setting righ and left boundary
    for(int j = 0; j <= jmax+1; j++)
      {
        geometry[0][j] = 0 ;
        geometry[imax+1][j] = 0;
      }
    // Setting inner Cell
    for(int i = 1; i <= imax; i++)
      { 
        for(int j = 1; j <= jmax; j++)
        {
          geometry[i][j] = 4;
        }
      }
    // Initialize Dummy Grid
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry);
    // Dummy Scenario Parameters
    #pragma region 
    double Re = 100;

    double dx = 0.1;
    double dy = 0.1;

    double GX = 0.0;
    double GY = 0.0;

    double alpha = 0.5;
    double omg = 1.7;
    double eps = 0.001;
    double dt = 0.08033235832; 
    
    int itermax = 100;
    #pragma endregion 

    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();

    // File and directory paths
    std::string logFileName;

    // Data Structure Declaration
    #pragma region 
    // Data structures for first time step
    matrix<double> U1;
    matrix<double> V1;

    // Data structures for first time step
    matrix<double> F2;
    matrix<double> G2;
    matrix<double> P2;
    matrix<double> U2;
    matrix<double> V2;
    #pragma endregion

    // Data Structure Initialization
    #pragma region 
    // Allocating memory for matrices
    U1.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V1.resize(imaxb,std::vector<double>(jmaxb,0.0));

    // Allocating memory for matrices
    F2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    P2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    U2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V2.resize(imaxb,std::vector<double>(jmaxb,0.0));
    #pragma endregion

    // Reading log files in data structures
    #pragma region 
    logFileName = ref_log_U_t + "1";
    read_matrix(logFileName.c_str(),U1,imaxb,jmaxb);
    logFileName = ref_log_V_t + "1";
    read_matrix(logFileName.c_str(),V1,imaxb,jmaxb);
    logFileName = ref_log_F_t + "22";
    read_matrix(logFileName.c_str(),F2,imaxb,jmaxb);
    logFileName = ref_log_G_t + "22";
    read_matrix(logFileName.c_str(),G2,imaxb,jmaxb);
    logFileName = ref_log_P_t + "22";
    read_matrix(logFileName.c_str(),P2,imaxb,jmaxb);
    #pragma endregion

    // Calculating ref files hasing
    #pragma region 
    logFileName = ref_log_U_t + "2";
    const std::string hashU2 = get_md5hash(logFileName);
    logFileName = ref_log_V_t + "2";
    const std::string hashV2 = get_md5hash(logFileName);
    #pragma endregion

    // Setting the Grid values to reference values
    #pragma region 
    domain.set_velocity(U1, velocity_type::U);
    domain.set_velocity(V1, velocity_type::V);
    domain.set_pressure(P2);
    #pragma endregion

    boundaryvalues(imax, jmax, domain,UT);
    calculate_uv(dt, dx, dy, imax, jmax, domain, F2, G2);

    domain.velocity(U2, velocity_type::U);
    domain.velocity(V2, velocity_type::V);


    write_matrix_less_precision(cmp_log_U_t, U2, imaxb, jmaxb);
    write_matrix_less_precision(cmp_log_V_t, V2, imaxb, jmaxb);


    const std::string hash_cmpU2 = get_md5hash(cmp_log_U_t);
    const std::string hash_cmpV2 = get_md5hash(cmp_log_V_t);

    REQUIRE (hashU2 == hash_cmpU2);
    std::cout<<hashU2<<" == "<<hash_cmpU2<<std::endl;
    REQUIRE (hashV2 == hash_cmpV2);
    std::cout<<hashV2<<" == "<<hash_cmpV2<<std::endl;


}

#pragma endregion


/*-------Boundary Type flag check*/
#pragma region 
/** Possible Cases For Cell Flag
       * F        f = 1
       * NS/BN    f = 34
       * NS/BS    f = 66
       * NS/BW    f = 130
       * NS/BE    f = 258
       * NS/BN/BW f = 162 
       * NS/BN/BE f = 290
       * NS/BS/BW f = 194
       * NS/BS/BE f = 322
       * FS/BN    f = 36
       * FS/BS    f = 68
       * FS/BW    f = 132
       * FS/BE    f = 260
       * FS/BN/BW f = 164
       * FS/BN/BE f = 292
       * FS/BS/BW f = 196
       * FS/BS/BE f = 324
       * OF/BN    f = 40
       * OF/BS    f = 72
       * OF/BW    f = 136
       * OF/BE    f = 264
       * OF/BN/BW f = 168
       * OF/BN/BE f = 296
       * OF/BS/BW f = 200
       * OF/BS/BE f = 328
       * IF/BN    f = 48
       * IF/BS    f = 80
       * IF/BW    f = 144
       * IF/BE    f = 272
       * IF/BN/BW f = 176
       * IF/BN/BE f = 304
       * IF/BS/BW f = 208
       * IF/BS/BE f = 336
       * 
       */

TEST_CASE("Test North Boundary", "[B_N] [!mayfail]")
{
    // NS/BN
    REQUIRE(B_N(static_cast<int>(boundary_type::NS_BN)) == true);
    // FS/BN
    REQUIRE(B_N(static_cast<int>(boundary_type::FS_BN)) == true);
    // OF/BN
    REQUIRE(B_N(static_cast<int>(boundary_type::OF_BN)) == true);
    // IF/BN
    REQUIRE(B_N(static_cast<int>(boundary_type::IF_BN)) == true);

}

TEST_CASE("Test South Boundary", "[B_S] [!mayfail]")
{
    // NS/BS
    REQUIRE(B_S(static_cast<int>(boundary_type::NS_BS)) == true);
    // FS/BS
    REQUIRE(B_S(static_cast<int>(boundary_type::FS_BS)) == true);
    // OF/BS
    REQUIRE(B_S(static_cast<int>(boundary_type::OF_BS)) == true);
    // IF/BS
    REQUIRE(B_S(static_cast<int>(boundary_type::IF_BS)) == true);

}

TEST_CASE("Test East Boundary", "[B_E] [!mayfail]")
{
    // NS/BE
    REQUIRE(B_E(static_cast<int>(boundary_type::NS_BE)) == true);
    // FS/BE
    REQUIRE(B_E(static_cast<int>(boundary_type::FS_BE)) == true);
    // OF/BE
    REQUIRE(B_E(static_cast<int>(boundary_type::OF_BE)) == true);
    // IF/BE
    REQUIRE(B_E(static_cast<int>(boundary_type::IF_BE)) == true);
    
}

TEST_CASE("Test West Boundary", "[B_W] [!mayfail]")
{
    // NS/BW
    REQUIRE(B_W(static_cast<int>(boundary_type::NS_BW)) == true);
    // FS/BW
    REQUIRE(B_W(static_cast<int>(boundary_type::FS_BW)) == true);
    // OF/BW
    REQUIRE(B_W(static_cast<int>(boundary_type::OF_BW)) == true);
    // IF/BW
    REQUIRE(B_W(static_cast<int>(boundary_type::IF_BW)) == true);
    
}

TEST_CASE("Test North West Boundary", "[B_NW] [!mayfail]")
{
    // NS/BN_BW
    REQUIRE(B_NW(static_cast<int>(boundary_type::NS_BN_BW)) == true);
    // FS/BN_BW
    REQUIRE(B_NW(static_cast<int>(boundary_type::FS_BN_BW)) == true);
    // OF/BN_BW
    REQUIRE(B_NW(static_cast<int>(boundary_type::OF_BN_BW)) == true);
    // IF/BN_BW
    REQUIRE(B_NW(static_cast<int>(boundary_type::IF_BN_BW)) == true);
    
}

TEST_CASE("Test North East Boundary", "[B_NE] [!mayfail]")
{
    // NS/BN
    REQUIRE(B_NE(static_cast<int>(boundary_type::NS_BN_BE)) == true);
    // FS/BN_BW
    REQUIRE(B_NE(static_cast<int>(boundary_type::FS_BN_BE)) == true);
    // OF/BN_BW
    REQUIRE(B_NE(static_cast<int>(boundary_type::OF_BN_BE)) == true);
    // IF/BN_BW
    REQUIRE(B_NE(static_cast<int>(boundary_type::IF_BN_BE)) == true);
    
}

TEST_CASE("Test South West Boundary", "[B_SW] [!mayfail]")
{
    // NS/BS_BW
    REQUIRE(B_SW(static_cast<int>(boundary_type::NS_BS_BW)) == true);
    // FS/BS_BW
    REQUIRE(B_SW(static_cast<int>(boundary_type::FS_BS_BW)) == true);
    // OF/BS_BW
    REQUIRE(B_SW(static_cast<int>(boundary_type::OF_BS_BW)) == true);
    // IF/BS_BW
    REQUIRE(B_SW(static_cast<int>(boundary_type::IF_BS_BW)) == true);
    
}

TEST_CASE("Test South East Boundary", "[B_SE] [!mayfail]")
{
    // NS/BS_BE
    REQUIRE(B_SE(static_cast<int>(boundary_type::NS_BS_BE)) == true);
    // FS/BS_BE
    REQUIRE(B_SE(static_cast<int>(boundary_type::FS_BS_BE)) == true);
    // OF/BS_BE
    REQUIRE(B_SE(static_cast<int>(boundary_type::OF_BS_BE)) == true);
    // IF/BS_BE
    REQUIRE(B_SE(static_cast<int>(boundary_type::IF_BS_BE)) == true);
    
}
#pragma endregion


/*----------------------unit Validation Tests using tolerance metric for element */

// Karman Vortex Validation Tests
#pragma region 

TEST_CASE("Test calculate dt Karman Vortex", "[calculate_dt_KV]")
{
    // Dummy Grid Parameters
    #pragma region 
    double xlength = 10;
    double ylength = 2;

    int imax = 98;
    int jmax = 18;

    double dt = 0.05;
    double dt_ref = 0.0;
    double t_end = 20.0;
    double tau = 0.5;

    double dt_Value = 2.0;

    int itermax = 500;
    double eps = 0.001;
    double omg = 1.7;
    double alpha = 0.9;

    double Re = 10000;
    double PR = 1;
    double GX = 0.0;
    double GY = 0.0;

    double UI = 1.0;
    double VI = 0.0;
    double PI = 0.0;
    double TI = 0.0;

    double UT = 0.0;

    double dx = 0.1;
    double dy = 0.1;

    double tol = 0.00000001;
    int** geometry;
    #pragma endregion
    
    //Getting log file names
    #pragma region
    // command line arrays
    int argc = 3;
    const char* argv[] = { "./sim","-k","-t"};
    // File paths string 
    std::string szFileName;
    std::string geoFileName;
    std::string problemName;
    std::vector<std::string> logFilesName;
    // Problem Flags
    bool temp_flag;

    logFilesName  = ParseTestingProblemFlags(argc,argv,szFileName,geoFileName,problemName,temp_flag);
    #pragma endregion
    // Geometry Initialization
    #pragma region 
    geometry = read_pgm((geoFileName).c_str());
    #pragma endregion

    // Initialzie Dummy Grid
    #pragma region    
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry) ;
    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();
    #pragma endregion

    // Data Structure Declaration
    #pragma region 
    matrix<double> U;
    matrix<double> V;
    matrix<double> P; 
    #pragma endregion

  
    // Data Structure Allocation
    #pragma region 
    // Input Data Structures
    U.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V.resize(imaxb,std::vector<double>(jmaxb,0.0));  
    #pragma endregion
     
    // Read Input matrices 
    #pragma region 
    // Reading dt
    std::ifstream fin(logFilesName[static_cast<int>(log_file::log_dt)]);
    fin >> dt_ref;

    // Input Data Structures  
    read_matrix(logFilesName[static_cast<int>(log_file::log_U_Input)],U,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_V_Input)],V,imaxb,jmaxb);

    #pragma endregion
   

    /* Real Scenario*/
    #pragma region 
  
    domain.set_velocity(U, velocity_type::U);
    domain.set_velocity(V, velocity_type::V);
    
    calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,PR,&temp_flag,domain);  
    #pragma endregion

    #pragma region 
    REQUIRE(fabs(dt - dt_ref) < tol );
    #pragma endregion
}

TEST_CASE( "Test F and G matrix Karman Vortex", "[calculate_fg_KV]" )
{
  // Dummy Grid Paramters
  #pragma region 
    
    double xlength = 10;
    double ylength = 2;

    int imax = 98;
    int jmax = 18;

    double dt = 0.05 ;
    double t_end = 20.0;
    double tau = 0.5;

    double dt_Value = 2.0;

    int itermax = 500;
    double eps = 0.001;
    double omg = 1.7;
    double alpha = 0.9;

    double Re = 10000;
    double GX = 0.0;
    double GY = 0.0;

    double UI = 1.0;
    double VI = 0.0;
    double PI = 0.0;
    double TI = 0.0;

    double UT = 0.0;

    double dx = 0.1;
    double dy = 0.1;

    double tol = 0.00000001;
    int** geometry;
    #pragma endregion
  
  //Getting log file names
  #pragma region
  // command line arrays
  int argc = 3;
  const char* argv[] = { "./sim","-k","-t"};
  // File paths string 
  std::string szFileName;
  std::string geoFileName;
  std::string problemName;
  std::vector<std::string> logFilesName;
  // Problem Flags
  bool temp_flag;

  logFilesName  = ParseTestingProblemFlags(argc,argv,szFileName,geoFileName,problemName,temp_flag);
  #pragma endregion

  // Geometry Initialization
  #pragma region 
  geometry = read_pgm((geoFileName).c_str());
  #pragma endregion
  
  // Initialize Dummy Grid
  #pragma region 
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry);

    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();
    #pragma endregion
  
  // Data Structure Declaration
  #pragma region 
    // Input Data Structures
    matrix<double> U;
    matrix<double> V;
    // Calculation Data Structures
    matrix<double> F_cal;
    matrix<double> G_cal;
    // Reference Data Structures
    matrix<double> F_ref;
    matrix<double> G_ref;
    #pragma endregion
  
  // Data Structure Allocation
  #pragma region 
    // Input Data Structures
    U.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    F_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    F_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    #pragma endregion
  
  // Reading log files in data structures
  #pragma region 
  // Reading dt
  std::ifstream fin(logFilesName[static_cast<int>(log_file::log_dt)]);
  fin >> dt;
  // Input Data Structures  
  read_matrix(logFilesName[static_cast<int>(log_file::log_U_Input)],U,imaxb,jmaxb);
  read_matrix(logFilesName[static_cast<int>(log_file::log_V_Input)],V,imaxb,jmaxb);
  // Reference Data Structures
  read_matrix(logFilesName[static_cast<int>(log_file::log_F)],F_ref,imaxb,jmaxb);
  read_matrix(logFilesName[static_cast<int>(log_file::log_G)],G_ref,imaxb,jmaxb);
  #pragma endregion
  
  // Dummy Scenario
  #pragma region
    
    domain.set_velocity(U, velocity_type::U);
    domain.set_velocity(V, velocity_type::V);
    
    calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, domain, F_cal, G_cal);

    #pragma endregion
  
  // Writing calculated matrices to log files
  #pragma region 
    write_matrix(cmp_log_F_t, F_cal, imaxb, jmaxb);
    write_matrix(cmp_log_G_t, G_cal, imaxb, jmaxb);
    #pragma endregion
  
  // Comparing Reference matrices to Calculated matrices
  #pragma region 

    for(int i = 0; i < imaxb; i++)
    {
      for(int j = 0; j < jmaxb; j++)
      {
        INFO("i = " << i << " j = " << j);
        REQUIRE(fabs(F_cal[i][j] - F_ref[i][j]) < tol);
        REQUIRE(fabs(G_cal[i][j] - G_ref[i][j]) < tol);
      }
    }

    #pragma endregion
  

}

TEST_CASE( "Test rhs of PPE matrix Karman Vortex", "[calculate_RS_KV]" )
{
  // Dummy Grid Parameters
    #pragma region 
    double xlength = 10;
    double ylength = 2;

    int imax = 98;
    int jmax = 18;

    double dt = 0.05 ;
    double t_end = 20.0;
    double tau = 0.5;

    double dt_Value = 2.0;

    int itermax = 500;
    double eps = 0.001;
    double omg = 1.7;
    double alpha = 0.9;

    double Re = 10000;
    double PR = 1;
    double GX = 0.0;
    double GY = 0.0;

    double UI = 1.0;
    double VI = 0.0;
    double PI = 0.0;
    double TI = 0.0;

    double UT = 0.0;

    double dx = 0.1;
    double dy = 0.1;

    double tol = 0.00000001;
    int** geometry;
    #pragma endregion
    
    //Getting log file names
    #pragma region
    // command line arrays
    int argc = 3;
    const char* argv[] = { "./sim","-k","-t"};
    // File paths string 
    std::string szFileName;
    std::string geoFileName;
    std::string problemName;
    std::vector<std::string> logFilesName;
    // Problem Flags
    bool temp_flag;
  
    logFilesName  = ParseTestingProblemFlags(argc,argv,szFileName,geoFileName,problemName,temp_flag);
    #pragma endregion

    // Geometry Initialization
    #pragma region 
    geometry = read_pgm((geoFileName).c_str());
    #pragma endregion

    // Initialzie Dummy Grid
    #pragma region    
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry) ;
    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();
    #pragma endregion

    
    
    // Data Structure Declaration
    #pragma region 
    // Input Data Structures
    matrix<double> F;
    matrix<double> G;
    // Calculation Data Structures
    matrix<double> RS_cal;
    // Reference Data Structures
    matrix<double> RS_ref;
    #pragma endregion

    // Data Structure Allocation
    #pragma region 
    // Input Data Structures
    F.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    RS_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    RS_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    #pragma endregion

    // Reading log files in data structures
    #pragma region 
    // Reading dt
    std::ifstream fin(logFilesName[static_cast<int>(log_file::log_dt)]);
    fin >> dt;    
    // Input Data Structures
    read_matrix(logFilesName[static_cast<int>(log_file::log_F)],F,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_G)],G,imaxb,jmaxb);

    // Reference Data Structures
    read_matrix(logFilesName[static_cast<int>(log_file::log_RS)],RS_ref,imaxb,jmaxb);
    #pragma endregion

    // Dummy Scenario
    #pragma region
    
    
    
    calculate_rs(dt, dx, dy, imax, jmax, F, G, RS_cal,domain);

    #pragma endregion
    // Writing calculated matrices to log files
    #pragma region 
    write_matrix(cmp_log_RS_t, RS_cal, imaxb, jmaxb);
    
    #pragma endregion
    // Comparing Reference matrices to Calculated matrices
    #pragma region 
    for(int i = 0; i < imaxb; i++)
    {
      for(int j = 0; j < jmaxb; j++)
      {
        if(RS_ref[i][j] == 0 )
        {
          // std::cout << "i = " << i << " j = " << j << std::endl;
            INFO("i = " << i << " j = " << j);
            REQUIRE(RS_cal[i][j] == RS_ref[i][j]);
        }
      }
    }

    for(int i = 0; i < imaxb; i++)
    {
      for(int j = 0; j < jmaxb; j++)
      {
        INFO("i = " << i << " j = " << j);
        REQUIRE(fabs(RS_cal[i][j] - RS_ref[i][j]) < tol);
        
      }
    }

    #pragma endregion
}

TEST_CASE( "Test Test SOR Karman Vortex", "[sor_KV]" )
{
  // Dummy Grid Paramters
    #pragma region 
    
    double xlength = 10;
    double ylength = 2;

    int imax = 98;
    int jmax = 18;

    double dt = 0.05;
    double t_end = 20.0;
    double tau = 0.5;

    double dt_Value = 2.0;

    int itermax = 500;
    double eps = 0.001;
    double omg = 1.7;
    double alpha = 0.9;

    double Re = 10000;
    double GX = 0.0;
    double GY = 0.0;

    double UI = 1.0;
    double VI = 0.0;
    double PI = 0.0;
    double TI = 0.0;

    double UT = 0.0;

    double dx = 0.1;
    double dy = 0.1;

    double tol = 0.00000001;
    int** geometry;
    #pragma endregion

    //Getting log file names
    #pragma region
    // command line arrays
    int argc = 3;
    const char* argv[] = { "./sim","-k","-t"};
    // File paths string 
    std::string szFileName;
    std::string geoFileName;
    std::string problemName;
    std::vector<std::string> logFilesName;
    // Problem Flags
    bool temp_flag;
  
    logFilesName  = ParseTestingProblemFlags(argc,argv,szFileName,geoFileName,problemName,temp_flag);
    #pragma endregion

    // Geometry Initialization
    #pragma region 
    geometry = read_pgm((geoFileName).c_str());
    #pragma endregion

    // Initialize Dummy Grid
    #pragma region 
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry);

    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();
    #pragma endregion

    
    
    // Data Structure Declaration
    #pragma region 
    // Input Data Structures
    matrix<double> U;
    matrix<double> V;
    matrix<double> RS;
    // Calculation Data Structures
    matrix<double> P_cal;
    // Reference Data Structures
    matrix<double> P_ref;
    #pragma endregion

    // Data Structure Allocation
    #pragma region 
    // Input Data Structures
    U.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V.resize(imaxb,std::vector<double>(jmaxb,0.0));
    RS.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    P_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    P_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    #pragma endregion

    // Reading log files in data structures
    #pragma region
    // Reading dt
    std::ifstream fin(logFilesName[static_cast<int>(log_file::log_dt)]);
    fin >> dt; 
    //Input Data    
    read_matrix(logFilesName[static_cast<int>(log_file::log_U_Input)],U,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_V_Input)],V,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_RS)],RS,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_P_Input)],P_cal,imaxb,jmaxb);

    // Reference Data
    read_matrix(logFilesName[static_cast<int>(log_file::log_P_Output)],P_ref,imaxb,jmaxb);

    #pragma endregion

    // Dummy Scenario
    #pragma region
    domain.set_velocity(U,velocity_type::U);
    domain.set_velocity(V,velocity_type::V);
    domain.set_pressure(P_cal);
    int it = 0;
    double res = 1.0;
    // Solve system using SOR
    while(it < itermax && res > eps){
      sor(omg,dx,dy,imax,jmax,domain,RS,&res);
      it++;
    }
    

    #pragma endregion
   
    // Writing calculated matrices to log files
    #pragma region 
    domain.pressure(P_cal);
    write_matrix(cmp_log_P_t, P_cal, imaxb, jmaxb);
    
    #pragma endregion
   
    // Comparing Reference matrices to Calculated matrices
    #pragma region 
   
    for(int i = 0; i < imaxb; i++)
    {
      for(int j = 0; j < jmaxb; j++)
      {
        INFO("i = " << i << " j = " << j);
        REQUIRE(fabs(P_cal[i][j] - P_ref[i][j]) < tol);
        
      }
    }

    #pragma endregion
}

TEST_CASE( "Test velocity (U and v) values  Karman Vortex", "[calculate_uv_KV]" )
{
  // Dummy Grid Paramters
    #pragma region 
    
    double xlength = 10;
    double ylength = 2;

    int imax = 98;
    int jmax = 18;

    double dt = 0.05 ;
    double t_end = 20.0;
    double tau = 0.5;

    double dt_Value = 2.0;

    int itermax = 500;
    double eps = 0.001;
    double omg = 1.7;
    double alpha = 0.9;

    double Re = 10000;
    double GX = 0.0;
    double GY = 0.0;

    double UI = 1.0;
    double VI = 0.0;
    double PI = 0.0;
    double TI = 0.0;

    double UT = 0.0;

    double dx = 0.1;
    double dy = 0.1;

    double tol = 0.00000001;
    int** geometry;
    #pragma endregion

  //Getting log file names
    #pragma region
    // command line arrays
    int argc = 3;
    const char* argv[] = { "./sim","-k","-t"};
    // File paths string 
    std::string szFileName;
    std::string geoFileName;
    std::string problemName;
    std::vector<std::string> logFilesName;
    // Problem Flags
    bool temp_flag;
  
    logFilesName  = ParseTestingProblemFlags(argc,argv,szFileName,geoFileName,problemName,temp_flag);
    #pragma endregion

    // Geometry Initialization
    #pragma region 
    geometry = read_pgm((geoFileName).c_str());
    #pragma endregion

    // Initialize Dummy Grid
    #pragma region 
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry);

    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();
    #pragma endregion

  // Data Structure Declaration
    #pragma region 
    // Input Data Structures
    matrix<double> F;
    matrix<double> G;
    matrix<double> P;
    // Calculation Data Structures
    matrix<double> U_cal;
    matrix<double> V_cal;
    // Reference Data Structures
    matrix<double> U_ref;
    matrix<double> V_ref;
    #pragma endregion

    // Data Structure Allocation
    #pragma region 
    // Input Data Structures
    F.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G.resize(imaxb,std::vector<double>(jmaxb,0.0));
    P.resize(imaxb,std::vector<double>(jmaxb,0.0));    
    // Calculation Data Structures
    U_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    U_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    #pragma endregion

    // Reading log files in data structures
    #pragma region 
    // Reading dt
    std::ifstream fin(logFilesName[static_cast<int>(log_file::log_dt)]);
    fin >> dt;
    // Input Data
    read_matrix(logFilesName[static_cast<int>(log_file::log_F)],F,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_G)],G,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_P_Output)],P,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_U_Input)],U_cal,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_V_Input)],V_cal,imaxb,jmaxb);

    // Reference Data Structures
    read_matrix(logFilesName[static_cast<int>(log_file::log_U_Output)],U_ref,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_V_Output)],V_ref,imaxb,jmaxb);
 
    #pragma endregion

    // Dummy Scenario
    #pragma region
    domain.set_velocity(U_cal,velocity_type::U);
    domain.set_velocity(V_cal,velocity_type::V);
    domain.set_pressure(P);
    calculate_uv(dt, dx, dy, imax, jmax, domain, F, G);

    #pragma endregion
    
    // Writing calculated matrices to log files
    #pragma region 
    domain.velocity(U_cal,velocity_type::U);
    domain.velocity(V_cal,velocity_type::V);
    write_matrix(cmp_log_U_t, U_cal , imaxb, jmaxb);
    write_matrix(cmp_log_V_t, V_cal , imaxb, jmaxb);
    #pragma endregion
    
    // Comparing Reference matrices to Calculated matrices
    #pragma region 
    // for(int i = 0; i < imaxb; i++)
    // {
    //   for(int j = 0; j < jmaxb; j++)
    //   {
    //     if(U_ref[i][j] == 0 )
    //     {
    //       // std::cout << "i = " << i << " j = " << j << std::endl;
    //         INFO("i = " << i << " j = " << j);
    //         REQUIRE(U_cal[i][j] == U_ref[i][j]);
    //     }
    //     if(V_ref[i][j] == 0 )
    //     {
    //       // std::cout << "i = " << i << " j = " << j << std::endl;
    //         INFO("i = " << i << " j = " << j);
    //         REQUIRE(V_cal[i][j] == V_ref[i][j]);
    //     }
        
    //   }
    // }

    for(int i = 1; i < imaxb-1; i++)
    {
      for(int j = 1; j < jmaxb-1; j++)
      {
        INFO("i = " << i << " j = " << j);
        REQUIRE(fabs(U_cal[i][j] - U_ref[i][j]) < tol);
        REQUIRE(fabs(V_cal[i][j] - V_ref[i][j]) < tol);
      }
    }

    #pragma endregion
}

#pragma endregion

// Flow Over A Step
#pragma region
TEST_CASE("Test calculate dt Flow Over A Step", "[calculate_dt_FS]")
{
    // Dummy Grid Paramters
    #pragma region 
    
    double xlength = 10;
    double ylength = 2;

    int imax = 98;
    int jmax = 18;

    double dt = 0.05;
    double dt_ref = 0.0;
    double t_end = 500.0;
    double tau = 0.5;

    double dt_Value = 10.0;

    int itermax = 500;
    double eps = 0.001;
    double omg = 1.7;
    double alpha = 0.9;

    double Re = 100;
    double PR = 1;
    double GX = 0.0;
    double GY = 0.0;

    double UI = 0.0;
    double VI = 0.0;
    double PI = 0.0;
    double TI = 0.0;

    double UT = 0.0;

    double dx = 0.1;
    double dy = 0.1;

    double tol = 0.00000001;
    int** geometry;
    #pragma endregion

    //Getting log file names
    #pragma region
    // command line arrays
    int argc = 3;
    const char* argv[] = { "./sim","-f","-t"};
    // File paths string 
    std::string szFileName;
    std::string geoFileName;
    std::string problemName;
    std::vector<std::string> logFilesName;
    // Problem Flags
    bool temp_flag;

    logFilesName  = ParseTestingProblemFlags(argc,argv,szFileName,geoFileName,problemName,temp_flag);
    #pragma endregion

    // Geometry Initialization
    #pragma region 
    geometry = read_pgm((geoFileName).c_str());
    #pragma endregion

    // Initialize Dummy Grid
    #pragma region 
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry);

    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();
    #pragma endregion

    // Data Structure Declaration
    #pragma region 
    matrix<double> U;
    matrix<double> V;
    matrix<double> P; 
    #pragma endregion

  
    // Data Structure Allocation
    #pragma region 
    // Input Data Structures
    U.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V.resize(imaxb,std::vector<double>(jmaxb,0.0));  
    #pragma endregion
     
    // Read Input matrices 
    #pragma region    
     // Reading dt
    std::ifstream fin(logFilesName[static_cast<int>(log_file::log_dt)]);
    fin >> dt_ref;

    // Input Data Structures  
    read_matrix(logFilesName[static_cast<int>(log_file::log_U_Input)],U,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_V_Input)],V,imaxb,jmaxb);

    #pragma endregion
   

    /* Real Scenario*/
    #pragma region 
  
    domain.set_velocity(U, velocity_type::U);
    domain.set_velocity(V, velocity_type::V);
    
    
    

    calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,PR,&temp_flag,domain);

    
    #pragma endregion

    #pragma region 
    REQUIRE(fabs(dt - dt_ref) < tol );
    #pragma endregion
}

TEST_CASE( "Test F and G matrix Flow Over A Step", "[calculate_fg_FS]" )
{
  // Dummy Grid Paramters
    #pragma region 
    
    double xlength = 10;
    double ylength = 2;

    int imax = 98;
    int jmax = 18;

    double dt = 0.05  ;
    double t_end = 500.0;
    double tau = 0.5;

    double dt_Value = 10.0;

    int itermax = 500;
    double eps = 0.001;
    double omg = 1.7;
    double alpha = 0.9;

    double Re = 100;
    double GX = 0.0;
    double GY = 0.0;

    double UI = 0.0;
    double VI = 0.0;
    double PI = 0.0;
    double TI = 0.0;

    double UT = 0.0;

    double dx = 0.1;
    double dy = 0.1;

    double tol = 0.00000001;
    int** geometry;
    #pragma endregion

    //Getting log file names
    #pragma region
    // command line arrays
    int argc = 3;
    const char* argv[] = { "./sim","-f","-t"};
    // File paths string 
    std::string szFileName;
    std::string geoFileName;
    std::string problemName;
    std::vector<std::string> logFilesName;
    // Problem Flags
    bool temp_flag;

    logFilesName  = ParseTestingProblemFlags(argc,argv,szFileName,geoFileName,problemName,temp_flag);
    #pragma endregion
   
    // Geometry Initialization
    #pragma region 
    geometry = read_pgm((geoFileName).c_str());
    #pragma endregion

    // Initialize Dummy Grid
    #pragma region 
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry);

    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();
    #pragma endregion

    
    // Data Structure Declaration
    #pragma region 
    // Input Data Structures
    matrix<double> U;
    matrix<double> V;
    // Calculation Data Structures
    matrix<double> F_cal;
    matrix<double> G_cal;
    // Reference Data Structures
    matrix<double> F_ref;
    matrix<double> G_ref;
    #pragma endregion

    // Data Structure Allocation
    #pragma region 
    // Input Data Structures
    U.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    F_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    

    // Allocating memory for matrices
    F_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    #pragma endregion

    // Reading log files in data structures
    #pragma region 
    // Reading dt
    std::ifstream fin(logFilesName[static_cast<int>(log_file::log_dt)]);
    fin >> dt;
    // Input Data Structures  
    read_matrix(logFilesName[static_cast<int>(log_file::log_U_Input)],U,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_V_Input)],V,imaxb,jmaxb);
    // Reference Data Structures
    read_matrix(logFilesName[static_cast<int>(log_file::log_F)],F_ref,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_G)],G_ref,imaxb,jmaxb);
    #pragma endregion

    // Dummy Scenario
    #pragma region
    
    domain.set_velocity(U, velocity_type::U);
    domain.set_velocity(V, velocity_type::V);
    
    calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, domain, F_cal, G_cal);

    #pragma endregion
    // Writing calculated matrices to log files
    #pragma region 
    write_matrix(cmp_log_F_t, F_cal, imaxb, jmaxb);
    write_matrix(cmp_log_G_t, G_cal, imaxb, jmaxb);
    #pragma endregion
    // Comparing Reference matrices to Calculated matrices
    #pragma region 

    for(int i = 0; i < imaxb; i++)
    {
      for(int j = 0; j < jmaxb; j++)
      {
        INFO("i = " << i << " j = " << j);
        REQUIRE(fabs(F_cal[i][j] - F_ref[i][j]) < tol);
        REQUIRE(fabs(G_cal[i][j] - G_ref[i][j]) < tol);
      }
    }

    #pragma endregion
} 

TEST_CASE( "Test rhs of PPE matrix Flow Over A Step", "[calculate_RS_FS]" )
{
  // Dummy Grid Paramters
    #pragma region 
    
    double xlength = 10;
    double ylength = 2;

    int imax = 98;
    int jmax = 18;

    double dt = 0.05  ;
    double t_end = 500.0;
    double tau = 0.5;

    double dt_Value = 10.0;

    int itermax = 500;
    double eps = 0.001;
    double omg = 1.7;
    double alpha = 0.9;

    double Re = 100;
    double GX = 0.0;
    double GY = 0.0;

    double UI = 0.0;
    double VI = 0.0;
    double PI = 0.0;
    double TI = 0.0;

    double UT = 0.0;

    double dx = 0.1;
    double dy = 0.1;

    double tol = 0.00000001;
    int** geometry;
    #pragma endregion

    //Getting log file names
    #pragma region
    // command line arrays
    int argc = 3;
    const char* argv[] = { "./sim","-f","-t"};
    // File paths string 
    std::string szFileName;
    std::string geoFileName;
    std::string problemName;
    std::vector<std::string> logFilesName;
    // Problem Flags
    bool temp_flag;

    logFilesName  = ParseTestingProblemFlags(argc,argv,szFileName,geoFileName,problemName,temp_flag);
    #pragma endregion
   
    // Geometry Initialization
    #pragma region 
    geometry = read_pgm((geoFileName).c_str());
    #pragma endregion
    // Initialize Dummy Grid
    #pragma region 
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry);

    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();
    #pragma endregion

    
    // Data Structure Declaration
    #pragma region 
    // Input Data Structures
    matrix<double> F;
    matrix<double> G;
    // Calculation Data Structures
    matrix<double> RS_cal;
    // Reference Data Structures
    matrix<double> RS_ref;
    #pragma endregion

    // Data Structure Allocation
    #pragma region 
    // Input Data Structures
    F.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    RS_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    RS_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    #pragma endregion

    // Reading log files in data structures
    #pragma region 
    // // File and directory paths
    // std::string logFileName;
    
    // // Input Data Structures
    // logFileName = Ref_Test_Data_t + "flog_F" + "1";
    // read_matrix(logFileName.c_str(),F,imaxb,jmaxb);
    // logFileName = Ref_Test_Data_t + "flog_G" + "1";
    // read_matrix(logFileName.c_str(),G,imaxb,jmaxb);

    // // Reference Data Structures
    // logFileName = Ref_Test_Data_t + "flog_RS" + "1";
    // read_matrix(logFileName.c_str(),RS_ref,imaxb,jmaxb);


    // Reading dt
    std::ifstream fin(logFilesName[static_cast<int>(log_file::log_dt)]);
    fin >> dt;    
    // Input Data Structures
    read_matrix(logFilesName[static_cast<int>(log_file::log_F)],F,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_G)],G,imaxb,jmaxb);

    // Reference Data Structures
    read_matrix(logFilesName[static_cast<int>(log_file::log_RS)],RS_ref,imaxb,jmaxb);
    #pragma endregion

    // Dummy Scenario
    #pragma region
    
    
    
    calculate_rs(dt, dx, dy, imax, jmax, F, G, RS_cal,domain);

    #pragma endregion
    // Writing calculated matrices to log files
    #pragma region 
    write_matrix(cmp_log_RS_t, RS_cal, imaxb, jmaxb);
    
    #pragma endregion
    // Comparing Reference matrices to Calculated matrices
    #pragma region 
    for(int i = 0; i < imaxb; i++)
    {
      for(int j = 0; j < jmaxb; j++)
      {
        if(RS_ref[i][j] == 0 )
        {
          // std::cout << "i = " << i << " j = " << j << std::endl;
            INFO("i = " << i << " j = " << j);
            REQUIRE(RS_cal[i][j] == RS_ref[i][j]);
        }
      }
    }

    for(int i = 0; i < imaxb; i++)
    {
      for(int j = 0; j < jmaxb; j++)
      {
        INFO("i = " << i << " j = " << j);
        REQUIRE(fabs(RS_cal[i][j] - RS_ref[i][j]) < tol);
        
      }
    }

    #pragma endregion
}

TEST_CASE( "Test Test SOR Flow Over A Step", "[sor_FS]" )
{
  // Dummy Grid Paramters
    #pragma region 
    
    double xlength = 10;
    double ylength = 2;

    int imax = 98;
    int jmax = 18;

    double dt = 0.05  ;
    double t_end = 500.0;
    double tau = 0.5;

    double dt_Value = 10.0;

    int itermax = 500;
    double eps = 0.001;
    double omg = 1.7;
    double alpha = 0.9;

    double Re = 100;
    double GX = 0.0;
    double GY = 0.0;

    double UI = 0.0;
    double VI = 0.0;
    double PI = 0.0;
    double TI = 0.0;

    double UT = 0.0;

    double dx = 0.1;
    double dy = 0.1;

    double tol = 0.00000001;
    int** geometry;
    #pragma endregion

   //Getting log file names
    #pragma region
    // command line arrays
    int argc = 3;
    const char* argv[] = { "./sim","-f","-t"};
    // File paths string 
    std::string szFileName;
    std::string geoFileName;
    std::string problemName;
    std::vector<std::string> logFilesName;
    // Problem Flags
    bool temp_flag;

    logFilesName  = ParseTestingProblemFlags(argc,argv,szFileName,geoFileName,problemName,temp_flag);
    #pragma endregion
   
    // Geometry Initialization
    #pragma region 
    geometry = read_pgm((geoFileName).c_str());
    #pragma endregion
    // Initialize Dummy Grid
    #pragma region 
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry);

    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();
    #pragma endregion

    
    
    // Data Structure Declaration
    #pragma region 
    // Input Data Structures
    matrix<double> U;
    matrix<double> V;
    matrix<double> RS;
    // Calculation Data Structures
    matrix<double> P_cal;
    // Reference Data Structures
    matrix<double> P_ref;
    #pragma endregion

    // Data Structure Allocation
    #pragma region 
    // Input Data Structures
    U.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V.resize(imaxb,std::vector<double>(jmaxb,0.0));
    RS.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    P_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    P_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    #pragma endregion

    // Reading log files in data structures
    #pragma region 
    // File and directory paths
    // std::string logFileName;
    
    // // Input Data Structures
    // logFileName = Ref_Test_Data_t + "flog_U_Input" + "1";
    // read_matrix(logFileName.c_str(),U,imaxb,jmaxb);
    // logFileName = Ref_Test_Data_t + "flog_V_Input" + "1";
    // read_matrix(logFileName.c_str(),V,imaxb,jmaxb);
    // logFileName = Ref_Test_Data_t + "flog_RS" + "1";
    // read_matrix(logFileName.c_str(),RS,imaxb,jmaxb);
    // logFileName = Ref_Test_Data_t + "flog_P_Input" + "1";
    // read_matrix(logFileName.c_str(),P_cal,imaxb,jmaxb);
    

    // // Reference Data
    // logFileName = Ref_Test_Data_t + "flog_P_Output" + "1";
    // read_matrix(logFileName.c_str(),P_ref,imaxb,jmaxb);

    // Reading dt
    std::ifstream fin(logFilesName[static_cast<int>(log_file::log_dt)]);
    fin >> dt; 
    //Input Data    
    read_matrix(logFilesName[static_cast<int>(log_file::log_U_Input)],U,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_V_Input)],V,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_RS)],RS,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_P_Input)],P_cal,imaxb,jmaxb);

    // Reference Data
    read_matrix(logFilesName[static_cast<int>(log_file::log_P_Output)],P_ref,imaxb,jmaxb);

    #pragma endregion

    // Dummy Scenario
    #pragma region
    domain.set_velocity(U,velocity_type::U);
    domain.set_velocity(V,velocity_type::V);
    domain.set_pressure(P_cal);
    int it = 0;
    double res = 1.0;
    // Solve system using SOR
    while(it < itermax && res > eps){
      sor(omg,dx,dy,imax,jmax,domain,RS,&res);
      it++;
    }
    

    #pragma endregion
   
    // Writing calculated matrices to log files
    #pragma region 
    domain.pressure(P_cal);
    write_matrix(cmp_log_P_t, P_cal, imaxb, jmaxb);
    
    #pragma endregion
   
    // Comparing Reference matrices to Calculated matrices
    #pragma region 
   
    for(int i = 0; i < imaxb; i++)
    {
      for(int j = 0; j < jmaxb; j++)
      {
        INFO("i = " << i << " j = " << j);
        REQUIRE(fabs(P_cal[i][j] - P_ref[i][j]) < tol);
        
      }
    }

    #pragma endregion
}

TEST_CASE( "Test velocity (U and v) values  Flow Over A Step", "[calculate_uv_FS]" )
{
  // Dummy Grid Paramters
    #pragma region 
    
    double xlength = 10;
    double ylength = 2;

    int imax = 98;
    int jmax = 18;

    double dt = 0.05  ;
    double t_end = 500.0;
    double tau = 0.5;

    double dt_Value = 10.0;

    int itermax = 500;
    double eps = 0.001;
    double omg = 1.7;
    double alpha = 0.9;

    double Re = 100;
    double GX = 0.0;
    double GY = 0.0;

    double UI = 0.0;
    double VI = 0.0;
    double PI = 0.0;
    double TI = 0.0;

    double UT = 0.0;

    double dx = 0.1;
    double dy = 0.1;

    double tol = 0.00000001;
    int** geometry;
    #pragma endregion

   //Getting log file names
    #pragma region
    // command line arrays
    int argc = 3;
    const char* argv[] = { "./sim","-f","-t"};
    // File paths string 
    std::string szFileName;
    std::string geoFileName;
    std::string problemName;
    std::vector<std::string> logFilesName;
    // Problem Flags
    bool temp_flag;

    logFilesName  = ParseTestingProblemFlags(argc,argv,szFileName,geoFileName,problemName,temp_flag);
    #pragma endregion
   
    // Geometry Initialization
    #pragma region 
    geometry = read_pgm((geoFileName).c_str());
    #pragma endregion
    // Initialize Dummy Grid
    #pragma region 
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry);

    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();
    #pragma endregion

  // Data Structure Declaration
    #pragma region 
    // Input Data Structures
    matrix<double> F;
    matrix<double> G;
    matrix<double> P;
    // Calculation Data Structures
    matrix<double> U_cal;
    matrix<double> V_cal;
    // Reference Data Structures
    matrix<double> U_ref;
    matrix<double> V_ref;
    #pragma endregion

    // Data Structure Allocation
    #pragma region 
    // Input Data Structures
    F.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G.resize(imaxb,std::vector<double>(jmaxb,0.0));
    P.resize(imaxb,std::vector<double>(jmaxb,0.0));    
    // Calculation Data Structures
    U_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    U_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    #pragma endregion

    // Reading log files in data structures
    #pragma region 
    // // File and directory paths
    // std::string logFileName;

    // // Input Data Structures
    // logFileName = Ref_Test_Data_t + "flog_F" + "1";
    // read_matrix(logFileName.c_str(),F,imaxb,jmaxb);
    // logFileName = Ref_Test_Data_t + "flog_G" + "1";
    // read_matrix(logFileName.c_str(),G,imaxb,jmaxb);
    // logFileName = Ref_Test_Data_t + "flog_P_Output" + "1";
    // read_matrix(logFileName.c_str(),P,imaxb,jmaxb);
    // logFileName = Ref_Test_Data_t + "flog_U_Input" + "1";
    // read_matrix(logFileName.c_str(),U_cal,imaxb,jmaxb);
    // logFileName = Ref_Test_Data_t + "flog_V_Input" + "1";
    // read_matrix(logFileName.c_str(),V_cal,imaxb,jmaxb);



    // // Reference Data Structures
    // logFileName = Ref_Test_Data_t + "flog_U_Output" + "1";
    // read_matrix(logFileName.c_str(),U_ref,imaxb,jmaxb);
    // logFileName = Ref_Test_Data_t + "flog_V_Output" + "1";
    // read_matrix(logFileName.c_str(),V_ref,imaxb,jmaxb);

    // Reading dt
    std::ifstream fin(logFilesName[static_cast<int>(log_file::log_dt)]);
    fin >> dt;
    // Input Data
    read_matrix(logFilesName[static_cast<int>(log_file::log_F)],F,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_G)],G,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_P_Output)],P,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_U_Input)],U_cal,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_V_Input)],V_cal,imaxb,jmaxb);

    // Reference Data Structures
    read_matrix(logFilesName[static_cast<int>(log_file::log_U_Output)],U_ref,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_V_Output)],V_ref,imaxb,jmaxb);
 
    #pragma endregion

    // Dummy Scenario
    #pragma region
    domain.set_velocity(U_cal,velocity_type::U);
    domain.set_velocity(V_cal,velocity_type::V);
    domain.set_pressure(P);
    calculate_uv(dt, dx, dy, imax, jmax, domain, F, G);

    #pragma endregion
    
    // Writing calculated matrices to log files
    #pragma region 
    domain.velocity(U_cal,velocity_type::U);
    domain.velocity(V_cal,velocity_type::V);
    write_matrix(cmp_log_U_t, U_cal , imaxb, jmaxb);
    write_matrix(cmp_log_V_t, V_cal , imaxb, jmaxb);
    #pragma endregion
    
    // Comparing Reference matrices to Calculated matrices
    #pragma region 
    // for(int i = 0; i < imaxb; i++)
    // {
    //   for(int j = 0; j < jmaxb; j++)
    //   {
    //     if(U_ref[i][j] == 0 )
    //     {
    //       // std::cout << "i = " << i << " j = " << j << std::endl;
    //         INFO("i = " << i << " j = " << j);
    //         REQUIRE(U_cal[i][j] == U_ref[i][j]);
    //     }
    //     if(V_ref[i][j] == 0 )
    //     {
    //       // std::cout << "i = " << i << " j = " << j << std::endl;
    //         INFO("i = " << i << " j = " << j);
    //         REQUIRE(V_cal[i][j] == V_ref[i][j]);
    //     }
        
    //   }
    // }

    for(int i = 1; i < imaxb-1; i++)
    {
      for(int j = 1; j < jmaxb-1; j++)
      {
        INFO("i = " << i << " j = " << j);
        REQUIRE(fabs(U_cal[i][j] - U_ref[i][j]) < tol);
        REQUIRE(fabs(V_cal[i][j] - V_ref[i][j]) < tol);
      }
    }

    #pragma endregion
}

#pragma endregion

// Natural Convection Validation Tests
#pragma region


TEST_CASE("Test calculate dt Natural Convection", "[calculate_dt_NC]")
{
    // Dummy Grid Parameters
    #pragma region 
    double xlength = 1;
    double ylength = 1;

    int imax = 48;
    int jmax = 48;

    double dt = 0.05 ;
    double dt_ref = 0.0;
    double t_end = 1000.0;
    double tau = 0.5;

    double dt_Value = 10.0;

    int itermax = 100;
    double eps = 0.00001;
    double omg = 1.7;
    double alpha = 0.5;
    double beta = 0.00021;

    double Re = 1000;
    double PR = 7;
    double GX = 0.0;
    double GY = -1.1;

    double UI = 0.0;
    double VI = 0.0;
    double PI = 0.0;
    double TI = 0.0;

    double UT = 0.0;

    double dx = 0.02;
    double dy = 0.02;

    double tol = 0.0000001;
    int** geometry;
    #pragma endregion
    

    //Getting log file names
    #pragma region
    // command line arrays
    int argc = 3;
    const char* argv[] = { "./sim","-n","-t"};
    // File paths string 
    std::string szFileName;
    std::string geoFileName;
    std::string problemName;
    std::vector<std::string> logFilesName;
    // Problem Flags
    bool temp_flag;

    logFilesName  = ParseTestingProblemFlags(argc,argv,szFileName,geoFileName,problemName,temp_flag);
    #pragma endregion

    // Geometry Initialization
    #pragma region 
    geometry = read_pgm((geoFileName).c_str());
    #pragma endregion

    // Initialzie Dummy Grid
    #pragma region    
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry) ;
    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();
    #pragma endregion

    // Data Structure Declaration
    #pragma region 
    matrix<double> U;
    matrix<double> V;
    matrix<double> P;
    matrix<double> T;
    #pragma endregion

  
    // Data Structure Allocation
    #pragma region 
    // Input Data Structures
    U.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V.resize(imaxb,std::vector<double>(jmaxb,0.0));
    T.resize(imaxb,std::vector<double>(jmaxb,0.0));
    #pragma endregion
     
    // Read Input matrices 
    #pragma region 
    // Reading dt
    std::ifstream fin(logFilesName[static_cast<int>(log_file::log_dt)]);
    fin >> dt_ref;

    // Input Data Structures  
    read_matrix(logFilesName[static_cast<int>(log_file::log_U_Input)],U,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_V_Input)],V,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_T_Input)],T,imaxb,jmaxb);

    #pragma endregion
   

    /* Real Scenario*/
    #pragma region 
  
    domain.set_velocity(U, velocity_type::U);
    domain.set_velocity(V, velocity_type::V);
    domain.set_temperature(T);    

    calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,PR,&temp_flag,domain);

    
    #pragma endregion

    #pragma region 
    REQUIRE(fabs(dt - dt_ref) < tol );
    #pragma endregion
}

TEST_CASE( "Test F and G matrix Natural Convection", "[calculate_fg_NC]" )
{
  // Dummy Grid Paramters
    #pragma region 
    
    double xlength = 1;
    double ylength = 1;

    int imax = 48;
    int jmax = 48;

    double dt = 0.05;
    double t_end = 1000.0;
    double tau = 0.5;

    double dt_Value = 10.0;

    int itermax = 100;
    double eps = 0.00001;
    double omg = 1.7;
    double alpha = 0.5;
    double beta = 0.00021;

    double Re = 1000;
    double PR = 7;
    double GX = 0.0;
    double GY = -1.1;

    double UI = 0.0;
    double VI = 0.0;
    double PI = 0.0;
    double TI = 0.0;

    double UT = 0.0;

    double dx = 0.02;
    double dy = 0.02;

    double tol = 0.0000001;
    int** geometry;
    #pragma endregion

    //Getting log file names
    #pragma region
    // command line arrays
    int argc = 3;
    const char* argv[] = { "./sim","-n","-t"};
    // File paths string 
    std::string szFileName;
    std::string geoFileName;
    std::string problemName;
    std::vector<std::string> logFilesName;
    // Problem Flags
    bool temp_flag;

    logFilesName  = ParseTestingProblemFlags(argc,argv,szFileName,geoFileName,problemName,temp_flag);
    #pragma endregion

    // Geometry Initialization
    #pragma region 
    geometry = read_pgm((geoFileName).c_str());
    #pragma endregion

    // Initialize Dummy Grid
    #pragma region 
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry);

    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();
    #pragma endregion

    #pragma region 
    // Input Data Structures
    matrix<double> U;
    matrix<double> V;
    matrix<double> T;
    // Calculation Data Structures
    matrix<double> F_cal;
    matrix<double> G_cal;
    // Reference Data Structures
    matrix<double> F_ref;
    matrix<double> G_ref;
    #pragma endregion

    // Data Structure Allocation
    #pragma region 
    // Input Data Structures
    U.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V.resize(imaxb,std::vector<double>(jmaxb,0.0));
    T.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    F_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    
    // Allocating memory for matrices
    F_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    
    #pragma endregion

    // Reading log files in data structures
    #pragma region 
    // Reading dt
    std::ifstream fin(logFilesName[static_cast<int>(log_file::log_dt)]);
    fin >> dt;
    // Input Data Structures  
    read_matrix(logFilesName[static_cast<int>(log_file::log_U_Input)],U,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_V_Input)],V,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_T_Input)],T,imaxb,jmaxb);
    // Reference Data Structures
    read_matrix(logFilesName[static_cast<int>(log_file::log_F)],F_ref,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_G)],G_ref,imaxb,jmaxb);
 
    #pragma endregion

    // Dummy Scenario
    #pragma region
    
    domain.set_velocity(U, velocity_type::U);
    domain.set_velocity(V, velocity_type::V);
    domain.set_temperature(T);

    calculate_fg_temp(Re,GX,GY,alpha,dt,dx,dy,imax,jmax,beta,domain,F_cal,G_cal);

    #pragma endregion

     #pragma region 
    write_matrix(cmp_log_F_t, F_cal, imaxb, jmaxb);
    write_matrix(cmp_log_G_t, G_cal, imaxb, jmaxb);
    #pragma endregion

    // Comparing Reference matrices to Calculated matrices
    #pragma region 
    

    for(int i = 0; i < imaxb; i++)
    {
      for(int j = 0; j < jmaxb; j++)
      {
        INFO("i = " << i << " j = " << j);
        REQUIRE(fabs(F_cal[i][j] - F_ref[i][j]) < tol);
        REQUIRE(fabs(G_cal[i][j] - G_ref[i][j]) < tol);
        
      }
    }

    #pragma endregion
}

TEST_CASE( "Test temperature values Natural Convection", "[calculate_temp_NC]" )
{
  // Dummy Grid Paramters
  #pragma region 
    
    double xlength = 1;
    double ylength = 1;

    int imax = 48;
    int jmax = 48;

    double dt = 0.0500000000 ;
    double t_end = 1000.0;
    double tau = 0.5;

    double dt_Value = 10.0;

    int itermax = 100;
    double eps = 0.00001;
    double omg = 1.7;
    double alpha = 0.5;
    double beta = 0.00021;

    double Re = 1000;
    double PR = 7;
    double GX = 0.0;
    double GY = -1.1;

    double UI = 0.0;
    double VI = 0.0;
    double PI = 0.0;
    double TI = 0.0;

    double UT = 0.0;

    double dx = 0.02;
    double dy = 0.02;

    double tol = 0.0000001;
    int** geometry;
    #pragma endregion
  
  //Getting log file names
    #pragma region
    // command line arrays
    int argc = 3;
    const char* argv[] = { "./sim","-n","-t"};
    // File paths string 
    std::string szFileName;
    std::string geoFileName;
    std::string problemName;
    std::vector<std::string> logFilesName;
    // Problem Flags
    bool temp_flag;

    logFilesName  = ParseTestingProblemFlags(argc,argv,szFileName,geoFileName,problemName,temp_flag);
    #pragma endregion

    // Geometry Initialization
    #pragma region 
    geometry = read_pgm((geoFileName).c_str());
    #pragma endregion
  // Initialize Dummy Grid
  #pragma region 
  Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry);
  int imaxb = domain.imaxb();
  int jmaxb = domain.jmaxb();
  #pragma endregion
  
  // Data Structure Declaration
  #pragma region 
  // Input Data Structures
  matrix<double> U;
  matrix<double> V;
  matrix<double> T;
  // Calculation Data Structures
  matrix<double> T_cal;
  // Reference Data Structures
  matrix<double> T_ref;
  #pragma endregion
  
  // Data Structure Allocation
  #pragma region 
    // Input Data Structures
    U.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V.resize(imaxb,std::vector<double>(jmaxb,0.0));
    T.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    T_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    
    // Allocating memory for matrices
    T_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    
    #pragma endregion
 
  // Reading log files in data structures
  #pragma region 
  // Reading dt
    std::ifstream fin(logFilesName[static_cast<int>(log_file::log_dt)]);
    fin >> dt;
    // Input Data Structures  
    read_matrix(logFilesName[static_cast<int>(log_file::log_U_Input)],U,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_V_Input)],V,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_T_Input)],T,imaxb,jmaxb);
    // Reference Data Structures
    read_matrix(logFilesName[static_cast<int>(log_file::log_T_Output)],T_ref,imaxb,jmaxb);
  #pragma endregion
  
  // Dummy Scenario
  #pragma region  
  domain.set_velocity(U, velocity_type::U);
  domain.set_velocity(V, velocity_type::V);
  domain.set_temperature(T);
  calculate_temp(dt,dx,dy,imax,jmax,Re,PR,alpha,0,domain);
  #pragma endregion
  
  // Write Output Data
  #pragma region 
  domain.temperature(T_cal);
  write_matrix(cmp_log_T_t, T_cal, imaxb, jmaxb);
  #pragma endregion
  
  // Comparing Reference matrices to Calculated matrices
  #pragma region 
    

    for(int i = 0; i < imaxb; i++)
    {
      for(int j = 0; j < jmaxb; j++)
      {
        INFO("i = " << i << " j = " << j);
        REQUIRE(fabs(T_cal[i][j] - T_ref[i][j]) < tol);
        
      }
    }

    #pragma endregion
}

TEST_CASE( "Test rhs of PPE matrix Natural Convection", "[calculate_RS_NC]" )
{
  // Dummy Grid Paramters
    #pragma region 
    
    double xlength = 1;
    double ylength = 1;

    int imax = 48;
    int jmax = 48;

    double dt = 0.0500000000 ;
    double t_end = 1000.0;
    double tau = 0.5;

    double dt_Value = 10.0;

    int itermax = 100;
    double eps = 0.00001;
    double omg = 1.7;
    double alpha = 0.5;
    double beta = 0.00021;

    double Re = 1000;
    double PR = 7;
    double GX = 0.0;
    double GY = -1.1;

    double UI = 0.0;
    double VI = 0.0;
    double PI = 0.0;
    double TI = 0.0;

    double UT = 0.0;

    double dx = 0.02;
    double dy = 0.02;

    double tol = 0.0000001;
    int** geometry;
    #pragma endregion

    //Getting log file names
    #pragma region
    // command line arrays
    int argc = 3;
    const char* argv[] = { "./sim","-n","-t"};
    // File paths string 
    std::string szFileName;
    std::string geoFileName;
    std::string problemName;
    std::vector<std::string> logFilesName;
    // Problem Flags
    bool temp_flag;

    logFilesName  = ParseTestingProblemFlags(argc,argv,szFileName,geoFileName,problemName,temp_flag);
    #pragma endregion

    // Geometry Initialization
    #pragma region 
    geometry = read_pgm((geoFileName).c_str());
    #pragma endregion
    
    // Initialize Dummy Grid
    #pragma region 
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry);

    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();
    #pragma endregion

    
    // Data Structure Declaration
    #pragma region 
    // Input Data Structures
    matrix<double> F;
    matrix<double> G;
    // Calculation Data Structures
    matrix<double> RS_cal;
    // Reference Data Structures
    matrix<double> RS_ref;
    #pragma endregion

    // Data Structure Allocation
    #pragma region 
    // Input Data Structures
    F.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    RS_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    RS_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    #pragma endregion

    // Reading log files in data structures
    #pragma region 
    // Redaing dt
    std::ifstream fin(logFilesName[static_cast<int>(log_file::log_dt)]);
    fin >> dt;    
    // Input Data Structures
    read_matrix(logFilesName[static_cast<int>(log_file::log_F)],F,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_G)],G,imaxb,jmaxb);

    // Reference Data
    read_matrix(logFilesName[static_cast<int>(log_file::log_RS)],RS_ref,imaxb,jmaxb);

    #pragma endregion

    // Dummy Scenario
    #pragma region
    
    
    
    calculate_rs(dt, dx, dy, imax, jmax, F, G, RS_cal,domain);

    #pragma endregion
    // Writing calculated matrices to log files
    #pragma region 
    write_matrix(cmp_log_RS_t, RS_cal, imaxb, jmaxb);
    
    #pragma endregion
    // Comparing Reference matrices to Calculated matrices
    #pragma region 
    for(int i = 0; i < imaxb; i++)
    {
      for(int j = 0; j < jmaxb; j++)
      {
        if(RS_ref[i][j] == 0 )
        {
          // std::cout << "i = " << i << " j = " << j << std::endl;
            INFO("i = " << i << " j = " << j);
            REQUIRE(RS_cal[i][j] == RS_ref[i][j]);
        }
      }
    }

    for(int i = 0; i < imaxb; i++)
    {
      for(int j = 0; j < jmaxb; j++)
      {
        INFO("i = " << i << " j = " << j);
        REQUIRE(fabs(RS_cal[i][j] - RS_ref[i][j]) < tol);
        
      }
    }

    #pragma endregion
}

TEST_CASE( "Test Test SOR Natural Convection", "[sor_NC]" )
{
  // Dummy Grid Paramters
    #pragma region 
    
    double xlength = 1;
    double ylength = 1;

    int imax = 48;
    int jmax = 48;

    double dt = 0.0500000000 ;
    double t_end = 1000.0;
    double tau = 0.5;

    double dt_Value = 10.0;

    int itermax = 100;
    double eps = 0.00001;
    double omg = 1.7;
    double alpha = 0.5;
    double beta = 0.00021;

    double Re = 1000;
    double PR = 7;
    double GX = 0.0;
    double GY = -1.1;

    double UI = 0.0;
    double VI = 0.0;
    double PI = 0.0;
    double TI = 0.0;

    double UT = 0.0;

    double dx = 0.02;
    double dy = 0.02;

    double tol = 0.0000001;
    int** geometry;
    #pragma endregion

    //Getting log file names
    #pragma region
    // command line arrays
    int argc = 3;
    const char* argv[] = { "./sim","-n","-t"};
    // File paths string 
    std::string szFileName;
    std::string geoFileName;
    std::string problemName;
    std::vector<std::string> logFilesName;
    // Problem Flags
    bool temp_flag;

    logFilesName  = ParseTestingProblemFlags(argc,argv,szFileName,geoFileName,problemName,temp_flag);
    #pragma endregion

    // Geometry Initialization
    #pragma region 
    geometry = read_pgm((geoFileName).c_str());
    #pragma endregion
    // Initialize Dummy Grid
    #pragma region 
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry);

    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();
    #pragma endregion

    
    // Data Structure Declaration
    #pragma region 
    // Input Data Structures
    matrix<double> U;
    matrix<double> V;
    matrix<double> RS;
    // Calculation Data Structures
    matrix<double> P_cal;
    // Reference Data Structures
    matrix<double> P_ref;
    #pragma endregion

    // Data Structure Allocation
    #pragma region 
    // Input Data Structures
    U.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V.resize(imaxb,std::vector<double>(jmaxb,0.0));
    RS.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    P_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    P_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    #pragma endregion

    // Reading log files in data structures
    #pragma region 
    // Reading dt
    std::ifstream fin(logFilesName[static_cast<int>(log_file::log_dt)]);
    fin >> dt; 
    //Input Data    
    read_matrix(logFilesName[static_cast<int>(log_file::log_U_Input)],U,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_V_Input)],V,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_RS)],RS,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_P_Input)],P_cal,imaxb,jmaxb);

    // Reference Data
    read_matrix(logFilesName[static_cast<int>(log_file::log_P_Output)],P_ref,imaxb,jmaxb);
    #pragma endregion

    // Dummy Scenario
    #pragma region
    domain.set_velocity(U,velocity_type::U);
    domain.set_velocity(V,velocity_type::V);
    domain.set_pressure(P_cal);
    int it = 0;
    double res = 1.0;
    // Solve system using SOR
    while(it < itermax && res > eps){
      sor(omg,dx,dy,imax,jmax,domain,RS,&res);
      it++;
    }
    

    #pragma endregion
   
    // Writing calculated matrices to log files
    #pragma region 
    domain.pressure(P_cal);
    write_matrix(cmp_log_P_t, P_cal, imaxb, jmaxb);
    
    #pragma endregion
   
    // Comparing Reference matrices to Calculated matrices
    #pragma region 
   
    for(int i = 0; i < imaxb; i++)
    {
      for(int j = 0; j < jmaxb; j++)
      {
        INFO("i = " << i << " j = " << j);
        REQUIRE(fabs(P_cal[i][j] - P_ref[i][j]) < tol);
        
      }
    }

    #pragma endregion
}

TEST_CASE( "Test velocity (U and v) values Natural Convection", "[calculate_uv_NC]" )
{
  // Dummy Grid Paramters
  #pragma region 
    
    double xlength = 1;
    double ylength = 1;

    int imax = 48;
    int jmax = 48;

    double dt = 0.0500000000 ;
    double t_end = 1000.0;
    double tau = 0.5;

    double dt_Value = 10.0;

    int itermax = 100;
    double eps = 0.00001;
    double omg = 1.7;
    double alpha = 0.5;
    double beta = 0.00021;

    double Re = 1000;
    double PR = 7;
    double GX = 0.0;
    double GY = -1.1;

    double UI = 0.0;
    double VI = 0.0;
    double PI = 0.0;
    double TI = 0.0;

    double UT = 0.0;

    double dx = 0.02;
    double dy = 0.02;

    double tol = 0.0000001;
    int** geometry;
    #pragma endregion
  
  //Getting log file names
    #pragma region
    // command line arrays
    int argc = 3;
    const char* argv[] = { "./sim","-n","-t"};
    // File paths string 
    std::string szFileName;
    std::string geoFileName;
    std::string problemName;
    std::vector<std::string> logFilesName;
    // Problem Flags
    bool temp_flag;

    logFilesName  = ParseTestingProblemFlags(argc,argv,szFileName,geoFileName,problemName,temp_flag);
    #pragma endregion

    // Geometry Initialization
    #pragma region 
    geometry = read_pgm((geoFileName).c_str());
    #pragma endregion
  // Initialize Dummy Grid
  #pragma region 
  Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry);
  int imaxb = domain.imaxb();
  int jmaxb = domain.jmaxb();
  #pragma endregion
  
  // Data Structure Declaration
    #pragma region 
    // Input Data Structures
    matrix<double> F;
    matrix<double> G;
    matrix<double> P;
    // Calculation Data Structures
    matrix<double> U_cal;
    matrix<double> V_cal;
    // Reference Data Structures
    matrix<double> U_ref;
    matrix<double> V_ref;
    #pragma endregion

    // Data Structure Allocation
    #pragma region 
    // Input Data Structures
    F.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G.resize(imaxb,std::vector<double>(jmaxb,0.0));
    P.resize(imaxb,std::vector<double>(jmaxb,0.0));    
    // Calculation Data Structures
    U_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    U_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    #pragma endregion

    // Reading log files in data structures
    #pragma region 
    // Reading dt
    std::ifstream fin(logFilesName[static_cast<int>(log_file::log_dt)]);
    fin >> dt;
    // Input Data
    read_matrix(logFilesName[static_cast<int>(log_file::log_F)],F,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_G)],G,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_P_Output)],P,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_U_Input)],U_cal,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_V_Input)],V_cal,imaxb,jmaxb);

    // Reference Data Structures
    read_matrix(logFilesName[static_cast<int>(log_file::log_U_Output)],U_ref,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_V_Output)],V_ref,imaxb,jmaxb);
 
    #pragma endregion

    // Dummy Scenario
    #pragma region
    domain.set_velocity(U_cal,velocity_type::U);
    domain.set_velocity(V_cal,velocity_type::V);
    domain.set_pressure(P);
    calculate_uv(dt, dx, dy, imax, jmax, domain, F, G);

    #pragma endregion
    
    // Writing calculated matrices to log files
    #pragma region 
    domain.velocity(U_cal,velocity_type::U);
    domain.velocity(V_cal,velocity_type::V);
    write_matrix(cmp_log_U_t, U_cal , imaxb, jmaxb);
    write_matrix(cmp_log_V_t, V_cal , imaxb, jmaxb);
    #pragma endregion
    
    // Comparing Reference matrices to Calculated matrices
    #pragma region 
    for(int i = 0; i < imaxb; i++)
    {
      for(int j = 0; j < jmaxb; j++)
      {
        if(U_ref[i][j] == 0 )
        {
          // std::cout << "i = " << i << " j = " << j << std::endl;
            INFO("i = " << i << " j = " << j);
            REQUIRE(U_cal[i][j] == U_ref[i][j]);
        }
        if(V_ref[i][j] == 0 )
        {
          // std::cout << "i = " << i << " j = " << j << std::endl;
            INFO("i = " << i << " j = " << j);
            REQUIRE(V_cal[i][j] == V_ref[i][j]);
        }
        
      }
    }

    for(int i = 1; i < imaxb-1; i++)
    {
      for(int j = 1; j < jmaxb-1; j++)
      {
        INFO("i = " << i << " j = " << j);
        REQUIRE(fabs(U_cal[i][j] - U_ref[i][j]) < tol);
        REQUIRE(fabs(V_cal[i][j] - V_ref[i][j]) < tol);
      }
    }

    #pragma endregion
}

#pragma endregion

// Fluid Trap New Data Set
#pragma region 
// TEST_CASE("Test calculate dt Fluid Trap", "[calculate_dt_FT]")
// {
//    // Dummy Grid Parameters
//     #pragma region 
//     double xlength = 2;
//     double ylength = 1;

//     int imax = 98;
//     int jmax = 48;

//     double dt = 0.5;
//     double dt_ref = 0.0;
//     double t_end = 2000.0;
//     double tau = 0.5;

//     double dt_Value = 10.0;

//     int itermax = 1000;
//     double eps = 0.00001;
//     double omg = 1.7;
//     double alpha = 0.5;
//     double beta = 0.00063;

//     double Re = 10000;
//     double PR = 7;
//     double GX = 0.0;
//     double GY = -9.81;

//     double UI = 0.0;
//     double VI = 0.0;
//     double PI = 0.0;
//     double TI = 0.0;

//     double UT = 0.0;

//     double dx = 0.02;
//     double dy = 0.02;

//     double tol = 0.0000001;
//     int** geometry;
//     #pragma endregion

//     //Getting log file names
//     #pragma region
//     // command line arrays
//     int argc = 3;
//     const char* argv[] = { "./sim","-ft","-t"};
//     // File paths string 
//     std::string szFileName;
//     std::string geoFileName;
//     std::string problemName;
//     std::vector<std::string> logFilesName;
//     // Problem Flags
//     bool temp_flag;

//     logFilesName  = ParseTestingProblemFlags(argc,argv,szFileName,geoFileName,problemName,temp_flag);
//     #pragma endregion

//     // Geometry Initialization
//     #pragma region 
//     geometry = read_pgm((geoFileName).c_str());
//     #pragma endregion

//     // Initialzie Dummy Grid
//     #pragma region    
//     Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry) ;
//     int imaxb = domain.imaxb();
//     int jmaxb = domain.jmaxb();
//     #pragma endregion

//     // Data Structure Declaration
//     #pragma region 
//     matrix<double> U;
//     matrix<double> V;
//     matrix<double> P;
//     matrix<double> T;
//     matrix<double> Test;
//     #pragma endregion

  
//     // Data Structure Allocation
//     #pragma region 
//     // Input Data Structures
//     U.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     V.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     T.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     #pragma endregion
     
//     // Read Input matrices 
//     #pragma region 
//     // Reading dt
//     std::ifstream fin(logFilesName[static_cast<int>(log_file::log_dt)]);
//     fin >> dt_ref;

//     // Input Data Structures  
//     read_matrix(logFilesName[static_cast<int>(log_file::log_U_Input)],U,imaxb,jmaxb);
//     read_matrix(logFilesName[static_cast<int>(log_file::log_V_Input)],V,imaxb,jmaxb);
//     read_matrix(logFilesName[static_cast<int>(log_file::log_T_Input)],T,imaxb,jmaxb);

//     #pragma endregion
   

//     /* Real Scenario*/
//     #pragma region 
  
//     domain.set_velocity(U, velocity_type::U);
//     domain.set_velocity(V, velocity_type::V);
//     domain.set_temperature(T);    

//     calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,PR,&temp_flag,domain);

    
//     #pragma endregion

//     #pragma region 
//     REQUIRE(fabs(dt - dt_ref) < tol );
//     #pragma endregion
// }

// TEST_CASE( "Test F and G matrix Fluid Trap", "[calculate_fg_FT]" )
// {
//   // Dummy Grid Parameters
//     #pragma region 
//     double xlength = 2;
//     double ylength = 1;

//     int imax = 98;
//     int jmax = 48;

//     double dt = 0.500000000 ;
//     double t_end = 2000.0;
//     double tau = 0.5;

//     double dt_Value = 10.0;

//     int itermax = 1000;
//     double eps = 0.00001;
//     double omg = 1.7;
//     double alpha = 0.5;
//     double beta = 0.00063;

//     double Re = 10000;
//     double PR = 7;
//     double GX = 0.0;
//     double GY = -9.81;

//     double UI = 0.0;
//     double VI = 0.0;
//     double PI = 0.0;
//     double TI = 0.0;

//     double UT = 0.0;

//     double dx = 0.02;
//     double dy = 0.02;

//     double tol = 0.000001;
//     int** geometry;
//     #pragma endregion
//     //Getting log file names
//     #pragma region
//     // command line arrays
//     int argc = 3;
//     const char* argv[] = { "./sim","-ft","-t"};
//     // File paths string 
//     std::string szFileName;
//     std::string geoFileName;
//     std::string problemName;
//     std::vector<std::string> logFilesName;
//     // Problem Flags
//     bool temp_flag;

//     logFilesName  = ParseTestingProblemFlags(argc,argv,szFileName,geoFileName,problemName,temp_flag);
//     #pragma endregion

//     // Geometry Initialization
//     #pragma region 
//     geometry = read_pgm((geoFileName).c_str());
//     #pragma endregion

//     // Initialize Dummy Grid
//     #pragma region 
//     Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry);

//     int imaxb = domain.imaxb();
//     int jmaxb = domain.jmaxb();
//     #pragma endregion

//     #pragma region 
//     // Input Data Structures
//     matrix<double> U;
//     matrix<double> V;
//     matrix<double> T;
//     // Calculation Data Structures
//     matrix<double> F_cal;
//     matrix<double> G_cal;
//     // Reference Data Structures
//     matrix<double> F_ref;
//     matrix<double> G_ref;
//     #pragma endregion

//     // Data Structure Allocation
//     #pragma region 
//     // Input Data Structures
//     U.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     V.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     T.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     // Allocating memory for matrices
//     F_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     G_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    
//     // Allocating memory for matrices
//     F_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     G_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    
//     #pragma endregion

//     // Reading log files in data structures
//     #pragma region 
//     // Reading dt
//     std::ifstream fin(logFilesName[static_cast<int>(log_file::log_dt)]);
//     fin >> dt;
//     // Input Data Structures  
//     read_matrix(logFilesName[static_cast<int>(log_file::log_U_Input)],U,imaxb,jmaxb);
//     read_matrix(logFilesName[static_cast<int>(log_file::log_V_Input)],V,imaxb,jmaxb);
//     read_matrix(logFilesName[static_cast<int>(log_file::log_T_Input)],T,imaxb,jmaxb);
//     // Reference Data Structures
//     read_matrix(logFilesName[static_cast<int>(log_file::log_F)],F_ref,imaxb,jmaxb);
//     read_matrix(logFilesName[static_cast<int>(log_file::log_G)],G_ref,imaxb,jmaxb);
 
//     #pragma endregion

//     // Dummy Scenario
//     #pragma region
    
//     domain.set_velocity(U, velocity_type::U);
//     domain.set_velocity(V, velocity_type::V);
//     domain.set_temperature(T);

//     calculate_fg_temp(Re,GX,GY,alpha,dt,dx,dy,imax,jmax,beta,domain,F_cal,G_cal);

//     #pragma endregion

//      #pragma region 
//     write_matrix(cmp_log_F_t, F_cal, imaxb, jmaxb);
//     write_matrix(cmp_log_G_t, G_cal, imaxb, jmaxb);
//     #pragma endregion

//     // Comparing Reference matrices to Calculated matrices
//     #pragma region 
    

//     for(int i = 0; i < imaxb; i++)
//     {
//       for(int j = 0; j < jmaxb; j++)
//       {
//         INFO("i = " << i << " j = " << j);
//         REQUIRE(fabs(F_cal[i][j] - F_ref[i][j]) < tol);
//         REQUIRE(fabs(G_cal[i][j] - G_ref[i][j]) < tol);
        
//       }
//     }

//     #pragma endregion
// }

// TEST_CASE( "Test temperature values Fluid Trap", "[calculate_temp_FT]" )
// {
//   // Dummy Grid Parameters
//     #pragma region 
//     double xlength = 2;
//     double ylength = 1;

//     int imax = 98;
//     int jmax = 48;

//     double dt = 0.500000000 ;
//     double t_end = 2000.0;
//     double tau = 0.5;

//     double dt_Value = 10.0;

//     int itermax = 1000;
//     double eps = 0.00001;
//     double omg = 1.7;
//     double alpha = 0.5;
//     double beta = 0.00063;

//     double Re = 10000;
//     double PR = 7;
//     double GX = 0.0;
//     double GY = -9.81;

//     double UI = 0.0;
//     double VI = 0.0;
//     double PI = 0.0;
//     double TI = 0.0;

//     double UT = 0.0;

//     double dx = 0.02;
//     double dy = 0.02;

//     double tol = 0.0000001;
//     int** geometry;
//     #pragma endregion
//   //Getting log file names
//     #pragma region
//     // command line arrays
//     int argc = 3;
//     const char* argv[] = { "./sim","-ft","-t"};
//     // File paths string 
//     std::string szFileName;
//     std::string geoFileName;
//     std::string problemName;
//     std::vector<std::string> logFilesName;
//     // Problem Flags
//     bool temp_flag;

//     logFilesName  = ParseTestingProblemFlags(argc,argv,szFileName,geoFileName,problemName,temp_flag);
//     #pragma endregion

//     // Geometry Initialization
//     #pragma region 
//     geometry = read_pgm((geoFileName).c_str());
//     #pragma endregion
//   // Initialize Dummy Grid
//   #pragma region 
//   Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry);
//   int imaxb = domain.imaxb();
//   int jmaxb = domain.jmaxb();
//   #pragma endregion
  
//   // Data Structure Declaration
//   #pragma region 
//   // Input Data Structures
//   matrix<double> U;
//   matrix<double> V;
//   matrix<double> T;
//   // Calculation Data Structures
//   matrix<double> T_cal;
//   // Reference Data Structures
//   matrix<double> T_ref;
//   #pragma endregion
  
//   // Data Structure Allocation
//   #pragma region 
//     // Input Data Structures
//     U.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     V.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     T.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     // Allocating memory for matrices
//     T_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    
//     // Allocating memory for matrices
//     T_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    
//     #pragma endregion
 
//   // Reading log files in data structures
//   #pragma region 
//   // Reading dt
//     std::ifstream fin(logFilesName[static_cast<int>(log_file::log_dt)]);
//     fin >> dt;
//     // Input Data Structures  
//     read_matrix(logFilesName[static_cast<int>(log_file::log_U_Input)],U,imaxb,jmaxb);
//     read_matrix(logFilesName[static_cast<int>(log_file::log_V_Input)],V,imaxb,jmaxb);
//     read_matrix(logFilesName[static_cast<int>(log_file::log_T_Input)],T,imaxb,jmaxb);
//     // Reference Data Structures
//     read_matrix(logFilesName[static_cast<int>(log_file::log_T_Output)],T_ref,imaxb,jmaxb);
//   #pragma endregion
  
//   // Dummy Scenario
//   #pragma region  
//   domain.set_velocity(U, velocity_type::U);
//   domain.set_velocity(V, velocity_type::V);
//   domain.set_temperature(T);
//   calculate_temp(dt,dx,dy,imax,jmax,Re,PR,alpha,0,domain);
//   #pragma endregion
  
//   // Write Output Data
//   #pragma region 
//   domain.temperature(T_cal);
//   write_matrix(cmp_log_T_t, T_cal, imaxb, jmaxb);
//   #pragma endregion
  
//   // Comparing Reference matrices to Calculated matrices
//   #pragma region 
    

//     for(int i = 0; i < imaxb; i++)
//     {
//       for(int j = 0; j < jmaxb; j++)
//       {
//         INFO("i = " << i << " j = " << j);
//         REQUIRE(fabs(T_cal[i][j] - T_ref[i][j]) < tol);
        
//       }
//     }

//     #pragma endregion
// }

// TEST_CASE( "Test rhs of PPE matrix Fluid Trap", "[calculate_RS_FT]" )
// {
//   // Dummy Grid Parameters
//     #pragma region 
//     double xlength = 2;
//     double ylength = 1;

//     int imax = 98;
//     int jmax = 48;

//     double dt = 0.500000000 ;
//     double t_end = 2000.0;
//     double tau = 0.5;

//     double dt_Value = 10.0;

//     int itermax = 1000;
//     double eps = 0.00001;
//     double omg = 1.7;
//     double alpha = 0.5;
//     double beta = 0.00063;

//     double Re = 10000;
//     double PR = 7;
//     double GX = 0.0;
//     double GY = -9.81;

//     double UI = 0.0;
//     double VI = 0.0;
//     double PI = 0.0;
//     double TI = 0.0;

//     double UT = 0.0;

//     double dx = 0.02;
//     double dy = 0.02;

//     double tol = 0.0000001;
//     int** geometry;
//     #pragma endregion
//     //Getting log file names
//     #pragma region
//     // command line arrays
//     int argc = 3;
//     const char* argv[] = { "./sim","-ft","-t"};
//     // File paths string 
//     std::string szFileName;
//     std::string geoFileName;
//     std::string problemName;
//     std::vector<std::string> logFilesName;
//     // Problem Flags
//     bool temp_flag;

//     logFilesName  = ParseTestingProblemFlags(argc,argv,szFileName,geoFileName,problemName,temp_flag);
//     #pragma endregion

//     // Geometry Initialization
//     #pragma region 
//     geometry = read_pgm((geoFileName).c_str());
//     #pragma endregion
    
//     // Initialize Dummy Grid
//     #pragma region 
//     Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry);

//     int imaxb = domain.imaxb();
//     int jmaxb = domain.jmaxb();
//     #pragma endregion

    
//     // Data Structure Declaration
//     #pragma region 
//     // Input Data Structures
//     matrix<double> F;
//     matrix<double> G;
//     // Calculation Data Structures
//     matrix<double> RS_cal;
//     // Reference Data Structures
//     matrix<double> RS_ref;
//     #pragma endregion

//     // Data Structure Allocation
//     #pragma region 
//     // Input Data Structures
//     F.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     G.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     // Allocating memory for matrices
//     RS_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     // Allocating memory for matrices
//     RS_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     #pragma endregion

//     // Reading log files in data structures
//     #pragma region 
//     // Redaing dt
//     std::ifstream fin(logFilesName[static_cast<int>(log_file::log_dt)]);
//     fin >> dt;    
//     // Input Data Structures
//     read_matrix(logFilesName[static_cast<int>(log_file::log_F)],F,imaxb,jmaxb);
//     read_matrix(logFilesName[static_cast<int>(log_file::log_G)],G,imaxb,jmaxb);

//     // Reference Data
//     read_matrix(logFilesName[static_cast<int>(log_file::log_RS)],RS_ref,imaxb,jmaxb);

//     #pragma endregion

//     // Dummy Scenario
//     #pragma region
    
    
    
//     calculate_rs(dt, dx, dy, imax, jmax, F, G, RS_cal,domain);

//     #pragma endregion
//     // Writing calculated matrices to log files
//     #pragma region 
//     write_matrix(cmp_log_RS_t, RS_cal, imaxb, jmaxb);
    
//     #pragma endregion
//     // Comparing Reference matrices to Calculated matrices
//     #pragma region 
//     for(int i = 0; i < imaxb; i++)
//     {
//       for(int j = 0; j < jmaxb; j++)
//       {
//         if(RS_ref[i][j] == 0 )
//         {
//           // std::cout << "i = " << i << " j = " << j << std::endl;
//             INFO("i = " << i << " j = " << j);
//             REQUIRE(RS_cal[i][j] == RS_ref[i][j]);
//         }
//       }
//     }

//     for(int i = 0; i < imaxb; i++)
//     {
//       for(int j = 0; j < jmaxb; j++)
//       {
//         INFO("i = " << i << " j = " << j);
//         REQUIRE(fabs(RS_cal[i][j] - RS_ref[i][j]) < tol);
        
//       }
//     }

//     #pragma endregion
// }

// TEST_CASE( "Test Test SOR Fluid Trap", "[sor_FT]" )
// {
//   // Dummy Grid Parameters
//     #pragma region 
//     double xlength = 2;
//     double ylength = 1;

//     int imax = 98;
//     int jmax = 48;

//     double dt = 0.500000000 ;
//     double t_end = 2000.0;
//     double tau = 0.5;

//     double dt_Value = 10.0;

//     int itermax = 1000;
//     double eps = 0.00001;
//     double omg = 1.7;
//     double alpha = 0.5;
//     double beta = 0.00063;

//     double Re = 10000;
//     double PR = 7;
//     double GX = 0.0;
//     double GY = -9.81;

//     double UI = 0.0;
//     double VI = 0.0;
//     double PI = 0.0;
//     double TI = 0.0;

//     double UT = 0.0;

//     double dx = 0.02;
//     double dy = 0.02;

//     double tol = 0.0000001;
//     int** geometry;
//     #pragma endregion
//     //Getting log file names
//     #pragma region
//     // command line arrays
//     int argc = 3;
//     const char* argv[] = { "./sim","-ft","-t"};
//     // File paths string 
//     std::string szFileName;
//     std::string geoFileName;
//     std::string problemName;
//     std::vector<std::string> logFilesName;
//     // Problem Flags
//     bool temp_flag;

//     logFilesName  = ParseTestingProblemFlags(argc,argv,szFileName,geoFileName,problemName,temp_flag);
//     #pragma endregion

//     // Geometry Initialization
//     #pragma region 
//     geometry = read_pgm((geoFileName).c_str());
//     #pragma endregion
//     // Initialize Dummy Grid
//     #pragma region 
//     Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry);

//     int imaxb = domain.imaxb();
//     int jmaxb = domain.jmaxb();
//     #pragma endregion

    
//     // Data Structure Declaration
//     #pragma region 
//     // Input Data Structures
//     matrix<double> U;
//     matrix<double> V;
//     matrix<double> RS;
//     // Calculation Data Structures
//     matrix<double> P_cal;
//     // Reference Data Structures
//     matrix<double> P_ref;
//     #pragma endregion

//     // Data Structure Allocation
//     #pragma region 
//     // Input Data Structures
//     U.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     V.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     RS.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     // Allocating memory for matrices
//     P_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     // Allocating memory for matrices
//     P_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     #pragma endregion

//     // Reading log files in data structures
//     #pragma region 
//     // Reading dt
//     std::ifstream fin(logFilesName[static_cast<int>(log_file::log_dt)]);
//     fin >> dt; 
//     //Input Data    
//     read_matrix(logFilesName[static_cast<int>(log_file::log_U_Input)],U,imaxb,jmaxb);
//     read_matrix(logFilesName[static_cast<int>(log_file::log_V_Input)],V,imaxb,jmaxb);
//     read_matrix(logFilesName[static_cast<int>(log_file::log_RS)],RS,imaxb,jmaxb);
//     read_matrix(logFilesName[static_cast<int>(log_file::log_P_Input)],P_cal,imaxb,jmaxb);

//     // Reference Data
//     read_matrix(logFilesName[static_cast<int>(log_file::log_P_Output)],P_ref,imaxb,jmaxb);
//     #pragma endregion

//     // Dummy Scenario
//     #pragma region
//     domain.set_velocity(U,velocity_type::U);
//     domain.set_velocity(V,velocity_type::V);
//     domain.set_pressure(P_cal);
//     int it = 0;
//     double res = 1.0;
//     // Solve system using SOR
//     while(it < itermax && res > eps){
//       sor(omg,dx,dy,imax,jmax,domain,RS,&res);
//       it++;
//     }
    

//     #pragma endregion
   
//     // Writing calculated matrices to log files
//     #pragma region 
//     domain.pressure(P_cal);
//     write_matrix(cmp_log_P_t, P_cal, imaxb, jmaxb);
    
//     #pragma endregion
   
//     // Comparing Reference matrices to Calculated matrices
//     #pragma region 
   
//     for(int i = 0; i < imaxb; i++)
//     {
//       for(int j = 0; j < jmaxb; j++)
//       {
//         INFO("i = " << i << " j = " << j);
//         REQUIRE(fabs(P_cal[i][j] - P_ref[i][j]) < tol);
        
//       }
//     }

//     #pragma endregion
// }

// TEST_CASE( "Test velocity (U and v) values Fluid Trap", "[calculate_uv_FT]" )
// {
//   // Dummy Grid Parameters
//     #pragma region 
//     double xlength = 2;
//     double ylength = 1;

//     int imax = 98;
//     int jmax = 48;

//     double dt = 0.500000000 ;
//     double t_end = 2000.0;
//     double tau = 0.5;

//     double dt_Value = 10.0;

//     int itermax = 1000;
//     double eps = 0.00001;
//     double omg = 1.7;
//     double alpha = 0.5;
//     double beta = 0.00063;

//     double Re = 10000;
//     double PR = 7;
//     double GX = 0.0;
//     double GY = -9.81;

//     double UI = 0.0;
//     double VI = 0.0;
//     double PI = 0.0;
//     double TI = 0.0;

//     double UT = 0.0;

//     double dx = 0.02;
//     double dy = 0.02;

//     double tol = 0.0000001;
//     int** geometry;
//     #pragma endregion
//   //Getting log file names
//     #pragma region
//     // command line arrays
//     int argc = 3;
//     const char* argv[] = { "./sim","-ft","-t"};
//     // File paths string 
//     std::string szFileName;
//     std::string geoFileName;
//     std::string problemName;
//     std::vector<std::string> logFilesName;
//     // Problem Flags
//     bool temp_flag;

//     logFilesName  = ParseTestingProblemFlags(argc,argv,szFileName,geoFileName,problemName,temp_flag);
//     #pragma endregion

//     // Geometry Initialization
//     #pragma region 
//     geometry = read_pgm((geoFileName).c_str());
//     #pragma endregion
//   // Initialize Dummy Grid
//   #pragma region 
//   Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry);
//   int imaxb = domain.imaxb();
//   int jmaxb = domain.jmaxb();
//   #pragma endregion
  
//   // Data Structure Declaration
//     #pragma region 
//     // Input Data Structures
//     matrix<double> F;
//     matrix<double> G;
//     matrix<double> P;
//     // Calculation Data Structures
//     matrix<double> U_cal;
//     matrix<double> V_cal;
//     // Reference Data Structures
//     matrix<double> U_ref;
//     matrix<double> V_ref;
//     #pragma endregion

//     // Data Structure Allocation
//     #pragma region 
//     // Input Data Structures
//     F.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     G.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     P.resize(imaxb,std::vector<double>(jmaxb,0.0));    
//     // Calculation Data Structures
//     U_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     V_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     // Allocating memory for matrices
//     U_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     V_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     #pragma endregion

//     // Reading log files in data structures
//     #pragma region 
//     // Reading dt
//     std::ifstream fin(logFilesName[static_cast<int>(log_file::log_dt)]);
//     fin >> dt;
//     // Input Data
//     read_matrix(logFilesName[static_cast<int>(log_file::log_F)],F,imaxb,jmaxb);
//     read_matrix(logFilesName[static_cast<int>(log_file::log_G)],G,imaxb,jmaxb);
//     read_matrix(logFilesName[static_cast<int>(log_file::log_P_Output)],P,imaxb,jmaxb);
//     read_matrix(logFilesName[static_cast<int>(log_file::log_U_Input)],U_cal,imaxb,jmaxb);
//     read_matrix(logFilesName[static_cast<int>(log_file::log_V_Input)],V_cal,imaxb,jmaxb);

//     // Reference Data Structures
//     read_matrix(logFilesName[static_cast<int>(log_file::log_U_Output)],U_ref,imaxb,jmaxb);
//     read_matrix(logFilesName[static_cast<int>(log_file::log_V_Output)],V_ref,imaxb,jmaxb);
 
//     #pragma endregion

//     // Dummy Scenario
//     #pragma region
//     domain.set_velocity(U_cal,velocity_type::U);
//     domain.set_velocity(V_cal,velocity_type::V);
//     domain.set_pressure(P);
//     calculate_uv(dt, dx, dy, imax, jmax, domain, F, G);

//     #pragma endregion
    
//     // Writing calculated matrices to log files
//     #pragma region 
//     domain.velocity(U_cal,velocity_type::U);
//     domain.velocity(V_cal,velocity_type::V);
//     write_matrix(cmp_log_U_t, U_cal , imaxb, jmaxb);
//     write_matrix(cmp_log_V_t, V_cal , imaxb, jmaxb);
//     #pragma endregion
    
//     // Comparing Reference matrices to Calculated matrices
//     #pragma region 
//     for(int i = 0; i < imaxb; i++)
//     {
//       for(int j = 0; j < jmaxb; j++)
//       {
//         if(U_ref[i][j] == 0 )
//         {
//           // std::cout << "i = " << i << " j = " << j << std::endl;
//             INFO("i = " << i << " j = " << j);
//             REQUIRE(U_cal[i][j] == U_ref[i][j]);
//         }
//         if(V_ref[i][j] == 0 )
//         {
//           // std::cout << "i = " << i << " j = " << j << std::endl;
//             INFO("i = " << i << " j = " << j);
//             REQUIRE(V_cal[i][j] == V_ref[i][j]);
//         }
        
//       }
//     }

//     for(int i = 1; i < imaxb-1; i++)
//     {
//       for(int j = 1; j < jmaxb-1; j++)
//       {
//         INFO("i = " << i << " j = " << j);
//         REQUIRE(fabs(U_cal[i][j] - U_ref[i][j]) < tol);
//         REQUIRE(fabs(V_cal[i][j] - V_ref[i][j]) < tol);
//       }
//     }

//     #pragma endregion
// }

#pragma endregion

// Fluid Trap Validation Tests
#pragma region 

// TEST_CASE("Test calculate dt Fluid Trap_", "[calculate_dt_FT_]")
// {
//     // Dummy Grid Parameters
//     #pragma region 
//     double xlength = 2;
//     double ylength = 1;

//     int imax = 98;
//     int jmax = 48;

//     double dt = 0.500000000 ;
//     double t_end = 2000.0;
//     double tau = 0.5;

//     double dt_Value = 10.0;

//     int itermax = 1000;
//     double eps = 0.00001;
//     double omg = 1.7;
//     double alpha = 0.5;
//     double beta = 0.00063;

//     double Re = 10000;
//     double PR = 7;
//     double GX = 0.0;
//     double GY = -9.81;

//     double UI = 0.0;
//     double VI = 0.0;
//     double PI = 0.0;
//     double TI = 0.0;

//     double UT = 0.0;

//     double dx = 0.02;
//     double dy = 0.02;

//     double tol = 0.0000001;
//     bool temp_flag = true;
//     int** geometry;
//     #pragma endregion
    

//     // Geometry Initialization
//     #pragma region 
//     geometry = read_pgm((std::string(Test_Data_t) + "Fluid_Trap" + ".pgm").c_str());
//     #pragma endregion

//     // Initialzie Dummy Grid
//     #pragma region    
//     Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry) ;
//     int imaxb = domain.imaxb();
//     int jmaxb = domain.jmaxb();
//     #pragma endregion

//     // Data Structure Declaration
//     #pragma region 
//     matrix<double> U;
//     matrix<double> V;
//     matrix<double> P;
//     matrix<double> T;
//     matrix<double> Test;
//     #pragma endregion

  
//     // Data Structure Allocation
//     #pragma region 
//     // Input Data Structures
//     U.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     V.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     T.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     #pragma endregion
     
//     #pragma region 
//     // File and directory paths
//     std::string logFileName;

//     // Input Data Structures
//     logFileName = Ref_Test_Data_t + "Tlog_U_Input" + "1";
//     read_matrix(logFileName.c_str(),U,imaxb,jmaxb);
//     logFileName = Ref_Test_Data_t + "Tlog_V_Input" + "1";
//     read_matrix(logFileName.c_str(),V,imaxb,jmaxb);
//     logFileName = Ref_Test_Data_t + "Tlog_T_Input" + "1";
//     read_matrix(logFileName.c_str(),T,imaxb,jmaxb);

//     #pragma endregion
   

//     /* Real Scenario*/
//     #pragma region 
  
//     domain.set_velocity(U, velocity_type::U);
//     domain.set_velocity(V, velocity_type::V);
//     domain.set_temperature(T);
//     domain.velocity(Test,velocity_type::U);
    

//     calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,PR,&temp_flag,domain);

    
//     #pragma endregion

//     #pragma region 
//     REQUIRE(dt == 0.5);
//     #pragma endregion
// }

// TEST_CASE( "Test F and G matrix Fluid Trap_", "[calculate_fg_FT_]" )
// {
//   // Dummy Grid Paramters
//     #pragma region 
    
//     double xlength = 2;
//     double ylength = 1;

//     int imax = 98;
//     int jmax = 48;

//     double dt = 0.500000000 ;
//     double t_end = 2000.0;
//     double tau = 0.5;

//     double dt_Value = 10.0;

//     int itermax = 1000;
//     double eps = 0.00001;
//     double omg = 1.7;
//     double alpha = 0.5;
//     double beta = 0.00063;

//     double Re = 10000;
//     double PR = 7;
//     double GX = 0.0;
//     double GY = -9.81;

//     double UI = 0.0;
//     double VI = 0.0;
//     double PI = 0.0;
//     double TI = 0.0;

//     double UT = 0.0;

//     double dx = 0.02;
//     double dy = 0.02;

//     double tol = 0.0000001;
//     int** geometry;
//     #pragma endregion

//     // Geometry Initialization
//     #pragma region 
//     geometry = read_pgm((std::string(Test_Data_t) + "Fluid_Trap" + ".pgm").c_str());
//     #pragma endregion

//     // Initialize Dummy Grid
//     #pragma region 
//     Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry);

//     int imaxb = domain.imaxb();
//     int jmaxb = domain.jmaxb();
//     #pragma endregion

//     #pragma region 
//     // Input Data Structures
//     matrix<double> U;
//     matrix<double> V;
//     matrix<double> T;
//     // Calculation Data Structures
//     matrix<double> F_cal;
//     matrix<double> G_cal;
//     // Reference Data Structures
//     matrix<double> F_ref;
//     matrix<double> G_ref;
//     #pragma endregion

//     // Data Structure Allocation
//     #pragma region 
//     // Input Data Structures
//     U.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     V.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     T.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     // Allocating memory for matrices
//     F_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     G_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    
//     // Allocating memory for matrices
//     F_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     G_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    
//     #pragma endregion

//     // Reading log files in data structures
//     #pragma region 
//     // File and directory paths
//     std::string logFileName;

//     // Input Data Structures
//     logFileName = Ref_Test_Data_t + "Tlog_U_Input" + "1";
//     read_matrix(logFileName.c_str(),U,imaxb,jmaxb);
//     logFileName = Ref_Test_Data_t + "Tlog_V_Input" + "1";
//     read_matrix(logFileName.c_str(),V,imaxb,jmaxb);
//     logFileName = Ref_Test_Data_t + "Tlog_T" + "1";
//     read_matrix(logFileName.c_str(),T,imaxb,jmaxb);

//     // Reference Data Structures
//     logFileName = Ref_Test_Data_t + "Tlog_F" + "1";
//     read_matrix(logFileName.c_str(),F_ref,imaxb,jmaxb);
//     logFileName = Ref_Test_Data_t + "Tlog_G" + "1";
//     read_matrix(logFileName.c_str(),G_ref,imaxb,jmaxb);
 
//     #pragma endregion

//     // Dummy Scenario
//     #pragma region
    
//     domain.set_velocity(U, velocity_type::U);
//     domain.set_velocity(V, velocity_type::V);
//     domain.set_temperature(T);

//     calculate_fg_temp(Re,GX,GY,alpha,dt,dx,dy,imax,jmax,beta,domain,F_cal,G_cal);

//     #pragma endregion

//      #pragma region 
//     write_matrix(cmp_log_F_t, F_cal, imaxb, jmaxb);
//     write_matrix(cmp_log_G_t, G_cal, imaxb, jmaxb);
//     #pragma endregion

//     // Comparing Reference matrices to Calculated matrices
//     #pragma region 
    

//     for(int i = 0; i < imaxb; i++)
//     {
//       for(int j = 0; j < jmaxb; j++)
//       {
//         INFO("i = " << i << " j = " << j);
//         REQUIRE(fabs(F_cal[i][j] - F_ref[i][j]) < tol);
//         REQUIRE(fabs(G_cal[i][j] - G_ref[i][j]) < tol);
        
//       }
//     }

//     #pragma endregion
// }

// TEST_CASE( "Test temperature values Fluid Trap_", "[calculate_temp_FT_]" )
// {
//   // Dummy Grid Paramters
//   #pragma region 
    
//     double xlength = 2;
//     double ylength = 1;

//     int imax = 98;
//     int jmax = 48;

//     double dt = 0.500000000 ;
//     double t_end = 2000.0;
//     double tau = 0.5;

//     double dt_Value = 10.0;

//     int itermax = 1000;
//     double eps = 0.00001;
//     double omg = 1.7;
//     double alpha = 0.5;
//     double beta = 0.00063;

//     double Re = 10000;
//     double PR = 7;
//     double GX = 0.0;
//     double GY = -9.81;

//     double UI = 0.0;
//     double VI = 0.0;
//     double PI = 0.0;
//     double TI = 0.0;

//     double UT = 0.0;

//     double dx = 0.02;
//     double dy = 0.02;

//     double tol = 0.0000001;
//     int** geometry;
//     #pragma endregion
  
//   // Geometry Initialization
//   #pragma region 
//     geometry = read_pgm((std::string(Test_Data_t) + "Fluid_Trap"  + ".pgm").c_str());
//     #pragma endregion
  
//   // Initialize Dummy Grid
//   #pragma region 
//   Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry);
//   int imaxb = domain.imaxb();
//   int jmaxb = domain.jmaxb();
//   #pragma endregion
  
//   // Data Structure Declaration
//   #pragma region 
//   // Input Data Structures
//   matrix<double> U;
//   matrix<double> V;
//   matrix<double> T;
//   // Calculation Data Structures
//   matrix<double> T_cal;
//   // Reference Data Structures
//   matrix<double> T_ref;
//   #pragma endregion
  
//   // Data Structure Allocation
//   #pragma region 
//     // Input Data Structures
//     U.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     V.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     T.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     // Allocating memory for matrices
//     T_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    
//     // Allocating memory for matrices
//     T_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    
//     #pragma endregion
 
//   // Reading log files in data structures
//   #pragma region 
//   // File and directory paths
//   std::string logFileName;
//   // Input Data Structures
//   logFileName = Ref_Test_Data_t + "Tlog_U_Input" + "1";
//   read_matrix(logFileName.c_str(),U,imaxb,jmaxb);
//   logFileName = Ref_Test_Data_t + "Tlog_V_Input" + "1";
//   read_matrix(logFileName.c_str(),V,imaxb,jmaxb);
//   logFileName = Ref_Test_Data_t + "Tlog_T_Input" + "1";
//   read_matrix(logFileName.c_str(),T,imaxb,jmaxb);
//   std::cout << "Flag of top wall" << domain.cell(1,jmax+1).flag() << std::endl;
//   // Reference Data Structures
//   logFileName = Ref_Test_Data_t + "Tlog_T" + "1";
//   read_matrix(logFileName.c_str(),T_ref,imaxb,jmaxb);
//   #pragma endregion
  
//   // Dummy Scenario
//   #pragma region  
//   domain.set_velocity(U, velocity_type::U);
//   domain.set_velocity(V, velocity_type::V);
//   domain.set_temperature(T);
//   calculate_temp(dt,dx,dy,imax,jmax,Re,PR,alpha,0,domain);
//   #pragma endregion
  
//   // Write Output Data
//   #pragma region 
//   domain.temperature(T_cal);
//   write_matrix(cmp_log_T_t, T_cal, imaxb, jmaxb);
//   #pragma endregion
  
//   // Comparing Reference matrices to Calculated matrices
//   #pragma region 
    

//     for(int i = 0; i < imaxb; i++)
//     {
//       for(int j = 0; j < jmaxb; j++)
//       {
//         INFO("i = " << i << " j = " << j);
//         REQUIRE(fabs(T_cal[i][j] - T_ref[i][j]) < tol);
        
//       }
//     }

//     #pragma endregion
// }

// TEST_CASE( "Test rhs of PPE matrix Fluid Trap_", "[calculate_RS_FT_]" )
// {
//   // Dummy Grid Paramters
//     #pragma region 
    
//     double xlength = 2;
//     double ylength = 1;

//     int imax = 98;
//     int jmax = 48;

//     double dt = 0.500000000 ;
//     double t_end = 2000.0;
//     double tau = 0.5;

//     double dt_Value = 10.0;

//     int itermax = 1000;
//     double eps = 0.00001;
//     double omg = 1.7;
//     double alpha = 0.5;
//     double beta = 0.00063;

//     double Re = 10000;
//     double PR = 7;
//     double GX = 0.0;
//     double GY = -9.81;

//     double UI = 0.0;
//     double VI = 0.0;
//     double PI = 0.0;
//     double TI = 0.0;

//     double UT = 0.0;

//     double dx = 0.02;
//     double dy = 0.02;

//     double tol = 0.0000001;
//     int** geometry;
//     #pragma endregion

//     // Geometry Initialization
//     #pragma region 
//     geometry = read_pgm((std::string(Test_Data_t) + "Fluid_Trap"  + ".pgm").c_str());
//     #pragma endregion

//     // Initialize Dummy Grid
//     #pragma region 
//     Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry);

//     int imaxb = domain.imaxb();
//     int jmaxb = domain.jmaxb();
//     #pragma endregion

    
//     // Data Structure Declaration
//     #pragma region 
//     // Input Data Structures
//     matrix<double> F;
//     matrix<double> G;
//     // Calculation Data Structures
//     matrix<double> RS_cal;
//     // Reference Data Structures
//     matrix<double> RS_ref;
//     #pragma endregion

//     // Data Structure Allocation
//     #pragma region 
//     // Input Data Structures
//     F.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     G.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     // Allocating memory for matrices
//     RS_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     // Allocating memory for matrices
//     RS_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     #pragma endregion

//     // Reading log files in data structures
//     #pragma region 
//     // File and directory paths
//     std::string logFileName;
    
//     // Input Data Structures
//     logFileName = Ref_Test_Data_t + "Tlog_F" + "1";
//     std::cout << "File name " << logFileName << std::endl;
//     read_matrix(logFileName.c_str(),F,imaxb,jmaxb);
//     logFileName = Ref_Test_Data_t + "Tlog_G" + "1";
//     read_matrix(logFileName.c_str(),G,imaxb,jmaxb);

//     // Reference Data Structures
//     logFileName = Ref_Test_Data_t + "Tlog_RS" + "1";
//     read_matrix(logFileName.c_str(),RS_ref,imaxb,jmaxb);

//     #pragma endregion

//     // Dummy Scenario
//     #pragma region
    
    
    
//     calculate_rs(dt, dx, dy, imax, jmax, F, G, RS_cal,domain);

//     #pragma endregion
//     // Writing calculated matrices to log files
//     #pragma region 
//     write_matrix(cmp_log_RS_t, RS_cal, imaxb, jmaxb);
    
//     #pragma endregion
//     // Comparing Reference matrices to Calculated matrices
//     #pragma region 
//     for(int i = 0; i < imaxb; i++)
//     {
//       for(int j = 0; j < jmaxb; j++)
//       {
//         if(RS_ref[i][j] == 0 )
//         {
//           // std::cout << "i = " << i << " j = " << j << std::endl;
//             INFO("i = " << i << " j = " << j);
//             REQUIRE(RS_cal[i][j] == RS_ref[i][j]);
//         }
//       }
//     }

//     for(int i = 0; i < imaxb; i++)
//     {
//       for(int j = 0; j < jmaxb; j++)
//       {
//         INFO("i = " << i << " j = " << j);
//         REQUIRE(fabs(RS_cal[i][j] - RS_ref[i][j]) < tol);
        
//       }
//     }

//     #pragma endregion
// }

// TEST_CASE( "Test Test SOR Fluid Trap_", "[sor_FT_]" )
// {
//   // Dummy Grid Paramters
//     #pragma region 
    
//     double xlength = 2;
//     double ylength = 1;

//     int imax = 98;
//     int jmax = 48;

//     double dt = 0.500000000 ;
//     double t_end = 2000.0;
//     double tau = 0.5;

//     double dt_Value = 10.0;

//     int itermax = 1000;
//     double eps = 0.00001;
//     double omg = 1.7;
//     double alpha = 0.5;
//     double beta = 0.00063;

//     double Re = 10000;
//     double PR = 7;
//     double GX = 0.0;
//     double GY = -9.81;

//     double UI = 0.0;
//     double VI = 0.0;
//     double PI = 0.0;
//     double TI = 0.0;

//     double UT = 0.0;

//     double dx = 0.02;
//     double dy = 0.02;

//     double tol = 0.00000001;
//     int** geometry;
//     #pragma endregion

//     // Geometry Initialization
//     #pragma region 
//     geometry = read_pgm((std::string(Test_Data_t) + "Fluid_Trap"  + ".pgm").c_str());
//     #pragma endregion

//     // Initialize Dummy Grid
//     #pragma region 
//     Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry);

//     int imaxb = domain.imaxb();
//     int jmaxb = domain.jmaxb();
//     #pragma endregion

    
//     // Data Structure Declaration
//     #pragma region 
//     // Input Data Structures
//     matrix<double> U;
//     matrix<double> V;
//     matrix<double> RS;
//     // Calculation Data Structures
//     matrix<double> P_cal;
//     // Reference Data Structures
//     matrix<double> P_ref;
//     #pragma endregion

//     // Data Structure Allocation
//     #pragma region 
//     // Input Data Structures
//     U.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     V.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     RS.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     // Allocating memory for matrices
//     P_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     // Allocating memory for matrices
//     P_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     #pragma endregion

//     // Reading log files in data structures
//     #pragma region 
//     // File and directory paths
//     std::string logFileName;
    
//     // Input Data Structures
//     logFileName = Ref_Test_Data_t + "Tlog_U_Input" + "1";
//     std::cout << "File name " << logFileName << std::endl;
//     read_matrix(logFileName.c_str(),U,imaxb,jmaxb);
//     logFileName = Ref_Test_Data_t + "Tlog_V_Input" + "1";
//     read_matrix(logFileName.c_str(),V,imaxb,jmaxb);
//     logFileName = Ref_Test_Data_t + "Tlog_RS" + "1";
//     read_matrix(logFileName.c_str(),RS,imaxb,jmaxb);
//     logFileName = Ref_Test_Data_t + "Tlog_P_Output" + "01";
//     read_matrix(logFileName.c_str(),P_cal,imaxb,jmaxb);
    

//     // Reference Data
//     logFileName = Ref_Test_Data_t + "Tlog_P_Output" + "11";
//     read_matrix(logFileName.c_str(),P_ref,imaxb,jmaxb);

//     #pragma endregion

//     // Dummy Scenario
//     #pragma region
//     domain.set_velocity(U,velocity_type::U);
//     domain.set_velocity(V,velocity_type::V);
//     domain.set_pressure(P_cal);
//     int it = 0;
//     double res = 1.0;
//     // Solve system using SOR
//     // while(it < itermax && res > eps){
//       std::cout << "res = " << res << std::endl;
//       sor(omg,dx,dy,imax,jmax,domain,RS,&res);
//       it++;
//     // }
    

//     #pragma endregion
   
//     // Writing calculated matrices to log files
//     #pragma region 
//     domain.pressure(P_cal);
//     std::cout << "P" << domain.cell(38,30).pressure() << std::endl;
//     write_matrix(cmp_log_P_t, P_cal, imaxb, jmaxb);
    
//     #pragma endregion
   
//     // Comparing Reference matrices to Calculated matrices
//     #pragma region 
   
//     for(int i = 0; i < imaxb; i++)
//     {
//       for(int j = 0; j < jmaxb; j++)
//       {
//         if(domain.cell(i,j).flag()&(1<<1))
//         {
//           // std::cout << "i = " << i << " j = " << j << std::endl;
//           INFO("i = " << i << " j = " << j);
//           INFO("flag: " << domain.cell(i,j).flag());
//           INFO("P_cal = " << P_cal[i][j] << " P_ref = " << P_ref[i][j]);
//           REQUIRE(fabs(P_cal[i][j] - P_ref[i][j]) < tol);
//         }
        
        
//       }
//     }

//     #pragma endregion
// }

// TEST_CASE( "Test velocity (U and v) values Fluid Trap_", "[calculate_uv_FT_]" )
// {
//   // Dummy Grid Paramters
//   #pragma region 
    
//     double xlength = 2;
//     double ylength = 1;

//     int imax = 98;
//     int jmax = 48;

//     double dt = 0.500000000 ;
//     double t_end = 2000.0;
//     double tau = 0.5;

//     double dt_Value = 10.0;

//     int itermax = 1000;
//     double eps = 0.00001;
//     double omg = 1.7;
//     double alpha = 0.5;
//     double beta = 0.00063;

//     double Re = 10000;
//     double PR = 7;
//     double GX = 0.0;
//     double GY = -9.81;

//     double UI = 0.0;
//     double VI = 0.0;
//     double PI = 0.0;
//     double TI = 0.0;

//     double UT = 0.0;

//     double dx = 0.02;
//     double dy = 0.02;

//     double tol = 0.0000001;
//     int** geometry;
//     #pragma endregion
  
//   // Geometry Initialization
//   #pragma region 
//   geometry = read_pgm((std::string(Test_Data_t) + "Fluid_Trap"  + ".pgm").c_str());
//   #pragma endregion
  
//   // Initialize Dummy Grid
//   #pragma region 
//   Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry);
//   int imaxb = domain.imaxb();
//   int jmaxb = domain.jmaxb();
//   #pragma endregion
  
//   // Data Structure Declaration
//     #pragma region 
//     // Input Data Structures
//     matrix<double> F;
//     matrix<double> G;
//     matrix<double> P;
//     // Calculation Data Structures
//     matrix<double> U_cal;
//     matrix<double> V_cal;
//     // Reference Data Structures
//     matrix<double> U_ref;
//     matrix<double> V_ref;
//     #pragma endregion

//     // Data Structure Allocation
//     #pragma region 
//     // Input Data Structures
//     F.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     G.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     P.resize(imaxb,std::vector<double>(jmaxb,0.0));    
//     // Calculation Data Structures
//     U_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     V_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     // Allocating memory for matrices
//     U_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     V_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     #pragma endregion

//     // Reading log files in data structures
//     #pragma region 
//     // File and directory paths
//     std::string logFileName;

//     // Input Data Structures
//     logFileName = Ref_Test_Data_t + "Tlog_F" + "1";
//     read_matrix(logFileName.c_str(),F,imaxb,jmaxb);
//     logFileName = Ref_Test_Data_t + "Tlog_G" + "1";
//     read_matrix(logFileName.c_str(),G,imaxb,jmaxb);
//     logFileName = Ref_Test_Data_t + "Tlog_P_Output" + "1";
//     read_matrix(logFileName.c_str(),P,imaxb,jmaxb);
//     logFileName = Ref_Test_Data_t + "Tlog_U_Input" + "1";
//     read_matrix(logFileName.c_str(),U_cal,imaxb,jmaxb);
//     logFileName = Ref_Test_Data_t + "Tlog_V_Input" + "1";
//     read_matrix(logFileName.c_str(),V_cal,imaxb,jmaxb);



//     // Reference Data Structures
//     logFileName = Ref_Test_Data_t + "Tlog_U_Output" + "1";
//     read_matrix(logFileName.c_str(),U_ref,imaxb,jmaxb);
//     logFileName = Ref_Test_Data_t + "Tlog_V_Output" + "1";
//     read_matrix(logFileName.c_str(),V_ref,imaxb,jmaxb);
 
//     #pragma endregion

//     // Dummy Scenario
//     #pragma region
//     domain.set_velocity(U_cal,velocity_type::U);
//     domain.set_velocity(V_cal,velocity_type::V);
//     domain.set_pressure(P);
//     calculate_uv(dt, dx, dy, imax, jmax, domain, F, G);

//     #pragma endregion
    
//     // Writing calculated matrices to log files
//     #pragma region 
//     domain.velocity(U_cal,velocity_type::U);
//     domain.velocity(V_cal,velocity_type::V);
//     write_matrix(cmp_log_U_t, U_cal , imaxb, jmaxb);
//     write_matrix(cmp_log_V_t, V_cal , imaxb, jmaxb);
//     #pragma endregion
    
//     // Comparing Reference matrices to Calculated matrices
//     #pragma region 
//     for(int i = 0; i < imaxb; i++)
//     {
//       for(int j = 0; j < jmaxb; j++)
//       {
//         if(U_ref[i][j] == 0 )
//         {
//           // std::cout << "i = " << i << " j = " << j << std::endl;
//             INFO("i = " << i << " j = " << j);
//             REQUIRE(U_cal[i][j] == U_ref[i][j]);
//         }
//         if(V_ref[i][j] == 0 )
//         {
//           // std::cout << "i = " << i << " j = " << j << std::endl;
//             INFO("i = " << i << " j = " << j);
//             REQUIRE(V_cal[i][j] == V_ref[i][j]);
//         }
        
//       }
//     }

//     for(int i = 1; i < imaxb-1; i++)
//     {
//       for(int j = 1; j < jmaxb-1; j++)
//       {
//         INFO("i = " << i << " j = " << j);
//         REQUIRE(fabs(U_cal[i][j] - U_ref[i][j]) < tol);
//         REQUIRE(fabs(V_cal[i][j] - V_ref[i][j]) < tol);
//       }
//     }

//     #pragma endregion
// }

#pragma endregion

/*----------------------Integration Tests---------------------------- */

// Karman Vortex Integration Tests


TEST_CASE("Test calculate dt & F and G matrix Karman Vortex", "[calculate_dt_fg_KV]")
{
    // Dummy Grid Parameters
    #pragma region 
    double xlength = 10;
    double ylength = 2;

    int imax = 98;
    int jmax = 18;

    double dt = 0.05;
    double dt_ref = 0.0;
    double t_end = 20.0;
    double tau = 0.5;

    double dt_Value = 2.0;

    int itermax = 500;
    double eps = 0.001;
    double omg = 1.7;
    double alpha = 0.9;

    double Re = 10000;
    double PR = 1;
    double GX = 0.0;
    double GY = 0.0;

    double UI = 1.0;
    double VI = 0.0;
    double PI = 0.0;
    double TI = 0.0;

    double UT = 0.0;

    double dx = 0.1;
    double dy = 0.1;

    double tol = 0.00000001;
    int** geometry;
    #pragma endregion
    
    //Getting log file names
    #pragma region
    // command line arrays
    int argc = 3;
    const char* argv[] = { "./sim","-k","-t"};
    // File paths string 
    std::string szFileName;
    std::string geoFileName;
    std::string problemName;
    std::vector<std::string> logFilesName;
    // Problem Flags
    bool temp_flag;

    logFilesName  = ParseTestingProblemFlags(argc,argv,szFileName,geoFileName,problemName,temp_flag);
    #pragma endregion
    // Geometry Initialization
    #pragma region 
    geometry = read_pgm((geoFileName).c_str());
    #pragma endregion

    // Initialzie Dummy Grid
    #pragma region    
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry) ;
    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();
    #pragma endregion

    // Data Structure Declaration
    #pragma region 
    // Input Data Structures
    matrix<double> U;
    matrix<double> V;
    // Calculation Data Structures
    matrix<double> F_cal;
    matrix<double> G_cal;
    // Reference Data Structures
    matrix<double> F_ref;
    matrix<double> G_ref;
    #pragma endregion
  
    // Data Structure Allocation
    #pragma region 
    // Input Data Structures
    U.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    F_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    F_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    #pragma endregion
  
    // Reading log files in data structures
    #pragma region 
    // Input Data Structures  
    read_matrix(logFilesName[static_cast<int>(log_file::log_U_Input)],U,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_V_Input)],V,imaxb,jmaxb);
    // Reference Data Structures
    read_matrix(logFilesName[static_cast<int>(log_file::log_F)],F_ref,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_G)],G_ref,imaxb,jmaxb);
    #pragma endregion
     
  
    /* Real Scenario*/
    #pragma region 
    domain.set_velocity(U, velocity_type::U);
    domain.set_velocity(V, velocity_type::V);
    
    calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,PR,&temp_flag,domain);  
    calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, domain, F_cal, G_cal);
    #pragma endregion

    // Writing calculated matrices to log files
    #pragma region 
    write_matrix(cmp_log_F_t, F_cal, imaxb, jmaxb);
    write_matrix(cmp_log_G_t, G_cal, imaxb, jmaxb);
    #pragma endregion
  
    // Comparing Reference matrices to Calculated matrices
    #pragma region 

    for(int i = 0; i < imaxb; i++)
    {
      for(int j = 0; j < jmaxb; j++)
      {
        INFO("i = " << i << " j = " << j);
        REQUIRE(fabs(F_cal[i][j] - F_ref[i][j]) < tol);
        REQUIRE(fabs(G_cal[i][j] - G_ref[i][j]) < tol);
      }
    }

    #pragma endregion
  
}

TEST_CASE("Test calculate dt & F and G and RS matrix Karman Vortex", "[calculate_dt_fg_RS_KV]")
{
    // Dummy Grid Parameters
    #pragma region 
    double xlength = 10;
    double ylength = 2;

    int imax = 98;
    int jmax = 18;

    double dt = 0.05;
    double dt_ref = 0.0;
    double t_end = 20.0;
    double tau = 0.5;

    double dt_Value = 2.0;

    int itermax = 500;
    double eps = 0.001;
    double omg = 1.7;
    double alpha = 0.9;

    double Re = 10000;
    double PR = 1;
    double GX = 0.0;
    double GY = 0.0;

    double UI = 1.0;
    double VI = 0.0;
    double PI = 0.0;
    double TI = 0.0;

    double UT = 0.0;

    double dx = 0.1;
    double dy = 0.1;

    double tol = 0.00000001;
    int** geometry;
    #pragma endregion
    
    //Getting log file names
    #pragma region
    // command line arrays
    int argc = 3;
    const char* argv[] = { "./sim","-k","-t"};
    // File paths string 
    std::string szFileName;
    std::string geoFileName;
    std::string problemName;
    std::vector<std::string> logFilesName;
    // Problem Flags
    bool temp_flag;

    logFilesName  = ParseTestingProblemFlags(argc,argv,szFileName,geoFileName,problemName,temp_flag);
    #pragma endregion
    // Geometry Initialization
    #pragma region 
    geometry = read_pgm((geoFileName).c_str());
    #pragma endregion

    // Initialzie Dummy Grid
    #pragma region    
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry) ;
    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();
    #pragma endregion

    // Data Structure Declaration
  #pragma region 
    // Input Data Structures
    matrix<double> U;
    matrix<double> V;
    // Calculation Data Structures
    matrix<double> F_cal;
    matrix<double> G_cal;
    // Calculation Data Structures
    matrix<double> RS_cal;
    // Reference Data Structures
    matrix<double> RS_ref;
    #pragma endregion
  
  // Data Structure Allocation
  #pragma region 
    // Input Data Structures
    U.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    F_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    RS_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    RS_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    #pragma endregion
  
  // Reading log files in data structures
  #pragma region 

  // Input Data Structures  
  read_matrix(logFilesName[static_cast<int>(log_file::log_U_Input)],U,imaxb,jmaxb);
  read_matrix(logFilesName[static_cast<int>(log_file::log_V_Input)],V,imaxb,jmaxb);
  // Reference Data Structures
  read_matrix(logFilesName[static_cast<int>(log_file::log_RS)],RS_ref,imaxb,jmaxb);
  #pragma endregion
     
  
   

    /* Real Scenario*/
    #pragma region 
  
    domain.set_velocity(U, velocity_type::U);
    domain.set_velocity(V, velocity_type::V);
    
    calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,PR,&temp_flag,domain);  
    calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, domain, F_cal, G_cal);
    calculate_rs(dt, dx, dy, imax, jmax, F_cal, G_cal, RS_cal,domain);
    #pragma endregion

    #pragma endregion
    // Writing calculated matrices to log files
    #pragma region 
    write_matrix(cmp_log_RS_t, RS_cal, imaxb, jmaxb);
    
    #pragma endregion
    // Comparing Reference matrices to Calculated matrices
    #pragma region 
    for(int i = 0; i < imaxb; i++)
    {
      for(int j = 0; j < jmaxb; j++)
      {
        if(RS_ref[i][j] == 0 )
        {
          
            INFO("i = " << i << " j = " << j);
            REQUIRE(RS_cal[i][j] == RS_ref[i][j]);
        }
      }
    }

    for(int i = 0; i < imaxb; i++)
    {
      for(int j = 0; j < jmaxb; j++)
      {
        INFO("i = " << i << " j = " << j);
        REQUIRE(fabs(RS_cal[i][j] - RS_ref[i][j]) < tol);
        
      }
    }

    #pragma endregion
}
 
TEST_CASE("Test calculate dt & F and G and RS matrix sor Karman Vortex", "[calculate_dt_fg_RS_sor_KV]")
{
    // Dummy Grid Parameters
    #pragma region 
    double xlength = 10;
    double ylength = 2;

    int imax = 98;
    int jmax = 18;

    double dt = 0.05;
    double dt_ref = 0.0;
    double t_end = 20.0;
    double tau = 0.5;

    double dt_Value = 2.0;

    int itermax = 500;
    double eps = 0.001;
    double omg = 1.7;
    double alpha = 0.9;

    double Re = 10000;
    double PR = 1;
    double GX = 0.0;
    double GY = 0.0;

    double UI = 1.0;
    double VI = 0.0;
    double PI = 0.0;
    double TI = 0.0;

    double UT = 0.0;

    double dx = 0.1;
    double dy = 0.1;

    double tol = 0.00000001;
    int** geometry;
    #pragma endregion
    
    //Getting log file names
    #pragma region
    // command line arrays
    int argc = 3;
    const char* argv[] = { "./sim","-k","-t"};
    // File paths string 
    std::string szFileName;
    std::string geoFileName;
    std::string problemName;
    std::vector<std::string> logFilesName;
    // Problem Flags
    bool temp_flag;

    logFilesName  = ParseTestingProblemFlags(argc,argv,szFileName,geoFileName,problemName,temp_flag);
    #pragma endregion
    // Geometry Initialization
    #pragma region 
    geometry = read_pgm((geoFileName).c_str());
    #pragma endregion

    // Initialzie Dummy Grid
    #pragma region    
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry) ;
    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();
    #pragma endregion

    // Data Structure Declaration
  #pragma region 
    // Input Data Structures
    matrix<double> U;
    matrix<double> V;
    // Calculation Data Structures
    matrix<double> F_cal;
    matrix<double> G_cal;
    // Calculation Data Structures
    matrix<double> RS_cal;
    matrix<double> P_cal;
    // Reference Data Structures
    matrix<double> P_ref;
    #pragma endregion
  
  // Data Structure Allocation
  #pragma region 
    // Input Data Structures
    U.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    F_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    RS_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    P_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    P_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    #pragma endregion
  
  // Reading log files in data structures
  #pragma region 

  // Input Data Structures  
  read_matrix(logFilesName[static_cast<int>(log_file::log_U_Input)],U,imaxb,jmaxb);
  read_matrix(logFilesName[static_cast<int>(log_file::log_V_Input)],V,imaxb,jmaxb);
  read_matrix(logFilesName[static_cast<int>(log_file::log_P_Input)],P_cal,imaxb,jmaxb);
  // Reference Data Structures
  read_matrix(logFilesName[static_cast<int>(log_file::log_P_Output)],P_ref,imaxb,jmaxb);
  #pragma endregion
     
  
   

    /* Real Scenario*/
    #pragma region 
  
    domain.set_velocity(U, velocity_type::U);
    domain.set_velocity(V, velocity_type::V);
    domain.set_pressure(P_cal);

    calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,PR,&temp_flag,domain);  
    calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, domain, F_cal, G_cal);
    calculate_rs(dt, dx, dy, imax, jmax, F_cal, G_cal, RS_cal,domain);
    
    int it = 0;
    double res = 1.0;
    // Solve system using SOR
    while(it < itermax && res > eps){
      sor(omg,dx,dy,imax,jmax,domain,RS_cal,&res);
      it++;
    }
    

    #pragma endregion

    #pragma endregion
   
    // Writing calculated matrices to log files
    #pragma region 
    domain.pressure(P_cal);
    write_matrix(cmp_log_P_t, P_cal, imaxb, jmaxb);
    
    #pragma endregion
   
    // Comparing Reference matrices to Calculated matrices
    #pragma region 
   
    for(int i = 0; i < imaxb; i++)
    {
      for(int j = 0; j < jmaxb; j++)
      {
        INFO("i = " << i << " j = " << j);
        REQUIRE(fabs(P_cal[i][j] - P_ref[i][j]) < tol);
        
      }
    }

    #pragma endregion
    
}

TEST_CASE("Test calculate dt & F and G and RS matrix sor uv Karman Vortex", "[calculate_dt_fg_RS_sor_uv_KV]")
{
    // Dummy Grid Parameters
    #pragma region 
    double xlength = 10;
    double ylength = 2;

    int imax = 98;
    int jmax = 18;

    double dt = 0.05;
    double dt_ref = 0.0;
    double t_end = 20.0;
    double tau = 0.5;

    double dt_Value = 2.0;

    int itermax = 500;
    double eps = 0.001;
    double omg = 1.7;
    double alpha = 0.9;

    double Re = 10000;
    double PR = 1;
    double GX = 0.0;
    double GY = 0.0;

    double UI = 1.0;
    double VI = 0.0;
    double PI = 0.0;
    double TI = 0.0;

    double UT = 0.0;

    double dx = 0.1;
    double dy = 0.1;

    double tol = 0.00000001;
    int** geometry;
    #pragma endregion
    
    //Getting log file names
    #pragma region
    // command line arrays
    int argc = 3;
    const char* argv[] = { "./sim","-k","-t"};
    // File paths string 
    std::string szFileName;
    std::string geoFileName;
    std::string problemName;
    std::vector<std::string> logFilesName;
    // Problem Flags
    bool temp_flag;

    logFilesName  = ParseTestingProblemFlags(argc,argv,szFileName,geoFileName,problemName,temp_flag);
    #pragma endregion
    // Geometry Initialization
    #pragma region 
    geometry = read_pgm((geoFileName).c_str());
    #pragma endregion

    // Initialzie Dummy Grid
    #pragma region    
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry) ;
    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();
    #pragma endregion

    // Data Structure Declaration
  #pragma region 
    // Input Data Structures
    matrix<double> U;
    matrix<double> V;
    // Calculation Data Structures
    matrix<double> F_cal;
    matrix<double> G_cal;
    // Calculation Data Structures
    matrix<double> RS_cal;
    matrix<double> P_cal;
    // Reference Data Structures
    matrix<double> U_ref;
    matrix<double> V_ref;
    #pragma endregion
  
  // Data Structure Allocation
  #pragma region 
    // Input Data Structures
    U.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    F_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    RS_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    P_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    U_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    #pragma endregion
  
  // Reading log files in data structures
  #pragma region 

  // Input Data Structures  
  read_matrix(logFilesName[static_cast<int>(log_file::log_U_Input)],U,imaxb,jmaxb);
  read_matrix(logFilesName[static_cast<int>(log_file::log_V_Input)],V,imaxb,jmaxb);
  read_matrix(logFilesName[static_cast<int>(log_file::log_P_Input)],P_cal,imaxb,jmaxb);
  // Reference Data Structures
    read_matrix(logFilesName[static_cast<int>(log_file::log_U_Output)],U_ref,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_V_Output)],V_ref,imaxb,jmaxb);
  #pragma endregion
     
  
   

    /* Real Scenario*/
    #pragma region 
  
    domain.set_velocity(U, velocity_type::U);
    domain.set_velocity(V, velocity_type::V);
    domain.set_pressure(P_cal);

    calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,PR,&temp_flag,domain);  
    calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, domain, F_cal, G_cal);
    calculate_rs(dt, dx, dy, imax, jmax, F_cal, G_cal, RS_cal,domain);
    
    int it = 0;
    double res = 1.0;
    // Solve system using SOR
    while(it < itermax && res > eps){
      sor(omg,dx,dy,imax,jmax,domain,RS_cal,&res);
      it++;
    }
    
    calculate_uv(dt, dx, dy, imax, jmax, domain, F_cal, G_cal);

    #pragma endregion

    #pragma endregion
   
    // Writing calculated matrices to log files
    #pragma region 
    domain.velocity(U,velocity_type::U);
    domain.velocity(V,velocity_type::V);
    write_matrix(cmp_log_U_t, U , imaxb, jmaxb);
    write_matrix(cmp_log_V_t, V , imaxb, jmaxb);
    #pragma endregion
    
    // Comparing Reference matrices to Calculated matrices
    #pragma region 
    

    for(int i = 1; i < imaxb-1; i++)
    {
      for(int j = 1; j < jmaxb-1; j++)
      {
        INFO("i = " << i << " j = " << j);
        REQUIRE(fabs(U[i][j] - U_ref[i][j]) < tol);
        REQUIRE(fabs(V[i][j] - V_ref[i][j]) < tol);
      }
    }

    #pragma endregion
    
}


// Flow Over A Step Integration Tests


TEST_CASE("Test calculate dt & F and G matrix Flow Over A Step", "[calculate_dt_fg_FS]")
{
    // Dummy Grid Paramters
    #pragma region 
    
    double xlength = 10;
    double ylength = 2;

    int imax = 98;
    int jmax = 18;

    double dt = 0.05;
    double dt_ref = 0.0;
    double t_end = 500.0;
    double tau = 0.5;

    double dt_Value = 10.0;

    int itermax = 500;
    double eps = 0.001;
    double omg = 1.7;
    double alpha = 0.9;

    double Re = 100;
    double PR = 1;
    double GX = 0.0;
    double GY = 0.0;

    double UI = 0.0;
    double VI = 0.0;
    double PI = 0.0;
    double TI = 0.0;

    double UT = 0.0;

    double dx = 0.1;
    double dy = 0.1;

    double tol = 0.00000001;
    int** geometry;
    #pragma endregion

    //Getting log file names
    #pragma region
    // command line arrays
    int argc = 3;
    const char* argv[] = { "./sim","-f","-t"};
    // File paths string 
    std::string szFileName;
    std::string geoFileName;
    std::string problemName;
    std::vector<std::string> logFilesName;
    // Problem Flags
    bool temp_flag;

    logFilesName  = ParseTestingProblemFlags(argc,argv,szFileName,geoFileName,problemName,temp_flag);
    #pragma endregion
    // Geometry Initialization
    #pragma region 
    geometry = read_pgm((geoFileName).c_str());
    #pragma endregion

    // Initialzie Dummy Grid
    #pragma region    
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry) ;
    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();
    #pragma endregion

    // Data Structure Declaration
    #pragma region 
    // Input Data Structures
    matrix<double> U;
    matrix<double> V;
    // Calculation Data Structures
    matrix<double> F_cal;
    matrix<double> G_cal;
    // Reference Data Structures
    matrix<double> F_ref;
    matrix<double> G_ref;
    #pragma endregion
  
  // Data Structure Allocation
  #pragma region 
    // Input Data Structures
    U.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V.resize(imaxb,std::vector<double>(jmaxb,0.0));
    
    // Allocating memory for matrices
    F_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    F_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    #pragma endregion
  
  // Reading log files in data structures
  #pragma region 

  // Input Data Structures  
  read_matrix(logFilesName[static_cast<int>(log_file::log_U_Input)],U,imaxb,jmaxb);
  read_matrix(logFilesName[static_cast<int>(log_file::log_V_Input)],V,imaxb,jmaxb);
  
  // Reference Data Structures
  read_matrix(logFilesName[static_cast<int>(log_file::log_F)],F_ref,imaxb,jmaxb);
  read_matrix(logFilesName[static_cast<int>(log_file::log_G)],G_ref,imaxb,jmaxb);
  #pragma endregion
     
  
   

    /* Real Scenario*/
    #pragma region 
  
    domain.set_velocity(U, velocity_type::U);
    domain.set_velocity(V, velocity_type::V);
    

    calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,PR,&temp_flag,domain); 
    
    calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, domain, F_cal, G_cal);
    #pragma endregion

    // Writing calculated matrices to log files
  #pragma region 
    write_matrix(cmp_log_F_t, F_cal, imaxb, jmaxb);
    write_matrix(cmp_log_G_t, G_cal, imaxb, jmaxb);
    #pragma endregion
  
  // Comparing Reference matrices to Calculated matrices
  #pragma region 

    for(int i = 0; i < imaxb; i++)
    {
      for(int j = 0; j < jmaxb; j++)
      {
        INFO("i = " << i << " j = " << j);
        REQUIRE(fabs(F_cal[i][j] - F_ref[i][j]) < tol);
        REQUIRE(fabs(G_cal[i][j] - G_ref[i][j]) < tol);
      }
    }

    #pragma endregion
  
}

TEST_CASE("Test calculate dt & F and G and RS matrix Flow Over A Step", "[calculate_dt_fg_RS_FS]")
{
    // Dummy Grid Paramters
    #pragma region 
    
    double xlength = 10;
    double ylength = 2;

    int imax = 98;
    int jmax = 18;

    double dt = 0.05;
    double dt_ref = 0.0;
    double t_end = 500.0;
    double tau = 0.5;

    double dt_Value = 10.0;

    int itermax = 500;
    double eps = 0.001;
    double omg = 1.7;
    double alpha = 0.9;

    double Re = 100;
    double PR = 1;
    double GX = 0.0;
    double GY = 0.0;

    double UI = 0.0;
    double VI = 0.0;
    double PI = 0.0;
    double TI = 0.0;

    double UT = 0.0;

    double dx = 0.1;
    double dy = 0.1;

    double tol = 0.00000001;
    int** geometry;
    #pragma endregion
    //Getting log file names
    #pragma region
    // command line arrays
    int argc = 3;
    const char* argv[] = { "./sim","-f","-t"};
    // File paths string 
    std::string szFileName;
    std::string geoFileName;
    std::string problemName;
    std::vector<std::string> logFilesName;
    // Problem Flags
    bool temp_flag;

    logFilesName  = ParseTestingProblemFlags(argc,argv,szFileName,geoFileName,problemName,temp_flag);
    #pragma endregion
    // Geometry Initialization
    #pragma region 
    geometry = read_pgm((geoFileName).c_str());
    #pragma endregion

    // Initialzie Dummy Grid
    #pragma region    
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry) ;
    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();
    #pragma endregion

    // Data Structure Declaration
  #pragma region 
    // Input Data Structures
    matrix<double> U;
    matrix<double> V;
    // Calculation Data Structures
    matrix<double> F_cal;
    matrix<double> G_cal;
    // Calculation Data Structures
    matrix<double> RS_cal;
    // Reference Data Structures
    matrix<double> RS_ref;
    #pragma endregion
  
  // Data Structure Allocation
  #pragma region 
    // Input Data Structures
    U.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    F_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    RS_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    RS_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    #pragma endregion
  
  // Reading log files in data structures
  #pragma region 

  // Input Data Structures  
  read_matrix(logFilesName[static_cast<int>(log_file::log_U_Input)],U,imaxb,jmaxb);
  read_matrix(logFilesName[static_cast<int>(log_file::log_V_Input)],V,imaxb,jmaxb);
  // Reference Data Structures
  read_matrix(logFilesName[static_cast<int>(log_file::log_RS)],RS_ref,imaxb,jmaxb);
  #pragma endregion
     
  
   

    /* Real Scenario*/
    #pragma region 
  
    domain.set_velocity(U, velocity_type::U);
    domain.set_velocity(V, velocity_type::V);
    
    calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,PR,&temp_flag,domain);  
    calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, domain, F_cal, G_cal);
    calculate_rs(dt, dx, dy, imax, jmax, F_cal, G_cal, RS_cal,domain);
    #pragma endregion

    #pragma endregion
    // Writing calculated matrices to log files
    #pragma region 
    write_matrix(cmp_log_RS_t, RS_cal, imaxb, jmaxb);
    
    #pragma endregion
    // Comparing Reference matrices to Calculated matrices
    #pragma region 
    for(int i = 0; i < imaxb; i++)
    {
      for(int j = 0; j < jmaxb; j++)
      {
        if(RS_ref[i][j] == 0 )
        {
          
            INFO("i = " << i << " j = " << j);
            REQUIRE(RS_cal[i][j] == RS_ref[i][j]);
        }
      }
    }

    for(int i = 0; i < imaxb; i++)
    {
      for(int j = 0; j < jmaxb; j++)
      {
        INFO("i = " << i << " j = " << j);
        REQUIRE(fabs(RS_cal[i][j] - RS_ref[i][j]) < tol);
        
      }
    }

    #pragma endregion
}

TEST_CASE("Test calculate dt & F and G and RS matrix sor Flow Over A Step", "[calculate_dt_fg_RS_sor_FS]")
{
    // Dummy Grid Paramters
    #pragma region 
    
    double xlength = 10;
    double ylength = 2;

    int imax = 98;
    int jmax = 18;

    double dt = 0.05;
    double dt_ref = 0.0;
    double t_end = 500.0;
    double tau = 0.5;

    double dt_Value = 10.0;

    int itermax = 500;
    double eps = 0.001;
    double omg = 1.7;
    double alpha = 0.9;

    double Re = 100;
    double PR = 1;
    double GX = 0.0;
    double GY = 0.0;

    double UI = 0.0;
    double VI = 0.0;
    double PI = 0.0;
    double TI = 0.0;

    double UT = 0.0;

    double dx = 0.1;
    double dy = 0.1;

    double tol = 0.00000001;
    int** geometry;
    #pragma endregion
    //Getting log file names
    #pragma region
    // command line arrays
    int argc = 3;
    const char* argv[] = { "./sim","-f","-t"};
    // File paths string 
    std::string szFileName;
    std::string geoFileName;
    std::string problemName;
    std::vector<std::string> logFilesName;
    // Problem Flags
    bool temp_flag;

    logFilesName  = ParseTestingProblemFlags(argc,argv,szFileName,geoFileName,problemName,temp_flag);
    #pragma endregion
    // Geometry Initialization
    #pragma region 
    geometry = read_pgm((geoFileName).c_str());
    #pragma endregion

    // Initialzie Dummy Grid
    #pragma region    
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry) ;
    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();
    #pragma endregion

    // Data Structure Declaration
  #pragma region 
    // Input Data Structures
    matrix<double> U;
    matrix<double> V;
    // Calculation Data Structures
    matrix<double> F_cal;
    matrix<double> G_cal;
    // Calculation Data Structures
    matrix<double> RS_cal;
    matrix<double> P_cal;
    // Reference Data Structures
    matrix<double> P_ref;
    #pragma endregion
  
  // Data Structure Allocation
  #pragma region 
    // Input Data Structures
    U.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    F_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    RS_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    P_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    P_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    #pragma endregion
  
  // Reading log files in data structures
  #pragma region 

  // Input Data Structures  
  read_matrix(logFilesName[static_cast<int>(log_file::log_U_Input)],U,imaxb,jmaxb);
  read_matrix(logFilesName[static_cast<int>(log_file::log_V_Input)],V,imaxb,jmaxb);
  read_matrix(logFilesName[static_cast<int>(log_file::log_P_Input)],P_cal,imaxb,jmaxb);
  // Reference Data Structures
  read_matrix(logFilesName[static_cast<int>(log_file::log_P_Output)],P_ref,imaxb,jmaxb);
  #pragma endregion
     
  
   

    /* Real Scenario*/
    #pragma region 
  
    domain.set_velocity(U, velocity_type::U);
    domain.set_velocity(V, velocity_type::V);
    domain.set_pressure(P_cal);

    calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,PR,&temp_flag,domain);  
    calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, domain, F_cal, G_cal);
    calculate_rs(dt, dx, dy, imax, jmax, F_cal, G_cal, RS_cal,domain);
    
    int it = 0;
    double res = 1.0;
    // Solve system using SOR
    while(it < itermax && res > eps){
      sor(omg,dx,dy,imax,jmax,domain,RS_cal,&res);
      it++;
    }
    

    #pragma endregion

    #pragma endregion
   
    // Writing calculated matrices to log files
    #pragma region 
    domain.pressure(P_cal);
    write_matrix(cmp_log_P_t, P_cal, imaxb, jmaxb);
    
    #pragma endregion
   
    // Comparing Reference matrices to Calculated matrices
    #pragma region 
   
    for(int i = 0; i < imaxb; i++)
    {
      for(int j = 0; j < jmaxb; j++)
      {
        INFO("i = " << i << " j = " << j);
        REQUIRE(fabs(P_cal[i][j] - P_ref[i][j]) < tol);
        
      }
    }

    #pragma endregion
    
}

TEST_CASE("Test calculate dt & F and G and RS matrix sor uv Flow Over A Step", "[calculate_dt_fg_RS_sor_uv_FS]")
{
   // Dummy Grid Paramters
    #pragma region 
    
    double xlength = 10;
    double ylength = 2;

    int imax = 98;
    int jmax = 18;

    double dt = 0.05;
    double dt_ref = 0.0;
    double t_end = 500.0;
    double tau = 0.5;

    double dt_Value = 10.0;

    int itermax = 500;
    double eps = 0.001;
    double omg = 1.7;
    double alpha = 0.9;

    double Re = 100;
    double PR = 1;
    double GX = 0.0;
    double GY = 0.0;

    double UI = 0.0;
    double VI = 0.0;
    double PI = 0.0;
    double TI = 0.0;

    double UT = 0.0;

    double dx = 0.1;
    double dy = 0.1;

    double tol = 0.00000001;
    int** geometry;
    #pragma endregion
    //Getting log file names
    #pragma region
    // command line arrays
    int argc = 3;
    const char* argv[] = { "./sim","-f","-t"};
    // File paths string 
    std::string szFileName;
    std::string geoFileName;
    std::string problemName;
    std::vector<std::string> logFilesName;
    // Problem Flags
    bool temp_flag;

    logFilesName  = ParseTestingProblemFlags(argc,argv,szFileName,geoFileName,problemName,temp_flag);
    #pragma endregion
    // Geometry Initialization
    #pragma region 
    geometry = read_pgm((geoFileName).c_str());
    #pragma endregion

    // Initialzie Dummy Grid
    #pragma region    
    Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry) ;
    int imaxb = domain.imaxb();
    int jmaxb = domain.jmaxb();
    #pragma endregion

    // Data Structure Declaration
  #pragma region 
    // Input Data Structures
    matrix<double> U;
    matrix<double> V;
    // Calculation Data Structures
    matrix<double> F_cal;
    matrix<double> G_cal;
    // Calculation Data Structures
    matrix<double> RS_cal;
    matrix<double> P_cal;
    // Reference Data Structures
    matrix<double> U_ref;
    matrix<double> V_ref;
    #pragma endregion
  
  // Data Structure Allocation
  #pragma region 
    // Input Data Structures
    U.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    F_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    G_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    RS_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    P_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
    // Allocating memory for matrices
    U_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    V_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
    #pragma endregion
  
  // Reading log files in data structures
  #pragma region 

  // Input Data Structures  
  read_matrix(logFilesName[static_cast<int>(log_file::log_U_Input)],U,imaxb,jmaxb);
  read_matrix(logFilesName[static_cast<int>(log_file::log_V_Input)],V,imaxb,jmaxb);
  read_matrix(logFilesName[static_cast<int>(log_file::log_P_Input)],P_cal,imaxb,jmaxb);
  // Reference Data Structures
    read_matrix(logFilesName[static_cast<int>(log_file::log_U_Output)],U_ref,imaxb,jmaxb);
    read_matrix(logFilesName[static_cast<int>(log_file::log_V_Output)],V_ref,imaxb,jmaxb);
  #pragma endregion
     
  
   

    /* Real Scenario*/
    #pragma region 
  
    domain.set_velocity(U, velocity_type::U);
    domain.set_velocity(V, velocity_type::V);
    domain.set_pressure(P_cal);

    calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,PR,&temp_flag,domain);  
    calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, domain, F_cal, G_cal);
    calculate_rs(dt, dx, dy, imax, jmax, F_cal, G_cal, RS_cal,domain);
    
    int it = 0;
    double res = 1.0;
    // Solve system using SOR
    while(it < itermax && res > eps){
      sor(omg,dx,dy,imax,jmax,domain,RS_cal,&res);
      it++;
    }
    
    calculate_uv(dt, dx, dy, imax, jmax, domain, F_cal, G_cal);

    #pragma endregion

    #pragma endregion
   
    // Writing calculated matrices to log files
    #pragma region 
    domain.velocity(U,velocity_type::U);
    domain.velocity(V,velocity_type::V);
    write_matrix(cmp_log_U_t, U , imaxb, jmaxb);
    write_matrix(cmp_log_V_t, V , imaxb, jmaxb);
    #pragma endregion
    
    // Comparing Reference matrices to Calculated matrices
    #pragma region 
    

    for(int i = 1; i < imaxb-1; i++)
    {
      for(int j = 1; j < jmaxb-1; j++)
      {
        INFO("i = " << i << " j = " << j);
        REQUIRE(fabs(U[i][j] - U_ref[i][j]) < tol);
        REQUIRE(fabs(V[i][j] - V_ref[i][j]) < tol);
      }
    }

    #pragma endregion
    
}


// Natural Convection  Integration Tests


// TEST_CASE("Test calculate dt & T and F and G matrix Natural Convection ", "[calculate_dt_fg_NC]")
// {
//     // Dummy Grid Paramters
//     #pragma region 
    
//     double xlength = 1;
//     double ylength = 1;

//     int imax = 48;
//     int jmax = 48;

//     double dt = 0.05;
//     double t_end = 1000.0;
//     double tau = 0.5;

//     double dt_Value = 10.0;

//     int itermax = 100;
//     double eps = 0.00001;
//     double omg = 1.7;
//     double alpha = 0.5;
//     double beta = 0.00021;

//     double Re = 1000;
//     double PR = 7;
//     double GX = 0.0;
//     double GY = -1.1;

//     double UI = 0.0;
//     double VI = 0.0;
//     double PI = 0.0;
//     double TI = 0.0;

//     double UT = 0.0;

//     double dx = 0.02;
//     double dy = 0.02;

//     double tol = 0.0000001;
//     int** geometry;
//     #pragma endregion

//     //Getting log file names
//     #pragma region
//     // command line arrays
//     int argc = 3;
//     const char* argv[] = { "./sim","-n","-t"};
//     // File paths string 
//     std::string szFileName;
//     std::string geoFileName;
//     std::string problemName;
//     std::vector<std::string> logFilesName;
//     // Problem Flags
//     bool temp_flag;

//     logFilesName  = ParseTestingProblemFlags(argc,argv,szFileName,geoFileName,problemName,temp_flag);
//     #pragma endregion
//     // Geometry Initialization
//     #pragma region 
//     geometry = read_pgm((geoFileName).c_str());
//     #pragma endregion

//     // Initialzie Dummy Grid
//     #pragma region    
//     Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry) ;
//     int imaxb = domain.imaxb();
//     int jmaxb = domain.jmaxb();
//     #pragma endregion

//     // Data Structure Declaration
//   #pragma region 
//     // Input Data Structures
//     matrix<double> U;
//     matrix<double> V;
//     matrix<double> T;
//     // Calculation Data Structures
//     matrix<double> F_cal;
//     matrix<double> G_cal;
//     // Reference Data Structures
//     matrix<double> F_ref;
//     matrix<double> G_ref;
//     #pragma endregion
  
//   // Data Structure Allocation
//   #pragma region 
//     // Input Data Structures
//     U.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     V.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     T.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     // Allocating memory for matrices
//     F_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     G_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     // Allocating memory for matrices
//     F_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     G_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     #pragma endregion
  
//   // Reading log files in data structures
//   #pragma region 

//   // Input Data Structures  
//   read_matrix(logFilesName[static_cast<int>(log_file::log_U_Input)],U,imaxb,jmaxb);
//   read_matrix(logFilesName[static_cast<int>(log_file::log_V_Input)],V,imaxb,jmaxb);
//   read_matrix(logFilesName[static_cast<int>(log_file::log_T_Input)],T,imaxb,jmaxb);
//   // Reference Data Structures
//   read_matrix(logFilesName[static_cast<int>(log_file::log_F)],F_ref,imaxb,jmaxb);
//   read_matrix(logFilesName[static_cast<int>(log_file::log_G)],G_ref,imaxb,jmaxb);
//   #pragma endregion
     
  
   

//     /* Real Scenario*/
//     #pragma region 
  
//     domain.set_velocity(U, velocity_type::U);
//     domain.set_velocity(V, velocity_type::V);
//     domain.set_temperature(T);    
//     calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,PR,&temp_flag,domain);  
//     calculate_temp(dt,dx,dy,imax,jmax,Re,PR,alpha,0,domain); 
//     calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, domain, F_cal, G_cal);
//     #pragma endregion

//     // Writing calculated matrices to log files
//   #pragma region 
//     write_matrix(cmp_log_F_t, F_cal, imaxb, jmaxb);
//     write_matrix(cmp_log_G_t, G_cal, imaxb, jmaxb);
//     #pragma endregion
  
//   // Comparing Reference matrices to Calculated matrices
//   #pragma region 

//     for(int i = 0; i < imaxb; i++)
//     {
//       for(int j = 0; j < jmaxb; j++)
//       {
//         INFO("i = " << i << " j = " << j);
//         REQUIRE(fabs(F_cal[i][j] - F_ref[i][j]) < tol);
//         REQUIRE(fabs(G_cal[i][j] - G_ref[i][j]) < tol);
//       }
//     }

//     #pragma endregion
  
// }

// TEST_CASE("Test calculate dt & T and F and G and RS matrix Natural Convection ", "[calculate_dt_fg_RS_NC]")
// {
//     // Dummy Grid Paramters
//     #pragma region 
    
//     double xlength = 1;
//     double ylength = 1;

//     int imax = 48;
//     int jmax = 48;

//     double dt = 0.05;
//     double t_end = 1000.0;
//     double tau = 0.5;

//     double dt_Value = 10.0;

//     int itermax = 100;
//     double eps = 0.00001;
//     double omg = 1.7;
//     double alpha = 0.5;
//     double beta = 0.00021;

//     double Re = 1000;
//     double PR = 7;
//     double GX = 0.0;
//     double GY = -1.1;

//     double UI = 0.0;
//     double VI = 0.0;
//     double PI = 0.0;
//     double TI = 0.0;

//     double UT = 0.0;

//     double dx = 0.02;
//     double dy = 0.02;

//     double tol = 0.0000001;
//     int** geometry;
//     #pragma endregion

//     //Getting log file names
//     #pragma region
//     // command line arrays
//     int argc = 3;
//     const char* argv[] = { "./sim","-n","-t"};
//     // File paths string 
//     std::string szFileName;
//     std::string geoFileName;
//     std::string problemName;
//     std::vector<std::string> logFilesName;
//     // Problem Flags
//     bool temp_flag;

//     logFilesName  = ParseTestingProblemFlags(argc,argv,szFileName,geoFileName,problemName,temp_flag);
//     #pragma endregion
//     // Geometry Initialization
//     #pragma region 
//     geometry = read_pgm((geoFileName).c_str());
//     #pragma endregion

//     // Initialzie Dummy Grid
//     #pragma region    
//     Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry) ;
//     int imaxb = domain.imaxb();
//     int jmaxb = domain.jmaxb();
//     #pragma endregion

//     // Data Structure Declaration
//   #pragma region 
//     // Input Data Structures
//     matrix<double> U;
//     matrix<double> V;
//     // Calculation Data Structures
//     matrix<double> F_cal;
//     matrix<double> G_cal;
//     // Calculation Data Structures
//     matrix<double> RS_cal;
//     // Reference Data Structures
//     matrix<double> RS_ref;
//     #pragma endregion
  
//   // Data Structure Allocation
//   #pragma region 
//     // Input Data Structures
//     U.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     V.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     // Allocating memory for matrices
//     F_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     G_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     RS_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     // Allocating memory for matrices
//     RS_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     #pragma endregion
  
//   // Reading log files in data structures
//   #pragma region 

//   // Input Data Structures  
//   read_matrix(logFilesName[static_cast<int>(log_file::log_U_Input)],U,imaxb,jmaxb);
//   read_matrix(logFilesName[static_cast<int>(log_file::log_V_Input)],V,imaxb,jmaxb);
//   // Reference Data Structures
//   read_matrix(logFilesName[static_cast<int>(log_file::log_RS)],RS_ref,imaxb,jmaxb);
//   #pragma endregion
     
  
   

//     /* Real Scenario*/
//     #pragma region 
  
//     domain.set_velocity(U, velocity_type::U);
//     domain.set_velocity(V, velocity_type::V);
    
//     calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,PR,&temp_flag,domain);  
//     calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, domain, F_cal, G_cal);
//     calculate_rs(dt, dx, dy, imax, jmax, F_cal, G_cal, RS_cal,domain);
//     #pragma endregion

//     #pragma endregion
//     // Writing calculated matrices to log files
//     #pragma region 
//     write_matrix(cmp_log_RS_t, RS_cal, imaxb, jmaxb);
    
//     #pragma endregion
//     // Comparing Reference matrices to Calculated matrices
//     #pragma region 
//     for(int i = 0; i < imaxb; i++)
//     {
//       for(int j = 0; j < jmaxb; j++)
//       {
//         if(RS_ref[i][j] == 0 )
//         {
          
//             INFO("i = " << i << " j = " << j);
//             REQUIRE(RS_cal[i][j] == RS_ref[i][j]);
//         }
//       }
//     }

//     for(int i = 0; i < imaxb; i++)
//     {
//       for(int j = 0; j < jmaxb; j++)
//       {
//         INFO("i = " << i << " j = " << j);
//         REQUIRE(fabs(RS_cal[i][j] - RS_ref[i][j]) < tol);
        
//       }
//     }

//     #pragma endregion
// }

// TEST_CASE("Test calculate dt & T and F and G and RS matrix sor Natural Convection ", "[calculate_dt_fg_RS_sor_NC]")
// {
//     // Dummy Grid Paramters
//     #pragma region 
    
//     double xlength = 1;
//     double ylength = 1;

//     int imax = 48;
//     int jmax = 48;

//     double dt = 0.05;
//     double t_end = 1000.0;
//     double tau = 0.5;

//     double dt_Value = 10.0;

//     int itermax = 100;
//     double eps = 0.00001;
//     double omg = 1.7;
//     double alpha = 0.5;
//     double beta = 0.00021;

//     double Re = 1000;
//     double PR = 7;
//     double GX = 0.0;
//     double GY = -1.1;

//     double UI = 0.0;
//     double VI = 0.0;
//     double PI = 0.0;
//     double TI = 0.0;

//     double UT = 0.0;

//     double dx = 0.02;
//     double dy = 0.02;

//     double tol = 0.0000001;
//     int** geometry;
//     #pragma endregion

//     //Getting log file names
//     #pragma region
//     // command line arrays
//     int argc = 3;
//     const char* argv[] = { "./sim","-n","-t"};
//     // File paths string 
//     std::string szFileName;
//     std::string geoFileName;
//     std::string problemName;
//     std::vector<std::string> logFilesName;
//     // Problem Flags
//     bool temp_flag;

//     logFilesName  = ParseTestingProblemFlags(argc,argv,szFileName,geoFileName,problemName,temp_flag);
//     #pragma endregion
//     // Geometry Initialization
//     #pragma region 
//     geometry = read_pgm((geoFileName).c_str());
//     #pragma endregion

//     // Initialzie Dummy Grid
//     #pragma region    
//     Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry) ;
//     int imaxb = domain.imaxb();
//     int jmaxb = domain.jmaxb();
//     #pragma endregion

//     // Data Structure Declaration
//   #pragma region 
//     // Input Data Structures
//     matrix<double> U;
//     matrix<double> V;
//     // Calculation Data Structures
//     matrix<double> F_cal;
//     matrix<double> G_cal;
//     // Calculation Data Structures
//     matrix<double> RS_cal;
//     matrix<double> P_cal;
//     // Reference Data Structures
//     matrix<double> P_ref;
//     #pragma endregion
  
//   // Data Structure Allocation
//   #pragma region 
//     // Input Data Structures
//     U.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     V.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     // Allocating memory for matrices
//     F_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     G_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     RS_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     P_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     // Allocating memory for matrices
//     P_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     #pragma endregion
  
//   // Reading log files in data structures
//   #pragma region 

//   // Input Data Structures  
//   read_matrix(logFilesName[static_cast<int>(log_file::log_U_Input)],U,imaxb,jmaxb);
//   read_matrix(logFilesName[static_cast<int>(log_file::log_V_Input)],V,imaxb,jmaxb);
//   read_matrix(logFilesName[static_cast<int>(log_file::log_P_Input)],P_cal,imaxb,jmaxb);
//   // Reference Data Structures
//   read_matrix(logFilesName[static_cast<int>(log_file::log_P_Output)],P_ref,imaxb,jmaxb);
//   #pragma endregion
     
  
   

//     /* Real Scenario*/
//     #pragma region 
  
//     domain.set_velocity(U, velocity_type::U);
//     domain.set_velocity(V, velocity_type::V);
//     domain.set_pressure(P_cal);

//     calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,PR,&temp_flag,domain);  
//     calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, domain, F_cal, G_cal);
//     calculate_rs(dt, dx, dy, imax, jmax, F_cal, G_cal, RS_cal,domain);
    
//     int it = 0;
//     double res = 1.0;
//     // Solve system using SOR
//     while(it < itermax && res > eps){
//       sor(omg,dx,dy,imax,jmax,domain,RS_cal,&res);
//       it++;
//     }
    

//     #pragma endregion

//     #pragma endregion
   
//     // Writing calculated matrices to log files
//     #pragma region 
//     domain.pressure(P_cal);
//     write_matrix(cmp_log_P_t, P_cal, imaxb, jmaxb);
    
//     #pragma endregion
   
//     // Comparing Reference matrices to Calculated matrices
//     #pragma region 
   
//     for(int i = 0; i < imaxb; i++)
//     {
//       for(int j = 0; j < jmaxb; j++)
//       {
//         INFO("i = " << i << " j = " << j);
//         REQUIRE(fabs(P_cal[i][j] - P_ref[i][j]) < tol);
        
//       }
//     }

//     #pragma endregion
    
// }

// TEST_CASE("Test calculate dt & T and F and G and RS matrix sor uv Natural Convection ", "[calculate_dt_fg_RS_sor_uv_NC]")
// {
//    // Dummy Grid Paramters
//     #pragma region 
    
//     double xlength = 1;
//     double ylength = 1;

//     int imax = 48;
//     int jmax = 48;

//     double dt = 0.05;
//     double t_end = 1000.0;
//     double tau = 0.5;

//     double dt_Value = 10.0;

//     int itermax = 100;
//     double eps = 0.00001;
//     double omg = 1.7;
//     double alpha = 0.5;
//     double beta = 0.00021;

//     double Re = 1000;
//     double PR = 7;
//     double GX = 0.0;
//     double GY = -1.1;

//     double UI = 0.0;
//     double VI = 0.0;
//     double PI = 0.0;
//     double TI = 0.0;

//     double UT = 0.0;

//     double dx = 0.02;
//     double dy = 0.02;

//     double tol = 0.0000001;
//     int** geometry;
//     #pragma endregion

//     //Getting log file names
//     #pragma region
//     // command line arrays
//     int argc = 3;
//     const char* argv[] = { "./sim","-n","-t"};
//     // File paths string 
//     std::string szFileName;
//     std::string geoFileName;
//     std::string problemName;
//     std::vector<std::string> logFilesName;
//     // Problem Flags
//     bool temp_flag;

//     logFilesName  = ParseTestingProblemFlags(argc,argv,szFileName,geoFileName,problemName,temp_flag);
//     #pragma endregion
//     // Geometry Initialization
//     #pragma region 
//     geometry = read_pgm((geoFileName).c_str());
//     #pragma endregion

//     // Initialzie Dummy Grid
//     #pragma region    
//     Grid domain = Grid(imax,jmax,boundary_size,PI,UI,VI,TI,geometry) ;
//     int imaxb = domain.imaxb();
//     int jmaxb = domain.jmaxb();
//     #pragma endregion

//     // Data Structure Declaration
//   #pragma region 
//     // Input Data Structures
//     matrix<double> U;
//     matrix<double> V;
//     // Calculation Data Structures
//     matrix<double> F_cal;
//     matrix<double> G_cal;
//     // Calculation Data Structures
//     matrix<double> RS_cal;
//     matrix<double> P_cal;
//     // Reference Data Structures
//     matrix<double> U_ref;
//     matrix<double> V_ref;
//     #pragma endregion
  
//   // Data Structure Allocation
//   #pragma region 
//     // Input Data Structures
//     U.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     V.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     // Allocating memory for matrices
//     F_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     G_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     RS_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     P_cal.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     // Allocating memory for matrices
//     U_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     V_ref.resize(imaxb,std::vector<double>(jmaxb,0.0));
//     #pragma endregion
  
//   // Reading log files in data structures
//   #pragma region 

//   // Input Data Structures  
//   read_matrix(logFilesName[static_cast<int>(log_file::log_U_Input)],U,imaxb,jmaxb);
//   read_matrix(logFilesName[static_cast<int>(log_file::log_V_Input)],V,imaxb,jmaxb);
//   read_matrix(logFilesName[static_cast<int>(log_file::log_P_Input)],P_cal,imaxb,jmaxb);
//   // Reference Data Structures
//     read_matrix(logFilesName[static_cast<int>(log_file::log_U_Output)],U_ref,imaxb,jmaxb);
//     read_matrix(logFilesName[static_cast<int>(log_file::log_V_Output)],V_ref,imaxb,jmaxb);
//   #pragma endregion
     
  
   

//     /* Real Scenario*/
//     #pragma region 
  
//     domain.set_velocity(U, velocity_type::U);
//     domain.set_velocity(V, velocity_type::V);
//     domain.set_pressure(P_cal);

//     calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,PR,&temp_flag,domain);  
//     calculate_fg(Re, GX, GY, alpha, dt, dx, dy, imax, jmax, domain, F_cal, G_cal);
//     calculate_rs(dt, dx, dy, imax, jmax, F_cal, G_cal, RS_cal,domain);
    
//     int it = 0;
//     double res = 1.0;
//     // Solve system using SOR
//     while(it < itermax && res > eps){
//       sor(omg,dx,dy,imax,jmax,domain,RS_cal,&res);
//       it++;
//     }
    
//     calculate_uv(dt, dx, dy, imax, jmax, domain, F_cal, G_cal);

//     #pragma endregion

//     #pragma endregion
   
//     // Writing calculated matrices to log files
//     #pragma region 
//     domain.velocity(U,velocity_type::U);
//     domain.velocity(V,velocity_type::V);
//     write_matrix(cmp_log_U_t, U , imaxb, jmaxb);
//     write_matrix(cmp_log_V_t, V , imaxb, jmaxb);
//     #pragma endregion
    
//     // Comparing Reference matrices to Calculated matrices
//     #pragma region 
    

//     for(int i = 1; i < imaxb-1; i++)
//     {
//       for(int j = 1; j < jmaxb-1; j++)
//       {
//         INFO("i = " << i << " j = " << j);
//         REQUIRE(fabs(U[i][j] - U_ref[i][j]) < tol);
//         REQUIRE(fabs(V[i][j] - V_ref[i][j]) < tol);
//       }
//     }

//     #pragma endregion
    
// }





