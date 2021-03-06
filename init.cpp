#include "init.hpp"
#include "helper.hpp"

int read_parameters(std::string szFileName, /* name of the file */
                    double *Re,             /* reynolds number   */
                    double *PR,             /* prandtl number   */
                    double *UI,             /* velocity x-direction */
                    double *VI,             /* velocity y-direction */
                    double *PI,             /* pressure */
                    double *TI,             /*temperature*/
                    double *GX,             /* gravitation x-direction */
                    double *GY,             /* gravitation y-direction */
                    double *t_end,          /* end time */
                    double *xlength,        /* length of the domain x-dir.*/
                    double *ylength,        /* length of the domain y-dir.*/
                    double *dt,             /* time step */
                    double *dx,             /* length of a cell x-dir. */
                    double *dy,             /* length of a cell y-dir. */
                    int *imax,              /* number of cells x-direction*/
                    int *jmax,              /* number of cells y-direction*/
                    double *alpha,          /* uppwind differencing factor*/
                    double *omg,            /* relaxation factor */
                    double *tau,            /* safety factor for time step*/
                    int *itermax,           /* max. number of iterations  */
                                            /* for pressure per time step */
                    double *eps,            /* accuracy bound for pressure*/
                    double *beta, double *dt_value, /* time for output */
                    double *UT) {
  // Reading Parameters
  get_file_double(szFileName, "UT", *UT);
  get_file_double(szFileName, "xlength", *xlength);
  get_file_double(szFileName, "ylength", *ylength);
  get_file_double(szFileName, "Re", *Re);
  get_file_double(szFileName, "PR", *PR);
  get_file_double(szFileName, "t_end", *t_end);
  get_file_double(szFileName, "dt", *dt);
  get_file_double(szFileName, "omg", *omg);
  get_file_double(szFileName, "eps", *eps);
  get_file_double(szFileName, "tau", *tau);
  get_file_double(szFileName, "alpha", *alpha);
  get_file_double(szFileName, "dt_value", *dt_value);
  get_file_double(szFileName, "UI", *UI);
  get_file_double(szFileName, "VI", *VI);
  get_file_double(szFileName, "GX", *GX);
  get_file_double(szFileName, "GY", *GY);
  get_file_double(szFileName, "PI", *PI);
  get_file_double(szFileName, "TI", *TI);

  get_file_double(szFileName, "beta", *beta);

  get_file_int(szFileName, "itermax", *itermax);
  get_file_int(szFileName, "imax", *imax);
  get_file_int(szFileName, "jmax", *jmax);

  *dx = *xlength / (double)(*imax);
  *dy = *ylength / (double)(*jmax);

  return 1;
}
