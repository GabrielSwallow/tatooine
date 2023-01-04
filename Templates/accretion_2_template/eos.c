#include "pluto.h"

/* **************************************************************** */
void SoundSpeed2 (double **u, double *cs2, double *h, int beg, int end,
                  int pos, Grid *grid)
/*
 *
 *    Define the square of the sound speed for different EOS
 *
 ****************************************************************** */
{
    int i, j, k;
    double *x1, *x2, *x3;

    x1 = grid[IDIR].x;
    x2 = grid[JDIR].x;
    x3 = grid[KDIR].x;

    i = g_i;
    j = g_j;
    k = g_k;

    #if BINARY_H == YES

    if (g_dir == IDIR)
    {
        x1 = (pos == FACE_CENTER ? grid[IDIR].xr : grid[IDIR].x);

        for (i = beg; i <= end; i++)
        {
            real K = 0.0;
            for (int l = 0; l < CENTRAL_OBJECT; l++)
            {
                real dx = x1[i]*cos(x2[j]) - g_nb.x[l];
                real dy = x1[i]*sin(x2[j]) - g_nb.y[l];

                real r3 = dx*dx + dy*dy;
                r3 *= sqrt(r3);

                K += CONST_G_CODE_UNITS*g_nb.m[l] / r3;
            }
            cs2[i] = g_isoSoundSpeed*g_isoSoundSpeed*x1[i]*x1[i]*K;
        }
    }
    else if (g_dir == JDIR)
    {
        x2 = (pos == FACE_CENTER ? grid[JDIR].xr : grid[JDIR].x);

        for (j = beg; j <= end; j++)
        {
            real K = 0.0;
            for (int l = 0; l < CENTRAL_OBJECT; l++)
            {
                real dx = x1[i]*cos(x2[j]) - g_nb.x[l];
                real dy = x1[i]*sin(x2[j]) - g_nb.y[l];

                real r3 = dx*dx + dy*dy;
                r3 *= sqrt(r3);

                K += CONST_G_CODE_UNITS*g_nb.m[l] / r3;
            }
            cs2[j] = g_isoSoundSpeed*g_isoSoundSpeed*x1[i]*x1[i]*K;
        }
    }
    #else
    if (g_dir == IDIR) 
    {
        x1 = (pos == FACE_CENTER ? grid[IDIR].xr : grid[IDIR].x);

        for (i = beg; i <= end; i++)
            cs2[i] = g_isoSoundSpeed*g_isoSoundSpeed
                    *CONST_G_CODE_UNITS*g_inputParam[M_CO]/x1[i];
    } 
    else if (g_dir == JDIR)
    {
        x2 = (pos == FACE_CENTER ? grid[JDIR].xr : grid[JDIR].x);

        for (j = beg; j <= end; j++) 
            cs2[j] = g_isoSoundSpeed*g_isoSoundSpeed
                    *CONST_G_CODE_UNITS*g_inputParam[M_CO]/x1[i];
    }
    #endif
}

/* *************************************************************** */
void Enthalpy (real **uprim, real *h, int beg, int end)
/*
 *
 *
 *
 ***************************************************************** */
{
  int i;
  double g_gammar;

  #if EOS == IDEAL
   g_gammar = g_gamma/(g_gamma - 1.0);
   for (i = beg; i <= end; i++){
     h[i] = g_gammar*uprim[i][PRS]/uprim[i][RHO];
   }
  #elif EOS == ISOTHERMAL 
   print (" Enthalpy not defined for isothermal EoS\n");
   QUIT_PLUTO(1);
  #endif
}

/* *************************************************************** */
void ENTROPY (real **v, real *s, int is, int ie)
/*
 *
 *
 *
 ***************************************************************** */
{
  int i;
  double rho;

  #if EOS == IDEAL
   for (i = is; i <= ie; i++){
     rho  = v[i][RHO];
     s[i] = v[i][PRS]/pow(rho,g_gamma);
   }
  #elif EOS == ISOTHERMAL || EOS == BAROTROPIC
   print (" Entropy not defined in isothermal or barotropic MHD\n");
   QUIT_PLUTO(1);
  #endif
}
