#include "pluto.h"

double nbodySmoothingSquared(int l, double *v, double x1, double x2, double x3)
{
    return G_SMOOTHING*G_SMOOTHING
          *g_inputParam[ASPECT_RATIO]*g_inputParam[ASPECT_RATIO]
          *x1*x1;
}
