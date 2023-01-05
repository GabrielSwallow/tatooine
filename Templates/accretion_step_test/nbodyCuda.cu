extern "C" {

#include "pluto.h"
#include "plutoCuda.cuh"

__device__ real cudaNbodySmoothingSquared(int l, uint3 dataIdx, dim3 dataDim)
{
    real R = cudaGrid.x1[dataIdx.x];
    return G_SMOOTHING*G_SMOOTHING
          *cudaInputParam[ASPECT_RATIO]
          *cudaInputParam[ASPECT_RATIO]
          *R*R;
}

} /* extern "C" */
