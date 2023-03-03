#ifndef EOS_CUDA_CUH
#define EOS_CUDA_CUH

__forceinline__ __device__ real calcSoundSpeed(real v[], int pos, int dir, 
                                               uint3 dataIdx)
{
    #if BINARY_H == YES
    real x1 = cudaGrid.x1[dataIdx.x];
    real x2 = cudaGrid.x2[dataIdx.y];

    if (dir == IDIR)
    {
        x1 = (pos == FACE_CENTER ? cudaGrid.x1_if[dataIdx.x+1] 
                                 : cudaGrid.x1[dataIdx.x]);
    }
    else if (dir == JDIR)
    {
        x2 = (pos == FACE_CENTER ? cudaGrid.x2_if[dataIdx.y+1]
                                 : cudaGrid.x2[dataIdx.y]);
    }

    real x = x1*cos(x2);
    real y = x1*sin(x2);

    real K = 0.0;
    for (int l = 0; l < CENTRAL_OBJECT; l++)
    {
        real dx = x - cudaNb.x[l]; 
        real dy = y - cudaNb.y[l]; 

        real r3 = dx*dx + dy*dy;
        r3 *= sqrt(r3);

        K += CONST_G_CODE_UNITS*cudaNb.m[l] / r3;
    }

    return cudaIsoSoundSpeed * x1 * sqrt(K);
    #else
    real x1 = cudaGrid.x1[dataIdx.x];

    if (dir == IDIR)
    {
        x1 = (pos == FACE_CENTER ? cudaGrid.x1_if[dataIdx.x+1] 
                                 : cudaGrid.x1[dataIdx.x]);
    }

    return cudaIsoSoundSpeed * sqrt(CONST_G_CODE_UNITS*cudaInputParam[M_CO]/x1);
    #endif
}

#endif /* EOS_CUDA_CUH */
