﻿using System;

namespace Burkardt.Uniform;

public static class Walk
{
    public static double[] uniform_walk ( int dim_num, int n, ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_WALK generates points on a uniform random walk.
        //
        //  Discussion:
        //
        //    The first point is at the origin.  Uniform random numbers are
        //    generated to determine the direction of the next step, which
        //    is always of length 1, and in coordinate direction.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 August 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the dimension of the space.
        //
        //    Input, int N, the number of points.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double UNIFORM_WALK[DIM_NUM*N], the points.
        //
    {
        int i;

        double[] x = new double[dim_num*n];

        int j = 0;
        for ( i = 0; i < dim_num; i++ )
        {
            x[i+j*dim_num] = 0.0;
        }

        for ( j = 1; j < n; j++ )
        {
            double dir = UniformRNG.r8_uniform_01 ( ref seed );
            dir = 2 * dim_num * ( dir - 0.5 );

            for ( i = 0; i < dim_num; i++ )
            {
                x[i+j*dim_num] = x[i+(j-1)*dim_num];
            }
            double arg = Math.Abs ( dir ) + 0.5;
            i = (int) arg;
            i = Math.Min ( i, dim_num );
            i = Math.Max ( i, 1 );
            i -= 1;

            switch (dir)
            {
                case < 0.0:
                    x[i+j*dim_num] -= 1.0;
                    break;
                default:
                    x[i+j*dim_num] += 1.0;
                    break;
            }

        }

        return x;
    }
}