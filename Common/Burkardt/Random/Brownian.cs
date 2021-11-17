using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace Burkardt.RandomNS;

public static partial class BRandom
{
    public static double[] brownian ( int dim_num, int n, ref typeMethods.r8NormalData data, ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BROWNIAN creates Brownian motion points.
        //
        //  Discussion:
        //
        //    A starting point is generated at the origin.  The next point
        //    is generated at a uniformly random angle and a (0,1) normally
        //    distributed distance from the previous point.
        //
        //    It is up to the user to rescale the data, if desired.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 August 2004
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
        //    Input, int &SEED, a seed for the random number generator.
        //
        //    Output, double BROWNIAN[DIM_NUM*N], the Brownian motion points.
        //
    {
        double[] direction;
        int i;
        int j;
        double r;
        double[] x;

        direction = new double[dim_num];
        x = new double[dim_num*n];
        //
        //  Initial point.
        //
        j = 0;
        for ( i = 0; i < dim_num; i++ )
        {
            x[i+j*dim_num] = 0.0;
        }
        //
        //  Generate angles and steps.
        //
        for ( j = 1; j < n; j++ )
        {
            r = typeMethods.r8_normal_01 ( ref data, ref seed );
            r = Math.Abs ( r );

            direction_uniform_nd ( dim_num, ref seed, ref direction );

            for ( i = 0; i < dim_num; i++ )
            {
                x[i+j*dim_num] = x[i+(j-1)*dim_num] + r * direction[i];
            }

        }
            
        return x;
    }
}