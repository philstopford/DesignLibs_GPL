using System;

namespace Burkardt.Uniform
{
    public static class Circle
    {
        public static double[] uniform_in_circle01_map ( int n, ref int seed )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    UNIFORM_IN_CIRCLE01_MAP maps uniform points into the unit circle.
            //
            //  Discussion:
            //
            //    The unit circle is centered at the origin and has radius 1.
            //
            //    This routine is valid for spatial dimension DIM_NUM = 2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 August 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, double UNIFORM_IN_CIRCLE01_MAP[DIM_NUM*N], the points.
            //
        {
            int DIM_NUM = 2;

            int j;
            double[] r;
            double[] t;
            double[] x;

            r = new double[n];
            t = new double[n];
            x = new double[DIM_NUM*n];

            for ( j = 0; j < n; j++ )
            {
                r[j] = UniformRNG.r8_uniform_01 ( ref seed );
                r[j] = Math.Sqrt ( r[j] );
            }

            for ( j = 0; j < n; j++ )
            {
                t[j] = 2.0 * Math.PI * UniformRNG.r8_uniform_01 ( ref seed );
            }

            for ( j = 0; j < n; j++ )
            {
                x[0+j*DIM_NUM] = r[j] * Math.Cos ( t[j] );
                x[1+j*DIM_NUM] = r[j] * Math.Sin ( t[j] );
            }

            return x;
        }
    }
}