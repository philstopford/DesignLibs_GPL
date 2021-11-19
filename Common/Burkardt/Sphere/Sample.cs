using System;
using Burkardt.Sampling;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.SphereNS;

public static class Sample
{
    public static GeometrySampleResult sample_sphere_uniform ( int m, int n, int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SAMPLE_SPHERE_UNIFORM samples points inside the unit sphere.
        //
        //  Discussion:
        //
        //    The sphere has center 0 and radius 1.
        //
        //    We first generate a point ON the sphere, and then distribute it
        //    IN the sphere.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Russell Cheng,
        //    Random Variate Generation,
        //    in Handbook of Simulation,
        //    edited by Jerry Banks,
        //    Wiley, 1998, pages 168.
        //
        //    Reuven Rubinstein,
        //    Monte Carlo Optimization, Simulation, and Sensitivity
        //    of Queueing Networks,
        //    Wiley, 1986, page 232.
        //
        //  Parameters:
        //
        //    Input, int M, the dimension of the space.
        //
        //    Input, int N, the number of points.
        //
        //    Input/output, int *SEED, a seed for the random number generator.
        //
        //    Output, double SAMPLE_SPHERE_UNIFORM[M*N], the points.
        //
    {
        GeometrySampleResult result = new();
        int j;

        double[] x = new double[m*n];
        double[] y = new double[m];

        double exponent = 1.0 / m;

        for ( j = 0; j < n; j++ )
        {
            //
            //  Fill a vector with normally distributed values.
            //
            typeMethods.r8vec_normal_01 ( m, ref seed, ref y );
            //
            //  Compute the length of the vector.
            //
            double norm = 0.0;
            int i;
            for ( i = 0; i < m; i++ )
            {
                norm += y[i] * y[i];
            }
            norm = Math.Sqrt ( norm );
            //
            //  Normalize the vector.
            //
            for ( i = 0; i < m; i++ )
            {
                y[i] /= norm;
            }
            //
            //  Now compute a value to map the point ON the sphere INTO the sphere.
            //
            double r = UniformRNG.r8_uniform_01 ( ref seed );
            r = Math.Pow ( r, exponent );

            for ( i = 0; i < m; i++ )
            {
                x[i+j*m] = r * y[i];
            }
        }

        result.result = x;
        result.seed = seed;

        return result;
    }
}