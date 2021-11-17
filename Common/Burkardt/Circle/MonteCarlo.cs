using System;
using Burkardt.Uniform;

namespace Burkardt.CircleNS;

public static class MonteCarlo
{
    public static double[] circle01_sample_ergodic ( int n, ref double angle )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE01_SAMPLE_ERGODIC samples the circumference of the unit circle in 2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 June 2017
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
        //    Krieger, 1992,
        //    ISBN: 0894647644,
        //    LC: QA298.R79.
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input/output, double ANGLE, the base angle, which could be anything
        //    in the range [0,2 PI).
        //
        //    Output, double X[2*N], the points.
        //
    {
        double[] c = { 0.0, 0.0 };
        const double r = 1.0;

        double golden_ratio = ( 1.0 + Math.Sqrt ( 5.0 ) ) / 2.0;

        double golden_angle = 2.0 * Math.PI / Math.Pow ( golden_ratio, 2 );

        double[] x = new double[2 * n];

        for (int j = 0; j < n; j++ )
        {
            x[0+j*2] = c[0] + r * Math.Cos ( angle );
            x[1+j*2] = c[1] + r * Math.Sin ( angle );
            angle += golden_angle % 2.0 * Math.PI ;
        }

        return x;
    }
        
    public static double[] circle01_sample_random ( int n, ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE01_SAMPLE_RANDOM samples the circumference of the unit circle in 2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 June 2017
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
        //    Krieger, 1992,
        //    ISBN: 0894647644,
        //    LC: QA298.R79.
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input/output, int &SEED, a seed for the random 
        //    number generator.
        //
        //    Output, double X[2*N], the points.
        //
    {
        double[] c = { 0.0, 0.0 };
        const double r = 1.0;

        double[] theta = UniformRNG.r8vec_uniform_01_new ( n, ref seed );

        double[] x = new double[2*n];

        for ( int j = 0; j < n; j++ )
        {
            x[0+j*2] = c[0] + r * Math.Cos ( 2.0 * Math.PI * theta[j] );
            x[1+j*2] = c[1] + r * Math.Sin ( 2.0 * Math.PI * theta[j] );
        }

        return x;
    }
        
}