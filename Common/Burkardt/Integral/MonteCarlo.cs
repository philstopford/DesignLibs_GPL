using System;
using Burkardt.Uniform;

namespace Burkardt.IntegralNS
{
    public static class MonteCarlo
    {
        public static double monte_carlo ( int dim_num, int m, Func< int, double[], double > f,
        ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MONTE_CARLO applies a Monte Carlo integration rule.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 November 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Ian Sloan, Stephen Joe,
        //    Lattice Methods for Multiple Integration,
        //    Oxford, 1994,
        //    ISBN: 0198534728,
        //    LC: QA311.S56
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int M, the number of points to use.
        //
        //    Input, double F ( int DIM_NUM, double X[] ), the name of the 
        //    user-supplied routine which evaluates the function.
        //
        //    Input/output, int *SEED, a seed for the random number generator.
        //
        //    Output, double MONTE_CARLO, the estimated integral.
        //
        {
            int j;
            double quad;
            double[] x = null;

            quad = 0.0;

            for ( j = 1; j <= m; j++ )
            {
                UniformRNG.r8vec_uniform_01 ( dim_num, ref seed, ref x );
                quad = quad + f ( dim_num, x );
            }

            quad = quad / ( double ) m;

            return quad;
        }
    }
}