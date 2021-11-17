using System;
using Burkardt.Uniform;

namespace Burkardt.IntegralNS;

public static class MonteCarlo
{
    public static double monte_carlo_nd ( Func<int, double[], double> func, int dim_num, 
            double[] a, double[] b, int eval_num, ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MONTE_CARLO_ND estimates a multidimensional integral using Monte Carlo.
        //
        //  Discussion:
        //
        //    Unlike the other routines, this routine requires the user to specify
        //    the number of function evaluations as an INPUT quantity.  
        //    No attempt at error estimation is made.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Philip Davis, Philip Rabinowitz,
        //    Methods of Numerical Integration,
        //    Second Edition,
        //    Dover, 2007,
        //    ISBN: 0486453391,
        //    LC: QA299.3.D28.
        //
        //  Parameters:
        //
        //    Input, double FUNC ( int dim_num, double x[] ), evaluates
        //    the function to be integrated.
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, double A[DIM_NUM], B[DIM_NUM], the integration limits.
        //
        //    Input, int EVAL_NUM, the number of function evaluations.
        //
        //    Input/output, int *SEED, a seed for the random number generator.
        //
        //    Output, double MONTE_CARLO_ND, the approximate value of the integral.
        //
    {
        int dim;
        int i;
        double result;
        double volume;
        double[] x;

        result = 0.0;

        for ( i = 0; i < eval_num; i++ )
        {
            x = UniformRNG.r8vec_uniform_01_new ( dim_num, ref seed );

            result += func ( dim_num, x );
        }

        volume = 1.0;
        for ( dim = 0; dim < dim_num; dim++ )
        {
            volume *= ( b[dim] - a[dim] );
        }

        result = result * volume / eval_num;

        return result;
    }
        
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
            quad += f ( dim_num, x );
        }

        quad /= m;

        return quad;
    }
}