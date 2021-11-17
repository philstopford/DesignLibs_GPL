using System;
using Burkardt.Uniform;

namespace Burkardt.HyperGeometry.Hypercube;

public static class MonteCarlo
{
        
    public static double hypercube01_monomial_integral ( int m, int[] e )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPERCUBE01_MONOMIAL_INTEGRAL: integrals in unit hypercube in M dimensions.
        //
        //  Discussion:
        //
        //    The integration region is 
        //
        //      0 <= X(1:M) <= 1.
        //
        //    The monomial is F(X) = product ( 1 <= I <= M ) X(I)^E(I).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 January 2014
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
        //    Academic Press, 1984, page 263.
        //
        //  Parameters:
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, int E[M], the exponents.  
        //    Each exponent must be nonnegative.
        //
        //    Output, double HYPERCUBE01_MONOMIAL_INTEGRAL, the integral.
        //
    {
        double integral;

        for (int i = 0; i < m; i++ )
        {
            switch (e[i])
            {
                case < 0:
                    Console.WriteLine();
                    Console.WriteLine("HYPERCUBE01_MONOMIAL_INTEGRAL - Fatal error!");
                    Console.WriteLine("  All exponents must be nonnegative.");
                    Console.WriteLine("  E[" + i + "] = " + e[i]);
                    return 0;
            }
        }

        integral = 1.0;
        for (int i = 0; i < m; i++ )
        {
            integral /= ( e[i] + 1 );
        }

        return integral;
    }
        
        
    public static double[] hypercube01_sample ( int m, int n, ref int seed )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPERCUBE01_SAMPLE samples the unit hypercube in 3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 January 2014
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
        //    Input, int M, the spatial dimension.
        //
        //    Input, int N, the number of points.
        //
        //    Input/output, int &SEED, a seed for the random 
        //    number generator.
        //
        //    Output, double X[M*N], the points.
        //
    {
        double[] x = UniformRNG.r8mat_uniform_01 ( m, n, ref seed );

        return x;
    }
        
        
    public static double hypercube01_volume ( int m )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPERCUBE01_VOLUME: volume of the unit hypercube in M dimensions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the spatial dimension.
        //
        //    Output, double HYPERCUBE01_VOLUME, the volume.
        //
    {
        double volume = 1.0;

        return volume;
    }
}