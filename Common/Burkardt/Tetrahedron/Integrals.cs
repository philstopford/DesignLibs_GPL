using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.TetrahedronNS
{
    public static class Integrals
    {
        public static double tet01_monomial_integral ( int dim_num, int[] expon )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TET01_MONOMIAL_INTEGRAL integrates a monomial over the unit tetrahedron.
            //
            //  Discussion:
            //
            //    This routine evaluates a monomial of the form
            //
            //      product ( 1 <= dim <= dim_num ) x(dim)^expon(dim)
            //
            //    where the exponents are nonnegative integers.  Note that
            //    if the combination 0^0 is encountered, it should be treated
            //    as 1.
            //
            //    Integral ( over unit tetrahedron ) x^l y^m z^n dx dy 
            //    = l! * m! * n! / ( l + m + n + 3 )!
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 July 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, int EXPON[DIM_NUM], the exponents.
            //
            //    Output, double TET01_MONOMIAL_INTEGRAL, the value of the integral of the
            //    monomial.
            //
        {
            int i;
            int k;
            double value;
            //
            //  The first computation ends with VALUE = 1.0;
            //
            value = 1.0;

            k = 0;

            for ( i = 1; i <= expon[0]; i++ )
            {
                k = k + 1;
                //  value = value * ( double ) ( i ) / ( double ) ( k );
            }

            for ( i = 1; i <= expon[1]; i++ )
            {
                k = k + 1;
                value = value * ( double ) ( i ) / ( double ) ( k );
            }

            for ( i = 1; i <= expon[2]; i++ )
            {
                k = k + 1;
                value = value * ( double ) ( i ) / ( double ) ( k );
            }

            k = k + 1;
            value = value / ( double ) ( k );

            k = k + 1;
            value = value / ( double ) ( k );

            k = k + 1;
            value = value / ( double ) ( k );

            return value;
        }
        
        public static double tetrahedron01_monomial_integral ( int[] e )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TETRAHEDRON01_MONOMIAL_INTEGRAL: integrals in the unit tetrahedron in 3D.
            //
            //  Discussion:
            //
            //    The monomial is F(X,Y,Z) = X^E(1) * Y^E(2) * Z^E(3).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    15 January 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int E[3], the exponents.  
            //    Each exponent must be nonnegative.
            //
            //    Output, double TETRAHEDRON01_MONOMIAL_INTEGRAL, the integral.
            //
        {
            int i;
            double integral;
            int j;
            int k;
            const int m = 3;

            for ( i = 0; i < m; i++ )
            {
                if ( e[i] < 0 )
                {
                    Console.WriteLine("");
                    Console.WriteLine("TETRAHEDRON01_MONOMIAL_INTEGRAL - Fatal error!");
                    Console.WriteLine("  All exponents must be nonnegative.");
                    return ( 1 );
                }
            }

            k = 0;
            integral = 1.0;

            for ( i = 0; i < m; i++ )
            {
                for ( j = 1; j <= e[i]; j++ )
                {
                    k = k + 1;
                    integral = integral * ( double ) ( j ) / ( double ) ( k );
                }
            }

            for ( i = 0; i < m; i++ )
            {
                k = k + 1;
                integral = integral / ( double ) ( k );
            }

            return integral;
        }
        
        public static double[] tetrahedron01_sample ( int n, ref int seed )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TETRAHEDRON01_SAMPLE samples the unit tetrahedron in 3D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    15 January 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
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
            //    Output, double TETRAHEDRON01_SAMPLE_01[3*N], the points.
            //
        {
            double[] e;
            double e_sum;
            int i;
            int j;
            const int m = 3;
            double[] x;

            x = new double[m*n];

            for ( j = 0; j < n; j++ )
            {
                e = UniformRNG.r8vec_uniform_01_new ( m + 1, ref seed );

                for ( i = 0; i < m + 1; i++ )
                {
                    e[i] = - Math.Log ( e[i] );
                }
                e_sum = typeMethods.r8vec_sum ( m + 1, e );

                for ( i = 0; i < m; i++ )
                {
                    x[i+j*m] = e[i] / e_sum;
                }
            }

            return x;
        }
        
        public static double tetrahedron01_volume ( )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TETRAHEDRON01_VOLUME returns the volume of the unit tetrahedron in 3D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    15 January 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, double TETRAHEDRON01_VOLUME, the volume.
            //
        {
            double volume;

            volume = 1.0 / 6.0;

            return volume;
        }
    }
}