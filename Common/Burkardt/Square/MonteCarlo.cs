using System;
using Burkardt.Uniform;

namespace Burkardt.Square
{
    public static class MonteCarlo
    {
        public static double square01_area()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SQUARE01_AREA: area of the unit square in 2D.
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
            //    Output, double SQUARE01_AREA, the area
            //
        {
            double area;

            area = 1.0;

            return area;
        }

        public static double square01_monomial_integral(int[] e)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SQUARE01_MONOMIAL_INTEGRAL: monomial integrals on the unit square in 2D.
            //
            //  Discussion:
            //
            //    The integration region is 
            //
            //      0 <= X <= 1,
            //      0 <= Y <= 1.
            //
            //    The monomial is F(X,Y) = X^E(1) * Y^E(2).
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
            //    Input, int E[2], the exponents.  
            //    Each exponent must be nonnegative.
            //
            //    Output, double SQUARE01_MONOMIAL_INTEGRAL, the integral.
            //
        {
            int i;
            double integral;
            const int m = 2;

            for (i = 0; i < m; i++)
            {
                if (e[i] < 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("SQUARE01_MONOMIAL_INTEGRAL - Fatal error!");
                    Console.WriteLine("  All exponents must be nonnegative.");
                    Console.WriteLine("  E[" + i + "] = " + e[i] + "");
                    return (1);
                }
            }

            integral = 1.0;
            for (i = 0; i < m; i++)
            {
                integral = integral / (double) (e[i] + 1);
            }

            return integral;
        }

        public static double[] square01_sample(int n, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SQUARE01_SAMPLE samples the unit square in 2D.
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
            //    Input, int N, the number of points.
            //
            //    Input/output, int &SEED, a seed for the random 
            //    number generator.
            //
            //    Output, double X[2*N], the points.
            //
        {
            const int m = 2;
            double[] x;

            x = UniformRNG.r8mat_uniform_01_new(m, n, ref seed);

            return x;
        }
    }
}