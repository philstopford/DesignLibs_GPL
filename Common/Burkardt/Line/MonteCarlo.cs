using System;
using Burkardt.Uniform;

namespace Burkardt.LineNS
{
    public static class MonteCarlo
    {
        public static double line01_length()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LINE01_LENGTH: length of the unit line in 1D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 January 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, double LINE01_LENGTH, the length.
            //
        {
            double length;

            length = 1.0;

            return length;
        }

        public static double line01_monomial_integral(int e)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LINE01_MONOMIAL_INTEGRAL returns monomial integrals on the unit line.
            //
            //  Discussion:
            //
            //    The integration region is 
            //
            //      0 <= X <= 1.
            //
            //    The monomial is F(X) = X^E.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 January 2014
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
            //    Input, int E, the exponent.  E must be nonnegative.
            //
            //    Output, double LINE01_MONOMIAL_INTEGRAL, the integral.
            //
        {
            double integral;

            if (e == -1)
            {
                Console.WriteLine("");
                Console.WriteLine("LINE01_MONOMIAL_INTEGRAL - Fatal error!");
                Console.WriteLine("  E = -1.");
                return 1;
            }

            integral = 1.0 / (double) (e + 1);

            return integral;
        }

        public static double[] line01_sample_ergodic(int n, ref double shift)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LINE01_SAMPLE_ERGODIC samples the unit line in 1D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 June 2017
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input/output, double &SHIFT, a value between 0 and 1.
            //
            //    Output, double X[N], the points.
            //
        {
            double golden;
            int j;
            double[] x;

            golden = (1.0 + Math.Sqrt(5.0)) / 2.0;

            x = new double[n];

            shift = (shift % 1.0);

            for (j = 0; j < n; j++)
            {
                x[j] = shift;
                shift = ((shift + golden) % 1.0);
            }

            return x;
        }

        public static double[] line01_sample_random(int n, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LINE01_SAMPLE_RANDOM samples the unit line in 1D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 January 2014
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
            //    Output, double X[N], the points.
            //
        {
            double[] x;

            x = UniformRNG.r8vec_uniform_01_new(n, ref seed);

            return x;
        }
    }
}