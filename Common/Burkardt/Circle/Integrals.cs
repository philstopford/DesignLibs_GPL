using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.CircleNS
{
    public static class Integrals
    {
        public static double circle01_length()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CIRCLE01_LENGTH: length of the circumference of the unit circle in 2D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    11 January 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, double CIRCLE01_LENGTH, the length.
            //
        {
            const double r = 1.0;

            double length = 2.0 * Math.PI * r;

            return length;
        }

        public static double circle01_monomial_integral(int[] e )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE01_MONOMIAL_INTEGRAL returns monomial integrals on the unit circle.
        //
        //  Discussion:
        //
        //    The integration region is 
        //
        //      X^2 + Y^2 = 1.
        //
        //    The monomial is F(X,Y) = X^E(1) * Y^E(2).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 January 2014
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
        //    Input, int E[2], the exponents of X and Y in the 
        //    monomial.  Each exponent must be nonnegative.
        //
        //    Output, double CIRCLE01_MONOMIAL_INTEGRAL, the integral.
        //
        {
            double arg;
            int i;
            double integral;

            if (e[0] < 0 || e[1] < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("CIRCLE01_MONOMIAL_INTEGRAL - Fatal error!");
                Console.WriteLine("  All exponents must be nonnegative.");
                Console.WriteLine("  E[0] = " + e[0] + "");
                Console.WriteLine("  E[1] = " + e[1] + "");
                return (1);
            }

            if ((e[0] % 2) == 1 || (e[1] % 2) == 1)
            {
                integral = 0.0;
            }
            else
            {
                integral = 2.0;

                for (i = 0; i < 2; i++)
                {
                    arg = 0.5 * (double) (e[i] + 1);
                    integral = integral * typeMethods.r8_gamma(arg);
                }

                arg = 0.5 * (double) (e[0] + e[1] + 2);
                integral = integral / typeMethods.r8_gamma(arg);
            }

            return integral;
        }

        public static double[] circle01_sample(int n, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CIRCLE01_SAMPLE samples the circumference of the unit circle in 2D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    11 January 2014
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
            double[] c = {
                0.0, 0.0
            }
            ;
            int j;
            const double r = 1.0;
            double[] theta;
            double[] x;

            theta = UniformRNG.r8vec_uniform_01_new(n, ref seed);

            x = new double[2 * n];

            for (j = 0; j < n; j++)
            {
                x[0 + j * 2] = c[0] + r *  Math.Cos(Math.PI * theta[j]);
                x[1 + j * 2] = c[1] + r * Math.Sin(Math.PI * theta[j]);
            }

            return x;
        }
    }
}