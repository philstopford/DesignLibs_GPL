using System;
using Burkardt.Types;

namespace Burkardt.SphereNS
{
    public static class Integrals
    {
        public static double sphere01_area()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE01_AREA returns the area of the unit sphere.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 January 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, double SPHERE01_AREA, the area.
            //
        {
            double area;
            const double r = 1.0;
            const double r8_pi = 3.141592653589793;

            area = 4.0 * r8_pi * r * r;

            return area;
        }

        public static double sphere01_monomial_integral(int[] e)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE01_MONOMIAL_INTEGRAL: integrals on the surface of the unit sphere in 3D.
            //
            //  Discussion:
            //
            //    The integration region is 
            //
            //      X^2 + Y^2 + Z^2 = 1.
            //
            //    The monomial is F(X,Y,Z) = X^E(1) * Y^E(2) * Z^E(3).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 September 2010
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
            //    Input, int E[3], the exponents of X, Y and Z in the 
            //    monomial.  Each exponent must be nonnegative.
            //
            //    Output, double SPHERE01_MONOMIAL_INTEGRAL, the integral.
            //
        {
            int i;
            double integral;
            double r8_pi = 3.141592653589793;

            if (e[0] < 0 || e[1] < 0 || e[2] < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("SPHERE01_MONOMIAL_INTEGRAL - Fatal error!");
                Console.WriteLine("  All exponents must be nonnegative.");
                Console.WriteLine("  E[0] = " + e[0] + "");
                Console.WriteLine("  E[1] = " + e[1] + "");
                Console.WriteLine("  E[2] = " + e[2] + "");
                return (1);
            }

            if (e[0] == 0 && e[1] == 0 && e[2] == 0)
            {
                integral = 2.0 * Math.Sqrt(r8_pi * r8_pi * r8_pi) / typeMethods.r8_gamma(1.5);
            }
            else if ((e[0] % 2) == 1 ||
                     (e[1] % 2) == 1 ||
                     (e[2] % 2) == 1)
            {
                integral = 0.0;
            }
            else
            {
                integral = 2.0;

                for (i = 0; i < 3; i++)
                {
                    integral = integral * typeMethods.r8_gamma(0.5 * (double)(e[i] + 1));
                }

                integral = integral
                           / typeMethods.r8_gamma(0.5 * (double)(e[0] + e[1] + e[2] + 3));

            }

            return integral;
        }

        public static double[] sphere01_sample(int n, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE01_SAMPLE uniformly samples the surface of the unit sphere in 3D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    25 September 2010
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
            //    Output, double X[3*N], the points.
            //
        {
            int i;
            int j;
            double norm;
            double[] x;
            typeMethods.r8vecNormalData data = new typeMethods.r8vecNormalData();

            x = typeMethods.r8mat_normal_01_new(3, n, ref data, ref seed);

            for (j = 0; j < n; j++)
            {
                //
                //  Compute the length of the vector.
                //
                norm = Math.Sqrt(Math.Pow(x[0 + j * 3], 2)
                                 + Math.Pow(x[1 + j * 3], 2)
                                 + Math.Pow(x[2 + j * 3], 2));
                //
                //  Normalize the vector.
                //
                for (i = 0; i < 3; i++)
                {
                    x[i + j * 3] = x[i + j * 3] / norm;
                }
            }

            return x;
        }
    }
}