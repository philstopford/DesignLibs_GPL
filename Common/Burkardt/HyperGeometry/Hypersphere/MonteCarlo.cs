using System;
using Burkardt.Types;

namespace Burkardt.HyperGeometry.Hypersphere
{
    public static class MonteCarlo
    {
        public static double hypersphere01_area(int m)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HYPERSPHERE01_AREA returns the area of the unit hypersphere.
            //
            //  Discussion:
            //
            //     M   Area
            //
            //     2    2        * PI
            //     3    4        * PI
            //     4  ( 2 /   1) * PI^2
            //     5  ( 8 /   3) * PI^2
            //     6  ( 1 /   1) * PI^3
            //     7  (16 /  15) * PI^3
            //     8  ( 1 /   3) * PI^4
            //     9  (32 / 105) * PI^4
            //    10  ( 1 /  12) * PI^5
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 January 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the spatial dimension.
            //
            //    Output, double HYPERSPHERE01_AREA, the area.
            //
        {
            double area;
            int i;
            int m_half;
            

            if ((m % 2) == 0)
            {
                m_half = m / 2;
                area = 2.0 * Math.Pow(Math.PI, m_half);
                for (i = 1; i <= m_half - 1; i++)
                {
                    area = area / (double) (i);
                }
            }
            else
            {
                m_half = (m - 1) / 2;
                area = Math.Pow(Math.PI, m_half) * Math.Pow(2.0, m);
                for (i = m_half + 1; i <= 2 * m_half; i++)
                {
                    area = area / (double) (i);
                }
            }

            return area;
        }

        public static double hypersphere01_monomial_integral(int m, int[] e)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HYPERSPHERE01_MONOMIAL_INTEGRAL: monomial integrals on unit hypersphere.
            //
            //  Discussion:
            //
            //    The integration region is 
            //
            //      sum ( 1 <= I <= M ) X(I)^2 = 1.
            //
            //    The monomial is F(X) = product ( 1 <= I <= M ) X(I)^E(I).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 January 2014
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
            //    Input, int E[M], the exponents.  Each exponent must be nonnegative.
            //
            //    Output, double HYPERSPHERE01_MONOMIAL_INTEGRAL, the integral.
            //
        {
            double arg;
            int i;
            double integral;

            for (i = 0; i < m; i++)
            {
                if (e[i] < 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("HYPERSPHERE01_MONOMIAL_INTEGRAL - Fatal error!");
                    Console.WriteLine("  All exponents must be nonnegative.");
                    Console.WriteLine("  E[" + i + "] = " + e[0] + "");
                    return (1);
                }
            }

            for (i = 0; i < m; i++)
            {
                if ((e[i] % 2) == 1)
                {
                    integral = 0.0;
                    return integral;
                }
            }

            integral = 2.0;

            for (i = 0; i < m; i++)
            {
                arg = 0.5 * (double) (e[i] + 1);
                integral = integral * typeMethods.r8_gamma(arg);
            }

            arg = 0.5 * (double) (typeMethods.i4vec_sum(m, e) + m);
            integral = integral / typeMethods.r8_gamma(arg);

            return integral;
        }

        public static double[] hypersphere01_sample(int m, int n, ref typeMethods.r8vecNormalData data, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HYPERSPHERE01_SAMPLE uniformly samples the surface of the unit hypersphere.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 January 2014
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
            int i;
            int j;
            double norm;
            double[] x;

            x = typeMethods.r8mat_normal_01_new(m, n, ref data, ref seed);

            for (j = 0; j < n; j++)
            {
                //
                //  Compute the length of the vector.
                //
                norm = 0.0;
                for (i = 0; i < m; i++)
                {
                    norm = norm + Math.Pow(x[i + j * m], 2);
                }

                norm = Math.Sqrt(norm);
                //
                //  Normalize the vector.
                //
                for (i = 0; i < m; i++)
                {
                    x[i + j * m] = x[i + j * m] / norm;
                }
            }

            return x;
        }
    }
}