using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.HyperGeometry.Hyperball;

public static class MonteCarlo
{
    public static double hyperball01_monomial_integral(int m, int[] e)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPERBALL01_MONOMIAL_INTEGRAL: monomial integrals in the unit hyperball.
        //
        //  Discussion:
        //
        //    The integration region is 
        //
        //      sum ( 1 <= I <= M ) X(I)^2 <= 1.
        //
        //    The monomial is F(X) = product ( 1 <= I <= M ) X(I)^E(I).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 January 2014
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
        //    Output, double HYPERBALL01_MONOMIAL_INTEGRAL, the integral.
        //
    {
        double arg;
        int i;
        double integral;
        const double r = 1.0;

        for (i = 0; i < m; i++)
        {
            switch (e[i])
            {
                case < 0:
                    Console.WriteLine("");
                    Console.WriteLine("HYPERBALL01_MONOMIAL_INTEGRAL - Fatal error!");
                    Console.WriteLine("  All exponents must be nonnegative.");
                    Console.WriteLine("  E[" + i + "] = " + e[i] + "");
                    return 1;
            }
        }

        for (i = 0; i < m; i++)
        {
            switch (e[i] % 2)
            {
                case 1:
                    integral = 0.0;
                    return integral;
            }
        }

        integral = 2.0;

        for (i = 0; i < m; i++)
        {
            arg = 0.5 * (e[i] + 1);
            integral *= typeMethods.r8_gamma(arg);
        }

        double s = typeMethods.i4vec_sum(m, e) + m;
        arg = 0.5 * s;

        integral /= typeMethods.r8_gamma(arg);
        /*
        The surface integral is now adjusted to give the volume integral.
        */
        integral = integral * Math.Pow(r, s) / s;

        return integral;
    }

    public static double[] hyperball01_sample(int m, int n, ref typeMethods.r8vecNormalData data, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPERBALL01_SAMPLE uniformly samples points from the unit hyperball.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 January 2014
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
        int j;

        double exponent = 1.0 / m;

        double[] x = typeMethods.r8mat_normal_01_new(m, n, ref data, ref seed);

        for (j = 0; j < n; j++)
        {
            //
            //  Compute the length of the vector.
            //
            double norm = 0.0;
            int i;
            for (i = 0; i < m; i++)
            {
                norm += Math.Pow(x[i + j * m], 2);
            }

            norm = Math.Sqrt(norm);
            //
            //  Normalize the vector.
            //
            for (i = 0; i < m; i++)
            {
                x[i + j * m] /= norm;
            }

            //
            //  Now compute a value to map the point ON the sphere INTO the sphere.
            //
            double r = UniformRNG.r8_uniform_01(ref seed);

            for (i = 0; i < m; i++)
            {
                x[i + j * m] = Math.Pow(r, exponent) * x[i + j * m];
            }
        }

        return x;
    }

    public static double[] hyperball01_indicator ( int dim_num, int point_num, double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPERBALL01_INDICATOR evaluates the unit hyperball indicator function.
        //
        //  Discussion:
        //
        //    F(X) = 1 if X is on or inside the unit sphere, and 0 elsewhere.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int POINT_NUM, the number of points to evaluate.
        //
        //    Input, double X[DIM_NUM*POINT_NUM], the points.
        //
        //    Output, double SPHERE_INDICATOR[POINT_NUM], the indicator value.
        //
    {
        int j;

        double[] value = new double[point_num];

        for ( j = 0; j < point_num; j++ )
        {
            double t = 0.0;
            int i;
            for ( i = 0; i < dim_num; i++ )
            {
                t += x[i+j*dim_num] * x[i+j*dim_num];
            }

            value[j] = t switch
            {
                <= 1.0 => 1.0,
                _ => 0.0
            };
        }
        return value;
    }
    public static double hyperball01_volume(int m)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPERBALL01_VOLUME returns the volume of the unit hyperball.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the spatial dimension.
        //
        //    Output, double HYPERBALL01_VOLUME, the volume.
        //
    {
        int i;
        int m_half;
            
        double volume;

        switch (m % 2)
        {
            case 0:
            {
                m_half = m / 2;
                volume = Math.Pow(Math.PI, m_half);
                for (i = 1; i <= m_half; i++)
                {
                    volume /= i;
                }

                break;
            }
            default:
            {
                m_half = (m - 1) / 2;
                volume = Math.Pow(Math.PI, m_half) * Math.Pow(2.0, m);
                for (i = m_half + 1; i <= 2 * m_half + 1; i++)
                {
                    volume /= i;
                }

                break;
            }
        }

        return volume;
    }
}