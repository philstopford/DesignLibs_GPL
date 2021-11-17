using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.SimplexNS;

public static class Integrals
{
    public static double simplex01_monomial_integral(int m, int[] e)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SIMPLEX01_MONOMIAL_INTEGRAL: integrals in the unit simplex in M dimensions.
        //
        //  Discussion:
        //
        //    The monomial is F(X) = product ( 1 <= I <= M ) X(I)^E(I).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, int E[M], the exponents.  
        //    Each exponent must be nonnegative.
        //
        //    Output, double SIMPLEX01_MONOMIAL_INTEGRAL, the integral.
        //
    {
        int i;
        double integral;
        int j;
        int k;

        for (i = 0; i < m; i++)
        {
            switch (e[i])
            {
                case < 0:
                    Console.WriteLine("");
                    Console.WriteLine("SIMPLEX01_MONOMIAL_INTEGRAL - Fatal error!");
                    Console.WriteLine("  All exponents must be nonnegative.");
                    return 1;
            }
        }

        k = 0;
        integral = 1.0;

        for (i = 0; i < m; i++)
        {
            for (j = 1; j <= e[i]; j++)
            {
                k += 1;
                integral = integral * j / k;
            }
        }

        for (i = 0; i < m; i++)
        {
            k += 1;
            integral /= k;
        }

        return integral;
    }

    public static double[] simplex01_sample(int m, int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SIMPLEX01_SAMPLE samples the unit simplex in M dimensions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 January 2015
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
        //    Input, int M, the spatial dimension.
        //
        //    Input, int N, the number of points.
        //
        //    Input/output, int &SEED, a seed for the random 
        //    number generator.
        //
        //    Output, double SIMPLEX01_SAMPLE_01[M*N], the points.
        //
    {
        double[] e;
        double e_sum;
        int i;
        int j;
        double[] x;

        x = new double[m * n];

        for (j = 0; j < n; j++)
        {
            e = UniformRNG.r8vec_uniform_01_new(m + 1, ref seed);

            for (i = 0; i < m + 1; i++)
            {
                e[i] = -Math.Log(e[i]);
            }

            e_sum = typeMethods.r8vec_sum(m + 1, e);

            for (i = 0; i < m; i++)
            {
                x[i + j * m] = e[i] / e_sum;
            }
        }

        return x;
    }

    public static double simplex01_volume(int m)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SIMPLEX01_VOLUME returns the volume of the unit simplex in M dimensions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the spatial dimension.
        //
        //    Output, double SIMPLEX01_VOLUME, the volume.
        //
    {
        int i;
        double volume;

        volume = 1.0;
        for (i = 1; i <= m; i++)
        {
            volume /= i;
        }

        return volume;
    }
}