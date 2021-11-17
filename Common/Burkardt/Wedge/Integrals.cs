using System;
using Burkardt.Uniform;

namespace Burkardt.Wedge;

public static class Integrals
{
    public static double wedge01_integral(int[] e)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WEDGE01_INTEGRAL returns the integral of a monomial in the unit wedge in 3D.
        //
        //  Discussion:
        //
        //    This routine returns the integral of
        //
        //      product ( 1 <= I <= 3 ) X(I)^E(I)
        //
        //    over the unit wedge.
        //
        //    The integration region is:
        //
        //      0 <= X
        //      0 <= Y
        //      X + Y <= 1
        //      -1 <= Z <= 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Arthur Stroud,
        //    Approximate Calculation of Multiple Integrals,
        //    Prentice Hall, 1971,
        //    ISBN: 0130438936,
        //    LC: QA311.S85.
        //
        //  Parameters:
        //
        //    Input, int E[3], the exponents.
        //
        //    Output, double WEDGE01_INTEGRAL, the integral of the monomial.
        //
    {
        int i;

        double value = 1.0;

        int k = e[0];

        for (i = 1; i <= e[1]; i++)
        {
            k += 1;
            value = value * i / k;
        }

        k += 1;
        value /= k;

        k += 1;
        value /= k;
        switch (e[2])
        {
            //
            //  Now account for integration in Z.
            //
            case -1:
                Console.WriteLine("");
                Console.WriteLine("WEDGE01_INTEGRAL - Fatal error!");
                Console.WriteLine("  E(3) = -1 is not a legal input.");
                return 1;
        }

        value = (e[2] % 2) switch
        {
            1 => 0.0,
            _ => value * 2.0 / (e[2] + 1)
        };

        return value;
    }

    public static double[] wedge01_sample(int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //   WEDGE01_SAMPLE samples points uniformly from the unit wedge in 3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 August 2014
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
        //    Output, double WEDGE01_SAMPLE[3*N], the points.
        //
    {
        int j;

        double[] x = new double[3 * n];

        for (j = 0; j < n; j++)
        {
            double[] e = UniformRNG.r8vec_uniform_01_new(4, ref seed);

            int i;
            for (i = 0; i < 3; i++)
            {
                e[i] = -Math.Log(e[i]);
            }

            double e_sum = 0.0;
            for (i = 0; i < 3; i++)
            {
                e_sum += e[i];
            }

            x[0 + j * 3] = e[0] / e_sum;
            x[1 + j * 3] = e[1] / e_sum;
            x[2 + j * 3] = 2.0 * e[3] - 1.0;
        }

        return x;
    }

    public static double wedge01_volume()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WEDGE01_VOLUME returns the volume of the unit wedge in 3D.
        //
        //  Discussion:
        //
        //    The unit wedge is:
        //
        //      0 <= X
        //      0 <= Y
        //      X + Y <= 1
        //      -1 <= Z <= 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double WEDGE01_VOLUME, the volume of the unit wedge.
        //
    {
        const double value = 1.0;

        return value;
    }
}