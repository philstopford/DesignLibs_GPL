using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.PyramidNS;

public static class MonteCarlo
{
    public static double pyramid01_integral(int[] expon )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PYRAMID01_INTEGRAL: monomial integral in a unit pyramid.
        //
        //  Discussion:
        //
        //    This function returns the value of the integral of X^ALPHA Y^BETA Z^GAMMA
        //    over the unit pyramid.
        //
        //    The integration region is:
        //
        //    - ( 1 - Z ) <= X <= 1 - Z
        //    - ( 1 - Z ) <= Y <= 1 - Z
        //              0 <= Z <= 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 March 2008
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
        //    Input, int EXPON[3], the exponents.
        //
        //    Output, double PYRAMID01_INTEGRAL, the integral of the monomial
        //    over the pyramid.
        //
    {
        double value = 0.0;

        switch (expon[0] % 2)
        {
            case 0 when expon[1] % 2 == 0:
            {
                int i_hi = 2 + expon[0] + expon[1];

                int i;
                for (i = 0; i <= i_hi; i++)
                {
                    value += typeMethods.r8_mop(i) * typeMethods.r8_choose(i_hi, i)
                             / (i + expon[2] + 1);
                }

                value = value
                    * 2.0 / (expon[0] + 1)
                    * 2.0 / (expon[1] + 1);
                break;
            }
        }

        return value;
    }

    public static double[] pyramid01_sample(int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PYRAMID01_SAMPLE: sample the unit pyramid.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of samples desired.
        //
        //    Input/output, int SEED, a seed for the random
        //    number generator.
        //
        //    Output, double PYRAMID01_SAMPLE[3*N], the sample values.
        //
    {
        int j;
        const int m = 3;
        const double one_third = 1.0 / 3.0;

        double[] x = UniformRNG.r8mat_uniform_01_new(m, n, ref seed);

        for (j = 0; j < n; j++)
        {
            x[2 + j * 3] = 1.0 - Math.Pow(x[2 + j * 3], one_third);
            x[1 + j * 3] = (1.0 - x[2 + j * 3]) * (2.0 * x[1 + j * 3] - 1.0);
            x[0 + j * 3] = (1.0 - x[2 + j * 3]) * (2.0 * x[0 + j * 3] - 1.0);
        }

        return x;
    }

    public static double pyramid01_volume()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PYRAMID01_VOLUME: volume of a unit pyramid with square base.
        //
        //  Discussion:
        //
        //    The volume of this unit pyramid is 4/3.
        //
        //    The integration region is:
        //
        //      - ( 1 - Z ) <= X <= 1 - Z
        //      - ( 1 - Z ) <= Y <= 1 - Z
        //                0 <= Z <= 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double PYRAMID01_VOLUME, the volume of the pyramid.
        //
    {
        const double volume = 4.0 / 3.0;

        return volume;
    }
}