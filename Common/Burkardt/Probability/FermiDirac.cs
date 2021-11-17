using System;
using Burkardt.Uniform;

namespace Burkardt.Probability;

public static class FermiDirac
{
    public static double fermi_dirac_sample(double u, double v, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FERMI_DIRAC_SAMPLE samples a (continuous) Fermi-Dirac distribution.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 April 2016
        //
        //  Author:
        //
        //    Original BASIC version by Frederick Ruckdeschel.
        //    C++ version by John Burkardt
        //
        //  Reference:
        //
        //    Frederick Ruckdeschel,
        //    BASIC Scientific Subroutines,
        //    Volume I,
        //    McGraw Hill, 1980,
        //    ISBN: 0-07-054201-5,
        //    LC: QA76.95.R82.
        //
        //  Parameters:
        //
        //    Input, double U, V, the parameters of the distribution.
        //    The value of U represents the halfway point for the distribution.
        //    Half the probability is to the left, and half to the right, of
        //    the value U.  The value of V controls the shape of the distribution.
        //    The ratio U/V determines the relative shape of the distribution.
        //    Values of U/V in excess of 100 will risk overflow.
        //
        //    Input/output, int &SEED, a seed for the random
        //    number generator.
        //
        //    Output, double FERMI_DIRAC_SAMPLE, a sample from the Fermi-Dirac distribution.
        //    Output values will be nonnegative, and roughly half of them should
        //    be less than or equal to U.
        //
    {
        int iter_max = 1000;
        double y1;

        double x = UniformRNG.r8_uniform_01(ref seed);
        double y = 1.0;
        double a = Math.Exp(4.0 * u / v);
        double b = (x - 1.0) * Math.Log(1.0 + a);

        int iter_num = 0;

        for (;;)
        {
            y1 = b + Math.Log(a + Math.Exp(y));

            if (Math.Abs(y - y1) < 0.001)
            {
                break;
            }

            y = y1;

            iter_num += 1;

            if (iter_max < iter_num)
            {
                break;
            }
        }

        double z = v * y1 / 4.0;

        return z;
    }
}