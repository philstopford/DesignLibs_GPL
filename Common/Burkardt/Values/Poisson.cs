namespace Burkardt.Values
{
    public static class Poisson
    {
        public static void poisson_cdf_values(ref int n_data, ref double a, ref int x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POISSON_CDF_VALUES returns some values of the Poisson CDF.
            //
            //  Discussion:
            //
            //    CDF(X)(A) is the probability of at most X successes in unit time,
            //    given that the expected mean number of successes is A.
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Needs["Statistics`DiscreteDistributions`]
            //      dist = PoissonDistribution [ a ]
            //      CDF [ dist, x ]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Milton Abramowitz, Irene Stegun,
            //    Handbook of Mathematical Functions,
            //    National Bureau of Standards, 1964,
            //    ISBN: 0-486-61272-4,
            //    LC: QA47.A34.
            //
            //    Stephen Wolfram,
            //    The Mathematica Book,
            //    Fourth Edition,
            //    Cambridge University Press, 1999,
            //    ISBN: 0-521-64314-7,
            //    LC: QA76.95.W65.
            //
            //    Daniel Zwillinger,
            //    CRC Standard Mathematical Tables and Formulae,
            //    30th Edition, CRC Press, 1996, pages 653-658.
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double A, the parameter of the function.
            //
            //    Output, int *X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 21;

            double[] a_vec =
            {
                0.02E+00,
                0.10E+00,
                0.10E+00,
                0.50E+00,
                0.50E+00,
                0.50E+00,
                1.00E+00,
                1.00E+00,
                1.00E+00,
                1.00E+00,
                2.00E+00,
                2.00E+00,
                2.00E+00,
                2.00E+00,
                5.00E+00,
                5.00E+00,
                5.00E+00,
                5.00E+00,
                5.00E+00,
                5.00E+00,
                5.00E+00
            };

            double[] fx_vec =
            {
                0.9801986733067553E+00,
                0.9048374180359596E+00,
                0.9953211598395555E+00,
                0.6065306597126334E+00,
                0.9097959895689501E+00,
                0.9856123220330293E+00,
                0.3678794411714423E+00,
                0.7357588823428846E+00,
                0.9196986029286058E+00,
                0.9810118431238462E+00,
                0.1353352832366127E+00,
                0.4060058497098381E+00,
                0.6766764161830635E+00,
                0.8571234604985470E+00,
                0.6737946999085467E-02,
                0.4042768199451280E-01,
                0.1246520194830811E+00,
                0.2650259152973617E+00,
                0.4404932850652124E+00,
                0.6159606548330631E+00,
                0.7621834629729387E+00
            };

            int[] x_vec =
            {
                0, 0, 1, 0,
                1, 2, 0, 1,
                2, 3, 0, 1,
                2, 3, 0, 1,
                2, 3, 4, 5,
                6
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                a = 0.0;
                x = 0;
                fx = 0.0;
            }
            else
            {
                a = a_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

    }
}