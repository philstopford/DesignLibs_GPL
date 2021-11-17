using System;

namespace Burkardt.AppliedStatistics;

public static partial class Algorithms
{
    public static void chi_square_cdf_values(ref int n_data, ref int a, ref double x, ref double fx)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHI_SQUARE_CDF_VALUES returns some values of the Chi-Square CDF.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Needs["Statistics`ContinuousDistributions`"]
        //      dist = ChiSquareDistribution [ df ]
        //      CDF [ dist, x ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 August 2004
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
        //  Parameters:
        //
        //    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, int *A, the parameter of the function.
        //
        //    Output, double *X, the argument of the function.
        //
        //    Output, double *FX, the value of the function.
        //
    {
        const int N_MAX = 21;

        int[] a_vec =  {
                1, 2, 1, 2,
                1, 2, 3, 4,
                1, 2, 3, 4,
                5, 3, 3, 3,
                3, 3, 10, 10,
                10
            }
            ;

        double[] fx_vec =  {
                0.7965567455405796E-01,
                0.4987520807317687E-02,
                0.1124629160182849E+00,
                0.9950166250831946E-02,
                0.4729107431344619E+00,
                0.1812692469220181E+00,
                0.5975750516063926E-01,
                0.1752309630642177E-01,
                0.6826894921370859E+00,
                0.3934693402873666E+00,
                0.1987480430987992E+00,
                0.9020401043104986E-01,
                0.3743422675270363E-01,
                0.4275932955291202E+00,
                0.6083748237289110E+00,
                0.7385358700508894E+00,
                0.8282028557032669E+00,
                0.8883897749052874E+00,
                0.1721156299558408E-03,
                0.3659846827343712E-02,
                0.1857593622214067E-01
            }
            ;

        double[] x_vec =  {
                0.01E+00,
                0.01E+00,
                0.02E+00,
                0.02E+00,
                0.40E+00,
                0.40E+00,
                0.40E+00,
                0.40E+00,
                1.00E+00,
                1.00E+00,
                1.00E+00,
                1.00E+00,
                1.00E+00,
                2.00E+00,
                3.00E+00,
                4.00E+00,
                5.00E+00,
                6.00E+00,
                1.00E+00,
                2.00E+00,
                3.00E+00
            }
            ;

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        n_data += 1;

        if (N_MAX < n_data)
        {
            n_data = 0;
            a = 0;
            x = 0.0;
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