namespace Burkardt.AppliedStatistics
{
    public static partial class Algorithms
    {

        public static void trigamma_values(ref int n_data, ref double x, ref double fx)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIGAMMA_VALUES returns some values of the TriGamma function.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      PolyGamma[1,x]
        //
        //    TriGamma(X) = d^2 ln ( Gamma ( X ) ) / d X^2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 September 2004
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
        //    Output, double *X, the argument of the function.
        //
        //    Output, double *FX, the value of the function.
        //
        {
            int N_MAX = 11;

            double[] fx_vec =  {
                0.1644934066848226E+01,
                0.1433299150792759E+01,
                0.1267377205423779E+01,
                0.1134253434996619E+01,
                0.1025356590529597E+01,
                0.9348022005446793E+00,
                0.8584318931245799E+00,
                0.7932328301639984E+00,
                0.7369741375017002E+00,
                0.6879720582426356E+00,
                0.6449340668482264E+00
            }
            ;

            double[] x_vec =  {
                1.0E+00,
                1.1E+00,
                1.2E+00,
                1.3E+00,
                1.4E+00,
                1.5E+00,
                1.6E+00,
                1.7E+00,
                1.8E+00,
                1.9E+00,
                2.0E+00
            }
            ;

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }
    }
}
