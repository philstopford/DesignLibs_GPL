namespace Burkardt.Values;

public static class AGM
{
    public static void agm_values(ref int n_data, ref double a, ref double b, ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    AGM_VALUES returns some values of the AGM.
        //
        //  Discussion:
        //
        //    The AGM is defined for nonnegative A and B.
        //
        //    The AGM of numbers A and B is defined by setting
        //
        //      A(0) = A,
        //      B(0) = B
        //
        //      A(N+1) = ( A(N) + B(N) ) / 2
        //      B(N+1) = sqrt ( A(N) * B(N) )
        //
        //    The two sequences both converge to AGM(A,B).
        //
        //    In Mathematica, the AGM can be evaluated by
        //
        //      ArithmeticGeometricMean [ a, b ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 February 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
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
        //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, ref double A, &B, the argument ofs the function.
        //
        //    Output, ref double FX, the value of the function.
        //
    {
        const int N_MAX = 15;

        double[] a_vec
                =
                {
                    22.0,
                    83.0,
                    42.0,
                    26.0,
                    4.0,
                    6.0,
                    40.0,
                    80.0,
                    90.0,
                    9.0,
                    53.0,
                    1.0,
                    1.0,
                    1.0,
                    1.5
                }
            ;
        double[] b_vec
                =
                {
                    96.0,
                    56.0,
                    7.0,
                    11.0,
                    63.0,
                    45.0,
                    75.0,
                    0.0,
                    35.0,
                    1.0,
                    53.0,
                    2.0,
                    4.0,
                    8.0,
                    8.0
                }
            ;
        double[]  fx_vec
                =
                {
                    52.274641198704240049,
                    68.836530059858524345,
                    20.659301196734009322,
                    17.696854873743648823,
                    23.867049721753300163,
                    20.717015982805991662,
                    56.127842255616681863,
                    0.000000000000000000,
                    59.269565081229636528,
                    3.9362355036495554780,
                    53.000000000000000000,
                    1.4567910310469068692,
                    2.2430285802876025701,
                    3.6157561775973627487,
                    4.0816924080221632670
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
            a = 0.0;
            b = 0.0;
            fx = 0.0;
        }
        else
        {
            a = a_vec[n_data - 1];
            b = b_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }
    }
}