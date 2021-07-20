namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double r8_rise(double x, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_RISE computes the rising factorial function [X]^N.
            //
            //  Discussion:
            //
            //    [X}^N = X * ( X + 1 ) * ( X + 2 ) * ... * ( X + N - 1 ).
            //
            //    Note that the number of ways of arranging N objects in M ordered
            //    boxes is [M}^N.  (Here, the ordering in each box matters).  Thus,
            //    2 objects in 2 boxes have the following 6 possible arrangements:
            //
            //      -/12, 1/2, 12/-, -/21, 2/1, 21/-.
            //
            //    Moreover, the number of non-decreasing maps from a set of
            //    N to a set of M ordered elements is [M]^N / N!.  Thus the set of
            //    nondecreasing maps from (1,2,3) to (a,b,c,d) is the 20 elements:
            //
            //      aaa, abb, acc, add, aab, abc, acd, aac, abd, aad
            //      bbb, bcc, bdd, bbc, bcd, bbd, ccc, cdd, ccd, ddd.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 May 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, the argument of the rising factorial function.
            //
            //    Input, int N, the order of the rising factorial function.
            //    If N = 0, RISE = 1, if N = 1, RISE = X.  Note that if N is
            //    negative, a "falling" factorial will be computed.
            //
            //    Output, double R8_RISE, the value of the rising factorial function.
            //
        {
            int i;
            double value;

            value = 1.0;

            if (0 < n)
            {
                for (i = 1; i <= n; i++)
                {
                    value = value * x;
                    x = x + 1.0;
                }
            }
            else if (n < 0)
            {
                for (i = -1; n <= i; i--)
                {
                    value = value * x;
                    x = x - 1.0;
                }
            }

            return value;
        }

        public static void r8_rise_values(ref int n_data, ref double x, ref int n, ref double f )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_RISE_VALUES returns some values of the rising factorial function.
        //
        //  Discussion:
        //
        //    Pochhammer(X,Y) = Gamma(X+Y) / Gamma(X)
        //
        //    For integer arguments, Pochhammer(M,N) = ( M + N - 1 )! / ( N - 1 )!
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Pochhammer[X,Y]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 December 2014
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
        //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, double &X, int &N, the arguments of the function.
        //
        //    Output, double &F, the value of the function.
        //
        {
            int N_MAX = 15;

            double[] f_vec =
            {
                1680.000000000000,
                1962.597656250000,
                2279.062500000000,
                2631.972656250000,
                3024.000000000000,
                1.000000000000000,
                7.500000000000000,
                63.75000000000000,
                605.6250000000000,
                6359.062500000000,
                73129.21875000000,
                914115.2343750000,
                1.234055566406250E+07,
                1.789380571289063E+08,
                2.773539885498047E+09
            }
            ;

            int[] n_vec =
            {
                4,
                4,
                4,
                4,
                4,
                0,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8,
                9
            }
            ;

            double[] x_vec =
            {
                5.00,
                5.25,
                5.50,
                5.75,
                6.00,
                7.50,
                7.50,
                7.50,
                7.50,
                7.50,
                7.50,
                7.50,
                7.50,
                7.50,
                7.50
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
                n = 0;
                f = 0.0;
            }
            else
            {
                x = x_vec[n_data - 1];
                n = n_vec[n_data - 1];
                f = f_vec[n_data - 1];
            }
        }

    }
}