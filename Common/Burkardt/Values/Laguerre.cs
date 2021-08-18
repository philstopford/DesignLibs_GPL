namespace Burkardt.Values
{
    public static class Laguerre
    {

        public static void laguerre_associated_values(ref int n_data, ref int n, ref int m, ref double x,
                ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LAGUERRE_ASSOCIATED_VALUES returns some values of the associated Laguerre polynomials.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      LaguerreL[n,m,x]
            //
            //    The associated Laguerre polynomials may be generalized so that the
            //    parameter M is allowed to take on arbitrary nonint *values.
            //    The resulting function is known as the generalized Laguerre function.
            //
            //    The polynomials satisfy the differential equation:
            //
            //      X * Y'' + (M+1-X) * Y' + (N-M) * Y = 0;
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 August 2004
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
            //    Output, ref int N, the order of the function.
            //
            //    Output, ref int M, the parameter.
            //
            //    Output, ref double X, the point where the function is evaluated.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 20;

            double[] fx_vec =
            {
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1500000000000000E+01,
                0.1625000000000000E+01,
                0.1479166666666667E+01,
                0.1148437500000000E+01,
                0.4586666666666667E+00,
                0.2878666666666667E+01,
                0.8098666666666667E+01,
                0.1711866666666667E+02,
                0.1045328776041667E+02,
                0.1329019368489583E+02,
                0.5622453647189670E+02,
                0.7484729341779436E+02,
                0.3238912982762806E+03,
                0.4426100000097533E+03,
                0.1936876572288250E+04
            };

            int[] m_vec =
            {
                0, 0, 0, 0,
                0, 1, 1, 1,
                1, 0, 1, 2,
                3, 2, 2, 3,
                3, 4, 4, 5
            };

            int[] n_vec =
            {
                1, 2, 3, 4,
                5, 1, 2, 3,
                4, 3, 3, 3,
                3, 4, 5, 6,
                7, 8, 9, 10
            };

            double[] x_vec =
            {
                0.00E+00,
                0.00E+00,
                0.00E+00,
                0.00E+00,
                0.00E+00,
                0.50E+00,
                0.50E+00,
                0.50E+00,
                0.50E+00,
                0.20E+00,
                0.20E+00,
                0.20E+00,
                0.20E+00,
                0.25E+00,
                0.25E+00,
                0.25E+00,
                0.25E+00,
                0.25E+00,
                0.25E+00,
                0.25E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                n = 0;
                m = 0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                n = n_vec[n_data - 1];
                m = m_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void laguerre_general_values(ref int n_data, ref int n, ref double a, ref double x,
                ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LAGUERRE_GENERAL_VALUES returns some values of the generalized Laguerre function.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      LaguerreL[n,a,x]
            //
            //    The functions satisfy the following differential equation:
            //
            //      X * Y'' + (ALPHA+1-X) * Y' + N * Y = 0;
            //
            //    Function values can be generated by the recursion:
            //
            //      L(0,ALPHA)(X) = 1
            //      L(1,ALPHA)(X) = 1+ALPHA-X
            //
            //      L(N,ALPHA)(X) = ( (2*N-1+ALPHA-X) * L(N-1,ALPHA)(X)
            //                     - (N-1+ALPHA) * L(N-2,ALPHA)(X) ) / N
            //
            //    The parameter ALPHA is required to be greater than -1.
            //
            //    For ALPHA = 0, the generalized Laguerre function L(N,ALPHA)(X)
            //    is equal to the Laguerre polynomial L(N)(X).
            //
            //    For ALPHA integral, the generalized Laguerre function
            //    L(N,ALPHA)(X) equals the associated Laguerre polynomial L(N,ALPHA)(X).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    21 August 2004
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
            //    Output, ref int N, the order of the function.
            //
            //    Output, ref double A, the parameter.
            //
            //    Output, ref double X, the point where the function is evaluated.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 20;

            double[] a_vec =
            {
                0.00E+00,
                0.25E+00,
                0.50E+00,
                0.75E+00,
                1.50E+00,
                2.50E+00,
                5.00E+00,
                1.20E+00,
                1.20E+00,
                1.20E+00,
                1.20E+00,
                1.20E+00,
                1.20E+00,
                5.20E+00,
                5.20E+00,
                5.20E+00,
                5.20E+00,
                5.20E+00,
                5.20E+00,
                5.20E+00
            };

            double[] fx_vec =
            {
                0.3726399739583333E-01,
                0.3494791666666667E+00,
                0.8710042317708333E+00,
                0.1672395833333333E+01,
                0.6657625325520833E+01,
                0.2395726725260417E+02,
                0.2031344319661458E+03,
                0.1284193996800000E+02,
                0.5359924801587302E+01,
                0.9204589064126984E+00,
                -0.1341585114857143E+01,
                -0.2119726307555556E+01,
                -0.1959193658349206E+01,
                0.1000000000000000E+01,
                0.5450000000000000E+01,
                0.1720125000000000E+02,
                0.4110393750000000E+02,
                0.8239745859375000E+02,
                0.1460179186171875E+03,
                0.2359204608298828E+03
            };

            int[] n_vec =
            {
                5,
                5,
                5,
                5,
                5,
                5,
                5,
                8,
                8,
                8,
                8,
                8,
                8,
                0,
                1,
                2,
                3,
                4,
                5,
                6
            };

            double[] x_vec =
            {
                0.25E+00,
                0.25E+00,
                0.25E+00,
                0.25E+00,
                0.25E+00,
                0.25E+00,
                0.25E+00,
                0.00E+00,
                0.20E+00,
                0.40E+00,
                0.60E+00,
                0.80E+00,
                1.00E+00,
                0.75E+00,
                0.75E+00,
                0.75E+00,
                0.75E+00,
                0.75E+00,
                0.75E+00,
                0.75E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                n = 0;
                a = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                n = n_vec[n_data - 1];
                a = a_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void laguerre_polynomial_values(ref int n_data, ref int n, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LAGUERRE_POLYNOMIAL_VALUES returns some values of the Laguerre polynomial.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      LaguerreL[n,x]
            //
            //  Differential equation:
            //
            //    X * Y'' + (1-X) * Y' + N * Y = 0;
            //
            //  First terms:
            //
            //      1
            //     -X    +  1
            //   (  X^2 -  4 X     +  2 ) / 2
            //   ( -X^3 +  9 X^2 -  18 X    +    6 ) / 6
            //   (  X^4 - 16 X^3 +  72 X^2 -   96 X +      24 ) / 24
            //   ( -X^5 + 25 X^4 - 200 X^3 +  600 X^2 -  600 x    +  120 ) / 120
            //   (  X^6 - 36 X^5 + 450 X^4 - 2400 X^3 + 5400 X^2 - 4320 X + 720 ) / 720
            //   ( -X^7 + 49 X^6 - 882 X^5 + 7350 X^4 - 29400 X^3
            //      + 52920 X^2 - 35280 X + 5040 ) / 5040
            //
            //  Recursion:
            //
            //    L(0)(X) = 1,
            //    L(1)(X) = 1-X,
            //    N * L(N)(X) = (2*N-1-X) * L(N-1)(X) - (N-1) * L(N-2)(X)
            //
            //  Orthogonality:
            //
            //    Integral ( 0 <= X < +oo ) exp ( - X ) * L(N)(X) * L(M)(X) dX
            //    = 0 if N /= M
            //    = 1 if N == M
            //
            //  Special values:
            //
            //    L(N)(0) = 1.
            //
            //  Relations:
            //
            //    L(N)(X) = (-1)^N / N! * exp ( x ) * (d/dx)^n ( exp ( - x ) * x^n )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 August 2004
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
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref int N, the order of the polynomial.
            //
            //    Output, ref double X, the point where the polynomial is evaluated.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 17;

            double[] fx_vec =
            {
                0.1000000000000000E+01,
                0.0000000000000000E+00,
                -0.5000000000000000E+00,
                -0.6666666666666667E+00,
                -0.6250000000000000E+00,
                -0.4666666666666667E+00,
                -0.2569444444444444E+00,
                -0.4047619047619048E-01,
                0.1539930555555556E+00,
                0.3097442680776014E+00,
                0.4189459325396825E+00,
                0.4801341790925124E+00,
                0.4962122235082305E+00,
                -0.4455729166666667E+00,
                0.8500000000000000E+00,
                -0.3166666666666667E+01,
                0.3433333333333333E+02
            };

            int[] n_vec =
            {
                0, 1, 2,
                3, 4, 5,
                6, 7, 8,
                9, 10, 11,
                12, 5, 5,
                5, 5
            };

            double[] x_vec =
            {
                1.0E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00,
                0.5E+00,
                3.0E+00,
                5.0E+00,
                1.0E+01
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                n = 0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                n = n_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

    }
}