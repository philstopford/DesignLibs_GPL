namespace Burkardt.Values
{
    public static class Hermite
    {

        public static void h_polynomial_values(ref int n_data, ref int n, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    H_POLYNOMIAL_VALUES: tabulated values of H(i,x).
            //
            //  Discussion:
            //
            //    H(i,x) is the physicist's Hermite polynomial of degree I.
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      HermiteH[n,x]
            //
            //  Differential equation:
            //
            //    Y'' - 2 X Y' + 2 N Y = 0;
            //
            //  First terms:
            //
            //      1
            //      2 X
            //      4 X^2     -  2
            //      8 X^3     - 12 X
            //     16 X^4     - 48 X^2     + 12
            //     32 X^5    - 160 X^3    + 120 X
            //     64 X^6    - 480 X^4    + 720 X^2    - 120
            //    128 X^7   - 1344 X^5   + 3360 X^3   - 1680 X
            //    256 X^8   - 3584 X^6  + 13440 X^4  - 13440 X^2   + 1680
            //    512 X^9   - 9216 X^7  + 48384 X^5  - 80640 X^3  + 30240 X
            //   1024 X^10 - 23040 X^8 + 161280 X^6 - 403200 X^4 + 302400 X^2 - 30240
            //
            //  Recursion:
            //
            //    H(0,X) = 1,
            //    H(1,X) = 2*X,
            //    H(N,X) = 2*X * H(N-1,X) - 2*(N-1) * H(N-2,X)
            //
            //  Norm:
            //
            //    Integral ( -oo < X < +oo ) exp ( - X^2 ) * H(N,X)^2 dX
            //    = sqrt ( PI ) * 2^N * N!
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 February 2012
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
            //    Output, int &N, the order of the polynomial.
            //
            //    Output, double &X, the point where the polynomial is evaluated.
            //
            //    Output, double &FX, the value of the function.
            //
        {
            int N_MAX = 18;

            double[] fx_vec =
            {
                0.1000000000000000E+01,
                0.1000000000000000E+02,
                0.9800000000000000E+02,
                0.9400000000000000E+03,
                0.8812000000000000E+04,
                0.8060000000000000E+05,
                0.7178800000000000E+06,
                0.6211600000000000E+07,
                0.5206568000000000E+08,
                0.4212712000000000E+09,
                0.3275529760000000E+10,
                0.2432987360000000E+11,
                0.1712370812800000E+12,
                0.0000000000000000E+00,
                0.4100000000000000E+02,
                -0.8000000000000000E+01,
                0.3816000000000000E+04,
                0.3041200000000000E+07
            };

            int[] n_vec =
            {
                0, 1, 2,
                3, 4, 5,
                6, 7, 8,
                9, 10, 11,
                12, 5, 5,
                5, 5, 5
            };

            double[] x_vec =
            {
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                0.0E+00,
                0.5E+00,
                1.0E+00,
                3.0E+00,
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

        public static void he_polynomial_values(ref int n_data, ref int n, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HE_POLYNOMIAL_VALUES: tabulated values of He(i,x).
            //
            //  Discussion:
            //
            //    He(i,x) represents the probabilist's Hermite polynomial.
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      He(n,x) = HermiteH[n,x/Sqrt[2]] / Sqrt [ 2^n ] 
            //
            //  First terms:
            //
            //   1
            //   X
            //   X^2  -  1
            //   X^3  -  3 X
            //   X^4  -  6 X^2 +   3
            //   X^5  - 10 X^3 +  15 X
            //   X^6  - 15 X^4 +  45 X^2 -   15
            //   X^7  - 21 X^5 + 105 X^3 -  105 X
            //   X^8  - 28 X^6 + 210 X^4 -  420 X^2 +  105
            //   X^9  - 36 X^7 + 378 X^5 - 1260 X^3 +  945 X
            //   X^10 - 45 X^8 + 630 X^6 - 3150 X^4 + 4725 X^2 - 945
            //
            //  Recursion:
            //
            //    He(0,X) = 1,
            //    He(1,X) = X,
            //    He(N,X) = X * He(N-1,X) - (N-1) * He(N-2,X)
            //
            //  Norm:
            //
            //    Integral ( -oo < X < +oo ) exp ( - 0.5 * X^2 ) * He(M,X) He(N,X) dX 
            //    = sqrt ( 2 * pi ) * N! * delta ( M, N )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 February 2012
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
            //    Output, int &N, the order of the polynomial.
            //
            //    Output, double &X, the point where the polynomial is evaluated.
            //
            //    Output, double &FX, the value of the function.
            //
        {
            int N_MAX = 18;

            double[] fx_vec
                    =
                    {
                        1.000000000000000E+00,
                        5.000000000000000E+00,
                        24.00000000000000E+00,
                        110.0000000000000E+00,
                        478.0000000000000E+00,
                        1950.000000000000E+00,
                        7360.000000000000E+00,
                        25100.00000000000E+00,
                        73980.00000000000E+00,
                        169100.0000000000E+00,
                        179680.0000000000E+00,
                        -792600.0000000000E+00,
                        -5939480.000000000E+00,
                        0.000000000000000E+00,
                        6.281250000000000E+00,
                        6.000000000000000E+00,
                        18.00000000000000E+00,
                        90150.00000000000E+00
                    }
                ;

            int[] n_vec
                    =
                    {
                        0, 1, 2,
                        3, 4, 5,
                        6, 7, 8,
                        9, 10, 11,
                        12, 5, 5,
                        5, 5, 5
                    }
                ;

            double[] x_vec
                    =
                    {
                        5.0E+00,
                        5.0E+00,
                        5.0E+00,
                        5.0E+00,
                        5.0E+00,
                        5.0E+00,
                        5.0E+00,
                        5.0E+00,
                        5.0E+00,
                        5.0E+00,
                        5.0E+00,
                        5.0E+00,
                        5.0E+00,
                        0.0E+00,
                        0.5E+00,
                        1.0E+00,
                        3.0E+00,
                        1.0E+01
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

                public static void hep_values(ref int n_data, ref int n, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HEP_VALUES returns values of the Hermite polynomials He(n,x).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 October 2014
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
            //    Output, int &N, the order of the function.
            //
            //    Output, double &X, the point where the function is evaluated.
            //
            //    Output, double &FX, the value of the function.
            //
        {
            int N_MAX = 18;

            double[] fx_vec =
            {
                1.000000000000000E+00,
                5.000000000000000E+00,
                24.00000000000000E+00,
                110.0000000000000E+00,
                478.0000000000000E+00,
                1950.000000000000E+00,
                7360.000000000000E+00,
                25100.00000000000E+00,
                73980.00000000000E+00,
                169100.0000000000E+00,
                179680.0000000000E+00,
                -792600.0000000000E+00,
                -5939480.000000000E+00,
                0.000000000000000E+00,
                6.281250000000000E+00,
                6.000000000000000E+00,
                18.00000000000000E+00,
                90150.00000000000E+00
            };

            int[] n_vec =
            {
                0, 1, 2,
                3, 4, 5,
                6, 7, 8,
                9, 10, 11,
                12, 5, 5,
                5, 5, 5
            };

            double[] x_vec =
            {
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                0.0E+00,
                0.5E+00,
                1.0E+00,
                3.0E+00,
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

        public static void hf_function_values(ref int n_data, ref int n, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HF_FUNCTION_VALUES: tabulated values of Hf(i,x).
            //
            //  Discussion:
            //
            //    Hf(i,x) represents the Hermite function of "degree" I.   
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Hf(n,x) = HermiteH[n,x] 
            //        * Exp [ -1/2 * x^2] / Sqrt [ 2^n * n! * Sqrt[Pi] ]
            //
            //    The Hermite functions are orthonormal:
            //
            //      Integral ( -oo < x < +oo ) Hf(m,x) Hf(n,x) dx = delta ( m, n )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 February 2012
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
            //    Output, int &N, the order of the polynomial.
            //
            //    Output, double &X, the point where the polynomial is evaluated.
            //
            //    Output, double &FX, the value of the function.
            //
        {
            int N_MAX = 23;

            double[] fx_vec =
            {
                0.7511255444649425E+00, 0.0000000000000000E+00, -0.5311259660135985E+00,
                0.0000000000000000E+00, 0.4599685791773266E+00, 0.0000000000000000E+00,
                0.4555806720113325E+00, 0.6442883651134752E+00, 0.3221441825567376E+00,
                -0.2630296236233334E+00, -0.4649750762925110E+00, -0.5881521185179581E-01,
                0.3905052515434106E+00, 0.2631861423064045E+00, -0.2336911435996523E+00,
                -0.3582973361472840E+00, 0.6146344487883041E-01, 0.3678312067984882E+00,
                0.9131969309166278E-01, 0.4385750950032321E+00, -0.2624689527931006E-01,
                0.5138426125477819E+00, 0.9355563118061758E-01
            };

            int[] n_vec =
            {
                0, 1, 2,
                3, 4, 5,
                0, 1, 2,
                3, 4, 5,
                6, 7, 8,
                9, 10, 11,
                12, 5, 5,
                5, 5
            };

            double[] x_vec =
            {
                0.0E+00, 0.0E+00, 0.0E+00,
                0.0E+00, 0.0E+00, 0.0E+00,
                1.0E+00, 1.0E+00, 1.0E+00,
                1.0E+00, 1.0E+00, 1.0E+00,
                1.0E+00, 1.0E+00, 1.0E+00,
                1.0E+00, 1.0E+00, 1.0E+00,
                1.0E+00, 0.5E+00, 2.0E+00,
                3.0E+00, 4.0E+00
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
        public static void hermite_function_values(ref int n_data, ref int n, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HERMITE_FUNCTION_VALUES returns some values of the Hermite function.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Hf(n,x) = HermiteH[n,x] 
            //        * Exp [ -1/2 * x^2] / Sqrt [ 2^n * n! * Sqrt[Pi] ]
            //
            //    The Hermite functions are orthonormal:
            //
            //      Integral ( -oo < x < +oo ) Hf(m,x) Hf(n,x) dx = delta ( m, n )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 February 2012
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
            int N_MAX = 23;

            double[] fx_vec =
            {
                0.7511255444649425E+00, 0.0000000000000000E+00, -0.5311259660135985E+00,
                0.0000000000000000E+00, 0.4599685791773266E+00, 0.0000000000000000E+00,
                0.4555806720113325E+00, 0.6442883651134752E+00, 0.3221441825567376E+00,
                -0.2630296236233334E+00, -0.4649750762925110E+00, -0.5881521185179581E-01,
                0.3905052515434106E+00, 0.2631861423064045E+00, -0.2336911435996523E+00,
                -0.3582973361472840E+00, 0.6146344487883041E-01, 0.3678312067984882E+00,
                0.9131969309166278E-01, 0.4385750950032321E+00, -0.2624689527931006E-01,
                0.5138426125477819E+00, 0.9355563118061758E-01
            };

            int[] n_vec =
            {
                0, 1, 2,
                3, 4, 5,
                0, 1, 2,
                3, 4, 5,
                6, 7, 8,
                9, 10, 11,
                12, 5, 5,
                5, 5
            };

            double[] x_vec =
            {
                0.0E+00, 0.0E+00, 0.0E+00,
                0.0E+00, 0.0E+00, 0.0E+00,
                1.0E+00, 1.0E+00, 1.0E+00,
                1.0E+00, 1.0E+00, 1.0E+00,
                1.0E+00, 1.0E+00, 1.0E+00,
                1.0E+00, 1.0E+00, 1.0E+00,
                1.0E+00, 0.5E+00, 2.0E+00,
                3.0E+00, 4.0E+00
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

        public static void hermite_poly_phys_values(ref int n_data, ref int n, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HERMITE_POLY_PHYS_VALUES returns some values of the physicist's Hermite polynomial.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      HermiteH[n,x]
            //
            //  Differential equation:
            //
            //    Y'' - 2 X Y' + 2 N Y = 0;
            //
            //  First terms:
            //
            //      1
            //      2 X
            //      4 X^2     -  2
            //      8 X^3     - 12 X
            //     16 X^4     - 48 X^2     + 12
            //     32 X^5    - 160 X^3    + 120 X
            //     64 X^6    - 480 X^4    + 720 X^2    - 120
            //    128 X^7   - 1344 X^5   + 3360 X^3   - 1680 X
            //    256 X^8   - 3584 X^6  + 13440 X^4  - 13440 X^2   + 1680
            //    512 X^9   - 9216 X^7  + 48384 X^5  - 80640 X^3  + 30240 X
            //   1024 X^10 - 23040 X^8 + 161280 X^6 - 403200 X^4 + 302400 X^2 - 30240
            //
            //  Recursion:
            //
            //    H(0,X) = 1,
            //    H(1,X) = 2*X,
            //    H(N,X) = 2*X * H(N-1,X) - 2*(N-1) * H(N-2,X)
            //
            //  Norm:
            //
            //    Integral ( -oo < X < +oo ) exp ( - X^2 ) * H(N,X)^2 dX
            //    = sqrt ( PI ) * 2^N * N!
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 February 2012
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
            int N_MAX = 18;

            double[] fx_vec =
            {
                0.1000000000000000E+01,
                0.1000000000000000E+02,
                0.9800000000000000E+02,
                0.9400000000000000E+03,
                0.8812000000000000E+04,
                0.8060000000000000E+05,
                0.7178800000000000E+06,
                0.6211600000000000E+07,
                0.5206568000000000E+08,
                0.4212712000000000E+09,
                0.3275529760000000E+10,
                0.2432987360000000E+11,
                0.1712370812800000E+12,
                0.0000000000000000E+00,
                0.4100000000000000E+02,
                -0.8000000000000000E+01,
                0.3816000000000000E+04,
                0.3041200000000000E+07
            };

            int[] n_vec =
            {
                0, 1, 2,
                3, 4, 5,
                6, 7, 8,
                9, 10, 11,
                12, 5, 5,
                5, 5, 5
            };

            double[] x_vec =
            {
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                0.0E+00,
                0.5E+00,
                1.0E+00,
                3.0E+00,
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

        public static void hermite_poly_prob_values(ref int n_data, ref int n, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HERMITE_POLY_PROB_VALUES: values of the probabilist's Hermite polynomial.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      He(n,x) = HermiteH[n,x/Sqrt[2]] / Sqrt [ 2^n ] 
            //
            //  First terms:
            //
            //   1
            //   X
            //   X^2  -  1
            //   X^3  -  3 X
            //   X^4  -  6 X^2 +   3
            //   X^5  - 10 X^3 +  15 X
            //   X^6  - 15 X^4 +  45 X^2 -   15
            //   X^7  - 21 X^5 + 105 X^3 -  105 X
            //   X^8  - 28 X^6 + 210 X^4 -  420 X^2 +  105
            //   X^9  - 36 X^7 + 378 X^5 - 1260 X^3 +  945 X
            //   X^10 - 45 X^8 + 630 X^6 - 3150 X^4 + 4725 X^2 - 945
            //
            //  Recursion:
            //
            //    He(0,X) = 1,
            //    He(1,X) = X,
            //    He(N,X) = X * He(N-1,X) - (N-1) * He(N-2,X)
            //
            //  Norm:
            //
            //    Integral ( -oo < X < +oo ) exp ( - 0.5 * X^2 ) * He(M,X) He(N,X) dX 
            //    = sqrt ( 2 * pi ) * N! * delta ( M, N )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 February 2012
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
            int N_MAX = 18;

            double[] fx_vec =
            {
                1.000000000000000E+00,
                5.000000000000000E+00,
                24.00000000000000E+00,
                110.0000000000000E+00,
                478.0000000000000E+00,
                1950.000000000000E+00,
                7360.000000000000E+00,
                25100.00000000000E+00,
                73980.00000000000E+00,
                169100.0000000000E+00,
                179680.0000000000E+00,
                -792600.0000000000E+00,
                -5939480.000000000E+00,
                0.000000000000000E+00,
                6.281250000000000E+00,
                6.000000000000000E+00,
                18.00000000000000E+00,
                90150.00000000000E+00
            };

            int[] n_vec =
            {
                0, 1, 2,
                3, 4, 5,
                6, 7, 8,
                9, 10, 11,
                12, 5, 5,
                5, 5, 5
            };

            double[] x_vec =
            {
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                0.0E+00,
                0.5E+00,
                1.0E+00,
                3.0E+00,
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