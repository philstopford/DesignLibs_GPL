namespace Burkardt.Values
{
    public static class Legendre
    {
        public static void lp_values(ref int n_data, ref int n, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LP_VALUES returns values of the Legendre polynomials P(n,x).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 September 2014
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
            int N_MAX = 22;

            double[] fx_vec =
            {
                0.1000000000000000E+01,
                0.2500000000000000E+00,
                -0.4062500000000000E+00,
                -0.3359375000000000E+00,
                0.1577148437500000E+00,
                0.3397216796875000E+00,
                0.2427673339843750E-01,
                -0.2799186706542969E+00,
                -0.1524540185928345E+00,
                0.1768244206905365E+00,
                0.2212002165615559E+00,
                0.0000000000000000E+00,
                -0.1475000000000000E+00,
                -0.2800000000000000E+00,
                -0.3825000000000000E+00,
                -0.4400000000000000E+00,
                -0.4375000000000000E+00,
                -0.3600000000000000E+00,
                -0.1925000000000000E+00,
                0.8000000000000000E-01,
                0.4725000000000000E+00,
                0.1000000000000000E+01
            };

            int[] n_vec =
            {
                0, 1, 2,
                3, 4, 5,
                6, 7, 8,
                9, 10, 3,
                3, 3, 3,
                3, 3, 3,
                3, 3, 3,
                3
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
                0.25E+00,
                0.25E+00,
                0.25E+00,
                0.25E+00,
                0.00E+00,
                0.10E+00,
                0.20E+00,
                0.30E+00,
                0.40E+00,
                0.50E+00,
                0.60E+00,
                0.70E+00,
                0.80E+00,
                0.90E+00,
                1.00E+00
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

        public static void p_polynomial_values(ref int n_data, ref int n, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P_POLYNOMIAL_VALUES: selected values of the Legendre polynomials P(n,x).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 March 2012
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
            int N_MAX = 22;

            double[] fx_vec =
            {
                0.1000000000000000E+01,
                0.2500000000000000E+00,
                -0.4062500000000000E+00,
                -0.3359375000000000E+00,
                0.1577148437500000E+00,
                0.3397216796875000E+00,
                0.2427673339843750E-01,
                -0.2799186706542969E+00,
                -0.1524540185928345E+00,
                0.1768244206905365E+00,
                0.2212002165615559E+00,
                0.0000000000000000E+00,
                -0.1475000000000000E+00,
                -0.2800000000000000E+00,
                -0.3825000000000000E+00,
                -0.4400000000000000E+00,
                -0.4375000000000000E+00,
                -0.3600000000000000E+00,
                -0.1925000000000000E+00,
                0.8000000000000000E-01,
                0.4725000000000000E+00,
                0.1000000000000000E+01
            };

            int[] n_vec =
            {
                0, 1, 2,
                3, 4, 5,
                6, 7, 8,
                9, 10, 3,
                3, 3, 3,
                3, 3, 3,
                3, 3, 3,
                3
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
                0.25E+00,
                0.25E+00,
                0.25E+00,
                0.25E+00,
                0.00E+00,
                0.10E+00,
                0.20E+00,
                0.30E+00,
                0.40E+00,
                0.50E+00,
                0.60E+00,
                0.70E+00,
                0.80E+00,
                0.90E+00,
                1.00E+00
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
        public static void pm_polynomial_values(ref int n_data, ref int n, ref int m, ref double x, ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PM_POLYNOMIAL_VALUES: selected values of Legendre polynomials Pm(n,m,x).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 March 2012
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
        //    Output, int &N, int &M, double &X,
        //    the arguments of the function.
        //
        //    Output, double &FX, the value of the function.
        //
        {
            int N_MAX = 20;

            double[] fx_vec =
            {
                0.0000000000000000E+00,
                -0.5000000000000000E+00,
                0.0000000000000000E+00,
                0.3750000000000000E+00,
                0.0000000000000000E+00,
                -0.8660254037844386E+00,
                -0.1299038105676658E+01,
                -0.3247595264191645E+00,
                0.1353164693413185E+01,
                -0.2800000000000000E+00,
                0.1175755076535925E+01,
                0.2880000000000000E+01,
                -0.1410906091843111E+02,
                -0.3955078125000000E+01,
                -0.9997558593750000E+01,
                0.8265311444100484E+02,
                0.2024442836815152E+02,
                -0.4237997531890869E+03,
                0.1638320624828339E+04,
                -0.2025687389227225E+05
            }
            ;

            int[] m_vec =
            {
                0, 0, 0, 0,
                0, 1, 1, 1,
                1, 0, 1, 2,
                3, 2, 2, 3,
                3, 4, 4, 5
            }
            ;

            int[] n_vec =
            {
                1, 2, 3, 4,
                5, 1, 2, 3,
                4, 3, 3, 3,
                3, 4, 5, 6,
                7, 8, 9, 10
            }
            ;

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
        public static void pmn_polynomial_values(ref int n_data, ref int n, ref int m, ref double x,
        ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PMN_POLYNOMIAL_VALUES: selected values of normalized Legendre polynomial Pmn(n,m,x).
        //
        //  Discussion:
        //
        //    In Mathematica, the unnormalized function can be evaluated by:
        //
        //      LegendreP [ n, m, x ]
        //
        //    The function is normalized by dividing by
        //
        //      sqrt ( 2 * ( n + m )! / ( 2 * n + 1 ) / ( n - m )! )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 March 2012
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
        //    Output, int &N, int &M, double &X,
        //    the arguments of the function.
        //
        //    Output, double &FX, the value of the function.
        //
        {
            int N_MAX = 21;

            double[] fx_vec =
            {
                0.7071067811865475E+00,
                0.6123724356957945E+00,
                -0.7500000000000000E+00,
                -0.1976423537605237E+00,
                -0.8385254915624211E+00,
                0.7261843774138907E+00,
                -0.8184875533567997E+00,
                -0.1753901900050285E+00,
                0.9606516343087123E+00,
                -0.6792832849776299E+00,
                -0.6131941618102092E+00,
                0.6418623720763665E+00,
                0.4716705890038619E+00,
                -0.1018924927466445E+01,
                0.6239615396237876E+00,
                0.2107022704608181E+00,
                0.8256314721961969E+00,
                -0.3982651281554632E+00,
                -0.7040399320721435E+00,
                0.1034723155272289E+01,
                -0.5667412129155530E+00
            }
            ;

            int[] m_vec =
            {
                0, 0, 1, 0,
                1, 2, 0, 1,
                2, 3, 0, 1,
                2, 3, 4, 0,
                1, 2, 3, 4,
                5
            }
            ;

            int[] n_vec =
            {
                0, 1, 1, 2,
                2, 2, 3, 3,
                3, 3, 4, 4,
                4, 4, 4, 5,
                5, 5, 5, 5,
                5
            }
            ;

            double[] x_vec =
            {
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50
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
        public static void pmns_polynomial_values(ref int n_data, ref int n, ref int m, ref double x,
        ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PMNS_POLYNOMIAL_VALUES: selected values of sphere-normalized Legendre polynomial Pmns(n,m,x).
        //
        //  Discussion:
        //
        //    In Mathematica, the unnormalized function can be evaluated by:
        //
        //      LegendreP [ n, m, x ]
        //
        //    The function is normalized for the sphere by dividing by
        //
        //      sqrt ( 4 * Math.PI * ( n + m )! / ( 2 * n + 1 ) / ( n - m )! )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 September 2010
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
        //    Output, int &N, int &M, double &X,
        //    the arguments of the function.
        //
        //    Output, double &FX, the value of the function.
        //
        {
            int N_MAX = 21;

            double[] fx_vec =
            {
                0.2820947917738781,
                0.2443012559514600,
                -0.2992067103010745,
                -0.07884789131313000,
                -0.3345232717786446,
                0.2897056515173922,
                -0.3265292910163510,
                -0.06997056236064664,
                0.3832445536624809,
                -0.2709948227475519,
                -0.2446290772414100,
                0.2560660384200185,
                0.1881693403754876,
                -0.4064922341213279,
                0.2489246395003027,
                0.08405804426339821,
                0.3293793022891428,
                -0.1588847984307093,
                -0.2808712959945307,
                0.4127948151484925,
                -0.2260970318780046
            }
            ;

            int[] m_vec =
            {
                0, 0, 1, 0,
                1, 2, 0, 1,
                2, 3, 0, 1,
                2, 3, 4, 0,
                1, 2, 3, 4,
                5
            }
            ;

            int[] n_vec =
            {
                0, 1, 1, 2,
                2, 2, 3, 3,
                3, 3, 4, 4,
                4, 4, 4, 5,
                5, 5, 5, 5,
                5
            }
            ;

            double[] x_vec =
            {
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50
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

        public static void pn_polynomial_values(ref int n_data, ref int n, ref double x,
        ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PN_POLYNOMIAL_VALUES: selected values of the normalized Legendre polynomials.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 March 2016
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
            int N_MAX = 22;

            double[] fx_vec =
            {
                0.7071067811865475,
                0.3061862178478972,
                -0.642337649721702,
                -0.6284815141846855,
                0.3345637065282053,
                0.7967179601799685,
                0.06189376866246124,
                -0.766588850921089,
                -0.4444760242953344,
                0.5450094674858101,
                0.7167706229835538,
                0.0000000000000000,
                -0.2759472322745781,
                -0.5238320341483518,
                -0.7155919752205163,
                -0.823164625090267,
                -0.8184875533567997,
                -0.6734983296193094,
                -0.360134523476992,
                0.1496662954709581,
                0.8839665576253438,
                1.870828693386971
            }
            ;

            int[] n_vec =
            {
                0, 1, 2,
                3, 4, 5,
                6, 7, 8,
                9, 10, 3,
                3, 3, 3,
                3, 3, 3,
                3, 3, 3,
                3
            }
            ;

            double[] x_vec =
            {
                0.25E+00,
                0.25E+00,
                0.25E+00,
                0.25E+00,
                0.25E+00,
                0.25E+00,
                0.25E+00,
                0.25E+00,
                0.25E+00,
                0.25E+00,
                0.25E+00,
                0.00E+00,
                0.10E+00,
                0.20E+00,
                0.30E+00,
                0.40E+00,
                0.50E+00,
                0.60E+00,
                0.70E+00,
                0.80E+00,
                0.90E+00,
                1.00E+00
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

        public static void legendre_associated_values(ref int n_data, ref int n, ref int m, ref double x,
                ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEGENDRE_ASSOCIATED_VALUES returns values of associated Legendre functions.
            //
            //  Discussion:
            //
            //    The function considered is the associated Legendre polynomial P^M_N(X).
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      LegendreP [ n, m, x ]
            //
            //  Differential equation:
            //
            //    (1-X*X) * Y'' - 2 * X * Y + ( N (N+1) - (M*M/(1-X*X)) * Y = 0;
            //
            //  First terms:
            //
            //    M = 0  ( = Legendre polynomials of first kind P(N)(X) )
            //
            //    P00 =    1
            //    P10 =    1 X
            //    P20 = (  3 X^2 -   1)/2
            //    P30 = (  5 X^3 -   3 X)/2
            //    P40 = ( 35 X^4 -  30 X^2 +   3)/8
            //    P50 = ( 63 X^5 -  70 X^3 +  15 X)/8
            //    P60 = (231 X^6 - 315 X^4 + 105 X^2 -  5)/16
            //    P70 = (429 X^7 - 693 X^5 + 315 X^3 - 35 X)/16
            //
            //    M = 1
            //
            //    P01 =   0
            //    P11 =   1 * SQRT(1-X*X)
            //    P21 =   3 * SQRT(1-X*X) * X
            //    P31 = 1.5 * SQRT(1-X*X) * (5*X*X-1)
            //    P41 = 2.5 * SQRT(1-X*X) * (7*X*X*X-3*X)
            //
            //    M = 2
            //
            //    P02 =   0
            //    P12 =   0
            //    P22 =   3 * (1-X*X)
            //    P32 =  15 * (1-X*X) * X
            //    P42 = 7.5 * (1-X*X) * (7*X*X-1)
            //
            //    M = 3
            //
            //    P03 =   0
            //    P13 =   0
            //    P23 =   0
            //    P33 =  15 * (1-X*X)^1.5
            //    P43 = 105 * (1-X*X)^1.5 * X
            //
            //    M = 4
            //
            //    P04 =   0
            //    P14 =   0
            //    P24 =   0
            //    P34 =   0
            //    P44 = 105 * (1-X*X)^2
            //
            //  Recursion:
            //
            //    if N < M:
            //      P(N,M) = 0;
            //    if N = M:
            //      P(N,M) = (2*M-1)!! * (1-X*X)^(M/2) where N!! means the product of
            //      all the odd integers less than or equal to N.
            //    if N = M+1:
            //      P(N,M) = X*(2*M+1)*P(M,M)
            //    if M+1 < N:
            //      P(N,M) = ( X*(2*N-1)*P(N-1,M) - (N+M-1)*P(N-2,M) )/(N-M)
            //
            //  Restrictions:
            //
            //    -1 <= X <= 1
            //     0 <= M <= N
            //
            //  Special values:
            //
            //    P(N,0)(X) = P(N)(X), that is, for M=0, the associated Legendre
            //    polynomial of the first kind equals the Legendre polynomial of the
            //    first kind.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    27 August 2004
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
            //    Output, ref int N, ref int M, ref double X,
            //    the arguments of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 20;

            double[] fx_vec =
            {
                0.0000000000000000E+00,
                -0.5000000000000000E+00,
                0.0000000000000000E+00,
                0.3750000000000000E+00,
                0.0000000000000000E+00,
                -0.8660254037844386E+00,
                -0.1299038105676658E+01,
                -0.3247595264191645E+00,
                0.1353164693413185E+01,
                -0.2800000000000000E+00,
                0.1175755076535925E+01,
                0.2880000000000000E+01,
                -0.1410906091843111E+02,
                -0.3955078125000000E+01,
                -0.9997558593750000E+01,
                0.8265311444100484E+02,
                0.2024442836815152E+02,
                -0.4237997531890869E+03,
                0.1638320624828339E+04,
                -0.2025687389227225E+05
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

        public static void legendre_associated_normalized_sphere_values(ref int n_data, ref int n, ref int m,
                ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEGENDRE_ASSOCIATED_NORMALIZED_SPHERE_VALUES: normalized associated Legendre.
            //
            //  Discussion:
            //
            //    The function considered is the associated Legendre polynomial P^M_N(X).
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      LegendreP [ n, m, x ]
            //
            //    The function is normalized for the sphere by dividing by
            //
            //      sqrt ( 4 * Math.PI * ( n + m )! / ( 2 * n + 1 ) / ( n - m )! )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 September 2010
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
            //    Output, ref int N, ref int M, ref double X,
            //    the arguments of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 21;

            double[] fx_vec =
            {
                0.2820947917738781,
                0.2443012559514600,
                -0.2992067103010745,
                -0.07884789131313000,
                -0.3345232717786446,
                0.2897056515173922,
                -0.3265292910163510,
                -0.06997056236064664,
                0.3832445536624809,
                -0.2709948227475519,
                -0.2446290772414100,
                0.2560660384200185,
                0.1881693403754876,
                -0.4064922341213279,
                0.2489246395003027,
                0.08405804426339821,
                0.3293793022891428,
                -0.1588847984307093,
                -0.2808712959945307,
                0.4127948151484925,
                -0.2260970318780046
            };

            int[] m_vec =
            {
                0, 0, 1, 0,
                1, 2, 0, 1,
                2, 3, 0, 1,
                2, 3, 4, 0,
                1, 2, 3, 4,
                5
            };

            int[] n_vec =
            {
                0, 1, 1, 2,
                2, 2, 3, 3,
                3, 3, 4, 4,
                4, 4, 4, 5,
                5, 5, 5, 5,
                5
            };

            double[] x_vec =
            {
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50
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

        public static void legendre_associated_normalized_values(ref int n_data, ref int n, ref int m,
                ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEGENDRE_ASSOCIATED_NORMALIZED_VALUES: normalized associated Legendre.
            //
            //  Discussion:
            //
            //    The function considered is the associated Legendre polynomial P^M_N(X).
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      LegendreP [ n, m, x ]
            //
            //    The function is normalized by dividing by
            //
            //      sqrt ( 2 * ( n + m )! / ( 2 * n + 1 ) / ( n - m )! )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 March 2012
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
            //    Output, ref int N, ref int M, ref double X,
            //    the arguments of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 21;

            double[] fx_vec =
            {
                0.7071067811865475E+00,
                0.6123724356957945E+00,
                -0.7500000000000000E+00,
                -0.1976423537605237E+00,
                -0.8385254915624211E+00,
                0.7261843774138907E+00,
                -0.8184875533567997E+00,
                -0.1753901900050285E+00,
                0.9606516343087123E+00,
                -0.6792832849776299E+00,
                -0.6131941618102092E+00,
                0.6418623720763665E+00,
                0.4716705890038619E+00,
                -0.1018924927466445E+01,
                0.6239615396237876E+00,
                0.2107022704608181E+00,
                0.8256314721961969E+00,
                -0.3982651281554632E+00,
                -0.7040399320721435E+00,
                0.1034723155272289E+01,
                -0.5667412129155530E+00
            };

            int[] m_vec =
            {
                0, 0, 1, 0,
                1, 2, 0, 1,
                2, 3, 0, 1,
                2, 3, 4, 0,
                1, 2, 3, 4,
                5
            };

            int[] n_vec =
            {
                0, 1, 1, 2,
                2, 2, 3, 3,
                3, 3, 4, 4,
                4, 4, 4, 5,
                5, 5, 5, 5,
                5
            };

            double[] x_vec =
            {
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50
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

        public static void legendre_function_q_values(ref int n_data, ref int n, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEGENDRE_FUNCTION_Q_VALUES returns values of the Legendre Q function.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      LegendreQ[n,x]
            //
            //  Differential equation:
            //
            //    (1-X*X) Y'' - 2 X Y' + N (N+1) = 0;
            //
            //  First terms:
            //
            //    Q(0)(X) = 0.5 * log((1+X)/(1-X))
            //    Q(1)(X) = Q(0)(X)*X - 1
            //    Q(2)(X) = Q(0)(X)*(3*X*X-1)/4 - 1.5*X
            //    Q(3)(X) = Q(0)(X)*(5*X*X*X-3*X)/4 - 2.5*X^2 + 2/3
            //    Q(4)(X) = Q(0)(X)*(35*X^4-30*X^2+3)/16 - 35/8 * X^3 + 55/24 * X
            //    Q(5)(X) = Q(0)(X)*(63*X^5-70*X^3+15*X)/16 - 63/8*X^4 + 49/8*X^2 - 8/15
            //
            //  Recursion:
            //
            //    Q(0) = 0.5 * log ( (1+X) / (1-X) )
            //    Q(1) = 0.5 * X * log ( (1+X) / (1-X) ) - 1.0
            //
            //    Q(N) = ( (2*N-1) * X * Q(N-1) - (N-1) * Q(N-2) ) / N
            //
            //  Restrictions:
            //
            //    -1 < X < 1
            //
            //  Special values:
            //
            //    Note that the Legendre function Q(N)(X) is equal to the
            //    associated Legendre function of the second kind,
            //    Q(N,M)(X) with M = 0.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 August 2004
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
            //    Output, ref int N, the order of the function.
            //
            //    Output, ref double X, the point where the function is evaluated.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 21;

            double[] fx_vec =
            {
                0.2554128118829953E+00,
                -0.9361467970292512E+00,
                -0.4787614548274669E+00,
                0.4246139251747229E+00,
                0.5448396833845414E+00,
                -0.9451328261673470E-01,
                -0.4973516573531213E+00,
                -0.1499018843853194E+00,
                0.3649161918783626E+00,
                0.3055676545072885E+00,
                -0.1832799367995643E+00,
                0.6666666666666667E+00,
                0.6268672028763330E+00,
                0.5099015515315237E+00,
                0.3232754180589764E+00,
                0.8026113738148187E-01,
                -0.1986547714794823E+00,
                -0.4828663183349136E+00,
                -0.7252886849144386E+00,
                -0.8454443502398846E+00,
                -0.6627096245052618E+00
            };

            int[] n_vec =
            {
                0, 1, 2,
                3, 4, 5,
                6, 7, 8,
                9, 10, 3,
                3, 3, 3,
                3, 3, 3,
                3, 3, 3
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
                0.25E+00,
                0.25E+00,
                0.25E+00,
                0.25E+00,
                0.00E+00,
                0.10E+00,
                0.20E+00,
                0.30E+00,
                0.40E+00,
                0.50E+00,
                0.60E+00,
                0.70E+00,
                0.80E+00,
                0.90E+00
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

        public static void legendre_normalized_polynomial_values(ref int n_data, ref int n, ref double x,
                ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEGENDRE_NORMALIZED_POLYNOMIAL_VALUES: the normalized Legendre polynomials.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 March 2016
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
            //    Output, ref int N, the order of the function.
            //
            //    Output, ref double X, the point where the function is evaluated.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 22;

            double[] fx_vec =
            {
                0.7071067811865475,
                0.3061862178478972,
                -0.642337649721702,
                -0.6284815141846855,
                0.3345637065282053,
                0.7967179601799685,
                0.06189376866246124,
                -0.766588850921089,
                -0.4444760242953344,
                0.5450094674858101,
                0.7167706229835538,
                0.0000000000000000,
                -0.2759472322745781,
                -0.5238320341483518,
                -0.7155919752205163,
                -0.823164625090267,
                -0.8184875533567997,
                -0.6734983296193094,
                -0.360134523476992,
                0.1496662954709581,
                0.8839665576253438,
                1.870828693386971
            };

            int[] n_vec =
            {
                0, 1, 2,
                3, 4, 5,
                6, 7, 8,
                9, 10, 3,
                3, 3, 3,
                3, 3, 3,
                3, 3, 3,
                3
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
                0.25E+00,
                0.25E+00,
                0.25E+00,
                0.25E+00,
                0.00E+00,
                0.10E+00,
                0.20E+00,
                0.30E+00,
                0.40E+00,
                0.50E+00,
                0.60E+00,
                0.70E+00,
                0.80E+00,
                0.90E+00,
                1.00E+00
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

        public static void legendre_polynomial_values(ref int n_data, ref int n, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEGENDRE_POLYNOMIAL_VALUES returns values of the Legendre polynomials.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      LegendreP [ n, x ]
            //
            //  Differential equation:
            //
            //    (1-X*X) * P(N)(X)'' - 2 * X * P(N)(X)' + N * (N+1) = 0;
            //
            //  First terms:
            //
            //    P( 0)(X) =       1
            //    P( 1)(X) =       1 X
            //    P( 2)(X) =  (    3 X^2 -       1)/2
            //    P( 3)(X) =  (    5 X^3 -     3 X)/2
            //    P( 4)(X) =  (   35 X^4 -    30 X^2 +     3)/8
            //    P( 5)(X) =  (   63 X^5 -    70 X^3 +    15 X)/8
            //    P( 6)(X) =  (  231 X^6 -   315 X^4 +   105 X^2 -     5)/16
            //    P( 7)(X) =  (  429 X^7 -   693 X^5 +   315 X^3 -    35 X)/16
            //    P( 8)(X) =  ( 6435 X^8 - 12012 X^6 +  6930 X^4 -  1260 X^2 +   35)/128
            //    P( 9)(X) =  (12155 X^9 - 25740 X^7 + 18018 X^5 -  4620 X^3 +  315 X)/128
            //    P(10)(X) =  (46189 X^10-109395 X^8 + 90090 X^6 - 30030 X^4 + 3465 X^2
            //                 -63 ) /256
            //
            //  Recursion:
            //
            //    P(0)(X) = 1
            //    P(1)(X) = X
            //    P(N)(X) = ( (2*N-1)*X*P(N-1)(X)-(N-1)*P(N-2)(X) ) / N
            //
            //    P'(0)(X) = 0;
            //    P'(1)(X) = 1
            //    P'(N)(X) = ( (2*N-1)*(P(N-1)(X)+X*P'(N-1)(X)-(N-1)*P'(N-2)(X) ) / N
            //
            //  Formula:
            //
            //    P(N)(X) = (1/2**N) * sum ( 0 <= M <= N/2 ) C(N,M) C(2N-2M,N) X^(N-2*M)
            //
            //  Orthogonality:
            //
            //    Integral ( -1 <= X <= 1 ) P(I)(X) * P(J)(X) dX
            //      = 0 if I =/= J
            //      = 2 / ( 2*I+1 ) if I = J.
            //
            //  Approximation:
            //
            //    A function F(X) defined on [-1,1] may be approximated by the series
            //
            //      C0*P(0)(X) + C1*P(1)(X) + ... + CN*P(N)(X)
            //
            //    where
            //
            //      C(I) = (2*I+1)/(2) * Integral ( -1 <= X <= 1 ) F(X) P(I)(X) dx.
            //
            //  Special values:
            //
            //    P(N)(1) = 1.
            //    P(N)(-1) = (-1)^N.
            //    | P(N)(X) | <= 1 in [-1,1].
            //
            //    P(N,0)(X) = P(N)(X), that is, for M=0, the associated Legendre
            //    function of the first kind and order N equals the Legendre polynomial
            //    of the first kind and order N.
            //
            //    The N zeroes of P(N)(X) are the abscissas used for Gauss-Legendre
            //    quadrature of the integral of a function F(X) with weight function 1
            //    over the interval [-1,1].
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 August 2004
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
            //    Output, ref int N, the order of the function.
            //
            //    Output, ref double X, the point where the function is evaluated.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 22;

            double[] fx_vec =
            {
                0.1000000000000000E+01,
                0.2500000000000000E+00,
                -0.4062500000000000E+00,
                -0.3359375000000000E+00,
                0.1577148437500000E+00,
                0.3397216796875000E+00,
                0.2427673339843750E-01,
                -0.2799186706542969E+00,
                -0.1524540185928345E+00,
                0.1768244206905365E+00,
                0.2212002165615559E+00,
                0.0000000000000000E+00,
                -0.1475000000000000E+00,
                -0.2800000000000000E+00,
                -0.3825000000000000E+00,
                -0.4400000000000000E+00,
                -0.4375000000000000E+00,
                -0.3600000000000000E+00,
                -0.1925000000000000E+00,
                0.8000000000000000E-01,
                0.4725000000000000E+00,
                0.1000000000000000E+01
            };

            int[] n_vec =
            {
                0, 1, 2,
                3, 4, 5,
                6, 7, 8,
                9, 10, 3,
                3, 3, 3,
                3, 3, 3,
                3, 3, 3,
                3
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
                0.25E+00,
                0.25E+00,
                0.25E+00,
                0.25E+00,
                0.00E+00,
                0.10E+00,
                0.20E+00,
                0.30E+00,
                0.40E+00,
                0.50E+00,
                0.60E+00,
                0.70E+00,
                0.80E+00,
                0.90E+00,
                1.00E+00
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

        public static void legendre_shifted_polynomial_values(ref int n_data, ref int n, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LEGENDRE_SHIFTED_POLYNOMIAL_VALUES: values of shifted Legendre polynomials.
            //
            //  Discussion:
            //
            //    If we denote the Legendre polynomial by P(n)(x), and the shifted 
            //    Legendre polynomial by P01(n)(x), then
            //
            //      P01(n)(x) = P(n)(2*x-1)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 March 2016
            //
            //  Author:
            //
            //    John Burkardt
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
            //    Output, ref double X, the point where the function is evaluated.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 22;

            double[] fx_vec =
            {
                0.1000000000000000E+01,
                0.2500000000000000E+00,
                -0.4062500000000000E+00,
                -0.3359375000000000E+00,
                0.1577148437500000E+00,
                0.3397216796875000E+00,
                0.2427673339843750E-01,
                -0.2799186706542969E+00,
                -0.1524540185928345E+00,
                0.1768244206905365E+00,
                0.2212002165615559E+00,
                0.0000000000000000E+00,
                -0.1475000000000000E+00,
                -0.2800000000000000E+00,
                -0.3825000000000000E+00,
                -0.4400000000000000E+00,
                -0.4375000000000000E+00,
                -0.3600000000000000E+00,
                -0.1925000000000000E+00,
                0.8000000000000000E-01,
                0.4725000000000000E+00,
                0.1000000000000000E+01
            };

            int[] n_vec =
            {
                0, 1, 2,
                3, 4, 5,
                6, 7, 8,
                9, 10, 3,
                3, 3, 3,
                3, 3, 3,
                3, 3, 3,
                3
            };

            double[] x_vec =
            {
                0.625E+00,
                0.625E+00,
                0.625E+00,
                0.625E+00,
                0.625E+00,
                0.625E+00,
                0.625E+00,
                0.625E+00,
                0.625E+00,
                0.625E+00,
                0.625E+00,
                0.50E+00,
                0.55E+00,
                0.60E+00,
                0.65E+00,
                0.70E+00,
                0.75E+00,
                0.80E+00,
                0.85E+00,
                0.90E+00,
                0.95E+00,
                1.00E+00
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