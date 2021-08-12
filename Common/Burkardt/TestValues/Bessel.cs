namespace Burkardt.TestValues
{
    public static class Bessel
    {
        public static void bessel_i0_int_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_I0_INT_VALUES returns some values of the Bessel I0 integral.
            //
            //  Discussion:
            //
            //    The function is defined by:
            //
            //      I0_INT(x) = Integral ( 0 <= t <= x ) I0(t) dt
            //
            //    The data was reported by McLeod.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    29 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Allan McLeod,
            //    Algorithm 757:
            //    MISCFUN: A software package to compute uncommon special functions,
            //    ACM Transactions on Mathematical Software,
            //    Volume 22, Number 3, September 1996, pages 288-301.
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 20;

            double[] fx_vec =
            {
                0.19531256208818052282E-02,
                -0.39062549670565734544E-02,
                0.62520348032546565850E-01,
                0.12516285581366971819E+00,
                -0.51051480879740303760E+00,
                0.10865210970235898158E+01,
                0.27750019054282535299E+01,
                -0.13775208868039716639E+02,
                0.46424372058106108576E+03,
                0.64111867658021584522E+07,
                -0.10414860803175857953E+08,
                0.44758598913855743089E+08,
                -0.11852985311558287888E+09,
                0.31430078220715992752E+09,
                -0.83440212900794309620E+09,
                0.22175367579074298261E+10,
                0.58991731842803636487E+10,
                -0.41857073244691522147E+11,
                0.79553885818472357663E+12,
                0.15089715082719201025E+17
            };

            double[] x_vec =
            {
                0.0019531250E+00,
                -0.0039062500E+00,
                0.0625000000E+00,
                0.1250000000E+00,
                -0.5000000000E+00,
                1.0000000000E+00,
                2.0000000000E+00,
                -4.0000000000E+00,
                8.0000000000E+00,
                18.0000000000E+00,
                -18.5000000000E+00,
                20.0000000000E+00,
                -21.0000000000E+00,
                22.0000000000E+00,
                -23.0000000000E+00,
                24.0000000000E+00,
                25.0000000000E+00,
                -27.0000000000E+00,
                30.0000000000E+00,
                40.0000000000E+00
            };

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

        public static void bessel_i0_spherical_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_I0_SPHERICAL_VALUES returns some values of the Spherical Bessel function i0.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Sqrt[Pi/(2*x)] * BesselI[1/2,x]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 January 2007
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
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 21;

            double[] fx_vec =
            {
                1.001667500198440E+00,
                1.006680012705470E+00,
                1.026880814507039E+00,
                1.061089303580402E+00,
                1.110132477734529E+00,
                1.175201193643801E+00,
                1.257884462843477E+00,
                1.360215358179667E+00,
                1.484729970750144E+00,
                1.634541271164267E+00,
                1.813430203923509E+00,
                2.025956895698133E+00,
                2.277595505698373E+00,
                2.574897010920645E+00,
                2.925685126512827E+00,
                3.339291642469967E+00,
                3.826838748926716E+00,
                4.401577467270101E+00,
                5.079293155726485E+00,
                5.878791279137455E+00,
                6.822479299281938E+00
            };

            double[] x_vec =
            {
                0.1E+00,
                0.2E+00,
                0.4E+00,
                0.6E+00,
                0.8E+00,
                1.0E+00,
                1.2E+00,
                1.4E+00,
                1.6E+00,
                1.8E+00,
                2.0E+00,
                2.2E+00,
                2.4E+00,
                2.6E+00,
                2.8E+00,
                3.0E+00,
                3.2E+00,
                3.4E+00,
                3.6E+00,
                3.8E+00,
                4.0E+00
            };

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

        public static void bessel_i0_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_I0_VALUES returns some values of the I0 Bessel function.
            //
            //  Discussion:
            //
            //    The modified Bessel functions In(Z) and Kn(Z) are solutions of
            //    the differential equation
            //
            //      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
            //
            //    The modified Bessel function I0(Z) corresponds to N = 0.
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      BesselI[0,x]
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
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 20;

            double[] fx_vec =
            {
                0.1000000000000000E+01,
                0.1010025027795146E+01,
                0.1040401782229341E+01,
                0.1092045364317340E+01,
                0.1166514922869803E+01,
                0.1266065877752008E+01,
                0.1393725584134064E+01,
                0.1553395099731217E+01,
                0.1749980639738909E+01,
                0.1989559356618051E+01,
                0.2279585302336067E+01,
                0.3289839144050123E+01,
                0.4880792585865024E+01,
                0.7378203432225480E+01,
                0.1130192195213633E+02,
                0.1748117185560928E+02,
                0.2723987182360445E+02,
                0.6723440697647798E+02,
                0.4275641157218048E+03,
                0.2815716628466254E+04
            };

            double[] x_vec =
            {
                0.00E+00,
                0.20E+00,
                0.40E+00,
                0.60E+00,
                0.80E+00,
                0.10E+01,
                0.12E+01,
                0.14E+01,
                0.16E+01,
                0.18E+01,
                0.20E+01,
                0.25E+01,
                0.30E+01,
                0.35E+01,
                0.40E+01,
                0.45E+01,
                0.50E+01,
                0.60E+01,
                0.80E+01,
                0.10E+02
            };

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

        public static void bessel_i1_spherical_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_I1_SPHERICAL_VALUES returns some values of the Spherical Bessel function i1.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Sqrt[Pi/(2*x)] * BesselI[3/2,x]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 January 2007
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
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 21;

            double[] fx_vec =
            {
                0.03336667857363341E+00,
                0.06693371456802954E+00,
                0.1354788933285401E+00,
                0.2072931911031093E+00,
                0.2841280857128948E+00,
                0.3678794411714423E+00,
                0.4606425870674146E+00,
                0.5647736480096238E+00,
                0.6829590627779635E+00,
                0.8182955028627777E+00,
                0.9743827435800610E+00,
                1.155432469636406E+00,
                1.366396525527973E+00,
                1.613118767572064E+00,
                1.902515460838681E+00,
                2.242790117769266E+00,
                2.643689828630357E+00,
                3.116811526884873E+00,
                3.675968313148932E+00,
                4.337627987747642E+00,
                5.121438384183637E+00
            };

            double[] x_vec =
            {
                0.1E+00,
                0.2E+00,
                0.4E+00,
                0.6E+00,
                0.8E+00,
                1.0E+00,
                1.2E+00,
                1.4E+00,
                1.6E+00,
                1.8E+00,
                2.0E+00,
                2.2E+00,
                2.4E+00,
                2.6E+00,
                2.8E+00,
                3.0E+00,
                3.2E+00,
                3.4E+00,
                3.6E+00,
                3.8E+00,
                4.0E+00
            };

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

        public static void bessel_i1_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_I1_VALUES returns some values of the I1 Bessel function.
            //
            //  Discussion:
            //
            //    The modified Bessel functions In(Z) and Kn(Z) are solutions of
            //    the differential equation
            //
            //      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      BesselI[1,x]
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
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 20;

            double[] fx_vec =
            {
                0.0000000000000000E+00,
                0.1005008340281251E+00,
                0.2040267557335706E+00,
                0.3137040256049221E+00,
                0.4328648026206398E+00,
                0.5651591039924850E+00,
                0.7146779415526431E+00,
                0.8860919814143274E+00,
                0.1084810635129880E+01,
                0.1317167230391899E+01,
                0.1590636854637329E+01,
                0.2516716245288698E+01,
                0.3953370217402609E+01,
                0.6205834922258365E+01,
                0.9759465153704450E+01,
                0.1538922275373592E+02,
                0.2433564214245053E+02,
                0.6134193677764024E+02,
                0.3998731367825601E+03,
                0.2670988303701255E+04
            };

            double[] x_vec =
            {
                0.00E+00,
                0.20E+00,
                0.40E+00,
                0.60E+00,
                0.80E+00,
                0.10E+01,
                0.12E+01,
                0.14E+01,
                0.16E+01,
                0.18E+01,
                0.20E+01,
                0.25E+01,
                0.30E+01,
                0.35E+01,
                0.40E+01,
                0.45E+01,
                0.50E+01,
                0.60E+01,
                0.80E+01,
                0.10E+02
            };

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

        public static void bessel_in_values(ref int n_data, ref int nu, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_IN_VALUES returns some values of the In Bessel function.
            //
            //  Discussion:
            //
            //    The modified Bessel functions In(Z) and Kn(Z) are solutions of
            //    the differential equation
            //
            //      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      BesselI[n,x]
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
            //    Output, ref int NU, the order of the function.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 28;

            double[] fx_vec =
            {
                0.5016687513894678E-02,
                0.1357476697670383E+00,
                0.6889484476987382E+00,
                0.1276466147819164E+01,
                0.2245212440929951E+01,
                0.1750561496662424E+02,
                0.2281518967726004E+04,
                0.3931278522104076E+08,
                0.2216842492433190E-01,
                0.2127399592398527E+00,
                0.1033115016915114E+02,
                0.1758380716610853E+04,
                0.2677764138883941E+21,
                0.2714631559569719E-03,
                0.9825679323131702E-02,
                0.2157974547322546E+01,
                0.7771882864032600E+03,
                0.2278548307911282E+21,
                0.2752948039836874E-09,
                0.3016963879350684E-06,
                0.4580044419176051E-02,
                0.2189170616372337E+02,
                0.1071597159477637E+21,
                0.3966835985819020E-24,
                0.4310560576109548E-18,
                0.5024239357971806E-10,
                0.1250799735644948E-03,
                0.5442008402752998E+19
            };

            int[] nu_vec =
            {
                2, 2, 2, 2,
                2, 2, 2, 2,
                3, 3, 3, 3,
                3, 5, 5, 5,
                5, 5, 10, 10,
                10, 10, 10, 20,
                20, 20, 20, 20
            };

            double[] x_vec =
            {
                0.2E+00,
                1.0E+00,
                2.0E+00,
                2.5E+00,
                3.0E+00,
                5.0E+00,
                10.0E+00,
                20.0E+00,
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00,
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00,
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00,
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                nu = 0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                nu = nu_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void bessel_ix_values(ref int n_data, ref double nu, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_IX_VALUES returns some values of the Ix Bessel function.
            //
            //  Discussion:
            //
            //    This set of data considers the less common case in which the
            //    index of the Bessel function In is actually not an integer.
            //    We may suggest this case by occasionally replacing the symbol
            //    "In" by "Ix".
            //
            //    The modified Bessel functions In(Z) and Kn(Z) are solutions of
            //    the differential equation
            //
            //      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      BesselI[n,x]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 March 2007
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
            //    Output, ref double NU, the order of the function.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 28;

            double[] fx_vec =
            {
                0.3592084175833614E+00,
                0.9376748882454876E+00,
                2.046236863089055E+00,
                3.053093538196718E+00,
                4.614822903407601E+00,
                26.47754749755907E+00,
                2778.784603874571E+00,
                4.327974627242893E+07,
                0.2935253263474798E+00,
                1.099473188633110E+00,
                21.18444226479414E+00,
                2500.906154942118E+00,
                2.866653715931464E+20,
                0.05709890920304825E+00,
                0.3970270801393905E+00,
                13.76688213868258E+00,
                2028.512757391936E+00,
                2.753157630035402E+20,
                0.4139416015642352E+00,
                1.340196758982897E+00,
                22.85715510364670E+00,
                2593.006763432002E+00,
                2.886630075077766E+20,
                0.03590910483251082E+00,
                0.2931108636266483E+00,
                11.99397010023068E+00,
                1894.575731562383E+00,
                2.716911375760483E+20
            };

            double[] nu_vec =
            {
                0.50E+00,
                0.50E+00,
                0.50E+00,
                0.50E+00,
                0.50E+00,
                0.50E+00,
                0.50E+00,
                0.50E+00,
                1.50E+00,
                1.50E+00,
                1.50E+00,
                1.50E+00,
                1.50E+00,
                2.50E+00,
                2.50E+00,
                2.50E+00,
                2.50E+00,
                2.50E+00,
                1.25E+00,
                1.25E+00,
                1.25E+00,
                1.25E+00,
                1.25E+00,
                2.75E+00,
                2.75E+00,
                2.75E+00,
                2.75E+00,
                2.75E+00
            };

            double[] x_vec =
            {
                0.2E+00,
                1.0E+00,
                2.0E+00,
                2.5E+00,
                3.0E+00,
                5.0E+00,
                10.0E+00,
                20.0E+00,
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00,
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00,
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00,
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                nu = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                nu = nu_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void bessel_j_spherical_values(ref int n_data, ref int n, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_J_SPHERICAL_VALUES returns values of the Spherical Bessel function j.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Sqrt[Pi/(2*x)] * BesselJ[n+1/2,x]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 January 2016
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
            //    Output, ref int N, the index of the function.

            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 22;

            double[] fx_vec =
            {
                0.8689780717709105E+00,
                0.2776712616989048E+00,
                0.05147914933043151E+00,
                0.006743927971987495E+00,
                0.0006838294584220406E+00,
                0.00005658597917091951E+00,
                3.955923765931341E-06,
                2.394450910776484E-07,
                1.277940110150618E-08,
                6.099572379372921E-10,
                2.633096568558721E-11,
                -0.05440211108893698E+00,
                0.07846694179875155E+00,
                0.07794219362856245E+00,
                -0.03949584498447032E+00,
                -0.1055892851176917E+00,
                -0.05553451162145218E+00,
                0.04450132233409427E+00,
                0.1133862306557747E+00,
                0.1255780236495678E+00,
                0.1000964095484906E+00,
                0.06460515449256426E+00
            };

            int[] n_vec =
            {
                0,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8,
                9,
                10,
                0,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8,
                9,
                10
            };

            double[] x_vec =
            {
                0.9050000000000000E+00,
                0.9050000000000000E+00,
                0.9050000000000000E+00,
                0.9050000000000000E+00,
                0.9050000000000000E+00,
                0.9050000000000000E+00,
                0.9050000000000000E+00,
                0.9050000000000000E+00,
                0.9050000000000000E+00,
                0.9050000000000000E+00,
                0.9050000000000000E+00,
                10.00000000000000E+00,
                10.00000000000000E+00,
                10.00000000000000E+00,
                10.00000000000000E+00,
                10.00000000000000E+00,
                10.00000000000000E+00,
                10.00000000000000E+00,
                10.00000000000000E+00,
                10.00000000000000E+00,
                10.00000000000000E+00,
                10.00000000000000E+00
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

        public static void bessel_j0_int_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_J0_INT_VALUES returns some values of the Bessel J0 integral.
            //
            //  Discussion:
            //
            //    The function is defined by:
            //
            //      J0_INT(x) = Integral ( 0 <= t <= x ) J0(t) dt
            //
            //    The data was reported by McLeod.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    29 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Allan McLeod,
            //    Algorithm 757:
            //    MISCFUN: A software package to compute uncommon special functions,
            //    ACM Transactions on Mathematical Software,
            //    Volume 22, Number 3, September 1996, pages 288-301.
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 20;

            double[] fx_vec =
            {
                0.97656242238978822427E-03,
                0.39062450329491108875E-02,
                -0.62479657927917933620E-01,
                0.12483733492120479139E+00,
                -0.48968050664604505505E+00,
                0.91973041008976023931E+00,
                -0.14257702931970265690E+01,
                0.10247341594606064818E+01,
                -0.12107468348304501655E+01,
                0.11008652032736190799E+01,
                -0.10060334829904124192E+01,
                0.81330572662485953519E+00,
                -0.10583788214211277585E+01,
                0.87101492116545875169E+00,
                -0.88424908882547488420E+00,
                0.11257761503599914603E+01,
                -0.90141212258183461184E+00,
                0.91441344369647797803E+00,
                -0.94482281938334394886E+00,
                0.92266255696016607257E+00
            };

            double[] x_vec =
            {
                0.0009765625E+00,
                0.0039062500E+00,
                -0.0625000000E+00,
                0.1250000000E+00,
                -0.5000000000E+00,
                1.0000000000E+00,
                -2.0000000000E+00,
                4.0000000000E+00,
                -8.0000000000E+00,
                16.0000000000E+00,
                -16.5000000000E+00,
                18.0000000000E+00,
                -20.0000000000E+00,
                25.0000000000E+00,
                -30.0000000000E+00,
                40.0000000000E+00,
                -50.0000000000E+00,
                75.0000000000E+00,
                -80.0000000000E+00,
                100.0000000000E+00
            };

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

        public static void bessel_j0_spherical_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_J0_SPHERICAL_VALUES returns some values of the Spherical Bessel function j0.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Sqrt[Pi/(2*x)] * BesselJ[1/2,x]
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
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 21;

            double[] fx_vec =
            {
                0.9983341664682815E+00,
                0.9933466539753061E+00,
                0.9735458557716262E+00,
                0.9410707889917256E+00,
                0.8966951136244035E+00,
                0.8414709848078965E+00,
                0.7766992383060220E+00,
                0.7038926642774716E+00,
                0.6247335019009407E+00,
                0.5410264615989973E+00,
                0.4546487134128408E+00,
                0.3674983653725410E+00,
                0.2814429918963129E+00,
                0.1982697583928709E+00,
                0.1196386250556803E+00,
                0.4704000268662241E-01,
                -0.1824191982111872E-01,
                -0.7515914765495039E-01,
                -0.1229223453596812E+00,
                -0.1610152344586103E+00,
                -0.1892006238269821E+00
            };

            double[] x_vec =
            {
                0.1E+00,
                0.2E+00,
                0.4E+00,
                0.6E+00,
                0.8E+00,
                1.0E+00,
                1.2E+00,
                1.4E+00,
                1.6E+00,
                1.8E+00,
                2.0E+00,
                2.2E+00,
                2.4E+00,
                2.6E+00,
                2.8E+00,
                3.0E+00,
                3.2E+00,
                3.4E+00,
                3.6E+00,
                3.8E+00,
                4.0E+00
            };

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

        public static void bessel_j0_zero_values(ref int n_data, ref int k, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_J0_ZERO_VALUES returns some values of zeroes of the J0 Bessel function.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 January 2016
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
            //    Output, ref double K, the index of the zero.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 30;

            double[] fx_vec =
            {
                2.40482555769577276862163187933E+00,
                5.52007811028631064959660411281E+00,
                8.65372791291101221695419871266E+00,
                11.7915344390142816137430449119E+00,
                14.9309177084877859477625939974E+00,
                18.0710639679109225431478829756E+00,
                21.2116366298792589590783933505E+00,
                24.3524715307493027370579447632E+00,
                27.4934791320402547958772882346E+00,
                30.6346064684319751175495789269E+00,
                33.7758202135735686842385463467E+00,
                36.9170983536640439797694930633E+00,
                40.0584257646282392947993073740E+00,
                43.1997917131767303575240727287E+00,
                46.3411883716618140186857888791E+00,
                49.4826098973978171736027615332E+00,
                52.6240518411149960292512853804E+00,
                55.7655107550199793116834927735E+00,
                58.9069839260809421328344066346E+00,
                62.0484691902271698828525002646E+00,
                65.1899648002068604406360337425E+00,
                68.3314693298567982709923038400E+00,
                71.4729816035937328250630738561E+00,
                74.6145006437018378838205404693E+00,
                77.7560256303880550377393718912E+00,
                80.8975558711376278637721434909E+00,
                84.0390907769381901578796383480E+00,
                87.1806298436411536512618050690E+00,
                90.3221726372104800557177667775E+00,
                93.4637187819447741711905915439E+00
            };

            int[] k_vec =
            {
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8,
                9,
                10,
                11,
                12,
                13,
                14,
                15,
                16,
                17,
                18,
                19,
                20,
                21,
                22,
                23,
                24,
                25,
                26,
                27,
                28,
                29,
                30
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                k = 0;
                fx = 0.0;
            }
            else
            {
                k = k_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void bessel_j0_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_J0_VALUES returns some values of the J0 Bessel function.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      BesselJ[0,x]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    10 August 2004
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
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 21;

            double[] fx_vec =
            {
                -0.1775967713143383E+00,
                -0.3971498098638474E+00,
                -0.2600519549019334E+00,
                0.2238907791412357E+00,
                0.7651976865579666E+00,
                0.1000000000000000E+01,
                0.7651976865579666E+00,
                0.2238907791412357E+00,
                -0.2600519549019334E+00,
                -0.3971498098638474E+00,
                -0.1775967713143383E+00,
                0.1506452572509969E+00,
                0.3000792705195556E+00,
                0.1716508071375539E+00,
                -0.9033361118287613E-01,
                -0.2459357644513483E+00,
                -0.1711903004071961E+00,
                0.4768931079683354E-01,
                0.2069261023770678E+00,
                0.1710734761104587E+00,
                -0.1422447282678077E-01
            };

            double[] x_vec =
            {
                -5.0E+00,
                -4.0E+00,
                -3.0E+00,
                -2.0E+00,
                -1.0E+00,
                0.0E+00,
                1.0E+00,
                2.0E+00,
                3.0E+00,
                4.0E+00,
                5.0E+00,
                6.0E+00,
                7.0E+00,
                8.0E+00,
                9.0E+00,
                10.0E+00,
                11.0E+00,
                12.0E+00,
                13.0E+00,
                14.0E+00,
                15.0E+00
            };

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

        public static void bessel_j1_spherical_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_J1_SPHERICAL_VALUES returns some values of the Spherical Bessel function j1.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Sqrt[Pi/(2*x)] * BesselJ[3/2,x]
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
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 21;

            double[] fx_vec =
            {
                0.3330001190255757E-01,
                0.6640038067032223E-01,
                0.1312121544218529E+00,
                0.1928919568034122E+00,
                0.2499855053465475E+00,
                0.3011686789397568E+00,
                0.3452845698577903E+00,
                0.3813753724123076E+00,
                0.4087081401263934E+00,
                0.4267936423844913E+00,
                0.4353977749799916E+00,
                0.4345452193763121E+00,
                0.4245152947656493E+00,
                0.4058301968314685E+00,
                0.3792360591872637E+00,
                0.3456774997623560E+00,
                0.3062665174917607E+00,
                0.2622467779189737E+00,
                0.2149544641595738E+00,
                0.1657769677515280E+00,
                0.1161107492591575E+00
            };

            double[] x_vec =
            {
                0.1E+00,
                0.2E+00,
                0.4E+00,
                0.6E+00,
                0.8E+00,
                1.0E+00,
                1.2E+00,
                1.4E+00,
                1.6E+00,
                1.8E+00,
                2.0E+00,
                2.2E+00,
                2.4E+00,
                2.6E+00,
                2.8E+00,
                3.0E+00,
                3.2E+00,
                3.4E+00,
                3.6E+00,
                3.8E+00,
                4.0E+00
            };

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

        public static void bessel_j1_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_J1_VALUES returns some values of the J1 Bessel function.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      BesselJ[1,x]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 August 2004
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
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 21;

            double[] fx_vec =
            {
                0.3275791375914652E+00,
                0.6604332802354914E-01,
                -0.3390589585259365E+00,
                -0.5767248077568734E+00,
                -0.4400505857449335E+00,
                0.0000000000000000E+00,
                0.4400505857449335E+00,
                0.5767248077568734E+00,
                0.3390589585259365E+00,
                -0.6604332802354914E-01,
                -0.3275791375914652E+00,
                -0.2766838581275656E+00,
                -0.4682823482345833E-02,
                0.2346363468539146E+00,
                0.2453117865733253E+00,
                0.4347274616886144E-01,
                -0.1767852989567215E+00,
                -0.2234471044906276E+00,
                -0.7031805212177837E-01,
                0.1333751546987933E+00,
                0.2051040386135228E+00
            };

            double[] x_vec =
            {
                -5.0E+00,
                -4.0E+00,
                -3.0E+00,
                -2.0E+00,
                -1.0E+00,
                0.0E+00,
                1.0E+00,
                2.0E+00,
                3.0E+00,
                4.0E+00,
                5.0E+00,
                6.0E+00,
                7.0E+00,
                8.0E+00,
                9.0E+00,
                10.0E+00,
                11.0E+00,
                12.0E+00,
                13.0E+00,
                14.0E+00,
                15.0E+00
            };

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

        public static void bessel_jn_values(ref int n_data, ref int nu, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_JN_VALUES returns some values of the Jn Bessel function.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      BesselJ[n,x]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    29 April 2001
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
            //    Output, ref int NU, the order of the function.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 20;

            double[] fx_vec =
            {
                0.1149034849319005E+00,
                0.3528340286156377E+00,
                0.4656511627775222E-01,
                0.2546303136851206E+00,
                -0.5971280079425882E-01,
                0.2497577302112344E-03,
                0.7039629755871685E-02,
                0.2611405461201701E+00,
                -0.2340615281867936E+00,
                -0.8140024769656964E-01,
                0.2630615123687453E-09,
                0.2515386282716737E-06,
                0.1467802647310474E-02,
                0.2074861066333589E+00,
                -0.1138478491494694E+00,
                0.3873503008524658E-24,
                0.3918972805090754E-18,
                0.2770330052128942E-10,
                0.1151336924781340E-04,
                -0.1167043527595797E+00
            };

            int[] nu_vec =
            {
                2, 2, 2, 2,
                2, 5, 5, 5,
                5, 5, 10, 10,
                10, 10, 10, 20,
                20, 20, 20, 20
            };

            double[] x_vec =
            {
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00,
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00,
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00,
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                nu = 0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                nu = nu_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void bessel_jx_values(ref int n_data, ref double nu, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_JX_VALUES returns some values of the Jx Bessel function.
            //
            //  Discussion:
            //
            //    This set of data considers the less common case in which the
            //    index of the Bessel function Jn is actually not an integer.
            //    We may suggest this case by occasionally replacing the symbol
            //    "Jn" by "Jx".
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      BesselJ[n,x]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 April 2007
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
            //    Output, ref double NU, the order of the function.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 28;

            double[] fx_vec =
            {
                0.3544507442114011E+00,
                0.6713967071418031E+00,
                0.5130161365618278E+00,
                0.3020049060623657E+00,
                0.06500818287737578E+00,
                -0.3421679847981618E+00,
                -0.1372637357550505E+00,
                0.1628807638550299E+00,
                0.2402978391234270E+00,
                0.4912937786871623E+00,
                -0.1696513061447408E+00,
                0.1979824927558931E+00,
                -0.1094768729883180E+00,
                0.04949681022847794E+00,
                0.2239245314689158E+00,
                0.2403772011113174E+00,
                0.1966584835818184E+00,
                0.02303721950962553E+00,
                0.3314145508558904E+00,
                0.5461734240402840E+00,
                -0.2616584152094124E+00,
                0.1296035513791289E+00,
                -0.1117432171933552E+00,
                0.03142623570527935E+00,
                0.1717922192746527E+00,
                0.3126634069544786E+00,
                0.1340289119304364E+00,
                0.06235967135106445E+00
            };
            double[] nu_vec =
            {
                0.50E+00,
                0.50E+00,
                0.50E+00,
                0.50E+00,
                0.50E+00,
                0.50E+00,
                0.50E+00,
                0.50E+00,
                1.50E+00,
                1.50E+00,
                1.50E+00,
                1.50E+00,
                1.50E+00,
                2.50E+00,
                2.50E+00,
                2.50E+00,
                2.50E+00,
                2.50E+00,
                1.25E+00,
                1.25E+00,
                1.25E+00,
                1.25E+00,
                1.25E+00,
                2.75E+00,
                2.75E+00,
                2.75E+00,
                2.75E+00,
                2.75E+00
            };

            double[] x_vec =
            {
                0.2E+00,
                1.0E+00,
                2.0E+00,
                2.5E+00,
                3.0E+00,
                5.0E+00,
                10.0E+00,
                20.0E+00,
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00,
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00,
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00,
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                nu = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                nu = nu_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void bessel_k0_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_K0_VALUES returns some values of the K0 Bessel function.
            //
            //  Discussion:
            //
            //    The modified Bessel functions In(Z) and Kn(Z) are solutions of
            //    the differential equation
            //
            //      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
            //
            //    The modified Bessel function K0(Z) corresponds to N = 0.
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      BesselK[0,x]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    15 August 2004
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
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 20;

            double[] fx_vec =
            {
                0.2427069024702017E+01,
                0.1752703855528146E+01,
                0.1114529134524434E+01,
                0.7775220919047293E+00,
                0.5653471052658957E+00,
                0.4210244382407083E+00,
                0.3185082202865936E+00,
                0.2436550611815419E+00,
                0.1879547519693323E+00,
                0.1459314004898280E+00,
                0.1138938727495334E+00,
                0.6234755320036619E-01,
                0.3473950438627925E-01,
                0.1959889717036849E-01,
                0.1115967608585302E-01,
                0.6399857243233975E-02,
                0.3691098334042594E-02,
                0.1243994328013123E-02,
                0.1464707052228154E-03,
                0.1778006231616765E-04
            };

            double[] x_vec =
            {
                0.1E+00,
                0.2E+00,
                0.4E+00,
                0.6E+00,
                0.8E+00,
                1.0E+00,
                1.2E+00,
                1.4E+00,
                1.6E+00,
                1.8E+00,
                2.0E+00,
                2.5E+00,
                3.0E+00,
                3.5E+00,
                4.0E+00,
                4.5E+00,
                5.0E+00,
                6.0E+00,
                8.0E+00,
                10.0E+00
            };

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

        public static void bessel_k0_int_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_K0_INT_VALUES returns some values of the Bessel K0 integral.
            //
            //  Discussion:
            //
            //    The function is defined by:
            //
            //      K0_INT(x) = Integral ( 0 <= t <= x ) K0(t) dt
            //
            //    The data was reported by McLeod.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    29 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Allan McLeod,
            //    Algorithm 757:
            //    MISCFUN: A software package to compute uncommon special functions,
            //    ACM Transactions on Mathematical Software,
            //    Volume 22, Number 3, September 1996, pages 288-301.
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 20;

            double[] fx_vec =
            {
                0.78587929563466784589E-02,
                0.26019991617330578111E-01,
                0.24311842237541167904E+00,
                0.39999633750480508861E+00,
                0.92710252093114907345E+00,
                0.12425098486237782662E+01,
                0.14736757343168286825E+01,
                0.15606495706051741364E+01,
                0.15673873907283660493E+01,
                0.15696345532693743714E+01,
                0.15701153443250786355E+01,
                0.15706574852894436220E+01,
                0.15707793116159788598E+01,
                0.15707942066465767196E+01,
                0.15707962315469192247E+01,
                0.15707963262340149876E+01,
                0.15707963267948756308E+01,
                0.15707963267948966192E+01,
                0.15707963267948966192E+01,
                0.15707963267948966192E+01
            };

            double[] x_vec =
            {
                0.0009765625E+00,
                0.0039062500E+00,
                0.0625000000E+00,
                0.1250000000E+00,
                0.5000000000E+00,
                1.0000000000E+00,
                2.0000000000E+00,
                4.0000000000E+00,
                5.0000000000E+00,
                6.0000000000E+00,
                6.5000000000E+00,
                8.0000000000E+00,
                10.0000000000E+00,
                12.0000000000E+00,
                15.0000000000E+00,
                20.0000000000E+00,
                30.0000000000E+00,
                50.0000000000E+00,
                80.0000000000E+00,
                100.0000000000E+00
            };

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

        public static void bessel_k1_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_K1_VALUES returns some values of the K1 Bessel function.
            //
            //  Discussion:
            //
            //    The modified Bessel functions In(Z) and Kn(Z) are solutions of
            //    the differential equation
            //
            //      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
            //
            //    The modified Bessel function K1(Z) corresponds to N = 1.
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      BesselK[1,x]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    15 August 2004
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
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 20;

            double[] fx_vec =
            {
                0.9853844780870606E+01,
                0.4775972543220472E+01,
                0.2184354424732687E+01,
                0.1302834939763502E+01,
                0.8617816344721803E+00,
                0.6019072301972346E+00,
                0.4345923910607150E+00,
                0.3208359022298758E+00,
                0.2406339113576119E+00,
                0.1826230998017470E+00,
                0.1398658818165224E+00,
                0.7389081634774706E-01,
                0.4015643112819418E-01,
                0.2223939292592383E-01,
                0.1248349888726843E-01,
                0.7078094908968090E-02,
                0.4044613445452164E-02,
                0.1343919717735509E-02,
                0.1553692118050011E-03,
                0.1864877345382558E-04
            };

            double[] x_vec =
            {
                0.1E+00,
                0.2E+00,
                0.4E+00,
                0.6E+00,
                0.8E+00,
                1.0E+00,
                1.2E+00,
                1.4E+00,
                1.6E+00,
                1.8E+00,
                2.0E+00,
                2.5E+00,
                3.0E+00,
                3.5E+00,
                4.0E+00,
                4.5E+00,
                5.0E+00,
                6.0E+00,
                8.0E+00,
                10.0E+00
            };

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

        public static void bessel_kn_values(ref int n_data, ref int nu, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_KN_VALUES returns some values of the Kn Bessel function.
            //
            //  Discussion:
            //
            //    The modified Bessel functions In(Z) and Kn(Z) are solutions of
            //    the differential equation
            //
            //      Z^2 * W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      BesselK[n,x]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 August 2004
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
            //    Output, ref int NU, the order of the function.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 28;

            double[] fx_vec =
            {
                0.4951242928773287E+02,
                0.1624838898635177E+01,
                0.2537597545660559E+00,
                0.1214602062785638E+00,
                0.6151045847174204E-01,
                0.5308943712223460E-02,
                0.2150981700693277E-04,
                0.6329543612292228E-09,
                0.7101262824737945E+01,
                0.6473853909486342E+00,
                0.8291768415230932E-02,
                0.2725270025659869E-04,
                0.3727936773826211E-22,
                0.3609605896012407E+03,
                0.9431049100596467E+01,
                0.3270627371203186E-01,
                0.5754184998531228E-04,
                0.4367182254100986E-22,
                0.1807132899010295E+09,
                0.1624824039795591E+06,
                0.9758562829177810E+01,
                0.1614255300390670E-02,
                0.9150988209987996E-22,
                0.6294369360424535E+23,
                0.5770856852700241E+17,
                0.4827000520621485E+09,
                0.1787442782077055E+03,
                0.1706148379722035E-20
            };

            int[] nu_vec =
            {
                2, 2, 2, 2,
                2, 2, 2, 2,
                3, 3, 3, 3,
                3, 5, 5, 5,
                5, 5, 10, 10,
                10, 10, 10, 20,
                20, 20, 20, 20
            };

            double[] x_vec =
            {
                0.2E+00,
                1.0E+00,
                2.0E+00,
                2.5E+00,
                3.0E+00,
                5.0E+00,
                10.0E+00,
                20.0E+00,
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00,
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00,
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00,
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                nu = 0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                nu = nu_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void bessel_kx_values(ref int n_data, ref double nu, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_KX_VALUES returns some values of the Kx Bessel function.
            //
            //  Discussion:
            //
            //    This set of data considers the less common case in which the
            //    index of the Bessel function Kn is actually not an integer.
            //    We may suggest this case by occasionally replacing the symbol
            //    "Kn" by "Kx".
            //
            //    The modified Bessel functions In(Z) and Kn(Z) are solutions of
            //    the differential equation
            //
            //      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      BesselK[n,x]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 April 2007
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
            //    Output, ref double NU, the order of the function.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 28;

            double[] fx_vec =
            {
                2.294489339798475E+00,
                0.4610685044478946E+00,
                0.1199377719680614E+00,
                0.06506594315400999E+00,
                0.03602598513176459E+00,
                0.003776613374642883E+00,
                0.00001799347809370518E+00,
                5.776373974707445E-10,
                0.9221370088957891E+00,
                0.1799066579520922E+00,
                0.004531936049571459E+00,
                0.00001979282590307570E+00,
                3.486992497366216E-23,
                3.227479531135262E+00,
                0.3897977588961997E+00,
                0.006495775004385758E+00,
                0.00002393132586462789E+00,
                3.627839645299048E-23,
                0.7311451879202114E+00,
                0.1567475478393932E+00,
                0.004257389528177461E+00,
                0.00001915541065869563E+00,
                3.463337593569306E-23,
                4.731184839919541E+00,
                0.4976876225514758E+00,
                0.007300864610941163E+00,
                0.00002546421294106458E+00,
                3.675275677913656E-23
            };

            double[] nu_vec =
            {
                0.50E+00,
                0.50E+00,
                0.50E+00,
                0.50E+00,
                0.50E+00,
                0.50E+00,
                0.50E+00,
                0.50E+00,
                1.50E+00,
                1.50E+00,
                1.50E+00,
                1.50E+00,
                1.50E+00,
                2.50E+00,
                2.50E+00,
                2.50E+00,
                2.50E+00,
                2.50E+00,
                1.25E+00,
                1.25E+00,
                1.25E+00,
                1.25E+00,
                1.25E+00,
                2.75E+00,
                2.75E+00,
                2.75E+00,
                2.75E+00,
                2.75E+00
            };

            double[] x_vec =
            {
                0.2E+00,
                1.0E+00,
                2.0E+00,
                2.5E+00,
                3.0E+00,
                5.0E+00,
                10.0E+00,
                20.0E+00,
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00,
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00,
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00,
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                nu = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                nu = nu_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void bessel_y0_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_Y0_VALUES returns some values of the Y0 Bessel function.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      BesselY[0,x]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 August 2004
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
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 16;

            double[] fx_vec =
            {
                -0.1534238651350367E+01,
                0.8825696421567696E-01,
                0.5103756726497451E+00,
                0.3768500100127904E+00,
                -0.1694073932506499E-01,
                -0.3085176252490338E+00,
                -0.2881946839815792E+00,
                -0.2594974396720926E-01,
                0.2235214893875662E+00,
                0.2499366982850247E+00,
                0.5567116728359939E-01,
                -0.1688473238920795E+00,
                -0.2252373126343614E+00,
                -0.7820786452787591E-01,
                0.1271925685821837E+00,
                0.2054642960389183E+00
            };


            double[] x_vec =
            {
                0.1E+00,
                1.0E+00,
                2.0E+00,
                3.0E+00,
                4.0E+00,
                5.0E+00,
                6.0E+00,
                7.0E+00,
                8.0E+00,
                9.0E+00,
                10.0E+00,
                11.0E+00,
                12.0E+00,
                13.0E+00,
                14.0E+00,
                15.0E+00
            };

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

        public static void bessel_y0_int_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_Y0_INT_VALUES returns some values of the Bessel Y0 integral.
            //
            //  Discussion:
            //
            //    The function is defined by:
            //
            //      Y0_INT(x) = Integral ( 0 <= t <= x ) Y0(t) dt
            //
            //    The data was reported by McLeod.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Allan McLeod,
            //    Algorithm 757:
            //    MISCFUN: A software package to compute uncommon special functions,
            //    ACM Transactions on Mathematical Software,
            //    Volume 22, Number 3, September 1996, pages 288-301.
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 20;

            double[] fx_vec =
            {
                -0.91442642860172110926E-02,
                -0.29682047390397591290E-01,
                -0.25391431276585388961E+00,
                -0.56179545591464028187E+00,
                -0.63706937660742309754E+00,
                -0.28219285008510084123E+00,
                0.38366964785312561103E+00,
                -0.12595061285798929390E+00,
                0.24129031832266684828E+00,
                0.17138069757627037938E+00,
                0.18958142627134083732E+00,
                0.17203846136449706946E+00,
                -0.16821597677215029611E+00,
                -0.93607927351428988679E-01,
                0.88229711948036648408E-01,
                -0.89324662736274161841E-02,
                -0.54814071000063488284E-01,
                -0.94958246003466381588E-01,
                -0.19598064853404969850E-01,
                -0.83084772357154773468E-02
            };

            double[] x_vec =
            {
                0.0019531250E+00,
                0.0078125000E+00,
                0.1250000000E+00,
                0.5000000000E+00,
                1.0000000000E+00,
                2.0000000000E+00,
                4.0000000000E+00,
                6.0000000000E+00,
                10.0000000000E+00,
                16.0000000000E+00,
                16.2500000000E+00,
                17.0000000000E+00,
                20.0000000000E+00,
                25.0000000000E+00,
                30.0000000000E+00,
                40.0000000000E+00,
                50.0000000000E+00,
                70.0000000000E+00,
                100.0000000000E+00,
                125.0000000000E+00
            };

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

        public static void bessel_y0_spherical_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_Y0_SPHERICAL_VALUES returns some values of the Spherical Bessel function y0.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Sqrt[Pi/(2*x)] * BesselY[1/2,x]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 August 2004
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
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 21;

            double[] fx_vec =
            {
                -0.9950041652780258E+01,
                -0.4900332889206208E+01,
                -0.2302652485007213E+01,
                -0.1375559358182797E+01,
                -0.8708833866839568E+00,
                -0.5403023058681397E+00,
                -0.3019647953972280E+00,
                -0.1214051020716007E+00,
                0.1824970143830545E-01,
                0.1262233859406039E+00,
                0.2080734182735712E+00,
                0.2675005078433390E+00,
                0.3072473814755190E+00,
                0.3295725974495951E+00,
                0.3365079788102351E+00,
                0.3299974988668152E+00,
                0.3119671174358603E+00,
                0.2843524095821944E+00,
                0.2490995600928186E+00,
                0.2081493978722149E+00,
                0.1634109052159030E+00
            };

            double[] x_vec =
            {
                0.1E+00,
                0.2E+00,
                0.4E+00,
                0.6E+00,
                0.8E+00,
                1.0E+00,
                1.2E+00,
                1.4E+00,
                1.6E+00,
                1.8E+00,
                2.0E+00,
                2.2E+00,
                2.4E+00,
                2.6E+00,
                2.8E+00,
                3.0E+00,
                3.2E+00,
                3.4E+00,
                3.6E+00,
                3.8E+00,
                4.0E+00
            };

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

        public static void bessel_y1_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_Y1_VALUES returns some values of the Y1 Bessel function.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      BesselY[1,x]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 August 2004
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
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 16;

            double[] fx_vec =
            {
                -0.6458951094702027E+01,
                -0.7812128213002887E+00,
                -0.1070324315409375E+00,
                0.3246744247918000E+00,
                0.3979257105571000E+00,
                0.1478631433912268E+00,
                -0.1750103443003983E+00,
                -0.3026672370241849E+00,
                -0.1580604617312475E+00,
                0.1043145751967159E+00,
                0.2490154242069539E+00,
                0.1637055374149429E+00,
                -0.5709921826089652E-01,
                -0.2100814084206935E+00,
                -0.1666448418561723E+00,
                0.2107362803687351E-01
            };

            double[] x_vec =
            {
                0.1E+00,
                1.0E+00,
                2.0E+00,
                3.0E+00,
                4.0E+00,
                5.0E+00,
                6.0E+00,
                7.0E+00,
                8.0E+00,
                9.0E+00,
                10.0E+00,
                11.0E+00,
                12.0E+00,
                13.0E+00,
                14.0E+00,
                15.0E+00
            };

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

        public static void bessel_y1_spherical_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_Y1_SPHERICAL_VALUES returns some values of the Spherical Bessel function y1.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Sqrt[Pi/(2*x)] * BesselY[3/2,x]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 August 2004
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
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 21;

            double[] fx_vec =
            {
                -0.1004987506942709E+03,
                -0.2549501110000635E+02,
                -0.6730177068289658E+01,
                -0.3233669719296388E+01,
                -0.1985299346979349E+01,
                -0.1381773290676036E+01,
                -0.1028336567803712E+01,
                -0.7906105943286149E+00,
                -0.6133274385019998E+00,
                -0.4709023582986618E+00,
                -0.3506120042760553E+00,
                -0.2459072254437506E+00,
                -0.1534232496148467E+00,
                -0.7151106706610352E-01,
                0.5427959479750482E-03,
                0.6295916360231598E-01,
                0.1157316440198251E+00,
                0.1587922092967723E+00,
                0.1921166676076864E+00,
                0.2157913917934037E+00,
                0.2300533501309578E+00
            };

            double[] x_vec =
            {
                0.1E+00,
                0.2E+00,
                0.4E+00,
                0.6E+00,
                0.8E+00,
                1.0E+00,
                1.2E+00,
                1.4E+00,
                1.6E+00,
                1.8E+00,
                2.0E+00,
                2.2E+00,
                2.4E+00,
                2.6E+00,
                2.8E+00,
                3.0E+00,
                3.2E+00,
                3.4E+00,
                3.6E+00,
                3.8E+00,
                4.0E+00
            };

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

        public static void bessel_yn_values(ref int n_data, ref int nu, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_YN_VALUES returns some values of the Yn Bessel function.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      BesselY[n,x]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 August 2004
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
            //    Output, ref int NU, the order of the function.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            //
            int N_MAX = 20;

            double[] fx_vec =
            {
                -0.1650682606816254E+01,
                -0.6174081041906827E+00,
                0.3676628826055245E+00,
                -0.5868082442208615E-02,
                0.9579316872759649E-01,
                -0.2604058666258122E+03,
                -0.9935989128481975E+01,
                -0.4536948224911019E+00,
                0.1354030476893623E+00,
                -0.7854841391308165E-01,
                -0.1216180142786892E+09,
                -0.1291845422080393E+06,
                -0.2512911009561010E+02,
                -0.3598141521834027E+00,
                0.5723897182053514E-02,
                -0.4113970314835505E+23,
                -0.4081651388998367E+17,
                -0.5933965296914321E+09,
                -0.1597483848269626E+04,
                0.1644263394811578E-01
            };

            int[] nu_vec =
            {
                2, 2, 2, 2,
                2, 5, 5, 5,
                5, 5, 10, 10,
                10, 10, 10, 20,
                20, 20, 20, 20
            };

            double[] x_vec =
            {
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00,
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00,
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00,
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                nu = 0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                nu = nu_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void bessel_yx_values(ref int n_data, ref double nu, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_YX_VALUES returns some values of the Yx Bessel function.
            //
            //  Discussion:
            //
            //    This set of data considers the less common case in which the
            //    index of the Bessel function Yn is actually not an integer.
            //    We may suggest this case by occasionally replacing the symbol
            //    "Yn" by "Yx".
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      BesselY[n,x]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 April 2007
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
            //    Output, ref double NU, the order of the function.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 28;

            double[] fx_vec =
            {
                -1.748560416961876E+00,
                -0.4310988680183761E+00,
                0.2347857104062485E+00,
                0.4042783022390569E+00,
                0.4560488207946332E+00,
                -0.1012177091851084E+00,
                0.2117088663313982E+00,
                -0.07280690478506185E+00,
                -1.102495575160179E+00,
                -0.3956232813587035E+00,
                0.3219244429611401E+00,
                0.1584346223881903E+00,
                0.02742813676191382E+00,
                -2.876387857462161E+00,
                -0.8282206324443037E+00,
                0.2943723749617925E+00,
                -0.1641784796149411E+00,
                0.1105304445562544E+00,
                -0.9319659251969881E+00,
                -0.2609445010948933E+00,
                0.2492796362185881E+00,
                0.2174410301416733E+00,
                -0.01578576650557229E+00,
                -4.023453301501028E+00,
                -0.9588998694752389E+00,
                0.2264260361047367E+00,
                -0.2193617736566760E+00,
                0.09413988344515077E+00
            };

            double[] nu_vec =
            {
                0.50E+00,
                0.50E+00,
                0.50E+00,
                0.50E+00,
                0.50E+00,
                0.50E+00,
                0.50E+00,
                0.50E+00,
                1.50E+00,
                1.50E+00,
                1.50E+00,
                1.50E+00,
                1.50E+00,
                2.50E+00,
                2.50E+00,
                2.50E+00,
                2.50E+00,
                2.50E+00,
                1.25E+00,
                1.25E+00,
                1.25E+00,
                1.25E+00,
                1.25E+00,
                2.75E+00,
                2.75E+00,
                2.75E+00,
                2.75E+00,
                2.75E+00
            };

            double[] x_vec =
            {
                0.2E+00,
                1.0E+00,
                2.0E+00,
                2.5E+00,
                3.0E+00,
                5.0E+00,
                10.0E+00,
                20.0E+00,
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00,
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00,
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00,
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                nu = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                nu = nu_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

    }
}