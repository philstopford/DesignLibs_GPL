namespace Burkardt.Values
{
    public static class Struve
    {
        public static void struve_h0_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    STRUVE_H0_VALUES returns some values of the Struve H0 function.
            //
            //  Discussion:
            //
            //    The function is defined by:
            //
            //      HO(x) = 2/pi * integral ( 0 <= t <= pi/2 ) sin ( x * cos ( t ) ) dt
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      StruveH[0,x]
            //
            //    The data was reported by McLeod.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 September 2004
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
            //    Allan McLeod,
            //    Algorithm 757:
            //    MISCFUN: A software package to compute uncommon special functions,
            //    ACM Transactions on Mathematical Software,
            //    Volume 22, Number 3, September 1996, pages 288-301.
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
                0.12433974658847434366E-02,
                -0.49735582423748415045E-02,
                0.39771469054536941564E-01,
                -0.15805246001653314198E+00,
                0.56865662704828795099E+00,
                0.66598399314899916605E+00,
                0.79085884950809589255E+00,
                -0.13501457342248639716E+00,
                0.20086479668164503137E+00,
                -0.11142097800261991552E+00,
                -0.17026804865989885869E+00,
                -0.13544931808186467594E+00,
                0.94393698081323450897E-01,
                -0.10182482016001510271E+00,
                0.96098421554162110012E-01,
                -0.85337674826118998952E-01,
                -0.76882290637052720045E-01,
                0.47663833591418256339E-01,
                -0.70878751689647343204E-01,
                0.65752908073352785368E-01
            };

            double[] x_vec =
            {
                0.0019531250E+00,
                -0.0078125000E+00,
                0.0625000000E+00,
                -0.2500000000E+00,
                1.0000000000E+00,
                1.2500000000E+00,
                2.0000000000E+00,
                -4.0000000000E+00,
                7.5000000000E+00,
                11.0000000000E+00,
                11.5000000000E+00,
                -16.0000000000E+00,
                20.0000000000E+00,
                25.0000000000E+00,
                -30.0000000000E+00,
                50.0000000000E+00,
                75.0000000000E+00,
                -80.0000000000E+00,
                100.0000000000E+00,
                -125.0000000000E+00
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

        public static void struve_h1_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    STRUVE_H1_VALUES returns some values of the Struve H1 function.
            //
            //  Discussion:
            //
            //    The function is defined by:
            //
            //      H1(x) = 2*x/pi * integral ( 0 <= t <= pi/2 )
            //        sin ( x * cos ( t ) )^2 * sin ( t ) dt
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      StruveH[1,x]
            //
            //    The data was reported by McLeod.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 September 2004
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
            //    Allan McLeod,
            //    Algorithm 757:
            //    MISCFUN: A software package to compute uncommon special functions,
            //    ACM Transactions on Mathematical Software,
            //    Volume 22, Number 3, September 1996, pages 288-301.
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
                0.80950369576367526071E-06,
                0.12952009724113229165E-04,
                0.82871615165407083021E-03,
                0.13207748375849572564E-01,
                0.19845733620194439894E+00,
                0.29853823231804706294E+00,
                0.64676372828356211712E+00,
                0.10697266613089193593E+01,
                0.38831308000420560970E+00,
                0.74854243745107710333E+00,
                0.84664854642567359993E+00,
                0.58385732464244384564E+00,
                0.80600584524215772824E+00,
                0.53880362132692947616E+00,
                0.72175037834698998506E+00,
                0.58007844794544189900E+00,
                0.60151910385440804463E+00,
                0.70611511147286827018E+00,
                0.61631110327201338454E+00,
                0.62778480765443656489E+00
            };

            double[] x_vec =
            {
                0.0019531250E+00,
                -0.0078125000E+00,
                0.0625000000E+00,
                -0.2500000000E+00,
                1.0000000000E+00,
                1.2500000000E+00,
                2.0000000000E+00,
                -4.0000000000E+00,
                7.5000000000E+00,
                11.0000000000E+00,
                11.5000000000E+00,
                -16.0000000000E+00,
                20.0000000000E+00,
                25.0000000000E+00,
                -30.0000000000E+00,
                50.0000000000E+00,
                75.0000000000E+00,
                -80.0000000000E+00,
                100.0000000000E+00,
                -125.0000000000E+00
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

        public static void struve_l0_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    STRUVE_L0_VALUES returns some values of the Struve L0 function.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      StruveL[0,x]
            //
            //    The data was reported by McLeod.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 September 2004
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
            //    Allan McLeod,
            //    Algorithm 757:
            //    MISCFUN: A software package to compute uncommon special functions,
            //    ACM Transactions on Mathematical Software,
            //    Volume 22, Number 3, September 1996, pages 288-301.
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
                0.12433985199262820188E-02,
                -0.19896526647882937004E-01,
                0.79715713253115014945E-01,
                -0.32724069939418078025E+00,
                0.71024318593789088874E+00,
                0.19374337579914456612E+01,
                -0.11131050203248583431E+02,
                0.16850062034703267148E+03,
                -0.28156522493745948555E+04,
                0.89344618796978400815E+06,
                0.11382025002851451057E+07,
                -0.23549701855860190304E+07,
                0.43558282527641046718E+08,
                0.49993516476037957165E+09,
                -0.57745606064408041689E+10,
                0.78167229782395624524E+12,
                -0.14894774793419899908E+17,
                0.29325537838493363267E+21,
                0.58940770556098011683E+25,
                -0.12015889579125463605E+30
            };

            double[] x_vec =
            {
                0.0019531250E+00,
                -0.0312500000E+00,
                0.1250000000E+00,
                -0.5000000000E+00,
                1.0000000000E+00,
                2.0000000000E+00,
                -4.0000000000E+00,
                7.0000000000E+00,
                -10.0000000000E+00,
                16.0000000000E+00,
                16.2500000000E+00,
                -17.0000000000E+00,
                20.0000000000E+00,
                22.5000000000E+00,
                -25.0000000000E+00,
                30.0000000000E+00,
                -40.0000000000E+00,
                50.0000000000E+00,
                60.0000000000E+00,
                -70.0000000000E+00
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

        public static void struve_l1_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    STRUVE_L1_VALUES returns some values of the Struve L1 function.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      StruveL[1,x]
            //
            //    The data was reported by McLeod.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 September 2004
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
            //    Allan McLeod,
            //    Algorithm 757:
            //    MISCFUN: A software package to compute uncommon special functions,
            //    ACM Transactions on Mathematical Software,
            //    Volume 22, Number 3, September 1996, pages 288-301.
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
                0.80950410749865126939E-06,
                0.20724649092571514607E-03,
                0.33191834066894516744E-02,
                0.53942182623522663292E-01,
                0.22676438105580863683E+00,
                0.11027597873677158176E+01,
                0.91692778117386847344E+01,
                0.15541656652426660966E+03,
                0.26703582852084829694E+04,
                0.86505880175304633906E+06,
                0.11026046613094942620E+07,
                0.22846209494153934787E+07,
                0.42454972750111979449E+08,
                0.48869614587997695539E+09,
                0.56578651292431051863E+10,
                0.76853203893832108948E+12,
                0.14707396163259352103E+17,
                0.29030785901035567967E+21,
                0.58447515883904682813E+25,
                0.11929750788892311875E+30
            };

            double[] x_vec =
            {
                0.0019531250E+00,
                -0.0312500000E+00,
                0.1250000000E+00,
                -0.5000000000E+00,
                1.0000000000E+00,
                2.0000000000E+00,
                -4.0000000000E+00,
                7.0000000000E+00,
                -10.0000000000E+00,
                16.0000000000E+00,
                16.2500000000E+00,
                -17.0000000000E+00,
                20.0000000000E+00,
                22.5000000000E+00,
                -25.0000000000E+00,
                30.0000000000E+00,
                -40.0000000000E+00,
                50.0000000000E+00,
                60.0000000000E+00,
                -70.0000000000E+00
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

    }
}