namespace Burkardt.Values
{
    public static class Transportation
    {
        public static void tran02_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRAN02_VALUES returns some values of the order 2 transportation function.
            //
            //  Discussion:
            //
            //    The function is defined by:
            //
            //      TRAN02(x) = integral ( 0 <= t <= x ) t^2 exp(t) / ( exp(t) - 1 )^2 dt
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
                0.19531247930394515480E-02,
                0.31249152314331109004E-01,
                0.12494577194783451032E+00,
                0.49655363615640595865E+00,
                0.97303256135517012845E+00,
                0.14121978695932525805E+01,
                0.18017185674405776809E+01,
                0.21350385339277043015E+01,
                0.24110500490169534620E+01,
                0.28066664045631179931E+01,
                0.28777421863296234131E+01,
                0.30391706043438554330E+01,
                0.31125074928667355940E+01,
                0.31656687817738577185E+01,
                0.32623520367816009184E+01,
                0.32843291144979517358E+01,
                0.32897895167775788137E+01,
                0.32898672226665499687E+01,
                0.32898681336064325400E+01,
                0.32898681336964528724E+01
            };

            double[] x_vec =
            {
                0.0019531250E+00,
                0.0312500000E+00,
                0.1250000000E+00,
                0.5000000000E+00,
                1.0000000000E+00,
                1.5000000000E+00,
                2.0000000000E+00,
                2.5000000000E+00,
                3.0000000000E+00,
                4.0000000000E+00,
                4.2500000000E+00,
                5.0000000000E+00,
                5.5000000000E+00,
                6.0000000000E+00,
                8.0000000000E+00,
                10.0000000000E+00,
                15.0000000000E+00,
                20.0000000000E+00,
                30.0000000000E+00,
                50.0000000000E+00
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

        public static void tran03_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRAN03_VALUES returns some values of the order 3 transportation function.
            //
            //  Discussion:
            //
            //    The function is defined by:
            //
            //      TRAN03(x) = integral ( 0 <= t <= x ) t^3 * exp(t) / ( exp(t) - 1 )^2 dt
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
                0.19073483296476379584E-05,
                0.48826138243180786081E-03,
                0.78074163848431205820E-02,
                0.12370868718812031049E+00,
                0.47984100657241749994E+00,
                0.10269431622039754738E+01,
                0.17063547219458658863E+01,
                0.24539217444475937661E+01,
                0.32106046629422467723E+01,
                0.45792174372291563703E+01,
                0.48722022832940370805E+01,
                0.56143866138422732286E+01,
                0.59984455864575470009E+01,
                0.63033953673480961120E+01,
                0.69579908688361166266E+01,
                0.71503227120085929750E+01,
                0.72110731475871876393E+01,
                0.72123221966388461839E+01,
                0.72123414161609465119E+01,
                0.72123414189575656868E+01
            };

            double[] x_vec =
            {
                0.0019531250E+00,
                0.0312500000E+00,
                0.1250000000E+00,
                0.5000000000E+00,
                1.0000000000E+00,
                1.5000000000E+00,
                2.0000000000E+00,
                2.5000000000E+00,
                3.0000000000E+00,
                4.0000000000E+00,
                4.2500000000E+00,
                5.0000000000E+00,
                5.5000000000E+00,
                6.0000000000E+00,
                8.0000000000E+00,
                10.0000000000E+00,
                15.0000000000E+00,
                20.0000000000E+00,
                30.0000000000E+00,
                50.0000000000E+00
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

        public static void tran04_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRAN04_VALUES returns some values of the order 4 transportation function.
            //
            //  Discussion:
            //
            //    The function is defined by:
            //
            //      TRAN04(x) = integral ( 0 <= t <= x ) t^4 * exp(t) / ( exp(t) - 1 )^2 dt
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
                0.24835263919461834041E-08,
                0.10172029353616724881E-04,
                0.65053332405940765479E-03,
                0.41150448004155727767E-01,
                0.31724404523442648241E+00,
                0.10079442901142373591E+01,
                0.22010881024333408363E+01,
                0.38846508619156545210E+01,
                0.59648223973714765245E+01,
                0.10731932392998622219E+02,
                0.11940028876819364777E+02,
                0.15359784316882182982E+02,
                0.17372587633093742893E+02,
                0.19122976016053166969E+02,
                0.23583979156921941515E+02,
                0.25273667677030441733E+02,
                0.25955198214572256372E+02,
                0.25975350935212241910E+02,
                0.25975757522084093747E+02,
                0.25975757609067315288E+02
            };

            double[] x_vec =
            {
                0.0019531250E+00,
                0.0312500000E+00,
                0.1250000000E+00,
                0.5000000000E+00,
                1.0000000000E+00,
                1.5000000000E+00,
                2.0000000000E+00,
                2.5000000000E+00,
                3.0000000000E+00,
                4.0000000000E+00,
                4.2500000000E+00,
                5.0000000000E+00,
                5.5000000000E+00,
                6.0000000000E+00,
                8.0000000000E+00,
                10.0000000000E+00,
                15.0000000000E+00,
                20.0000000000E+00,
                30.0000000000E+00,
                50.0000000000E+00
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

        public static void tran05_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRAN05_VALUES returns some values of the order 5 transportation function.
            //
            //  Discussion:
            //
            //    The function is defined by:
            //
            //      TRAN05(x) = integral ( 0 <= t <= x ) t^5 * exp(t) / ( exp(t) - 1 )^2 dt
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
                0.36379780361036116971E-11,
                0.23840564453948442379E-06,
                0.60982205372226969189E-04,
                0.15410004586376649337E-01,
                0.23661587923909478926E+00,
                0.11198756851307629651E+01,
                0.32292901663684049171E+01,
                0.70362973105160654056E+01,
                0.12770557691044159511E+02,
                0.29488339015245845447E+02,
                0.34471340540362254586E+02,
                0.50263092218175187785E+02,
                0.60819909101127165207E+02,
                0.70873334429213460498E+02,
                0.10147781242977788097E+03,
                0.11638074540242071077E+03,
                0.12409623901262967878E+03,
                0.12442270155632550228E+03,
                0.12443132790838589548E+03,
                0.12443133061720432435E+03
            };

            double[] x_vec =
            {
                0.0019531250E+00,
                0.0312500000E+00,
                0.1250000000E+00,
                0.5000000000E+00,
                1.0000000000E+00,
                1.5000000000E+00,
                2.0000000000E+00,
                2.5000000000E+00,
                3.0000000000E+00,
                4.0000000000E+00,
                4.2500000000E+00,
                5.0000000000E+00,
                5.5000000000E+00,
                6.0000000000E+00,
                8.0000000000E+00,
                10.0000000000E+00,
                15.0000000000E+00,
                20.0000000000E+00,
                30.0000000000E+00,
                50.0000000000E+00
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

        public static void tran06_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRAN06_VALUES returns some values of the order 6 transportation function.
            //
            //  Discussion:
            //
            //    The function is defined by:
            //
            //      TRAN06(x) = integral ( 0 <= t <= x ) t^6 * exp(t) / ( exp(t) - 1 )^2 dt
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
                0.56843405953641209574E-14,
                0.59601180165247401484E-08,
                0.60978424397580572815E-05,
                0.61578909866319494394E-02,
                0.18854360275680840514E+00,
                0.13319251347921659134E+01,
                0.50857202271697616755E+01,
                0.13729222365466557122E+02,
                0.29579592481641441292E+02,
                0.88600835706899853768E+02,
                0.10916037113373004909E+03,
                0.18224323749575359518E+03,
                0.23765383125586756031E+03,
                0.29543246745959381136E+03,
                0.50681244381280455592E+03,
                0.63878231134946125623E+03,
                0.72699203556994876111E+03,
                0.73230331643146851717E+03,
                0.73248692015882096369E+03,
                0.73248700462879996604E+03
            };

            double[] x_vec =
            {
                0.0019531250E+00,
                0.0312500000E+00,
                0.1250000000E+00,
                0.5000000000E+00,
                1.0000000000E+00,
                1.5000000000E+00,
                2.0000000000E+00,
                2.5000000000E+00,
                3.0000000000E+00,
                4.0000000000E+00,
                4.2500000000E+00,
                5.0000000000E+00,
                5.5000000000E+00,
                6.0000000000E+00,
                8.0000000000E+00,
                10.0000000000E+00,
                15.0000000000E+00,
                20.0000000000E+00,
                30.0000000000E+00,
                50.0000000000E+00
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

        public static void tran07_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRAN07_VALUES returns some values of the order 7 transportation function.
            //
            //  Discussion:
            //
            //    The function is defined by:
            //
            //      TRAN07(x) = integral ( 0 <= t <= x ) t^7 * exp(t) / ( exp(t) - 1 )^2 dt
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
                0.92518563327283409427E-17,
                0.15521095556949867541E-09,
                0.63516238373841716290E-06,
                0.25638801246626135714E-02,
                0.15665328993811649746E+00,
                0.16538225039181097423E+01,
                0.83763085709508211054E+01,
                0.28078570717830763747E+02,
                0.72009676046751991365E+02,
                0.28174905701691911450E+03,
                0.36660227975327792529E+03,
                0.70556067982603601123E+03,
                0.99661927562755629434E+03,
                0.13288914430417403901E+04,
                0.27987640273169129925E+04,
                0.39721376409416504325E+04,
                0.49913492839319899726E+04,
                0.50781562639825019000E+04,
                0.50820777202028708434E+04,
                0.50820803580047164618E+04
            };

            double[] x_vec =
            {
                0.0019531250E+00,
                0.0312500000E+00,
                0.1250000000E+00,
                0.5000000000E+00,
                1.0000000000E+00,
                1.5000000000E+00,
                2.0000000000E+00,
                2.5000000000E+00,
                3.0000000000E+00,
                4.0000000000E+00,
                4.2500000000E+00,
                5.0000000000E+00,
                5.5000000000E+00,
                6.0000000000E+00,
                8.0000000000E+00,
                10.0000000000E+00,
                15.0000000000E+00,
                20.0000000000E+00,
                30.0000000000E+00,
                50.0000000000E+00
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

        public static void tran08_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRAN08_VALUES returns some values of the order 8 transportation function.
            //
            //  Discussion:
            //
            //    The function is defined by:
            //
            //      TRAN08(x) = integral ( 0 <= t <= x ) t^8 * exp(t) / ( exp(t) - 1 )^2 dt
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
                0.15488598634539359463E-19,
                0.41574269117845953797E-11,
                0.68050651245227411689E-07,
                0.10981703519563009836E-02,
                0.13396432776187883834E+00,
                0.21153387806998617182E+01,
                0.14227877028750735641E+02,
                0.59312061431647843226E+02,
                0.18139614577043147745E+03,
                0.93148001928992220863E+03,
                0.12817928112604611804E+04,
                0.28572838386329242218E+04,
                0.43872971687877730010E+04,
                0.62993229139406657611E+04,
                0.16589426277154888511E+05,
                0.27064780798797398935E+05,
                0.38974556062543661284E+05,
                0.40400240716905025786E+05,
                0.40484316504120655568E+05,
                0.40484399001892184901E+05
            };

            double[] x_vec =
            {
                0.0019531250E+00,
                0.0312500000E+00,
                0.1250000000E+00,
                0.5000000000E+00,
                1.0000000000E+00,
                1.5000000000E+00,
                2.0000000000E+00,
                2.5000000000E+00,
                3.0000000000E+00,
                4.0000000000E+00,
                4.2500000000E+00,
                5.0000000000E+00,
                5.5000000000E+00,
                6.0000000000E+00,
                8.0000000000E+00,
                10.0000000000E+00,
                15.0000000000E+00,
                20.0000000000E+00,
                30.0000000000E+00,
                50.0000000000E+00
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

        public static void tran09_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRAN09_VALUES returns some values of the order 9 transportation function.
            //
            //  Discussion:
            //
            //    The function is defined by:
            //
            //      TRAN09(x) = integral ( 0 <= t <= x ) t^9 * exp(t) / ( exp(t) - 1 )^2 dt
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
                0.26469772870084897671E-22,
                0.11367943653594246210E-12,
                0.74428246255329800255E-08,
                0.48022728485415366194E-03,
                0.11700243014358676725E+00,
                0.27648973910899914391E+01,
                0.24716631405829192997E+02,
                0.12827119828849828583E+03,
                0.46842894800662208986E+03,
                0.31673967371627895718E+04,
                0.46140886546630195390E+04,
                0.11952718545392302185E+05,
                0.20001612666477027728E+05,
                0.31011073271851366554E+05,
                0.10352949905541130133E+06,
                0.19743173017140591390E+06,
                0.33826030414658460679E+06,
                0.36179607036750755227E+06,
                0.36360622124777561525E+06,
                0.36360880558827162725E+06
            };

            double[] x_vec =
            {
                0.0019531250E+00,
                0.0312500000E+00,
                0.1250000000E+00,
                0.5000000000E+00,
                1.0000000000E+00,
                1.5000000000E+00,
                2.0000000000E+00,
                2.5000000000E+00,
                3.0000000000E+00,
                4.0000000000E+00,
                4.2500000000E+00,
                5.0000000000E+00,
                5.5000000000E+00,
                6.0000000000E+00,
                8.0000000000E+00,
                10.0000000000E+00,
                15.0000000000E+00,
                20.0000000000E+00,
                30.0000000000E+00,
                50.0000000000E+00
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