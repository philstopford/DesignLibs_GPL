namespace Burkardt.Values;

public static class Log
{

    public static void log_values(ref int n_data, ref double x, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_VALUES returns some values of the natural logarithm function.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Log[x]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 June 2007
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
        const int N_MAX = 20;

        double[] fx_vec =
        {
            -11.512925464970228420E+00,
            -4.6051701859880913680E+00,
            -2.3025850929940456840E+00,
            -1.6094379124341003746E+00,
            -1.2039728043259359926E+00,
            -0.91629073187415506518E+00,
            -0.69314718055994530942E+00,
            -0.51082562376599068321E+00,
            -0.35667494393873237891E+00,
            -0.22314355131420975577E+00,
            -0.10536051565782630123E+00,
            0.00000000000000000000E+00,
            0.69314718055994530942E+00,
            1.0986122886681096914E+00,
            1.1447298858494001741E+00,
            1.6094379124341003746E+00,
            2.3025850929940456840E+00,
            2.9957322735539909934E+00,
            4.6051701859880913680E+00,
            18.631401766168018033E+00
        };

        double[] x_vec =
        {
            1.0E-05,
            1.0E-02,
            0.1E+00,
            0.2E+00,
            0.3E+00,
            0.4E+00,
            0.5E+00,
            0.6E+00,
            0.7E+00,
            0.8E+00,
            0.9E+00,
            1.0E+00,
            2.0E+00,
            3.0E+00,
            3.1415926535897932385E+00,
            5.0E+00,
            10.0E+00,
            20.0E+00,
            100.0E+00,
            123456789.0E+00
        };

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        n_data += 1;

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

    public static void log_normal_cdf_values(ref int n_data, ref double mu, ref double sigma,
            ref double x, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_NORMAL_CDF_VALUES returns some values of the Log Normal CDF.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Needs["Statistics`ContinuousDistributions`"]
        //      dist = LogNormalDistribution [ mu, sigma ]
        //      CDF [ dist, x ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 August 2004
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
        //    Output, ref double MU, the mean of the distribution.
        //
        //    Output, ref double SIGMA, the standard deviation of the distribution.
        //
        //    Output, ref double X, the argument of the function.
        //
        //    Output, ref double FX, the value of the function.
        //
    {
        const int N_MAX = 12;

        double[] fx_vec =
        {
            0.2275013194817921E-01,
            0.2697049307349095E+00,
            0.5781741008028732E+00,
            0.7801170895122241E+00,
            0.4390310097476894E+00,
            0.4592655190218048E+00,
            0.4694258497695908E+00,
            0.4755320473858733E+00,
            0.3261051056816658E+00,
            0.1708799040927608E+00,
            0.7343256357952060E-01,
            0.2554673736161761E-01
        };

        double[] mu_vec =
        {
            0.1000000000000000E+01,
            0.1000000000000000E+01,
            0.1000000000000000E+01,
            0.1000000000000000E+01,
            0.1000000000000000E+01,
            0.1000000000000000E+01,
            0.1000000000000000E+01,
            0.1000000000000000E+01,
            0.2000000000000000E+01,
            0.3000000000000000E+01,
            0.4000000000000000E+01,
            0.5000000000000000E+01
        };

        double[] sigma_vec =
        {
            0.5000000000000000E+00,
            0.5000000000000000E+00,
            0.5000000000000000E+00,
            0.5000000000000000E+00,
            0.2000000000000000E+01,
            0.3000000000000000E+01,
            0.4000000000000000E+01,
            0.5000000000000000E+01,
            0.2000000000000000E+01,
            0.2000000000000000E+01,
            0.2000000000000000E+01,
            0.2000000000000000E+01
        };

        double[] x_vec =
        {
            0.1000000000000000E+01,
            0.2000000000000000E+01,
            0.3000000000000000E+01,
            0.4000000000000000E+01,
            0.2000000000000000E+01,
            0.2000000000000000E+01,
            0.2000000000000000E+01,
            0.2000000000000000E+01,
            0.3000000000000000E+01,
            0.3000000000000000E+01,
            0.3000000000000000E+01,
            0.3000000000000000E+01
        };

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        n_data += 1;

        if (N_MAX < n_data)
        {
            n_data = 0;
            mu = 0.0;
            sigma = 0.0;
            x = 0.0;
            fx = 0.0;
        }
        else
        {
            mu = mu_vec[n_data - 1];
            sigma = sigma_vec[n_data - 1];
            x = x_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }
    }

    public static void log_series_cdf_values(ref int n_data, ref double t, ref int n, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_SERIES_CDF_VALUES returns some values of the log series CDF.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Needs["Statistics`DiscreteDistributions`]
        //      dist = LogSeriesDistribution [ t ]
        //      CDF [ dist, n ]
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
        //    Output, ref double T, the parameter of the function.
        //
        //    Output, ref int N, the argument of the function.
        //
        //    Output, ref double FX, the value of the function.
        //
    {
        const int N_MAX = 29;

        double[] fx_vec =
        {
            0.9491221581029903E+00,
            0.9433541128559735E+00,
            0.9361094611773272E+00,
            0.9267370278044118E+00,
            0.9141358246245129E+00,
            0.8962840235449100E+00,
            0.8690148741955517E+00,
            0.8221011541254772E+00,
            0.7213475204444817E+00,
            0.6068261510845583E+00,
            0.5410106403333613E+00,
            0.4970679476476894E+00,
            0.4650921887927060E+00,
            0.4404842934597863E+00,
            0.4207860535926143E+00,
            0.4045507673897055E+00,
            0.3908650337129266E+00,
            0.2149757685421097E+00,
            0.0000000000000000E+00,
            0.2149757685421097E+00,
            0.3213887739704539E+00,
            0.3916213575531612E+00,
            0.4437690508633213E+00,
            0.4850700239649681E+00,
            0.5191433267738267E+00,
            0.5480569580144867E+00,
            0.5731033910767085E+00,
            0.5951442521714636E+00,
            0.6147826594068904E+00
        };

        int[] n_vec =
        {
            1, 1, 1, 1, 1,
            1, 1, 1, 1, 1,
            1, 1, 1, 1, 1,
            1, 1, 1, 0, 1,
            2, 3, 4, 5, 6,
            7, 8, 9, 10
        };

        double[] t_vec =
        {
            0.1000000000000000E+00,
            0.1111111111111111E+00,
            0.1250000000000000E+00,
            0.1428571428571429E+00,
            0.1666666666666667E+00,
            0.2000000000000000E+00,
            0.2500000000000000E+00,
            0.3333333333333333E+00,
            0.5000000000000000E+00,
            0.6666666666666667E+00,
            0.7500000000000000E+00,
            0.8000000000000000E+00,
            0.8333333333333333E+00,
            0.8571485714857149E+00,
            0.8750000000000000E+00,
            0.8888888888888889E+00,
            0.9000000000000000E+00,
            0.9900000000000000E+00,
            0.9900000000000000E+00,
            0.9900000000000000E+00,
            0.9900000000000000E+00,
            0.9900000000000000E+00,
            0.9900000000000000E+00,
            0.9900000000000000E+00,
            0.9900000000000000E+00,
            0.9900000000000000E+00,
            0.9900000000000000E+00,
            0.9900000000000000E+00,
            0.9900000000000000E+00
        };

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        n_data += 1;

        if (N_MAX < n_data)
        {
            n_data = 0;
            t = 0.0;
            n = 0;
            fx = 0.0;
        }
        else
        {
            t = t_vec[n_data - 1];
            n = n_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }
    }

    public static void log10_values(ref int n_data, ref double x, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG10_VALUES returns some values of the logarithm base 10 function.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Log[x]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 March 2010
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
        const int N_MAX = 20;

        double[] fx_vec =
        {
            -5.0000000000000000000,
            -2.0000000000000000000,
            -1.0000000000000000000,
            -0.69897000433601880479,
            -0.52287874528033756270,
            -0.39794000867203760957,
            -0.30102999566398119521,
            -0.22184874961635636749,
            -0.15490195998574316929,
            -0.096910013008056414359,
            -0.045757490560675125410,
            0.000000000000000000000,
            0.30102999566398119521,
            0.47712125471966243730,
            0.49714987269413385435,
            0.69897000433601880479,
            1.0000000000000000000,
            1.3010299956639811952,
            2.0000000000000000000,
            8.0915149771692704475
        };

        double[] x_vec =
        {
            1.0E-05,
            1.0E-02,
            0.1E+00,
            0.2E+00,
            0.3E+00,
            0.4E+00,
            0.5E+00,
            0.6E+00,
            0.7E+00,
            0.8E+00,
            0.9E+00,
            1.0E+00,
            2.0E+00,
            3.0E+00,
            3.1415926535897932385E+00,
            5.0E+00,
            10.0E+00,
            20.0E+00,
            100.0E+00,
            123456789.0E+00
        };

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        n_data += 1;

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

    public static void logarithmic_integral_values(ref int n_data, ref double x, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOGARITHMIC_INTEGRAL_VALUES returns values of the logarithmic integral LI(X).
        //
        //  Discussion:
        //
        //    The logarithmic integral is defined as:
        //
        //      LI(X) = integral ( 0 <= T <= Z ) dT / log ( T )
        //
        //    The principal value of the integral is taken.  There is a
        //    branch cut discontinuity in the complex plane from -oo to +1.
        //
        //    Abramowitz and Stegun assume 1 < X.
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      LogIntegral[x]
        //
        //    There is a simple relationship with the exponential integral EI:
        //
        //      LI(X) = EI(LN(X))
        //
        //    The function LI(X) provides a good approximation to PI(X),
        //    the number of primes less than or equal to X.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 August 2004
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
        const int N_MAX = 28;

        double[] fx_vec =
        {
            0.0000000000000000E+00,
            -0.3238978959329102E-01,
            -0.8512648672879405E-01,
            -0.1574149028946895E+00,
            -0.2529494192126213E+00,
            -0.3786710430610880E+00,
            -0.5468514142104170E+00,
            -0.7809468775455607E+00,
            -0.1134011957382327E+01,
            -0.1775800683423525E+01,
            -0.2443622553873225E+01,
            -0.3124190050507211E+01,
            -0.2872935510329120E+01,
            -0.2164282524138207E+01,
            -0.1440351296279408E+01,
            -0.6864884538258716E+00,
            0.1250649863152964E+00,
            0.1045163780117493E+01,
            0.2967585095039051E+01,
            0.5253718299558931E+01,
            0.8519716463711059E+01,
            0.1360509217709172E+02,
            0.2193466832805100E+02,
            0.3604254831722944E+02,
            0.6051306533791733E+02,
            0.1037211171690373E+03,
            0.1810780396816945E+03,
            0.3211144156746837E+03
        };

        double[] x_vec =
        {
            0.000000E+00,
            0.100000E+00,
            0.200000E+00,
            0.300000E+00,
            0.400000E+00,
            0.500000E+00,
            0.600000E+00,
            0.700000E+00,
            0.800000E+00,
            0.900000E+00,
            0.950000E+00,
            0.975000E+00,
            0.103125E+01,
            0.106250E+01,
            0.112500E+01,
            0.125000E+01,
            0.150000E+01,
            0.200000E+01,
            0.400000E+01,
            0.800000E+01,
            0.160000E+02,
            0.320000E+02,
            0.640000E+02,
            0.128000E+03,
            0.256000E+03,
            0.512000E+03,
            0.102400E+04,
            0.204800E+04
        };

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        n_data += 1;

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


    public static void polylogarithm_values(ref int n_data, ref int n, ref double z, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYLOGARITHM_VALUES returns some values of the polylogarithm.
        //
        //  Discussion:
        //
        //    The polylogarithm of n and z is defined as
        //
        //      f[n,z] = Sum ( 1 <= k < +oo ) z^k / k^n
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      PolyLog[n,z]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 September 2004
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
        //    Output, ref int N, the exponent of the denominator.
        //
        //    Output, ref double Z, the base of the numerator.
        //
        //    Output, ref double FX, the value of the function.
        //
    {
        const int N_MAX = 12;

        double[] fx_vec =
        {
            0.1644934066848226E+01,
            0.1202056903159594E+01,
            0.1000994575127818E+01,
            0.5822405264650125E+00,
            0.5372131936080402E+00,
            0.5002463206060068E+00,
            0.3662132299770635E+00,
            0.3488278611548401E+00,
            0.3334424797228716E+00,
            0.1026177910993911E+00,
            0.1012886844792230E+00,
            0.1000097826564961E+00
        };

        int[] n_vec =
        {
            2, 3, 10, 2, 3, 10, 2, 3, 10, 2, 3, 10
        };

        double[] z_vec =
        {
            0.1000000000000000E+01,
            0.1000000000000000E+01,
            0.1000000000000000E+01,
            0.5000000000000000E+00,
            0.5000000000000000E+00,
            0.5000000000000000E+00,
            0.3333333333333333E+00,
            0.3333333333333333E+00,
            0.3333333333333333E+00,
            0.1000000000000000E+00,
            0.1000000000000000E+00,
            0.1000000000000000E+00
        };

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        n_data += 1;

        if (N_MAX < n_data)
        {
            n_data = 0;
            n = 0;
            z = 0.0;
            fx = 0.0;
        }
        else
        {
            n = n_vec[n_data - 1];
            z = z_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }
    }
}