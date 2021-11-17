namespace Burkardt.Values;

public static class Cosine
{
    public static void ci_values(ref int n_data, ref double x, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CI_VALUES returns some values of the cosine integral function.
        //
        //  Discussion:
        //
        //    The cosine integral is defined by
        //
        //      CI(X) = - integral ( X <= T < +oo ) ( cos ( T ) ) / T  dT
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      CosIntegral[x]
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
        const int N_MAX = 16;

        double[] fx_vec =
        {
            -0.1777840788066129E+00,
            -0.2227070695927976E-01,
            0.1005147070088978E+00,
            0.1982786159524672E+00,
            0.2760678304677729E+00,
            0.3374039229009681E+00,
            0.4204591828942405E+00,
            0.4620065850946773E+00,
            0.4717325169318778E+00,
            0.4568111294183369E+00,
            0.4229808287748650E+00,
            0.2858711963653835E+00,
            0.1196297860080003E+00,
            -0.3212854851248112E-01,
            -0.1409816978869304E+00,
            -0.1934911221017388E+00
        };

        double[] x_vec =
        {
            0.5E+00,
            0.6E+00,
            0.7E+00,
            0.8E+00,
            0.9E+00,
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
            4.5E+00
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

    public static void cin_values(ref int n_data, ref double x, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIN_VALUES returns some values of the alternate cosine integral function.
        //
        //  Discussion:
        //
        //    The alternate cosine integral is defined by
        //
        //      CIN(X) = gamma + log(X) + integral ( 0 <= T <= X ) ( cos ( T ) - 1 ) / T  dT
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      EulerGamma + Log[x] - CosIntegral[x]
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
        //    Output, ref double X, the argument of the function.
        //
        //    Output, ref double FX, the value of the function.
        //
    {
        const int N_MAX = 16;

        double[] fx_vec =
        {
            0.6185256314820045E-01,
            0.8866074809482194E-01,
            0.1200260139539026E+00,
            0.1557934976348559E+00,
            0.1957873187759337E+00,
            0.2398117420005647E+00,
            0.3390780388012470E+00,
            0.4516813164280685E+00,
            0.5754867772153906E+00,
            0.7081912003853150E+00,
            0.8473820166866132E+00,
            0.1207635200410304E+01,
            0.1556198167561642E+01,
            0.1862107181909382E+01,
            0.2104491723908354E+01,
            0.2274784183779546E+01
        };

        double[] x_vec =
        {
            0.5E+00,
            0.6E+00,
            0.7E+00,
            0.8E+00,
            0.9E+00,
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
            4.5E+00
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

    public static void cinh_values(ref int n_data, ref double x, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CINH_VALUES returns some values of the alternate hyperbolic cosine integral.
        //
        //  Discussion:
        //
        //    The alternate hyperbolic cosine integral is defined by
        //
        //      CINH(X) = integral ( 0 <= T < X ) ( cosh ( T ) - 1 ) / T  dT
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Integrate [ ( Cosh[t] - 1 ) / t, { t, 0, x } ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 March 2010
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
        //    Input/output, ref int N_DATA.  The user sets N_DATA to 0
        //    before the first call.  On each call, the routine increments N_DATA by 1,
        //    and returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, ref double X, the argument of the function.
        //
        //    Output, ref double FX, the value of the function.
        //
    {
        const int N_MAX = 17;

        double[] fx_vec =
        {
            0.00000000000000000,
            0.06315467070191883,
            0.09136085223843649,
            0.1250284547325902,
            0.1643278712460683,
            0.2094587379417273,
            0.2606512760786754,
            0.3823047024751071,
            0.5318061742668980,
            0.7122865135136963,
            0.9275748842583805,
            1.182304077185436,
            2.030919091578478,
            3.284564141195967,
            5.129213294250493,
            7.850037532801762,
            11.88451858691463
        };
        double[] x_vec =
        {
            0.0,
            0.5,
            0.6,
            0.7,
            0.8,
            0.9,
            1.0,
            1.2,
            1.4,
            1.6,
            1.8,
            2.0,
            2.5,
            3.0,
            3.5,
            4.0,
            4.5
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


    public static void cos_values(ref int n_data, ref double x, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COS_VALUES returns some values of the cosine function.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Cos[x]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 June 2007
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
        const int N_MAX = 13;

        double[] fx_vec =
        {
            1.0000000000000000000,
            0.96592582628906828675,
            0.87758256189037271612,
            0.86602540378443864676,
            0.70710678118654752440,
            0.54030230586813971740,
            0.50000000000000000000,
            0.00000000000000000000,
            -0.41614683654714238700,
            -0.98999249660044545727,
            -1.0000000000000000000,
            -0.65364362086361191464,
            0.28366218546322626447
        };

        double[] x_vec =
        {
            0.0000000000000000000,
            0.26179938779914943654,
            0.50000000000000000000,
            0.52359877559829887308,
            0.78539816339744830962,
            1.0000000000000000000,
            1.0471975511965977462,
            1.5707963267948966192,
            2.0000000000000000000,
            3.0000000000000000000,
            3.1415926535897932385,
            4.0000000000000000000,
            5.0000000000000000000
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

    public static void cos_degree_values(ref int n_data, ref double x, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COS_DEGREE_VALUES: values of the cosine function for degree arguments.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Cos[x Degree]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 March 2010
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
        const int N_MAX = 22;

        double[] fx_vec =
        {
            0.99619469809174553230,
            1.0000000000000000000,
            0.99984769515639123916,
            0.99939082701909573001,
            0.99862953475457387378,
            0.99756405025982424761,
            0.99619469809174553230,
            0.98480775301220805937,
            0.96592582628906828675,
            0.86602540378443864676,
            0.70710678118654752440,
            0.50000000000000000000,
            0.25881904510252076235,
            0.087155742747658173558,
            0.069756473744125300776,
            0.052335956242943832722,
            0.034899496702500971646,
            0.017452406437283512819,
            0.000000000000000000000,
            -0.017452406437283512819,
            -0.25881904510252076235,
            -1.0000000000000000000
        };
        double[] x_vec =
        {
            -5.0,
            0.0,
            1.0,
            2.0,
            3.0,
            4.0,
            5.0,
            10.0,
            15.0,
            30.0,
            45.0,
            60.0,
            75.0,
            85.0,
            86.0,
            87.0,
            88.0,
            89.0,
            90.0,
            91.0,
            105.0,
            180.0
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

    public static void cos_power_int_values(ref int n_data, ref double a, ref double b, ref int n,
            ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COS_POWER_INT_VALUES returns some values of the sine power integral.
        //
        //  Discussion:
        //
        //    The function has the form
        //
        //      COS_POWER_INT(A,B,N) = integral ( A <= T <= B ) ( cos(T) )^N dt
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Integrate [ ( Cos[x] )^n, { x, a, b } ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 March 2012
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
        //    Output, ref double A, &B, the limits of integration.
        //
        //    Output, ref int N, the power.
        //
        //    Output, ref double FX, the value of the function.
        //
    {
        const int N_MAX = 11;

        double[] a_vec =
        {
            0.00E+00,
            0.00E+00,
            0.00E+00,
            0.00E+00,
            0.00E+00,
            0.00E+00,
            0.00E+00,
            0.00E+00,
            0.00E+00,
            0.00E+00,
            0.00E+00
        };

        double[] b_vec =
        {
            3.141592653589793,
            3.141592653589793,
            3.141592653589793,
            3.141592653589793,
            3.141592653589793,
            3.141592653589793,
            3.141592653589793,
            3.141592653589793,
            3.141592653589793,
            3.141592653589793,
            3.141592653589793
        };

        double[] fx_vec =
        {
            3.141592653589793,
            0.0,
            1.570796326794897,
            0.0,
            1.178097245096172,
            0.0,
            0.9817477042468104,
            0.0,
            0.8590292412159591,
            0.0,
            0.7731263170943632
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
            10
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
            a = 0.0;
            b = 0.0;
            n = 0;
            fx = 0.0;
        }
        else
        {
            a = a_vec[n_data - 1];
            b = b_vec[n_data - 1];
            n = n_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }
    }

    public static void cosh_values(ref int n_data, ref double x, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COSH_VALUES returns some values of the hyperbolic cosine function.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Cosh[x]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2007
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
        const int N_MAX = 18;

        double[] fx_vec =
        {
            74.209948524787844444,
            1.5430806348152437785,
            1.0000000000000000000,
            1.0050041680558035990,
            1.0200667556190758463,
            1.0453385141288604850,
            1.0810723718384548093,
            1.1276259652063807852,
            1.1854652182422677038,
            1.2551690056309430182,
            1.3374349463048445980,
            1.4330863854487743878,
            1.5430806348152437785,
            3.7621956910836314596,
            10.067661995777765842,
            27.308232836016486629,
            74.209948524787844444,
            11013.232920103323140
        };

        double[] x_vec =
        {
            -5.0,
            -1.0,
            0.0,
            0.1,
            0.2,
            0.3,
            0.4,
            0.5,
            0.6,
            0.7,
            0.8,
            0.9,
            1.0,
            2.0,
            3.0,
            4.0,
            5.0,
            10.0
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

}