namespace Burkardt.Values
{
    public static class arc
    {
        public static void arccos_values(ref int n_data, ref double x, ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ARCCOS_VALUES returns some values of the arc cosine function.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      ArcCos[x]
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
            int N_MAX = 12;

            double[] fx_vec
             = 
            {
                1.6709637479564564156,
                1.5707963267948966192,
                1.4706289056333368229,
                1.3694384060045658278,
                1.2661036727794991113,
                1.1592794807274085998,
                1.0471975511965977462,
                0.92729521800161223243,
                0.79539883018414355549,
                0.64350110879328438680,
                0.45102681179626243254,
                0.00000000000000000000
            }
            ;

            double[] x_vec
             = 
            {
                -0.1,
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
                1.0
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
                fx = 0.0;
            }
            else
            {
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void arccosh_values(ref int n_data, ref double x, ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ARCCOSH_VALUES returns some values of the hyperbolic arc cosine function.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      ArcCosh[x]
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
            int N_MAX = 15;

            double[] fx_vec
             = 
            {
                0.0000000000000000000,
                0.14130376948564857735,
                0.44356825438511518913,
                0.62236250371477866781,
                0.75643291085695958624,
                0.86701472649056510395,
                0.96242365011920689500,
                1.3169578969248167086,
                1.7627471740390860505,
                1.8115262724608531070,
                2.0634370688955605467,
                2.2924316695611776878,
                2.9932228461263808979,
                5.2982923656104845907,
                7.6009022095419886114
            }
            ;

            double[] x_vec
             = 
            {
                1.0,
                1.01,
                1.1,
                1.2,
                1.3,
                1.4,
                1.5,
                2.0,
                3.0,
                3.1415926535897932385,
                4.0,
                5.0,
                10.0,
                100.0,
                1000.0
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
                fx = 0.0;
            }
            else
            {
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void arcsin_values(ref int n_data, ref double x, ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ARCSIN_VALUES returns some values of the arc sine function.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      ArcSin[x]
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
            int N_MAX = 12;

            double[] fx_vec
             = 
            {
                -0.10016742116155979635,
                0.00000000000000000000,
                0.10016742116155979635,
                0.20135792079033079146,
                0.30469265401539750797,
                0.41151684606748801938,
                0.52359877559829887308,
                0.64350110879328438680,
                0.77539749661075306374,
                0.92729521800161223243,
                1.1197695149986341867,
                1.5707963267948966192
            }
            ;

            double[] x_vec
             = 
            {
                -0.1,
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
                1.0
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
                fx = 0.0;
            }
            else
            {
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void arcsinh_values(ref int n_data, ref double x, ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ARCSINH_VALUES returns some values of the hyperbolic arc sine function.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      ArcSinh[x]
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
            int N_MAX = 20;

            double[] fx_vec
             = 
            {
                -2.3124383412727526203,
                -0.88137358701954302523,
                0.00000000000000000000,
                0.099834078899207563327,
                0.19869011034924140647,
                0.29567304756342243910,
                0.39003531977071527608,
                0.48121182505960344750,
                0.56882489873224753010,
                0.65266656608235578681,
                0.73266825604541086415,
                0.80886693565278246251,
                0.88137358701954302523,
                1.4436354751788103425,
                1.8184464592320668235,
                2.0947125472611012942,
                2.3124383412727526203,
                2.9982229502979697388,
                5.2983423656105887574,
                7.6009027095419886115
            }
            ;

            double[] x_vec
             = 
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
                10.0,
                100.0,
                1000.0
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
                fx = 0.0;
            }
            else
            {
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void arctan_values(ref int n_data, ref double x, ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ARCTAN_VALUES returns some values of the arc tangent function.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      ArcTan[x]
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
            int N_MAX = 11;

            double[] f_vec
             = 
            {
                0.00000000000000000000,
                0.24497866312686415417,
                0.32175055439664219340,
                0.46364760900080611621,
                0.78539816339744830962,
                1.1071487177940905030,
                1.2490457723982544258,
                1.3258176636680324651,
                1.3734007669450158609,
                1.4711276743037345919,
                1.5208379310729538578
            }
            ;

            double[] x_vec
             = 
            {
                0.00000000000000000000,
                0.25000000000000000000,
                0.33333333333333333333,
                0.50000000000000000000,
                1.0000000000000000000,
                2.0000000000000000000,
                3.0000000000000000000,
                4.0000000000000000000,
                5.0000000000000000000,
                10.000000000000000000,
                20.000000000000000000
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
                fx = 0.0;
            }
            else
            {
                x = x_vec[n_data - 1];
                fx = f_vec[n_data - 1];
            }
        }

        public static void arctan_int_values(ref int n_data, ref double x, ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ARCTAN_INT_VALUES returns some values of the inverse tangent integral.
        //
        //  Discussion:
        //
        //    The function is defined by:
        //
        //      ARCTAN_INT(x) = Integral ( 0 <= t <= x ) arctan ( t ) / t dt
        //
        //    The data was reported by McLeod.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 August 2004
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

            double[] f_vec
             = 
            {
                0.19531241721588483191E-02,
                -0.39062433772980711281E-02,
                0.78124470192576499535E-02,
                0.15624576181996527280E-01,
                -0.31246610349485401551E-01,
                0.62472911335014397321E-01,
                0.12478419717389654039E+00,
                -0.24830175098230686908E+00,
                0.48722235829452235711E+00,
                0.91596559417721901505E+00,
                0.12749694484943800618E+01,
                -0.15760154034463234224E+01,
                0.24258878412859089996E+01,
                0.33911633326292997361E+01,
                0.44176450919422186583E+01,
                -0.47556713749547247774E+01,
                0.50961912150934111303E+01,
                0.53759175735714876256E+01,
                -0.61649904785027487422E+01,
                0.72437843013083534973E+01
            }
            ;

            double[] x_vec
             = 
            {
                0.0019531250E+00,
                -0.0039062500E+00,
                0.0078125000E+00,
                0.0156250000E+00,
                -0.0312500000E+00,
                0.0625000000E+00,
                0.1250000000E+00,
                -0.2500000000E+00,
                0.5000000000E+00,
                1.0000000000E+00,
                1.5000000000E+00,
                -2.0000000000E+00,
                4.0000000000E+00,
                8.0000000000E+00,
                16.0000000000E+00,
                -20.0000000000E+00,
                25.0000000000E+00,
                30.0000000000E+00,
                -50.0000000000E+00,
                100.0000000000E+00
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
                fx = 0.0;
            }
            else
            {
                x = x_vec[n_data - 1];
                fx = f_vec[n_data - 1];
            }
        }

        public static void arctan2_values(ref int n_data, ref double x, ref double y, ref double fxy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ARCTAN2_VALUES: arc tangent function of two arguments.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      ArcTan[x,y]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 March 2010
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
        //    Output, ref double X, &Y, the arguments of the function.
        //
        //    Output, ref double FXY, the value of the function.
        //
        {
            int N_MAX = 19;

            double[] f_vec
             = 
            {
                -1.5707963267948966192,
                -1.0471975511965977462,
                -0.52359877559829887308,
                0.00000000000000000000,
                0.52359877559829887308,
                1.0471975511965977462,
                1.5707963267948966192,
                2.0943951023931954923,
                2.6179938779914943654,
                3.1415926535897932385,
                -2.6179938779914943654,
                -2.0943951023931954923,
                -1.5707963267948966192,
                -1.0471975511965977462,
                -0.52359877559829887308,
                0.00000000000000000000,
                0.52359877559829887308,
                1.0471975511965977462,
                1.5707963267948966192
            }
            ;

            double[] x_vec
             = 
            {
                0.00000000000000000000,
                0.50000000000000000000,
                0.86602540378443864676,
                1.00000000000000000000,
                0.86602540378443864676,
                0.50000000000000000000,
                0.00000000000000000000,
                -0.50000000000000000000,
                -0.86602540378443864676,
                -1.00000000000000000000,
                -0.86602540378443864676,
                -0.50000000000000000000,
                0.00000000000000000000,
                0.50000000000000000000,
                0.86602540378443864676,
                1.00000000000000000000,
                0.86602540378443864676,
                0.50000000000000000000,
                0.00000000000000000000
            }
            ;

            double[] y_vec
             = 
            {
                -1.00000000000000000000,
                -0.86602540378443864676,
                -0.50000000000000000000,
                0.00000000000000000000,
                0.50000000000000000000,
                0.86602540378443864676,
                1.00000000000000000000,
                0.86602540378443864676,
                0.50000000000000000000,
                0.00000000000000000000,
                -0.50000000000000000000,
                -0.86602540378443864676,
                -1.00000000000000000000,
                -0.86602540378443864676,
                -0.50000000000000000000,
                0.00000000000000000000,
                0.50000000000000000000,
                0.86602540378443864676,
                1.00000000000000000000
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
                y = 0.0;
                fxy = 0.0;
            }
            else
            {
                x = x_vec[n_data - 1];
                y = y_vec[n_data - 1];
                fxy = f_vec[n_data - 1];
            }
        }

        public static void arctanh_values(ref int n_data, ref double x, ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ARCTANH_VALUES returns some values of the hyperbolic arc tangent function.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      ArcTanh[x]
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
            int N_MAX = 15;

            double[] fx_vec
             = 
            {
                -0.54930614433405484570,
                0.00000000000000000000,
                0.0010000003333335333335,
                0.10033534773107558064,
                0.20273255405408219099,
                0.30951960420311171547,
                0.42364893019360180686,
                0.54930614433405484570,
                0.69314718055994530942,
                0.86730052769405319443,
                1.0986122886681096914,
                1.4722194895832202300,
                2.6466524123622461977,
                3.8002011672502000318,
                7.2543286192620472067
            }
            ;

            double[] x_vec
             = 
            {
                -0.5,
                0.0,
                0.001,
                0.1,
                0.2,
                0.3,
                0.4,
                0.5,
                0.6,
                0.7,
                0.8,
                0.9,
                0.99,
                0.999,
                0.999999
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