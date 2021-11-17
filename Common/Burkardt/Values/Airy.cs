using System.Numerics;

namespace Burkardt.Values;

public static class Airy
{
    public static void airy_ai_values(ref int n_data, ref double x, ref double ai)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    AIRY_AI_VALUES returns some values of the Airy Ai(x) function.
        //
        //  Discussion:
        //
        //    The Airy functions Ai(X) and Bi(X) are a pair of linearly independent
        //    solutions of the differential equation:
        //
        //      W'' - X * W = 0;
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      AiryAi[x]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 August 2004
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
        //    Output, ref double AI, the value of the Airy AI function.
        //
    {
        const int N_MAX = 11;

        double[] ai_vec
                =
                {
                    0.3550280538878172E+00,
                    0.3292031299435381E+00,
                    0.3037031542863820E+00,
                    0.2788064819550049E+00,
                    0.2547423542956763E+00,
                    0.2316936064808335E+00,
                    0.2098000616663795E+00,
                    0.1891624003981501E+00,
                    0.1698463174443649E+00,
                    0.1518868036405444E+00,
                    0.1352924163128814E+00
                }
            ;

        double[] x_vec
                =
                {
                    0.0E+00,
                    0.1E+00,
                    0.2E+00,
                    0.3E+00,
                    0.4E+00,
                    0.5E+00,
                    0.6E+00,
                    0.7E+00,
                    0.8E+00,
                    0.9E+00,
                    1.0E+00
                }
            ;

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
            ai = 0.0;
        }
        else
        {
            x = x_vec[n_data - 1];
            ai = ai_vec[n_data - 1];
        }
    }

    public static void airy_ai_int_values(ref int n_data, ref double x, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    AIRY_AI_INT_VALUES returns some values of the integral of the Airy function.
        //
        //  Discussion:
        //
        //    The function is defined by:
        //
        //      AIRY_AI_INT(x) = Integral ( 0 <= t <= x ) Ai(t) dt
        //
        //    The data was reported by McLeod.
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
        const int N_MAX = 20;

        double[] fx_vec
                =
                {
                    -0.75228838916610124300E+00,
                    -0.57348350185854889466E+00,
                    -0.76569840313421291743E+00,
                    -0.65181015505382467421E+00,
                    -0.55881974894471876922E+00,
                    -0.56902352870716815309E+00,
                    -0.47800749642926168100E+00,
                    -0.46567398346706861416E+00,
                    -0.96783140945618013679E-01,
                    -0.34683049857035607494E-03,
                    0.34658366917927930790E-03,
                    0.27657581846051227124E-02,
                    0.14595330491185717833E+00,
                    0.23631734191710977960E+00,
                    0.33289264538612212697E+00,
                    0.33318759129779422976E+00,
                    0.33332945170523851439E+00,
                    0.33333331724248357420E+00,
                    0.33333333329916901594E+00,
                    0.33333333333329380187E+00
                }
            ;

        double[] x_vec
                =
                {
                    -12.0000000000E+00,
                    -11.0000000000E+00,
                    -10.0000000000E+00,
                    -9.5000000000E+00,
                    -9.0000000000E+00,
                    -6.5000000000E+00,
                    -4.0000000000E+00,
                    -1.0000000000E+00,
                    -0.2500000000E+00,
                    -0.0009765625E+00,
                    0.0009765625E+00,
                    0.0078125000E+00,
                    0.5000000000E+00,
                    1.0000000000E+00,
                    4.0000000000E+00,
                    4.5000000000E+00,
                    6.0000000000E+00,
                    8.0000000000E+00,
                    10.0000000000E+00,
                    12.0000000000E+00
                }
            ;

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

    public static void airy_ai_prime_values(ref int n_data, ref double x, ref double aip)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    AIRY_AI_PRIME_VALUES returns some values of the Airy function Ai'(x).
        //
        //  Discussion:
        //
        //    The Airy functions Ai(X) and Bi(X) are a pair of linearly independent
        //    solutions of the differential equation:
        //
        //      W'' - X * W = 0;
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      AiryAiPrime[x]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 August 2004
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
        //    Output, ref double AIP, the derivative of the Airy AI function.
        //
    {
        const int N_MAX = 11;

        double[] aip_vec
                =
                {
                    -0.2588194037928068E+00,
                    -0.2571304219075862E+00,
                    -0.2524054702856195E+00,
                    -0.2451463642190548E+00,
                    -0.2358320344192082E+00,
                    -0.2249105326646839E+00,
                    -0.2127932593891585E+00,
                    -0.1998511915822805E+00,
                    -0.1864128638072717E+00,
                    -0.1727638434616347E+00,
                    -0.1591474412967932E+00
                }
            ;

        double[] x_vec
                =
                {
                    0.0E+00,
                    0.1E+00,
                    0.2E+00,
                    0.3E+00,
                    0.4E+00,
                    0.5E+00,
                    0.6E+00,
                    0.7E+00,
                    0.8E+00,
                    0.9E+00,
                    1.0E+00
                }
            ;

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
            aip = 0.0;
        }
        else
        {
            x = x_vec[n_data - 1];
            aip = aip_vec[n_data - 1];
        }
    }

    public static void airy_bi_values(ref int n_data, ref double x, ref double bi)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    AIRY_BI_VALUES returns some values of the Airy Bi(x) function.
        //
        //  Discussion:
        //
        //    The Airy functions Ai(X) and Bi(X) are a pair of linearly independent
        //    solutions of the differential equation:
        //
        //      W'' - X * W = 0;
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      AiryBi[x]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 August 2004
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
        //    Output, ref double BI, the value of the Airy BI function.
        //
    {
        const int N_MAX = 11;

        double[] bi_vec
                =
                {
                    0.6149266274460007E+00,
                    0.6598616901941892E+00,
                    0.7054642029186612E+00,
                    0.7524855850873156E+00,
                    0.8017730000135972E+00,
                    0.8542770431031555E+00,
                    0.9110633416949405E+00,
                    0.9733286558781659E+00,
                    0.1042422171231561E+01,
                    0.1119872813134447E+01,
                    0.1207423594952871E+01
                }
            ;

        double[] x_vec
                =
                {
                    0.0E+00,
                    0.1E+00,
                    0.2E+00,
                    0.3E+00,
                    0.4E+00,
                    0.5E+00,
                    0.6E+00,
                    0.7E+00,
                    0.8E+00,
                    0.9E+00,
                    1.0E+00
                }
            ;

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
            bi = 0.0;
        }
        else
        {
            x = x_vec[n_data - 1];
            bi = bi_vec[n_data - 1];
        }
    }

    public static void airy_bi_int_values(ref int n_data, ref double x, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    AIRY_BI_INT_VALUES returns some values of the integral of the Airy function.
        //
        //  Discussion:
        //
        //    The function is defined by:
        //
        //      AIRY_BI_INT(x) = Integral ( 0 <= t <= x ) Bi(t) dt
        //
        //    The data was reported by McLeod.
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
        const int N_MAX = 20;

        double[] fx_vec
                =
                {
                    0.17660819031554631869E-01,
                    -0.15040424806140020451E-01,
                    0.14756446293227661920E-01,
                    -0.11847304264848446271E+00,
                    -0.64916741266165856037E-01,
                    0.97260832464381044540E-01,
                    0.50760058495287539119E-01,
                    -0.37300500963429492179E+00,
                    -0.13962988442666578531E+00,
                    -0.12001735266723296160E-02,
                    0.12018836117890354598E-02,
                    0.36533846550952011043E+00,
                    0.87276911673800812196E+00,
                    0.48219475263803429675E+02,
                    0.44006525804904178439E+06,
                    0.17608153976228301458E+07,
                    0.73779211705220007228E+07,
                    0.14780980310740671617E+09,
                    0.97037614223613433849E+11,
                    0.11632737638809878460E+15
                }
            ;

        double[] x_vec
                =
                {
                    -12.0000000000E+00,
                    -10.0000000000E+00,
                    -8.0000000000E+00,
                    -7.5000000000E+00,
                    -7.0000000000E+00,
                    -6.5000000000E+00,
                    -4.0000000000E+00,
                    -1.0000000000E+00,
                    -0.2500000000E+00,
                    -0.0019531250E+00,
                    0.0019531250E+00,
                    0.5000000000E+00,
                    1.0000000000E+00,
                    4.0000000000E+00,
                    8.0000000000E+00,
                    8.5000000000E+00,
                    9.0000000000E+00,
                    10.0000000000E+00,
                    12.0000000000E+00,
                    14.0000000000E+00
                }
            ;

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

    public static void airy_bi_prime_values(ref int n_data, ref double x, ref double bip)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    AIRY_BI_PRIME_VALUES returns some values of the Airy function Bi'(x).
        //
        //  Discussion:
        //
        //    The Airy functions Ai(X) and Bi(X) are a pair of linearly independent
        //    solutions of the differential equation:
        //
        //      W'' - X * W = 0;
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      AiryBiPrime[x]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 August 2004
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
        //    Output, ref double BIP, the derivative of the Airy BI function.
        //
    {
        const int N_MAX = 11;

        double[] bip_vec
                =
                {
                    0.4482883573538264E+00,
                    0.4515126311496465E+00,
                    0.4617892843621509E+00,
                    0.4800490287524480E+00,
                    0.5072816760506224E+00,
                    0.5445725641405923E+00,
                    0.5931444786342857E+00,
                    0.6544059191721400E+00,
                    0.7300069016152518E+00,
                    0.8219038903072090E+00,
                    0.9324359333927756E+00
                }
            ;

        double[] x_vec
                =
                {
                    0.0E+00,
                    0.1E+00,
                    0.2E+00,
                    0.3E+00,
                    0.4E+00,
                    0.5E+00,
                    0.6E+00,
                    0.7E+00,
                    0.8E+00,
                    0.9E+00,
                    1.0E+00
                }
            ;

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
            bip = 0.0;
        }
        else
        {
            x = x_vec[n_data - 1];
            bip = bip_vec[n_data - 1];
        }
    }

    public static void airy_cai_values(ref int n_data, ref Complex x,
            ref Complex cai )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    AIRY_CAI_VALUES returns some values of the Airy Ai(x) for complex argument.
        //
        //  Discussion:
        //
        //    The Airy functions Ai(X) and Bi(X) are a pair of linearly independent
        //    solutions of the differential equation:
        //
        //      W'' - X * W = 0;
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      AiryAi[x]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 April 2007
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
        //    Output, complex <double> &X, the argument of the function.
        //
        //    Output, complex <double> &CAI, the value of the Airy AI function.
        //
    {
        const int N_MAX = 10;

        Complex[] cai_vec
                =

                {
                    new(0.1352924163128814, +0.0000000000000000),
                    new(0.1433824486882056, -0.1092193342707378),
                    new(0.2215404472324631, -0.2588711788891803),
                    new(0.4763929771766866, -0.3036484220291284),
                    new(0.5983692170633874, -0.08154602160771214),
                    new(0.5355608832923521, +0.00000000000000000),
                    new(0.5983692170633874, +0.08154602160771214),
                    new(0.4763929771766866, +0.3036484220291284),
                    new(0.2215404472324631, +0.2588711788891803),
                    new(0.1433824486882056, +0.1092193342707378)
                }
            ;

        Complex[] x_vec
                =

                {
                    new(1.0000000000000000, +0.0000000000000000),
                    new(0.8090169943749474, +0.5877852522924731),
                    new(0.3090169943749474, +0.9510565162951536),
                    new(-0.3090169943749474, +0.9510565162951536),
                    new(-0.8090169943749474, +0.5877852522924731),
                    new(-1.0000000000000000, +0.0000000000000000),
                    new(-0.8090169943749474, -0.5877852522924731),
                    new(-0.3090169943749474, -0.9510565162951536),
                    new(0.3090169943749474, -0.9510565162951536),
                    new(0.8090169943749474, -0.5877852522924731)
                }
            ;

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
            cai = 0.0;
        }
        else
        {
            x = x_vec[n_data - 1];
            cai = cai_vec[n_data - 1];
        }
    }

    public static void airy_cbi_values(ref int n_data, ref Complex x,
            ref Complex cbi )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    AIRY_CBI_VALUES returns some values of the Airy Bi(x) for complex argument.
        //
        //  Discussion:
        //
        //    The Airy functions Ai(X) and Bi(X) are a pair of linearly independent
        //    solutions of the differential equation:
        //
        //      W'' - X * W = 0;
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      AiryAi[x]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 April 2007
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
        //    Output, complex <double> &X, the argument of the function.
        //
        //    Output, complex <double> &CBI, the value of the Airy BI function.
        //
    {
        const int N_MAX = 10;

        Complex[] cbi_vec
                =

                {
                    new(1.207423594952871, +0.0000000000000000),
                    new(0.9127160108293936, +0.3800456133135556),
                    new(0.6824453575635721, +0.3343047153635002),
                    new(0.5726265660086474, +0.3988641086982559),
                    new(0.2511841251049547, +0.3401447690712719),
                    new(0.1039973894969446, +0.0000000000000000),
                    new(0.2511841251049547, -0.3401447690712719),
                    new(0.5726265660086474, -0.3988641086982559),
                    new(0.6824453575635721, -0.3343047153635002),
                    new(0.9127160108293936, -0.3800456133135556)
                }
            ;

        Complex[] x_vec
                =

                {
                    new(1.0000000000000000, +0.0000000000000000),
                    new(0.8090169943749474, +0.5877852522924731),
                    new(0.3090169943749474, +0.9510565162951536),
                    new(-0.3090169943749474, +0.9510565162951536),
                    new(-0.8090169943749474, +0.5877852522924731),
                    new(-1.0000000000000000, +0.0000000000000000),
                    new(-0.8090169943749474, -0.5877852522924731),
                    new(-0.3090169943749474, -0.9510565162951536),
                    new(0.3090169943749474, -0.9510565162951536),
                    new(0.8090169943749474, -0.5877852522924731)
                }
            ;

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
            cbi = 0.0;
        }
        else
        {
            x = x_vec[n_data - 1];
            cbi = cbi_vec[n_data - 1];
        }
    }

    public static void airy_gi_values(ref int n_data, ref double x, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    AIRY_GI_VALUES returns some values of the Airy Gi function.
        //
        //  Discussion:
        //
        //    The function is defined by:
        //
        //      AIRY_GI(x) = Integral ( 0 <= t < +oo ) sin ( x*t+t^3/3) dt / pi
        //
        //    The data was reported by McLeod.
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
        const int N_MAX = 20;

        double[] fx_vec
                =
                {
                    0.20468308070040542435E+00,
                    0.18374662832557904078E+00,
                    -0.11667221729601528265E+00,
                    0.31466934902729557596E+00,
                    -0.37089040722426257729E+00,
                    -0.25293059772424019694E+00,
                    0.28967410658692701936E+00,
                    -0.34644836492634090590E+00,
                    0.28076035913873049496E+00,
                    0.21814994508094865815E+00,
                    0.20526679000810503329E+00,
                    0.22123695363784773258E+00,
                    0.23521843981043793760E+00,
                    0.82834303363768729338E-01,
                    0.45757385490989281893E-01,
                    0.44150012014605159922E-01,
                    0.39951133719508907541E-01,
                    0.35467706833949671483E-01,
                    0.31896005100679587981E-01,
                    0.26556892713512410405E-01
                }
            ;

        double[] x_vec
                =
                {
                    -0.0019531250E+00,
                    -0.1250000000E+00,
                    -1.0000000000E+00,
                    -4.0000000000E+00,
                    -8.0000000000E+00,
                    -8.2500000000E+00,
                    -9.0000000000E+00,
                    -10.0000000000E+00,
                    -11.0000000000E+00,
                    -13.0000000000E+00,
                    0.0019531250E+00,
                    0.1250000000E+00,
                    1.0000000000E+00,
                    4.0000000000E+00,
                    7.0000000000E+00,
                    7.2500000000E+00,
                    8.0000000000E+00,
                    9.0000000000E+00,
                    10.0000000000E+00,
                    12.0000000000E+00
                }
            ;

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

    public static void airy_hi_values(ref int n_data, ref double x, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    AIRY_HI_VALUES returns some values of the Airy Hi function.
        //
        //  Discussion:
        //
        //    The function is defined by:
        //
        //      AIRY_HI(x) = Integral ( 0 <= t < +oo ) exp(x*t-t^3/3) dt / pi
        //
        //    The data was reported by McLeod.
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
        const int N_MAX = 20;

        double[] fx_vec
                =
                {
                    0.40936798278458884024E+00,
                    0.37495291608048868619E+00,
                    0.22066960679295989454E+00,
                    0.77565356679703713590E-01,
                    0.39638826473124717315E-01,
                    0.38450072575004151871E-01,
                    0.35273216868317898556E-01,
                    0.31768535282502272742E-01,
                    0.28894408288051391369E-01,
                    0.24463284011678541180E-01,
                    0.41053540139998941517E+00,
                    0.44993502381204990817E+00,
                    0.97220515514243332184E+00,
                    0.83764237105104371193E+02,
                    0.80327744952044756016E+05,
                    0.15514138847749108298E+06,
                    0.11995859641733262114E+07,
                    0.21472868855967642259E+08,
                    0.45564115351632913590E+09,
                    0.32980722582904761929E+12
                }
            ;

        double[] x_vec
                =
                {
                    -0.0019531250E+00,
                    -0.1250000000E+00,
                    -1.0000000000E+00,
                    -4.0000000000E+00,
                    -8.0000000000E+00,
                    -8.2500000000E+00,
                    -9.0000000000E+00,
                    -10.0000000000E+00,
                    -11.0000000000E+00,
                    -13.0000000000E+00,
                    0.0019531250E+00,
                    0.1250000000E+00,
                    1.0000000000E+00,
                    4.0000000000E+00,
                    7.0000000000E+00,
                    7.2500000000E+00,
                    8.0000000000E+00,
                    9.0000000000E+00,
                    10.0000000000E+00,
                    12.0000000000E+00
                }
            ;

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