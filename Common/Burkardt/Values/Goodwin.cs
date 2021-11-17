namespace Burkardt.Values;

public static class Goodwin
{
    public static void goodwin_values(ref int n_data, ref double x, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GOODWIN_VALUES returns some values of the Goodwin and Staton function.
        //
        //  Discussion:
        //
        //    The function is defined by:
        //
        //      GOODWIN(x) = Integral ( 0 <= t < +oo ) exp ( -t^2 ) / ( t + x ) dt
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
        const int N_MAX = 20;

        double[] fx_vec =
        {
            0.59531540040441651584E+01,
            0.45769601268624494109E+01,
            0.32288921331902217638E+01,
            0.19746110873568719362E+01,
            0.96356046208697728563E+00,
            0.60513365250334458174E+00,
            0.51305506459532198016E+00,
            0.44598602820946133091E+00,
            0.37344458206879749357E+00,
            0.35433592884953063055E+00,
            0.33712156518881920994E+00,
            0.29436170729362979176E+00,
            0.25193499644897222840E+00,
            0.22028778222123939276E+00,
            0.19575258237698917033E+00,
            0.17616303166670699424E+00,
            0.16015469479664778673E+00,
            0.14096116876193391066E+00,
            0.13554987191049066274E+00,
            0.11751605060085098084E+00
        };

        double[] x_vec =
        {
            0.0019531250E+00,
            0.0078125000E+00,
            0.0312500000E+00,
            0.1250000000E+00,
            0.5000000000E+00,
            1.0000000000E+00,
            1.2500000000E+00,
            1.5000000000E+00,
            1.8750000000E+00,
            2.0000000000E+00,
            2.1250000000E+00,
            2.5000000000E+00,
            3.0000000000E+00,
            3.5000000000E+00,
            4.0000000000E+00,
            4.5000000000E+00,
            5.0000000000E+00,
            5.7500000000E+00,
            6.0000000000E+00,
            7.0000000000E+00
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