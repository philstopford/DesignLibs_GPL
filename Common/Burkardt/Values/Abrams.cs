namespace Burkardt.Values;

public class Abrams
{
    public static void abram0_values(ref int n_data, ref double x, ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ABRAM0_VALUES returns some values of the Abramowitz0 function.
        //
        //  Discussion:
        //
        //    The function is defined by:
        //
        //      ABRAM0(X) = integral ( 0 <= T < +oo ) exp ( -T * T - X / T ) dT
        //
        //    The data was reported by McLeod.
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
        //    Allan McLeod,
        //    Algorithm 757:
        //    MISCFUN: A software package to compute uncommon special functions,
        //    ACM Transactions on Mathematical Software,
        //    Volume 22, Number 3, September 1996, pages 288-301.
        //
        //  Parameters:
        //
        //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, double &X, the argument of the function.
        //
        //    Output, double &FX, the value of the function.
        //
    {
        const int N_MAX = 20;

        double[] fx_vec =
            {
                0.87377726306985360531E+00,
                0.84721859650456925922E+00,
                0.77288934483988301615E+00,
                0.59684345853450151603E+00,
                0.29871735283675888392E+00,
                0.15004596450516388138E+00,
                0.11114662419157955096E+00,
                0.83909567153151897766E-01,
                0.56552321717943417515E-01,
                0.49876496603033790206E-01,
                0.44100889219762791328E-01,
                0.19738535180254062496E-01,
                0.86193088287161479900E-02,
                0.40224788162540127227E-02,
                0.19718658458164884826E-02,
                0.10045868340133538505E-02,
                0.15726917263304498649E-03,
                0.10352666912350263437E-04,
                0.91229759190956745069E-06,
                0.25628287737952698742E-09
            }
            ;
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
                3.0000000000E+00,
                4.0000000000E+00,
                5.0000000000E+00,
                6.0000000000E+00,
                7.0000000000E+00,
                10.0000000000E+00,
                15.0000000000E+00,
                20.0000000000E+00,
                40.0000000000E+00
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

    public static void abram1_values(ref int n_data, ref double x, ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ABRAM1_VALUES returns some values of the Abramowitz1 function.
        //
        //  Discussion:
        //
        //    The function is defined by:
        //
        //      ABRAM1(x) = integral ( 0 <= t < oo ) t * exp ( -t^2 - x / t ) dt
        //
        //    The data was reported by McLeod.
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
        //    Allan McLeod,
        //    Algorithm 757:
        //    MISCFUN: A software package to compute uncommon special functions,
        //    ACM Transactions on Mathematical Software,
        //    Volume 22, Number 3, September 1996, pages 288-301.
        //
        //  Parameters:
        //
        //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, double &X, the argument of the function.
        //
        //    Output, double &FX, the value of the function.
        //
    {
        const int N_MAX = 20;

        double[] fx_vec =
            {
                0.49828219848799921792E+00,
                0.49324391773047288556E+00,
                0.47431612784691234649E+00,
                0.41095983258760410149E+00,
                0.25317617388227035867E+00,
                0.14656338138597777543E+00,
                0.11421547056018366587E+00,
                0.90026307383483764795E-01,
                0.64088214170742303375E-01,
                0.57446614314166191085E-01,
                0.51581624564800730959E-01,
                0.25263719555776416016E-01,
                0.11930803330196594536E-01,
                0.59270542280915272465E-02,
                0.30609215358017829567E-02,
                0.16307382136979552833E-02,
                0.28371851916959455295E-03,
                0.21122150121323238154E-04,
                0.20344578892601627337E-05,
                0.71116517236209642290E-09
            }
            ;

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
                3.0000000000E+00,
                4.0000000000E+00,
                5.0000000000E+00,
                6.0000000000E+00,
                7.0000000000E+00,
                10.0000000000E+00,
                15.0000000000E+00,
                20.0000000000E+00,
                40.0000000000E+00
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

    public static void abram2_values(ref int n_data, ref double x, ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ABRAM2_VALUES returns some values of the Abramowitz2 function.
        //
        //  Discussion:
        //
        //    The function is defined by:
        //
        //      ABRAM2(x) = Integral ( 0 <= t < +oo ) t^2 * exp( -t^2 - x / t ) dt
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
        //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, double &X, the argument of the function.
        //
        //    Output, double &FX, the value of the function.
        //
    {
        const int N_MAX = 20;

        double[] fx_vec =
            {
                0.44213858162107913430E+00,
                0.43923379545684026308E+00,
                0.42789857297092602234E+00,
                0.38652825661854504406E+00,
                0.26538204413231368110E+00,
                0.16848734838334595000E+00,
                0.13609200032513227112E+00,
                0.11070330027727917352E+00,
                0.82126019995530382267E-01,
                0.74538781999594581763E-01,
                0.67732034377612811390E-01,
                0.35641808698811851022E-01,
                0.17956589956618269083E-01,
                0.94058737143575370625E-02,
                0.50809356204299213556E-02,
                0.28149565414209719359E-02,
                0.53808696422559303431E-03,
                0.44821756380146327259E-04,
                0.46890678427324100410E-05,
                0.20161544850996420504E-08
            }
            ;

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
                3.0000000000E+00,
                4.0000000000E+00,
                5.0000000000E+00,
                6.0000000000E+00,
                7.0000000000E+00,
                10.0000000000E+00,
                15.0000000000E+00,
                20.0000000000E+00,
                40.0000000000E+00
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