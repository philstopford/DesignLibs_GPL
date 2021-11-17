namespace Burkardt.Elliptic;

public class FK
{
    public static double evaluate(double k)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPTIC_FK evaluates the complete elliptic integral F(K).
        //
        //  Discussion:
        //
        //    The value is computed using Carlson elliptic integrals:
        //
        //      F(k) = RF ( 0, 1-k^2, 1 ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 June 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double K, the argument.
        //
        //    Output, double ELLIPTIC_FK, the function value.
        //
    {
        double errtol;
        int ierr = 0;
        double value = 0;
        double x;
        double y;
        double z;

        x = 0.0;
        y = (1.0 - k) * (1.0 + k);
        z = 1.0;
        errtol = 1.0E-03;

        value = Integral.rf(x, y, z, errtol, ref ierr);

        return value;
    }

    public static void values(ref int n_data, ref double x, ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPTIC_FK_VALUES returns values of the complete elliptic integral F(K).
        //
        //  Discussion:
        //
        //    This is one form of what is sometimes called the complete elliptic
        //    integral of the first kind.
        //
        //    The function is defined by the formula:
        //
        //      F(K) = integral ( 0 <= T <= PI/2 )
        //        dT / sqrt ( 1 - K^2 * sin ( T )^2 )
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      EllipticK[k^2]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 May 2018
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
        //    Output, double &X, the argument of the function.
        //
        //    Output, double &FX, the value of the function.
        //
    {
        const int N_MAX = 20;

        double[] fx_vec =
            {
                1.570796326794897E+00,
                1.591003453790792E+00,
                1.612441348720219E+00,
                1.635256732264580E+00,
                1.659623598610528E+00,
                1.685750354812596E+00,
                1.713889448178791E+00,
                1.744350597225613E+00,
                1.777519371491253E+00,
                1.813883936816983E+00,
                1.854074677301372E+00,
                1.898924910271554E+00,
                1.949567749806026E+00,
                2.007598398424376E+00,
                2.075363135292469E+00,
                2.156515647499643E+00,
                2.257205326820854E+00,
                2.389016486325580E+00,
                2.578092113348173E+00,
                2.908337248444552E+00
            }
            ;

        double[] x_vec =
            {
                0.0000000000000000E+00,
                0.2236067977499790E+00,
                0.3162277660168379E+00,
                0.3872983346207417E+00,
                0.4472135954999579E+00,
                0.5000000000000000E+00,
                0.5477225575051661E+00,
                0.5916079783099616E+00,
                0.6324555320336759E+00,
                0.6708203932499369E+00,
                0.7071067811865476E+00,
                0.7416198487095663E+00,
                0.7745966692414834E+00,
                0.8062257748298550E+00,
                0.8366600265340756E+00,
                0.8660254037844386E+00,
                0.8944271909999159E+00,
                0.9219544457292888E+00,
                0.9486832980505138E+00,
                0.9746794344808963E+00
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