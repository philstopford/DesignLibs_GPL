using System;

namespace Burkardt.Elliptic;

public static class FA_inc
{
    public static double evaluate(double phi, double a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPTIC_INC_FA evaluates the incomplete elliptic integral F(PHI,A).
        //
        //  Discussion:
        //
        //    The value is computed using Carlson elliptic integrals:
        //
        //      k = sin ( a * Math.PI / 180 )
        //      F(phi,k) = sin(phi) * RF ( cos^2 ( phi ), 1-k^2 sin^2 ( phi ), 1 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 June 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double PHI, A, the arguments.
        //    0 <= PHI <= PI/2.
        //    0 <= sin^2 ( A * Math.PI / 180 ) * sin^2(PHI) <= 1.
        //
        //    Output, double ELLIPTIC_INC_FA, the function value.
        //
    {
        int ierr = 0;

        double k = Math.Sin(a * Math.PI / 180.0);

        double cp = Math.Cos(phi);
        double sp = Math.Sin(phi);
        double x = cp * cp;
        double y = (1.0 - k * sp) * (1.0 + k * sp);
        const double z = 1.0;
        const double errtol = 1.0E-03;

        double value = Integral.rf(x, y, z, errtol, ref ierr);

        if (ierr != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("ELLIPTIC_INC_FA - Fatal error!");
            Console.WriteLine("  RF returned IERR = " + ierr + "");
            return 1;
        }

        value = sp * value;

        return value;
    }

    public static void values(ref int n_data, ref double phi, ref double a, ref double fa )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPTIC_INC_FA_VALUES: values of the incomplete elliptic integral F(PHI,A).
        //
        //  Discussion:
        //
        //    This is the incomplete elliptic integral of the first kind.
        //
        //      F(PHI,A) = integral ( 0 <= T <= PHI ) 
        //        dT / sqrt ( 1 - sin^2 ( A ) * sin^2 ( T ) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 June 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    US Department of Commerce, 1964.
        //
        //    Stephen Wolfram,
        //    The Mathematica Book,
        //    Fourth Edition,
        //    Wolfram Media / Cambridge University Press, 1999.
        //
        //  Parameters:
        //
        //    Input/output, int &N_DATA.  The user sets N_DATA to 0 
        //    before the first call.  On each call, the routine increments N_DATA by 1, 
        //    and returns the corresponding data when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, double &PHI, &A, the arguments.
        //
        //    Output, double &FA, the function value.
        //
    {
        const int N_MAX = 20;

        double[] a_vec =
            {
                123.0821233267548,
                11.26931745051486,
                -94.88806452075445,
                -99.71407853545323,
                57.05881039324191,
                -19.71363287074183,
                56.31230299738043,
                -91.55605346417718,
                -27.00654574696468,
                -169.2293728595904,
                61.96859564803047,
                -158.7324398933148,
                105.0883958999383,
                -48.95883872360177,
                -42.58568835110901,
                11.65603284687828,
                -8.398113719173338,
                17.69362213019626,
                73.8803420626852,
                -69.82492339645128
            }
            ;

        double[] fa_vec =
            {
                0.3478806460316299,
                1.313180577009584,
                0.7037956689264326,
                0.4157626844675118,
                0.06888475483285136,
                0.09697816754845832,
                0.6605394722518515,
                1.82758346036751,
                1.482258783392487,
                0.1485295339221232,
                1.753800062701494,
                0.193528896465351,
                0.4199100508706138,
                0.1790836490491233,
                1.446048832279763,
                1.094097652100984,
                1.358947908427035,
                1.46400078231538,
                0.3009092014525799,
                0.6621341112075102
            }
            ;

        double[] phi_vec =
            {
                0.3430906586047127,
                1.302990057703935,
                0.6523628380743488,
                0.4046022501376546,
                0.06884642871852312,
                0.0969609046794745,
                0.630370432896175,
                1.252375418911598,
                1.409796082144801,
                0.1485105463502483,
                1.349466184634646,
                0.1933711786970301,
                0.4088829927466769,
                0.1785430666405224,
                1.292588374416351,
                1.087095515757691,
                1.352794600489329,
                1.432530166308616,
                0.2968093345769761,
                0.6235880396594726
            }
            ;

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        if (N_MAX <= n_data)
        {
            n_data = 0;
            a = 0.0;
            fa = 0.0;
            phi = 0.0;
        }
        else
        {
            a = a_vec[n_data];
            fa = fa_vec[n_data];
            phi = phi_vec[n_data];
            n_data += 1;
        }
    }
}