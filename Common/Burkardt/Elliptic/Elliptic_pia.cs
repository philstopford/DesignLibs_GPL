using System;

namespace Burkardt.Elliptic;

public static class PIA
{
    public static double evaluate(double n, double a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPTIC_PIA evaluates the complete elliptic integral Pi(N,A).
        //
        //  Discussion:
        //
        //    This is one form of what is sometimes called the complete elliptic
        //    integral of the third kind.
        //
        //    The double is defined by the formula:
        //
        //      Pi(N,A) = integral ( 0 <= T <= PI/2 )
        //        dT / (1 - N sin^2(T) ) sqrt ( 1 - sin^2(A) * sin ( T )^2 )
        //
        //    In MATLAB, the double can be evaluated by:
        //
        //      ellipticPi(n,(sin(a*pi/180)^2)
        //
        //    The value is computed using Carlson elliptic integrals:
        //
        //      k = sin ( a * Math.PI / 180 )
        //      Pi(n,k) = RF ( 0, 1 - k^2, 1 ) + 1/3 n RJ ( 0, 1 - k^2, 1, 1 - n )
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
        //    Input, double N, A, the arguments.
        //
        //    Output, double ELLIPTIC_PIA, the function value.
        //
    {
        double errtol;
        int ierr = 0;
        double k;
        double p;
            
        double value = 0;
        double x;
        double y;
        double z;

        k = Math.Sin(a * Math.PI / 180.0);
        x = 0.0;
        y = (1.0 - k) * (1.0 + k);
        z = 1.0;
        p = 1.0 - n;
        errtol = 1.0E-03;

        value = Integral.rf(x, y, z, errtol, ref ierr)
                + n * Integral.rj(x, y, z, p, errtol, ref ierr) / 3.0;

        return value;
    }

    public static void values(ref int n_data, ref double n, ref double a, ref double pia )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPTIC_PIA_VALUES returns values of the complete elliptic integral Pi(A).
        //
        //  Discussion:
        //
        //    This is one form of what is sometimes called the complete elliptic
        //    integral of the third kind.
        //
        //    The function is defined by the formula:
        //
        //      Pi(N,A) = integral ( 0 <= T <= PI/2 )
        //        dT / (1 - N sin^2(T) ) sqrt ( 1 - sin^2(A) * sin ( T )^2 )
        //
        //    In MATLAB, the function can be evaluated by:
        //
        //      ellipticPi(n,(sin(A*pi/180))^2)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 May 2018
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
        //    Output, double &N, &A, the arguments.
        //
        //    Output, double &PIA, the value of the function.
        //
    {
        const int N_MAX = 20;

        double[] a_vec =
            {
                30.00000000000000,
                45.00000000000000,
                60.00000000000000,
                77.07903361841643,
                30.00000000000000,
                45.00000000000000,
                60.00000000000000,
                77.07903361841643,
                30.00000000000000,
                45.00000000000000,
                60.00000000000000,
                77.07903361841643,
                30.00000000000000,
                45.00000000000000,
                60.00000000000000,
                77.07903361841643,
                30.00000000000000,
                45.00000000000000,
                60.00000000000000,
                77.07903361841643
            }
            ;

        double[] n_vec =
            {
                -10.0,
                -10.0,
                -10.0,
                -10.0,
                -3.0,
                -3.0,
                -3.0,
                -3.0,
                -1.0,
                -1.0,
                -1.0,
                -1.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.5,
                0.5,
                0.5,
                0.5
            }
            ;

        double[] pia_vec =
            {
                0.4892245275965397,
                0.5106765677902629,
                0.5460409271920561,
                0.6237325893535237,
                0.823045542660675,
                0.8760028274011437,
                0.9660073560143946,
                1.171952391481798,
                1.177446843000566,
                1.273127366749682,
                1.440034318657551,
                1.836472172302591,
                1.685750354812596,
                1.854074677301372,
                2.156515647499643,
                2.908337248444552,
                2.413671504201195,
                2.701287762095351,
                3.234773471249465,
                4.633308147279891
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
            a = 0.0;
            n = 0.0;
            pia = 0.0;
        }
        else
        {
            a = a_vec[n_data - 1];
            n = n_vec[n_data - 1];
            pia = pia_vec[n_data - 1];
        }
    }
}