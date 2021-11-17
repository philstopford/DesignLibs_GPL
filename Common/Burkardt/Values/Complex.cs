using System.Numerics;

namespace Burkardt.Values;

public class Cmplex
{
    public static void c8_log_values(ref int n_data, ref Complex z, ref Complex fz )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8_LOG_VALUES: the logarithm of a complex value.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 January 2019
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    David Collens,
        //    Algorithm 243: Logarithm of a Complex Number,
        //    Communications of the Association for Computing Machinery,
        //    Volume 7, Number 11, November 1964, page 660.
        //
        //  Parameters:
        //
        //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, complex <double> &Z, the argument of the function.
        //
        //    Output, complex <double> &FZ, the value of the function.
        //
    {
        const int N_MAX = 12;

        Complex[] fz_vec =
            {
                new(1.039720770839918, -2.356194490192345),
                new(0.804718956217050, +2.677945044588987),
                new(0.346573590279973, -2.356194490192345),
                new(0.000000000000000, +3.141592653589793),
                new(0.693147180559945, -1.570796326794897),
                new(0.000000000000000, -1.570796326794897),
                new(0.000000000000000, +1.570796326794897),
                new(0.693147180559945, +1.570796326794897),
                new(0.346573590279973, -0.785398163397448),
                new(0.000000000000000, +0.000000000000000),
                new(1.039720770839918, -0.785398163397448),
                new(0.804718956217050, +0.463647609000806)
            }
            ;
        Complex[] z_vec =
            {
                new(-2.0, -2.0),
                new(-2.0, +1.0),
                new(-1.0, -1.0),
                new(-1.0, +0.0),
                new(0.0, -2.0),
                new(0.0, -1.0),
                new(0.0, +1.0),
                new(0.0, +2.0),
                new(1.0, -1.0),
                new(1.0, +0.0),
                new(2.0, -2.0),
                new(2.0, +1.0)
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
            z = new Complex(0.0, 0.0);
            fz = new Complex(0.0, 0.0);
        }
        else
        {
            z = z_vec[n_data - 1];
            fz = fz_vec[n_data - 1];
        }
    }

}