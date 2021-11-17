using System;

namespace Burkardt.Elliptic;

public static class Jacobi_sn
{
    public static double evaluate(double u, double m)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    JACOBI_SN evaluates the Jacobi elliptic function SN(U,M).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 June 2018
        //
        //  Author:
        //
        //    Original ALGOL version by Roland Bulirsch.
        //    C++ version by John Burkardt
        //
        //  Reference:
        //
        //    Roland Bulirsch,
        //    Numerical calculation of elliptic integrals and elliptic functions,
        //    Numerische Mathematik,
        //    Volume 7, Number 1, 1965, pages 78-90.
        //
        //  Parameters:
        //
        //    Input, double U, M, the arguments.
        //
        //    Output, double JACOBI_SN, the function value.
        //
    {
        double cn = 0;
        double dn = 0;
        double sn = 0;

        Jacobi.sncndn(u, m, ref sn, ref cn, ref dn);

        return sn;
    }

    public static void values(ref int n_data, ref double u, ref double a, ref double k,
            ref double m, ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    jacobi_sn_values returns values of the Jacobi elliptic function SN(U,M).
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      JacobiSN[ u, m ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 November 2020
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
        //  Input:
        //
        //    int &N_DATA.  The user sets N_DATA to 0 before the first call.  
        //
        //  Output:
        //
        //    int &N_DATA.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    double &U the argument of the function.
        //
        //    double &A, &K, &M, the parameters of the function.
        //
        //    double &FX, the value of the function.
        //
    {
        const int N_MAX = 20;

        double[] m_vec =
            {
                0.0E+00,
                0.0E+00,
                0.0E+00,
                0.0E+00,
                0.0E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00
            }
            ;

        double[] fx_vec =
            {
                0.9983341664682815E-01,
                0.1986693307950612E+00,
                0.4794255386042030E+00,
                0.8414709848078965E+00,
                0.9092974268256817E+00,
                0.9975068547462484E-01,
                0.1980217429819704E+00,
                0.4707504736556573E+00,
                0.8030018248956439E+00,
                0.9946623253580177E+00,
                0.9966799462495582E-01,
                0.1973753202249040E+00,
                0.4621171572600098E+00,
                0.7615941559557649E+00,
                0.9640275800758169E+00,
                0.9993292997390670E+00,
                -0.1973753202249040E+00,
                -0.4621171572600098E+00,
                -0.7615941559557649E+00,
                -0.9640275800758169E+00
            }
            ;

        double[] u_vec =
            {
                0.1E+00,
                0.2E+00,
                0.5E+00,
                1.0E+00,
                2.0E+00,
                0.1E+00,
                0.2E+00,
                0.5E+00,
                1.0E+00,
                2.0E+00,
                0.1E+00,
                0.2E+00,
                0.5E+00,
                1.0E+00,
                2.0E+00,
                4.0E+00,
                -0.2E+00,
                -0.5E+00,
                -1.0E+00,
                -2.0E+00
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
            k = 0.0;
            m = 0.0;
            u = 0.0;
            fx = 0.0;
        }
        else
        {
            m = m_vec[n_data - 1];
            k = Math.Sqrt(m);
            a = Math.Asin(k);
            u = u_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }
    }
}