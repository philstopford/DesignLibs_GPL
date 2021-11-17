using System;

namespace Burkardt.Elliptic;

public static class Jacobi_cn
{
    public static double evaluate(double u, double m)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    JACOBI_CN evaluates the Jacobi elliptic function CN(U,M).
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
        //    Output, double JACOBI_CN, the function value.
        //
    {
        double cn = 0;
        double dn = 0;
        double sn = 0;

        Jacobi.sncndn(u, m, ref sn, ref cn, ref dn);

        return cn;
    }

    public static void values(ref int n_data, ref double u, ref double a, ref double k,
            ref double m, ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    jacobi_cn_values returns values of the Jacobi elliptic function CN(U,M).
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      JacobiCN[ u, m ]
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
                0.9950041652780258E+00,
                0.9800665778412416E+00,
                0.8775825618903727E+00,
                0.5403023058681397E+00,
                -0.4161468365471424E+00,
                0.9950124626090582E+00,
                0.9801976276784098E+00,
                0.8822663948904403E+00,
                0.5959765676721407E+00,
                -0.1031836155277618E+00,
                0.9950207489532265E+00,
                0.9803279976447253E+00,
                0.8868188839700739E+00,
                0.6480542736638854E+00,
                0.2658022288340797E+00,
                0.3661899347368653E-01,
                0.9803279976447253E+00,
                0.8868188839700739E+00,
                0.6480542736638854E+00,
                0.2658022288340797E+00
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