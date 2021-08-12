using System;

namespace Burkardt.TestValues
{
    public static class Jacobi
    {

        public static void jacobi_cn_values(ref int n_data, ref double u, ref double a, ref double k,
                ref double m, ref double fx)

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
            //    ref int N_DATA.  The user sets N_DATA to 0 before the first call.  
            //
            //  Output:
            //
            //    ref int N_DATA.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    ref double U the argument of the function.
            //
            //    ref double A, &K, &M, the parameters of the function.
            //
            //    ref double FX, the value of the function.
            //
        {
            int N_MAX = 20;

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
            };

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
            };

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
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

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

        public static void jacobi_dn_values(ref int n_data, ref double u, ref double a, ref double k,
                ref double m, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    jacobi_dn_values returns values of the Jacobi elliptic function DN(U,M).
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      JacobiDN[ u, m ]
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
            //    ref int N_DATA.  The user sets N_DATA to 0 before the first call.  
            //
            //  Output:
            //
            //    ref int N_DATA.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    ref double U the argument of the function.
            //
            //    ref double A, &K, &M, the parameters of the function.
            //
            //    ref double FX, the value of the function.
            //
        {
            int N_MAX = 20;

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
            };

            double[] fx_vec =
            {
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.9975093485144243E+00,
                0.9901483195224800E+00,
                0.9429724257773857E+00,
                0.8231610016315963E+00,
                0.7108610477840873E+00,
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
            };

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
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

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

        public static void jacobi_poly_values(ref int n_data, ref int n, ref double a, ref double b, ref double x,
                ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    JACOBI_POLY_VALUES returns some values of the Jacobi polynomial.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      JacobiP[ n, a, b, x ]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 May 2018
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
            //    Output, ref int N, the degree of the polynomial.
            //
            //    Output, ref double A, &B, parameters of the function.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 26;

            double[] a_vec =
            {
                0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 2.0,
                3.0, 4.0, 5.0, 0.0,
                0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0,
                0.0, 0.0
            };

            double[] b_vec =
            {
                1.0, 1.0, 1.0, 1.0,
                1.0, 1.0, 1.0, 1.0,
                1.0, 1.0, 1.0, 2.0,
                3.0, 4.0, 5.0, 1.0,
                1.0, 1.0, 1.0, 1.0,
                1.0, 1.0, 1.0, 1.0,
                1.0, 1.0
            };

            double[] fx_vec =
            {
                0.1000000000000000E+01,
                0.2500000000000000E+00,
                -0.3750000000000000E+00,
                -0.4843750000000000E+00,
                -0.1328125000000000E+00,
                0.2753906250000000E+00,
                -0.1640625000000000E+00,
                -0.1174804687500000E+01,
                -0.2361328125000000E+01,
                -0.2616210937500000E+01,
                0.1171875000000000E+00,
                0.4218750000000000E+00,
                0.5048828125000000E+00,
                0.5097656250000000E+00,
                0.4306640625000000E+00,
                -0.6000000000000000E+01,
                0.3862000000000000E-01,
                0.8118400000000000E+00,
                0.3666000000000000E-01,
                -0.4851200000000000E+00,
                -0.3125000000000000E+00,
                0.1891200000000000E+00,
                0.4023400000000000E+00,
                0.1216000000000000E-01,
                -0.4396200000000000E+00,
                0.1000000000000000E+01
            };

            int[] n_vec =
            {
                0, 1, 2, 3,
                4, 5, 5, 5,
                5, 5, 5, 5,
                5, 5, 5, 5,
                5, 5, 5, 5,
                5, 5, 5, 5,
                5, 5
            };

            double[] x_vec =
            {
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                -1.0E+00,
                -0.8E+00,
                -0.6E+00,
                -0.4E+00,
                -0.2E+00,
                0.0E+00,
                0.2E+00,
                0.4E+00,
                0.6E+00,
                0.8E+00,
                1.0E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                n = 0;
                a = 0.0;
                b = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                n = n_vec[n_data - 1];
                a = a_vec[n_data - 1];
                b = b_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void jacobi_sn_values(ref int n_data, ref double u, ref double a, ref double k,
                ref double m, ref double fx)

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
            //    ref int N_DATA.  The user sets N_DATA to 0 before the first call.  
            //
            //  Output:
            //
            //    ref int N_DATA.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    ref double U the argument of the function.
            //
            //    ref double A, &K, &M, the parameters of the function.
            //
            //    ref double FX, the value of the function.
            //
        {
            int N_MAX = 20;

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
            };

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
            };

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
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

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
}