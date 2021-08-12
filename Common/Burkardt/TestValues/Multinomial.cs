namespace Burkardt.TestValues
{
    public static class Multinomial
    {

        public static void multinomial_pdf_sizes(ref int n_data, ref int m)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MULTINOMIAL_PDF_SIZES returns sizes of some multinomial PDF data.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    31 July 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input/output, integer &N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, integer &M, the size of the given problem.
            //
        {
            int N_MAX = 10;

            int[] m_vec =
            {
                2, 2, 2, 3, 5,
                5, 5, 5, 5, 5
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                m = 0;
            }
            else
            {
                m = m_vec[n_data - 1];
            }
        }

        public static void multinomial_pdf_values(ref int n_data, int m, ref int n, ref double[] p, ref int[] x,
                ref double pdf)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MULTINOMIAL_PDF_VALUES returns some values of the multinomial PDF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    31 July 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Input, int M, the number of outcomes.
            //
            //    Output, ref int N, the number of trials.
            //
            //    Output, double P[M], the probability of each outcome on one trial.
            //
            //    Output, int X[M], the number of times each outcome occurred in
            //    N trials.
            //
            //    Output, ref double PDF, the probability of X.
            //
        {
            int N_MAX = 10;

            int i;
            int[] n_vec =
            {
                3, 4, 3, 3, 3,
                3, 3, 3, 3, 3
            };

            int[] offset =
            {
                0, 2, 4, 6, 9,
                14, 19, 24, 29, 34
            };

            double[] p_vec =
            {
                0.7, 0.3,
                0.7, 0.3,
                0.5, 0.5,
                0.6, 0.0, 0.4,
                0.6, 0.1, 0.1, 0.1, 0.1,
                0.6, 0.1, 0.1, 0.1, 0.1,
                0.6, 0.1, 0.1, 0.1, 0.1,
                0.6, 0.1, 0.1, 0.1, 0.1,
                0.6, 0.1, 0.1, 0.1, 0.1,
                0.6, 0.1, 0.1, 0.1, 0.1
            };

            double[] pdf_vec =
            {
                0.441,
                0.2646,
                0.375,
                0.0,
                0.216,
                0.108,
                0.018,
                0.036,
                0.001,
                0.006
            };

            int[] x_vec =
            {
                2, 1,
                2, 2,
                2, 1,
                1, 1, 1,
                3, 0, 0, 0, 0,
                2, 1, 0, 0, 0,
                1, 0, 2, 0, 0,
                1, 0, 0, 1, 1,
                0, 0, 0, 3, 0,
                0, 1, 1, 1, 0
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
                p = null;
                x = null;
                pdf = 0.0;
            }
            else
            {
                n = n_vec[n_data - 1];
                for (i = 0; i < m; i++)
                {
                    p[i] = p_vec[i + offset[n_data - 1]];
                }

                for (i = 0; i < m; i++)
                {
                    x[i] = x_vec[i + offset[n_data - 1]];
                }

                pdf = pdf_vec[n_data - 1];
            }
        }

    }
}