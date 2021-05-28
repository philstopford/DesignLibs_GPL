using System;
using Burkardt.Types;

namespace Burkardt.Probability
{
    public static class FProb
    {
        public static double f_cdf(double x, int m, int n)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F_CDF evaluates the F central CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 October 2004
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
        //  Parameters:
        //
        //    Input, double X, the argument of the CDF.
        //
        //    Input, int M, N, the parameters of the PDF.
        //    1 <= M,
        //    1 <= N.
        //
        //    Output, double F_CDF, the value of the CDF.
        //
        {
            double cdf;

            if (x <= 0.0)
            {
                cdf = 0.0;
            }
            else
            {
                double arg1 = 0.5 * (double) (n);
                double arg2 = 0.5 * (double) (m);

                double arg3 = (double) (n) / ((double) (n) + (double) (m) * x);

                cdf = 1.0 - Beta.beta_inc(arg1, arg2, arg3);
            }

            return cdf;
        }

        public static void f_cdf_values(ref int n_data, ref int a, ref int b, ref double x, ref double fx )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F_CDF_VALUES returns some values of the F CDF test function.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Needs["Statistics`ContinuousDistributions`"]
        //      dist = FRatioDistribution [ dfn, dfd ]
        //      CDF [ dist, x ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 August 2004
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
        //    Output, int &A, int &B, the parameters of the function.
        //
        //    Output, double &X, the argument of the function.
        //
        //    Output, double &FX, the value of the function.
        //
        {
            int N_MAX = 20;

            int[] a_vec =
            {
                1,
                1,
                5,
                1,
                2,
                4,
                1,
                6,
                8,
                1,
                3,
                6,
                1,
                1,
                1,
                1,
                2,
                3,
                4,
                5
            };

            int[] b_vec =
            {
                1,
                5,
                1,
                5,
                10,
                20,
                5,
                6,
                16,
                5,
                10,
                12,
                5,
                5,
                5,
                5,
                5,
                5,
                5,
                5
            };

            double[] fx_vec =
            {
                0.5000000000000000E+00,
                0.4999714850534485E+00,
                0.4996034370170990E+00,
                0.7496993658293228E+00,
                0.7504656462757382E+00,
                0.7514156325324275E+00,
                0.8999867031372156E+00,
                0.8997127554259699E+00,
                0.9002845660853669E+00,
                0.9500248817817622E+00,
                0.9500574946122442E+00,
                0.9501926400000000E+00,
                0.9750133887312993E+00,
                0.9900022327445249E+00,
                0.9949977837872073E+00,
                0.9989999621122122E+00,
                0.5687988496283079E+00,
                0.5351452100063650E+00,
                0.5143428032407864E+00,
                0.5000000000000000E+00
            };

            double[] x_vec =
            {
                1.00E+00,
                0.528E+00,
                1.89E+00,
                1.69E+00,
                1.60E+00,
                1.47E+00,
                4.06E+00,
                3.05E+00,
                2.09E+00,
                6.61E+00,
                3.71E+00,
                3.00E+00,
                10.01E+00,
                16.26E+00,
                22.78E+00,
                47.18E+00,
                1.00E+00,
                1.00E+00,
                1.00E+00,
                1.00E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                a = 0;
                b = 0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                a = a_vec[n_data - 1];
                b = b_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static bool f_check(int m, int n)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F_CHECK checks the parameters of the F PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the parameters of the PDF.
        //    1 <= M,
        //    1 <= N.
        //
        //    Output, bool F_CHECK, is true if the parameters are legal.
        //
        {
            if (m <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("F_CHECK - Warning!");
                Console.WriteLine("  M <= 0.");
                return false;
            }

            if (n <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("F_CHECK - Warning!");
                Console.WriteLine("  N <= 0.");
                return false;
            }

            return true;
        }

        public static double f_mean(int m, int n)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F_MEAN returns the mean of the F central PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the parameters of the PDF.
        //    1 <= M,
        //    1 <= N.
        //    Note, however, that the mean is not defined unless 3 <= N.
        //
        //    Output, double F_MEAN, the mean of the PDF.
        //
        {
            if (n < 3)
            {
                Console.WriteLine(" ");
                Console.WriteLine("F_MEAN - Fatal error!");
                Console.WriteLine("  The mean is not defined for N < 3.");
                return (1);
            }

            double mean = (double) (n) / (double) (n - 2);

            return mean;
        }

        public static double f_pdf(double x, int m, int n)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F_PDF evaluates the F central PDF.
        //
        //  Discussion:
        //
        //    PDF(M,N;X) = M^(M/2) * X^((M-2)/2)
        //      / ( Beta(M/2,N/2) * N^(M/2) * ( 1 + (M/N) * X )^((M+N)/2)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //    0.0 <= X
        //
        //    Input, int M, N, the parameters of the PDF.
        //    1 <= M,
        //    1 <= N.
        //
        //    Output, double F_PDF, the value of the PDF.
        //
        {
            double pdf;

            if (x < 0.0)
            {
                pdf = 0.0;
            }
            else
            {
                double a = (double) (m);
                double b = (double) (n);

                double top = Math.Sqrt(Math.Pow(a, m) * Math.Pow(b, n) * Math.Pow(x, m - 2));
                double bot1 = typeMethods.r8_beta(a / 2.0, b / 2.0);
                double bot2 = Math.Sqrt(Math.Pow(b + a * x, m + n));

                pdf = top / (bot1 * bot2);
            }

            return pdf;
        }

        public static double f_sample(int m, int n, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F_SAMPLE samples the F central PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the parameters of the PDF.
        //    1 <= M,
        //    1 <= N.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double F_SAMPLE, a sample of the PDF.
        //
        {
            double a = (double) (m);
            double xm = Chi.chi_square_sample(a, ref seed);

            a = (double) (n);
            double xn = Chi.chi_square_sample(a, ref seed);

            double x = (double) (n) * xm / ((double) (m) * xn);

            return x;
        }

        public static double f_variance(int m, int n)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F_VARIANCE returns the variance of the F central PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the parameters of the PDF.
        //    1 <= M,
        //    1 <= N.
        //    Note, however, that the variance is not defined unless 5 <= N.
        //
        //    Output, double F_VARIANCE, the variance of the PDF.
        //
        {
            if (n < 5)
            {
                Console.WriteLine(" ");
                Console.WriteLine("F_VARIANCE - Fatal error!");
                Console.WriteLine("  The variance is not defined for N < 5.");
                return (1);
            }

            double variance = (double) (2 * n * n * (m + n - 2)) /
                              (double) (m * (n - 2) * (n - 2) * (n - 4));

            return variance;
        }

        public static void f_noncentral_cdf_values(ref int n_data, ref int n1, ref int n2, ref double lambda,
        ref double x, ref double fx )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F_NONCENTRAL_CDF_VALUES returns some values of the F CDF test function.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Needs["Statistics`ContinuousDistributions`"]
        //      dist = NoncentralFRatioDistribution [ n1, n2, lambda ]
        //      CDF [ dist, x ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 August 2004
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
        //    Output, int &N1, int &N2, the numerator and denominator
        //    degrees of freedom.
        //
        //    Output, double &LAMBDA, the noncentrality parameter.
        //
        //    Output, double &X, the argument of the function.
        //
        //    Output, double &FX, the value of the function.
        //
        {
            int N_MAX = 22;

            double[] fx_vec =
            {
                0.5000000000000000E+00,
                0.6367825323508774E+00,
                0.5840916116305482E+00,
                0.3234431872392788E+00,
                0.4501187879813550E+00,
                0.6078881441188312E+00,
                0.7059275551414605E+00,
                0.7721782003263727E+00,
                0.8191049017635072E+00,
                0.3170348430749965E+00,
                0.4327218008454471E+00,
                0.4502696915707327E+00,
                0.4261881186594096E+00,
                0.6753687206341544E+00,
                0.4229108778879005E+00,
                0.6927667261228938E+00,
                0.3632174676491226E+00,
                0.4210054012695865E+00,
                0.4266672258818927E+00,
                0.4464016600524644E+00,
                0.8445888579504827E+00,
                0.4339300273343604E+00
            };

            double[] lambda_vec =
            {
                0.00E+00,
                0.00E+00,
                0.25E+00,
                1.00E+00,
                1.00E+00,
                1.00E+00,
                1.00E+00,
                1.00E+00,
                1.00E+00,
                2.00E+00,
                1.00E+00,
                1.00E+00,
                1.00E+00,
                2.00E+00,
                1.00E+00,
                1.00E+00,
                0.00E+00,
                1.00E+00,
                1.00E+00,
                1.00E+00,
                1.00E+00,
                1.00E+00
            };

            int[] n1_vec =
            {
                1, 1, 1, 1,
                1, 1, 1, 1,
                1, 1, 2, 2,
                3, 3, 4, 4,
                5, 5, 6, 6,
                8, 16
            };

            int[] n2_vec =
            {
                1, 5, 5, 5,
                5, 5, 5, 5,
                5, 5, 5, 10,
                5, 5, 5, 5,
                1, 5, 6, 12,
                16, 8
            };

            double[] x_vec =
            {
                1.00E+00,
                1.00E+00,
                1.00E+00,
                0.50E+00,
                1.00E+00,
                2.00E+00,
                3.00E+00,
                4.00E+00,
                5.00E+00,
                1.00E+00,
                1.00E+00,
                1.00E+00,
                1.00E+00,
                1.00E+00,
                1.00E+00,
                2.00E+00,
                1.00E+00,
                1.00E+00,
                1.00E+00,
                1.00E+00,
                2.00E+00,
                2.00E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                n1 = 0;
                n2 = 0;
                lambda = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                n1 = n1_vec[n_data - 1];
                n2 = n2_vec[n_data - 1];
                lambda = lambda_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static bool f_noncentral_check(double a, int m, int n)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F_NONCENTRAL_CHECK checks the parameters of the F noncentral PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, a parameter of the PDF.
        //
        //    Input, int M, N, the parameters of the PDF.
        //    1 <= M,
        //    1 <= N.
        //
        //    Output, bool F_NONCENTRAL_CHECK, is true if the parameters are legal.
        //
        {
            if (a <= 0.0)
            {
                Console.WriteLine(" ");
                Console.WriteLine("F_NONCENTRAL_CHECK - Warning!");
                Console.WriteLine("  A <= 0.");
                return false;
            }

            if (m <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("F_NONCENTRAL_CHECK - Warning!");
                Console.WriteLine("  M <= 0.");
                return false;
            }

            if (n <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("F_NONCENTRAL_CHECK - Warning!");
                Console.WriteLine("  N <= 0.");
                return false;
            }

            return true;
        }

        public static double f_noncentral_mean(double a, int m, int n)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F_NONCENTRAL_MEAN returns the mean of the F noncentral PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 February 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, a parameter of the PDF.
        //
        //    Input, int M, N, parameters of the PDF.
        //    1 <= M,
        //    1 <= N.
        //    Note, however, that the mean is not defined unless 3 <= N.
        //
        //    Output, double F_NONCENTAL_MEAN, the mean of the PDF.
        //
        {
            if (n < 3)
            {
                Console.WriteLine(" ");
                Console.WriteLine("F_NONCENTRAL_MEAN - Fatal error!");
                Console.WriteLine("  The mean is not defined for N < 3.");
                return (1);
            }

            double mean = ((double) (m) + a) * (double) (n)
                          / ((double) (m) * (double) (n - 2));

            return mean;
        }

        public static double f_noncentral_variance(double a, int m, int n)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F_NONCENTRAL_VARIANCE returns the variance of the F noncentral PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, a parameter of the PDF.
        //
        //    Input, int M, N, parameters of the PDF.
        //    1 <= M,
        //    1 <= N.
        //    Note, however, that the variance is not defined unless 5 <= N.
        //
        //    Output, double F_NONCENTRAL_VARIANCE, the variance of the PDF.
        //
        {
            if (n < 5)
            {
                Console.WriteLine(" ");
                Console.WriteLine("F_NONCENTRAL_VARIANCE - Fatal error!");
                Console.WriteLine("  The variance is not defined for N < 5.");
                return (1);
            }

            double mr = (double) (m);
            double nr = (double) (n);

            double variance = ((mr + a) * (mr + a) + 2.0 * (mr + a) * nr * nr)
                              / ((nr - 2.0) * (nr - 4.0) * mr * mr)
                              - (mr + a) * (mr + a) * nr * nr
                              / ((nr - 2.0) * (nr - 2.0) * mr * mr);

            return variance;
        }
    }
}