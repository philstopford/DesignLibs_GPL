using System;
using Burkardt.Uniform;

namespace Burkardt.Probability
{
    public static class Exponential
    {
        public static double exponential_01_cdf(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXPONENTIAL_01_CDF evaluates the Exponential 01 CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //
        //    Output, double EXPONENTIAL_01_CDF, the value of the CDF.
        //
        {
            double cdf;

            if (x <= 0.0)
            {
                cdf = 0.0;
            }
            else
            {
                cdf = 1.0 - Math.Exp(-x);
            }

            return cdf;
        }

        public static double exponential_01_cdf_inv(double cdf)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXPONENTIAL_01_CDF_INV inverts the Exponential 01 CDF.
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
        //  Parameters:
        //
        //    Input, double CDF, the value of the CDF.
        //    0.0 <= CDF <= 1.0.
        //
        //    Output, double EXPONENTIAL_CDF_INV, the corresponding argument.
        //
        {
            double x;

            if (cdf < 0.0 || 1.0 < cdf)
            {
                Console.WriteLine(" ");
                Console.WriteLine("EXPONENTIAL_01_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return (1);
            }

            x = -Math.Log(1.0 - cdf);

            return x;
        }

        public static double exponential_01_mean()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXPONENTIAL_01_MEAN returns the mean of the Exponential 01 PDF.
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
        //  Parameters:
        //
        //    Output, double EXPONENTIAL_MEAN, the mean of the PDF.
        //
        {
            double mean = 1.0;

            return mean;
        }

        public static double exponential_01_pdf(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXPONENTIAL_01_PDF evaluates the Exponential 01 PDF.
        //
        //  Discussion:
        //
        //    PDF(X) = EXP ( - X )
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
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //    0.0 <= X
        //
        //    Output, double PDF, the value of the PDF.
        //
        {
            double pdf;

            if (x < 0.0)
            {
                pdf = 0.0;
            }
            else
            {
                pdf = Math.Exp(-x);
            }

            return pdf;
        }

        public static double exponential_01_sample(ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXPONENTIAL_01_SAMPLE samples the Exponential PDF with parameter 1.
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
        //  Parameters:
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double EXPONENTIAL_01_SAMPLE, a sample of the PDF.
        //
        {
            double cdf = UniformRNG.r8_uniform_01(ref seed);

            double x = -Math.Log(1.0 - cdf);

            return x;
        }

        public static double exponential_01_variance()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXPONENTIAL_01_VARIANCE returns the variance of the Exponential 01 PDF.
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
        //  Parameters:
        //
        //    Output, double VARIANCE, the variance of the PDF.
        //
        {
            double variance = 1.0;

            return variance;
        }

        public static double exponential_cdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXPONENTIAL_CDF evaluates the Exponential CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //
        //    Input, double A, B, the parameter of the PDF.
        //    0.0 < B.
        //
        //    Output, double EXPONENTIAL_CDF, the value of the CDF.
        //
        {
            double cdf;

            if (x <= a)
            {
                cdf = 0.0;
            }
            else
            {
                cdf = 1.0 - Math.Exp((a - x) / b);
            }

            return cdf;
        }

        public static double exponential_cdf_inv(double cdf, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXPONENTIAL_CDF_INV inverts the Exponential CDF.
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
        //  Parameters:
        //
        //    Input, double CDF, the value of the CDF.
        //    0.0 <= CDF <= 1.0.
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < B.
        //
        //    Output, double EXPONENTIAL_CDF_INV, the corresponding argument.
        //
        {
            if (cdf < 0.0 || 1.0 < cdf)
            {
                Console.WriteLine(" ");
                Console.WriteLine("EXPONENTIAL_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return (1);
            }

            double x = a - b * Math.Log(1.0 - cdf);

            return x;
        }

        public static void exponential_cdf_values(ref int n_data, ref double lambda, ref double x,
            ref double fx )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXPONENTIAL_CDF_VALUES returns some values of the Exponential CDF.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Needs["Statistics`ContinuousDistributions`"]
        //      dist = ExponentialDistribution [ lambda ]
        //      CDF [ dist, x ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 August 2004
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
        //    Output, double &LAMBDA, the parameter of the distribution.
        //
        //    Output, double &X, the argument of the function.
        //
        //    Output, double &FX, the value of the function.
        //
        {
            int N_MAX = 9;

            double[] fx_vec =
            {
                0.3934693402873666E+00,
                0.6321205588285577E+00,
                0.7768698398515702E+00,
                0.8646647167633873E+00,
                0.8646647167633873E+00,
                0.9816843611112658E+00,
                0.9975212478233336E+00,
                0.9996645373720975E+00,
                0.9999546000702375E+00
            } ;

            double[] lambda_vec =
            {
                0.5000000000000000E+00,
                0.5000000000000000E+00,
                0.5000000000000000E+00,
                0.5000000000000000E+00,
                0.1000000000000000E+01,
                0.2000000000000000E+01,
                0.3000000000000000E+01,
                0.4000000000000000E+01,
                0.5000000000000000E+01
            } ;

            double[] x_vec =
            {
                0.1000000000000000E+01,
                0.2000000000000000E+01,
                0.3000000000000000E+01,
                0.4000000000000000E+01,
                0.2000000000000000E+01,
                0.2000000000000000E+01,
                0.2000000000000000E+01,
                0.2000000000000000E+01,
                0.2000000000000000E+01
            } ;

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                lambda = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                lambda = lambda_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static bool exponential_check(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXPONENTIAL_CHECK checks the parameters of the Exponential CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameter of the PDF.
        //    0.0 < B.
        //
        //    Output, bool EXPONENTIAL_CHECK, is true if the parameters are legal.
        //
        {
            if (b <= 0.0)
            {
                Console.WriteLine(" ");
                Console.WriteLine("EXPONENTIAL_CHECK - Warning!");
                Console.WriteLine("  B <= 0.0");
                return false;
            }

            return true;
        }

        public static double exponential_mean(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXPONENTIAL_MEAN returns the mean of the Exponential PDF.
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
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < B.
        //
        //    Output, double EXPONENTIAL_MEAN, the mean of the PDF.
        //
        {
            double mean = a + b;

            return mean;
        }

        public static double exponential_pdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXPONENTIAL_PDF evaluates the Exponential PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B;X) = ( 1 / B ) * EXP ( ( A - X ) / B )
        //
        //    The time interval between two Poisson events is a random
        //    variable with the Exponential PDF.  The parameter B is the
        //    average interval between events.
        //
        //    In another context, the Exponential PDF is related to
        //    the Boltzmann distribution, which describes the relative
        //    probability of finding a system, which is in thermal equilibrium
        //    at absolute temperature T, in a given state having energy E.
        //    The relative probability is
        //
        //      Boltzmann_Relative_Probability(E,T) = exp ( - E / ( k * T ) ),
        //
        //    where k is the Boltzmann constant,
        //
        //      k = 1.38 * 10**(-23) joules / degree Kelvin
        //
        //    and normalization requires a determination of the possible
        //    energy states of the system.
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
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //    A <= X
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < B.
        //
        //    Output, double EXPONENTIAL_PDF, the value of the PDF.
        //
        {
            double pdf;

            if (x < a)
            {
                pdf = 0.0;
            }
            else
            {
                pdf = (1.0 / b) * Math.Exp((a - x) / b);
            }

            return pdf;
        }

        public static double exponential_sample(double a, double b, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXPONENTIAL_SAMPLE samples the Exponential PDF.
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
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < B.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double EXPONENTIAL_SAMPLE, a sample of the PDF.
        //
        {
            double cdf = UniformRNG.r8_uniform_01(ref seed);

            double x = exponential_cdf_inv(cdf, a, b);

            return x;
        }

        public static double exponential_variance(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXPONENTIAL_VARIANCE returns the variance of the Exponential PDF.
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
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < B.
        //
        //    Output, double EXPONENTIAL_VARIANCE, the variance of the PDF.
        //
        {
            double variance = b * b;

            return variance;
        }
    }
}