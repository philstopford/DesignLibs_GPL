using System;
using Burkardt.Uniform;

namespace Burkardt.Probability
{
    public static class Rayleigh
    {
        public static double rayleigh_cdf(double x, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RAYLEIGH_CDF evaluates the Rayleigh CDF.
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
        //    Input, double X, the argument of the CDF.
        //    0.0 <= X.
        //
        //    Input, double A, the parameter of the PDF.
        //    0.0 < A.
        //
        //    Output, double RAYLEIGH_CDF, the value of the CDF.
        //
        {
            double cdf;

            if (x < 0.0)
            {
                cdf = 0.0;
            }
            else
            {
                cdf = 1.0 - Math.Exp(-x * x / (2.0 * a * a));
            }

            return cdf;
        }

        public static double rayleigh_cdf_inv(double cdf, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RAYLEIGH_CDF_INV inverts the Rayleigh CDF.
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
        //    Input, double CDF, the value of the CDF.
        //    0.0 <= CDF <= 1.0.
        //
        //    Input, double A, the parameter of the PDF.
        //    0.0 < A.
        //
        //    Output, double RAYLEIGH_CDF_INV, the corresponding argument.
        //
        {
            double x;

            if (cdf < 0.0 || 1.0 < cdf)
            {
                Console.WriteLine(" ");
                Console.WriteLine("RAYLEIGH_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return (1);
            }

            x = Math.Sqrt(-2.0 * a * a * Math.Log(1.0 - cdf));

            return x;
        }

        public static void rayleigh_cdf_values(ref int n_data, ref double sigma, ref double x, ref double fx )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RAYLEIGH_CDF_VALUES returns some values of the Rayleigh CDF.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Needs["Statistics`ContinuousDistributions`"]
        //      dist = RayleighDistribution [ sigma ]
        //      CDF [ dist, x ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 September 2004
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
        //    Output, double &SIGMA, the shape parameter of the distribution.
        //
        //    Output, double &X, the argument of the function.
        //
        //    Output, double &FX, the value of the function.
        //
        {
            int N_MAX = 9;

            double[] fx_vec =
            {
                0.8646647167633873E+00,
                0.9996645373720975E+00,
                0.9999999847700203E+00,
                0.999999999999987E+00,
                0.8646647167633873E+00,
                0.3934693402873666E+00,
                0.1992625970831920E+00,
                0.1175030974154046E+00,
                0.7688365361336422E-01
            }
            ;

            double[] sigma_vec =
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
            }
            ;

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
            }
            ;

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                sigma = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                sigma = sigma_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static bool rayleigh_check(double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RAYLEIGH_CHECK checks the parameter of the Rayleigh PDF.
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
        //    Input, double A, the parameter of the PDF.
        //    0.0 < A.
        //
        //    Output, bool RAYLEIGH_CHECK, is true if the parameter is legal.
        //
        {
            if (a <= 0.0)
            {
                Console.WriteLine(" ");
                Console.WriteLine("RAYLEIGH_CHECK - Warning!");
                Console.WriteLine("  A <= 0.");
                return false;
            }

            return true;
        }

        public static double rayleigh_mean(double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RAYLEIGH_MEAN returns the mean of the Rayleigh PDF.
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
        //    Input, double A, the parameter of the PDF.
        //    0.0 < A.
        //
        //    Output, double RAYLEIGH_MEAN, the mean of the PDF.
        //
        {
            double mean;
            const double r8_pi = 3.14159265358979323;

            mean = a * Math.Sqrt(0.5 * r8_pi);

            return mean;
        }

        public static double rayleigh_pdf(double x, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RAYLEIGH_PDF evaluates the Rayleigh PDF.
        //
        //  Discussion:
        //
        //    PDF(A;X) = ( X / A^2 ) * EXP ( - X^2 / ( 2 * A^2 ) )
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
        //    Input, double A, the parameter of the PDF.
        //    0 < A.
        //
        //    Output, double RAYLEIGH_PDF, the value of the PDF.
        //
        {
            double pdf;

            if (x < 0.0)
            {
                pdf = 0.0;
            }
            else
            {
                pdf = (x / (a * a)) * Math.Exp(-x * x / (2.0 * a * a));
            }

            return pdf;
        }

        public static double rayleigh_sample(double a, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RAYLEIGH_SAMPLE samples the Rayleigh PDF.
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
        //    Input, double A, the parameter of the PDF.
        //    0.0 < A.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double RAYLEIGH_SAMPLE, a sample of the PDF.
        //
        {
            double cdf = UniformRNG.r8_uniform_01(ref seed);

            double x = rayleigh_cdf_inv(cdf, a);

            return x;
        }

        public static double rayleigh_variance(double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RAYLEIGH_VARIANCE returns the variance of the Rayleigh PDF.
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
        //    Input, double A, the parameters of the PDF.
        //    0.0 < A.
        //
        //    Output, double RAYLEIGH_VARIANCE, the variance of the PDF.
        //
        {
            const double r8_pi = 3.14159265358979323;

            double variance = 2.0 * a * a * (1.0 - 0.25 * r8_pi);

            return variance;
        }
    }
}