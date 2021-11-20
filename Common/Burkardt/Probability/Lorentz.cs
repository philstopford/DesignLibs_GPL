using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Probability;

public static class Lorentz
{
    public static double lorentz_cdf(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LORENTZ_CDF evaluates the Lorentz CDF.
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
        //    Input, double X, the argument of the CDF.
        //
        //    Output, double LORENTZ_CDF, the value of the CDF.
        //
    {
            

        double cdf = 0.5 + Math.Atan(x) / Math.PI;

        return cdf;
    }

    public static double lorentz_cdf_inv(double cdf)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LORENTZ_CDF_INV inverts the Lorentz CDF.
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
        //    Output, double LORENTZ_CDF_INV, the corresponding argument.
        //
    {
            

        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("LORENTZ_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
            default:
            {
                double x = Math.Tan(Math.PI * (cdf - 0.5));

                return x;
            }
        }
    }

    public static double lorentz_mean()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LORENTZ_MEAN returns the mean of the Lorentz PDF.
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
        //    Output, double LORENTZ_MEAN, the mean of the PDF.
        //
    {
        const double mean = 0.0;

        return mean;
    }

    public static double lorentz_pdf(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LORENTZ_PDF evaluates the Lorentz PDF.
        //
        //  Discussion:
        //
        //    PDF(X) = 1 / ( PI * ( 1 + X^2 ) )
        //
        //    The chief interest of the Lorentz PDF is that it is easily
        //    inverted, and can be used to dominate other PDF's in an
        //    acceptance/rejection method.
        //
        //    LORENTZ_PDF(X) = CAUCHY_PDF(0,1;X)
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
        //
        //    Output, double LORENTZ_PDF, the value of the PDF.
        //
    {
            

        double pdf = 1.0 / (Math.PI * (1.0 + x * x));

        return pdf;
    }

    public static double lorentz_sample(ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LORENTZ_SAMPLE samples the Lorentz PDF.
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
        //    Output, double LORENTZ_SAMPLE, a sample of the PDF.
        //
    {
        double cdf = UniformRNG.r8_uniform_01(ref seed);

        double x = lorentz_cdf_inv(cdf);

        return x;
    }

    public static double lorentz_variance()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LORENTZ_VARIANCE returns the variance of the Lorentz PDF.
        //
        //  Discussion:
        //
        //    The variance of the Lorentz PDF is not well defined.  This routine
        //    is made available for completeness only, and simply returns
        //    a "very large" number.
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
        //    Output, double LORENTZ_VARIANCE, the mean of the PDF.
        //
    {
        return typeMethods.r8_huge();
    }
}