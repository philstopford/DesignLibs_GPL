using System;

namespace Burkardt.Probability;

public static class Dipole
{
    public static double dipole_cdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIPOLE_CDF evaluates the Dipole CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the CDF.
        //
        //    Input, double A, B, the parameters of the PDF.
        //    A is arbitrary, but represents an angle, so only 0 <= A <= 2 * PI
        //    is interesting, and -1.0 <= B <= 1.0.
        //
        //    Output, double DIPOLE_CDF, the value of the CDF.
        //
    {
        double cdf = 0.5 + 1.0 / Math.PI * Math.Atan(x) + b * b * (x * Math.Cos(2.0 * a)
                                                                   - Math.Sin(2.0 * a)) / (Math.PI * (1.0 + x * x));

        return cdf;
    }

    public static double dipole_cdf_inv(double cdf, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIPOLE_CDF_INV inverts the Dipole CDF.
        //
        //  Discussion:
        //
        //    A simple bisection method is used.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double CDF, the value of the CDF.
        //
        //    Input, double A, B, the parameters of the PDF.
        //    -1.0 <= B <= 1.0.
        //
        //    Output, double DIPOLE_CDF_INV, the corresponding argument of the CDF.
        //
    {
        double cdf1;
        double cdf2;
        const int it_max = 100;
        const double r8_huge = 1.0E+30;
        const double tol = 0.0001;
        double x;

        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("DIPOLE_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
            case 0.0:
                x = -r8_huge;
                return x;
            case 1.0:
                x = r8_huge;
                return x;
        }

        //
        //  Seek X1 < X < X2.
        //
        double x1 = -1.0;

        for (;;)
        {
            cdf1 = dipole_cdf(x1, a, b);

            if (cdf1 <= cdf)
            {
                break;
            }

            x1 = 2.0 * x1;
        }

        double x2 = 1.0;

        for (;;)
        {
            cdf2 = dipole_cdf(x2, a, b);

            if (cdf <= cdf2)
            {
                break;
            }

            x2 = 2.0 * x2;
        }

        //
        //  Now use bisection.
        //
        int it = 0;

        for (;;)
        {
            it += 1;

            double x3 = 0.5 * (x1 + x2);
            double cdf3 = dipole_cdf(x3, a, b);

            if (Math.Abs(cdf3 - cdf) < tol)
            {
                x = x3;
                break;
            }

            if (it_max < it)
            {
                Console.WriteLine(" ");
                Console.WriteLine("DIPOLE_CDF_INV - Fatal error!");
                Console.WriteLine("  Iteration limit exceeded.");
                return 1;
            }

            if (cdf3 <= cdf && cdf1 <= cdf || cdf <= cdf3 && cdf <= cdf1)
            {
                x1 = x3;
                cdf1 = cdf3;
            }
            else
            {
                x2 = x3;
                cdf2 = cdf3;
            }

        }

        return x;
    }

    public static bool dipole_check(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIPOLE_CHECK checks the parameters of the Dipole CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    A is arbitrary, but represents an angle, so only 0 <= A <= 2 * PI
        //    is interesting, and -1.0 <= B <= 1.0.
        //
        //    Output, bool DIPOLE_CHECK, is true if the parameters are legal.
        //
    {
        switch (b)
        {
            case < -1.0:
            case > 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("DIPOLE_CHECK - Warning!");
                Console.WriteLine("  -1.0 <= B <= 1.0 is required.");
                Console.WriteLine("  The input B = " + b + "");
                return false;
            default:
                return true;
        }
    }

    public static double dipole_pdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIPOLE_PDF evaluates the Dipole PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B;X) =
        //        1 / ( PI * ( 1 + X^2 ) )
        //      + B^2 * ( ( 1 - X^2 ) * cos ( 2 * A ) + 2.0 * X * sin ( 2 * A ) )
        //      / ( PI * ( 1 + X )^2 )
        //
        //    Densities of this kind commonly occur in the analysis of resonant
        //    scattering of elementary particles.
        //
        //    DIPOLE_PDF(A,0;X) = CAUCHY_PDF(A;X)
        //    A = 0, B = 1 yields the single channel dipole distribution.
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
        //  Reference:
        //
        //    Robert Knop,
        //    Algorithm 441,
        //    ACM Transactions on Mathematical Software.
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //
        //    Input, double A, B, the parameters of the PDF.
        //    A is arbitrary, but represents an angle, so only 0 <= A <= 2 * PI
        //      is interesting,
        //    and -1.0 <= B <= 1.0.
        //
        //    Output, double DIPOLE_PDF, the value of the PDF.
        //
    {
        double pdf = 1.0 / (Math.PI * (1.0 + x * x))
                     + b * b * ((1.0 - x * x) * Math.Cos(2.0 * a)
                                + 2.0 * x * Math.Sin(2.0 * x))
                     / (Math.PI * (1.0 + x * x) * (1.0 + x * x));

        return pdf;
    }

    public static double dipole_sample(double a, double b, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIPOLE_SAMPLE samples the Dipole PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Robert Knop,
        //    Algorithm 441,
        //    ACM Transactions on Mathematical Software.
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    A is arbitrary, but represents an angle, so only 0 <= A <= 2 * PI
        //      is interesting,
        //    and -1.0 <= B <= 1.0.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double DIPOLE_SAMPLE, a sample of the PDF.
        //
    {
        //
        //  Find (X1,X2) at random in a circle.
        //
        double a2 = b * Math.Sin(a);
        double b2 = b * Math.Cos(a);
        double c2 = 1.0;

        double[] xc = Disk.disk_sample(a2, b2, c2, ref seed);
        //
        //  The dipole variate is the ratio X1 / X2.
        //
        double x = xc[0] / xc[1];

        return x;
    }
}