using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Probability;

public static class VonMises
{
    public static double von_mises_cdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VON_MISES_CDF evaluates the von Mises CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 November 2006
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Geoffrey Hill
        //    C++ version by John Burkardt
        //
        //  Reference:
        //
        //    Geoffrey Hill,
        //    ACM TOMS Algorithm 518,
        //    Incomplete Bessel Function I0: The von Mises Distribution,
        //    ACM Transactions on Mathematical Software,
        //    Volume 3, Number 3, September 1977, pages 279-284.
        //
        //    Kanti Mardia, Peter Jupp,
        //    Directional Statistics,
        //    Wiley, 2000, QA276.M335
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the CDF.
        //    A - PI <= X <= A + PI.
        //
        //    Input, double A, B, the parameters of the PDF.
        //    -PI <= A <= PI,
        //    0.0 < B.
        //
        //    Output, double VON_MISES_CDF, the value of the CDF.
        //
    {
        const double a1 = 12.0;
        const double a2 = 0.8;
        const double a3 = 8.0;
        const double a4 = 1.0;
        double c;
        const double c1 = 56.0;
        double cdf;
        const double ck = 10.5;

        double r;
        double s;
        double v;
        switch (x - a)
        {
            //
            //  We expect -PI <= X - A <= PI.
            //
            case <= -Math.PI:
                cdf = 0.0;
                return cdf;
            case >= Math.PI:
                cdf = 1.0;
                return cdf;
        }

        //
        //  Convert the angle (X - A) modulo 2 PI to the range ( 0, 2 * PI ).
        //
        double z = b;

        double u = typeMethods.r8_modp(x - a + Math.PI, 2.0 * Math.PI);

        switch (u)
        {
            case < 0.0:
                u += 2.0 * Math.PI;
                break;
        }

        double y = u - Math.PI;
        //
        //  For small B, sum IP terms by backwards recursion.
        //
        if (z <= ck)
        {
            v = 0.0;

            switch (z)
            {
                case > 0.0:
                {
                    int ip = (int) (z * a2 - a3 / (z + a4) + a1);
                    double p = ip;
                    s = Math.Sin(y);
                    c = Math.Cos(y);
                    y = p * y;
                    double sn = Math.Sin(y);
                    double cn = Math.Cos(y);
                    r = 0.0;
                    z = 2.0 / z;

                    int n;
                    for (n = 2; n <= ip; n++)
                    {
                        p -= 1.0;
                        y = sn;
                        sn = sn * c - cn * s;
                        cn = cn * c + y * s;
                        r = 1.0 / (p * z + r);
                        v = (sn / p + v) * r;
                    }

                    break;
                }
            }

            cdf = (u * 0.5 + v) / Math.PI;
        }
        //
        //  For large B, compute the normal approximation and left tail.
        //
        else
        {
            c = 24.0 * z;
            v = c - c1;
            r = Math.Sqrt((54.0 / (347.0 / v + 26.0 - c) - 6.0 + c) / 12.0);
            z = Math.Sin(0.5 * y) * r;
            s = 2.0 * z * z;
            v = v - s + 3.0;
            y = (c - s - s - 16.0) / 3.0;
            y = ((s + 1.75) * s + 83.5) / v - y;
            double arg = z * (1.0 - s / y / y);
            double erfx = typeMethods.r8_error_f(arg);
            cdf = 0.5 * erfx + 0.5;
        }

        cdf = Math.Max(cdf, 0.0);
        cdf = Math.Min(cdf, 1.0);

        return cdf;
    }

    public static double von_mises_cdf_inv(double cdf, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VON_MISES_CDF_INV inverts the von Mises CDF.
        //
        //  Discussion:
        //
        //    A simple bisection method is used on the interval [ A - PI, A + PI ].
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
        //
        //    Input, double A, B, the parameters of the PDF.
        //    -PI <= A <= PI,
        //    0.0 < B.
        //
        //    Output, double VON_MISES_CDF_INV, the corresponding argument of the CDF.
        //    A - PI <= X <= A + PI.
        //
    {
        const int it_max = 100;
            
        const double tol = 0.0001;
        double x;

        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("VON_MISES_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
            case 0.0:
                x = a - Math.PI;
                return x;
            case 1.0:
                x = a + Math.PI;
                return x;
        }

        double x1 = a - Math.PI;
        double cdf1 = 0.0;

        double x2 = a + Math.PI;
        //
        //  Now use bisection.
        //
        int it = 0;

        for (;;)
        {
            it += 1;

            double x3 = 0.5 * (x1 + x2);
            double cdf3 = von_mises_cdf(x3, a, b);

            if (Math.Abs(cdf3 - cdf) < tol)
            {
                x = x3;
                break;
            }

            if (it_max < it)
            {
                Console.WriteLine(" ");
                Console.WriteLine("VON_MISES_CDF_INV - Fatal error!");
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
            }
        }

        return x;
    }

    public static void von_mises_cdf_values(ref int n_data, ref double a, ref double b, ref double x,
            ref double fx )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VON_MISES_CDF_VALUES returns some values of the von Mises CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 December 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Kanti Mardia, Peter Jupp,
        //    Directional Statistics,
        //    Wiley, 2000, QA276.M335
        //
        //  Parameters:
        //
        //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, double &A, &B, the parameters of the function.
        //
        //    Output, double &X, the argument of the function.
        //
        //    Output, double &FX, the value of the function.
        //
    {
        const int N_MAX = 23;

        double[] a_vec =
        {
            0.0E+00,
            0.0E+00,
            0.0E+00,
            0.0E+00,
            0.0E+00,
            0.1E+01,
            0.1E+01,
            0.1E+01,
            0.1E+01,
            0.1E+01,
            0.1E+01,
            -0.2E+01,
            -0.1E+01,
            0.0E+01,
            0.1E+01,
            0.2E+01,
            0.3E+01,
            0.0E+00,
            0.0E+00,
            0.0E+00,
            0.0E+00,
            0.0E+00,
            0.0E+00
        };

        double[] b_vec =
        {
            0.1E+01,
            0.1E+01,
            0.1E+01,
            0.1E+01,
            0.1E+01,
            0.2E+01,
            0.2E+01,
            0.2E+01,
            0.2E+01,
            0.2E+01,
            0.2E+01,
            0.3E+01,
            0.3E+01,
            0.3E+01,
            0.3E+01,
            0.3E+01,
            0.3E+01,
            0.0E+00,
            0.1E+01,
            0.2E+01,
            0.3E+01,
            0.4E+01,
            0.5E+01
        };

        double[] fx_vec =
        {
            0.2535089956281180E-01,
            0.1097539041177346E+00,
            0.5000000000000000E+00,
            0.8043381312498558E+00,
            0.9417460124555197E+00,
            0.5000000000000000E+00,
            0.6018204118446155E+00,
            0.6959356933122230E+00,
            0.7765935901304593E+00,
            0.8410725934916615E+00,
            0.8895777369550366E+00,
            0.9960322705517925E+00,
            0.9404336090170247E+00,
            0.5000000000000000E+00,
            0.5956639098297530E-01,
            0.3967729448207649E-02,
            0.2321953958111930E-03,
            0.6250000000000000E+00,
            0.7438406999109122E+00,
            0.8369224904294019E+00,
            0.8941711407897124E+00,
            0.9291058600568743E+00,
            0.9514289900655436E+00
        };

        double[] x_vec =
        {
            -0.2617993977991494E+01,
            -0.1570796326794897E+01,
            0.0000000000000000E+00,
            0.1047197551196598E+01,
            0.2094395102393195E+01,
            0.1000000000000000E+01,
            0.1200000000000000E+01,
            0.1400000000000000E+01,
            0.1600000000000000E+01,
            0.1800000000000000E+01,
            0.2000000000000000E+01,
            0.0000000000000000E+00,
            0.0000000000000000E+00,
            0.0000000000000000E+00,
            0.0000000000000000E+00,
            0.0000000000000000E+00,
            0.0000000000000000E+00,
            0.7853981633974483E+00,
            0.7853981633974483E+00,
            0.7853981633974483E+00,
            0.7853981633974483E+00,
            0.7853981633974483E+00,
            0.7853981633974483E+00
        };

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
            b = 0.0;
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

    public static bool von_mises_check(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VON_MISES_CHECK checks the parameters of the von Mises PDF.
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
        //    -PI <= A <= PI,
        //    0.0 < B.
        //
        //    Output, bool VON_MISES_CHECK, is true if the parameters are legal.
        //
    {
            

        switch (a)
        {
            case < -Math.PI:
            case > Math.PI:
                Console.WriteLine(" ");
                Console.WriteLine("VON_MISES_CHECK - Warning!");
                Console.WriteLine("  A < -PI or PI < A.");
                return false;
        }

        switch (b)
        {
            case <= 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("VON_MISES_CHECK - Warning!");
                Console.WriteLine("  B <= 0.0");
                return false;
            default:
                return true;
        }
    }

    public static double von_mises_circular_variance(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VON_MISES_CIRCULAR_VARIANCE returns the circular variance of the von Mises PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 December 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    -PI <= A <= PI,
        //    0.0 < B.
        //
        //    Output, double VON_MISES_CIRCULAR_VARIANCE, the circular variance of the PDF.
        //
    {
        double value = 1.0 - Bessel.bessel_i1(b) / Bessel.bessel_i0(b);

        return value;
    }

    public static double von_mises_mean(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VON_MISES_MEAN returns the mean of the von Mises PDF.
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
        //    Input, double A, B, the parameters of the PDF.
        //    -PI <= A <= PI,
        //    0.0 < B.
        //
        //    Output, double VON_MISES_MEAN, the mean of the PDF.
        //
    {
        return a;
    }

    public static double von_mises_pdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VON_MISES_PDF evaluates the von Mises PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B;X) = EXP ( B * COS ( X - A ) ) / ( 2 * PI * I0(B) )
        //
        //    where:
        //
        //      I0(*) is the modified Bessel function of the first
        //      kind of order 0.
        //
        //    The von Mises distribution for points on the unit circle is
        //    analogous to the normal distribution of points on a line.
        //    The variable X is interpreted as a deviation from the angle A,
        //    with B controlling the amount of dispersion.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Jerry Banks, editor,
        //    Handbook of Simulation,
        //    Engineering and Management Press Books, 1998, page 160.
        //
        //    D J Best, N I Fisher,
        //    Efficient Simulation of the von Mises Distribution,
        //    Applied Statistics,
        //    Volume 28, Number 2, pages 152-157.
        //
        //    Kanti Mardia, Peter Jupp,
        //    Directional Statistics,
        //    Wiley, 2000, QA276.M335
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //    A - PI <= X <= A + PI.
        //
        //    Input, double A, B, the parameters of the PDF.
        //    -PI <= A <= PI,
        //    0.0 < B.
        //
        //    Output, double VON_MISES_PDF, the value of the PDF.
        //
    {
        double pdf;
            

        if (x < a - Math.PI)
        {
            pdf = 0.0;
        }
        else if (x <= a + Math.PI)
        {
            pdf = Math.Exp(b * Math.Cos(x - a)) / (2.0 * Math.PI * Bessel.bessel_i0(b));
        }
        else
        {
            pdf = 0.0;
        }

        return pdf;
    }

    public static double von_mises_sample(double a, double b, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VON_MISES_SAMPLE samples the von Mises PDF.
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
        //  Reference:
        //
        //    D J Best, N I Fisher,
        //    Efficient Simulation of the von Mises Distribution,
        //    Applied Statistics,
        //    Volume 28, Number 2, pages 152-157.
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    -PI <= A <= PI,
        //    0.0 < B.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double VON_MISES_SAMPLE, a sample of the PDF.
        //
    {
        double f;

        double tau = 1.0 + Math.Sqrt(1.0 + 4.0 * b * b);
        double rho = (tau - Math.Sqrt(2.0 * tau)) / (2.0 * b);
        double r = (1.0 + rho * rho) / (2.0 * rho);

        for (;;)
        {
            double u1 = UniformRNG.r8_uniform_01(ref seed);
            double z = Math.Cos(Math.PI * u1);
            f = (1.0 + r * z) / (r + z);
            double c = b * (r - f);

            double u2 = UniformRNG.r8_uniform_01(ref seed);

            if (u2 < c * (2.0 - c))
            {
                break;
            }

            if (c <= Math.Log(c / u2) + 1.0)
            {
                break;
            }

        }

        double u3 = UniformRNG.r8_uniform_01(ref seed);

        double x = a + typeMethods.r8_sign(u3 - 0.5) * Math.Acos(f);

        return x;
    }
}