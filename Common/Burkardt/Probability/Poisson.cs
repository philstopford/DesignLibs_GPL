using System;
using Burkardt.Types;

namespace Burkardt.Probability;

public static class Poisson
{
    public static double poisson_cdf(int k, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POISSON_CDF evaluates the Poisson CDF.
        //
        //  Discussion:
        //
        //    CDF(K,A) is the probability that the number of events observed
        //    in a unit time period will be no greater than K, given that the
        //    expected number of events in a unit time period is A.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 January 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int K, the argument of the CDF.
        //
        //    Input, double A, the parameter of the PDF.
        //    0 < A.
        //
        //    Output, double POISSON_CDF, the value of the CDF.
        //
    {
        double cdf;
        switch (a)
        {
            //
            //  Check.
            //
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("POISSON_CDF - Fatal error!");
                Console.WriteLine("  A <= 0.");
                return 0;
        }

        switch (k)
        {
            //
            //  Special cases.
            //
            case < 0:
                cdf = 0.0;
                return cdf;
        }

        //
        //  General case.
        //
        double next = Math.Exp(-a);
        cdf = next;

        for (int i = 1; i <= k; i++)
        {
            double last = next;
            next = last * a / i;
            cdf += next;
        }

        return cdf;
    }

    public static int poisson_cdf_inv(double cdf, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POISSON_CDF_INV inverts the Poisson CDF.
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
        //    Input, double CDF, a value of the CDF.
        //    0 <= CDF < 1.
        //
        //    Input, double A, the parameter of the PDF.
        //    0.0 < A.
        //
        //    Output, int POISSON_CDF_INV, the corresponding argument.
        //
    {
        int x = 0;
        int xmax = 100;

        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("POISSON_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
        }

        //
        //  Now simply start at X = 0, and find the first value for which
        //  CDF(X-1) <= CDF <= CDF(X).
        //
        double sum2 = 0.0;

        for (int i = 0; i <= xmax; i++)
        {
            double sumold = sum2;

            double newval = 0;
            switch (i)
            {
                case 0:
                    newval = Math.Exp(-a);
                    sum2 = newval;
                    break;
                default:
                {
                    double last = newval;
                    newval = last * a / i;
                    sum2 += newval;
                    break;
                }
            }

            if (sumold <= cdf && cdf <= sum2)
            {
                x = i;
                return x;
            }
        }

        switch (x)
        {
            case > 100:
                Console.WriteLine(" ");
                Console.WriteLine("POISSON_CDF_INV - Warning!");
                Console.WriteLine("  Exceeded XMAX = " + xmax + "");
                x = xmax;
                break;
        }
            
        return x;
    }

    public static void poisson_cdf_values(ref int n_data, ref double a, ref int x, ref double fx )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POISSON_CDF_VALUES returns some values of the Poisson CDF.
        //
        //  Discussion:
        //
        //    CDF(X)(A) is the probability of at most X successes in unit time,
        //    given that the expected mean number of successes is A.
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Needs["Statistics`DiscreteDistributions`]
        //      dist = PoissonDistribution [ a ]
        //      CDF [ dist, x ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 August 2004
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
        //    Daniel Zwillinger,
        //    CRC Standard Mathematical Tables and Formulae,
        //    30th Edition, CRC Press, 1996, pages 653-658.
        //
        //  Parameters:
        //
        //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, double &A, the parameter of the function.
        //
        //    Output, int *X, the argument of the function.
        //
        //    Output, double &FX, the value of the function.
        //
    {
        const int N_MAX = 21;

        double[] a_vec =
        {
            0.02E+00,
            0.10E+00,
            0.10E+00,
            0.50E+00,
            0.50E+00,
            0.50E+00,
            1.00E+00,
            1.00E+00,
            1.00E+00,
            1.00E+00,
            2.00E+00,
            2.00E+00,
            2.00E+00,
            2.00E+00,
            5.00E+00,
            5.00E+00,
            5.00E+00,
            5.00E+00,
            5.00E+00,
            5.00E+00,
            5.00E+00
        };

        double[] fx_vec =
        {
            0.9801986733067553E+00,
            0.9048374180359596E+00,
            0.9953211598395555E+00,
            0.6065306597126334E+00,
            0.9097959895689501E+00,
            0.9856123220330293E+00,
            0.3678794411714423E+00,
            0.7357588823428846E+00,
            0.9196986029286058E+00,
            0.9810118431238462E+00,
            0.1353352832366127E+00,
            0.4060058497098381E+00,
            0.6766764161830635E+00,
            0.8571234604985470E+00,
            0.6737946999085467E-02,
            0.4042768199451280E-01,
            0.1246520194830811E+00,
            0.2650259152973617E+00,
            0.4404932850652124E+00,
            0.6159606548330631E+00,
            0.7621834629729387E+00
        };

        int[] x_vec =
        {
            0, 0, 1, 0,
            1, 2, 0, 1,
            2, 3, 0, 1,
            2, 3, 0, 1,
            2, 3, 4, 5,
            6
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
            x = 0;
            fx = 0.0;
        }
        else
        {
            a = a_vec[n_data - 1];
            x = x_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }
    }

    public static bool poisson_check(double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POISSON_CHECK checks the parameter of the Poisson PDF.
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
        //    Output, bool POISSON_CHECK, is true if the parameters are legal.
        //
    {
        switch (a)
        {
            case <= 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("POISSON_CHECK - Warning!");
                Console.WriteLine("  A <= 0.");
                return false;
            default:
                return true;
        }
    }

    public static double poisson_kernel(double r, int n, double[] c, double[] x, double[] y )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POISSON_KERNEL evaluates the Poisson kernel.
        //
        //  Discussion:
        //
        //    P(X,Y) = ( R^2 - |X-C|^2 ) / ( R * A * |X-Y|^N )
        //
        //    where the N-dimensional ball has radius R and center C,
        //    and A is the area of the unit sphere.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 November 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the ball.
        //
        //    Input, int N, the spatial dimension.
        //
        //    Input, double C(N), the center of the ball.
        //
        //    Input, double X(N), a point inside the ball.
        //
        //    Input, double Y(N), a point on the surface of the ball.
        //
        //    Output, double POISSON_KERNEL, the Poisson kernel function P(X,Y).
        //
    {
        double area;
        double b;
        double p;
        double t;
        double xc_diff_norm;
        double xy_diff_norm;

        xc_diff_norm = typeMethods.r8vec_diff_norm(n, x, c);
        xy_diff_norm = typeMethods.r8vec_diff_norm(n, x, y);
        area = Misc.sphere_unit_area_nd(n);

        t = (r + xc_diff_norm) * (r - xc_diff_norm);
        b = r * area * Math.Pow(xy_diff_norm, n);
        p = t / b;

        return p;
    }

    public static double poisson_mean(double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POISSON_MEAN returns the mean of the Poisson PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 February 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, the parameter of the PDF.
        //    0 < A.
        //
        //    Output, double POISSON_MEAN, the mean of the PDF.
        //
    {
        switch (a)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("POISSON_MEAN - Fatal error!");
                Console.WriteLine("  A <= 0.");
                return 0;
            default:
                return a;
        }
    }

    public static double poisson_pdf(int k, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POISSON_PDF evaluates the Poisson PDF.
        //
        //  Discussion:
        //
        //    PDF(K,A) is the probability that the number of events observed
        //    in a unit time period will be K, given the expected number
        //    of events in a unit time.
        //
        //    The parameter A is the expected number of events per unit time.
        //
        //    The Poisson PDF is a discrete version of the exponential PDF.
        //
        //    The time interval between two Poisson events is a random
        //    variable with the exponential PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 February 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int K, the argument of the PDF.
        //
        //    Input, double A, the parameter of the PDF.
        //    0 < A.
        //
        //    Output, double POISSON_PDF, the value of the PDF.
        //
    {
        double pdf;
        switch (a)
        {
            //
            //  Check.
            //
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("POISSON_PDF - Fatal error!");
                Console.WriteLine("  A <= 0.");
                pdf = 0.0;
                return pdf;
            default:
                pdf = Math.Exp(-a) * Math.Pow(a, k) / typeMethods.r8_factorial(k);

                return pdf;
        }
    }

    public static int poisson_sample(double a, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POISSON_SAMPLE samples the Poisson PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 February 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, the parameter of the PDF.
        //    0 < A.
        //
        //    Input/output, int &SEED, the random number generator seed.
        //
        //    Output, int POISSON_SAMPLE, a sample of the PDF.
        //
    {
        int KMAX = 100;
        switch (a)
        {
            //
            //  Check.
            //
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("POISSON_SAMPLE - Fatal error!");
                Console.WriteLine("  A <= 0.");
                return 0;
        }

        //
        //  Pick a random value of CDF.
        //
        double cdf = Uniform.uniform_01_sample(ref seed);
        //
        //  Now simply start at K = 0, and find the first value for which
        //  CDF(K-1) <= CDF <= CDF(K).
        //
        double sum = 0.0;
        double next = 0;

        for (int i = 0; i <= KMAX; i++)
        {
            double sumold = sum;

            switch (i)
            {
                case 0:
                    next = Math.Exp(-a);
                    sum = next;
                    break;
                default:
                {
                    double last = next;
                    next = last * a / i;
                    sum += next;
                    break;
                }
            }

            if (sumold <= cdf && cdf <= sum)
            {
                return i;
            }

        }

        Console.WriteLine("");
        Console.WriteLine("POISSON_SAMPLE - Warning!");
        Console.WriteLine("  Exceeded KMAX = " + KMAX + "");

        return KMAX;
    }

    public static double poisson_variance(double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POISSON_VARIANCE returns the variance of the Poisson PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 February 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, the parameter of the PDF.
        //    0 < A.
        //
        //    Output, double POISSON_VARIANCE, the variance of the PDF.
        //
    {
        switch (a)
        {
            //
            //  Check.
            //
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("POISSON_VARIANCE - Fatal error!");
                Console.WriteLine("  A <= 0.");
                return 0.0;
            default:
                return a;
        }
    }
}