﻿using System;

namespace Burkardt.AppliedStatistics;

public static partial class Algorithms
{
    public static void student_noncentral_cdf_values(ref int n_data, ref int df, ref double lambda,
            ref double x, ref double fx)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    STUDENT_NONCENTRAL_CDF_VALUES returns values of the noncentral Student CDF.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Needs["Statistics`ContinuousDistributions`"]
        //      dist = NoncentralStudentTDistribution [ df, lambda ]
        //      CDF [ dist, x ]
        //
        //    Mathematica seems to have some difficulty computing this function
        //    to the desired number of digits.
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
        //    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, int *DF, double *LAMBDA, the parameters of the
        //    function.
        //
        //    Output, double *X, the argument of the function.
        //
        //    Output, double *FX, the value of the function.
        //
    {
        const int N_MAX = 30;

        int[] df_vec =  {
                1, 2, 3,
                1, 2, 3,
                1, 2, 3,
                1, 2, 3,
                1, 2, 3,
                15, 20, 25,
                1, 2, 3,
                10, 10, 10,
                10, 10, 10,
                10, 10, 10
            }
            ;

        double[] fx_vec =  {
                0.8975836176504333E+00,
                0.9522670169E+00,
                0.9711655571887813E+00,
                0.8231218864E+00,
                0.9049021510E+00,
                0.9363471834E+00,
                0.7301025986E+00,
                0.8335594263E+00,
                0.8774010255E+00,
                0.5248571617E+00,
                0.6293856597E+00,
                0.6800271741E+00,
                0.20590131975E+00,
                0.2112148916E+00,
                0.2074730718E+00,
                0.9981130072E+00,
                0.9994873850E+00,
                0.9998391562E+00,
                0.168610566972E+00,
                0.16967950985E+00,
                0.1701041003E+00,
                0.9247683363E+00,
                0.7483139269E+00,
                0.4659802096E+00,
                0.9761872541E+00,
                0.8979689357E+00,
                0.7181904627E+00,
                0.9923658945E+00,
                0.9610341649E+00,
                0.8688007350E+00
            }
            ;

        double[] lambda_vec =  {
                0.0E+00,
                0.0E+00,
                0.0E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00,
                2.0E+00,
                2.0E+00,
                2.0E+00,
                4.0E+00,
                4.0E+00,
                4.0E+00,
                7.0E+00,
                7.0E+00,
                7.0E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00,
                2.0E+00,
                3.0E+00,
                4.0E+00,
                2.0E+00,
                3.0E+00,
                4.0E+00,
                2.0E+00,
                3.0E+00,
                4.0E+00
            }
            ;

        double[] x_vec =  {
                3.00E+00,
                3.00E+00,
                3.00E+00,
                3.00E+00,
                3.00E+00,
                3.00E+00,
                3.00E+00,
                3.00E+00,
                3.00E+00,
                3.00E+00,
                3.00E+00,
                3.00E+00,
                3.00E+00,
                3.00E+00,
                3.00E+00,
                15.00E+00,
                15.00E+00,
                15.00E+00,
                0.05E+00,
                0.05E+00,
                0.05E+00,
                4.00E+00,
                4.00E+00,
                4.00E+00,
                5.00E+00,
                5.00E+00,
                5.00E+00,
                6.00E+00,
                6.00E+00,
                6.00E+00
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
            df = 0;
            lambda = 0.0;
            x = 0.0;
            fx = 0.0;
        }
        else
        {
            df = df_vec[n_data - 1];
            lambda = lambda_vec[n_data - 1];
            x = x_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }
    }


    public static double tnc(double t, double df, double delta, ref int ifault)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TNC computes the tail of the noncentral T distribution.
        //
        //  Discussion:
        //
        //    This routine computes the cumulative probability at T of the
        //    non-central T-distribution with DF degrees of freedom (which may
        //    be fractional) and non-centrality parameter DELTA.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 January 2008
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Russell Lenth.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Russell Lenth,
        //    Algorithm AS 243:
        //    Cumulative Distribution Function of the Non-Central T Distribution,
        //    Applied Statistics,
        //    Volume 38, Number 1, 1989, pages 185-189.
        //
        //    William Guenther,
        //    Evaluation of probabilities for the noncentral distributions and
        //    difference of two T-variables with a desk calculator,
        //    Journal of Statistical Computation and Simulation,
        //    Volume 6, Number 3-4, 1978, pages 199-206.
        //
        //  Parameters:
        //
        //    Input, double T, the point whose cumulative probability
        //    is desired.
        //
        //    Input, double DF, the number of degrees of freedom.
        //
        //    Input, double DELTA, the noncentrality parameter.
        //
        //    Output, int *IFAULT, error flag.
        //    0, no error.
        //    nonzero, an error occcurred.
        //
        //    Output, double TNC, the tail of the noncentral
        //    T distribution.
        //
    {
        const double alnrpi = 0.57236494292470008707;
        const double errmax = 1.0E-10;
        const int itrmax = 100;
        const double r2pi = 0.79788456080286535588;

        double value = 0.0;

        switch (df)
        {
            case <= 0.0:
                ifault = 2;
                return value;
        }

        ifault = 0;

        double tt = t;
        double del = delta;
        bool negdel = false;

        switch (t)
        {
            case < 0.0:
                negdel = true;
                tt = -tt;
                del = -del;
                break;
        }

        //
        //  Initialize twin series.
        //
        double en = 1.0;
        double x = t * t / (t * t + df);

        switch (x)
        {
            case <= 0.0:
            {
                ifault = 0;
                value += alnorm(del, true);

                value = negdel switch
                {
                    true => 1.0 - value,
                    _ => value
                };

                return value;
            }
        }

        double lambda = del * del;
        double p = 0.5 * Math.Exp(-0.5 * lambda);
        double q = r2pi * p * del;
        double s = 0.5 - p;
        double a = 0.5;
        double b = 0.5 * df;
        double rxb = Math.Pow(1.0 - x, b);
        double albeta = alnrpi + Helpers.LogGamma(b) - Helpers.LogGamma(a + b);
        double xodd = betain(x, a, b, albeta, ref ifault);
        double godd = 2.0 * rxb * Math.Exp(a * Math.Log(x) - albeta);
        double xeven = 1.0 - rxb;
        double geven = b * x * rxb;
        value = p * xodd + q * xeven;
        //
        //  Repeat until convergence.
        //
        for (;;)
        {
            a += 1.0;
            xodd -= godd;
            xeven -= geven;
            godd = godd * x * (a + b - 1.0) / a;
            geven = geven * x * (a + b - 0.5) / (a + 0.5);
            p = p * lambda / (2.0 * en);
            q = q * lambda / (2.0 * en + 1.0);
            s -= p;
            en += 1.0;
            value = value + p * xodd + q * xeven;
            double errbd = 2.0 * s * (xodd - godd);

            if (errbd <= errmax)
            {
                ifault = 0;
                break;
            }

            if (!(itrmax < en))
            {
                continue;
            }

            ifault = 1;
            break;
        }

        value += alnorm(del, true);

        value = negdel switch
        {
            true => 1.0 - value,
            _ => value
        };

        return value;
    }

}