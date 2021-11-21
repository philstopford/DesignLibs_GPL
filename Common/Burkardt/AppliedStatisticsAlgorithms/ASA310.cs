using System;

namespace Burkardt.AppliedStatistics;

public static partial class Algorithms
{
    public static double ncbeta(double a, double b, double lambda, double x, double errmax,
            ref int ifault)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NCBETA computes the noncentral Beta CDF.
        //
        //  Discussion:
        //
        //    Three corrections needed to be made to the text of this routine.
        //    They are noted in the comments below.
        //
        //    Two of these corrections were errors in transcription made when
        //    producing the online copy distributed by APSTAT.
        //
        //    One error, an error of omission, occurred in the original printed
        //    copy of the routine, and was carried over into the online copy.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 February 2008
        //
        //  Author:
        //
        //    Original FORTRAN77 version by R Chattamvelli, R Shanmugam.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    R Chattamvelli, R Shanmugam,
        //    Algorithm AS 310:
        //    Computing the Non-central Beta Distribution Function,
        //    Applied Statistics,
        //    Volume 46, Number 1, 1997, pages 146-156.
        //
        //  Parameters:
        //
        //    Input, double A, B, the shape parameters.
        //    0 <= A, 0 <= B.
        //
        //    Input, double LAMBDA, the noncentrality parameter.
        //    0 <= LAMBDA.
        //
        //    Input, double X, the value at which the CDF is desired.
        //
        //    Input, double ERRMAX, the precision tolerance.
        //
        //    Output, int *IFAULT, error flag.
        //    0, no error occurred.
        //    1, X is 0 or 1.
        //    2, X < 0 or 1 < X.
        //    3, A, B or LAMBDA is less than 0.
        //
        //    Output, double NCBETA, the value of the noncentral Beta CDF.
        //
    {
        ifault = 0;
        double value = x;
        switch (lambda)
        {
            //
            //  Check parameters.
            //
            case <= 0.0:
                ifault = 3;
                return value;
        }

        switch (a)
        {
            case <= 0.0:
                ifault = 3;
                return value;
        }

        switch (b)
        {
            case <= 0.0:
                ifault = 3;
                return value;
        }

        switch (x)
        {
            case <= 0.0:
                value = 0.0;
                return value;
            case >= 1.0:
                value = 1.0;
                return value;
        }

        double c = 0.5 * lambda;
        switch (lambda)
        {
            //
            //  AS 226 as it stands is sufficient in this situation.
            //
            case < 54.0:
                value = betanc(x, a, b, lambda, ref ifault);
                return value;
        }

        int m = (int) (c + 0.5);
        double mr = m;
        int iterlo = m - (int) (5.0 * Math.Sqrt(mr));
        int iterhi = m + (int) (5.0 * Math.Sqrt(mr));
        double t = -c + mr * Math.Log(c) - Helpers.LogGamma(mr + 1.0);
        double q = Math.Exp(t);
        double r = q;
        double psum = q;

        double beta = Helpers.LogGamma(a + mr)
                      + Helpers.LogGamma(b)
                      - Helpers.LogGamma(a + mr + b);

        double s1 = (a + mr) * Math.Log(x)
            + b * Math.Log(1.0 - x) - Math.Log(a + mr) - beta;
        double gx = Math.Exp(s1);
        double fx = gx;
        double temp = betain(x, a + mr, b, beta, ref ifault);
        double ftemp = temp;
        //
        //  The online copy of AS 310 has "SUM = Q - TEMP" which is incorrect.
        //
        double sum = q * temp;
        int iter1 = m;
        //
        //  The first set of iterations starts from M and goes downwards
        //
        for (;;)
        {
            if (iter1 < iterlo)
            {
                break;
            }

            if (q < errmax)
            {
                break;
            }

            //
            //  The online copy of AS 310 has "Q = Q - ITER1 / C" which is incorrect.
            //
            q = q * iter1 / c;
            gx = (a + iter1) / (x * (a + b + iter1 - 1.0)) * gx;
            iter1 -= 1;
            temp += gx;
            psum += q;
            sum += q * temp;
        }

        double t0 = Helpers.LogGamma(a + b)
                    - Helpers.LogGamma(a + 1.0)
                    - Helpers.LogGamma(b);

        double s0 = a * Math.Log(x) + b * Math.Log(1.0 - x);
        //
        //  Both the online copy of AS 310 and the text printed in the reference
        //  did not initialize the variable S to zero, which is incorrect.
        //  JVB, 12 January 2008.
        //
        double s = 0.0;
        for (int i = 1; i <= iter1; i++)
        {
            int j = i - 1;
            s += Math.Exp(t0 + s0 + j * Math.Log(x));
            double t1 = Math.Log(a + b + j) - Math.Log(a + 1.0 + j) + t0;
            t0 = t1;
        }

        //
        //  Compute the first part of error bound.
        //
        double errbd = (1.0 - gammad(c, iter1, ref ifault)) * (temp + s);

        q = r;
        temp = ftemp;
        gx = fx;
        int iter2 = m;

        for (;;)
        {
            double ebd = errbd + (1.0 - psum) * temp;

            if (ebd < errmax || iterhi <= iter2)
            {
                break;
            }

            iter2 += 1;
            q = q * c / iter2;
            psum += q;
            temp -= gx;
            gx = x * (a + b + iter2 - 1.0) / (a + iter2) * gx;
            sum += q * temp;
        }

        value = sum;

        return value;
    }
}