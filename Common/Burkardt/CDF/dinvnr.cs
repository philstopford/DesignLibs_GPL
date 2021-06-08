using System;

namespace Burkardt.CDFLib
{
    public static partial class CDF
    {
        public static double dinvnr(double p, double q)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DINVNR computes the inverse of the normal distribution.
            //
            //  Discussion:
            //
            //    Returns X such that CUMNOR(X)  =   P,  i.e., the  integral from -
            //    infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU is P
            //
            //    The rational function on page 95 of Kennedy and Gentle is used as a start
            //    value for the Newton method of finding roots.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 February 2021
            //
            //  Author:
            //
            //    Barry Brown, James Lovato, Kathy Russell.
            //
            //  Reference:
            //
            //    Kennedy and Gentle,
            //    Statistical Computing,
            //    Marcel Dekker, NY, 1980,
            //    QA276.4  K46
            //
            //  Parameters:
            //
            //    Input, double *P, *Q, the probability, and the complementary
            //    probability.
            //
            //    Output, double DINVNR, the argument X for which the
            //    Normal CDF has the value P.
            //
        {
            double maxit = 100;
            double eps = (1.0e-13);
            double r2pi = 0.3989422804014326e0;
            double nhalf = (-0.5e0);
            double dennor(double x)
            {
                return (r2pi * Math.Exp(nhalf * (x) * (x)));
            }

            double ccum;
            double cum;
            double dinvnr;
            double dx;
            int i;
            double pp;
            bool qporq;
            double strtx;
            double xcur;
            //
            //  FIND MINIMUM OF P AND Q
            //
            qporq = p <= q;
            if (!qporq) goto S10;
            pp = p;
            goto S20;
            S10:
            pp = q;
            S20:
            //
            //  INITIALIZATION STEP
            //
            strtx = stvaln(pp);
            xcur = strtx;
            //
            //  NEWTON INTERATIONS
            //
            for (i = 1; i <= maxit; i++)
            {
                cumnor(xcur, cum, ccum);
                dx = (cum - pp) / dennor(xcur);
                xcur = xcur - dx;
                if (Math.Abs(dx / xcur) < eps) goto S40;
            }

            dinvnr = strtx;
            //
            //  IF WE GET HERE, NEWTON HAS FAILED
            //
            if (!qporq) dinvnr = -dinvnr;
            return dinvnr;
            S40:
            //
            //  IF WE GET HERE, NEWTON HAS SUCCEDED
            //
            dinvnr = xcur;
            if (!qporq) dinvnr = -dinvnr;
            return dinvnr;
        }
    }
}