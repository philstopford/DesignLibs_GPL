using System;

namespace Burkardt.CDFLib
{
    public static partial class CDF
    {
        public static void cumnor(double arg, ref double result, ref double ccum)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CUMNOR computes the cumulative normal distribution.
            //
            //  Discussion:
            //
            //    This function evaluates the normal distribution function:
            //
            //                              / x
            //                     1       |       -t*t/2
            //          P(x) = ----------- |      e       dt
            //                 sqrt(2 pi)  |
            //                             /-oo
            //
            //    This transportable program uses rational functions that
            //    theoretically approximate the normal distribution function to
            //    at least 18 significant decimal digits.  The accuracy achieved
            //    depends on the arithmetic system, the compiler, the intrinsic
            //    functions, and proper selection of the machine-dependent
            //    constants.
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
            //    William Cody
            //
            //  Reference:
            //
            //    William Cody,
            //    Rational Chebyshev approximations for the error function,
            //    Mathematics of Computation,
            //    1969, pages 631-637.
            //
            //    William Cody,
            //    Algorithm 715:
            //    SPECFUN - A Portable FORTRAN Package of Special Function Routines
            //      and Test Drivers,
            //    ACM Transactions on Mathematical Software,
            //    Volume 19, 1993, pages 22-32.
            //
            //  Parameters:
            //
            //    Input, double *ARG, the upper limit of integration.
            //
            //    Output, double *CUM, *CCUM, the Normal density CDF and
            //    complementary CDF.
            //
            //  Local:
            //
            //    double EPS, the argument below which anorm(x)
            //    may be represented by 0.5D+00 and above which  x*x  will not underflow.
            //    A conservative value is the largest machine number X
            //    such that   1.0D+00 + X = 1.0D+00   to machine precision.
            //
        {
            double[] a =  {
                2.2352520354606839287e00,1.6102823106855587881e02,1.0676894854603709582e03,
                1.8154981253343561249e04,6.5682337918207449113e-2
            }
            ;
            double[] b =  {
                4.7202581904688241870e01,9.7609855173777669322e02,1.0260932208618978205e04,
                4.5507789335026729956e04
            }
            ;
            double[] c =  {
                3.9894151208813466764e-1,8.8831497943883759412e00,9.3506656132177855979e01,
                5.9727027639480026226e02,2.4945375852903726711e03,6.8481904505362823326e03,
                1.1602651437647350124e04,9.8427148383839780218e03,1.0765576773720192317e-8
            }
            ;
            double[] d =  {
                2.2266688044328115691e01,2.3538790178262499861e02,1.5193775994075548050e03,
                6.4855582982667607550e03,1.8615571640885098091e04,3.4900952721145977266e04,
                3.8912003286093271411e04,1.9685429676859990727e04
            }
            ;
            double half = 0.5e0;
            double[] p =  {
                2.1589853405795699e-1,1.274011611602473639e-1,2.2235277870649807e-2,
                1.421619193227893466e-3,2.9112874951168792e-5,2.307344176494017303e-2
            }
            ;
            double one = 1.0e0;
            double[] q =  {
                1.28426009614491121e00,4.68238212480865118e-1,6.59881378689285515e-2,
                3.78239633202758244e-3,7.29751555083966205e-5
            }
            ;
            double sixten = 1.60e0;
            double sqrpi = 3.9894228040143267794e-1;
            double thrsh = 0.66291e0;
            double root32 = 5.656854248e0;
            double zero = 0.0e0;
            int K1 = 1;
            int K2 = 2;
            int i;
            double del, eps, temp, x, xden, xnum, y, xsq, min;
            //
            //  Machine dependent constants
            //
            eps = dpmpar(K1) * 0.5e0;
            min = dpmpar(K2);
            x = arg;
            y = Math.Abs(x);
            if (y <= thrsh)
            {
                //
                //  Evaluate  anorm  for  |X| <= 0.66291
                //
                xsq = zero;
                if (y > eps) xsq = x * x;
                xnum = a[4] * xsq;
                xden = xsq;
                for (i = 0; i < 3; i++)
                {
                    xnum = (xnum + a[i]) * xsq;
                    xden = (xden + b[i]) * xsq;
                }

                result = x * (xnum + a[3]) / (xden + b[3]);
                temp = result;
                result = half + temp;
                ccum = half - temp;
            }
            //
            //  Evaluate  anorm  for 0.66291 <= |X| <= sqrt(32)
            //
            else if (y <= root32)
            {
                xnum = c[8] * y;
                xden = y;
                for (i = 0; i < 7; i++)
                {
                    xnum = (xnum + c[i]) * y;
                    xden = (xden + d[i]) * y;
                }

                result = (xnum + c[7]) / (xden + d[7]);
                xsq = Math.Truncate(y * sixten) / sixten;
                del = (y - xsq) * (y + xsq);
                result = Math.Exp(-(xsq * xsq * half)) * Math.Exp(-(del * half)) * result;
                ccum = one - result;
                if (x > zero)
                {
                    temp = result;
                    result = ccum;
                    ccum = temp;
                }
            }
            //
            //  Evaluate  anorm  for |X| > sqrt(32)
            //
            else
            {
                result = zero;
                xsq = one / (x * x);
                xnum = p[5] * xsq;
                xden = xsq;
                for (i = 0; i < 4; i++)
                {
                    xnum = (xnum + p[i]) * xsq;
                    xden = (xden + q[i]) * xsq;
                }

                result = xsq * (xnum + p[4]) / (xden + q[4]);
                result = (sqrpi - result) / y;
                xsq = Math.Truncate(x * sixten) / sixten;
                del = (x - xsq) * (x + xsq);
                result = Math.Exp(-(xsq * xsq * half)) * Math.Exp(-(del * half)) * result;
                ccum = one - result;
                if (x > zero)
                {
                    temp = result;
                    result = ccum;
                    ccum = temp;
                }
            }

            if (result < min) result = 0.0e0;
            //
            //  Fix up for negative argument, erf, etc.
            //
            if (ccum < min) ccum = 0.0e0;
        }
    }
}