using System;
using Burkardt.Uniform;

namespace Burkardt.AppliedStatistics
{
    public static partial class Algorithms
    {
        public static double alnorm ( double x, bool upper )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ALNORM computes the cumulative density of the standard normal distribution.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 January 2008
        //
        //  Author:
        //
        //    Original FORTRAN77 version by David Hill.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    David Hill,
        //    Algorithm AS 66:
        //    The Normal Integral,
        //    Applied Statistics,
        //    Volume 22, Number 3, 1973, pages 424-427.
        //
        //  Parameters:
        //
        //    Input, double X, is one endpoint of the semi-infinite interval
        //    over which the integration takes place.
        //
        //    Input, bool UPPER, determines whether the upper or lower
        //    interval is to be integrated:
        //    .TRUE.  => integrate from X to + Infinity;
        //    .FALSE. => integrate from - Infinity to X.
        //
        //    Output, double ALNORM, the integral of the standard normal
        //    distribution over the desired interval.
        //
        {
            double a1 = 5.75885480458;
            double a2 = 2.62433121679;
            double a3 = 5.92885724438;
            double b1 = -29.8213557807;
            double b2 = 48.6959930692;
            double c1 = -0.000000038052;
            double c2 = 0.000398064794;
            double c3 = -0.151679116635;
            double c4 = 4.8385912808;
            double c5 = 0.742380924027;
            double c6 = 3.99019417011;
            double con = 1.28;
            double d1 = 1.00000615302;
            double d2 = 1.98615381364;
            double d3 = 5.29330324926;
            double d4 = -15.1508972451;
            double d5 = 30.789933034;
            double ltone = 7.0;
            double p = 0.398942280444;
            double q = 0.39990348504;
            double r = 0.398942280385;
            bool up;
            double utzero = 18.66;
            double value;
            double y;
            double z;

            up = upper;
            z = x;

            if ( z < 0.0 )
            {
                up = !up;
                z = - z;
            }

            if ( ltone < z && ( ( !up ) || utzero < z ) )
            {
                if ( up )
                {
                    value = 0.0;
                }
                else
                {
                    value = 1.0;
                }
                return value;
            }

            y = 0.5 * z * z;

            if ( z <= con )
            {
                value = 0.5 - z * ( p - q * y 
                / ( y + a1 + b1 
                / ( y + a2 + b2 
                / ( y + a3 ))));
            }
            else
            {
                value = r * Math.Exp ( - y ) 
                / ( z + c1 + d1 
                / ( z + c2 + d2 
                / ( z + c3 + d3 
                / ( z + c4 + d4 
                / ( z + c5 + d5 
                / ( z + c6 ))))));
            }

            if ( !up )
            {
                value = 1.0 - value;
            }

            return value;
        }
        
        public static double prncst ( double st, int idf, double d, ref int ifault )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PRNCST computes the lower tail of noncentral T distribution.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 January 2008
        //
        //  Author:
        //
        //    Original FORTRAN77 version by BE Cooper.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    BE Cooper,
        //    Algorithm AS 5:
        //    The Integral of the Non-Central T-Distribution,
        //    Applied Statistics,
        //    Volume 17, Number 2, 1968, page 193.
        //
        //  Parameters:
        //
        //    Input, double ST, the argument.
        //
        //    Input, int IDF, the number of degrees of freedom.
        //
        //    Input, double D, the noncentrality parameter.
        //
        //    Output, int *IFAULT, error flag.
        //    0, no error occurred.
        //    nonzero, an error occurred.
        //
        //    Output, double PRNCST, the value of the lower tail of
        //    the noncentral T distribution.
        //
        //  Local Parameters:
        //
        //    Local, double G1, 1.0 / sqrt(2.0 * pi)
        //
        //    Local, double G2, 1.0 / (2.0 * pi)
        //
        //    Local, double G3, sqrt(2.0 * pi)
        //
        {
            double a;
            double ak;
            double b;
            double da;
            double drb;
            double emin = 12.5;
            double f;
            double fk;
            double fkm1;
            double fmkm1;
            double fmkm2;
            double g1 = 0.3989422804;
            double g2 = 0.1591549431;
            double g3 = 2.5066282746;
            int ioe;
            int k;
            double rb;
            double sum;
            double value;

            f = ( double ) ( idf );
            //
            //  For very large IDF, use the normal approximation.
            //
            if ( 100 < idf )
            {
                ifault = 1;

                a = Math.Sqrt ( 0.5 * f ) * Math.Exp ( Helpers.LogGamma ( 0.5 * ( f - 1.0 ) ) 
                - Helpers.LogGamma ( 0.5 * f ) ) * d;

                value = alnorm ( ( st - a ) / Math.Sqrt ( f * ( 1.0 + d * d ) 
                / ( f - 2.0 ) - a * a ), false );
                return value;
            }

            ifault = 0;
            ioe = ( idf % 2 );
            a = st / Math.Sqrt ( f );
            b = f / ( f + st * st );
            rb = Math.Sqrt ( b );
            da = d * a;
            drb = d * rb;

            if ( idf == 1 )
            {
                value = alnorm ( drb, true ) + 2.0 * tfn ( drb, a );
                return value;
            }

            sum = 0.0;

            if ( Math.Abs ( drb ) < emin )
            {
                fmkm2 = a * rb * Math.Exp ( - 0.5 * drb * drb ) 
                * alnorm ( a * drb, false ) * g1;
            }
            else
            {
                fmkm2 = 0.0;
            }

            fmkm1 = b * da * fmkm2;

            if ( Math.Abs ( d ) < emin )
            {
                fmkm1 = fmkm1 + b * a * g2 * Math.Exp ( - 0.5 * d * d );
            }

            if ( ioe == 0 )
            {
                sum = fmkm2;
            }
            else
            {
                sum = fmkm1;
            }

            ak = 1.0;
            fk = 2.0;

            for ( k = 2; k <= idf - 2; k = k + 2 )
            {
                fkm1 = fk - 1.0;
                fmkm2 = b * ( da * ak * fmkm1 + fmkm2 ) * fkm1 / fk;
                ak = 1.0 / ( ak * fkm1 );
                fmkm1 = b * ( da * ak * fmkm2 + fmkm1 ) * fk / ( fk + 1.0 );

                if ( ioe == 0 )
                {
                    sum = sum + fmkm2;
                }
                else
                {
                    sum = sum + fmkm1;
                }
                ak = 1.0 / ( ak * fk );
                fk = fk + 2.0;
            }

            if ( ioe == 0 )
            {
                value = alnorm ( d, true ) + sum * g3;
            }
            else
            {
                value = alnorm ( drb, true ) + 2.0 * ( sum + tfn ( drb, a ) );
            }

            return value;
        }
        
        
        public static double tfn ( double x, double fx )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TFN calculates the T-function of Owen.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 January 2008
        //
        //  Author:
        //
        //    Original FORTRAN77 version by JC Young, Christoph Minder.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    MA Porter, DJ Winstanley,
        //    Remark AS R30:
        //    A Remark on Algorithm AS76:
        //    An Integral Useful in Calculating Noncentral T and Bivariate
        //    Normal Probabilities,
        //    Applied Statistics,
        //    Volume 28, Number 1, 1979, page 113.
        //
        //    JC Young, Christoph Minder,
        //    Algorithm AS 76: 
        //    An Algorithm Useful in Calculating Non-Central T and 
        //    Bivariate Normal Distributions,
        //    Applied Statistics,
        //    Volume 23, Number 3, 1974, pages 455-457.
        //
        //  Parameters:
        //
        //    Input, double X, FX, the parameters of the function.
        //
        //    Output, double TFN, the value of the T-function.
        //
        {
            int NG = 5;

            double fxs;
            int i;
            double[] r = {
            0.1477621, 
            0.1346334, 
            0.1095432, 
            0.0747257, 
            0.0333357 };
            double r1;
            double r2;
            double rt;
            double tp = 0.159155;
            double tv1 = 1.0E-35;
            double tv2 = 15.0;
            double tv3 = 15.0;
            double tv4 = 1.0E-05;
            double[] u = {
            0.0744372, 
            0.2166977, 
            0.3397048,
            0.4325317, 
            0.4869533 };
            double value;
            double x1;
            double x2;
            double xs;
            //
            //  Test for X near zero.
            //
            if ( Math.Abs ( x ) < tv1 )
            {
            value = tp * Math.Atan ( fx );
            return value;
            }
            //
            //  Test for large values of abs(X).
            //
            if ( tv2 < Math.Abs ( x ) )
            {
            value = 0.0;
            return value;
            }
            //
            //  Test for FX near zero.
            //
            if ( Math.Abs ( fx ) < tv1 )
            {
            value = 0.0;
            return value;
            }
            //
            //  Test whether abs ( FX ) is so large that it must be truncated.
            //
            xs = - 0.5 * x * x;
            x2 = fx;
            fxs = fx * fx;
            //
            //  Computation of truncation point by Newton iteration.
            //
            if ( tv3 <= Math.Log ( 1.0 + fxs ) - xs * fxs )
            {
            x1 = 0.5 * fx;
            fxs = 0.25 * fxs;

            for ( ; ; )
            {
            rt = fxs + 1.0;

            x2 = x1 + ( xs * fxs + tv3 - Math.Log ( rt ) ) 
            / ( 2.0 * x1 * ( 1.0 / rt - xs ) );

            fxs = x2 * x2;

            if ( Math.Abs ( x2 - x1 ) < tv4 )
            {
            break;
            }
            x1 = x2;
            }
            }
            //
            //  Gaussian quadrature.
            //
            rt = 0.0;
            for ( i = 0; i < NG; i++ )
            {
            r1 = 1.0 + fxs * Math.Pow ( 0.5 + u[i], 2 );
            r2 = 1.0 + fxs * Math.Pow ( 0.5 - u[i], 2 );

            rt = rt + r[i] * ( Math.Exp ( xs * r1 ) / r1 + Math.Exp ( xs * r2 ) / r2 );
            }

            value = rt * x2 * tp;

            return value;
        }

        public static void rnorm ( ref int seed, ref double u1, ref double u2 )
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    RNORM returns two independent standard random normal deviates.
            //
            //  Discussion:
            //
            //    This routine sets U1 and U2 to two independent standardized 
            //    random normal deviates.   This is a version of the 
            //    method given in Knuth.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 April 2014
            //
            //  Author:
            //
            //    Original FORTRAN77 version by William Smith, Ronald Hocking.
            //    This C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Donald Knuth,
            //    The Art of Computer Programming,
            //    Volume 2, Seminumerical Algorithms,
            //    Third Edition,
            //    Addison Wesley, 1997,
            //    ISBN: 0201896842,
            //    LC: QA76.6.K64.
            //
            //  Parameters:
            //
            //    Input/output, int &SEED, a seed for the random 
            //    number generator.
            //
            //    Output, double &U1, &U2, two standard random normal deviates.
            //
        {
            for ( ; ; )
            {
                double x = UniformRNG.r8_uniform_01 ( ref seed );
                double y = UniformRNG.r8_uniform_01 ( ref seed );
                x = 2.0 * x - 1.0;
                y = 2.0 * y - 1.0;
                double s = x * x + y * y;

                if ( s <= 1.0 )
                {
                    s = Math.Sqrt ( - 2.0 * Math.Log ( s ) / s );
                    u1 = x * s;
                    u2 = y * s;
                    break;
                }
            }
        }
        
    }
}