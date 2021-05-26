using System;

namespace Burkardt.AppliedStatistics
{
    public static partial class Algorithms
    {

        public static void normp(double z, ref double p, ref double q, ref double pdf)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMP computes the cumulative density of the standard normal distribution.
        //
        //  Discussion:
        //
        //    This is algorithm 5666 from Hart, et al.
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
        //    Original FORTRAN77 version by Alan Miller.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, 
        //    Charles Mesztenyi, John Rice, Henry Thacher, 
        //    Christoph Witzgall,
        //    Computer Approximations,
        //    Wiley, 1968,
        //    LC: QA297.C64.
        //
        //  Parameters:
        //
        //    Input, double Z, divides the real line into two 
        //    semi-infinite intervals, over each of which the standard normal 
        //    distribution is to be integrated.
        //
        //    Output, double *P, *Q, the integrals of the standard normal
        //    distribution over the intervals ( - Infinity, Z] and 
        //    [Z, + Infinity ), respectively.
        //
        //    Output, double *PDF, the value of the standard normal distribution
        //    at Z.
        //
        {
            double cutoff = 7.071;
            double p0 = 220.2068679123761;
            double p1 = 221.2135961699311;
            double p2 = 112.0792914978709;
            double p3 = 33.91286607838300;
            double p4 = 6.373962203531650;
            double p5 = 0.7003830644436881;
            double p6 = 0.03526249659989109;
            double q0 = 440.4137358247522;
            double q1 = 793.8265125199484;
            double q2 = 637.3336333788311;
            double q3 = 296.5642487796737;
            double q4 = 86.78073220294608;
            double q5 = 16.06417757920695;
            double q6 = 1.755667163182642;
            double q7 = 0.08838834764831844;
            double root2pi = 2.506628274631001;

            double zabs = Math.Abs(z);
            //
            //  37 < |Z|.
            //
            if (37.0 < zabs)
            {
                pdf = 0.0;
                p = 0.0;
            }
            //
            //  |Z| <= 37.
            //
            else
            {
                double expntl = Math.Exp(-0.5 * zabs * zabs);
                pdf = expntl / root2pi;
                //
                //  |Z| < CUTOFF = 10 / sqrt(2).
                //
                if (zabs < cutoff)
                {
                    p = expntl * ((((((
                                           p6 * zabs
                                           + p5) * zabs
                                       + p4) * zabs
                                      + p3) * zabs
                                     + p2) * zabs
                                    + p1) * zabs
                                   + p0) / (((((((
                                                     q7 * zabs
                                                     + q6) * zabs
                                                 + q5) * zabs
                                                + q4) * zabs
                                               + q3) * zabs
                                              + q2) * zabs
                                             + q1) * zabs
                                            + q0);
                }
                //
                //  CUTOFF <= |Z|.
                //
                else
                {
                    p = pdf / (
                        zabs + 1.0 / (
                            zabs + 2.0 / (
                                zabs + 3.0 / (
                                    zabs + 4.0 / (
                                        zabs + 0.65)))));
                }
            }

            if (z < 0.0)
            {
                q = 1.0 - p;
            }
            else
            {
                q = p;
                p = 1.0 - q;
            }
        }

        public static void nprob(double z, ref double p, ref double q, ref double pdf)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NPROB computes the cumulative density of the standard normal distribution.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 June 2013
        //
        //  Author:
        //
        //    Original FORTRAN77 version by AG Adams.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    AG Adams,
        //    Algorithm 39:
        //    Areas Under the Normal Curve,
        //    Computer Journal,
        //    Volume 12, Number 2, May 1969, pages 197-198.
        //
        //  Parameters:
        //
        //    Input, double Z, divides the real line into 
        //    two semi-infinite intervals, over each of which the standard normal 
        //    distribution is to be integrated.
        //
        //    Output, double *P, *Q, the integrals of the standard normal
        //    distribution over the intervals ( - Infinity, Z] and 
        //    [Z, + Infinity ), respectively.
        //
        //    Output, double *PDF, the value of the standard normal
        //    distribution at Z.
        //
        {
            double a0 = 0.5;
            double a1 = 0.398942280444;
            double a2 = 0.399903438504;
            double a3 = 5.75885480458;
            double a4 = 29.8213557808;
            double a5 = 2.62433121679;
            double a6 = 48.6959930692;
            double a7 = 5.92885724438;
            double b0 = 0.398942280385;
            double b1 = 0.000000038052;
            double b2 = 1.00000615302;
            double b3 = 0.000398064794;
            double b4 = 1.98615381364;
            double b5 = 0.151679116635;
            double b6 = 5.29330324926;
            double b7 = 4.8385912808;
            double b8 = 15.1508972451;
            double b9 = 0.742380924027;
            double b10 = 30.789933034;
            double b11 = 3.99019417011;
            double y;
            double zabs;

            zabs = Math.Abs(z);
            //
            //  |Z| between 0 and 1.28
            //
            if (zabs <= 1.28)
            {
                y = a0 * z * z;
                pdf = Math.Exp(-y) * b0;

                q = a0 - zabs * (a1 - a2 * y
                    / (y + a3 - a4
                        / (y + a5 + a6
                            / (y + a7))));
            }
            //
            //  |Z| between 1.28 and 12.7
            //
            else if (zabs <= 12.7)
            {
                y = a0 * z * z;
                pdf = Math.Exp(-y) * b0;

                q = pdf
                     / (zabs - b1 + b2
                         / (zabs + b3 + b4
                             / (zabs - b5 + b6
                                 / (zabs + b7 - b8
                                     / (zabs + b9 + b10
                                         / (zabs + b11))))));
            }
            //
            //  Z far out in tail.
            //
            else
            {
                q = 0.0;
                pdf = 0.0;
            }

            if (z < 0.0)
            {
                p = q;
                q = 1.0 - p;
            }
            else
            {
                p = 1.0 - q;
            }
        }


    }
}