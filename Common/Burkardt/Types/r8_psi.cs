using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double r8_psi(double xx)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_PSI evaluates the psi function.
        //
        //  Discussion:
        //
        //    This routine evaluates the logarithmic derivative of the
        //    Gamma function,
        //
        //      PSI(X) = d/dX ( GAMMA(X) ) / GAMMA(X)
        //             = d/dX LN ( GAMMA(X) )
        //
        //    for real X, where either
        //
        //      - XMAX1 < X < - XMIN, and X is not a negative integer,
        //
        //    or
        //
        //      XMIN < X.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 June 2013
        //
        //  Author:
        //
        //    Original FORTRAN77 version by William Cody, Anthony Strecok, Henry Thacher.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    William Cody, Anthony Strecok, Henry Thacher,
        //    Chebyshev Approximations for the Psi Function,
        //    Mathematics of Computation,
        //    Volume 27, Number 121, January 1973, pages 123-127.
        //
        //  Parameters:
        //
        //    Input, double XX, the argument of the psi function.
        //
        //    Output, double R8_PSI, the value of the psi function.  
        //    PSI is assigned the value 0 when the psi function is undefined.
        //
        {
            double den;
            double fourth = 0.25;
            double[] p1 =
            {
                4.5104681245762934160E-03,
                5.4932855833000385356E+00,
                3.7646693175929276856E+02,
                7.9525490849151998065E+03,
                7.1451595818951933210E+04,
                3.0655976301987365674E+05,
                6.3606997788964458797E+05,
                5.8041312783537569993E+05,
                1.6585695029761022321E+05
            };
            double[] p2 =
            {
                -2.7103228277757834192E+00,
                -1.5166271776896121383E+01,
                -1.9784554148719218667E+01,
                -8.8100958828312219821E+00,
                -1.4479614616899842986E+00,
                -7.3689600332394549911E-02,
                -6.5135387732718171306E-21
            };
            double piov4 = 0.78539816339744830962;
            double[] q1 =
            {
                9.6141654774222358525E+01,
                2.6287715790581193330E+03,
                2.9862497022250277920E+04,
                1.6206566091533671639E+05,
                4.3487880712768329037E+05,
                5.4256384537269993733E+05,
                2.4242185002017985252E+05,
                6.4155223783576225996E-08
            };
            double[] q2 =
            {
                4.4992760373789365846E+01,
                2.0240955312679931159E+02,
                2.4736979003315290057E+02,
                1.0742543875702278326E+02,
                1.7463965060678569906E+01,
                8.8427520398873480342E-01
            };
            double sgn;
            double three = 3.0;
            double upper;
            double value;
            double w;
            double x;
            double x01 = 187.0E+00;
            double x01d = 128.0E+00;
            double x02 = 6.9464496836234126266E-04;
            double xlarge = 2.04E+15;
            double xmax1 = 3.60E+16;
            double xmin1 = 5.89E-39;
            double xsmall = 2.05E-09;
            double z;

            x = xx;
            w = Math.Abs(x);
            double aug = 0.0;
            //
            //  Check for valid arguments, then branch to appropriate algorithm.
            //
            if (xmax1 <= -x || w < xmin1)
            {
                if (0.0 < x)
                {
                    value = -typeMethods.r8_huge();
                }
                else
                {
                    value = typeMethods.r8_huge();
                }

                return value;
            }

            //
            //  X < 0.5, use reflection formula: psi(1-x) = psi(x) + Math.PI * cot(pi*x)
            //  Use 1/X for PI*COTAN(PI*X)  when  XMIN1 < |X| <= XSMALL.
            //
            if (x < 0.5)
            {
                if (w <= xsmall)
                {
                    aug = -1.0 / x;
                }
                //
                //  Argument reduction for cotangent.
                //
                else
                {
                    if (x < 0.0)
                    {
                        sgn = piov4;
                    }
                    else
                    {
                        sgn = -piov4;
                    }

                    w = w - (double) ((int) (w));
                    int nq = (int) (w * 4.0);
                    w = 4.0 * (w - (double) (nq) * fourth);
                    //
                    //  W is now related to the fractional part of 4.0 * X.
                    //  Adjust argument to correspond to values in the first
                    //  quadrant and determine the sign.
                    //
                    int n = nq / 2;

                    if (n + n != nq)
                    {
                        w = 1.0 - w;
                    }

                    z = piov4 * w;

                    if ((n % 2) != 0)
                    {
                        sgn = -sgn;
                    }

                    //
                    //  Determine the final value for  -pi * cotan(pi*x).
                    //
                    n = (nq + 1) / 2;
                    if ((n % 2) == 0)
                    {
                        //
                        //  Check for singularity.
                        //
                        if (z == 0.0)
                        {
                            if (0.0 < x)
                            {
                                value = -typeMethods.r8_huge();
                            }
                            else
                            {
                                value = typeMethods.r8_huge();
                            }

                            return value;
                        }

                        aug = sgn * (4.0 / Math.Tan(z));
                    }
                    else
                    {
                        aug = sgn * (4.0 * Math.Tan(z));
                    }
                }

                x = 1.0 - x;
            }

            //
            //  0.5 <= X <= 3.0.
            //
            if (x <= three)
            {
                den = x;
                upper = p1[0] * x;
                for (int i = 0; i < 7; i++)
                {
                    den = (den + q1[i]) * x;
                    upper = (upper + p1[i + 1]) * x;
                }

                den = (upper + p1[8]) / (den + q1[7]);
                x = (x - x01 / x01d) - x02;
                value = den * x + aug;
                return value;
            }

            //
            //  3.0 < X.
            //
            if (x < xlarge)
            {
                w = 1.0 / (x * x);
                den = w;
                upper = p2[0] * w;
                for (int i = 0; i < 5; i++)
                {
                    den = (den + q2[i]) * w;
                    upper = (upper + p2[i + 1]) * w;
                }

                aug = (upper + p2[6]) / (den + q2[5]) - 0.5 / x + aug;
            }

            value = aug + Math.Log(x);

            return value;
        }
        
    }
}