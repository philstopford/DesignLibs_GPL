using System;

namespace Burkardt.CDFLib;

public static class bivariatenormal
{
    public class BivnorData
    {
        public const int idig = 15;
    }
    public static double bivnor( ref BivnorData data, double ah, double ak, double r)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BIVNOR computes the bivariate normal CDF.
        //
        //  Discussion:
        //
        //    BIVNOR computes the probability for two normal variates X and Y
        //    whose correlation is R, that AH <= X and AK <= Y.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 April 2012
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Thomas Donnelly.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Thomas Donnelly,
        //    Algorithm 462: Bivariate Normal Distribution,
        //    Communications of the ACM,
        //    October 1973, Volume 16, Number 10, page 638.
        //
        //  Parameters:
        //
        //    Input, double AH, AK, the lower limits of integration.
        //
        //    Input, double R, the correlation between X and Y.
        //
        //    Output, double BIVNOR, the bivariate normal CDF.
        //
        //  Local Parameters:
        //
        //    Local, int IDIG, the number of significant digits
        //    to the right of the decimal point desired in the answer.
        //
    {
        double con;
        double gw = 0;
        int is_ = 0;
        const double twopi = Math.PI * 2;
        double wh = 0;
        double wk = 0;

        double b = 0.0;

        double gh = Gauss.gauss(-ah) / 2.0;
        double gk = Gauss.gauss(-ak) / 2.0;

        switch (r)
        {
            case 0.0:
                b = 4.00 * gh * gk;
                b = Math.Max(b, 0.0);
                b = Math.Min(b, 1.0);
                return b;
        }

        double rr = (1.0 + r) * (1.0 - r);

        switch (rr)
        {
            case < 0.0:
                Console.WriteLine("");
                Console.WriteLine("BIVNOR - Fatal error!");
                Console.WriteLine("  1 < |R|.");
                break;
            case 0.0:
            {
                switch (r)
                {
                    case < 0.0:
                    {
                        b = (ah + ak) switch
                        {
                            < 0.0 => 2.0 * (gh + gk) - 1.0,
                            _ => b
                        };

                        break;
                    }
                    default:
                    {
                        b = (ah - ak) switch
                        {
                            < 0.0 => 2.0 * gk,
                            _ => 2.0 * gh
                        };

                        break;
                    }
                }

                b = Math.Max(b, 0.0);
                b = Math.Min(b, 1.0);
                return b;
            }
        }

        double sqr = Math.Sqrt(rr);

        switch (BivnorData.idig)
        {
            case 15:
                con = twopi * 1.0E-15 / 2.0;
                break;
            default:
            {
                con = twopi / 2.0;
                int i;
                for (i = 1; i <= BivnorData.idig; i++)
                {
                    con /= 10.0;
                }

                break;
            }
        }

        switch (ah)
        {
            //
            //  (0,0)
            //
            case 0.0 when ak == 0.0:
                b = 0.25 + Math.Asin(r) / twopi;
                b = Math.Max(b, 0.0);
                b = Math.Min(b, 1.0);
                return b;
            //
            //  (0,nonzero)
            //
            case 0.0 when ak != 0.0:
                b = gk;
                wh = -ak;
                wk = (ah / ak - r) / sqr;
                gw = 2.0 * gk;
                is_ = 1;
                break;
            //
            default:
            {
                if (ah != 0.0 && ak == 0.0)
                {
                    b = gh;
                    wh = -ah;
                    wk = (ak / ah - r) / sqr;
                    gw = 2.0 * gh;
                    is_ = -1;
                }
                //
                //  (nonzero,nonzero)
                //
                else if (ah != 0.0 && ak != 0.0)
                {
                    b = gh + gk;
                    switch (ah * ak)
                    {
                        case < 0.0:
                            b -= 0.5;
                            break;
                    }

                    wh = -ah;
                    wk = (ak / ah - r) / sqr;
                    gw = 2.0 * gh;
                    is_ = -1;
                }

                break;
            }
        }

        for (;;)
        {
            double sgn = -1.0;

            if (wk != 0.0)
            {
                double t = 0.0;
                switch (Math.Abs(wk))
                {
                    case 1.0:
                        t = wk * gw * (1.0 - gw) / 2.0;
                        b += sgn * t;
                        break;
                    default:
                    {
                        switch (Math.Abs(wk))
                        {
                            case > 1.0:
                            {
                                sgn = -sgn;
                                wh *= wk;
                                double g2 = Gauss.gauss(wh);
                                wk = 1.0 / wk;

                                switch (wk)
                                {
                                    case < 0.0:
                                        b += 0.5;
                                        break;
                                }

                                b = b - (gw + g2) / 2.0 + gw * g2;
                                break;
                            }
                        }

                        double h2 = wh * wh;
                        double a2 = wk * wk;
                        double h4 = h2 / 2.0;
                        double ex = Math.Exp(-h4);
                        double w2 = h4 * ex;
                        double ap = 1.0;
                        double s2 = ap - ex;
                        double sp = ap;
                        double s1 = 0.0;
                        double sn = s1;
                        double conex = Math.Abs(con / wk);

                        for (;;)
                        {
                            double cn = ap * s2 / (sn + sp);
                            s1 += cn;

                            if (Math.Abs(cn) <= conex)
                            {
                                break;
                            }

                            sn = sp;
                            sp += 1.0;
                            s2 -= w2;
                            w2 = w2 * h4 / sp;
                            ap = -ap * a2;
                        }

                        t = (Math.Atan(wk) - wk * s1) / twopi;
                        b += sgn * t;
                        break;
                    }
                }
            }

            if (0 <= is_)
            {
                break;
            }

            if (ak == 0.0)
            {
                break;
            }

            wh = -ak;
            wk = (ah / ak - r) / sqr;
            gw = 2.0 * gk;
            is_ = 1;
        }

        b = Math.Max(b, 0.0);
        b = Math.Min(b, 1.0);

        return b;
    }
}