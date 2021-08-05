using System;

namespace Burkardt.CDFLib
{
    public static class bivariatenormal
    {
        public class BivnorData
        {
            public int idig = 15;

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
            double a2;
            double ap;
            double b;
            double cn;
            double con;
            double conex;
            double ex;
            double g2;
            double gh;
            double gk;
            double gw = 0;
            double h2;
            double h4;
            int i;
            int is_ = 0;
            double rr;
            double s1;
            double s2;
            double sgn;
            double sn;
            double sp;
            double sqr;
            double t;
            double twopi = Math.PI * 2;
            double w2;
            double wh = 0;
            double wk = 0;

            b = 0.0;

            gh = Gauss.gauss(-ah) / 2.0;
            gk = Gauss.gauss(-ak) / 2.0;

            if (r == 0.0)
            {
                b = 4.00 * gh * gk;
                b = Math.Max(b, 0.0);
                b = Math.Min(b, 1.0);
                return b;
            }

            rr = (1.0 + r) * (1.0 - r);

            if (rr < 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("BIVNOR - Fatal error!");
                Console.WriteLine("  1 < |R|.");
            }

            if (rr == 0.0)
            {
                if (r < 0.0)
                {
                    if (ah + ak < 0.0)
                    {
                        b = 2.0 * (gh + gk) - 1.0;
                    }
                }
                else
                {
                    if (ah - ak < 0.0)
                    {
                        b = 2.0 * gk;
                    }
                    else
                    {
                        b = 2.0 * gh;
                    }
                }

                b = Math.Max(b, 0.0);
                b = Math.Min(b, 1.0);
                return b;
            }

            sqr = Math.Sqrt(rr);

            if (data.idig == 15)
            {
                con = twopi * 1.0E-15 / 2.0;
            }
            else
            {
                con = twopi / 2.0;
                for (i = 1; i <= data.idig; i++)
                {
                    con = con / 10.0;
                }
            }

            //
            //  (0,0)
            //
            if (ah == 0.0 && ak == 0.0)
            {
                b = 0.25 + Math.Asin(r) / twopi;
                b = Math.Max(b, 0.0);
                b = Math.Min(b, 1.0);
                return b;
            }

            //
            //  (0,nonzero)
            //
            if (ah == 0.0 && ak != 0.0)
            {
                b = gk;
                wh = -ak;
                wk = (ah / ak - r) / sqr;
                gw = 2.0 * gk;
                is_ = 1;
            }
            //
            //  (nonzero,0)
            //
            else if (ah != 0.0 && ak == 0.0)
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
                if (ah * ak < 0.0)
                {
                    b = b - 0.5;
                }

                wh = -ah;
                wk = (ak / ah - r) / sqr;
                gw = 2.0 * gh;
                is_ = -1;
            }

            for (;;)
            {
                sgn = -1.0;
                t = 0.0;

                if (wk != 0.0)
                {
                    if (Math.Abs(wk) == 1.0)
                    {
                        t = wk * gw * (1.0 - gw) / 2.0;
                        b = b + sgn * t;
                    }
                    else
                    {
                        if (1.0 < Math.Abs(wk))
                        {
                            sgn = -sgn;
                            wh = wh * wk;
                            g2 = Gauss.gauss(wh);
                            wk = 1.0 / wk;

                            if (wk < 0.0)
                            {
                                b = b + 0.5;
                            }

                            b = b - (gw + g2) / 2.0 + gw * g2;
                        }

                        h2 = wh * wh;
                        a2 = wk * wk;
                        h4 = h2 / 2.0;
                        ex = Math.Exp(-h4);
                        w2 = h4 * ex;
                        ap = 1.0;
                        s2 = ap - ex;
                        sp = ap;
                        s1 = 0.0;
                        sn = s1;
                        conex = Math.Abs(con / wk);

                        for (;;)
                        {
                            cn = ap * s2 / (sn + sp);
                            s1 = s1 + cn;

                            if (Math.Abs(cn) <= conex)
                            {
                                break;
                            }

                            sn = sp;
                            sp = sp + 1.0;
                            s2 = s2 - w2;
                            w2 = w2 * h4 / sp;
                            ap = -ap * a2;
                        }

                        t = (Math.Atan(wk) - wk * s1) / twopi;
                        b = b + sgn * t;
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
}