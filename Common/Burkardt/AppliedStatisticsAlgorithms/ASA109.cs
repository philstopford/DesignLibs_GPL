using System;

namespace Burkardt.AppliedStatistics
{
    public static partial class Algorithms
    {
        public static double xinbta(double p, double q, double beta, double alpha, ref int ifault)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XINBTA computes inverse of the incomplete Beta function.
        //
        //  Discussion:
        //
        //    The accuracy exponent SAE was loosened from -37 to -30, because
        //    the code would not otherwise accept the results of an iteration
        //    with p = 0.3, q = 3.0, alpha = 0.2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 September 2014
        //
        //  Author:
        //
        //    Original FORTRAN77 version by GW Cran, KJ Martin, GE Thomas.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    GW Cran, KJ Martin, GE Thomas,
        //    Remark AS R19 and Algorithm AS 109:
        //    A Remark on Algorithms AS 63: The Incomplete Beta Integral
        //    and AS 64: Inverse of the Incomplete Beta Integeral,
        //    Applied Statistics,
        //    Volume 26, Number 1, 1977, pages 111-114.
        //
        //  Parameters:
        //
        //    Input, double P, Q, the parameters of the incomplete
        //    Beta function.
        //
        //    Input, double BETA, the logarithm of the value of
        //    the complete Beta function.
        //
        //    Input, double ALPHA, the value of the incomplete Beta
        //    function.  0 <= ALPHA <= 1.
        //
        //    Output, int &IFAULT, error flag.
        //    0, no error occurred.
        //    nonzero, an error occurred.
        //
        //    Output, double XINBTA, the argument of the incomplete
        //    Beta function which produces the value ALPHA.
        //
        //  Local Parameters:
        //
        //    Local, double SAE, requests an accuracy of about 10^SAE.
        //
        {
            double a;
            double acu;
            double adj;
            double fpu;
            double g;
            double h;
            int iex;
            bool indx;
            double pp;
            double prev;
            double qq;
            double r;
            double s;
            double sae = -30.0;
            double sq;
            double t;
            double tx;
            double value;
            double w;
            double xin;
            double y;
            double yprev;

            fpu = Math.Pow(10.0, sae);

            ifault = 0;
            value = alpha;
            //
            //  Test for admissibility of parameters.
            //
            if (p <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("XINBTA - Fatal error!");
                Console.WriteLine("  P <= 0.0.");
                ifault = 1;
                return 1;
            }

            if (q <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("XINBTA - Fatal error!");
                Console.WriteLine("  Q <= 0.0.");
                ifault = 1;
                return 1;
            }

            if (alpha < 0.0 || 1.0 < alpha)
            {
                Console.WriteLine("");
                Console.WriteLine("XINBTA - Fatal error!");
                Console.WriteLine("  ALPHA not between 0 and 1.");
                ifault = 2;
                return 1;
            }

            //
            //  If the answer is easy to determine, return immediately.
            //
            if (alpha == 0.0)
            {
                value = 0.0;
                return value;
            }

            if (alpha == 1.0)
            {
                value = 1.0;
                return value;
            }

            //
            //  Change tail if necessary.
            //
            if (0.5 < alpha)
            {
                a = 1.0 - alpha;
                pp = q;
                qq = p;
                indx = true;
            }
            else
            {
                a = alpha;
                pp = p;
                qq = q;
                indx = false;
            }

            //
            //  Calculate the initial approximation.
            //
            r = Math.Sqrt(-Math.Log(a * a));

            y = r - (2.30753 + 0.27061 * r)
                / (1.0 + (0.99229 + 0.04481 * r) * r);

            if (1.0 < pp && 1.0 < qq)
            {
                r = (y * y - 3.0) / 6.0;
                s = 1.0 / (pp + pp - 1.0);
                t = 1.0 / (qq + qq - 1.0);
                h = 2.0 / (s + t);
                w = y * Math.Sqrt(h + r) / h - (t - s)
                    * (r + 5.0 / 6.0 - 2.0 / (3.0 * h));
                value = pp / (pp + qq * Math.Exp(w + w));
            }
            else
            {
                r = qq + qq;
                t = 1.0 / (9.0 * qq);
                t = r * Math.Pow(1.0 - t + y * Math.Sqrt(t), 3);

                if (t <= 0.0)
                {
                    value = 1.0 - Math.Exp((Math.Log((1.0 - a) * qq) + beta) / qq);
                }
                else
                {
                    t = (4.0 * pp + r - 2.0) / t;

                    if (t <= 1.0)
                    {
                        value = Math.Exp((Math.Log(a * pp) + beta) / pp);
                    }
                    else
                    {
                        value = 1.0 - 2.0 / (t + 1.0);
                    }
                }
            }

            //
            //  Solve for X by a modified Newton-Raphson method,
            //  using the function BETAIN.
            //
            r = 1.0 - pp;
            t = 1.0 - qq;
            yprev = 0.0;
            sq = 1.0;
            prev = 1.0;

            if (value < 0.0001)
            {
                value = 0.0001;
            }

            if (0.9999 < value)
            {
                value = 0.9999;
            }

            iex = (int)Math.Max(-5.0 / pp / pp - 1.0 / Math.Pow(a, 0.2) - 13.0, sae);

            acu = Math.Pow(10.0, iex);
            //
            //  Iteration loop.
            //
            for (;;)
            {
                y = betain(value, pp, qq, beta, ref ifault);

                if (ifault != 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("XINBTA - Fatal error!");
                    Console.WriteLine("  BETAIN returned IFAULT = " + ifault + "");
                    ifault = 1;
                    return 1;
                }

                xin = value;
                y = (y - a) * Math.Exp(beta + r * Math.Log(xin) + t * Math.Log(1.0 - xin));

                if (y * yprev <= 0.0)
                {
                    prev = Math.Max(sq, fpu);
                }

                g = 1.0;

                for (;;)
                {
                    for (;;)
                    {
                        adj = g * y;
                        sq = adj * adj;

                        if (sq < prev)
                        {
                            tx = value - adj;

                            if (0.0 <= tx && tx <= 1.0)
                            {
                                break;
                            }
                        }

                        g = g / 3.0;
                    }

                    //
                    //  Check whether the current estimate is acceptable.
                    //  The change "VALUE = TX" was suggested by Ivan Ukhov.
                    //
                    if (prev <= acu || y * y <= acu)
                    {
                        value = tx;
                        if (indx)
                        {
                            value = 1.0 - value;
                        }

                        return value;
                    }

                    if (tx != 0.0 && tx != 1.0)
                    {
                        break;
                    }
                }
            }
        }
    }
}
