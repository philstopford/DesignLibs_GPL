using System;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double r8_hyper_2f1(double a, double b, double c, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_HYPER_2F1 evaluates the hypergeometric function F(A,B,C,X).
        //
        //  Discussion:
        //
        //    A minor bug was corrected.  The HW variable, used in several places as
        //    the "old" value of a quantity being iteratively improved, was not
        //    being initialized.  JVB, 11 February 2008.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 July 2009
        //
        //  Author:
        //
        //    Original FORTRAN77 original by Shanjie Zhang, Jianming Jin.
        //    C++ version by John Burkardt.
        //
        //    The original FORTRAN77 version of this routine is copyrighted by
        //    Shanjie Zhang and Jianming Jin.  However, they give permission to
        //    incorporate this routine into a user program provided that the copyright
        //    is acknowledged.
        //
        //  Reference:
        //
        //    Shanjie Zhang, Jianming Jin,
        //    Computation of Special Functions,
        //    Wiley, 1996,
        //    ISBN: 0-471-11963-6,
        //    LC: QA351.C45
        //
        //  Parameters:
        //
        //    Input, double A, B, C, X, the arguments of the function.
        //    C must not be equal to a nonpositive integer.
        //    X < 1.
        //
        //    Output, double R8_HYPER_2F1, the value of the function.
        //
    {
        double c0;
        const double el = 0.5772156649015329;
        double gc;
        double gca;
        double gcab;
        double gcb;
        double hf = 0;
        double hw;
        int k = 0;
        int nm = 0;
        double r;

        bool l0 = Math.Abs(c - (int) c) <= typeMethods.r8_epsilon() && c < 0.0;
        bool l1 = 1.0 - x < 1.0E-15 && c - a - b <= 0.0;
        bool l2 = Math.Abs(a - (int) a) <= typeMethods.r8_epsilon() && a < 0.0;
        bool l3 = Math.Abs(b - (int) b) <= typeMethods.r8_epsilon() && b < 0.0;
        bool l4 = Math.Abs(c - a - (int) (c - a)) <= typeMethods.r8_epsilon() && c - a <= 0.0;
        bool l5 = Math.Abs(c - b - (int) (c - b)) <= typeMethods.r8_epsilon() && c - b <= 0.0;

        switch (l0)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("R8_HYPER_2F1 - Fatal error!");
                Console.WriteLine("  The hypergeometric series is divergent.");
                Console.WriteLine("  C is integral and negative.");
                Console.WriteLine("  C = " + c + "");
                return 1;
        }

        switch (l1)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("R8_HYPER_2F1 - Fatal error!");
                Console.WriteLine("  The hypergeometric series is divergent.");
                Console.WriteLine("  1 - X < 0, C - A - B <= 0");
                Console.WriteLine("  A = " + a + "");
                Console.WriteLine("  B = " + b + "");
                Console.WriteLine("  C = " + c + "");
                Console.WriteLine("  X = " + x + "");
                return 1;
        }

        double eps = x switch
        {
            > 0.95 => 1.0E-08,
            _ => 1.0E-15
        };

        if (x == 0.0 || a == 0.0 || b == 0.0)
        {
            hf = 1.0;
            return hf;
        }

        if (Math.Abs(1.0 - x - eps) <= typeMethods.r8_epsilon() && 0.0 < c - a - b)
        {
            gc = Helpers.Gamma(c);
            gcab = Helpers.Gamma(c - a - b);
            gca = Helpers.Gamma(c - a);
            gcb = Helpers.Gamma(c - b);
            hf = gc * gcab / (gca * gcb);
            return hf;
        }
        if (1.0 + x <= eps && Math.Abs(c - a + b - 1.0) <= eps)
        {
            double g0 = Math.Sqrt(Math.PI) * Math.Pow(2.0, -a);
            double g1 = Helpers.Gamma(c);
            double g2 = Helpers.Gamma(1.0 + a / 2.0 - b);
            double g3 = Helpers.Gamma(0.5 + 0.5 * a);
            hf = g0 * g1 / (g2 * g3);
            return hf;
        }
        if (l2 || l3)
        {
            nm = l3 switch
            {
                true => (int) Math.Abs(b),
                _ => l2 switch
                {
                    true => (int) Math.Abs(a),
                    _ => nm
                }
            };

            hf = 1.0;
            r = 1.0;

            for (k = 1; k <= nm; k++)
            {
                r = r * (a + k - 1.0) * (b + k - 1.0)
                    / (k * (c + k - 1.0)) * x;
                hf += r;
            }

            return hf;
        }
        if (l4 || l5)
        {
            nm = l5 switch
            {
                true => (int) Math.Abs(c - b),
                _ => l4 switch
                {
                    true => (int) Math.Abs(c - a),
                    _ => nm
                }
            };

            hf = 1.0;
            r = 1.0;
            for (k = 1; k <= nm; k++)
            {
                r = r * (c - a + k - 1.0) * (c - b + k - 1.0)
                    / (k * (c + k - 1.0)) * x;
                hf += r;
            }

            hf = Math.Pow(1.0 - x, c - a - b) * hf;
            return hf;
        }

        double aa = a;
        double bb = b;
        double x1 = x;

        switch (x)
        {
            case < 0.0:
            {
                x /= x - 1.0;
                if (a < c && b < a && 0.0 < b)
                {
                    a = bb;
                    b = aa;
                }

                b = c - b;
                break;
            }
        }

        switch (x)
        {
            case >= 0.75:
            {
                double gm = 0.0;

                double c1;
                double gb;
                double ga;
                double r0;
                double r1;
                switch (Math.Abs(c - a - b - (int) (c - a - b)))
                {
                    case < 1.0E-15:
                    {
                        int m = (int) (c - a - b);
                        ga = Helpers.Gamma(a);
                        gb = Helpers.Gamma(b);
                        gc = Helpers.Gamma(c);
                        double gam = Helpers.Gamma(a + m);
                        double gbm = Helpers.Gamma(b + m);

                        double pa = r8_psi(a);
                        double pb = r8_psi(b);

                        if (m != 0)
                        {
                            gm = 1.0;
                        }

                        int j;
                        for (j = 1; j <= Math.Abs(m) - 1; j++)
                        {
                            gm *= j;
                        }

                        double rm = 1.0;
                        for (j = 1; j <= Math.Abs(m); j++)
                        {
                            rm *= j;
                        }

                        double f0 = 1.0;
                        r0 = 1.0;
                        r1 = 1.0;
                        double sp0 = 0.0;
                        double sp = 0.0;

                        double f1;
                        double sm;
                        double rp;
                        switch (m)
                        {
                            case >= 0:
                            {
                                c0 = gm * gc / (gam * gbm);
                                c1 = -gc * Math.Pow(x - 1.0, m) / (ga * gb * rm);

                                for (k = 1; k <= m - 1; k++)
                                {
                                    r0 = r0 * (a + k - 1.0) * (b + k - 1.0)
                                        / (k * (k - m)) * (1.0 - x);
                                    f0 += r0;
                                }

                                for (k = 1; k <= m; k++)
                                {
                                    sp0 = sp0 + 1.0 / (a + k - 1.0) + 1.0 / (b + k - 1.0)
                                          - 1.0 / k;
                                }

                                f1 = pa + pb + sp0 + 2.0 * el + Math.Log(1.0 - x);
                                hw = f1;

                                for (k = 1; k <= 250; k++)
                                {
                                    sp = sp + (1.0 - a) / (k * (a + k - 1.0))
                                            + (1.0 - b) / (k * (b + k - 1.0));

                                    sm = 0.0;
                                    for (j = 1; j <= m; j++)
                                    {
                                        sm = sm + (1.0 - a)
                                                / ((j + k) * (a + j + k - 1.0))
                                                + 1.0 / (b + j + k - 1.0);
                                    }

                                    rp = pa + pb + 2.0 * el + sp + sm + Math.Log(1.0 - x);

                                    r1 = r1 * (a + m + k - 1.0) * (b + m + k - 1.0)
                                        / (k * (m + k)) * (1.0 - x);

                                    f1 += r1 * rp;

                                    if (Math.Abs(f1 - hw) < Math.Abs(f1) * eps)
                                    {
                                        break;
                                    }

                                    hw = f1;
                                }

                                hf = f0 * c0 + f1 * c1;
                                break;
                            }
                            case < 0:
                            {
                                m = -m;
                                c0 = gm * gc / (ga * gb * Math.Pow(1.0 - x, m));
                                c1 = -Math.Pow(-1.0, m) * gc / (gam * gbm * rm);

                                for (k = 1; k <= m - 1; k++)
                                {
                                    r0 = r0 * (a - m + k - 1.0) * (b - m + k - 1.0)
                                        / (k * (k - m)) * (1.0 - x);
                                    f0 += r0;
                                }

                                for (k = 1; k <= m; k++)
                                {
                                    sp0 += 1.0 / k;
                                }

                                f1 = pa + pb - sp0 + 2.0 * el + Math.Log(1.0 - x);
                                hw = f1;

                                for (k = 1; k <= 250; k++)
                                {
                                    sp = sp + (1.0 - a)
                                            / (k * (a + k - 1.0))
                                            + (1.0 - b) / (k * (b + k - 1.0));

                                    sm = 0.0;
                                    for (j = 1; j <= m; j++)
                                    {
                                        sm += 1.0 / (j + k);
                                    }

                                    rp = pa + pb + 2.0 * el + sp - sm + Math.Log(1.0 - x);

                                    r1 = r1 * (a + k - 1.0) * (b + k - 1.0)
                                        / (k * (m + k)) * (1.0 - x);

                                    f1 += r1 * rp;

                                    if (Math.Abs(f1 - hw) < Math.Abs(f1) * eps)
                                    {
                                        break;
                                    }

                                    hw = f1;
                                }

                                hf = f0 * c0 + f1 * c1;
                                break;
                            }
                        }

                        break;
                    }
                    default:
                    {
                        ga = Helpers.Gamma(a);
                        gb = Helpers.Gamma(b);
                        gc = Helpers.Gamma(c);
                        gca = Helpers.Gamma(c - a);
                        gcb = Helpers.Gamma(c - b);
                        gcab = Helpers.Gamma(c - a - b);
                        double gabc = Helpers.Gamma(a + b - c);
                        c0 = gc * gcab / (gca * gcb);
                        c1 = gc * gabc / (ga * gb) * Math.Pow(1.0 - x, c - a - b);
                        hf = 0.0;
                        hw = hf;
                        r0 = c0;
                        r1 = c1;

                        for (k = 1; k <= 250; k++)
                        {
                            r0 = r0 * (a + k - 1.0) * (b + k - 1.0)
                                / (k * (a + b - c + k)) * (1.0 - x);

                            r1 = r1 * (c - a + k - 1.0) * (c - b + k - 1.0)
                                / (k * (c - a - b + k)) * (1.0 - x);

                            hf = hf + r0 + r1;

                            if (Math.Abs(hf - hw) < Math.Abs(hf) * eps)
                            {
                                break;
                            }

                            hw = hf;
                        }

                        hf = hf + c0 + c1;
                        break;
                    }
                }

                break;
            }
            default:
            {
                double a0 = 1.0;

                if (a < c && c < 2.0 * a && b < c && c < 2.0 * b)
                {
                    a0 = Math.Pow(1.0 - x, c - a - b);
                    a = c - a;
                    b = c - b;
                }

                hf = 1.0;
                hw = hf;
                r = 1.0;

                for (k = 1; k <= 250; k++)
                {
                    r = r * (a + k - 1.0) * (b + k - 1.0)
                        / (k * (c + k - 1.0)) * x;

                    hf += r;

                    if (Math.Abs(hf - hw) <= Math.Abs(hf) * eps)
                    {
                        break;
                    }

                    hw = hf;
                }

                hf = a0 * hf;
                break;
            }
        }

        switch (x1)
        {
            case < 0.0:
                x = x1;
                c0 = 1.0 / Math.Pow(1.0 - x, aa);
                hf = c0 * hf;
                break;
        }

        a = aa;
        b = bb;

        switch (k)
        {
            case > 120:
                Console.WriteLine("");
                Console.WriteLine("R8_HYPER_2F1 - Warning!");
                Console.WriteLine("  A large number of iterations were needed.");
                Console.WriteLine("  The accuracy of the results should be checked.");
                break;
        }

        return hf;
    }
}