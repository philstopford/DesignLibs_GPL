using System;
using Burkardt.Quadrature;
using Burkardt.Weight;

namespace TOMS655Test;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TOMS655_TEST.
        //
        //  Discussion:
        //
        //    TOMS655_TEST tests the TOMS655 library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = 0;
        double alpha = 0;
        double b = 0;
        double beta = 0;
        int kind;
        int nt;

        Console.WriteLine("");
        Console.WriteLine("TOMS655_TEST");
        Console.WriteLine("  Test the TOMS655 library.");

        test01();
        test02();
        test03();
        test04();
        test05();
        test06();
        test07();
        test08();
        test09();
        //
        //  Compute 15 points of an example of each rule.
        //
        for (kind = 1; kind <= 9; kind++)
        {
            nt = 15;
            switch (kind)
            {
                case 8:
                    alpha = 1.0;
                    beta = -alpha - 2 * nt - 2;
                    break;
                default:
                    alpha = 0.0;
                    beta = 0.0;
                    break;
            }

            test10(nt, kind, alpha, beta);
        }

        //
        //  Compute 15 points of an example of each rule using nondefault A, B.
        //
        for (kind = 1; kind <= 9; kind++)
        {
            nt = 15;

            switch (kind)
            {
                case 1:
                    alpha = 0.0;
                    beta = 0.0;
                    a = 0.0;
                    b = 1.0;
                    break;
                case 2:
                    alpha = 0.0;
                    beta = 0.0;
                    a = 0.0;
                    b = 1.0;
                    break;
                case 3:
                    alpha = 1.0;
                    beta = 0.0;
                    a = 0.0;
                    b = 1.0;
                    break;
                case 4:
                    alpha = 1.5;
                    beta = 0.5;
                    a = 0.0;
                    b = 1.0;
                    break;
                case 5:
                    alpha = 1.0;
                    beta = 0.0;
                    a = 1.0;
                    b = 1.0;
                    break;
                case 6:
                    alpha = 1.0;
                    beta = 0.0;
                    a = 0.0;
                    b = 0.5;
                    break;
                case 7:
                    alpha = 1.0;
                    beta = 0.0;
                    a = 0.0;
                    b = 1.0;
                    break;
                case 8:
                    alpha = 1.0;
                    beta = -alpha - 2 * nt - 2;
                    a = 0.0;
                    b = 1.0;
                    break;
                case 9:
                    alpha = 0.0;
                    beta = 0.0;
                    a = 0.0;
                    b = 1.0;
                    break;
            }

            cgqf_test(nt, kind, alpha, beta, a, b);
        }

        wm_test();

        Console.WriteLine("");
        Console.WriteLine("TOMS655_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests CIQFS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 January 2010
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
        //    C++ version by John Burkardt.
        //
    {
        double alpha;
        double beta;
        int i;
        int key;
        int kind;
        int lu;
        int[] mlt;
        int[] ndx;
        int nt;
        int nwts;
        double pi = 3.14159265358979323846264338327950;
        double[] t;

        Console.WriteLine("  ----------------------------------------");
        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Test CIQFS.");
        //
        //  Number of knots.
        //
        nt = 5;
        //
        //  Set the knots in the default interval [-1,+1].
        //
        t = new double[nt];

        for (i = 1; i <= nt; i++)
        {
            t[i - 1] = Math.Cos((2 * i - 1) * pi / (2 * nt));
        }

        //
        //  Set the knot multiplicities.
        //
        mlt = new int[nt];
        for (i = 0; i < nt; i++)
        {
            mlt[i] = 2;
        }

        //
        //  Set the size of the weights array.
        //
        nwts = 0;
        for (i = 0; i < nt; i++)
        {
            nwts += mlt[i];
        }

        //
        //  Because KEY = 1, NDX will be set up for us.
        //
        ndx = new int[nt];
        //
        //  KEY = 1 indicates that the WTS array should hold the weights
        //  in the usual order.
        //
        key = 1;
        //
        //  Request Legendre weight function.
        //
        kind = 1;
        //
        //  ALPHA, BETA not used in Legendre weight function but set anyway.
        //
        alpha = 0.0;
        beta = 0.0;
        //
        //  LU controls printing.
        //  A positive value requests that we compute and print weights, and
        //  conduct a moments check.
        //
        lu = 6;
        //
        //  This call returns the WTS array.
        //
        CIQFS.ciqfs(nt, t, mlt, nwts, ref ndx, key, kind, alpha, beta, lu);
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests CIQFS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 February 2010
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
        //    C++ version by John Burkardt.
        //
    {
        double a;
        double alpha;
        double b;
        double beta;
        int i;
        int key;
        int kind;
        int lu;
        int[] mlt;
        int[] ndx;
        int nt;
        int nwts;
        double[] t;
        double[] wts;

        Console.WriteLine("  ----------------------------------------");
        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Test CIQF, CIQFS, CGQF and CGQFS");
        Console.WriteLine("  with all classical weight functions.");
        //
        //  Try all weight functions.
        //
        for (kind = 1; kind <= 9; kind++)
        {
            //
            //  Number of knots.
            //
            nt = 5;
            //
            //  Set parameters ALPHA and BETA.
            //
            alpha = 0.5;
            if (kind != 8)
            {
                beta = 2.0;
            }
            else
            {
                beta = -16.0;
            }

            //
            //  Set A and B.
            //
            a = -0.5;
            b = 2.0;
            //
            //  Have CGQF compute the knots and weights.
            //
            lu = 6;
            t = new double[nt];
            wts = new double[nt];

            Console.WriteLine("");
            Console.WriteLine("  Knots and weights of Gauss quadrature formula");
            Console.WriteLine("  computed by CGQF.");
            CGQF.cgqf(nt, kind, alpha, beta, a, b, lu, ref t, ref wts);
            //
            //  Now compute the weights for the same knots by CIQF.
            //
            //  Set the knot multiplicities.
            //
            mlt = new int[nt];
            for (i = 0; i < nt; i++)
            {
                mlt[i] = 2;
            }

            //
            //  Set the size of the weights array.
            //
            nwts = 0;
            for (i = 0; i < nt; i++)
            {
                nwts += mlt[i];
            }

            //
            //  We need to deallocate and reallocate WTS because it is now of
            //  dimension NWTS rather than NT.
            //
            wts = new double[nwts];
            //
            //  Because KEY = 1, NDX will be set up for us.
            //
            ndx = new int[nt];
            //
            //  KEY = 1 indicates that the WTS array should hold the weights
            //  in the usual order.
            //
            key = 1;
            //
            //  LU controls printing.
            //  A positive value requests that we compute and print weights, and
            //  conduct a moments check.
            //
            lu = 6;

            Console.WriteLine("");
            Console.WriteLine("  Weights of Gauss quadrature formula computed from the");
            Console.WriteLine("  knots by CIQF.");

            wts = CIQF.ciqf(nt, t, mlt, nwts, ref ndx, key, kind, alpha, beta, a, b, lu);

        }
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests CEIQFS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 January 2010
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
        //    C++ version by John Burkardt.
        //
    {
        double alpha;
        double beta;
        int i;
        int kind;
        int[] mlt;
        int nt;
        double pi = 3.14159265358979323846264338327950;
        double qfsum;
        double qfsx;
        double[] t;

        Console.WriteLine("  ----------------------------------------");
        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  Test CEIQFS.");
        //
        //  Number of knots.
        //
        nt = 5;
        //
        //  Set the knots in the default interval [-1,+1].
        //
        t = new double[nt];

        for (i = 1; i <= nt; i++)
        {
            t[i - 1] = Math.Cos((2 * i - 1) * pi / (2 * nt));
        }

        //
        //  Set the knot multiplicities.
        //
        mlt = new int[nt];
        for (i = 0; i < nt; i++)
        {
            mlt[i] = 2;
        }

        //
        //  Set KIND to the Legendre weight function.
        //
        kind = 1;
        //
        //  ALPHA, BETA not used in Legendre weight function but set anyway.
        //
        alpha = 0.0;
        beta = 0.0;
        //
        //  Call CEIQFS to set up the quadrature formula and evaluate it on F.
        //
        qfsum = CEIQFS.ceiqfs(nt, t, mlt, kind, alpha, beta, f);

        Console.WriteLine("");
        Console.WriteLine("  Integral of sin(x) on -1, 1 by Fejer type rule");
        Console.WriteLine("  with " + nt + " points of multiplicity 2.");
        Console.WriteLine("  Quadrature formula:" + qfsum.ToString("0.################").PadLeft(24) + "");

        qfsx = Math.Cos(-1.0) - Math.Cos(1.0);
        Console.WriteLine("  Exact value       :" + qfsx.ToString("0.################").PadLeft(24) + "");
        Console.WriteLine("  Error             :" + Math.Abs(qfsum - qfsx) + "");
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests CEIQF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 January 2010
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
        //    C++ version by John Burkardt.
        //
    {
        double a;
        double alpha;
        double b;
        double beta;
        int i;
        int kind;
        int[] mlt;
        int nt;
        double pi = 3.14159265358979323846264338327950;
        double qfsum;
        double qfsx;
        double[] t;

        Console.WriteLine("  ----------------------------------------");
        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  Test CEIQF.");
        //
        //  Number of knots.
        //
        nt = 5;
        //
        //  Set the knots in the default interval [-1,+1].
        //
        t = new double[nt];

        for (i = 1; i <= nt; i++)
        {
            t[i - 1] = Math.Cos((2 * i - 1) * pi / (2 * nt));
        }

        //
        //  Set the knot multiplicities.
        //
        mlt = new int[nt];
        for (i = 0; i < nt; i++)
        {
            mlt[i] = 2;
        }

        //
        //  Set KIND to the Legendre weight function.
        //
        kind = 1;
        //
        //  ALPHA, BETA not used in Legendre weight function but set anyway.
        //
        alpha = 0.0;
        beta = 0.0;
        //
        //  Set nonstandard interval A, B.
        //
        a = -0.5;
        b = 2.0;
        //
        //  Shift knots from [-1,1] to [A,B].
        //
        for (i = 0; i < nt; i++)
        {
            t[i] = ((b - a) * t[i] + (a + b)) / 2.0;
        }

        //
        //  Call CEIQF to set up the quadrature formula and evaluate it on F.
        //
        qfsum = CEIQF.ceiqf(nt, t, mlt, kind, alpha, beta, a, b, f);

        Console.WriteLine("");
        Console.WriteLine("  Integral of sin(x) from " + a + " to " + b + "");
        Console.WriteLine("  by Fejer type rule with " + nt + " points");
        Console.WriteLine("  of multiplicity 2.");
        Console.WriteLine("  Quadrature formula:" + qfsum.ToString("0.################").PadLeft(24) + "");

        qfsx = Math.Cos(a) - Math.Cos(b);
        Console.WriteLine("  Exact value       :" + qfsx.ToString("0.################").PadLeft(24) + "");
        Console.WriteLine("  Error             :" + Math.Abs(qfsum - qfsx) + "");

    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests CLIQFS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 January 2010
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
        //    C++ version by John Burkardt.
        //
    {
        double alpha;
        double beta;
        int i;
        int kind;
        int lu;
        int nt;
        double pi = 3.14159265358979323846264338327950;
        double[] t;

        Console.WriteLine("  ----------------------------------------");
        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  Test CLIQFS.");
        //
        //  Number of knots.
        //
        nt = 5;
        //
        //  Set the knots in the default interval [-1,+1].
        //
        t = new double[nt];

        for (i = 1; i <= nt; i++)
        {
            t[i - 1] = Math.Cos((2 * i - 1) * pi / (2 * nt));
        }

        //
        //  Request Legendre weight function.
        //
        kind = 1;
        //
        //  ALPHA, BETA not used in Legendre weight function but set anyway.
        //
        alpha = 0.0;
        beta = 0.0;
        //
        //  LU controls printing.
        //  A positive value requests that we compute and print weights, and
        //  conduct a moments check.
        //
        lu = 6;
        //
        //  This call returns the WTS array.
        //
        CLIQFS.cliqfs(nt, t, kind, alpha, beta, lu);
    }

    private static void test06()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 tests CLIQF and EIQFS..
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 January 2010
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
        //    C++ version by John Burkardt.
        //
    {
        double a;
        double alpha;
        double b;
        double beta;
        int i;
        int kind;
        int lu;
        int nt;
        double pi = 3.14159265358979323846264338327950;
        double qfsum;
        double qfsx;
        double[] t;
        double[] wts;

        Console.WriteLine("  ----------------------------------------");
        Console.WriteLine("");
        Console.WriteLine("TEST06");
        Console.WriteLine("  Test CLIQF and EIQFS.");
        //
        //  Number of knots.
        //
        nt = 5;
        //
        //  Set the knots in the default interval [-1,+1].
        //
        t = new double[nt];

        for (i = 1; i <= nt; i++)
        {
            t[i - 1] = Math.Cos((2 * i - 1) * pi / (2 * nt));
        }

        //
        //  Set KIND to the Legendre weight function.
        //
        kind = 1;
        //
        //  ALPHA, BETA not used in Legendre weight function but set anyway.
        //
        alpha = 0.0;
        beta = 0.0;
        //
        //  Set nonstandard interval A, B.
        //
        a = -0.5;
        b = 2.0;
        //
        //  Shift knots from [-1,1] to [A,B].
        //
        for (i = 0; i < nt; i++)
        {
            t[i] = ((b - a) * t[i] + (a + b)) / 2.0;
        }

        //
        //  LU controls printout.
        //
        lu = 6;
        //
        //  Call CLIQF to set up the quadrature formula.
        //
        wts = CLIQF.cliqf(nt, t, kind, alpha, beta, a, b, lu);
        //
        //  Call EIQFS to evaluate the quadrature formula.
        //
        qfsum = EIQFS.eiqfs(nt, t, wts, f);

        Console.WriteLine("");
        Console.WriteLine("  Integral of sin(x) from " + a + " to " + b + "");
        Console.WriteLine("  by Fejer type rule with " + nt + " points");
        Console.WriteLine("  of multiplicity 1.");
        Console.WriteLine("  Quadrature formula:" + qfsum.ToString("0.################").PadLeft(24) + "");

        qfsx = Math.Cos(a) - Math.Cos(b);
        Console.WriteLine("  Exact value       :" + qfsx.ToString("0.################").PadLeft(24) + "");
        Console.WriteLine("  Error             :" + Math.Abs(qfsum - qfsx) + "");
    }

    private static void test07()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST07 tests CEGQF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 January 2010
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
        //    C++ version by John Burkardt.
        //
    {
        double a;
        double alpha;
        double b;
        double beta;
        int kind;
        int nt;
        double qfsum;
        double qfsx;

        Console.WriteLine("  ----------------------------------------");
        Console.WriteLine("");
        Console.WriteLine("TEST07");
        Console.WriteLine("  Test CEGQF.");
        //
        //  Number of knots.
        //
        nt = 12;
        //
        //  Request exponential weight function.
        //
        kind = 7;
        //
        //  Set ALPHA and BETA.
        //
        alpha = 1.0;
        beta = 0.0;
        //
        //  Set interval [A,B].
        //
        a = -0.5;
        b = 2.0;
        //
        //  Call CEGQF to compute and evaluate the Gauss quadrature formula.
        //
        qfsum = CEGQF.cegqf(nt, kind, alpha, beta, a, b, f);

        Console.WriteLine("");
        Console.WriteLine("  Integral of x*sin(x) from " + a + " to " + b + "");
        Console.WriteLine("  by Gauss-exponential rule with " + nt + " points");
        Console.WriteLine("  Quadrature formula:" + qfsum.ToString("0.################").PadLeft(24) + "");

        qfsx = (b - a) * 0.5 * (Math.Cos(a) - Math.Cos(b))
            + Math.Sin(b) + Math.Sin(a) - 2.0 * Math.Sin((a + b) / 2.0);

        Console.WriteLine("  Exact value       :" + qfsx.ToString("0.################").PadLeft(24) + "");
        Console.WriteLine("  Error             :" + Math.Abs(qfsum - qfsx) + "");
    }

    private static void test08()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST08 tests CEGQFS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 January 2010
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
        //    C++ version by John Burkardt.
        //
    {
        double alpha;
        double beta;
        int kind;
        int nt;
        double qfsum;
        double qfsx;

        Console.WriteLine("  ----------------------------------------");
        Console.WriteLine("");
        Console.WriteLine("TEST08");
        Console.WriteLine("  Test CEGQFS.");
        //
        //  Number of knots.
        //
        nt = 12;
        //
        //  Request exponential weight function.
        //
        kind = 7;
        //
        //  Set ALPHA and BETA.
        //
        alpha = 1.0;
        beta = 0.0;
        //
        //  Call CEGQFS to compute and evaluate the Gauss quadrature formula.
        //
        qfsum = CEGQFS.cegqfs(nt, kind, alpha, beta, f);

        Console.WriteLine("");
        Console.WriteLine("  Integral of x*sin(x) from -1 to +1");
        Console.WriteLine("  by Gauss-exponential rule with " + nt + " points.");
        Console.WriteLine("  Quadrature formula: " + qfsum.ToString("0.################").PadLeft(24) + "");

        qfsx = Math.Cos(-1.0) - Math.Cos(+1.0);

        Console.WriteLine("  Exact value       :" + qfsx.ToString("0.################").PadLeft(24) + "");
        Console.WriteLine("  Error             :" + Math.Abs(qfsum - qfsx) + "");
    }

    private static void test09()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST09 calls CGQFS to compute and print generalized Gauss-Hermite rules.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 January 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double alpha;
        double beta;
        int io;
        int kind;
        int nt;
        double[] t;
        double[] wts;

        Console.WriteLine("");
        Console.WriteLine("TEST09");
        Console.WriteLine("  Call CGQFS to compute generalized Hermite rules.");

        nt = 15;
        kind = 6;
        alpha = 1.0;
        beta = 0.0;
        io = -6;
        t = new double[nt];
        wts = new double[nt];

        Console.WriteLine("");
        Console.WriteLine("  NT = " + nt + "");
        Console.WriteLine("  ALPHA = " + alpha + "");

        CGQFS.cgqfs(nt, kind, alpha, beta, io, ref t, ref wts);
    }

    private static void test10(int nt, int kind, double alpha, double beta)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST10 calls CDGQF to compute a quadrature formula.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        double[] t;
        double[] wts;

        Console.WriteLine("");
        Console.WriteLine("TEST10");
        Console.WriteLine("  Call CDGQF to compute a quadrature formula.");
        Console.WriteLine("");
        Console.WriteLine("  KIND = " + kind + "");
        Console.WriteLine("  ALPHA = " + alpha + "");
        Console.WriteLine("  BETA  = " + beta + "");

        t = new double[nt];
        wts = new double[nt];

        CDGQF.cdgqf(nt, kind, alpha, beta, ref t, ref wts);

        Console.WriteLine("");
        Console.WriteLine(" Index     Abscissas                 Weights");
        Console.WriteLine("");
        for (i = 0; i < nt; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + t[i].ToString("0.################").PadLeft(24)
                                   + "  " + wts[i].ToString("0.################").PadLeft(24) + "");
        }
    }

    private static void cgqf_test(int nt, int kind, double alpha, double beta, double a, double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CGQF_TEST calls CGQF to compute a quadrature formula.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 February 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int lo;
        double[] t;
        double[] wts;

        Console.WriteLine("");
        Console.WriteLine("CGQF_TEST");
        Console.WriteLine("  CGQF computes a quadrature formula with nondefault");
        Console.WriteLine("  values of parameters A and B.");
        Console.WriteLine("");
        Console.WriteLine("  KIND =  " + kind + "");
        Console.WriteLine("  ALPHA = " + alpha + "");
        Console.WriteLine("  BETA  = " + beta + "");
        Console.WriteLine("  A =     " + a + "");
        Console.WriteLine("  B  =    " + b + "");

        lo = 0;
        t = new double[nt];
        wts = new double[nt];

        CGQF.cgqf(nt, kind, alpha, beta, a, b, lo, ref t, ref wts);

        Console.WriteLine("");
        Console.WriteLine(" Index     Abscissas                 Weights");
        Console.WriteLine("");
        for (i = 0; i < nt; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + t[i].ToString("0.################").PadLeft(24)
                                   + "  " + wts[i].ToString("0.################").PadLeft(24) + "");
        }
    }

    private static void wm_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WM_TEST calls WM_TESTER with various parameter values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 November 2015
        //
        //  Author:
        //
        //    John Burkardt.
        //
    {
        double alpha;
        double beta;
        int kind;
        int m;

        m = 5;
        kind = 1;
        alpha = 0.0;
        beta = 0.0;
        wm_tester(m, kind, alpha, beta);

        m = 5;
        kind = 2;
        alpha = 0.0;
        beta = 0.0;
        wm_tester(m, kind, alpha, beta);

        m = 5;
        kind = 3;
        alpha = 0.5;
        beta = 0.0;
        wm_tester(m, kind, alpha, beta);

        m = 5;
        kind = 4;
        alpha = 0.25;
        beta = 0.75;
        wm_tester(m, kind, alpha, beta);

        m = 5;
        kind = 5;
        alpha = 2.0;
        beta = 0.0;
        wm_tester(m, kind, alpha, beta);

        m = 5;
        kind = 6;
        alpha = 1.0;
        beta = 0.0;
        wm_tester(m, kind, alpha, beta);

        m = 5;
        kind = 7;
        alpha = 2.0;
        beta = 0.0;
        wm_tester(m, kind, alpha, beta);

        m = 5;
        kind = 8;
        alpha = -0.5;
        beta = -6.0;
        wm_tester(m, kind, alpha, beta);

        m = 5;
        kind = 9;
        alpha = 0.0;
        beta = 0.0;
        wm_tester(m, kind, alpha, beta);
    }

    private static void wm_tester(int m, int kind, double alpha, double beta)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WM_TESTER tests WM.
        //
        //  Discussion:
        //
        //    Moment(K) = Integral ( A <= X <= B ) X^(K-1) * W(X) dx
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 November 2015
        //
        //  Author:
        //
        //    John Burkardt.
        //
        //  Reference:
        //
        //    Sylvan Elhay, Jaroslav Kautsky,
        //    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
        //    Interpolatory Quadrature,
        //    ACM Transactions on Mathematical Software,
        //    Volume 13, Number 4, December 1987, pages 399-415.
        //
        //  Parameters:
        //
        //    Input, int M, the number of moments to evaluate.
        //
        //    Input, int KIND, the rule.
        //    1, Legendre,             (a,b)       1.0
        //    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
        //    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
        //    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
        //    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
        //    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
        //    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
        //    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
        //    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
        //
        //    Input, double ALPHA, the value of Alpha, if needed.
        //
        //    Input, double BETA, the value of Beta, if needed.
        //
    {
        int i;
        double[] w;

        w = WM.wm(m, kind, alpha, beta);

        Console.WriteLine("");
        Console.WriteLine("WM_TESTER:");
        Console.WriteLine("  WM_TEST computes moments for rule " + kind + "");
        Console.WriteLine("  with ALPHA = " + alpha + ", BETA = " + beta + "");
        Console.WriteLine("");
        Console.WriteLine("  Order          Moment");
        Console.WriteLine("");
        for (i = 0; i < m; i++)
        {
            Console.WriteLine("   "
                              + "  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                              + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static double f(double x, int i)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F returns values of the integrand or its derivatives.
        //
        //  Discussion:
        //
        //    This function is an example of an integrand function.
        //
        //    The package can generate quadrature formulas that use derivative
        //    information as well as function values.  Therefore, this routine is
        //    set up to provide derivatives of any order as well as the function
        //    value.  In an actual application, the highest derivative needed
        //    is of order one less than the highest knot multiplicity.
        //
        //    In other words, in the usual case where knots are not repeated,
        //    this routine only needs to return function values, not any derivatives.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 January 2010
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
        //    C++ version by John Burkardt.
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Input, int I, the order of the derivative of F to
        //    be evaluated.
        //
        //    Output, double F, the value of the I-th derivative of F at X.
        //
    {
        int l;
        double value = 0;

        l = i % 4;

        value = l switch
        {
            0 => Math.Sin(x),
            1 => Math.Cos(x),
            2 => -Math.Sin(x),
            3 => -Math.Cos(x),
            _ => value
        };

        return value;
    }
}