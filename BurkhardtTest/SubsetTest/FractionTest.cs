using System;
using Burkardt.Function;
using Burkardt.RationalNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace SubsetTestNS;

public static class FractionTest
{
    public static void cfrac_to_rat_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CFRAC_TO_RAT_TEST tests CFRAC_TO_RAT.
        //
        //  Discussion:
        //
        //    Compute the continued fraction form of 4096/15625.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 October 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int M = 10;

        int[] a = new int[M];
        int bot = 15625;
        bool error = false;
        int i;
        int n = 0;
        int[] p = new int[M];
        int[] q = new int[M];
        int top = 4096;

        Console.WriteLine("");
        Console.WriteLine("CFRAC_TO_RAT_TEST");
        Console.WriteLine("  CFRAC_TO_RAT continued fraction => fraction.");
        Console.WriteLine("");
        Console.WriteLine("  Regular fraction is " + top + "/" + bot + "");

        Rational.rat_to_cfrac(top, bot, M, ref n, ref a, ref error);

        typeMethods.i4vec1_print(n, a, "  Continued fraction coefficients:");

        Fraction.cfrac_to_rat(n, a, ref p, ref q);

        Console.WriteLine("");
        Console.WriteLine("  The continued fraction convergents.");
        Console.WriteLine("  The last row contains the value of the continued");
        Console.WriteLine("  fraction, written as a common fraction.");
        Console.WriteLine("");
        Console.WriteLine("  I, P(I), Q(I), P(I)/Q(I)");
        Console.WriteLine("");

        for (i = 0; i < n; i++)
        {
            Console.WriteLine(i.ToString().PadLeft(3) + "  "
                                                      + p[i].ToString().PadLeft(6) + "  "
                                                      + q[i].ToString().PadLeft(6) + "  "
                                                      + (p[i] / (double)q[i]).ToString().PadLeft(14) + "");
        }

    }

    public static void cfrac_to_rfrac_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CFRAC_TO_RFRAC_TEST tests CFRAC_TO_RFRAC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 October 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int MAXM = 10;

        bool error = false;
        double[] g = new double[2 * MAXM];
        double[] h = new double[2 * MAXM];
        int i;
        int m;
        double[] p = new double[MAXM];
        double[] q = new double[MAXM + 1];
        string cout;

        m = 3;

        p[0] = 1.0;
        p[1] = 1.0;
        p[2] = 2.0;

        q[0] = 1.0;
        q[1] = 3.0;
        q[2] = 1.0;
        q[3] = 1.0;

        Console.WriteLine("");
        Console.WriteLine("CFRAC_TO_RFRAC_TEST");
        Console.WriteLine("  CFRAC_TO_RFRAC: continued fraction to ratio;");

        Console.WriteLine("");
        Console.WriteLine("  Rational polynomial fraction coefficients:");
        Console.WriteLine("");

        cout = "  P:  ";
        for (i = 0; i < m; i++)
        {
            cout += p[i].ToString().PadLeft(12);
        }

        Console.WriteLine(cout);

        cout = "  Q:  ";
        for (i = 0; i < m + 1; i++)
        {
            cout += q[i].ToString().PadLeft(12);
        }

        Console.WriteLine(cout);

        Fraction.rfrac_to_cfrac(m, p, q, ref h, ref error);

        typeMethods.r8vec_print(2 * m, h, "  Continued fraction coefficients:");

        for (i = 0; i < 2 * m; i++)
        {
            g[i] = 1.0;
        }

        Fraction.cfrac_to_rfrac(2 * m, g, h, ref p, ref q);

        Console.WriteLine("");
        Console.WriteLine("  Recovered rational polynomial:");
        Console.WriteLine("");

        cout = "  P:  ";
        for (i = 0; i < m; i++)
        {
            cout += p[i].ToString().PadLeft(12);
        }

        Console.WriteLine(cout);

        cout = "  Q:  ";
        for (i = 0; i < m + 1; i++)
        {
            cout += q[i].ToString().PadLeft(12);
        }

        Console.WriteLine(cout);
    }

    public static void jfrac_to_rfrac_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    JFRAC_TO_RFRAC_TEST tests JFRAC_TO_RFRAC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int MAXM = 10;

        int i;
        int m;
        double[] p = new double[MAXM];
        double[] q = new double[MAXM];
        double[] r = new double[MAXM];
        double[] s = new double[MAXM];
        int seed;
        string cout = "";

        //
        //  Generate the data, but force Q(M+1) to be 1.  
        //  That will make it easier to see that the two operations are inverses
        //  of each other.  JFRAC_TO_RFRAC is free to scale its output, and chooses
        //  a scaling in which Q(M+1) is 1.
        //
        seed = 123456789;
        m = 6;
        UniformRNG.r8vec_uniform_01(m, ref seed, ref p);
        UniformRNG.r8vec_uniform_01(m + 1, ref seed, ref q);

        for (i = 0; i < m; i++)
        {
            q[i] /= q[m];
        }

        q[m] = 1.0;

        Console.WriteLine("");
        Console.WriteLine("JFRAC_TO_RFRAC_TEST");
        Console.WriteLine("  JFRAC_TO_RFRAC converts a J fraction");
        Console.WriteLine("  to a rational polynomial fraction.");
        Console.WriteLine("");
        Console.WriteLine("  The original rational polynomial coefficients:");
        Console.WriteLine("");

        for (i = 0; i < m; i++)
        {
            cout += p[i].ToString().PadLeft(14) + "  ";
        }

        Console.WriteLine(cout);

        cout = "";

        for (i = 0; i < m + 1; i++)
        {
            cout += q[i].ToString().PadLeft(14) + "  ";
        }

        Console.WriteLine(cout);

        cout = "";

        Fraction.rfrac_to_jfrac(m, p, q, ref r, ref s);

        Console.WriteLine("");
        Console.WriteLine("  The J fraction coefficients:");
        Console.WriteLine("");

        for (i = 0; i < m; i++)
        {
            cout += r[i].ToString().PadLeft(14) + "  ";
        }

        Console.WriteLine(cout);

        cout = "";

        for (i = 0; i < m; i++)
        {
            cout += s[i].ToString().PadLeft(14) + "  ";
        }

        Console.WriteLine(cout);

        cout = "";

        Fraction.jfrac_to_rfrac(m, r, s, ref p, ref q);

        Console.WriteLine("");
        Console.WriteLine("  The recovered rational polynomial:");
        Console.WriteLine("");

        for (i = 0; i < m; i++)
        {
            cout += p[i].ToString().PadLeft(14) + "  ";
        }

        Console.WriteLine(cout);

        cout = "";

        for (i = 0; i < m + 1; i++)
        {
            cout += q[i].ToString().PadLeft(14) + "  ";
        }

        Console.WriteLine(cout);

    }

    public static void rfrac_to_cfrac_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RFRAC_TO_CFRAC_TEST tests RFRAC_TO_CFRAC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 October 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int MAXM = 10;

        bool error = false;
        double[] g = new double[2 * MAXM];
        double[] h = new double[2 * MAXM];
        int i;
        int m;
        double[] p = new double[MAXM];
        double[] q = new double[MAXM + 1];

        m = 3;

        p[0] = 1.0;
        p[1] = 1.0;
        p[2] = 2.0;

        q[0] = 1.0;
        q[1] = 3.0;
        q[2] = 1.0;
        q[3] = 1.0;

        Console.WriteLine("");
        Console.WriteLine("RFRAC_TO_CFRAC_TEST");
        Console.WriteLine("  RFRAC_TO_CFRAC: ratio to continued fration.");

        Console.WriteLine("");
        Console.WriteLine("  Rational polynomial fraction coefficients:");
        Console.WriteLine("");

        string cout = "  P:  ";
        for (i = 0; i < m; i++)
        {
            cout += p[i].ToString().PadLeft(12);
        }

        Console.WriteLine(cout);

        cout = "  Q:  ";
        for (i = 0; i < m + 1; i++)
        {
            cout += q[i].ToString().PadLeft(12);
        }

        Console.WriteLine(cout);

        Fraction.rfrac_to_cfrac(m, p, q, ref h, ref error);

        typeMethods.r8vec_print(2 * m, h, "  Continued fraction coefficients:");

        for (i = 0; i < 2 * m; i++)
        {
            g[i] = 1.0;
        }

        Fraction.cfrac_to_rfrac(2 * m, g, h, ref p, ref q);

        Console.WriteLine("");
        Console.WriteLine("  Recovered rational polynomial:");
        Console.WriteLine("");

        cout = "  P:  ";
        for (i = 0; i < m; i++)
        {
            cout += p[i].ToString().PadLeft(12);
        }

        Console.WriteLine(cout);

        cout = "  Q:  ";
        for (i = 0; i < m + 1; i++)
        {
            cout += q[i].ToString().PadLeft(12);
        }

        Console.WriteLine(cout);
    }

    public static void rfrac_to_jfrac_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RFRAC_TO_JFRAC_TEST tests RFRAC_TO_JFRAC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int MAXM = 10;

        int i;
        int m;
        double[] p = new double[MAXM];
        double[] q = new double[MAXM];
        double[] r = new double[MAXM];
        double[] s = new double[MAXM];
        int seed;
        //
        //  Generate the data, but force Q(M+1) to be 1.  
        //  That will make it easier to see that the two operations are inverses
        //  of each other.  JFRAC_TO_RFRAC is free to scale its output, and chooses
        //  a scaling in which Q(M+1) is 1.
        //
        seed = 123456789;
        m = 6;
        UniformRNG.r8vec_uniform_01(m, ref seed, ref p);
        UniformRNG.r8vec_uniform_01(m + 1, ref seed, ref q);

        for (i = 0; i < m; i++)
        {
            q[i] /= q[m];
        }

        q[m] = 1.0;

        Console.WriteLine("");
        Console.WriteLine("RFRAC_TO_JFRAC_TEST");
        Console.WriteLine("  RFRAC_TO_JFRAC converts a rational polynomial");
        Console.WriteLine("  fraction to a J fraction.");
        Console.WriteLine("");
        Console.WriteLine("  The original rational polynomial coefficients:");
        Console.WriteLine("");

        string cout = "";
        for (i = 0; i < m; i++)
        {
            cout += p[i].ToString().PadLeft(14) + "  ";
        }

        Console.WriteLine(cout);

        cout = "";
        for (i = 0; i < m + 1; i++)
        {
            cout += q[i].ToString().PadLeft(14) + "  ";
        }

        Console.WriteLine(cout);

        Fraction.rfrac_to_jfrac(m, p, q, ref r, ref s);

        Console.WriteLine("");
        Console.WriteLine("  The J fraction coefficients:");
        Console.WriteLine("");

        cout = "";
        for (i = 0; i < m; i++)
        {
            cout += r[i].ToString().PadLeft(14) + "  ";
        }

        Console.WriteLine(cout);
        cout = "";
        for (i = 0; i < m; i++)
        {
            cout += s[i].ToString().PadLeft(14) + "  ";
        }

        Console.WriteLine(cout);

        Fraction.jfrac_to_rfrac(m, r, s, ref p, ref q);

        Console.WriteLine("");
        Console.WriteLine("  The recovered rational polynomial:");
        Console.WriteLine("");

        cout = "";
        for (i = 0; i < m; i++)
        {
            cout += p[i].ToString().PadLeft(14) + "  ";
        }

        Console.WriteLine(cout);

        cout = "";
        for (i = 0; i < m + 1; i++)
        {
            cout += q[i].ToString().PadLeft(14) + "  ";
        }

        Console.WriteLine(cout);
    }


}