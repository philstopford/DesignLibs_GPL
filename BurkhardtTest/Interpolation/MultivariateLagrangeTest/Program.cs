using System;
using Burkardt.Interpolation;
using Burkardt.MonomialNS;
using Burkardt.Types;
using Burkardt.PolynomialNS;

namespace MultivariateLagrangeTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for LAGRANGE_ND_TEST.
        //
        //  Discussion:
        //
        //    LAGRANGE_ND_TEST tests the LAGRANGE_ND library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int option;

        Console.WriteLine("");
        Console.WriteLine("LAGRANGE_ND_TEST");
        Console.WriteLine("  Test the LAGRANGE_ND library.");

        test01();
        test02();
        test03();
        test04();
        test05();
        test06();
        test07();
        test08();
        test09();
        test10();

        option = 0;
        test11(option);

        option = 1;
        test11(option);
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("LAGRANGE_ND_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests MONO_BETWEEN_ENUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int d;
        int n1;
        int n2;
        int v;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  MONO_BETWEEN_ENUM can enumerate the number of monomials");
        Console.WriteLine("  in D variables, of total degree between N1 and N2.");

        d = 3;
        Console.WriteLine("");
        Console.WriteLine("  Using spatial dimension D = " + d + "");
        Console.WriteLine("");
        string cout = "   N2:";
        for (n2 = 0; n2 <= 8; n2++)
        {
            cout += "  " + n2.ToString().PadLeft(4);
        }

        Console.WriteLine(cout);
        Console.WriteLine("  N1 +------------------------------------------------------");
        for (n1 = 0; n1 <= 8; n1++)
        {
            cout = "  " + n1.ToString().PadLeft(2) + " |";
            for (n2 = 0; n2 <= 8; n2++)
            {
                v = Monomial.mono_between_enum(d, n1, n2);
                cout += "  " + v.ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
        }
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests MONO_TOTAL_ENUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int d;
        int n;
        int v;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  MONO_TOTAL_ENUM can enumerate the number of monomials");
        Console.WriteLine("  in D variables, of total degree N.");

        Console.WriteLine("");
        string cout = "    N:";
        for (n = 0; n <= 8; n++)
        {
            cout += n.ToString().PadLeft(4);
        }

        Console.WriteLine(cout);
        Console.WriteLine("   D +------------------------------------------------------");
        for (d = 1; d <= 8; d++)
        {
            cout = "  " + d.ToString().PadLeft(2) + " |";
            for (n = 0; n <= 8; n++)
            {
                v = Monomial.mono_total_enum(d, n);
                cout += "  " + v.ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
        }
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests MONO_UPTO_ENUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int d;
        int n;
        int v;

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  MONO_UPTO_ENUM can enumerate the number of monomials");
        Console.WriteLine("  in D variables, of total degree 0 up to N.");

        Console.WriteLine("");
        string cout = "    N:";
        for (n = 0; n <= 8; n++)
        {
            cout += "  " + n.ToString().PadLeft(4);
        }

        Console.WriteLine(cout);
        Console.WriteLine("   D +------------------------------------------------------");
        for (d = 1; d <= 8; d++)
        {
            cout = "  " + d.ToString().PadLeft(2) + " |";
            for (n = 0; n <= 8; n++)
            {
                v = Monomial.mono_upto_enum(d, n);
                cout += " " + v.ToString().PadLeft(5);
            }

            Console.WriteLine(cout);
        }
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests MONO_BETWEEN_NEXT_GRLEX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int d = 3;
        int i;
        int j;
        int n1;
        int n2;
        int[] x = new int[3];

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  MONO_BETWEEN_NEXT_GRLEX can list the monomials");
        Console.WriteLine("  in D variables, of total degree N between N1 and N2,");
        Console.WriteLine("  one at a time.");
        Console.WriteLine("");
        Console.WriteLine("  We start the process with (0,0,...,0,N1).");
        Console.WriteLine("  The process ends with (N2,0,...,0,0)");

        n1 = 2;
        n2 = 3;

        Console.WriteLine("");
        Console.WriteLine("  Let D =  " + d + "");
        Console.WriteLine("      N1 = " + n1 + "");
        Console.WriteLine("      N2 = " + n2 + "");
        Console.WriteLine("");

        x[0] = 0;
        x[1] = 0;
        x[2] = n1;

        j = 1;

        for (;;)
        {
            string cout = "  " + j.ToString().PadLeft(2) + ":";
            for (i = 0; i < d; i++)
            {
                cout += "  " + x[i].ToString().PadLeft(1);
            }

            Console.WriteLine(cout);

            if (x[0] == n2)
            {
                break;
            }

            Monomial.mono_between_next_grlex(d, n1, n2, ref x);
            j += 1;
        }
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests LAGRANGE_ND_COMPLETE in 1D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] c;
        int d;
        int[] e;
        double error;
        int i;
        int j;
        string label;
        int n;
        int nd;
        int o;
        double[] pc;
        int[] pe;
        int[] po;
        int r;
        double[] v;
        double[] value;
        double[] xd =  {
                0.0, 1.0, 2.0, 3.0, 4.0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  LAGRANGE_COMPLETE determines");
        Console.WriteLine("  the Lagrange interpolating polynomials L(x)");
        Console.WriteLine("  for ND points in D dimensions, assuming that");
        Console.WriteLine("  the number of points exactly coincides with");
        Console.WriteLine("  R = Pi(D,N), the number of monomials of degree N or less");
        Console.WriteLine("");
        Console.WriteLine("  As a special demonstration, this code runs in 1D");

        nd = 5;

        d = 1;
        n = 4;
        r = Monomial.mono_upto_enum(d, n);

        pc = new double[nd * r];
        pe = new int[nd * r];
        po = new int[nd];

        c = new double[r];
        e = new int[r];

        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension D = " + d + "");
        Console.WriteLine("  Maximum degree N = " + n + "");
        Console.WriteLine("  Number of monomials R = " + r + "");
        Console.WriteLine("  Number of data points ND = " + nd + "");

        typeMethods.r8mat_transpose_print(d, nd, xd, "  Data points XD:");

        MultivariateLagrange.lagrange_complete(d, n, r, nd, xd, ref po, ref pc, ref pe);
        //
        //  Print the polynomials.
        //
        Console.WriteLine("");
        Console.WriteLine("  Lagrange polynomials for XD data points:");
        Console.WriteLine("");

        for (i = 0; i < nd; i++)
        {
            o = po[i];
            for (j = 0; j < o; j++)
            {
                c[j] = pc[i + j * nd];
                e[j] = pe[i + j * nd];
            }

            label =  "  P(" + i + ")(x) =";
            Polynomial.polynomial_print(d, o, c, e, label);
        }

        //
        //  Evaluate the polynomials at XD.
        //
        value = new double[nd * nd];

        for (i = 0; i < nd; i++)
        {
            o = po[i];
            for (j = 0; j < o; j++)
            {
                c[j] = pc[i + j * nd];
                e[j] = pe[i + j * nd];
            }

            v = Polynomial.polynomial_value(d, o, c, e, nd, xd);

            for (j = 0; j < nd; j++)
            {
                value[i + j * nd] = v[j];
            }
        }

        error = typeMethods.r8mat_is_identity(nd, value);
        Console.WriteLine("");
        Console.WriteLine("  Frobenius norm of Lagrange matrix error = " + error + "");
    }

    private static void test06()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 tests LAGRANGE_ND_COMPLETE in 2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] c;
        int d;
        int[] e;
        double error;
        int i;
        int j;
        string label;
        int n;
        int nd;
        int o;
        double[] pc;
        int[] pe;
        int[] po;
        int r;
        double[] v;
        double[] value;
        double[] xd =  {
                0.0, 0.0,
                1.0, 0.0,
                2.0, 0.0,
                0.0, 1.0,
                1.0, 1.0,
                0.0, 2.0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("TEST06");
        Console.WriteLine("  LAGRANGE_COMPLETE determines");
        Console.WriteLine("  the Lagrange interpolating polynomials L(x)");
        Console.WriteLine("  for ND points in D dimensions, assuming that");
        Console.WriteLine("  the number of points exactly coincides with");
        Console.WriteLine("  R = Pi(D,N), the number of monomials of degree N or less");
        Console.WriteLine("");
        Console.WriteLine("  The data points are the grid nodes of a triangle.");

        nd = 6;

        d = 2;
        n = 2;
        r = Monomial.mono_upto_enum(d, n);

        pc = new double[nd * r];
        pe = new int[nd * r];
        po = new int[nd];

        c = new double[r];
        e = new int[r];

        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension D = " + d + "");
        Console.WriteLine("  Maximum degree N = " + n + "");
        Console.WriteLine("  Number of monomials R = " + r + "");
        Console.WriteLine("  Number of data points ND = " + nd + "");

        typeMethods.r8mat_transpose_print(d, nd, xd, "  Data points XD:");

        MultivariateLagrange.lagrange_complete(d, n, r, nd, xd, ref po, ref pc, ref pe);
        //
        //  Print the polynomials.
        //
        Console.WriteLine("");
        Console.WriteLine("  Lagrange polynomials for XD data points:");
        Console.WriteLine("");

        for (i = 0; i < nd; i++)
        {
            o = po[i];
            for (j = 0; j < o; j++)
            {
                c[j] = pc[i + j * nd];
                e[j] = pe[i + j * nd];
            }

            label = "  P(" + i + ")(x) =";
            Polynomial.polynomial_print(d, o, c, e, label);
        }

        //
        //  Evaluate the polynomials at XD.
        //
        value = new double[nd * nd];

        for (i = 0; i < nd; i++)
        {
            o = po[i];
            for (j = 0; j < o; j++)
            {
                c[j] = pc[i + j * nd];
                e[j] = pe[i + j * nd];
            }

            v = Polynomial.polynomial_value(d, o, c, e, nd, xd);

            for (j = 0; j < nd; j++)
            {
                value[i + j * nd] = v[j];
            }
        }

        error = typeMethods.r8mat_is_identity(nd, value);
        Console.WriteLine("");
        Console.WriteLine("  Frobenius norm of Lagrange matrix error = " + error + "");
    }

    private static void test07()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST07 tests LAGRANGE_ND_COMPLETE in 3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] c;
        int d;
        int[] e;
        double error;
        int i;
        int j;
        string label;
        int n;
        int nd;
        int o;
        double[] pc;
        int[] pe;
        int[] po;
        int r;
        double[] v;
        double[] value;
        double[] xd =  {
                0.0, 0.0, 0.0,
                1.0, 0.0, 0.0,
                2.0, 0.0, 0.0,
                0.0, 1.0, 0.0,
                1.0, 1.0, 0.0,
                0.0, 2.0, 0.0,
                0.0, 0.0, 1.0,
                1.0, 0.0, 1.0,
                0.0, 1.0, 1.0,
                0.0, 0.0, 2.0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("TEST07");
        Console.WriteLine("  LAGRANGE_COMPLETE determines");
        Console.WriteLine("  the Lagrange interpolating polynomials L(x)");
        Console.WriteLine("  for ND points in D dimensions, assuming that");
        Console.WriteLine("  the number of points exactly coincides with");
        Console.WriteLine("  R = Pi(D,N), the number of monomials of degree N or less");
        Console.WriteLine("");
        Console.WriteLine("  The data points are the grid nodes of a tetrahedron.");

        nd = 10;

        d = 3;
        n = 2;
        r = Monomial.mono_upto_enum(d, n);

        pc = new double[nd * r];
        pe = new int[nd * r];
        po = new int[nd];

        c = new double[r];
        e = new int[r];

        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension D = " + d + "");
        Console.WriteLine("  Maximum degree N = " + n + "");
        Console.WriteLine("  Number of monomials R = " + r + "");
        Console.WriteLine("  Number of data points ND = " + nd + "");

        typeMethods.r8mat_transpose_print(d, nd, xd, "  Data points XD:");

        MultivariateLagrange.lagrange_complete(d, n, r, nd, xd, ref po, ref pc, ref pe);
        //
        //  Print the polynomials.
        //
        Console.WriteLine("");
        Console.WriteLine("  Lagrange polynomials for XD data points:");
        Console.WriteLine("");

        for (i = 0; i < nd; i++)
        {
            o = po[i];
            for (j = 0; j < o; j++)
            {
                c[j] = pc[i + j * nd];
                e[j] = pe[i + j * nd];
            }

            label = "  P(" + i + ")(x) =";
            Polynomial.polynomial_print(d, o, c, e, label);
        }

        //
        //  Evaluate the polynomials at XD.
        //
        value = new double[nd * nd];

        for (i = 0; i < nd; i++)
        {
            o = po[i];
            for (j = 0; j < o; j++)
            {
                c[j] = pc[i + j * nd];
                e[j] = pe[i + j * nd];
            }

            v = Polynomial.polynomial_value(d, o, c, e, nd, xd);

            for (j = 0; j < nd; j++)
            {
                value[i + j * nd] = v[j];
            }
        }

        error = typeMethods.r8mat_is_identity(nd, value);
        Console.WriteLine("");
        Console.WriteLine("  Frobenius norm of Lagrange matrix error = " + error + "");
    }

    private static void test08()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST08 tests LAGRANGE_PARTIAL in 2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] c;
        int d;
        int[] e;
        double error;
        int i;
        int j;
        string label;
        int n;
        int nd;
        int o;
        double[] pc;
        int[] pe;
        int[] po;
        int r;
        double[] v;
        double[] value;
        double[] xd =  {
                0.0, 0.0,
                -1.0, 0.0,
                1.0, 0.0,
                0.0, -1.0,
                0.0, 1.0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("TEST08");
        Console.WriteLine("  LAGRANGE_PARTIAL determines");
        Console.WriteLine("  the Lagrange interpolating polynomials L(x)");
        Console.WriteLine("  for ND points in D dimensions, assuming that");
        Console.WriteLine("  the number of points is less than or equal to");
        Console.WriteLine("  R = Pi(D,N), the number of monomials of degree N or less");
        Console.WriteLine("");
        Console.WriteLine("  For this example, the data points are the same as those");
        Console.WriteLine("  used by the level 1 Clenshaw Curtis sparse grid in 2D.");

        nd = 5;

        d = 2;
        n = 2;
        r = Monomial.mono_upto_enum(d, n);

        pc = new double[nd * r];
        pe = new int[nd * r];
        po = new int[nd];

        c = new double[r];
        e = new int[r];

        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension D = " + d + "");
        Console.WriteLine("  Maximum degree N = " + n + "");
        Console.WriteLine("  Number of monomials R = " + r + "");
        Console.WriteLine("  Number of data points ND = " + nd + "");

        typeMethods.r8mat_transpose_print(d, nd, xd, "  Data points XD:");

        MultivariateLagrange.lagrange_partial(d, n, r, nd, xd, ref po, ref pc, ref pe);
        //
        //
        //  Print the polynomials.
        //
        Console.WriteLine("");
        Console.WriteLine("  Lagrange polynomials for XD data points:");
        Console.WriteLine("");

        for (i = 0; i < nd; i++)
        {
            o = po[i];
            for (j = 0; j < o; j++)
            {
                c[j] = pc[i + j * nd];
                e[j] = pe[i + j * nd];
            }

            label = "  P(" + i + ")(x) =";
            Polynomial.polynomial_print(d, o, c, e, label);
        }

        //
        //  Evaluate the polynomials at XD.
        //
        value = new double[nd * nd];

        for (i = 0; i < nd; i++)
        {
            o = po[i];
            for (j = 0; j < o; j++)
            {
                c[j] = pc[i + j * nd];
                e[j] = pe[i + j * nd];
            }

            v = Polynomial.polynomial_value(d, o, c, e, nd, xd);

            for (j = 0; j < nd; j++)
            {
                value[i + j * nd] = v[j];
            }

                
        }

        error = typeMethods.r8mat_is_identity(nd, value);
        Console.WriteLine("");
        Console.WriteLine("  Frobenius norm of Lagrange matrix error = " + error + "");
    }
    //****************************************************************************80

    private static void test09()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST09 tests LAGRANGE_PARTIAL in 3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] c;
        int d;
        int[] e;
        double error;
        int i;
        int j;
        string label;
        int n;
        int nd;
        int o;
        double[] pc;
        int[] pe;
        int[] po;
        int r;
        double[] v;
        double[] value;
        double[] xd =  {
                0.0, 0.0, 0.0,
                -1.0, 0.0, 0.0,
                1.0, 0.0, 0.0,
                0.0, -1.0, 0.0,
                0.0, 1.0, 0.0,
                0.0, 0.0, -1.0,
                0.0, 0.0, 1.0,
                -0.707106781187, 0.0, 0.0,
                0.707106781187, 0.0, 0.0,
                -1.0, -1.0, 0.0,
                1.0, -1.0, 0.0,
                -1.0, 1.0, 0.0,
                1.0, 1.0, 0.0,
                0.0, -0.707106781187, 0.0,
                0.0, 0.707106781187, 0.0,
                -1.0, 0.0, -1.0,
                1.0, 0.0, -1.0,
                -1.0, 0.0, 1.0,
                1.0, 0.0, 1.0,
                0.0, -1.0, -1.0,
                0.0, 1.0, -1.0,
                0.0, -1.0, 1.0,
                0.0, 1.0, 1.0,
                0.0, 0.0, -0.707106781187,
                0.0, 0.0, 0.707106781187
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("TEST09");
        Console.WriteLine("  LAGRANGE_PARTIAL determines");
        Console.WriteLine("  the Lagrange interpolating polynomials L(x)");
        Console.WriteLine("  for ND points in D dimensions, assuming that");
        Console.WriteLine("  the number of points is less than or equal to");
        Console.WriteLine("  R = Pi(D,N), the number of monomials of degree N or less");
        Console.WriteLine("");
        Console.WriteLine("  For this example, the data points are the same as those");
        Console.WriteLine("  used by the level 2 Clenshaw Curtis sparse grid in 3D.");

        nd = 25;

        d = 3;
        n = 4;
        r = Monomial.mono_upto_enum(d, n);

        pc = new double[nd * r];
        pe = new int[nd * r];
        po = new int[nd];

        c = new double[r];
        e = new int[r];

        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension D = " + d + "");
        Console.WriteLine("  Maximum degree N = " + n + "");
        Console.WriteLine("  Number of monomials R = " + r + "");
        Console.WriteLine("  Number of data points ND = " + nd + "");

        typeMethods.r8mat_transpose_print(d, nd, xd, "  Data points XD:");

        MultivariateLagrange.lagrange_partial(d, n, r, nd, xd, ref po, ref pc, ref pe);
        //
        //
        //  Print the polynomials.
        //
        Console.WriteLine("");
        Console.WriteLine("  Lagrange polynomials for XD data points:");
        Console.WriteLine("");

        for (i = 0; i < nd; i++)
        {
            o = po[i];
            for (j = 0; j < o; j++)
            {
                c[j] = pc[i + j * nd];
                e[j] = pe[i + j * nd];
            }

            label = "  P(" + i + ")(x) =";
            Polynomial.polynomial_print(d, o, c, e, label);
        }

        //
        //  Evaluate the polynomials at XD.
        //
        value = new double[nd * nd];

        for (i = 0; i < nd; i++)
        {
            o = po[i];
            for (j = 0; j < o; j++)
            {
                c[j] = pc[i + j * nd];
                e[j] = pe[i + j * nd];
            }

            v = Polynomial.polynomial_value(d, o, c, e, nd, xd);

            for (j = 0; j < nd; j++)
            {
                value[i + j * nd] = v[j];
            }

                
        }

        error = typeMethods.r8mat_is_identity(nd, value);
        Console.WriteLine("");
        Console.WriteLine("  Frobenius norm of Lagrange matrix error = " + error + "");
    }

    private static void test10()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST10 tests LAGRANGE_PARTIAL2 in 2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] c;
        int d;
        int[] e;
        double error;
        double f;
        int i;
        int j;
        int k;
        string label;
        int n;
        int nd;
        int ni;
        int o;
        double[] pc;
        double[] pd;
        int[] pe;
        int pn;
        int[] po;
        int r;
        double[] v;
        double[] value;
        double[] xd =  {
                0.0, 0.0,
                -1.0, 0.0,
                1.0, 0.0,
                0.0, -1.0,
                0.0, 1.0,
                -1.0, 1.0,
                1.0, 1.0,
                -1.0, -1.0,
                1.0, -1.0,
                -0.5, 0.0,
                0.0, -0.5,
                0.0, +0.5,
                0.5, 0.0
            }
            ;
        double[] xyi;
        double[] zi;

        Console.WriteLine("");
        Console.WriteLine("TEST10");
        Console.WriteLine("  LAGRANGE_PARTIAL2 determines");
        Console.WriteLine("  the Lagrange interpolating polynomials L(x)");
        Console.WriteLine("  for ND points in D dimensions, assuming that");
        Console.WriteLine("  the number of points is less than or equal to");
        Console.WriteLine("  R = Pi(D,N), the number of monomials of degree N or less");
        Console.WriteLine("");
        Console.WriteLine("  For this example, the data points are the same as those");
        Console.WriteLine("  used by the level 2 Clenshaw Curtis sparse grid in 2D.");

        nd = 13;
        ni = 11 * 11;

        d = 2;
        n = 4;
        r = Monomial.mono_upto_enum(d, n);

        pc = new double[nd * r];
        pe = new int[nd * r];
        po = new int[nd];

        c = new double[r];
        e = new int[r];

        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension D = " + d + "");
        Console.WriteLine("  Maximum degree N = " + n + "");
        Console.WriteLine("  Number of monomials R = " + r + "");
        Console.WriteLine("  Number of data points ND = " + nd + "");

        typeMethods.r8mat_transpose_print(d, nd, xd, "  Data points XD:");

        MultivariateLagrange.lagrange_partial2(d, n, r, nd, xd, ref po, ref pc, ref pe);
        //
        //  Print the polynomials.
        //
        Console.WriteLine("");
        Console.WriteLine("  Lagrange polynomials for XD data points:");
        Console.WriteLine("");

        for (i = 0; i < nd; i++)
        {
            o = po[i];
            for (j = 0; j < o; j++)
            {
                c[j] = pc[i + j * nd];
                e[j] = pe[i + j * nd];
            }

            label = "  P(" + i + ")(x) =";
            Polynomial.polynomial_print(d, o, c, e, label);
        }

        //
        //  Evaluate the polynomials at XD.
        //
        value = new double[nd * nd];

        for (i = 0; i < nd; i++)
        {
            o = po[i];
            for (j = 0; j < o; j++)
            {
                c[j] = pc[i + j * nd];
                e[j] = pe[i + j * nd];
            }

            v = Polynomial.polynomial_value(d, o, c, e, nd, xd);

            for (j = 0; j < nd; j++)
            {
                value[i + j * nd] = v[j];
            }

                
        }

        error = typeMethods.r8mat_is_identity(nd, value);
        Console.WriteLine("");
        Console.WriteLine("  Frobenius norm of Lagrange matrix error = " + error + "");
        //
        //  Evaluate a function at the data points.
        //
        pd = new double[nd];
        for (i = 0; i < nd; i++)
        {
            pd[i] = Math.Sin(xd[0 + i * 2]) * Math.Cos(xd[1 + i * 2]);
        }

        //
        //  Compare exact function and interpolant at a grid of points.
        //
        xyi = new double[2 * ni];

        k = 0;
        for (j = 1; j <= 11; j++)
        {
            for (i = 1; i <= 11; i++)
            {
                xyi[0 + k * 2] = ((11 - i) * -1.0
                                  + (i - 1) * +1.0)
                                 / (11 - 1);
                xyi[1 + k * 2] = ((11 - j) * -1.0
                                  + (j - 1) * +1.0)
                                 / (11 - 1);
                k += 1;
            }
        }

        pn = nd;
        zi = MultivariateLagrange.interpolant_value(d, r, pn, po, pc, pe, pd, ni, xyi);

        error = 0.0;
        for (k = 0; k < ni; k++)
        {
            f = Math.Sin(xyi[0 + k * 2]) * Math.Cos(xyi[1 + k * 2]);
            if (error < Math.Abs(zi[k] - f))
            {
                error = Math.Abs(zi[k] - f);
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  Maximum absolute interpolant error on 11x11 grid = " + error + "");
    }

    private static void test11(int option)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGRANGE_ND_TEST11 tests LAGRANGE_PARTIAL3 in 2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 February 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int OPTION, determines the initial basis:
        //    0, use monomials, 1, x, y, x^2, xy, y^2, x^3, ...
        //    1, use Legendre products, 1, y, x, (3y^2-1)/2, xy, (3^x^2-1), (5y^3-3)/2,...
        //
    {
        double[] c;
        int d;
        int[] e;
        double error;
        double f;
        int i;
        int j;
        int k;
        string label;
        int n;
        int n2 = 0;
        int nd;
        int ni;
        int o;
        double[] pc = new double[1];
        double[] pd;
        int[] pe = new int[1];
        int pn;
        int[] po;
        int r;
        double[] v;
        double[] value;
        double[] xd =  {
                0.0000000000000000, 0.0000000000000000,
                -1.0000000000000000, 0.0000000000000000,
                1.0000000000000000, 0.0000000000000000,
                0.0000000000000000, -1.0000000000000000,
                0.0000000000000000, 1.0000000000000000,
                -0.7071067811865475, 0.0000000000000000,
                0.7071067811865476, 0.0000000000000000,
                -1.0000000000000000, -1.0000000000000000,
                1.0000000000000000, -1.0000000000000000,
                -1.0000000000000000, 1.0000000000000000,
                1.0000000000000000, 1.0000000000000000,
                0.0000000000000000, -0.7071067811865475,
                0.0000000000000000, 0.7071067811865476,
                -0.9238795325112867, 0.0000000000000000,
                -0.3826834323650897, 0.0000000000000000,
                0.3826834323650898, 0.0000000000000000,
                0.9238795325112867, 0.0000000000000000,
                -0.7071067811865475, -1.0000000000000000,
                0.7071067811865476, -1.0000000000000000,
                -0.7071067811865475, 1.0000000000000000,
                0.7071067811865476, 1.0000000000000000,
                -1.0000000000000000, -0.7071067811865475,
                1.0000000000000000, -0.7071067811865475,
                -1.0000000000000000, 0.7071067811865476,
                1.0000000000000000, 0.7071067811865476,
                0.0000000000000000, -0.9238795325112867,
                0.0000000000000000, -0.3826834323650897,
                0.0000000000000000, 0.3826834323650898,
                0.0000000000000000, 0.9238795325112867,
                -0.9807852804032304, 0.0000000000000000,
                -0.8314696123025453, 0.0000000000000000,
                -0.5555702330196020, 0.0000000000000000,
                -0.1950903220161282, 0.0000000000000000,
                0.1950903220161283, 0.0000000000000000,
                0.5555702330196023, 0.0000000000000000,
                0.8314696123025452, 0.0000000000000000,
                0.9807852804032304, 0.0000000000000000,
                -0.9238795325112867, -1.0000000000000000,
                -0.3826834323650897, -1.0000000000000000,
                0.3826834323650898, -1.0000000000000000,
                0.9238795325112867, -1.0000000000000000,
                -0.9238795325112867, 1.0000000000000000,
                -0.3826834323650897, 1.0000000000000000,
                0.3826834323650898, 1.0000000000000000,
                0.9238795325112867, 1.0000000000000000,
                -0.7071067811865475, -0.7071067811865475,
                0.7071067811865476, -0.7071067811865475,
                -0.7071067811865475, 0.7071067811865476,
                0.7071067811865476, 0.7071067811865476,
                -1.0000000000000000, -0.9238795325112867,
                1.0000000000000000, -0.9238795325112867,
                -1.0000000000000000, -0.3826834323650897,
                1.0000000000000000, -0.3826834323650897,
                -1.0000000000000000, 0.3826834323650898,
                1.0000000000000000, 0.3826834323650898,
                -1.0000000000000000, 0.9238795325112867,
                1.0000000000000000, 0.9238795325112867,
                0.0000000000000000, -0.9807852804032304,
                0.0000000000000000, -0.8314696123025453,
                0.0000000000000000, -0.5555702330196020,
                0.0000000000000000, -0.1950903220161282,
                0.0000000000000000, 0.1950903220161283,
                0.0000000000000000, 0.5555702330196023,
                0.0000000000000000, 0.8314696123025452,
                0.0000000000000000, 0.9807852804032304
            }
            ;
        double[] xyi;
        double[] zi;

        Console.WriteLine("");
        Console.WriteLine("LAGRANGE_ND_TEST11");
        Console.WriteLine("  LAGRANGE_PARTIAL3 determines");
        Console.WriteLine("  the Lagrange interpolating polynomials L(x)");
        Console.WriteLine("  for ND points in D dimensions, assuming that");
        Console.WriteLine("  the number of points is less than or equal to");
        Console.WriteLine("  R = Pi(D,N), the number of monomials of degree N or less");
        Console.WriteLine("");
        Console.WriteLine("  If LAGRANGE_PARTIAL3 determines that the problem is not");
        Console.WriteLine("  well-posed for the given value of N, it increases N");
        Console.WriteLine("  until a suitable value is found.");
        Console.WriteLine("");
        Console.WriteLine("  For this example, the data points are the same as those");
        Console.WriteLine("  used by the level 2 Clenshaw Curtis sparse grid in 2D.");

        nd = 65;
        ni = 11 * 11;

        d = 2;
        n = 10;
        po = new int[nd];

        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension D = " + d + "");
        Console.WriteLine("  Maximum degree N = " + n + "");
        Console.WriteLine("  Number of data points ND = " + nd + "");
        Console.WriteLine("  Monomial/Legendre option OPTION = " + option + "");

        typeMethods.r8mat_transpose_print(d, nd, xd, "  Data points XD:");

        MultivariateLagrange.lagrange_partial3(d, n, nd, xd, option, ref po, ref pc, ref pe, ref n2);

        if (n < n2)
        {
            Console.WriteLine("");
            Console.WriteLine("  LAGRANGE_PARTIAL3 increased N to " + n2 + "");
        }

        r = Monomial.mono_upto_enum(d, n2);
        Console.WriteLine("  Number of monomials R = " + r + "");
        c = new double[r];
        e = new int[r];
        //
        //  Print the polynomials.
        //
        Console.WriteLine("");
        Console.WriteLine("  (First 2) Lagrange polynomials for XD data points:");
        Console.WriteLine("");

        //for ( i = 0; i < nd; i++ )
        for (i = 0; i < 2; i++)
        {
            o = po[i];
            for (j = 0; j < o; j++)
            {
                c[j] = pc[i + j * nd];
                e[j] = pe[i + j * nd];
            }

            label = "  P(" + i + ")(x) =";
            Polynomial.polynomial_print(d, o, c, e, label);
        }

        //
        //  Evaluate the polynomials at XD.
        //
        value = new double[nd * nd];

        for (i = 0; i < nd; i++)
        {
            o = po[i];
            for (j = 0; j < o; j++)
            {
                c[j] = pc[i + j * nd];
                e[j] = pe[i + j * nd];
            }

            v = Polynomial.polynomial_value(d, o, c, e, nd, xd);

            for (j = 0; j < nd; j++)
            {
                value[i + j * nd] = v[j];
            }

                
        }

        error = typeMethods.r8mat_is_identity(nd, value);
        Console.WriteLine("");
        Console.WriteLine("  Frobenius norm of Lagrange matrix error = " + error + "");
        //
        //  Evaluate a function at the data points.
        //
        pd = new double[nd];
        for (i = 0; i < nd; i++)
        {
            pd[i] = Math.Exp(xd[0 + i * 2] * xd[1 + i * 2]);
        }

        //
        //  Compare exact function and interpolant at a grid of points.
        //
        xyi = new double[2 * ni];

        k = 0;
        for (j = 1; j <= 11; j++)
        {
            for (i = 1; i <= 11; i++)
            {
                xyi[0 + k * 2] = ((11 - i) * -1.0
                                  + (i - 1) * +1.0)
                                 / (11 - 1);
                xyi[1 + k * 2] = ((11 - j) * -1.0
                                  + (j - 1) * +1.0)
                                 / (11 - 1);
                k += 1;
            }
        }

        pn = nd;
        zi = MultivariateLagrange.interpolant_value(d, r, pn, po, pc, pe, pd, ni, xyi);

        error = 0.0;
        for (k = 0; k < ni; k++)
        {
            f = Math.Exp(xyi[0 + k * 2] * xyi[1 + k * 2]);
            if (error < Math.Abs(zi[k] - f))
            {
                error = Math.Abs(zi[k] - f);
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  Maximum absolute interpolant error on 11x11 grid = " + error + "");
    }
}