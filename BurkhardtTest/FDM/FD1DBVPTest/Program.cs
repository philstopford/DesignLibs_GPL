using System;
using Burkardt.FDM;
using Burkardt.Types;

namespace FDM1DBVPTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for FD1D_BVP_TEST.
        //
        //  Discussion:
        //
        //    FD1D_BVP_TEST tests the FD1D_BVP library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("FD1D_BVP_TEST");
        Console.WriteLine("  Test the FD1D_BVP library.");

        fd1d_bvp_test01();
        fd1d_bvp_test02();
        fd1d_bvp_test03();
        fd1d_bvp_test04();
        fd1d_bvp_test05();

        Console.WriteLine("");
        Console.WriteLine("FD1D_BVP_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void fd1d_bvp_test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FD1D_BVP_TEST01 carries out test case #1.
        //
        //  Discussion:
        //
        //    Use A1, C1, F1, EXACT1.
        //
        //    Repeat using a nonuniform mesh.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        string filename;
        int i;
        int n = 21;
        double[] u;
        double[] u2;
        double[] uexact;
        double[] x;
        double x1 = 0.0;
        double x2 = 1.0;

        Console.WriteLine("");
        Console.WriteLine("FD1D_BVP_TEST01");
        Console.WriteLine("  A1(X)  = 1.0");
        Console.WriteLine("  A1'(X) = 0.0");
        Console.WriteLine("  C1(X)  = 0.0");
        Console.WriteLine("  F1(X)  = X * ( X + 3 ) * exp ( X )");
        Console.WriteLine("  U1(X)  = X * ( 1 - X ) * exp ( X )");
        Console.WriteLine("");
        Console.WriteLine("  Number of nodes = " + n + "");
        Console.WriteLine("  X1 = " + x1 + "");
        Console.WriteLine("  X2 = " + x2 + "");

        x = typeMethods.r8vec_even(n, x1, x2);

        u = FDM_1D_BVP.fd1d_bvp(n, a1, a1prime, c1, f1, x);

        uexact = exact1(n, x);

        Console.WriteLine("");
        Console.WriteLine("     I         X        U             Uexact         Error");
        Console.WriteLine("");

        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + uexact[i].ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + Math.Abs(u[i] - uexact[i]).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Repeat, using a nonuniform mesh.");

        x = typeMethods.r8vec_even(n, x1, x2);

        for (i = 0; i < n; i++)
        {
            x[i] = Math.Sqrt(x[i]);
        }

        u = FDM_1D_BVP.fd1d_bvp(n, a1, a1prime, c1, f1, x);

        uexact = exact1(n, x);

        Console.WriteLine("");
        Console.WriteLine("     I         X        U             Uexact         Error");
        Console.WriteLine("");

        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + uexact[i].ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + Math.Abs(u[i] - uexact[i]).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        //
        //  Write the data to files.
        //
        filename = "fd1d_bvp_test01_nodes.txt";
        typeMethods.r8mat_write(filename, n, 1, x);

        u2 = new double[n * 2];
        for (i = 0; i < n; i++)
        {
            u2[i + 0 * n] = u[i];
        }

        for (i = 0; i < n; i++)
        {
            u2[i + 1 * n] = uexact[i];
        }

        filename = "fd1d_bvp_test01_values.txt";
        typeMethods.r8mat_write(filename, n, 2, u2);
    }

    private static void fd1d_bvp_test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FD1D_BVP_TEST02 carries out test case #2.
        //
        //  Discussion:
        //
        //    Use A1, C2, F2, EXACT1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        string filename;
        int i;
        int n = 11;
        double[] u;
        double[] u2;
        double[] uexact;
        double[] x;
        double x1 = 0.0;
        double x2 = 1.0;

        Console.WriteLine("");
        Console.WriteLine("FD1D_BVP_TEST02");
        Console.WriteLine("  A1(X)  = 1.0");
        Console.WriteLine("  A1''(X) = 0.0");
        Console.WriteLine("  C2(X)  = 2.0");
        Console.WriteLine("  F2(X)  = X * ( 5 - X ) * exp ( X )");
        Console.WriteLine("  U1(X)  = X * ( 1 - X ) * exp ( X )");
        Console.WriteLine("");
        Console.WriteLine("  Number of nodes = " + n + "");
        Console.WriteLine("  X1 = " + x1 + "");
        Console.WriteLine("  X2 = " + x2 + "");

        x = typeMethods.r8vec_even(n, x1, x2);

        u = FDM_1D_BVP.fd1d_bvp(n, a1, a1prime, c2, f2, x);

        uexact = exact1(n, x);

        Console.WriteLine("");
        Console.WriteLine("     I         X        U             Uexact         Error");
        Console.WriteLine("");

        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + uexact[i].ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + Math.Abs(u[i] - uexact[i]).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        //
        //  Write the data to files.
        //
        filename = "fd1d_bvp_test02_nodes.txt";
        typeMethods.r8mat_write(filename, n, 1, x);

        u2 = new double[n * 2];
        for (i = 0; i < n; i++)
        {
            u2[i + 0 * n] = u[i];
        }

        for (i = 0; i < n; i++)
        {
            u2[i + 1 * n] = uexact[i];
        }

        filename = "fd1d_bvp_test02_values.txt";
        typeMethods.r8mat_write(filename, n, 2, u2);
    }

    private static void fd1d_bvp_test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FD1D_BVP_TEST03 carries out test case #3.
        //
        //  Discussion:
        //
        //    Use A1, C3, F3, EXACT1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        string filename;
        int i;
        int n = 11;
        double[] u;
        double[] u2;
        double[] uexact;
        double[] x;
        double x1 = 0.0;
        double x2 = 1.0;

        Console.WriteLine("");
        Console.WriteLine("FD1D_BVP_TEST03");
        Console.WriteLine("  A1(X)  = 1.0");
        Console.WriteLine("  A1''(X) = 0.0");
        Console.WriteLine("  C3(X)  = 2.0 * X");
        Console.WriteLine("  F3(X)  = - X * ( 2 * X * X - 3 * X - 3 ) * exp ( X )");
        Console.WriteLine("  U1(X)  = X * ( 1 - X ) * exp ( X )");
        Console.WriteLine("");
        Console.WriteLine("  Number of nodes = " + n + "");
        Console.WriteLine("  X1 = " + x1 + "");
        Console.WriteLine("  X2 = " + x2 + "");

        x = typeMethods.r8vec_even(n, x1, x2);

        u = FDM_1D_BVP.fd1d_bvp(n, a1, a1prime, c3, f3, x);

        uexact = exact1(n, x);

        Console.WriteLine("");
        Console.WriteLine("     I         X        U             Uexact         Error");
        Console.WriteLine("");

        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + uexact[i].ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + Math.Abs(u[i] - uexact[i]).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        //
        //  Write the data to files.
        //
        filename = "fd1d_bvp_test03_nodes.txt";
        typeMethods.r8mat_write(filename, n, 1, x);

        u2 = new double[n * 2];
        for (i = 0; i < n; i++)
        {
            u2[i + 0 * n] = u[i];
        }

        for (i = 0; i < n; i++)
        {
            u2[i + 1 * n] = uexact[i];
        }

        filename = "fd1d_bvp_test03_values.txt";
        typeMethods.r8mat_write(filename, n, 2, u2);

    }

    private static void fd1d_bvp_test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FD1D_BVP_TEST04 carries out test case #4.
        //
        //  Discussion:
        //
        //    Use A2, C1, F4, EXACT1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        string filename;
        int i;
        int n = 11;
        double[] u;
        double[] u2;
        double[] uexact;
        double[] x;
        double x1 = 0.0;
        double x2 = 1.0;

        Console.WriteLine("");
        Console.WriteLine("FD1D_BVP_TEST04");
        Console.WriteLine("  A2(X)  = 1.0 + X * X");
        Console.WriteLine("  A2''(X) = 2.0 * X");
        Console.WriteLine("  C1(X)  = 0.0");
        Console.WriteLine("  F4(X)  = ( X + 3 X^2 + 5 X^3 + X^4 ) * exp ( X )");
        Console.WriteLine("  U1(X)  = X * ( 1 - X ) * exp ( X )");
        Console.WriteLine("");
        Console.WriteLine("  Number of nodes = " + n + "");
        Console.WriteLine("  X1 = " + x1 + "");
        Console.WriteLine("  X2 = " + x2 + "");

        x = typeMethods.r8vec_even(n, x1, x2);

        u = FDM_1D_BVP.fd1d_bvp(n, a2, a2prime, c1, f4, x);

        uexact = exact1(n, x);

        Console.WriteLine("");
        Console.WriteLine("     I         X        U             Uexact         Error");
        Console.WriteLine("");

        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + uexact[i].ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + Math.Abs(u[i] - uexact[i]).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        //
        //  Write the data to files.
        //
        filename = "fd1d_bvp_test04_nodes.txt";
        typeMethods.r8mat_write(filename, n, 1, x);

        u2 = new double[n * 2];
        for (i = 0; i < n; i++)
        {
            u2[i + 0 * n] = u[i];
        }

        for (i = 0; i < n; i++)
        {
            u2[i + 1 * n] = uexact[i];
        }

        filename = "fd1d_bvp_test04_values.txt";
        typeMethods.r8mat_write(filename, n, 2, u2);

    }

    private static void fd1d_bvp_test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FD1D_BVP_TEST05 carries out test case #5.
        //
        //  Discussion:
        //
        //    Use A3, C1, F5, EXACT1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        string filename;
        int i;
        int n = 11;
        double[] u;
        double[] u2;
        double[] uexact;
        double[] x;
        double x1 = 0.0;
        double x2 = 1.0;

        Console.WriteLine("");
        Console.WriteLine("FD1D_BVP_TEST05");
        Console.WriteLine("  A3(X)  = 1.0 + X * X for X <= 1/3");
        Console.WriteLine("         = 7/9 + X     for      1/3 < X");
        Console.WriteLine("  A3''(X) = 2.0 * X     for X <= 1/3");
        Console.WriteLine("           1           for      1/3 < X");
        Console.WriteLine("  C1(X)  = 0.0");
        Console.WriteLine("  F5(X)  = ( X + 3 X^2 + 5 X^3 + X^4 ) * exp ( X )");
        Console.WriteLine("                       for X <= 1/3");
        Console.WriteLine("         = ( - 1 + 10/3 X + 43/9 X^2 + X^3 ) * exp ( X )");
        Console.WriteLine("                       for      1/3 <= X");
        Console.WriteLine("  U1(X)  = X * ( 1 - X ) * exp ( X )");
        Console.WriteLine("");
        Console.WriteLine("  Number of nodes = " + n + "");
        Console.WriteLine("  X1 = " + x1 + "");
        Console.WriteLine("  X2 = " + x2 + "");

        x = typeMethods.r8vec_even(n, x1, x2);

        u = FDM_1D_BVP.fd1d_bvp(n, a3, a3prime, c1, f5, x);

        uexact = exact1(n, x);

        Console.WriteLine("");
        Console.WriteLine("     I         X        U             Uexact         Error");
        Console.WriteLine("");

        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + uexact[i].ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + Math.Abs(u[i] - uexact[i]).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        //
        //  Write the data to files.
        //
        filename = "fd1d_bvp_test05_nodes.txt";
        typeMethods.r8mat_write(filename, n, 1, x);

        u2 = new double[n * 2];
        for (i = 0; i < n; i++)
        {
            u2[i + 0 * n] = u[i];
        }

        for (i = 0; i < n; i++)
        {
            u2[i + 1 * n] = uexact[i];
        }

        filename = "fd1d_bvp_test05_values.txt";
        typeMethods.r8mat_write(filename, n, 2, u2);

    }

    private static double a1(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    A1 evaluates A function #1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double A1, the value of A(X).
        //
    {
        double value = 1.0;

        return value;
    }

    private static double a1prime(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    A1PRIME evaluates A' function #1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double A1PRIME, the value of A'(X).
        //
    {
        double value = 0.0;

        return value;
    }

    private static double a2(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    A2 evaluates A function #2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double A2, the value of A(X).
        //
    {
        double value = 1.0 + x * x;

        return value;
    }

    private static double a2prime(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    A2PRIME evaluates A' function #2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double A2PRIME, the value of A'(X).
        //
    {
        double value = 2.0 * x;

        return value;
    }

    private static double a3(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    A3 evaluates A function #3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double A3, the value of A(X).
        //
    {
        double value = x switch
        {
            <= 1.0 / 3.0 => 1.0 + x * x,
            _ => x + 7.0 / 9.0
        };

        return value;
    }

    private static double a3prime(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    A3PRIME evaluates A' function #3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double A3PRIME, the value of A'(X).
        //
    {
        double value = x switch
        {
            <= 1.0 / 3.0 => 2.0 * x,
            _ => 1.0
        };

        return value;
    }

    private static double c1(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C1 evaluates C function #1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double C1, the value of C(X).
        //
    {
        double value = 0.0;

        return value;
    }

    private static double c2(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C2 evaluates C function #2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double C2, the value of C(X).
        //
    {
        double value = 2.0;

        return value;
    }

    private static double c3(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C3 evaluates C function #3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double C3, the value of C(X).
        //
    {
        double value = 2.0 * x;

        return value;
    }

    private static double f1(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F1 evaluates right hand side function #1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double F1, the value of F(X).
        //
    {
        double value = x * (x + 3.0) * Math.Exp(x);

        return value;
    }

    private static double f2(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F2 evaluates right hand side function #2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double F2, the value of F(X).
        //
    {
        double value = x * (5.0 - x) * Math.Exp(x);

        return value;
    }

    private static double f3(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F3 evaluates right hand side function #3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double F3, the value of F(X).
        //
    {
        double value = -x * (2.0 * x * x - 3.0 * x - 3.0) * Math.Exp(x);

        return value;
    }

    private static double f4(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F4 evaluates right hand side function #4.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double F4, the value of F(X).
        //
    {
        double value = (x + 3.0 * x * x + 5.0 * x * x * x + x * x * x * x) * Math.Exp(x);

        return value;
    }

    private static double f5(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F5 evaluates right hand side function #5.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double F5, the value of F(X).
        //
    {
        double value = x switch
        {
            <= 1.0 / 3.0 => (x + 3.0 * x * x + 5.0 * x * x * x + x * x * x * x) * Math.Exp(x),
            _ => (-1.0 + 10.0 / 3.0 * x + 43.0 / 9.0 * x * x + x * x * x) * Math.Exp(x)
        };

        return value;
    }

    private static double[] exact1(int n, double[] x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT1 evaluates exact solution #1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input, double X[N], the evaluation points.
        //
        //    Output, double EXACT1[N], the values of U(X(1:N)).
        //
    {
        double[] uexact = new double[n];

        for (int i = 0; i < n; i++)
        {
            uexact[i] = x[i] * (1.0 - x[i]) * Math.Exp(x[i]);
        }

        return uexact;
    }

    private static double[] exact2(int n, double[] x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT2 returns exact solution #2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input, double X[N], the evaluation points.
        //
        //    Output, double EXACT2[N], the values of U(X(1:N)).
        //
    {
        double[] uexact = new double[n];

        for (int i = 0; i < n; i++)
        {
            uexact[i] = x[i] switch
            {
                <= 2.0 / 3.0 => x[i] * (1.0 - x[i]) * Math.Exp(x[i]),
                _ => x[i] * (1.0 - x[i]) * Math.Exp(2.0 / 3.0)
            };
        }

        return uexact;
    }

    private double[] exact3(int n, double[] x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT3 returns exact solution #3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input, double X[N], the evaluation points.
        //
        //    Output, double EXACT3[N], the values of U(X(1:N)).
        //
    {
        double[] uexact = new double[n];

        for (int i = 0; i < n; i++)
        {
            uexact[i] = x[i] switch
            {
                <= 2.0 / 3.0 => x[i] * (1.0 - x[i]) * Math.Exp(x[i]),
                _ => x[i] * (1.0 - x[i])
            };
        }

        return uexact;
    }
}