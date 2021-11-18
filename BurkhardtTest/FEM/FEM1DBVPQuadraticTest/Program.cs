using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.FEM;
using Burkardt.Types;

namespace FEM1DBVPQuadraticTest;

internal class Program
{
    private static void Main(string[] args)
    {
        Console.WriteLine("");
        Console.WriteLine("FEM1D_BVP_QUADRATIC_TEST");
        Console.WriteLine("  Test the FEM1D_BVP_QUADRATIC library.");

        test00();
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
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("FEM1D_BVP_QUADRATIC_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test00()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST00 carries out test case #0.
        //
        //  Discussion:
        //
        //    - uxx + u = x  for 0 < x < 1
        //    u(0) = u(1) = 0
        //
        //    exact  = x - sinh(x) / sinh(1)
        //    exact' = 1 - cosh(x) / sinh(1)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Dianne O'Leary,
        //    Scientific Computing with Case Studies,
        //    SIAM, 2008,
        //    ISBN13: 978-0-898716-66-5,
        //    LC: QA401.O44.
        //
    {
        int i;
        int n = 11;
        double e1;
        double e2;
        double h1s;
        double mx;
        double[] u;
        double uexact;
        double[] x;
        double x_first;
        double x_last;

        Console.WriteLine("");
        Console.WriteLine("TEST00");
        Console.WriteLine("  Solve -( A(x) U'(x) )' + C(x) U(x) = F(x)");
        Console.WriteLine("  for 0 < x < 1, with U(0) = U(1) = 0.");
        Console.WriteLine("  A(X)  = 1.0");
        Console.WriteLine("  C(X)  = 1.0");
        Console.WriteLine("  F(X)  = X");
        Console.WriteLine("  U(X)  = X - SINH(X) / SINH(1)");
        Console.WriteLine("");
        Console.WriteLine("  Number of nodes = " + n + "");
        //
        //  Geometry definitions.
        //
        x_first = 0.0;
        x_last = 1.0;
        x = typeMethods.r8vec_linspace_new(n, x_first, x_last);

        u = FEM_1D_BVP_Quadratic.fem1d_bvp_quadratic(n, a00, c00, f00, x);

        Console.WriteLine("");
        Console.WriteLine("     I    X         U         Uexact    Error");
        Console.WriteLine("");

        for (i = 0; i < n; i++)
        {
            uexact = exact00(x[i]);
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + uexact.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + Math.Abs(u[i] - uexact).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        e1 = FEM_Error.l1_error(n, x, u, exact00);
        e2 = FEM_Error.l2_error_quadratic(n, x, u, exact00);
        h1s = FEM_Error.h1s_error_quadratic(n, x, u, exact_ux00);
        mx = FEM_Error.max_error_quadratic(n, x, u, exact00);
        Console.WriteLine("");
        Console.WriteLine("  l1 norm of error  = " + e1 + "");
        Console.WriteLine("  L2 norm of error  = " + e2 + "");
        Console.WriteLine("  Seminorm of error = " + h1s + "");
        Console.WriteLine("  Max norm of error = " + mx + "");

    }

    private static double a00(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    A00 evaluates A function #0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double A00, the value of A(X).
        //
    {
        double value = 0;

        value = 1.0;

        return value;
    }

    private static double c00(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C00 evaluates C function #0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double C00, the value of C(X).
        //
    {
        double value = 0;

        value = 1.0;

        return value;
    }

    private static double exact00(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT00 evaluates exact solution #00.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT00, the value of U(X).
        //
    {
        double value = 0;

        value = x - Math.Sinh(x) / Math.Sinh(1.0);

        return value;
    }

    private static double exact_ux00(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT_UX00 evaluates the derivative of exact solution #00.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT_UX00, the value of dUdX(X).
        //
    {
        double value = 0;

        value = 1.0 - Math.Cosh(x) / Math.Sinh(1.0);

        return value;
    }

    private static double f00(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F00 evaluates right hand side function #00.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double F00, the value of F(X).
        //
    {
        double value = 0;

        value = x;

        return value;
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 carries out test case #1.
        //
        //  Discussion:
        //
        //    Use A1, C1, F1, EXACT1, EXACT_UX1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Dianne O'Leary,
        //    Scientific Computing with Case Studies,
        //    SIAM, 2008,
        //    ISBN13: 978-0-898716-66-5,
        //    LC: QA401.O44.
        //
    {
        int i;
        int n = 11;
        double e1;
        double e2;
        double h1s;
        double mx;
        double[] u;
        double uexact;
        double[] x;
        double x_first;
        double x_last;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Solve -( A(x) U'(x) )' + C(x) U(x) = F(x)");
        Console.WriteLine("  for 0 < x < 1, with U(0) = U(1) = 0.");
        Console.WriteLine("  A1(X)  = 1.0");
        Console.WriteLine("  C1(X)  = 0.0");
        Console.WriteLine("  F1(X)  = X * ( X + 3 ) * Math.Exp ( X )");
        Console.WriteLine("  U1(X)  = X * ( 1 - X ) * Math.Exp ( X )");
        Console.WriteLine("");
        Console.WriteLine("  Number of nodes = " + n + "");
        //
        //  Geometry definitions.
        //
        x_first = 0.0;
        x_last = 1.0;
        x = typeMethods.r8vec_linspace_new(n, x_first, x_last);

        u = FEM_1D_BVP_Quadratic.fem1d_bvp_quadratic(n, a1, c1, f1, x);

        Console.WriteLine("");
        Console.WriteLine("     I    X         U         Uexact    Error");
        Console.WriteLine("");

        for (i = 0; i < n; i++)
        {
            uexact = exact1(x[i]);
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + uexact.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + Math.Abs(u[i] - uexact).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        e1 = FEM_Error.l1_error(n, x, u, exact1);
        e2 = FEM_Error.l2_error_quadratic(n, x, u, exact1);
        h1s = FEM_Error.h1s_error_quadratic(n, x, u, exact_ux1);
        mx = FEM_Error.max_error_quadratic(n, x, u, exact1);
        Console.WriteLine("");
        Console.WriteLine("  l1 norm of error  = " + e1 + "");
        Console.WriteLine("  L2 norm of error  = " + e2 + "");
        Console.WriteLine("  Seminorm of error = " + h1s + "");
        Console.WriteLine("  Max norm of error = " + mx + "");

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
        //    14 June 2014
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
        double value = 0;

        value = 1.0;

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
        //    20 August 2010
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
        double value = 0;

        value = 0.0;

        return value;
    }

    private static double exact1(double x)

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
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT1, the value of U(X).
        //
    {
        double value = 0;

        value = x * (1.0 - x) * Math.Exp(x);

        return value;
    }

    private static double exact_ux1(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT_UX1 evaluates the derivative of exact solution #1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT_UX1, the value of dUdX(X).
        //
    {
        double value = 0;

        value = (1.0 - x - x * x) * Math.Exp(x);

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
        //    20 August 2010
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
        double value = 0;

        value = x * (x + 3.0) * Math.Exp(x);

        return value;
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 carries out test case #2.
        //
        //  Discussion:
        //
        //    Use A2, C2, F2, EXACT2, EXACT_UX2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Dianne O'Leary,
        //    Scientific Computing with Case Studies,
        //    SIAM, 2008,
        //    ISBN13: 978-0-898716-66-5,
        //    LC: QA401.O44.
        //
    {
        int i;
        int n = 11;
        double e1;
        double e2;
        double h1s;
        double mx;
        double[] u;
        double uexact;
        double[] x;
        double x_first;
        double x_last;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Solve -( A(x) U'(x) )' + C(x) U(x) = F(x)");
        Console.WriteLine("  for 0 < x < 1, with U(0) = U(1) = 0.");
        Console.WriteLine("  A2(X)  = 1.0");
        Console.WriteLine("  C2(X)  = 2.0");
        Console.WriteLine("  F2(X)  = X * ( 5 - X ) * Math.Exp ( X )");
        Console.WriteLine("  U2(X)  = X * ( 1 - X ) * Math.Exp ( X )");
        Console.WriteLine("");
        Console.WriteLine("  Number of nodes = " + n + "");
        //
        //  Geometry definitions.
        //
        x_first = 0.0;
        x_last = 1.0;
        x = typeMethods.r8vec_linspace_new(n, x_first, x_last);

        u = FEM_1D_BVP_Quadratic.fem1d_bvp_quadratic(n, a2, c2, f2, x);

        Console.WriteLine("");
        Console.WriteLine("     I    X         U         Uexact    Error");
        Console.WriteLine("");

        for (i = 0; i < n; i++)
        {
            uexact = exact2(x[i]);
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + uexact.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + Math.Abs(u[i] - uexact).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        e1 = FEM_Error.l1_error(n, x, u, exact2);
        e2 = FEM_Error.l2_error_quadratic(n, x, u, exact2);
        h1s = FEM_Error.h1s_error_quadratic(n, x, u, exact_ux2);
        mx = FEM_Error.max_error_quadratic(n, x, u, exact2);
        Console.WriteLine("");
        Console.WriteLine("  l1 norm of error  = " + e1 + "");
        Console.WriteLine("  L2 norm of error  = " + e2 + "");
        Console.WriteLine("  Seminorm of error = " + h1s + "");
        Console.WriteLine("  Max norm of error = " + mx + "");

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
        //    14 June 2014
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
        double value = 0;

        value = 1.0;

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
        //    20 August 2010
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
        double value = 0;

        value = 2.0;

        return value;
    }

    private static double exact2(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT2 evaluates exact solution #2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT2, the value of U(X).
        //
    {
        double value = 0;

        value = x * (1.0 - x) * Math.Exp(x);

        return value;
    }

    private static double exact_ux2(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT_UX2 evaluates the derivative of exact solution #2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT_UX2, the value of dUdX(X).
        //
    {
        double value = 0;

        value = (1.0 - x - x * x) * Math.Exp(x);

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
        //    20 August 2010
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
        double value = 0;

        value = x * (5.0 - x) * Math.Exp(x);

        return value;
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 carries out test case #3.
        //
        //  Discussion:
        //
        //    Use A3, C3, F3, EXACT3, EXACT_UX3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Dianne O'Leary,
        //    Scientific Computing with Case Studies,
        //    SIAM, 2008,
        //    ISBN13: 978-0-898716-66-5,
        //    LC: QA401.O44.
        //
    {
        int i;
        int n = 11;
        double e1;
        double e2;
        double h1s;
        double mx;
        double[] u;
        double uexact;
        double[] x;
        double x_first;
        double x_last;

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  Solve -( A(x) U'(x) )' + C(x) U(x) = F(x)");
        Console.WriteLine("  for 0 < x < 1, with U(0) = U(1) = 0.");
        Console.WriteLine("  A3(X)  = 1.0");
        Console.WriteLine("  C3(X)  = 2.0 * X");
        Console.WriteLine("  F3(X)  = - X * ( 2 * X * X - 3 * X - 3 ) * Math.Exp ( X )");
        Console.WriteLine("  U3(X)  = X * ( 1 - X ) * Math.Exp ( X )");
        Console.WriteLine("");
        Console.WriteLine("  Number of nodes = " + n + "");
        //
        //  Geometry definitions.
        //
        x_first = 0.0;
        x_last = 1.0;
        x = typeMethods.r8vec_linspace_new(n, x_first, x_last);

        u = FEM_1D_BVP_Quadratic.fem1d_bvp_quadratic(n, a3, c3, f3, x);

        Console.WriteLine("");
        Console.WriteLine("     I    X         U         Uexact    Error");
        Console.WriteLine("");

        for (i = 0; i < n; i++)
        {
            uexact = exact3(x[i]);
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + uexact.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + Math.Abs(u[i] - uexact).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        e1 = FEM_Error.l1_error(n, x, u, exact3);
        e2 = FEM_Error.l2_error_quadratic(n, x, u, exact3);
        h1s = FEM_Error.h1s_error_quadratic(n, x, u, exact_ux3);
        mx = FEM_Error.max_error_quadratic(n, x, u, exact3);
        Console.WriteLine("");
        Console.WriteLine("  l1 norm of error  = " + e1 + "");
        Console.WriteLine("  L2 norm of error  = " + e2 + "");
        Console.WriteLine("  Seminorm of error = " + h1s + "");
        Console.WriteLine("  Max norm of error = " + mx + "");

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
        //    14 June 2014
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
        double value = 0;

        value = 1.0;

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
        //    20 August 2010
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
        double value = 0;

        value = 2.0 * x;

        return value;
    }

    private static double exact3(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT3 evaluates exact solution #3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT3, the value of U(X).
        //
    {
        double value = 0;

        value = x * (1.0 - x) * Math.Exp(x);

        return value;
    }

    private static double exact_ux3(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT_UX3 evaluates the derivative of exact solution #3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT_UX3, the value of dUdX(X).
        //
    {
        double value = 0;

        value = (1.0 - x - x * x) * Math.Exp(x);

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
        //    20 August 2010
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
        double value = 0;

        value = -x * (2.0 * x * x - 3.0 * x - 3.0) * Math.Exp(x);

        return value;
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 carries out test case #4.
        //
        //  Discussion:
        //
        //    Use A4, C4, F4, EXACT4, EXACT_UX4.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Dianne O'Leary,
        //    Scientific Computing with Case Studies,
        //    SIAM, 2008,
        //    ISBN13: 978-0-898716-66-5,
        //    LC: QA401.O44.
        //
    {
        int i;
        int n = 11;
        double e1;
        double e2;
        double h1s;
        double mx;
        double[] u;
        double uexact;
        double[] x;
        double x_first;
        double x_last;

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  Solve -( A(x) U'(x) )' + C(x) U(x) = F(x)");
        Console.WriteLine("  for 0 < x < 1, with U(0) = U(1) = 0.");
        Console.WriteLine("  A4(X)  = 1.0 + X * X");
        Console.WriteLine("  C4(X)  = 0.0");
        Console.WriteLine("  F4(X)  = ( X + 3 X^2 + 5 X^3 + X^4 ) * Math.Exp ( X )");
        Console.WriteLine("  U4(X)  = X * ( 1 - X ) * Math.Exp ( X )");
        Console.WriteLine("");
        Console.WriteLine("  Number of nodes = " + n + "");
        //
        //  Geometry definitions.
        //
        x_first = 0.0;
        x_last = 1.0;
        x = typeMethods.r8vec_linspace_new(n, x_first, x_last);

        u = FEM_1D_BVP_Quadratic.fem1d_bvp_quadratic(n, a4, c4, f4, x);

        Console.WriteLine("");
        Console.WriteLine("     I    X         U         Uexact    Error");
        Console.WriteLine("");

        for (i = 0; i < n; i++)
        {
            uexact = exact4(x[i]);
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + uexact.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + Math.Abs(u[i] - uexact).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        e1 = FEM_Error.l1_error(n, x, u, exact4);
        e2 = FEM_Error.l2_error_quadratic(n, x, u, exact4);
        h1s = FEM_Error.h1s_error_quadratic(n, x, u, exact_ux4);
        mx = FEM_Error.max_error_quadratic(n, x, u, exact4);
        Console.WriteLine("");
        Console.WriteLine("  l1 norm of error  = " + e1 + "");
        Console.WriteLine("  L2 norm of error  = " + e2 + "");
        Console.WriteLine("  Seminorm of error = " + h1s + "");
        Console.WriteLine("  Max norm of error = " + mx + "");

    }

    private static double a4(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    A4 evaluates A function #4.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double A4, the value of A(X).
        //
    {
        double value = 0;

        value = 1.0 + x * x;

        return value;
    }

    private static double c4(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C4 evaluates C function #4.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double C4, the value of C(X).
        //
    {
        double value = 0;

        value = 0.0;

        return value;
    }

    private static double exact4(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT4 evaluates exact solution #4.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT4, the value of U(X).
        //
    {
        double value = 0;

        value = x * (1.0 - x) * Math.Exp(x);

        return value;
    }

    private static double exact_ux4(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT_UX4 evaluates the derivative of exact solution #4.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT_UX4, the value of dUdX(X).
        //
    {
        double value = 0;

        value = (1.0 - x - x * x) * Math.Exp(x);

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
        //    20 August 2010
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
        double value = 0;

        value = (x + 3.0 * x * x + 5.0 * x * x * x + x * x * x * x) * Math.Exp(x);

        return value;
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 carries out test case #5.
        //
        //  Discussion:
        //
        //    Use A5, C5, F5, EXACT5, EXACT_UX5.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Dianne O'Leary,
        //    Scientific Computing with Case Studies,
        //    SIAM, 2008,
        //    ISBN13: 978-0-898716-66-5,
        //    LC: QA401.O44.
        //
    {
        int i;
        int n = 11;
        double e1;
        double e2;
        double h1s;
        double mx;
        double[] u;
        double uexact;
        double[] x;
        double x_first;
        double x_last;

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  Solve -( A(x) U'(x) )' + C(x) U(x) = F(x)");
        Console.WriteLine("  for 0 < x < 1, with U(0) = U(1) = 0.");
        Console.WriteLine("  A5(X)  = 1.0 + X * X for X <= 1/3");
        Console.WriteLine("         = 7/9 + X     for      1/3 < X");
        Console.WriteLine("  C5(X)  = 0.0");
        Console.WriteLine("  F5(X)  = ( X + 3 X^2 + 5 X^3 + X^4 ) * Math.Exp ( X )");
        Console.WriteLine("                       for X <= 1/3");
        Console.WriteLine("         = ( - 1 + 10/3 X + 43/9 X^2 + X^3 ) .* Math.Exp ( X )");
        Console.WriteLine("                       for      1/3 <= X");
        Console.WriteLine("  U5(X)  = X * ( 1 - X ) * Math.Exp ( X )");
        Console.WriteLine("");
        Console.WriteLine("  Number of nodes = " + n + "");
        //
        //  Geometry definitions.
        //
        x_first = 0.0;
        x_last = 1.0;
        x = typeMethods.r8vec_linspace_new(n, x_first, x_last);

        u = FEM_1D_BVP_Quadratic.fem1d_bvp_quadratic(n, a5, c5, f5, x);

        Console.WriteLine("");
        Console.WriteLine("     I    X         U         Uexact    Error");
        Console.WriteLine("");

        for (i = 0; i < n; i++)
        {
            uexact = exact5(x[i]);
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + uexact.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + Math.Abs(u[i] - uexact).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        e1 = FEM_Error.l1_error(n, x, u, exact5);
        e2 = FEM_Error.l2_error_quadratic(n, x, u, exact5);
        h1s = FEM_Error.h1s_error_quadratic(n, x, u, exact_ux5);
        mx = FEM_Error.max_error_quadratic(n, x, u, exact5);
        Console.WriteLine("");
        Console.WriteLine("  l1 norm of error  = " + e1 + "");
        Console.WriteLine("  L2 norm of error  = " + e2 + "");
        Console.WriteLine("  Seminorm of error = " + h1s + "");
        Console.WriteLine("  Max norm of error = " + mx + "");

    }

    private static double a5(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    A5 evaluates A function #5.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double A5, the value of A(X).
        //
    {
        double value = x switch
        {
            <= 1.0 / 3.0 => 1.0 + x * x,
            _ => x + 7.0 / 9.0
        };

        return value;
    }

    private static double c5(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C5 evaluates C function #5.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double C5, the value of C(X).
        //
    {
        double value = 0;

        value = 0.0;

        return value;
    }

    private static double exact5(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT5 evaluates exact solution #5.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT5, the value of U(X).
        //
    {
        double value = 0;

        value = x * (1.0 - x) * Math.Exp(x);

        return value;
    }

    private static double exact_ux5(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT_UX5 evaluates the derivative of exact solution #5.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT_UX5, the value of dUdX(X).
        //
    {
        double value = 0;

        value = (1.0 - x - x * x) * Math.Exp(x);

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
        //    20 August 2010
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

    private static void test06()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 does an error analysis.
        //
        //  Discussion:
        //
        //    Use A6, C6, F6, EXACT6, EXACT_UX6.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Dianne O'Leary,
        //    Scientific Computing with Case Studies,
        //    SIAM, 2008,
        //    ISBN13: 978-0-898716-66-5,
        //    LC: QA401.O44.
        //
    {
        int i;
        int n;
        double e1;
        double e2;
        double h1s;
        double mx;
        double[] u;
        double[] x;
        double x_first;
        double x_last;

        Console.WriteLine("");
        Console.WriteLine("TEST06");
        Console.WriteLine("  Solve -( A(x) U'(x) )' + C(x) U(x) = F(x)");
        Console.WriteLine("  for 0 < x < 1, with U(0) = U(1) = 0.");
        Console.WriteLine("  A6(X)  = 1.0");
        Console.WriteLine("  C6(X)  = 0.0");
        Console.WriteLine("  F6(X)  = pi*pi*sin(pi*X)");
        Console.WriteLine("  U6(X)  = sin(pi*x)");
        Console.WriteLine("");
        Console.WriteLine("  Compute L2 norm and seminorm of error for various N.");
        Console.WriteLine("");
        Console.WriteLine("     N        l1 error      L2 error      Seminorm error  Maxnorm error");
        Console.WriteLine("");

        n = 11;
        for (i = 0; i <= 4; i++)
        {
            //
            //  Geometry definitions.
            //
            x_first = 0.0;
            x_last = 1.0;
            x = typeMethods.r8vec_linspace_new(n, x_first, x_last);

            u = FEM_1D_BVP_Quadratic.fem1d_bvp_quadratic(n, a6, c6, f6, x);

            e1 = FEM_Error.l1_error(n, x, u, exact6);
            e2 = FEM_Error.l2_error_quadratic(n, x, u, exact6);
            h1s = FEM_Error.h1s_error_quadratic(n, x, u, exact_ux6);
            mx = FEM_Error.max_error_quadratic(n, x, u, exact6);

            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + e1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + e2.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + h1s.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + mx.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            n = 2 * (n - 1) + 1;
        }
    }

    private static double a6(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    A6 evaluates A function #6.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double A6, the value of A(X).
        //
    {
        double value = 0;

        value = 1.0;

        return value;
    }

    private static double c6(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C6 evaluates C function #6.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double C6, the value of C(X).
        //
    {
        double value = 0;

        value = 0.0;

        return value;
    }

    private static double exact6(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT6 returns exact solution #6.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT6, the value of U(X).
        //
    {
        double value = 0;

        value = Math.Sin(Math.PI * x);

        return value;
    }

    private static double exact_ux6(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT_UX6 returns the derivative of exact solution #6.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT_UX6, the value of U(X).
        //
    {
        double value = 0;

        value = Math.PI * Math.Cos(Math.PI * x);

        return value;
    }

    private static double f6(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F6 evaluates right hand side function #6.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double F6, the value of F(X).
        //
    {
        double value = 0;

        value = Math.PI * Math.PI * Math.Sin(Math.PI * x);

        return value;
    }

    private static void test07()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST07 does an error analysis.
        //
        //  Discussion:
        //
        //    Use A7, C7, F7, EXACT7, EXACT_UX7.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Eric Becker, Graham Carey, John Oden,
        //    Finite Elements, An Introduction, Volume I,
        //    Prentice-Hall, 1981, page 123-124,
        //    ISBN: 0133170578,
        //    LC: TA347.F5.B4.
        //
    {
        int i;
        int n;
        double e1;
        double e2;
        double h1s;
        double mx;
        double[] u;
        double[] x;
        double x_first;
        double x_last;

        Console.WriteLine("");
        Console.WriteLine("TEST07");
        Console.WriteLine("  Solve -( A(x) U'(x) )' + C(x) U(x) = F(x)");
        Console.WriteLine("  for 0 < x < 1, with U(0) = U(1) = 0.");
        Console.WriteLine("  Becker/Carey/Oden example");
        Console.WriteLine("");
        Console.WriteLine("  Compute L2 norm and seminorm of error for various N.");
        Console.WriteLine("");
        Console.WriteLine("     N        l1 error      L2 error      Seminorm error  Maxnorm error");
        Console.WriteLine("");

        n = 11;
        for (i = 0; i <= 4; i++)
        {
            //
            //  Geometry definitions.
            //
            x_first = 0.0;
            x_last = 1.0;
            x = typeMethods.r8vec_linspace_new(n, x_first, x_last);

            u = FEM_1D_BVP_Quadratic.fem1d_bvp_quadratic(n, a7, c7, f7, x);

            e1 = FEM_Error.l1_error(n, x, u, exact7);
            e2 = FEM_Error.l2_error_quadratic(n, x, u, exact7);
            h1s = FEM_Error.h1s_error_quadratic(n, x, u, exact_ux7);
            mx = FEM_Error.max_error_quadratic(n, x, u, exact7);

            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + e1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + e2.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + mx.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + h1s.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            n = 2 * (n - 1) + 1;
        }

    }

    private static double a7(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    A7 evaluates A function #7.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double A7, the value of A(X).
        //
    {
        double alpha;
        double value = 0;
        double x0;

        alpha = 30.0;
        x0 = 1.0 / 3.0;
        value = 1.0 / alpha + alpha * Math.Pow(x - x0, 2);

        return value;
    }

    private static double c7(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C7 evaluates C function #7.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double C7, the value of C(X).
        //
    {
        double value = 0;

        value = 0.0;

        return value;
    }

    private static double exact7(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT7 returns exact solution #7.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT7, the value of U(X).
        //
    {
        double alpha;
        double value = 0;
        double x0;

        alpha = 30.0;
        x0 = 1.0 / 3.0;
        value = (1.0 - x)
                * (Math.Atan(alpha * (x - x0)) + Math.Atan(alpha * x0));

        return value;
    }

    private static double exact_ux7(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT_UX7 returns the derivative of exact solution #7.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT_UX7, the value of U(X).
        //
    {
        double alpha;
        double value = 0;
        double x0;

        alpha = 30.0;
        x0 = 1.0 / 3.0;
        value = -Math.Atan(alpha * (x - x0)) - Math.Atan(alpha * x0)
                + (1.0 - x) * alpha / (1.0 + alpha * alpha * Math.Pow(x - x0, 2));


        return value;
    }

    private static double f7(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F7 evaluates right hand side function #7.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double F7, the value of F(X).
        //
    {
        double alpha;
        double value = 0;
        double x0;

        alpha = 30.0;
        x0 = 1.0 / 3.0;
        value = 2.0 * (1.0 + alpha * (x - x0) *
            (Math.Atan(alpha * (x - x0)) + Math.Atan(alpha * x0)));

        return value;
    }

    private static void test08()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST08 carries out test case #8.
        //
        //  Discussion:
        //
        //    Use A8, C8, F8, EXACT8, EXACT_UX8.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Dianne O'Leary,
        //    Scientific Computing with Case Studies,
        //    SIAM, 2008,
        //    ISBN13: 978-0-898716-66-5,
        //    LC: QA401.O44.
        //
    {
        int i;
        int n = 11;
        double e1;
        double e2;
        double h1s;
        double mx;
        double[] u;
        double uexact;
        double[] x;
        double x_first;
        double x_last;

        Console.WriteLine("");
        Console.WriteLine("TEST08");
        Console.WriteLine("  Solve -( A(x) U'(x) )' + C(x) U(x) = F(x)");
        Console.WriteLine("  for 0 < x < 1, with U(0) = U(1) = 0.");
        Console.WriteLine("  A8(X)  = 1.0");
        Console.WriteLine("  C8(X)  = 0.0");
        Console.WriteLine("  F8(X) = X * ( X + 3 ) * Math.Exp ( X ),   X <= 2/3");
        Console.WriteLine("        = 2 * Math.Exp ( 2/3),                   2/3 < X");
        Console.WriteLine("  U8(X) = X * ( 1 - X ) * Math.Exp ( X ),   X <= 2/3");
        Console.WriteLine("        = X * ( 1 - X ) * Math.Exp ( 2/3 ),      2/3 < X");
        Console.WriteLine("");
        Console.WriteLine("  Number of nodes = " + n + "");
        //
        //  Geometry definitions.
        //
        x_first = 0.0;
        x_last = 1.0;
        x = typeMethods.r8vec_linspace_new(n, x_first, x_last);

        u = FEM_1D_BVP_Quadratic.fem1d_bvp_quadratic(n, a8, c8, f8, x);

        Console.WriteLine("");
        Console.WriteLine("     I    X         U         Uexact    Error");
        Console.WriteLine("");

        for (i = 0; i < n; i++)
        {
            uexact = exact8(x[i]);
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + uexact.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + Math.Abs(u[i] - uexact).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        e1 = FEM_Error.l1_error(n, x, u, exact8);
        e2 = FEM_Error.l2_error_quadratic(n, x, u, exact8);
        h1s = FEM_Error.h1s_error_quadratic(n, x, u, exact_ux8);
        mx = FEM_Error.max_error_quadratic(n, x, u, exact8);
        Console.WriteLine("");
        Console.WriteLine("  l1 norm of error  = " + e1 + "");
        Console.WriteLine("  L2 norm of error  = " + e2 + "");
        Console.WriteLine("  Seminorm of error = " + h1s + "");
        Console.WriteLine("  Max norm of error = " + mx + "");
    }

    private static double a8(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    A8 evaluates A function #8.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double A8, the value of A(X).
        //
    {
        double value = 0;

        value = 1.0;

        return value;
    }

    private static double c8(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C8 evaluates C function #8.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double C8, the value of C(X).
        //
    {
        double value = 0;

        value = 0.0;

        return value;
    }

    private static double exact8(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT8 evaluates exact solution #8.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT8, the value of U(X).
        //
    {
        double value = x switch
        {
            <= 2.0 / 3.0 => x * (1.0 - x) * Math.Exp(x),
            _ => x * (1.0 - x) * Math.Exp(2.0 / 3.0)
        };

        return value;
    }

    private static double exact_ux8(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT_UX8 evaluates the derivative of exact solution #8.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT_UX8, the value of dUdX(X).
        //
    {
        double value = x switch
        {
            <= 2.0 / 3.0 => (1.0 - x - x * x) * Math.Exp(x),
            _ => (1.0 - 2.0 * x) * Math.Exp(2.0 / 3.0)
        };

        return value;
    }

    private static double f8(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F8 evaluates right hand side function #8.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double F8, the value of F(X).
        //
    {
        double value = x switch
        {
            <= 2.0 / 3.0 => x * (x + 3.0) * Math.Exp(x),
            _ => 2.0 * Math.Exp(2.0 / 3.0)
        };

        return value;
    }

    private static void test09()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST09 carries out test case #9.
        //
        //  Discussion:
        //
        //    Use A9, C9, F9, EXACT9, EXACT_UX9.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Dianne O'Leary,
        //    Scientific Computing with Case Studies,
        //    SIAM, 2008,
        //    ISBN13: 978-0-898716-66-5,
        //    LC: QA401.O44.
        //
    {
        int i;
        int n = 11;
        double e1;
        double e2;
        double h1s;
        double mx;
        double[] u;
        double uexact;
        double[] x;
        double x_first;
        double x_last;

        Console.WriteLine("");
        Console.WriteLine("TEST09");
        Console.WriteLine("  Solve -( A(x) U'(x) )' + C(x) U(x) = F(x)");
        Console.WriteLine("  for 0 < x < 1, with U(0) = U(1) = 0.");
        Console.WriteLine("  A9(X)  = 1.0");
        Console.WriteLine("  C9(X)  = 0.0");
        Console.WriteLine("  F9(X) = X * ( X + 3 ) * Math.Exp ( X ),   X <= 2/3");
        Console.WriteLine("        = 2 * Math.Exp ( 2/3),                   2/3 < X");
        Console.WriteLine("  U9(X) = X * ( 1 - X ) * Math.Exp ( X ),   X <= 2/3");
        Console.WriteLine("        = X * ( 1 - X ) * Math.Exp ( 2/3 ),      2/3 < X");
        Console.WriteLine("");
        Console.WriteLine("  Number of nodes = " + n + "");
        //
        //  Geometry definitions.
        //
        x_first = 0.0;
        x_last = 1.0;
        x = typeMethods.r8vec_linspace_new(n, x_first, x_last);

        u = FEM_1D_BVP_Quadratic.fem1d_bvp_quadratic(n, a9, c9, f9, x);

        Console.WriteLine("");
        Console.WriteLine("     I    X         U         Uexact    Error");
        Console.WriteLine("");

        for (i = 0; i < n; i++)
        {
            uexact = exact9(x[i]);
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + uexact.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + Math.Abs(u[i] - uexact).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        e1 = FEM_Error.l1_error(n, x, u, exact9);
        e2 = FEM_Error.l2_error_quadratic(n, x, u, exact9);
        h1s = FEM_Error.h1s_error_quadratic(n, x, u, exact_ux9);
        mx = FEM_Error.max_error_quadratic(n, x, u, exact9);
        Console.WriteLine("");
        Console.WriteLine("  l1 norm of error  = " + e1 + "");
        Console.WriteLine("  L2 norm of error  = " + e2 + "");
        Console.WriteLine("  Seminorm of error = " + h1s + "");
        Console.WriteLine("  Max norm of error = " + mx + "");
    }

    private static double a9(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    A9 evaluates A function #9.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double A8, the value of A(X).
        //
    {
        double value = 0;

        value = 1.0;

        return value;
    }

    private static double c9(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C9 evaluates C function #9.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double C9, the value of C(X).
        //
    {
        double value = 0;

        value = 0.0;

        return value;
    }

    private static double exact9(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT9 evaluates exact solution #9.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT9, the value of U(X).
        //
    {
        double value = x switch
        {
            <= 2.0 / 3.0 => x * (1.0 - x) * Math.Exp(x),
            _ => x * (1.0 - x)
        };

        return value;
    }

    private static double exact_ux9(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT_UX9 evaluates the derivative of exact solution #9.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT_UX9, the value of dUdX(X).
        //
    {
        double value = x switch
        {
            <= 2.0 / 3.0 => (1.0 - x - x * x) * Math.Exp(x),
            _ => 1.0 - 2.0 * x
        };

        return value;
    }

    private static double f9(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F9 evaluates right hand side function #9.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double F9, the value of F(X).
        //
    {
        double value = x switch
        {
            <= 2.0 / 3.0 => x * (x + 3.0) * Math.Exp(x),
            _ => 2.0
        };

        return value;
    }

    private static void test10()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST10 tests FEM1D_BVP_QUADRATIC.
        //
        //  Discussion:
        //
        //    We want to compute errors and do convergence rates for the 
        //    following problem:
        //
        //    - uxx + u = x  for 0 < x < 1
        //    u(0) = u(1) = 0
        //
        //    exact  = x - sinh(x) / sinh(1)
        //    exact' = 1 - cosh(x) / sinh(1)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Dianne O'Leary,
        //    Scientific Computing with Case Studies,
        //    SIAM, 2008,
        //    ISBN13: 978-0-898716-66-5,
        //    LC: QA401.O44.
        //
    {
        string command_filename;
        List<string> command_unit = new();
        string data_filename;
        List<string> data_unit = new();
        int e_log;
        int e_log_max = 6;
        double[] h_plot;
        double h1;
        double[] h1_plot;
        double l2;
        double[] l2_plot;
        double mx;
        double[] mx_plot;
        int n;
        int ne;
        int ne1;
        int ne2;
        int[] ne_plot;
        string output_filename;
        double r;
        double[] u;
        double[] x;
        double x_hi;
        double x_lo;

        Console.WriteLine("");
        Console.WriteLine("TEST10");
        Console.WriteLine("  Solve -( A(x) U'(x) )' + C(x) U(x) = F(x)");
        Console.WriteLine("  for 0 < x < 1, with U(0) = U(1) = 0.");
        Console.WriteLine("  A(X)  = 1.0");
        Console.WriteLine("  C(X)  = 1.0");
        Console.WriteLine("  F(X)  = X");
        Console.WriteLine("  U(X)  = X - SINH(X) / SINH(1)");
        Console.WriteLine("");
        Console.WriteLine(" log(E)    E         L2error         H1error        Maxerror");
        Console.WriteLine("");

        h_plot = new double[e_log_max + 1];
        h1_plot = new double[e_log_max + 1];
        l2_plot = new double[e_log_max + 1];
        mx_plot = new double[e_log_max + 1];
        ne_plot = new int[e_log_max + 1];

        for (e_log = 0; e_log <= e_log_max; e_log++)
        {
            ne = (int) Math.Pow(2, e_log + 1);

            n = ne + 1;
            x_lo = 0.0;
            x_hi = 1.0;
            x = typeMethods.r8vec_linspace_new(n, x_lo, x_hi);

            u = FEM_1D_BVP_Quadratic.fem1d_bvp_quadratic(n, a10, c10, f10, x);

            ne_plot[e_log] = ne;

            h_plot[e_log] = (x_hi - x_lo) / ne;

            l2 = FEM_Error.l2_error_quadratic(n, x, u, exact10);
            l2_plot[e_log] = l2;

            h1 = FEM_Error.h1s_error_quadratic(n, x, u, exact_ux10);
            h1_plot[e_log] = h1;

            mx = FEM_Error.max_error_quadratic(n, x, u, exact10);
            mx_plot[e_log] = mx;

            Console.WriteLine("  " + e_log.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + ne.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + l2.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + h1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + mx.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

        }

        Console.WriteLine("");
        Console.WriteLine(" log(E1)  E1 / E2          L2rate          H1rate         Maxrate");
        Console.WriteLine("");

        for (e_log = 0; e_log < e_log_max; e_log++)
        {
            ne1 = ne_plot[e_log];
            ne2 = ne_plot[e_log + 1];
            r = ne2 / (double) ne1;
            l2 = l2_plot[e_log] / l2_plot[e_log + 1];
            l2 = Math.Log(l2) / Math.Log(r);
            h1 = h1_plot[e_log] / h1_plot[e_log + 1];
            h1 = Math.Log(h1) / Math.Log(r);
            mx = mx_plot[e_log] / mx_plot[e_log + 1];
            mx = Math.Log(mx) / Math.Log(r);
            Console.WriteLine("  " + e_log.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + ne1.ToString(CultureInfo.InvariantCulture).PadLeft(4) + " /" + ne2.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + l2.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + h1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + mx.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        //
        //  Create the data file.
        //
        data_filename = "data.txt";
        for (e_log = 0; e_log <= e_log_max; e_log++)
        {
            data_unit.Add("  " + ne_plot[e_log]
                               + "  " + l2_plot[e_log]
                               + "  " + h1_plot[e_log]
                               + "  " + mx_plot[e_log] + "");
        }

        File.WriteAllLines(data_filename, data_unit);
        data_unit.Clear();
        Console.WriteLine("");
        Console.WriteLine("  Created graphics data file \"" + data_filename + "\".");
        //
        //  Plot the L2 error as a function of NE.
        //
        command_filename = "commands_l2.txt";

        output_filename = "l2.png";

        command_unit.Add("# " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("# Usage:");
        command_unit.Add("#  gnuplot < " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("set term png");
        command_unit.Add("set output '" + output_filename + "'");
        command_unit.Add("set xlabel '<---NE--->'");
        command_unit.Add("set ylabel '<---L2(NE)--->'");
        command_unit.Add("set title 'L2 error versus number of elements NE'");
        command_unit.Add("set logscale xy");
        command_unit.Add("set size ratio -1");
        command_unit.Add("set grid");
        command_unit.Add("set style data lines");
        command_unit.Add("plot '" + data_filename + "' using 1:2 with points pt 7 ps 2 lc rgb 'blue',\\");
        command_unit.Add("     '" + data_filename + "' using 1:2 lw 3 linecolor rgb 'red'");

        File.WriteAllLines(command_filename, command_unit);
        command_unit.Clear();
        Console.WriteLine("  Created graphics command file \"" + command_filename + "\".");
        //
        //  Plot the H1 error as a function of NE.
        //
        command_filename = "commands_h1.txt";

        output_filename = "h1.png";

        command_unit.Add("# " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("# Usage:");
        command_unit.Add("#  gnuplot < " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("set term png");
        command_unit.Add("set output '" + output_filename + "'");
        command_unit.Add("set xlabel '<---NE--->'");
        command_unit.Add("set ylabel '<---H1(NE)--->'");
        command_unit.Add("set title 'H1 error versus number of elements NE'");
        command_unit.Add("set logscale xy");
        command_unit.Add("set size ratio -1");
        command_unit.Add("set grid");
        command_unit.Add("set style data lines");
        command_unit.Add("plot '" + data_filename + "' using 1:3 with points pt 7 ps 2 lc rgb 'blue',\\");
        command_unit.Add("     '" + data_filename + "' using 1:3 lw 3 linecolor rgb 'red'");

        File.WriteAllLines(command_filename, command_unit);
        command_unit.Clear();
        Console.WriteLine("  Created graphics command file \"" + command_filename + "\".");
        //
        //  Plot the MX error as a function of NE.
        //
        command_filename = "commands_mx.txt";

        output_filename = "mx.png";

        command_unit.Add("# " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("# Usage:");
        command_unit.Add("#  gnuplot < " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("set term png");
        command_unit.Add("set output '" + output_filename + "'");
        command_unit.Add("set xlabel '<---NE--->'");
        command_unit.Add("set ylabel '<---MX(NE)--->'");
        command_unit.Add("set title 'Max error versus number of elements NE'");
        command_unit.Add("set logscale xy");
        command_unit.Add("set size ratio -1");
        command_unit.Add("set grid");
        command_unit.Add("set style data lines");
        command_unit.Add("plot '" + data_filename + "' using 1:4 with points pt 7 ps 2 lc rgb 'blue',\\");
        command_unit.Add("     '" + data_filename + "' using 1:4 lw 3 linecolor rgb 'red'");

        File.WriteAllLines(command_filename, command_unit);
        command_unit.Clear();
        Console.WriteLine("  Created graphics command file \"" + command_filename + "\".");

    }

    private static double a10(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    A10 evaluates A function #0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double A10, the value of A(X).
        //
    {
        double value = 0;

        value = 1.0;

        return value;
    }

    private static double c10(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C10 evaluates C function #10.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double C10, the value of C(X).
        //
    {
        double value = 0;

        value = 1.0;

        return value;
    }

    private static double exact10(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT10 evaluates exact solution #00.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT10, the value of U(X).
        //
    {
        double value = 0;

        value = x - Math.Sinh(x) / Math.Sinh(1.0);

        return value;
    }

    private static double exact_ux10(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT_UX10 evaluates the derivative of exact solution #10.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT_UX10, the value of dUdX(X).
        //
    {
        double value = 0;

        value = 1.0 - Math.Cosh(x) / Math.Sinh(1.0);

        return value;
    }

    private static double f10(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F10 evaluates right hand side function #10.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double F10, the value of F(X).
        //
    {
        double value = 0;

        value = x;

        return value;
    }
}