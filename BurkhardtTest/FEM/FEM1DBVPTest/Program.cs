using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using Burkardt.FEM;
using Burkardt.Types;

namespace FEM1DBVPTest;

internal static class Program
{
    private static void Main()
    {
//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FEM1D_BVP_LINEAR_TEST.
//
//  Location:
//
//    http://people.sc.fsu.edu/~jburkardt/cpp_src/fem1d_bvp_linear/fem1d_bvp_linear_prb.cpp
//
//  Discussion:
//
//    FEM1D_BVP_LINEAR_TEST tests the FEM1D_BVP_LINEAR library.
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
        {
            Console.WriteLine("");
            Console.WriteLine("FEM1D_BVP_LINEAR_TEST");
            Console.WriteLine("  Test the FEM1D_BVP_LINEAR library.");

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

            Console.WriteLine("");
            Console.WriteLine("FEM1D_BVP_LINEAR_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");

        }

        static void test00()
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
//    10 July 2015
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
            const int n = 11;

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
            double x_first = 0.0;
            double x_last = 1.0;
            double[] x = typeMethods.r8vec_linspace_new(n, x_first, x_last);

            double[] u = FEM_1D_BVP.fem1d_bvp_linear(n, FEM_Test_Methods.a00, FEM_Test_Methods.c00, FEM_Test_Methods.f00, x);

            Console.WriteLine("");
            Console.WriteLine("     I    X         U         Uexact    Error");
            Console.WriteLine("");

            for (i = 0; i < n; i++)
            {
                double uexact = FEM_Test_Methods.exact00(x[i]);
                Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                       + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + uexact.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + Math.Abs(u[i] - uexact).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }

            double e1 = FEM_Error.l1_error(n, x, u, FEM_Test_Methods.exact00);
            double e2 = FEM_Error.l2_error_linear(n, x, u, FEM_Test_Methods.exact00);
            double h1s = FEM_Error.h1s_error_linear(n, x, u, FEM_Test_Methods.exact_ux00);
            double mx = FEM_Error.max_error_linear(n, x, u, FEM_Test_Methods.exact00);
            Console.WriteLine("");
            Console.WriteLine("  l1 norm of error  = " + e1 + "");
            Console.WriteLine("  L2 norm of error  = " + e2 + "");
            Console.WriteLine("  Seminorm of error = " + h1s + "");
            Console.WriteLine("  Max norm of error = " + mx + "");

        }


        static void test01()
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
            const int n = 11;

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  Solve -( A(x) U'(x) )' + C(x) U(x) = F(x)");
            Console.WriteLine("  for 0 < x < 1, with U(0) = U(1) = 0.");
            Console.WriteLine("  A1(X)  = 1.0");
            Console.WriteLine("  C1(X)  = 0.0");
            Console.WriteLine("  F1(X)  = X * ( X + 3 ) * exp ( X )");
            Console.WriteLine("  U1(X)  = X * ( 1 - X ) * exp ( X )");
            Console.WriteLine("");
            Console.WriteLine("  Number of nodes = " + n + "");
//
//  Geometry definitions.
//
            double x_first = 0.0;
            double x_last = 1.0;
            double[] x = typeMethods.r8vec_linspace_new(n, x_first, x_last);

            double[] u = FEM_1D_BVP.fem1d_bvp_linear(n, FEM_Test_Methods.a1, FEM_Test_Methods.c1, FEM_Test_Methods.f1, x);

            Console.WriteLine("");
            Console.WriteLine("     I    X         U         Uexact    Error");
            Console.WriteLine("");

            for (i = 0; i < n; i++)
            {
                double uexact = FEM_Test_Methods.exact1(x[i]);
                Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                       + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + uexact.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + Math.Abs(u[i] - uexact).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }

            double e1 = FEM_Error.l1_error(n, x, u, FEM_Test_Methods.exact1);
            double e2 = FEM_Error.l2_error_linear(n, x, u, FEM_Test_Methods.exact1);
            double h1s = FEM_Error.h1s_error_linear(n, x, u, FEM_Test_Methods.exact_ux1);
            double mx = FEM_Error.max_error_linear(n, x, u, FEM_Test_Methods.exact1);
            Console.WriteLine("");
            Console.WriteLine("  l1 norm of error  = " + e1 + "");
            Console.WriteLine("  L2 norm of error  = " + e2 + "");
            Console.WriteLine("  Seminorm of error = " + h1s + "");
            Console.WriteLine("  Max norm of error = " + mx + "");

        }


        static void test02()
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
            const int n = 11;

            Console.WriteLine("");
            Console.WriteLine("TEST02");
            Console.WriteLine("  Solve -( A(x) U'(x) )' + C(x) U(x) = F(x)");
            Console.WriteLine("  for 0 < x < 1, with U(0) = U(1) = 0.");
            Console.WriteLine("  A2(X)  = 1.0");
            Console.WriteLine("  C2(X)  = 2.0");
            Console.WriteLine("  F2(X)  = X * ( 5 - X ) * exp ( X )");
            Console.WriteLine("  U2(X)  = X * ( 1 - X ) * exp ( X )");
            Console.WriteLine("");
            Console.WriteLine("  Number of nodes = " + n + "");
//
//  Geometry definitions.
//
            const double x_first = 0.0;
            const double x_last = 1.0;
            double[] x = typeMethods.r8vec_linspace_new(n, x_first, x_last);

            double[] u = FEM_1D_BVP.fem1d_bvp_linear(n, FEM_Test_Methods.a2, FEM_Test_Methods.c2, FEM_Test_Methods.f2, x);

            Console.WriteLine("");
            Console.WriteLine("     I    X         U         Uexact    Error");
            Console.WriteLine("");

            for (i = 0; i < n; i++)
            {
                double uexact = FEM_Test_Methods.exact2(x[i]);
                Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                       + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + uexact.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + Math.Abs(u[i] - uexact).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }

            double e1 = FEM_Error.l1_error(n, x, u, FEM_Test_Methods.exact2);
            double e2 = FEM_Error.l2_error_linear(n, x, u, FEM_Test_Methods.exact2);
            double h1s = FEM_Error.h1s_error_linear(n, x, u, FEM_Test_Methods.exact_ux2);
            double mx = FEM_Error.max_error_linear(n, x, u, FEM_Test_Methods.exact2);
            Console.WriteLine("");
            Console.WriteLine("  l1 norm of error  = " + e1 + "");
            Console.WriteLine("  L2 norm of error  = " + e2 + "");
            Console.WriteLine("  Seminorm of error = " + h1s + "");
            Console.WriteLine("  Max norm of error = " + mx + "");
        }


        static void test03()
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
            const int n = 11;

            Console.WriteLine("");
            Console.WriteLine("TEST03");
            Console.WriteLine("  Solve -( A(x) U'(x) )' + C(x) U(x) = F(x)");
            Console.WriteLine("  for 0 < x < 1, with U(0) = U(1) = 0.");
            Console.WriteLine("  A3(X)  = 1.0");
            Console.WriteLine("  C3(X)  = 2.0 * X");
            Console.WriteLine("  F3(X)  = - X * ( 2 * X * X - 3 * X - 3 ) * exp ( X )");
            Console.WriteLine("  U3(X)  = X * ( 1 - X ) * exp ( X )");
            Console.WriteLine("");
            Console.WriteLine("  Number of nodes = " + n + "");
//
//  Geometry definitions.
//
            const double x_first = 0.0;
            const double x_last = 1.0;
            double[] x = typeMethods.r8vec_linspace_new(n, x_first, x_last);

            double[] u = FEM_1D_BVP.fem1d_bvp_linear(n, FEM_Test_Methods.a3, FEM_Test_Methods.c3, FEM_Test_Methods.f3, x);

            Console.WriteLine("");
            Console.WriteLine("     I    X         U         Uexact    Error");
            Console.WriteLine("");

            for (i = 0; i < n; i++)
            {
                double uexact = FEM_Test_Methods.exact3(x[i]);
                Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                       + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + uexact.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + Math.Abs(u[i] - uexact).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }

            double e1 = FEM_Error.l1_error(n, x, u, FEM_Test_Methods.exact3);
            double e2 = FEM_Error.l2_error_linear(n, x, u, FEM_Test_Methods.exact3);
            double h1s = FEM_Error.h1s_error_linear(n, x, u, FEM_Test_Methods.exact_ux3);
            double mx = FEM_Error.max_error_linear(n, x, u, FEM_Test_Methods.exact3);
            Console.WriteLine("");
            Console.WriteLine("  l1 norm of error  = " + e1 + "");
            Console.WriteLine("  L2 norm of error  = " + e2 + "");
            Console.WriteLine("  Seminorm of error = " + h1s + "");
            Console.WriteLine("  Max norm of error = " + mx + "");
        }

        static void test04()
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
            const int n = 11;

            Console.WriteLine("");
            Console.WriteLine("TEST04");
            Console.WriteLine("  Solve -( A(x) U'(x) )' + C(x) U(x) = F(x)");
            Console.WriteLine("  for 0 < x < 1, with U(0) = U(1) = 0.");
            Console.WriteLine("  A4(X)  = 1.0 + X * X");
            Console.WriteLine("  C4(X)  = 0.0");
            Console.WriteLine("  F4(X)  = ( X + 3 X^2 + 5 X^3 + X^4 ) * exp ( X )");
            Console.WriteLine("  U4(X)  = X * ( 1 - X ) * exp ( X )");
            Console.WriteLine("");
            Console.WriteLine("  Number of nodes = " + n + "");
//
//  Geometry definitions.
//
            double x_first = 0.0;
            double x_last = 1.0;
            double[] x = typeMethods.r8vec_linspace_new(n, x_first, x_last);

            double[] u = FEM_1D_BVP.fem1d_bvp_linear(n, FEM_Test_Methods.a4, FEM_Test_Methods.c4, FEM_Test_Methods.f4, x);

            Console.WriteLine("");
            Console.WriteLine("     I    X         U         Uexact    Error");
            Console.WriteLine("");

            for (i = 0; i < n; i++)
            {
                double uexact = FEM_Test_Methods.exact4(x[i]);
                Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                       + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + uexact.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + Math.Abs(u[i] - uexact).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }

            double e1 = FEM_Error.l1_error(n, x, u, FEM_Test_Methods.exact4);
            double e2 = FEM_Error.l2_error_linear(n, x, u, FEM_Test_Methods.exact4);
            double h1s = FEM_Error.h1s_error_linear(n, x, u, FEM_Test_Methods.exact_ux4);
            double mx = FEM_Error.max_error_linear(n, x, u, FEM_Test_Methods.exact4);
            Console.WriteLine("");
            Console.WriteLine("  l1 norm of error  = " + e1 + "");
            Console.WriteLine("  L2 norm of error  = " + e2 + "");
            Console.WriteLine("  Seminorm of error = " + h1s + "");
            Console.WriteLine("  Max norm of error = " + mx + "");
        }


        static void test05()
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
            const int n = 11;

            Console.WriteLine("");
            Console.WriteLine("TEST05");
            Console.WriteLine("  Solve -( A(x) U'(x) )' + C(x) U(x) = F(x)");
            Console.WriteLine("  for 0 < x < 1, with U(0) = U(1) = 0.");
            Console.WriteLine("  A5(X)  = 1.0 + X * X for X <= 1/3");
            Console.WriteLine("         = 7/9 + X     for      1/3 < X");
            Console.WriteLine("  C5(X)  = 0.0");
            Console.WriteLine("  F5(X)  = ( X + 3 X^2 + 5 X^3 + X^4 ) * exp ( X )");
            Console.WriteLine("                       for X <= 1/3");
            Console.WriteLine("         = ( - 1 + 10/3 X + 43/9 X^2 + X^3 ) .* exp ( X )");
            Console.WriteLine("                       for      1/3 <= X");
            Console.WriteLine("  U5(X)  = X * ( 1 - X ) * exp ( X )");
            Console.WriteLine("");
            Console.WriteLine("  Number of nodes = " + n + "");
//
//  Geometry definitions.
//
            double x_first = 0.0;
            double x_last = 1.0;
            double[] x = typeMethods.r8vec_linspace_new(n, x_first, x_last);

            double[] u = FEM_1D_BVP.fem1d_bvp_linear(n, FEM_Test_Methods.a5, FEM_Test_Methods.c5, FEM_Test_Methods.f5, x);

            Console.WriteLine("");
            Console.WriteLine("     I    X         U         Uexact    Error");
            Console.WriteLine("");

            for (i = 0; i < n; i++)
            {
                double uexact = FEM_Test_Methods.exact5(x[i]);
                Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                       + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + uexact.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + Math.Abs(u[i] - uexact).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }

            double e1 = FEM_Error.l1_error(n, x, u, FEM_Test_Methods.exact5);
            double e2 = FEM_Error.l2_error_linear(n, x, u, FEM_Test_Methods.exact5);
            double h1s = FEM_Error.h1s_error_linear(n, x, u, FEM_Test_Methods.exact_ux5);
            double mx = FEM_Error.max_error_linear(n, x, u, FEM_Test_Methods.exact5);
            Console.WriteLine("");
            Console.WriteLine("  l1 norm of error  = " + e1 + "");
            Console.WriteLine("  L2 norm of error  = " + e2 + "");
            Console.WriteLine("  Seminorm of error = " + h1s + "");
            Console.WriteLine("  Max norm of error = " + mx + "");
        }

        static void test06()
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

            int n = 11;
            for (i = 0; i <= 4; i++)
            {
//
//  Geometry definitions.
//
                double x_first = 0.0;
                double x_last = 1.0;
                double[] x = typeMethods.r8vec_linspace_new(n, x_first, x_last);

                double[] u = FEM_1D_BVP.fem1d_bvp_linear(n, FEM_Test_Methods.a6, FEM_Test_Methods.c6, FEM_Test_Methods.f6,
                    x);

                double e1 = FEM_Error.l1_error(n, x, u, FEM_Test_Methods.exact6);
                double e2 = FEM_Error.l2_error_linear(n, x, u, FEM_Test_Methods.exact6);
                double h1s = FEM_Error.h1s_error_linear(n, x, u, FEM_Test_Methods.exact_ux6);
                double mx = FEM_Error.max_error_linear(n, x, u, FEM_Test_Methods.exact6);

                Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                       + "  " + e1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + e2.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + h1s.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + mx.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

                n = 2 * (n - 1) + 1;

            }
        }


        static void test07()

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

            int n = 11;
            for (i = 0; i <= 4; i++)
            {
//
//  Geometry definitions.
//
                double x_first = 0.0;
                double x_last = 1.0;
                double[] x = typeMethods.r8vec_linspace_new(n, x_first, x_last);

                double[] u = FEM_1D_BVP.fem1d_bvp_linear(n, FEM_Test_Methods.a7, FEM_Test_Methods.c7, FEM_Test_Methods.f7,
                    x);

                double e1 = FEM_Error.l1_error(n, x, u, FEM_Test_Methods.exact7);
                double e2 = FEM_Error.l2_error_linear(n, x, u, FEM_Test_Methods.exact7);
                double h1s = FEM_Error.h1s_error_linear(n, x, u, FEM_Test_Methods.exact_ux7);
                double mx = FEM_Error.max_error_linear(n, x, u, FEM_Test_Methods.exact7);

                Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                       + "  " + e1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + e2.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + h1s.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + mx.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

                n = 2 * (n - 1) + 1;

            }
        }

        static void test08()
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
            const int n = 11;

            Console.WriteLine("");
            Console.WriteLine("TEST08");
            Console.WriteLine("  Solve -( A(x) U'(x) )' + C(x) U(x) = F(x)");
            Console.WriteLine("  for 0 < x < 1, with U(0) = U(1) = 0.");
            Console.WriteLine("  A8(X)  = 1.0");
            Console.WriteLine("  C8(X)  = 0.0");
            Console.WriteLine("  F8(X) = X * ( X + 3 ) * exp ( X ),   X <= 2/3");
            Console.WriteLine("        = 2 * exp ( 2/3),                   2/3 < X");
            Console.WriteLine("  U8(X) = X * ( 1 - X ) * exp ( X ),   X <= 2/3");
            Console.WriteLine("        = X * ( 1 - X ) * exp ( 2/3 ),      2/3 < X");
            Console.WriteLine("");
            Console.WriteLine("  Number of nodes = " + n + "");
//
//  Geometry definitions.
//
            double x_first = 0.0;
            double x_last = 1.0;
            double[] x = typeMethods.r8vec_linspace_new(n, x_first, x_last);

            double[] u = FEM_1D_BVP.fem1d_bvp_linear(n, FEM_Test_Methods.a8, FEM_Test_Methods.c8, FEM_Test_Methods.f8, x);

            Console.WriteLine("");
            Console.WriteLine("     I    X         U         Uexact    Error");
            Console.WriteLine("");

            for (i = 0; i < n; i++)
            {
                double uexact = FEM_Test_Methods.exact8(x[i]);
                Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                       + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + uexact.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + Math.Abs(u[i] - uexact).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }

            double e1 = FEM_Error.l1_error(n, x, u, FEM_Test_Methods.exact8);
            double e2 = FEM_Error.l2_error_linear(n, x, u, FEM_Test_Methods.exact8);
            double h1s = FEM_Error.h1s_error_linear(n, x, u, FEM_Test_Methods.exact_ux8);
            double mx = FEM_Error.max_error_linear(n, x, u, FEM_Test_Methods.exact8);
            Console.WriteLine("");
            Console.WriteLine("  l1 norm of error  = " + e1 + "");
            Console.WriteLine("  L2 norm of error  = " + e2 + "");
            Console.WriteLine("  Seminorm of error = " + h1s + "");
            Console.WriteLine("  Max norm of error = " + mx + "");
        }

        static void test09()
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
            const int n = 11;

            Console.WriteLine("");
            Console.WriteLine("TEST09");
            Console.WriteLine("  Solve -( A(x) U'(x) )' + C(x) U(x) = F(x)");
            Console.WriteLine("  for 0 < x < 1, with U(0) = U(1) = 0.");
            Console.WriteLine("  A9(X)  = 1.0");
            Console.WriteLine("  C9(X)  = 0.0");
            Console.WriteLine("  F9(X) = X * ( X + 3 ) * exp ( X ),   X <= 2/3");
            Console.WriteLine("        = 2 * exp ( 2/3),                   2/3 < X");
            Console.WriteLine("  U9(X) = X * ( 1 - X ) * exp ( X ),   X <= 2/3");
            Console.WriteLine("        = X * ( 1 - X ) * exp ( 2/3 ),      2/3 < X");
            Console.WriteLine("");
            Console.WriteLine("  Number of nodes = " + n + "");
//
//  Geometry definitions.
//
            double x_first = 0.0;
            double x_last = 1.0;
            double[] x = typeMethods.r8vec_linspace_new(n, x_first, x_last);

            double[] u = FEM_1D_BVP.fem1d_bvp_linear(n, FEM_Test_Methods.a9, FEM_Test_Methods.c9, FEM_Test_Methods.f9, x);

            Console.WriteLine("");
            Console.WriteLine("     I    X         U         Uexact    Error");
            Console.WriteLine("");

            for (i = 0; i < n; i++)
            {
                double uexact = FEM_Test_Methods.exact9(x[i]);
                Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                       + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + uexact.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + Math.Abs(u[i] - uexact).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }

            double e1 = FEM_Error.l1_error(n, x, u, FEM_Test_Methods.exact9);
            double e2 = FEM_Error.l2_error_linear(n, x, u, FEM_Test_Methods.exact9);
            double h1s = FEM_Error.h1s_error_linear(n, x, u, FEM_Test_Methods.exact_ux9);
            double mx = FEM_Error.max_error_linear(n, x, u, FEM_Test_Methods.exact9);
            Console.WriteLine("");
            Console.WriteLine("  l1 norm of error  = " + e1 + "");
            Console.WriteLine("  L2 norm of error  = " + e2 + "");
            Console.WriteLine("  Seminorm of error = " + h1s + "");
            Console.WriteLine("  Max norm of error = " + mx + "");
        }

        static void test10()
//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests FEM1D_BVP_LINEAR.
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
            int e_log;
            const int e_log_max = 6;
            double h1;
            double l2;
            double mx;

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

            double[] h_plot = new double[e_log_max + 1];
            double[] h1_plot = new double[e_log_max + 1];
            double[] l2_plot = new double[e_log_max + 1];
            double[] mx_plot = new double[e_log_max + 1];
            int[] ne_plot = new int[e_log_max + 1];

            for (e_log = 0; e_log <= e_log_max; e_log++)
            {
                int ne = (int) Math.Pow(2, e_log);

                int n = ne + 1;
                const double x_lo = 0.0;
                const double x_hi = 1.0;
                double[] x = typeMethods.r8vec_linspace_new(n, x_lo, x_hi);

                double[] u = FEM_1D_BVP.fem1d_bvp_linear(n, FEM_Test_Methods.a10, FEM_Test_Methods.c10, FEM_Test_Methods.f10,
                    x);

                ne_plot[e_log] = ne;

                h_plot[e_log] = (x_hi - x_lo) / ne;

                l2 = FEM_Error.l2_error_linear(n, x, u, FEM_Test_Methods.exact10);
                l2_plot[e_log] = l2;

                h1 = FEM_Error.h1s_error_linear(n, x, u, FEM_Test_Methods.exact_ux10);
                h1_plot[e_log] = h1;

                mx = FEM_Error.max_error_linear(n, x, u, FEM_Test_Methods.exact10);
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
                int ne1 = ne_plot[e_log];
                int ne2 = ne_plot[e_log + 1];
                double r = ne2 / (double) ne1;
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
            string data_filename = "data.txt";
            List<string> data_unit = new();
            for (e_log = 0; e_log <= e_log_max; e_log++)
            {
                data_unit.Add("  " + ne_plot[e_log]
                                   + "  " + l2_plot[e_log]
                                   + "  " + h1_plot[e_log]
                                   + "  " + mx_plot[e_log] + "");
            }

            File.WriteAllLines(data_filename, data_unit);

            Console.WriteLine("");
            Console.WriteLine("  Created graphics data file \"" + data_filename + "\".");
//
//  Plot the L2 error as a function of NE.
//
            string command_filename = "commands_l2.txt";

            List<string> command_unit = new();

            string output_filename = "l2.png";

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

            Console.WriteLine("  Created graphics command file \"" + command_filename + "\".");
//
//  Plot the H1 error as a function of NE.
//
            command_filename = "commands_h1.txt";

            command_unit.Clear();

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

            Console.WriteLine("  Created graphics command file \"" + command_filename + "\".");
//
//  Plot the MX error as a function of NE.
//
            command_filename = "commands_mx.txt";

            command_unit.Clear();

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
            Console.WriteLine("  Created graphics command file \"" + command_filename + "\".");
        }

    }
}