using System;
using Burkardt.FEM;
using Burkardt.Types;

namespace FEM1DBVPQuadraticTest
{
    class Program
    {
        static void Main(string[] args)
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

            u = FEM_1D_BVP_Quadratic.fem1d_bvp_quadratic(n, FEM_Test_Methods.a00, FEM_Test_Methods.c00, FEM_Test_Methods.f00, x);

            Console.WriteLine("");
            Console.WriteLine("     I    X         U         Uexact    Error");
            Console.WriteLine("");

            for (i = 0; i < n; i++)
            {
                uexact = FEM_Test_Methods.exact00(x[i]);
                Console.WriteLine("  " + i.ToString().PadLeft(4)
                    + "  " + x[i].ToString().PadLeft(8)
                    + "  " + u[i].ToString().PadLeft(14)
                    + "  " + uexact.ToString().PadLeft(14)
                    + "  " + Math.Abs(u[i] - uexact).ToString().PadLeft(14) + "");
            }

            e1 = FEM_Error.l1_error(n, x, u, FEM_Test_Methods.exact00);
            e2 = FEM_Error.l2_error_quadratic(n, x, u,FEM_Test_Methods. exact00);
            h1s = FEM_Error.h1s_error_quadratic(n, x, u, FEM_Test_Methods.exact_ux00);
            mx = FEM_Error.max_error_quadratic(n, x, u, FEM_Test_Methods.exact00);
            Console.WriteLine("");
            Console.WriteLine("  l1 norm of error  = " + e1 + "");
            Console.WriteLine("  L2 norm of error  = " + e2 + "");
            Console.WriteLine("  Seminorm of error = " + h1s + "");
            Console.WriteLine("  Max norm of error = " + mx + "");

        }
    }
}