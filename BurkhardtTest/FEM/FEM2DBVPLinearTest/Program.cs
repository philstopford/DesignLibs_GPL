using System;
using Burkardt.FEM;
using Burkardt.Types;

namespace FEM2DBVPLinearTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for FEM2D_BVP_LINEAR_TEST.
        //
        //  Discussion:
        //
        //    FEM2D_BVP_LINEAR_TEST tests the FEM2D_BVP_LINEAR library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("FEM2D_BVP_LINEAR_TEST");
        Console.WriteLine("  Test the FEM2D_BVP_LINEAR library.");

        test01();

        Console.WriteLine("");
        Console.WriteLine("FEM2D_BVP_LINEAR_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
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
        //    Use A1, C1, F1, EXACT1, EXACT_UX1, EXACT_UX2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double e1;
        double e2;
        double h1s;
        int i;
        int j;
        int nx = 3;
        int ny = 3;
        double[] u;
        double uexact;
        double[] x;
        double x_first;
        double x_last;
        double[] y;
        double y_first;
        double y_last;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Solve - del ( A del U ) + C U = F ");
        Console.WriteLine("  on the unit square with zero boundary conditions.");
        Console.WriteLine("  A1(X,Y) = 1.0");
        Console.WriteLine("  C1(X,Y) = 0.0");
        Console.WriteLine("  F1(X,Y) = 2*X*(1-X)+2*Y*(1-Y)");
        Console.WriteLine("  U1(X,Y) = X * ( 1 - X ) * Y * ( 1 - Y )");
        Console.WriteLine("");
        Console.WriteLine("  Number of X grid values NX = " + nx + "");
        Console.WriteLine("  Number of Y grid values NY = " + ny + "");
        //
        //  Geometry definitions.
        //
        x_first = 0.0;
        x_last = 1.0;
        x = typeMethods.r8vec_even_new(nx, x_first, x_last);

        y_first = 0.0;
        y_last = 1.0;
        y = typeMethods.r8vec_even_new(ny, y_first, y_last);

        u = FEM_2D_BVP_Linear.fem2d_bvp_linear(nx, ny, a1, c1, f1, x, y);

        Console.WriteLine("");
        Console.WriteLine("     I     J    X         Y         U         Uexact    Error");
        Console.WriteLine("");

        for (j = 0; j < ny; j++)
        {
            for (i = 0; i < nx; i++)
            {
                uexact = exact1(x[i], y[j]);
                Console.WriteLine(i.ToString().PadLeft(6) + "  "
                                                          + j.ToString().PadLeft(4) + "  "
                                                          + x[i].ToString().PadLeft(8) + "  "
                                                          + y[j].ToString().PadLeft(8) + "  "
                                                          + u[i + j * nx].ToString().PadLeft(14) + "  "
                                                          + uexact.ToString().PadLeft(14) + "  "
                                                          + Math.Abs(u[i + j * nx] - uexact).ToString().PadLeft(14) + "");
            }
        }

        e1 = FEM_2D_BVP_Linear.fem2d_l1_error(nx, ny, x, y, u, exact1);
        e2 = FEM_2D_BVP_Linear.fem2d_l2_error_linear(nx, ny, x, y, u, exact1);
        h1s = FEM_2D_BVP_Linear.fem2d_h1s_error_linear(nx, ny, x, y, u, exact_ux1, exact_uy1);
        Console.WriteLine("");
        Console.WriteLine("  l1 norm of error  = " + e1 + "");
        Console.WriteLine("  L2 norm of error  = " + e2 + "");
        Console.WriteLine("  Seminorm of error = " + h1s + "");

    }

    private static double a1(double x, double y)

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
        //    20 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, Y, the evaluation point.
        //
        //    Output, double A1, the value of A(X,Y).
        //
    {
        double value = 1.0;

        return value;
    }

    private static double c1(double x, double y)

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
        //    20 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, Y, the evaluation point.
        //
        //    Output, double C1, the value of C(X,Y).
        //
    {
        const double value = 0.0;

        return value;
    }

    private static double exact1(double x, double y)

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
        //    20 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, Y, the evaluation point.
        //
        //    Output, double EXACT1, the value of U(X,Y).
        //
    {
        double value = x * (1.0 - x) * y * (1.0 - y);

        return value;
    }

    private static double exact_ux1(double x, double y)

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
        //    20 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, Y, the evaluation point.
        //
        //    Output, double EXACT_UX1, the value of dUdX(X,Y).
        //
    {
        double value = (1.0 - 2.0 * x) * (y - y * y);

        return value;
    }

    private static double exact_uy1(double x, double y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT_UY1 evaluates the derivative of exact solution #1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, Y, the evaluation point.
        //
        //    Output, double EXACT_UY1, the value of dUdY(X,Y).
        //
    {
        double value = (x - x * x) * (1.0 - 2.0 * y);

        return value;
    }

    private static double f1(double x, double y)

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
        //    20 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, Y, the evaluation point.
        //
        //    Output, double F1, the value of F(X,Y).
        //
    {
        double value = 2.0 * x * (1.0 - x)
                       + 2.0 * y * (1.0 - y);

        return value;
    }
}