using System;
using System.Globalization;
using Burkardt.FEM;
using Burkardt.Types;

namespace FEM2DBVPQuadraticTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for FEM2D_BVP_QUADRATIC_TEST.
        //
        //  Discussion:
        //
        //    FEM2D_BVP_QUADRATIC_TEST tests the FEM2D_BVP_QUADRATIC library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("FEM2D_BVP_QUADRATIC_TEST");
        Console.WriteLine("  Test the FEM2D_BVP_QUADRATIC library.");

        test01();

        Console.WriteLine("");
        Console.WriteLine("FEM2D_BVP_QUADRATIC_TEST");
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
        int j;
        const int nx = 3;
        const int ny = 3;

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
        const double x_first = 0.0;
        const double x_last = 1.0;
        double[] x = typeMethods.r8vec_even_new(nx, x_first, x_last);

        const double y_first = 0.0;
        const double y_last = 1.0;
        double[] y = typeMethods.r8vec_even_new(ny, y_first, y_last);

        double[] u = FEM_2D_BVP_Quadratic.fem2d_bvp_quadratic(nx, ny, a1, c1, f1, x, y);

        Console.WriteLine("");
        Console.WriteLine("     I     J    X         Y         U         Uexact    Error");
        Console.WriteLine("");

        for (j = 0; j < ny; j++)
        {
            int i;
            for (i = 0; i < nx; i++)
            {
                double uexact = exact1(x[i], y[j]);
                Console.WriteLine(i.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                                                                                      + j.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                                                                                      + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                                                                                      + y[j].ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                                                                                      + u[i + j * nx].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "  "
                                                                                      + uexact.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "  "
                                                                                      + Math.Abs(u[i + j * nx] - uexact).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }
        }

        double e1 = FEM_2D_BVP_Quadratic.fem2d_l1_error(nx, ny, x, y, u, exact1);
        double e2 = FEM_2D_BVP_Quadratic.fem2d_l2_error_quadratic(nx, ny, x, y, u, exact1);
        double h1s = FEM_2D_BVP_Quadratic.fem2d_h1s_error_quadratic(nx, ny, x, y, u, exact_ux1, exact_uy1);
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
        const double value = 1.0;

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