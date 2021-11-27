using System;
using System.Globalization;
using Burkardt.FEM;
using Burkardt.MatrixNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkadt.FEM2DBVPSereneTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for FEM2D_BVP_SERENE_TEST.
        //
        //  Discussion:
        //
        //    FEM2D_BVP_SERENE_TEST tests the FEM2D_BVP_SERENE library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("FEM2D_BVP_SERENE_TEST");

        Console.WriteLine("  Test the FEM2D_BVP_SERENE library.");

        test01();
        test02();
        test03();

        Console.WriteLine("");
        Console.WriteLine("FEM2D_BVP_SERENE_TEST");
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
        //    Use A1, C1, F1, EXACT1, EXACT_UX1, EXACT_UX1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int nx = 5;
        const int ny = 5;

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
        double x_first = 0.0;
        double x_last = 1.0;
        double[] x = typeMethods.r8vec_linspace_new(nx, x_first, x_last);

        const double y_first = 0.0;
        const double y_last = 1.0;
        double[] y = typeMethods.r8vec_linspace_new(ny, y_first, y_last);

        const bool show11 = false;
        double[] u = FEM_2D_BVP_Serene.fem2d_bvp_serene(nx, ny, a1, c1, f1, x, y, show11);

        switch (nx * ny)
        {
            case <= 25:
            {
                Console.WriteLine("");
                Console.WriteLine("     I     J    X         Y               U               Uexact     Error");
                Console.WriteLine("");

                int k = 0;

                int j;
                for (j = 0; j < ny; j++)
                {
                    int inc = (j % 2) switch
                    {
                        0 => 1,
                        _ => 2
                    };

                    int i;
                    for (i = 0; i < nx; i += inc)
                    {
                        double uexact = exact1(x[i], y[j]);
                        Console.WriteLine(i.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                                                                  + j.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                                                                  + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                                                                  + y[j].ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                                                                  + u[k].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "  "
                                                                  + uexact.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "  "
                                                                  + Math.Abs(u[k] - uexact).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                        k += 1;
                    }
                }

                break;
            }
        }

        double e1 = FEM_2D_BVP_Serene.fem2d_l1_error_serene(nx, ny, x, y, u, exact1);
        double e2 = FEM_2D_BVP_Serene.fem2d_l2_error_serene(nx, ny, x, y, u, exact1);
        double h1s = FEM_2D_BVP_Serene.fem2d_h1s_error_serene(nx, ny, x, y, u, exact_ux1, exact_uy1);

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
        //    07 July 2014
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
        //    07 July 2014
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
        //    07 July 2014
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
        //    07 July 2014
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
        //    07 July 2014
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
        //    07 July 2014
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

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 checks the basis functions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        double xq;
        double[] xx =  {
                2.0, 1.0, 0.0, 0.0, 0.0, 1.0, 2.0, 2.0
            }
            ;
        double yq;
        double[] yy =  {
                5.0, 5.0, 5.0, 4.0, 3.0, 3.0, 3.0, 4.0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Basis function checks.");
        //
        //  Check that V is identity matrix at nodes.
        //
        Console.WriteLine("");
        Console.WriteLine("  The matrix Aij = V(j)(X(i),Y(i)) should be the identity.");
        Console.WriteLine("");

        double xw = 0.0;
        double ys = 3.0;
        double xe = 2.0;
        double yn = 5.0;

        for (i = 0; i < 8; i++)
        {
            xq = xx[i];
            yq = yy[i];
            FEM_2D_BVP_Serene.basis_serene(xq, yq, xw, ys, xe, yn, xx, yy);
            int j;
            for (j = 0; j < 8; j++)
            {
            }

            Console.WriteLine("");
        }

        //
        //  Check that VX and VY sum to zero anywhere.
        //
        Console.WriteLine("");
        Console.WriteLine("  The vectors dVdX(1:8)(X,Y) and dVdY(1:8)(X,Y)");
        Console.WriteLine("  should both sum to zero for any (X,Y).");

        int seed = 123456789;
        xq = 2.0 * UniformRNG.r8_uniform_01(ref seed);
        yq = 3.0 + 2.0 * UniformRNG.r8_uniform_01(ref seed);
        xw = 0.0;
        ys = 3.0;
        xe = 2.0;
        yn = 5.0;

        double[] vx = FEM_2D_BVP_Serene.basis_dx_serene(xq, yq, xw, ys, xe, yn, xx, yy);
        double[] vy = FEM_2D_BVP_Serene.basis_dy_serene(xq, yq, xw, ys, xe, yn, xx, yy);

        Console.WriteLine("");
        Console.WriteLine("  Random evaluation point is (" + xq + "," + yq + ")");
        Console.WriteLine("");
        Console.WriteLine("              dVdX        dVdY");
        Console.WriteLine("");
        for (i = 0; i < 8; i++)
        {
            Console.WriteLine(i.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                                                      + vx[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                                                      + vy[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Sum:  "
                          + typeMethods.r8vec_sum(8, vx).ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                          + typeMethods.r8vec_sum(8, vy).ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
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
        //    Use A3, C3, F3, EXACT3, EXACT_UX3, EXACT_UX3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        const int nx = 5;
        const int ny = 5;

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  Solve - del ( A del U ) + C U = F ");
        Console.WriteLine("  on the unit square with zero boundary conditions.");
        Console.WriteLine("  A1(X,Y) = 0.0");
        Console.WriteLine("  C1(X,Y) = 1.0");
        Console.WriteLine("  F1(X,Y) = X * ( 1 - X ) * Y * ( 1 - Y )");
        Console.WriteLine("  U1(X,Y) = X * ( 1 - X ) * Y * ( 1 - Y )");
        Console.WriteLine("");
        Console.WriteLine("  This example is contrived so that the system matrix");
        Console.WriteLine("  is the WATHEN matrix.");
        Console.WriteLine("");
        Console.WriteLine("  Number of X grid values NX = " + nx + "");
        Console.WriteLine("  Number of Y grid values NY = " + ny + "");
        //
        //  Geometry definitions.
        //
        const double x_first = 0.0;
        const double x_last = 1.0;
        double[] x = typeMethods.r8vec_linspace_new(nx, x_first, x_last);

        const double y_first = 0.0;
        const double y_last = 1.0;
        double[] y = typeMethods.r8vec_linspace_new(ny, y_first, y_last);

        double[] u = FEM_2D_BVP_Serene.fem2d_bvp_serene(nx, ny, a3, c3, f3, x, y, true);

        switch (nx * ny)
        {
            case <= 25:
            {
                Console.WriteLine("");
                Console.WriteLine("     I     J    X         Y               U               Uexact     Error");
                Console.WriteLine("");

                int k = 0;

                for (j = 0; j < ny; j++)
                {
                    int inc = (j % 2) switch
                    {
                        0 => 1,
                        _ => 2
                    };

                    for (i = 0; i < nx; i += inc)
                    {
                        double uexact = exact3(x[i], y[j]);
                        Console.WriteLine(i.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                                                                  + j.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                                                                  + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                                                                  + y[j].ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                                                                  + u[k].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "  "
                                                                  + uexact.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "  "
                                                                  + Math.Abs(u[k] - uexact).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                        k += 1;
                    }
                }

                break;
            }
        }

        double e1 = FEM_2D_BVP_Serene.fem2d_l1_error_serene(nx, ny, x, y, u, exact3);
        double e2 = FEM_2D_BVP_Serene.fem2d_l2_error_serene(nx, ny, x, y, u, exact3);
        double h1s = FEM_2D_BVP_Serene.fem2d_h1s_error_serene(nx, ny, x, y, u, exact_ux3, exact_uy3);

        Console.WriteLine("");
        Console.WriteLine("  l1 norm of error  = " + e1 + "");
        Console.WriteLine("  L2 norm of error  = " + e2 + "");
        Console.WriteLine("  Seminorm of error = " + h1s + "");

        //
        //  Pull out the Wathen matrix from MATLAB.
        //  It will have been multiplied by a random scale factor.
        //  While my numbering scheme is
        //    3  2  1
        //    4     8
        //    5  6  7
        //  the numbering scheme used here is 
        //    1  2  3
        //    4     5
        //    6  7  8
        //    
        double[] amat = WathenMatrix.wathen(1, 1, 8);

        double scale = 0.5 * amat[0 + 2 * 8];
        for (j = 0; j < 8; j++)
        {
            for (i = 0; i < 8; i++)
            {
                amat[i + j * 8] /= scale;
            }
        }

        typeMethods.r8mat_print(8, 8, amat, "  Wathen matrix:");
    }

    private static double a3(double x, double y)

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
        //    07 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, Y, the evaluation point.
        //
        //    Output, double A3, the value of A(X,Y).
        //
    {
        double value = 0.0;

        return value;
    }

    private static double c3(double x, double y)

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
        //    07 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, Y, the evaluation point.
        //
        //    Output, double C3, the value of C(X,Y).
        //
    {
        double value = 0;

        value = 1.0;

        return value;
    }

    private static double exact3(double x, double y)

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
        //    07 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, Y, the evaluation point.
        //
        //    Output, double EXACT3, the value of U(X,Y).
        //
    {
        double value = x * (1.0 - x) * y * (1.0 - y);

        return value;
    }

    private static double exact_ux3(double x, double y)

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
        //    07 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, Y, the evaluation point.
        //
        //    Output, double EXACT_UX3, the value of dUdX(X,Y).
        //
    {
        double value = (1.0 - 2.0 * x) * (y - y * y);

        return value;
    }

    private static double exact_uy3(double x, double y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT_UY3 evaluates the derivative of exact solution #3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, Y, the evaluation point.
        //
        //    Output, double EXACT_UY3, the value of dUdY(X,Y).
        //
    {
        double value = (x - x * x) * (1.0 - 2.0 * y);

        return value;
    }

    private static double f3(double x, double y)

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
        //    07 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, Y, the evaluation point.
        //
        //    Output, double F3, the value of F(X,Y).
        //
    {
        double value = x * (1.0 - x) * y * (1.0 - y);

        return value;
    }
}