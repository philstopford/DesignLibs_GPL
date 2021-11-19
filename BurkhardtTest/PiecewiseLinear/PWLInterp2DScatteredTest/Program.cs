using System;
using Burkardt.PiecewiseLinear;
using Burkardt.TriangulationNS;
using InterpTest;

namespace PWLInterp2DScatteredTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for PWL_INTERP_2D_SCATTERED_TEST.
        //
        //  Discussion:
        //
        //    PWL_INTERP_2D_SCATTERED_TEST tests the PWL_INTERP_2D_SCATTERED library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int prob;
        int prob_num;

        Console.WriteLine("");
        Console.WriteLine("PWL_INTERP_2D_SCATTERED_TEST:");
        Console.WriteLine("  Test the PWL_INTERP_2D_SCATTERED library.");
        Console.WriteLine("  The R8LIB library is needed.");
        Console.WriteLine("  This test also needs the TEST_INTERP_2D library.");

        test01();
        test02();
        //
        //  Numerical tests.
        //
        prob_num = Data_2D.f00_num();

        for (prob = 1; prob <= prob_num; prob++)
        {
            test03(prob);
        }

        Console.WriteLine("");
        Console.WriteLine("PWL_INTERP_2D_SCATTERED_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests R8TRIS2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] element_neighbor = new int[3 * 2 * 9];
        int element_num = 0;
        int node_num = 9;
        double[] node_xy =
        {
            0.0, 0.0,
            0.0, 1.0,
            0.2, 0.5,
            0.3, 0.6,
            0.4, 0.5,
            0.6, 0.4,
            0.6, 0.5,
            1.0, 0.0,
            1.0, 1.0
        };
        int[] triangle = new int[3 * 2 * 9];

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  R8TRIS2 computes the Delaunay triangulation of");
        Console.WriteLine("  a set of nodes in 2D.");
        //
        //  Set up the Delaunay triangulation.
        //
        Delauney.r8tris2(node_num, ref node_xy, ref element_num, ref triangle, ref element_neighbor);

        Print.triangulation_order3_print(node_num, element_num, node_xy,
            triangle, element_neighbor);
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests PWL_INTERP_2D_SCATTERED_VALUE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] element_neighbor = new int[3 * 2 * 9];
        int element_num = 0;
        int i;
        int j;
        int k;
        int ni = 25;
        int node_num = 9;
        double[] node_xy =
        {
            0.0, 0.0,
            0.0, 1.0,
            0.2, 0.5,
            0.3, 0.6,
            0.4, 0.5,
            0.6, 0.4,
            0.6, 0.5,
            1.0, 0.0,
            1.0, 1.0
        };
        int[] triangle = new int[3 * 2 * 9];
        double x;
        double[] xyi = new double[2 * 25];
        double y;
        double[] zd = new double[9];
        double ze;
        double[] zi;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  PWL_INTERP_2D_SCATTERED_VALUE evaluates a");
        Console.WriteLine("  piecewise linear interpolant to scattered data.");
        //
        //  Set up the Delaunay triangulation.
        //
        Delauney.r8tris2(node_num, ref node_xy, ref element_num, ref triangle, ref element_neighbor);

        for (j = 0; j < element_num; j++)
        {
            for (i = 0; i < 3; i++)
            {
                switch (element_neighbor[i + j * 3])
                {
                    case > 0:
                        element_neighbor[i + j * 3] -= 1;
                        break;
                }
            }
        }

        Print.triangulation_order3_print(node_num, element_num, node_xy,
            triangle, element_neighbor);
        //
        //  Define the Z data.
        //
        for (i = 0; i < node_num; i++)
        {
            x = node_xy[0 + i * 2];
            y = node_xy[1 + i * 2];
            zd[i] = x + 2.0 * y;
        }

        //
        //  Define the interpolation points.
        //
        k = 0;
        for (i = 0; i <= 4; i++)
        {
            for (j = 0; j <= 4; j++)
            {
                xyi[0 + k * 2] = (i - 1) / 4.0;
                xyi[1 + k * 2] = (j - 1) / 4.0;
                k += 1;
            }
        }

        //
        //  Evaluate the interpolant.
        //
        zi = Interp2D.pwl_interp_2d_scattered_value(node_num, node_xy, zd, element_num,
            triangle, element_neighbor, ni, xyi);

        Console.WriteLine("");
        Console.WriteLine("     K      Xi(K)       Yi(K)       Zi(K)       Z(X,Y)");
        Console.WriteLine("");
        for (k = 0; k < ni; k++)
        {
            ze = xyi[0 + k * 2] + 2.0 * xyi[1 + k * 2];
            Console.WriteLine("  " + k.ToString().PadLeft(4)
                                   + "  " + xyi[0 + k * 2].ToString().PadLeft(10)
                                   + "  " + xyi[1 + k * 2].ToString().PadLeft(10)
                                   + "  " + zi[k].ToString().PadLeft(10)
                                   + "  " + ze.ToString().PadLeft(10) + "");
        }
    }

    private static void test03(int prob)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests PWL_INTERP_2D_SCATTERED_VALUE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] element_neighbor;
        int element_num = 0;
        int g;
        int i;
        int j;
        int k;
        int nd;
        int ni = 25;
        double rms;
        int[] triangle;
        double[] xd;
        double[] xi = new double[25];
        double[] xyd;
        double[] xyi = new double[2 * 25];
        double[] yd;
        double[] yi = new double[25];
        double[] zd;
        double[] ze;
        double[] zi;

        g = 2;
        nd = Data_2D.g00_size(g);

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  PWL_INTERP_2D_SCATTERED_VALUE evaluates a");
        Console.WriteLine("  piecewise linear interpolant to scattered data.");
        Console.WriteLine("  Here, we use grid number " + g + "");
        Console.WriteLine("  with " + nd + " scattered points in the unit square");
        Console.WriteLine("  on problem " + prob + "");
        //
        //  Get the data points and evaluate the function there.
        //
        xd = new double[nd];
        yd = new double[nd];

        Data_2D.g00_xy(g, nd, ref xd, ref yd);

        zd = new double[nd];
        Data_2D.f00_f0(prob, nd, xd, yd, ref zd);

        xyd = new double[2 * nd];

        for (i = 0; i < nd; i++)
        {
            xyd[0 + i * 2] = xd[i];
            xyd[1 + i * 2] = yd[i];
        }

        //
        //  Set up the Delaunay triangulation.
        //
        element_neighbor = new int[3 * 2 * nd];
        triangle = new int[3 * 2 * nd];

        Delauney.r8tris2(nd, ref xyd, ref element_num, ref triangle, ref element_neighbor);

        for (j = 0; j < element_num; j++)
        {
            for (i = 0; i < 3; i++)
            {
                switch (element_neighbor[i + j * 3])
                {
                    case > 0:
                        element_neighbor[i + j * 3] -= 1;
                        break;
                }
            }
        }

        //
        //  Define the interpolation points.
        //
        k = 0;
        for (i = 1; i <= 5; i++)
        {
            for (j = 1; j <= 5; j++)
            {
                xyi[0 + k * 2] = (2 * i - 1) / 10.0;
                xyi[1 + k * 2] = (2 * j - 1) / 10.0;
                k += 1;
            }
        }

        for (k = 0; k < ni; k++)
        {
            xi[k] = xyi[0 + k * 2];
            yi[k] = xyi[1 + k * 2];
        }

        ze = new double[ni];
        Data_2D.f00_f0(prob, ni, xi, yi, ref ze);
        //
        //  Evaluate the interpolant.
        //
        zi = Interp2D.pwl_interp_2d_scattered_value(nd, xyd, zd, element_num,
            triangle, element_neighbor, ni, xyi);

        rms = 0.0;
        for (k = 0; k < ni; k++)
        {
            rms += Math.Pow(zi[k] - ze[k], 2);
        }

        rms = Math.Sqrt(rms / ni);

        Console.WriteLine("");
        Console.WriteLine("  RMS error is " + rms + "");

        Console.WriteLine("");
        Console.WriteLine("     K      Xi(K)       Yi(K)       Zi(K)       Z(X,Y)");
        Console.WriteLine("");

        for (k = 0; k < ni; k++)
        {
            Console.WriteLine("  " + k.ToString().PadLeft(4)
                                   + "  " + xyi[0 + k * 2].ToString().PadLeft(10)
                                   + "  " + xyi[1 + k * 2].ToString().PadLeft(10)
                                   + "  " + zi[k].ToString().PadLeft(10)
                                   + "  " + ze[k].ToString().PadLeft(10) + "");
        }
    }
}