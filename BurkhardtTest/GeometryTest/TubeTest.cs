using System;
using System.Linq;
using Burkardt.Tube;
using Burkardt.Types;

namespace GeometryTest;

public static class TubeTest
{
    public static void tube_2d_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TUBE_2D_TEST tests TUBE_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int DIM_NUM = 2;
        const int TEST_NUM = 4;

        double[] dist_test = { 0.5, 0.5, 1.0, 1.0 };
        int[] n_test = { 4, 5, 5, 5 };
        double[] p_test = {
            0.0,  0.0,
            4.0,  3.0,
            4.0,  0.0,
            0.0,  0.0,
            0.0,  0.0,
            2.0,  0.0,
            2.0,  1.0,
            0.0,  1.0,
            0.0,  0.0,
            10.0, 20.0,
            20.0, 20.0,
            10.0, 10.0,
            20.0, 10.0,
            10.0, 20.0,
            0.0,  0.0,
            10.0,  0.0,
            10.0, 10.0,
            10.0,  0.0,
            0.0,  0.0 };
        int test;

        Console.WriteLine("");
        Console.WriteLine("TUBE_2D_TEST");
        Console.WriteLine("  TUBE_2D computes corners of a tube of radius");
        Console.WriteLine("  DIST surrounding a sequence of points.");

        for ( test = 0; test < TEST_NUM; test++ )
        {
            int n = n_test[test];
            double dist = dist_test[test];

            int nlo = 0;
            int j;
            for ( j = 0; j < test; j++ )
            {
                nlo += n_test[j];
            }

            double[] p = p_test.Skip( + DIM_NUM * nlo).ToArray();

            Console.WriteLine("");
            Console.WriteLine("  Test " + test + "");
            Console.WriteLine("  Number of points N = " + n + "");
            Console.WriteLine("  Tube radius DIST = " + dist + "");

            typeMethods.r8mat_transpose_print ( DIM_NUM, n, p, "  Points to surround:" );

            double[] p1 = new double[DIM_NUM*n];
            double[] p2 = new double[DIM_NUM*n];

            Geometry.tube_2d ( dist, n, p, ref p1, ref p2 );

            typeMethods.r8mat_transpose_print ( DIM_NUM, n, p1, "  P1:" );
            typeMethods.r8mat_transpose_print ( DIM_NUM, n, p2, "  P2:" );
        }
    }


}