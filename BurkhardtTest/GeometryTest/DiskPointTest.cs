using System;
using System.Linq;
using Burkardt.Types;

namespace GeometryTest;

public static class DiskPointTest
{
    public static void disk_point_dist_3d_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DISK_POINT_DIST_3D_TEST tests DISK_POINT_DIST_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;
        int TEST_NUM = 5;

        double[] axis = { 0.0, 1.0, 1.0 };
        double dist;
        double[] dist_test = { 2.0, 0.0, 0.0, 8.0, 10.0 };
        double[] p;
        double[] pc = { 0.0, 1.4142135, 1.4142135 };
        double[] p_test = {
            0.0,  0.0,        0.0,
            0.0,  0.70710677, 2.1213202,
            2.0,  1.4142135,  1.4142135,
            10.0,  1.4142135,  1.4142135,
            10.0,  5.6568542,  5.6568542 };
        double r = 2.0;
        int test;

        Console.WriteLine("");
        Console.WriteLine("DISK_POINT_DIST_3D_TEST");
        Console.WriteLine("  DISK_POINT_DIST_3D finds the distance from");
        Console.WriteLine("  a disk to a point in 3D.");

        Console.WriteLine("");
        Console.WriteLine("  Disk radius = " + r + "");
        typeMethods.r8vec_print ( DIM_NUM, pc, "  Disk center: " );
        typeMethods.r8vec_print ( DIM_NUM, axis, "  Disk axis: " );

        for ( test = 0; test < TEST_NUM; test++ )
        {
            p = p_test.Skip( + test * DIM_NUM).ToArray();

            typeMethods.r8vec_print ( DIM_NUM, p, "  Point: " );

            dist = Burkardt.Disk.Geometry.disk_point_dist_3d ( pc, r, axis, p );

            Console.WriteLine("");
            Console.WriteLine("  Distance = " + dist
                                              + "  Expected = " + dist_test[test] + "");
        }
    }

}