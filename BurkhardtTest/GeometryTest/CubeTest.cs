using System;
using Burkardt.Cube;
using Burkardt.Geometry;

namespace GeometryTest;

public static class CubeTest
{
    public static void test020()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST020 tests CUBE_SIZE_3D and CUBE_SHAPE_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int DIM_NUM = 3;

        int edge_num = 0;
        int face_num = 0;
        int face_order_max = 0;
        int point_num = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST020");
        Console.WriteLine("  For the cube,");
        Console.WriteLine("  CUBE_SIZE_3D returns dimension information;");
        Console.WriteLine("  CUBE_SHAPE_3D returns face and order information.");
        Console.WriteLine("  SHAPE_PRINT_3D prints this information.");
        //
        //  Get the sizes.
        //
        Geometry.cube_size_3d(ref point_num, ref edge_num, ref face_num, ref face_order_max);

        Console.WriteLine("");
        Console.WriteLine("    Number of vertices: " + point_num + "");
        Console.WriteLine("    Number of edges   : " + edge_num + "");
        Console.WriteLine("    Number of faces   : " + face_num + "");
        Console.WriteLine("    Maximum face order: " + face_order_max + "");
        //
        //  Make room for the data.
        //
        int[] face_order = new int[face_num];
        int[] face_point = new int[face_order_max * face_num];
        double[] point_coord = new double[DIM_NUM * point_num];
        //
        //  Get the data.
        //
        Geometry.cube_shape_3d(point_num, face_num, face_order_max, ref point_coord,
            ref face_order, ref face_point);
        //
        //  Print the data.
        //
        Shape.shape_print_3d(point_num, face_num, face_order_max,
            point_coord, face_order, face_point);

    }

    public static void cube01_volume_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CUBE01_VOLUME_TEST tests CUBE01_VOLUME.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 January 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("CUBE01_VOLUME_TEST");
        Console.WriteLine("  CUBE01_VOLUME returns the volume of the unit cube.");

        double volume = Geometry.cube01_volume();

        Console.WriteLine("");
        Console.WriteLine("  Volume = " + volume + "");
    }

}