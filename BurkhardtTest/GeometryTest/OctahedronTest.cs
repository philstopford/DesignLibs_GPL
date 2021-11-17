using System;
using Burkardt.Octahedron;

namespace GeometryTest;

public static class OctahedronTest
{
    public static void test0475 ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0475 tests OCTAHEDRON_SIZE_3D, OCTAHEDRON_SHAPE_3D, SHAPE_PRINT_3D.
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
        int DIM_NUM = 3;

        int edge_num = 0;
        int face_num = 0;
        int[] face_order;
        int face_order_max = 0;
        int[] face_point;
        int point_num = 0;
        double[] point_coord;

        Console.WriteLine("");
        Console.WriteLine("TEST0475");
        Console.WriteLine("  For the octahedron,");
        Console.WriteLine("  OCTAHEDRON_SIZE_3D returns dimension information;");
        Console.WriteLine("  OCTAHEDRON_SHAPE_3D returns face and order information.");
        Console.WriteLine("  SHAPE_PRINT_3D prints this information.");
        //
        //  Get the data sizes.
        //
        Geometry.octahedron_size_3d ( ref point_num, ref edge_num, ref face_num, ref face_order_max );

        Console.WriteLine("");
        Console.WriteLine("    Number of vertices: " + point_num + "");
        Console.WriteLine("    Number of edges   : " + edge_num + "");
        Console.WriteLine("    Number of faces   : " + face_num + "");
        Console.WriteLine("    Maximum face order: " + face_order_max + "");
        //
        //  Make room for the data.
        //
        face_order = new int[face_num];
        face_point = new int[face_order_max*face_num];
        point_coord = new double[DIM_NUM*point_num];
        //
        //  Get the data.
        //
        Geometry.octahedron_shape_3d ( point_num, face_num, face_order_max, ref point_coord,
            ref face_order, ref face_point );
        //
        //  Print the data.
        //
        Burkardt.Geometry.Shape.shape_print_3d ( point_num, face_num, face_order_max,
            point_coord, face_order, face_point );

    }

}