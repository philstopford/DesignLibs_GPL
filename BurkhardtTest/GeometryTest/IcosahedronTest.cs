using System;
using Burkardt.Icosahedron;

namespace GeometryTest;

public static class IcosahedronTest
{
    public static void test0325 ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0325 tests ICOS_SIZE, ICOS_SHAPE, SHAPE_PRINT_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;

        int edge_num = 0;
        int[] edge_point;
        int face_num = 0;
        int[] face_order;
        int face_order_max = 0;
        int[] face_point;
        int point_num = 0;
        double[] point_coord;

        Console.WriteLine("");
        Console.WriteLine("TEST0325");
        Console.WriteLine("  For the icosahedron,");
        Console.WriteLine("  ICOS_SIZE returns dimension information;");
        Console.WriteLine("  ICOS_SHAPE returns face and order information.");
        Console.WriteLine("  SHAPE_PRINT_3D prints this information.");
        //
        //  Get the data sizes.
        //
        Geometry.icos_size ( ref point_num, ref edge_num, ref face_num, ref face_order_max );

        Console.WriteLine("");
        Console.WriteLine("    Number of vertices: " + point_num + "");
        Console.WriteLine("    Number of edges   : " + edge_num + "");
        Console.WriteLine("    Number of faces   : " + face_num + "");
        Console.WriteLine("    Maximum face order: " + face_order_max + "");
        //
        //  Make room for the data.
        //
        point_coord = new double[DIM_NUM*point_num];
        edge_point = new int[2*edge_num];
        face_order = new int[face_num];
        face_point = new int[face_order_max*face_num];
        //
        //  Get the data.
        //
        Geometry.icos_shape ( point_num, edge_num, face_num, face_order_max, ref point_coord,
            ref edge_point, ref face_order, ref face_point );
        //
        //  Print the data.
        //
        Burkardt.Geometry.Shape.shape_print_3d ( point_num, face_num, face_order_max,
            point_coord, face_order, face_point );
    }

    public static void test179 ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST179 tests SOCCER_SIZE_3D, SOCCER_SHAPE_3D, SHAPE_PRINT_3D.
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
        Console.WriteLine("TEST179");
        Console.WriteLine("  For the truncated icosahedron, or soccer ball,");
        Console.WriteLine("  SOCCER_SIZE_3D returns dimension information;");
        Console.WriteLine("  SOCCER_SHAPE_3D returns face and order information.");
        Console.WriteLine("  SHAPE_PRINT_3D prints this information.");
        //
        //  Get the data sizes.
        //
        Geometry.soccer_size_3d ( ref point_num, ref edge_num, ref face_num, ref face_order_max );

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
        Geometry.soccer_shape_3d ( point_num, face_num, face_order_max, ref point_coord,
            ref face_order, ref face_point );
        //
        //  Print the data.
        //
        Burkardt.Geometry.Shape.shape_print_3d ( point_num, face_num, face_order_max,
            point_coord, face_order, face_point );
            
    }
}