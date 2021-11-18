using System;

namespace GeometryTest;

public static class DualTest
{
    public static void dual_size_3d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DUAL_SIZE_3D_TEST tests DUAL_SIZE_3D;
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

        int edge_num1 = 0;
        int edge_num2 = 0;
        int[] edge_point1;
        int face_num1 = 0;
        int face_num2 = 0;
        int[] face_order1;
        int face_order_max1 = 0;
        int face_order_max2 = 0;
        int[] face_point1;
        int point_num1 = 0;
        int point_num2 = 0;
        double[] point_coord1;

        Console.WriteLine("");
        Console.WriteLine("DUAL_SIZE_3D_TEST");
        Console.WriteLine("  DUAL_SIZE_3D finds the sizes of the dual of a");
        Console.WriteLine("  polyhedron;");
        //
        //  Get the CUBE shape.
        //
        Burkardt.Cube.Geometry.cube_size_3d(ref point_num1, ref edge_num1, ref face_num1, ref face_order_max1);

        Console.WriteLine("");
        Console.WriteLine("  The cube:");
        Console.WriteLine("    Number of vertices: " + point_num1 + "");
        Console.WriteLine("    Number of edges   : " + edge_num1 + "");
        Console.WriteLine("    Number of faces   : " + face_num1 + "");
        Console.WriteLine("    Maximum face order: " + face_order_max1 + "");

        face_order1 = new int[face_num1];
        face_point1 = new int[face_order_max1 * face_num1];
        point_coord1 = new double[DIM_NUM * point_num1];

        Burkardt.Cube.Geometry.cube_shape_3d(point_num1, face_num1, face_order_max1, ref point_coord1,
            ref face_order1, ref face_point1);

        Burkardt.Geometry.Dual.dual_size_3d(point_num1, edge_num1, face_num1, face_order_max1,
            point_coord1, face_order1, face_point1, ref point_num2, ref edge_num2,
            ref face_num2, ref face_order_max2);

        Console.WriteLine("");
        Console.WriteLine("  The dual of the cube:");
        Console.WriteLine("    Number of vertices: " + point_num2 + "");
        Console.WriteLine("    Number of edges   : " + edge_num2 + "");
        Console.WriteLine("    Number of faces   : " + face_num2 + "");
        Console.WriteLine("    Maximum face order: " + face_order_max2 + "");

        //
        //  Get the DODECAHEDRON shape.
        //
        Burkardt.Dodecahedron.Geometry.dodec_size_3d(ref point_num1, ref edge_num1, ref face_num1,
            ref face_order_max1);

        Console.WriteLine("");
        Console.WriteLine("  The dodecahedron:");
        Console.WriteLine("    Number of vertices: " + point_num1 + "");
        Console.WriteLine("    Number of edges   : " + edge_num1 + "");
        Console.WriteLine("    Number of faces   : " + face_num1 + "");
        Console.WriteLine("    Maximum face order: " + face_order_max1 + "");

        face_order1 = new int[face_num1];
        face_point1 = new int[face_order_max1 * face_num1];
        point_coord1 = new double[DIM_NUM * point_num1];

        Burkardt.Dodecahedron.Geometry.dodec_shape_3d(point_num1, face_num1, face_order_max1, ref point_coord1,
            ref face_order1, ref face_point1);

        Burkardt.Geometry.Dual.dual_size_3d(point_num1, edge_num1, face_num1, face_order_max1,
            point_coord1, face_order1, face_point1, ref point_num2, ref edge_num2,
            ref face_num2, ref face_order_max2);

        Console.WriteLine("");
        Console.WriteLine("  The dual of the dodecahedron:");
        Console.WriteLine("    Number of vertices: " + point_num2 + "");
        Console.WriteLine("    Number of edges   : " + edge_num2 + "");
        Console.WriteLine("    Number of faces   : " + face_num2 + "");
        Console.WriteLine("    Maximum face order: " + face_order_max2 + "");

        //
        //  Get the ICOSAHEDRON shape.
        //
        Burkardt.Icosahedron.Geometry.icos_size(ref point_num1, ref edge_num1, ref face_num1, ref face_order_max1);

        Console.WriteLine("");
        Console.WriteLine("  The icosahedron:");
        Console.WriteLine("    Number of vertices: " + point_num1 + "");
        Console.WriteLine("    Number of edges   : " + edge_num1 + "");
        Console.WriteLine("    Number of faces   : " + face_num1 + "");
        Console.WriteLine("    Maximum face order: " + face_order_max1 + "");

        edge_point1 = new int[2 * edge_num1];
        face_order1 = new int[face_num1];
        face_point1 = new int[face_order_max1 * face_num1];
        point_coord1 = new double[DIM_NUM * point_num1];

        Burkardt.Icosahedron.Geometry.icos_shape(point_num1, edge_num1, face_num1, face_order_max1,
            ref point_coord1, ref edge_point1, ref face_order1, ref face_point1);

        Burkardt.Geometry.Dual.dual_size_3d(point_num1, edge_num1, face_num1, face_order_max1,
            point_coord1, face_order1, face_point1, ref point_num2, ref edge_num2,
            ref face_num2, ref face_order_max2);

        Console.WriteLine("");
        Console.WriteLine("  The dual of the icosahedron:");
        Console.WriteLine("    Number of vertices: " + point_num2 + "");
        Console.WriteLine("    Number of edges   : " + edge_num2 + "");
        Console.WriteLine("    Number of faces   : " + face_num2 + "");
        Console.WriteLine("    Maximum face order: " + face_order_max2 + "");

        //
        //  Get the octahedron shape.
        //
        Burkardt.Octahedron.Geometry.octahedron_size_3d(ref point_num1, ref edge_num1, ref face_num1,
            ref face_order_max1);

        Console.WriteLine("");
        Console.WriteLine("  The octahedron:");
        Console.WriteLine("    Number of vertices: " + point_num1 + "");
        Console.WriteLine("    Number of edges   : " + edge_num1 + "");
        Console.WriteLine("    Number of faces   : " + face_num1 + "");
        Console.WriteLine("    Maximum face order: " + face_order_max1 + "");

        face_order1 = new int[face_num1];
        face_point1 = new int[face_order_max1 * face_num1];
        point_coord1 = new double[DIM_NUM * point_num1];

        Burkardt.Octahedron.Geometry.octahedron_shape_3d(point_num1, face_num1, face_order_max1,
            ref point_coord1, ref face_order1, ref face_point1);

        Burkardt.Geometry.Dual.dual_size_3d(point_num1, edge_num1, face_num1, face_order_max1,
            point_coord1, face_order1, face_point1, ref point_num2, ref edge_num2,
            ref face_num2, ref face_order_max2);

        Console.WriteLine("");
        Console.WriteLine("  The dual of the octahedron:");
        Console.WriteLine("    Number of vertices: " + point_num2 + "");
        Console.WriteLine("    Number of edges   : " + edge_num2 + "");
        Console.WriteLine("    Number of faces   : " + face_num2 + "");
        Console.WriteLine("    Maximum face order: " + face_order_max2 + "");

        //
        //  Get the soccer ball shape.
        //
        Burkardt.Icosahedron.Geometry.soccer_size_3d(ref point_num1, ref edge_num1, ref face_num1,
            ref face_order_max1);

        Console.WriteLine("");
        Console.WriteLine("  The soccer ball:");
        Console.WriteLine("    Number of vertices: " + point_num1 + "");
        Console.WriteLine("    Number of edges   : " + edge_num1 + "");
        Console.WriteLine("    Number of faces   : " + face_num1 + "");
        Console.WriteLine("    Maximum face order: " + face_order_max1 + "");

        face_order1 = new int[face_num1];
        face_point1 = new int[face_order_max1 * face_num1];
        point_coord1 = new double[DIM_NUM * point_num1];

        Burkardt.Icosahedron.Geometry.soccer_shape_3d(point_num1, face_num1, face_order_max1,
            ref point_coord1, ref face_order1, ref face_point1);

        Burkardt.Geometry.Dual.dual_size_3d(point_num1, edge_num1, face_num1, face_order_max1,
            point_coord1, face_order1, face_point1, ref point_num2, ref edge_num2,
            ref face_num2, ref face_order_max2);

        Console.WriteLine("");
        Console.WriteLine("  The dual of the soccer ball:");
        Console.WriteLine("    Number of vertices: " + point_num2 + "");
        Console.WriteLine("    Number of edges   : " + edge_num2 + "");
        Console.WriteLine("    Number of faces   : " + face_num2 + "");
        Console.WriteLine("    Maximum face order: " + face_order_max2 + "");

        //
        //  Get the tetrahedron shape.
        //
        Burkardt.TetrahedronNS.Geometry.tetrahedron_size_3d(ref point_num1, ref edge_num1, ref face_num1,
            ref face_order_max1);

        Console.WriteLine("");
        Console.WriteLine("  The tetrahedron:");
        Console.WriteLine("    Number of vertices: " + point_num1 + "");
        Console.WriteLine("    Number of edges   : " + edge_num1 + "");
        Console.WriteLine("    Number of faces   : " + face_num1 + "");
        Console.WriteLine("    Maximum face order: " + face_order_max1 + "");

        face_order1 = new int[face_num1];
        face_point1 = new int[face_order_max1 * face_num1];
        point_coord1 = new double[DIM_NUM * point_num1];

        Burkardt.TetrahedronNS.Geometry.tetrahedron_shape_3d(point_num1, face_num1, face_order_max1,
            ref point_coord1, ref face_order1, ref face_point1);

        Burkardt.Geometry.Dual.dual_size_3d(point_num1, edge_num1, face_num1, face_order_max1,
            point_coord1, face_order1, face_point1, ref point_num2, ref edge_num2,
            ref face_num2, ref face_order_max2);

        Console.WriteLine("");
        Console.WriteLine("  The dual of the tetrahedron:");
        Console.WriteLine("    Number of vertices: " + point_num2 + "");
        Console.WriteLine("    Number of edges   : " + edge_num2 + "");
        Console.WriteLine("    Number of faces   : " + face_num2 + "");
        Console.WriteLine("    Maximum face order: " + face_order_max2 + "");


    }

    public static void dual_shape_3d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DUAL_SHAPE_3D_TEST tests DUAL_SHAPE_3D;
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

        int edge_num1 = 0;
        int edge_num2 = 0;
        int face_num1 = 0;
        int face_num2 = 0;
        int[] face_order1;
        int[] face_order2;
        int face_order_max1 = 0;
        int face_order_max2 = 0;
        int[] face_point1;
        int[] face_point2;
        int point_num1 = 0;
        int point_num2 = 0;
        double[] point_coord1;
        double[] point_coord2;

        Console.WriteLine("");
        Console.WriteLine("DUAL_SHAPE_3D_TEST");
        Console.WriteLine("  DUAL_SHAPE_3D finds the dual of a polyhedron.");
        //
        //  Get the dodecahedron shape.
        //
        Console.WriteLine("");
        Console.WriteLine("  The dodecahedron:");

        Burkardt.Dodecahedron.Geometry.dodec_size_3d(ref point_num1, ref edge_num1, ref face_num1,
            ref face_order_max1);

        Console.WriteLine("");
        Console.WriteLine("    Number of vertices: " + point_num1 + "");
        Console.WriteLine("    Number of edges   : " + edge_num1 + "");
        Console.WriteLine("    Number of faces   : " + face_num1 + "");
        Console.WriteLine("    Maximum face order: " + face_order_max1 + "");

        face_order1 = new int[face_num1];
        face_point1 = new int [face_order_max1 * face_num1];
        point_coord1 = new double[DIM_NUM * point_num1];

        Burkardt.Dodecahedron.Geometry.dodec_shape_3d(point_num1, face_num1, face_order_max1, ref point_coord1,
            ref face_order1, ref face_point1);

        Burkardt.Geometry.Shape.shape_print_3d(point_num1, face_num1, face_order_max1,
            point_coord1, face_order1, face_point1);
        //
        //  Get the dual.
        //
        Console.WriteLine("");
        Console.WriteLine("  The dual of the dodecahedron:");

        Burkardt.Geometry.Dual.dual_size_3d(point_num1, edge_num1, face_num1, face_order_max1,
            point_coord1, face_order1, face_point1, ref point_num2, ref edge_num2,
            ref face_num2, ref face_order_max2);

        Console.WriteLine("");
        Console.WriteLine("    Number of vertices: " + point_num2 + "");
        Console.WriteLine("    Number of edges   : " + edge_num2 + "");
        Console.WriteLine("    Number of faces   : " + face_num2 + "");
        Console.WriteLine("    Maximum face order: " + face_order_max2 + "");

        face_order2 = new int[face_num2];
        face_point2 = new int[face_order_max2 * face_num2];
        point_coord2 = new double[DIM_NUM * point_num2];

        Burkardt.Geometry.Dual.dual_shape_3d(point_num1, face_num1, face_order_max1, point_coord1,
            face_order1, face_point1, point_num2, face_num2, face_order_max2,
            ref point_coord2, ref face_order2, ref face_point2);

        Burkardt.Geometry.Shape.shape_print_3d(point_num2, face_num2, face_order_max2,
            point_coord2, face_order2, face_point2);

    }

}