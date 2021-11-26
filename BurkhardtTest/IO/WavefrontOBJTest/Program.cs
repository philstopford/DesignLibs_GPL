using System;
using Burkardt.IO;

namespace WavefrontOBJTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for OBJ_IO_TEST.
        //
        //  Discussion:
        //
        //    OBJ_IO_TEST tests the OBJ_IO library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 January 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("OBJ_IO_TEST:");
        Console.WriteLine("  Test the OBJ_IO library.");

        string filename = "cube.obj";
        test01(filename);

        filename = "cube.obj";
        test02 ( filename );

        filename = "cube_normals.obj";
        test03(filename);

        filename = "cube_no_normals.obj";
        //test04 ( filename );

        Console.WriteLine("");
        Console.WriteLine("OBJ_IO_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01(string filename)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests OBJ_SIZE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 March 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int face_num = 0;
        int node_num = 0;
        int normal_num = 0;
        int order_max = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  OBJ_SIZE determines the size of various objects");
        Console.WriteLine("  in an OBJ file.");

        WavefrontOBJ.obj_size(filename, ref node_num, ref face_num, ref normal_num, ref order_max);

        WavefrontOBJ.obj_size_print(filename, node_num, face_num, normal_num, order_max);
    }

    private static void test02(string input_file_name)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests OBJ_READ.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 March 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int face_num = 0;
        int node_num = 0;
        int normal_num = 0;
        int order_max = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  OBJ_READ reads an object in an OBJ file.");

        WavefrontOBJ.obj_size(input_file_name, ref node_num, ref face_num, ref normal_num, ref order_max);

        int[] face_node = new int[order_max * face_num];
        int[] face_order = new int[face_num];
        double[] node_xyz = new double[3 * node_num];
        double[] normal_vector = new double[3 * normal_num];

        //obj_read ( input_file_name, node_num, face_num, normal_num, 
        //  order_max, node_xyz, face_order, face_node, normal_vector, vertex_normal );

        WavefrontOBJ.obj_face_node_print(face_num, order_max, face_order, face_node);
        WavefrontOBJ.obj_normal_vector_print(normal_num, normal_vector);
        WavefrontOBJ.obj_node_xyz_print(node_num, node_xyz);
    }

    private static void test03(string output_filename)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests OBJ_WRITE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 January 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] face_node =
        {
            1, 3, 2,
            2, 3, 4,
            1, 6, 5,
            1, 2, 6,
            3, 7, 4,
            4, 7, 8,
            5, 6, 8,
            5, 8, 7,
            1, 5, 7,
            1, 7, 3,
            2, 4, 6,
            6, 4, 8
        };
        const int face_num = 12;
        int[] face_order =
        {
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3
        };
        const int node_num = 8;
        double[] node_xyz =
        {
            0.0, 0.0, 0.0,
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            1.0, 1.0, 0.0,
            0.0, 0.0, 1.0,
            1.0, 0.0, 1.0,
            0.0, 1.0, 1.0,
            1.0, 1.0, 1.0
        };
        const int normal_num = 6;
        double[] normal_vector =
        {
            0.0, 0.0, 1.0,
            0.0, 0.0, -1.0,
            0.0, 1.0, 0.0,
            0.0, -1.0, 0.0,
            1.0, 0.0, 0.0,
            -1.0, 0.0, 0.0
        };
        const int order_max = 3;
        int[] vertex_normal =
        {
            2, 2, 2,
            2, 2, 2,
            4, 4, 4,
            4, 4, 4,
            3, 3, 3,
            3, 3, 3,
            1, 1, 1,
            1, 1, 1,
            6, 6, 6,
            6, 6, 6,
            5, 5, 5,
            5, 5, 5
        };

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  OBJ_WRITE writes an ASCII OBJ file.");

        WavefrontOBJ.obj_write(output_filename, node_num, face_num, normal_num,
            order_max, node_xyz, face_order, face_node, normal_vector, vertex_normal);

        Console.WriteLine("");
        Console.WriteLine("  Graphics data was written to the OBJ file \""
                          + output_filename + "\".");
    }
}