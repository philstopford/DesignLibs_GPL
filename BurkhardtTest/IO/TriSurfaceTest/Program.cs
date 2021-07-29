using System;
using Burkardt;

namespace TriSurfaceTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for TRI_SURFACE_IO_TEST.
            //
            //  Discussion:
            //
            //    TRI_SURFACE_IO_TEST tests the TRI_SURFACE_IO library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    28 September 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("TRI_SURFACE_IO_TEST:");
            Console.WriteLine("  Test the TRI_SURFACE_IO library.");

            test01("sphere_nodes.txt", "sphere_triangles.txt");
            test02("sphere_nodes.txt", "sphere_triangles.txt");
            test03("cube_nodes.txt", "cube_triangles.txt");

            Console.WriteLine("");
            Console.WriteLine("TRI_SURFACE_IO_TEST:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");

        }

        static void test01(string node_file_name, string triangle_file_name)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 tests TRI_SURFACE_SIZE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    28 September 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int dim_num = 0;
            int node_num = 0;
            int order_num = 0;
            int triangle_num = 0;

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  TRI_SURFACE_SIZE determines the size of various objects");
            Console.WriteLine("  in a TRI_SURFACE file.");

            TriSurface.tri_surface_size(node_file_name, triangle_file_name, ref dim_num, ref node_num,
                ref order_num, ref triangle_num);

            TriSurface.tri_surface_size_print(node_file_name, triangle_file_name, dim_num,
                node_num, order_num, triangle_num);

        }

        static void test02(string node_file_name, string triangle_file_name)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 tests TRI_SURFACE_READ.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 September 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int dim_num = 0;
            int node_num = 0;
            double[] node_xyz = null;
            int order_num = 0;
            int triangle_num = 0;
            int[] triangle_node = null;

            Console.WriteLine("");
            Console.WriteLine("TEST02");
            Console.WriteLine("  TRI_SURFACE_READ reads data from a TRI_SURFACE file.");

            TriSurface.tri_surface_size(node_file_name, triangle_file_name, ref dim_num, ref node_num,
                ref order_num, ref triangle_num);

            TriSurface.tri_surface_size_print(node_file_name, triangle_file_name, dim_num,
                node_num, order_num, triangle_num);

            TriSurface.tri_surface_read(node_file_name, triangle_file_name, dim_num, node_num,
                order_num, triangle_num, ref node_xyz, ref triangle_node);

            TriSurface.tri_surface_print(node_file_name, triangle_file_name, dim_num, node_num,
                order_num, triangle_num, node_xyz, triangle_node);
        }

        static void test03(string node_file_name, string triangle_file_name)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST03 tests TRI_SURFACE_WRITE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 September 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 3;
            int NODE_NUM = 8;
            int ORDER_NUM = 3;
            int TRIANGLE_NUM = 12;

            int dim_num = DIM_NUM;
            int node_num = NODE_NUM;
            int order_num = 3;
            int triangle_num = 12;

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
            int[] triangle_node =
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

            Console.WriteLine("");
            Console.WriteLine("TEST03");
            Console.WriteLine("  TRI_SURFACE_WRITE writes TRI_SURFACE data to two files.");

            TriSurface.tri_surface_write(node_file_name, triangle_file_name, dim_num, node_num,
                order_num, triangle_num, node_xyz, triangle_node);

            Console.WriteLine("");
            Console.WriteLine("  Graphics data was written to:");
            Console.WriteLine("    Node file:     \"" + node_file_name + "\".");
            Console.WriteLine("    Triangle file: \"" + triangle_file_name + "\".");
        }
    }
}