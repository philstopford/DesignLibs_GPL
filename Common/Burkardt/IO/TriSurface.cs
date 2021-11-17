using System;
using Burkardt.Table;
using Burkardt.Types;

namespace Burkardt.IO;

public static class TriSurface
{
    public static void tri_surface_print(string node_file_name, string triangle_file_name,
            int dim_num, int node_num, int order_num, int triangle_num,
            double[] node_xyz, int[] triangle_node)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRI_SURFACE_PRINT prints graphics information from a pair of TRI_SURFACE files.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 September 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string NODE_FILE_NAME, the name of the node file.
        //
        //    Input, string TRIANGLE_FILE_NAME, the name of the triangle file.
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int NODE_NUM, the number of points.
        //
        //    Input, int ORDER_NUM, the order of the triangles.
        //
        //    Input, int TRIANGLE_NUM, the number of triangles.
        //
        //    Input, double NODE_XYZ[DIM_NUM*NODE_NUM], the node coordinates.
        //
        //    Input, int TRIANGLE_NODE[ORDER_NUM*TRIANGLE_NUM], 
        //    the nodes that form the triangles.
        //
    {
        typeMethods.r8mat_transpose_print(dim_num, node_num, node_xyz, "  Node coordinates");

        typeMethods.i4mat_transpose_print(order_num, triangle_num, triangle_node,
            "  Triangle nodes");
    }

    public static void tri_surface_read(string node_file_name, string triangle_file_name,
            int dim_num, int node_num, int order_num, int triangle_num,
            ref double[] node_xyz, ref int[] triangle_node)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRI_SURFACE_READ reads graphics information from a pair of TRI_SURFACE files.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 September 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string NODE_FILE_NAME, the name of the node file.
        //
        //    Input, string TRIANGLE_FILE_NAME, the name of the triangle file.
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int NODE_NUM, the number of points.
        //
        //    Input, int ORDER_NUM, the order of the triangles.
        //
        //    Input, int TRIANGLE_NUM, the number of triangles.
        //
        //    Output, double NODE_XYZ[DIM_NUM*NODE_NUM], the node coordinates.
        //
        //    Output, int TRIANGLE_NODE[ORDER_NUM*TRIANGLE_NUM], 
        //    the nodes that form the triangles.
        //
    {
        node_xyz = typeMethods.r8mat_data_read(node_file_name, dim_num, node_num);

        triangle_node = typeMethods.i4mat_data_read(triangle_file_name, order_num,
            triangle_num);
    }

    public static void tri_surface_size(string node_file_name, string triangle_file_name,
            ref int dim_num, ref int node_num, ref int order_num, ref int triangle_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRI_SURFACE_SIZE determines the size of a TRI_SURFACE object.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 September 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string NODE_FILE_NAME, the name of the node file.
        //
        //    Input, string TRIANGLE_FILE_NAME, the name of the triangle file.
        //
        //    Output, int *DIM_NUM, the spatial dimension.
        //
        //    Output, int *NODE_NUM, the number of nodes.
        //
        //    Output, int *ORDER_NUM, the order of the triangles.
        //
        //    Output, int *TRIANGLE_NUM, the number of triangles.
        //
    {
        dim_num = -1;
        node_num = -1;
        triangle_num = -1;

        TableHeader h = typeMethods.r8mat_header_read(node_file_name);
        dim_num = h.m;
        node_num = h.n;

        switch (dim_num)
        {
            case < 2:
            case > 3:
                Console.WriteLine("");
                Console.WriteLine("TRI_SURFACE_SIZE - Warning!");
                Console.WriteLine("  The spatial dimension DIM_NUM = " + dim_num + "");
                Console.WriteLine("  This seems an unlikely value.");
                break;
        }

        h = typeMethods.i4mat_header_read(triangle_file_name);
        order_num = h.m;
        triangle_num = h.n;

        if (order_num != 3 && order_num != 6)
        {
            Console.WriteLine("");
            Console.WriteLine("TRI_SURFACE_SIZE - Fatal error!");
            Console.WriteLine("  The order of the triangles seems to be " + order_num + "");
            Console.WriteLine("  Only the values 3 and 6 are acceptable.");
        }
    }

    public static void tri_surface_size_print(string node_file_name, string triangle_file_name,
            int dim_num, int node_num, int order_num, int triangle_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRI_SURFACE_SIZE_PRINT prints sizes associated with a TRI_SURFACE file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 September 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string NODE_FILE_NAME, the name of the node file.
        //
        //    Input, string TRIANGLE_FILE_NAME, the name of the triangle file.
        //
        //    Input, int DIM_NUM, the number of spatial dimensions.
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, int ORDER_NUM, the order of the triangles.
        //
        //    Input, int TRIANGLE_NUM, the number of triangles.
        //
    {
        Console.WriteLine("");
        Console.WriteLine("TRI_SURFACE_SIZE_PRINT:");
        Console.WriteLine("");
        Console.WriteLine("  Node file     \"" + node_file_name + "\".");
        Console.WriteLine("  Triangle file \"" + triangle_file_name + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension  = " + dim_num + "");
        Console.WriteLine("  Nodes              = " + node_num + "");
        Console.WriteLine("  Triangle order     = " + order_num + "");
        Console.WriteLine("  Triangles          = " + triangle_num + "");
    }

    public static void tri_surface_write(string node_file_name, string triangle_file_name,
            int dim_num, int node_num, int order_num, int triangle_num,
            double[] node_xyz, int[] triangle_node)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRI_SURFACE_WRITE writes graphics information to a pair of TRI_SURFACE files.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 September 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string NODE_FILE_NAME, the name of the node file.
        //
        //    Input, string TRIANGLE_FILE_NAME, the name of the triangle file.
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int NODE_NUM, the number of points.
        //
        //    Input, int ORDER_NUM, the order of the triangles.
        //
        //    Input, int TRIANGLE_NUM, the number of triangles.
        //
        //    Input, double NODE_XYZ[DIM_NUM*NODE_NUM], the node coordinates.
        //
        //    Input, int TRIANGLE_NODE[ORDER_NUM*TRIANGLE_NUM], 
        //    the nodes that form the triangles.
        //
    {
        typeMethods.r8mat_write(node_file_name, dim_num, node_num, node_xyz);

        typeMethods.i4mat_write(triangle_file_name, order_num, triangle_num, triangle_node);

    }
}