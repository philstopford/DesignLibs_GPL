using System;
using Burkardt.IO;
using Burkardt.Types;

namespace MeditToFEM;

internal static class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for MEDIT_TO_FEM.
        //
        //  Discussion:
        //
        //    MEDIT_TO_FEM converts mesh data from MEDIT to FEM format.
        //
        //  Usage:
        //
        //    medit_to_fem prefix
        //
        //    where 'prefix' is the common filename prefix:
        //
        //    * 'prefix'.mesh is the MEDIT mesh file.
        //    * 'prefix'_nodes.txt will contain the node coordinates.
        //    * 'prefix'_elements.txt will contain the element node connectivity.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 November 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int dim = 0;
        int edges = 0;
        int element_num;
        int element_order;
        int hexahedrons = 0;
        string prefix;
        int quadrilaterals = 0;
        int tetrahedrons = 0;
        int triangles = 0;
        int vertices = 0;

        Console.WriteLine("");
        Console.WriteLine("MEDIT_TO_FEM");
        Console.WriteLine("  Read a mesh description created by the MEDIT program:");
        Console.WriteLine("  * 'prefix'.mesh, the MEDIT mesh file.");
        Console.WriteLine("  Write out two simple FEM files listing nodes and elements.");
        Console.WriteLine("  * 'prefix'_nodes.txt, node coordinates.");
        Console.WriteLine("  * 'prefix'_elements.txt, element connectivity.");
        //
        //  Get the filename prefix.
        //
        try
        {
            prefix = args[0];
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  Please enter the filename prefix.");

            prefix = Console.ReadLine();
        }

        //
        //  Create the filenames.
        //
        string medit_filename = prefix + ".mesh";
        string fem_node_filename = prefix + "_nodes.txt";
        string fem_element_filename = prefix + "_elements.txt";
        //
        //  Read MEDIT sizes.
        //
        Mesh.mesh_size_read(medit_filename, ref dim, ref vertices, ref edges, ref triangles,
            ref quadrilaterals, ref tetrahedrons, ref hexahedrons);
        //
        //  Report sizes.
        //
        Console.WriteLine("");
        Console.WriteLine("  Number of dimensions = " + dim + "");
        Console.WriteLine("  Number of vertices = " + vertices + "");
        Console.WriteLine("  Number of edges = " + edges + "");
        Console.WriteLine("  Number of triangles = " + triangles + "");
        Console.WriteLine("  Number of quadrilaterals = " + quadrilaterals + "");
        Console.WriteLine("  Number of tetrahedrons = " + tetrahedrons + "");
        Console.WriteLine("  Number of hexahedrons = " + hexahedrons + "");
        //
        //  Allocate memory.
        //
        int[] edge_label = new int[edges];
        int[] edge_vertex = new int[2 * edges];
        int[] hexahedron_label = new int[hexahedrons];
        int[] hexahedron_vertex = new int[8 * hexahedrons];
        int[] quadrilateral_label = new int[quadrilaterals];
        int[] quadrilateral_vertex = new int[4 * quadrilaterals];
        int[] tetrahedron_label = new int[tetrahedrons];
        int[] tetrahedron_vertex = new int[4 * tetrahedrons];
        int[] triangle_label = new int[triangles];
        int[] triangle_vertex = new int[3 * triangles];
        double[] vertex_coordinate = new double[dim * vertices];
        int[] vertex_label = new int[vertices];
        //
        //  Read MEDIT data.
        //
        Mesh.mesh_data_read(medit_filename, dim, vertices, edges, triangles,
            quadrilaterals, tetrahedrons, hexahedrons, ref vertex_coordinate,
            ref vertex_label, ref edge_vertex, ref edge_label, ref triangle_vertex, ref triangle_label,
            ref quadrilateral_vertex, ref quadrilateral_label, ref tetrahedron_vertex,
            ref tetrahedron_label, ref hexahedron_vertex, ref hexahedron_label);
        //
        //  Choose the FEM data.
        //
        //  We need to assume that there is only one element type.
        //  If there are elements of multiple dimension, take the highest.
        //
        //m = dim;
        //node_num = vertices;
        typeMethods.r8mat_write(fem_node_filename, dim, vertices, vertex_coordinate);

        Console.WriteLine("");
        Console.WriteLine("  Created node coordinate file '" + fem_node_filename + "'");

        switch (hexahedrons)
        {
            case > 0 when dim == 3:
                element_order = 8;
                element_num = hexahedrons;
                typeMethods.i4mat_write(fem_element_filename, element_order, element_num,
                    hexahedron_vertex);
                break;
            default:
            {
                switch (tetrahedrons)
                {
                    case > 0 when dim == 3:
                        element_order = 4;
                        element_num = tetrahedrons;
                        typeMethods.i4mat_write(fem_element_filename, element_order, element_num,
                            tetrahedron_vertex);
                        break;
                    default:
                    {
                        switch (quadrilaterals)
                        {
                            case > 0 when dim == 2:
                                element_order = 4;
                                element_num = quadrilaterals;
                                typeMethods.i4mat_write(fem_element_filename, element_order, element_num,
                                    quadrilateral_vertex);
                                break;
                            default:
                            {
                                switch (triangles)
                                {
                                    case > 0 when dim == 2:
                                        element_order = 3;
                                        element_num = triangles;
                                        typeMethods.i4mat_write(fem_element_filename, element_order, element_num,
                                            triangle_vertex);
                                        break;
                                    default:
                                        Console.WriteLine("");
                                        Console.WriteLine("MEDIT_TO_FEM - Fatal error!");
                                        Console.WriteLine("  Unexpected values for spatial dimension");
                                        Console.WriteLine("  and number of nonzero objects.");
                                        return;
                                }

                                break;
                            }
                        }

                        break;
                    }
                }

                break;
            }
        }

        Console.WriteLine("  Created element connectivity file '" + fem_element_filename + "'");

        Console.WriteLine("");
        Console.WriteLine("MEDIT_TO_FEM:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}