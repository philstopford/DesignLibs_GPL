using System;
using Burkardt.IO;

namespace MeshIOTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN tests the MESH_IO library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    25 October 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            string filename;

            Console.WriteLine("");
            Console.WriteLine("MESH_IO_TEST");
            Console.WriteLine("  Test the MESH_IO library.");
            //
            //  Create the file hexahexa_2x2x2.mesh
            //
            test01();
            //
            //  Read and print the sizes of file hexahexa_2x2x2.mesh.
            //
            filename = "hexahexa_2x2x2.mesh";
            test03(filename);
            //
            //  Create the file cyl248.mesh
            //
            test02();
            //
            //  Read and print the sizes of file cyl248.mesh.
            //
            filename = "cyl248.mesh";
            test03(filename);
            //
            //  Read and print the data in file cyl248.mesh.
            //
            filename = "cyl248.mesh";
            test04(filename);
            //
            //  Terminate.
            //
            Console.WriteLine("");
            Console.WriteLine("MESH_IO_TEST");
            Console.WriteLine("  Normal end of execution.");

            Console.WriteLine("");
        }

        static void test01()

            //****************************************************************************80
            //
            //// TEST01 creates a MESH dataset and writes it to a file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 October 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int dim = 0;
            int[] edge_label;
            int[] edge_vertex;
            int edges = 0;
            string filename;
            int[] hexahedron_label;
            int[] hexahedron_vertex;
            int hexahedrons = 0;
            int[] quadrilateral_label;
            int[] quadrilateral_vertex;
            int quadrilaterals = 0;
            int[] tetrahedron_label;
            int[] tetrahedron_vertex;
            int tetrahedrons = 0;
            int[] triangle_label;
            int[] triangle_vertex;
            int triangles = 0;
            double[] vertex_coordinate;
            int[] vertex_label;
            int vertices = 0;

            Console.WriteLine("");
            Console.WriteLine("TEST01:");
            Console.WriteLine("  Create a hexahedral mesh and write it to a file.");
            //
            //  Get sizes.
            //
            Mesh.hexahexa_2x2x2_size(ref dim, ref vertices, ref edges, ref triangles, ref quadrilaterals,
                ref tetrahedrons, ref hexahedrons);
            //
            //  Allocate memory.
            //
            edge_label = new int[edges];
            edge_vertex = new int[2 * edges];
            hexahedron_label = new int[hexahedrons];
            hexahedron_vertex = new int[8 * hexahedrons];
            quadrilateral_label = new int[quadrilaterals];
            quadrilateral_vertex = new int[4 * quadrilaterals];
            tetrahedron_label = new int[tetrahedrons];
            tetrahedron_vertex = new int[4 * tetrahedrons];
            triangle_label = new int[triangles];
            triangle_vertex = new int[3 * triangles];
            vertex_coordinate = new double[3 * vertices];
            vertex_label = new int[vertices];
            //
            //  Get the data.
            //
            Mesh.hexahexa_2x2x2_data(dim, vertices, edges, triangles,
                quadrilaterals, tetrahedrons, hexahedrons, ref vertex_coordinate, ref vertex_label,
                ref edge_vertex, ref edge_label, ref triangle_vertex, ref triangle_label,
                ref quadrilateral_vertex, ref quadrilateral_label, ref tetrahedron_vertex,
                ref tetrahedron_label, ref hexahedron_vertex, ref hexahedron_label);
            //
            //  Write the data.
            //
            filename = "hexahexa_2x2x2.mesh";

            Mesh.mesh_write(filename, dim, vertices, edges, triangles,
                quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate,
                vertex_label, edge_vertex, edge_label, triangle_vertex, triangle_label,
                quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex,
                tetrahedron_label, hexahedron_vertex, hexahedron_label);

            Console.WriteLine("");
            Console.WriteLine("  Created the file \"" + filename + "\".");
        }

        static void test02()

            //****************************************************************************80
            //
            //// TEST02 creates a MESH dataset and writes it to a file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 October 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int dim = 0;
            int[] edge_label;
            int[] edge_vertex;
            int edges = 0;
            string filename;
            int[] hexahedron_label;
            int[] hexahedron_vertex;
            int hexahedrons = 0;
            int[] quadrilateral_label;
            int[] quadrilateral_vertex;
            int quadrilaterals = 0;
            int[] tetrahedron_label;
            int[] tetrahedron_vertex;
            int tetrahedrons = 0;
            int[] triangle_label;
            int[] triangle_vertex;
            int triangles = 0;
            double[] vertex_coordinate;
            int[] vertex_label;
            int vertices = 0;

            Console.WriteLine("");
            Console.WriteLine("TEST02:");
            Console.WriteLine("  Create a tetrahedral mesh and write it to a file.");
            //
            //  Get sizes.
            //
            Mesh.cyl248_size(ref dim, ref vertices, ref edges, ref triangles, ref quadrilaterals,
                ref tetrahedrons, ref hexahedrons);
            //
            //  Allocate memory.
            //
            edge_label = new int[edges];
            edge_vertex = new int[2 * edges];
            hexahedron_label = new int[hexahedrons];
            hexahedron_vertex = new int[8 * hexahedrons];
            quadrilateral_label = new int[quadrilaterals];
            quadrilateral_vertex = new int[4 * quadrilaterals];
            tetrahedron_label = new int[tetrahedrons];
            tetrahedron_vertex = new int[4 * tetrahedrons];
            triangle_label = new int[triangles];
            triangle_vertex = new int[3 * triangles];
            vertex_coordinate = new double[3 * vertices];
            vertex_label = new int[vertices];
            //
            //  Get the data.
            //
            Mesh.cyl248_data(dim, vertices, edges, triangles,
                quadrilaterals, tetrahedrons, hexahedrons, ref vertex_coordinate, ref vertex_label,
                ref edge_vertex, ref edge_label, ref triangle_vertex, ref triangle_label,
                ref quadrilateral_vertex, ref quadrilateral_label, ref tetrahedron_vertex,
                ref tetrahedron_label, ref hexahedron_vertex, ref hexahedron_label);
            //
            //  Write the data.
            //
            filename = "cyl248.mesh";

            Mesh.mesh_write(filename, dim, vertices, edges, triangles,
                quadrilaterals, tetrahedrons, hexahedrons, vertex_coordinate,
                vertex_label, edge_vertex, edge_label, triangle_vertex, triangle_label,
                quadrilateral_vertex, quadrilateral_label, tetrahedron_vertex,
                tetrahedron_label, hexahedron_vertex, hexahedron_label);

            Console.WriteLine("");
            Console.WriteLine("  Created the file \"" + filename + "\".");
        }

        static void test03(string filename)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MESH_IO_TEST03 reads and prints the sizes in a MESH dataset.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    25 October 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int dim = 0;
            int edges = 0;
            int hexahedrons = 0;
            int quadrilaterals = 0;
            int tetrahedrons = 0;
            int triangles = 0;
            int vertices = 0;

            Console.WriteLine("");
            Console.WriteLine("MESH_IO_TEST03");
            Console.WriteLine("  Read a mesh file and print its sizes.");
            //
            //  Read sizes.
            //
            Mesh.mesh_size_read(filename, ref dim, ref vertices, ref edges, ref triangles,
                ref quadrilaterals, ref tetrahedrons, ref hexahedrons);
            //
            //  Print sizes.
            //
            Console.WriteLine("");
            Console.WriteLine("  Header information for \"" + filename + "\"");

            Mesh.mesh_size_print(dim, vertices, edges, triangles, quadrilaterals,
                tetrahedrons, hexahedrons);

        }

        static void test04(string filename)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MESH_IO_TEST04 reads a MESH dataset and prints its data.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    25 October 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int dim = 0;
            int[] edge_label;
            int[] edge_vertex;
            int edges = 0;
            int[] hexahedron_label;
            int[] hexahedron_vertex;
            int hexahedrons = 0;
            int[] quadrilateral_label;
            int[] quadrilateral_vertex;
            int quadrilaterals = 0;
            int[] tetrahedron_label;
            int[] tetrahedron_vertex;
            int tetrahedrons = 0;
            int[] triangle_label;
            int[] triangle_vertex;
            int triangles = 0;
            double[] vertex_coordinate;
            int[] vertex_label;
            int vertices = 0;

            Console.WriteLine("");
            Console.WriteLine("MESH_IO_TEST04");
            Console.WriteLine("  Read a mesh file and print its data.");
            //
            //  Read sizes.
            //
            Mesh.mesh_size_read(filename, ref dim, ref vertices, ref edges, ref triangles,
                ref quadrilaterals, ref tetrahedrons, ref hexahedrons);
            //
            //  Allocate memory.
            //
            edge_label = new int[edges];
            edge_vertex = new int[2 * edges];
            hexahedron_label = new int[hexahedrons];
            hexahedron_vertex = new int[8 * hexahedrons];
            quadrilateral_label = new int[quadrilaterals];
            quadrilateral_vertex = new int[4 * quadrilaterals];
            tetrahedron_label = new int[tetrahedrons];
            tetrahedron_vertex = new int[4 * tetrahedrons];
            triangle_label = new int[triangles];
            triangle_vertex = new int[3 * triangles];
            vertex_coordinate = new double[3 * vertices];
            vertex_label = new int[vertices];
            //
            //  Read the data.
            //
            Mesh.mesh_data_read(filename, dim, vertices, edges, triangles,
                quadrilaterals, tetrahedrons, hexahedrons, ref vertex_coordinate,
                ref vertex_label, ref edge_vertex, ref edge_label, ref triangle_vertex, ref triangle_label,
                ref quadrilateral_vertex, ref quadrilateral_label, ref tetrahedron_vertex,
                ref tetrahedron_label, ref hexahedron_vertex, ref hexahedron_label);
            //
            //  Print the data.
            //
            Console.WriteLine("");
            Console.WriteLine("  Data for file \"" + filename + "\".");

            Mesh.mesh_data_print(dim, vertices, edges, triangles, quadrilaterals,
                tetrahedrons, hexahedrons, vertex_coordinate, vertex_label, edge_vertex,
                edge_label, triangle_vertex, triangle_label, quadrilateral_vertex,
                quadrilateral_label, tetrahedron_vertex, tetrahedron_label,
                hexahedron_vertex, hexahedron_label);
        }
    }
}