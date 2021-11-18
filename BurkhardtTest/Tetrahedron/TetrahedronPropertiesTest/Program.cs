using System;
using Burkardt.Table;
using Burkardt.TetrahedronNS;
using Burkardt.Types;

namespace TetrahedronPropertiesTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TETRAHEDRON_PROPERTIES.
        //
        //
        //  Discussion:
        //
        //    TETRAHEDRON_PROPERTIES reports properties of a tetrahedron.
        //
        //  Usage:
        //
        //    tetrahedron_properties filename
        //
        //    where "filename" is a file containing the coordinates of the vertices.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] ab = new double[3];
        double[] ac = new double[3];
        double[] ad = new double[3];
        double[] bc = new double[3];
        double[] bd = new double[3];
        double[] cd = new double[3];
        double[] centroid;
        double[] circum_center = new double[3];
        double circum_radius = 0;
        double[] dihedral_angles;
        int dim_num;
        double[] edge_length;
        double[] face_angles = new double[3 * 4];
        double[] face_areas = new double[4];
        int i;
        int j;
        double[] in_center = new double[3];
        double in_radius = 0;
        string node_filename;
        int node_num;
        double[] node_xyz;
        double quality1;
        double quality2;
        double quality3;
        double quality4;
        const double r8_pi = 3.141592653589793;
        double[] solid_angles;
        double volume;

        Console.WriteLine("");

        try
        {
            node_filename = args[0];
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  Please enter the name of the node coordinate file.");
            node_filename = Console.ReadLine();
        }

        //
        //  Read the node data.
        //
        TableHeader h = typeMethods.r8mat_header_read(node_filename);
        dim_num = h.m;
        node_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Read the header of \"" + node_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM = " + dim_num + "");
        Console.WriteLine("  Number of points NODE_NUM = " + node_num + "");

        Console.WriteLine("");
        Console.WriteLine("TETRAHEDRON_PROPERTIES:");
        Console.WriteLine("  C++ version:");
        Console.WriteLine("  Determine properties of a tetrahedron.");

        if (dim_num != 3)
        {
            Console.WriteLine("");
            Console.WriteLine("TETRAHEDRON_PROPERTIES - Fatal error!");
            Console.WriteLine("  Dataset must have spatial dimension 3.");
            return;
        }

        if (node_num != 4)
        {
            Console.WriteLine("");
            Console.WriteLine("TETRAHEDRON_PROPERTIES - Fatal error!");
            Console.WriteLine("  Dataset must have 4 nodes.");
            return;
        }

        node_xyz = typeMethods.r8mat_data_read(node_filename, dim_num, node_num);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + node_filename + "\".");

        typeMethods.r8mat_transpose_print(dim_num, node_num, node_xyz, "  Node coordinates:");
        //
        //  CIRCUMSPHERE
        //
        Properties.tetrahedron_circumsphere(node_xyz, ref circum_radius, ref circum_center);

        Console.WriteLine("");
        Console.WriteLine("  CIRCUM_RADIUS = " + circum_radius + "");
        Console.WriteLine("  CIRCUM_CENTER: "
                          + "  " + circum_center[0].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + circum_center[1].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + circum_center[2].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        //
        //  CENTROID
        //
        centroid = Properties.tetrahedron_centroid(node_xyz);

        Console.WriteLine("");
        Console.WriteLine("  CENTROID: "
                          + "  " + centroid[0].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + centroid[1].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + centroid[2].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        //
        //  DIHEDRAL ANGLES
        //
        dihedral_angles = Properties.tetrahedron_dihedral_angles(node_xyz);

        typeMethods.r8vec_print(6, dihedral_angles, "  DIHEDRAL_ANGLES (radians)");

        for (i = 0; i < 6; i++)
        {
            dihedral_angles[i] = dihedral_angles[i] * 180.0 / r8_pi;
        }

        typeMethods.r8vec_print(6, dihedral_angles, "  DIHEDRAL_ANGLES (degrees)");
        //
        //  EDGES
        //
        Properties.tetrahedron_edges(node_xyz, ref ab, ref ac, ref ad, ref bc, ref bd, ref cd);

        Console.WriteLine("");
        typeMethods.r8vec_transpose_print(3, ab, "  EDGE AB:");
        typeMethods.r8vec_transpose_print(3, ac, "  EDGE AC:");
        typeMethods.r8vec_transpose_print(3, ad, "  EDGE AD:");
        typeMethods.r8vec_transpose_print(3, bc, "  EDGE BC:");
        typeMethods.r8vec_transpose_print(3, bd, "  EDGE BD:");
        typeMethods.r8vec_transpose_print(3, cd, "  EDGE CD:");
        //
        //  EDGE LENGTHS
        //
        edge_length = Properties.tetrahedron_edge_length(node_xyz);

        typeMethods.r8vec_print(6, edge_length, "  EDGE_LENGTHS");
        //
        //  FACE ANGLES
        //
        Properties.tetrahedron_face_angles(node_xyz, ref face_angles);

        typeMethods.r8mat_transpose_print(3, 4, face_angles, "  FACE_ANGLES (radians)");

        for (j = 0; j < 4; j++)
        {
            for (i = 0; i < 3; i++)
            {
                face_angles[i + j * 3] = face_angles[i + j * 3] * 180.0 / r8_pi;
            }
        }

        typeMethods.r8mat_transpose_print(3, 4, face_angles, "  FACE_ANGLES (degrees)");
        //
        //  FACE AREAS
        //
        Properties.tetrahedron_face_areas(node_xyz, ref face_areas);

        typeMethods.r8vec_print(4, face_areas, "  FACE_AREAS");
        //
        //  INSPHERE
        //
        Properties.tetrahedron_insphere(node_xyz, ref in_radius, ref in_center);

        Console.WriteLine("");
        Console.WriteLine("  IN_RADIUS = " + in_radius + "");
        Console.WriteLine("  IN_CENTER: "
                          + "  " + in_center[0].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + in_center[1].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                          + "  " + in_center[2].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        //
        //  QUALITY1
        //
        quality1 = Properties.tetrahedron_quality1(node_xyz);

        Console.WriteLine("");
        Console.WriteLine("  QUALITY1 = " + quality1 + "");
        //
        //  QUALITY2
        //
        quality2 = Properties.tetrahedron_quality2(node_xyz);

        Console.WriteLine("");
        Console.WriteLine("  QUALITY2 = " + quality2 + "");
        //
        //  QUALITY3
        //
        quality3 = Properties.tetrahedron_quality3(node_xyz);

        Console.WriteLine("");
        Console.WriteLine("  QUALITY3 = " + quality3 + "");
        //
        //  QUALITY4
        //
        quality4 = Properties.tetrahedron_quality4(node_xyz);

        Console.WriteLine("");
        Console.WriteLine("  QUALITY4 = " + quality4 + "");
        //
        //  SOLID ANGLES
        //
        solid_angles = Properties.tetrahedron_solid_angles(node_xyz);

        typeMethods.r8vec_print(4, solid_angles, "  SOLID_ANGLES (steradians)");
        //
        //  VOLUME
        //
        volume = Tetrahedron.tetrahedron_volume(node_xyz);

        Console.WriteLine("");
        Console.WriteLine("  VOLUME = " + volume + "");

        Console.WriteLine("");
        Console.WriteLine("TETRAHEDRON_PROPERTIES:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}