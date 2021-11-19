using System;
using Burkardt.Table;
using Burkardt.TriangleNS;
using Burkardt.TriangulationNS;
using Burkardt.Types;

namespace TableVoronoiTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TABLE_VORONOI.
        //
        //  Discussion:
        //
        //    TABLE_VORONOI uses GEOMPACK to get some Voronoi diagram information.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 May 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Max Gunzburger and John Burkardt,
        //    Uniformity Measures for Point Samples in Hypercubes.
        //
        //  Parameters:
        //
        //    Command line argument FILE_NAME, is the name of a file in the TABLE format
        //    containing the coordinates of a point set to be analyzed.
        //
    {
        int i;
        string file_name;

        Console.WriteLine("");

        Console.WriteLine("");
        Console.WriteLine("TABLE_VORONOI");
        Console.WriteLine("");
        Console.WriteLine("  This program is given the coordinates of a set of");
        Console.WriteLine("  points in the plane, calls GEOMPACK to determine the");
        Console.WriteLine("  Delaunay triangulation of those points, and then");
        Console.WriteLine("  digests that data to produce information defining");
        Console.WriteLine("  the Voronoi diagram.");
        Console.WriteLine("");
        Console.WriteLine("  The input file contains the following data:");
        Console.WriteLine("");
        Console.WriteLine("    G_NUM:    the number of generators;");
        Console.WriteLine("    G_XY:     the (X,Y) coordinates of the generators.");
        Console.WriteLine("");
        Console.WriteLine("  The computed Voronoi information includes:");
        Console.WriteLine("");
        Console.WriteLine("    G_DEGREE: the degree of each Voronoi cell;");
        Console.WriteLine("    G_START:  the index of the first Voronoi vertex;");
        Console.WriteLine("    G_FACE:   the list of all Voronoi vertices;");
        Console.WriteLine("");
        Console.WriteLine("    V_NUM:    the number of (finite) Voronoi vertices;");
        Console.WriteLine("    V_XY:     the (X,Y) coordinates of the Voronoi vertices;");
        Console.WriteLine("");
        Console.WriteLine("    I_NUM:    the number of Voronoi vertices at infinity;");
        Console.WriteLine("    I_XY:     the directions associated with the Voronoi");
        Console.WriteLine("              vertices at infinity.");
        //
        //  If the input file was not specified, get it now.
        //
        try
        {
            for (i = 1; i < args.Length; i++)
            {
                handle_file(args[i]);
            }
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("TABLE_VORONOI:");
            Console.WriteLine("  Please enter the name of a file to be analyzed.");

            file_name = Console.ReadLine();

            handle_file(file_name);

        }

        Console.WriteLine("");
        Console.WriteLine("TABLE_VORONOI:");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");
    }

    private static void handle_file(string file_name)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HANDLE_FILE handles a single file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 February 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, character FILE_NAME, the name of a file whose data
        //    is to be read and processed.
        //
    {
        int NDIM = 2;

        int[] g_degree;
        int[] g_face;
        int g_num;
        int[] g_start;
        double[] g_xy;
        int i_num = 0;
        double[] i_xy;
        int m;
        int v_num = 0;
        double[] v_xy;

        Console.WriteLine("");
        Console.WriteLine("HANDLE_FILE");
        Console.WriteLine("  Read the TABLE file \"" + file_name + "\".");

        TableHeader h = typeMethods.dtable_header_read(file_name);
        m = h.m;
        g_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  DTABLE_HEADER_READ has read the header.");
        Console.WriteLine("");
        Console.WriteLine("  The spatial dimension of the data M = " + m + "");
        Console.WriteLine("  The number of generators, G_NUM = " + g_num + "");

        if (m < NDIM)
        {
            Console.WriteLine("");
            Console.WriteLine("HANDLE - Fatal error!");
            Console.WriteLine("  The input spatial dimension is less than 2.");
            return;
        }

        if (NDIM < m)
        {
            Console.WriteLine("");
            Console.WriteLine("HANDLE - Warning:");
            Console.WriteLine("  The input spatial dimension is greater than 2.");
            Console.WriteLine("  Only the first two coordinates will be considered.");
        }

        i_xy = new double[NDIM * g_num];
        g_degree = new int[g_num];
        g_face = new int[6 * g_num];
        g_start = new int[g_num];
        v_xy = new double[NDIM * 2 * g_num];

        g_xy = typeMethods.dtable_data_read(file_name, NDIM, g_num);

        Console.WriteLine("");
        Console.WriteLine("  DTABLE_DATA_READ has read the data.");

        typeMethods.r8mat_transpose_print(NDIM, g_num, g_xy, "  The generators");

        voronoi_data(g_num, g_xy, ref g_degree, ref g_start, ref g_face, ref v_num, ref v_xy,
            ref i_num, ref i_xy);

        Console.WriteLine("");
        Console.WriteLine("  V_NUM: Number of Voronoi vertices = " + v_num + "");

        typeMethods.r8mat_transpose_print(NDIM, v_num, v_xy, "  Voronoi vertices:");

        Console.WriteLine("");
        Console.WriteLine("  I_NUM: Number of Voronoi vertices at infinity = " +
                          i_num + "");

        Console.WriteLine("");
        Console.WriteLine("  I_XY, the directions of the Voronoi vertices at infinity:");
        Console.WriteLine("");

        typeMethods.r8mat_transpose_print(NDIM, i_num, i_xy, "  Directions at infinity:");

    }

    private static void voronoi_data(int g_num, double[] g_xy, ref int[] g_degree, ref int[] g_start,
            ref int[] g_face, ref int v_num, ref double[] v_xy, ref int i_num, ref double[] i_xy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VORONOI_DATA returns data defining the Voronoi diagram.
        //
        //  Discussion:
        //
        //    The routine first determines the Delaunay triangulation.
        //
        //    The Voronoi diagram is then determined from this information.
        //
        //    In particular, the circumcenter of each Delaunay triangle
        //    is a vertex of a Voronoi polygon.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 February 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int G_NUM, the number of generators.
        //
        //    Input, double G_XY[2*G_NUM], the point coordinates.
        //
        //    Output, int G_DEGREE[G_NUM], the degree of each Voronoi cell.
        //
        //    Output, int G_START[G_NUM], the index in G_FACE of the first
        //    vertex at which to begin a traversal of the boundary of the 
        //    cell associated with each point.
        //
        //    Output, int G_FACE[6*G_NUM], the sequence of vertices to be used
        //    in a traversal of the boundary of the cell associated with each point.
        //
        //    Output, int *V_NUM, the number of vertices of the Voronoi diagram.
        //
        //    Output, double V_XY[2*V_NUM], the coordinates of the vertices
        //    of the Voronoi diagram.
        //
        //    Output, int *I_NUM, the number of vertices at infinity of the 
        //    Voronoi diagram.
        //
        //    Output, double I_XY[2*I_NUM], the direction of the
        //    vertices at infinity.
        //
    {
        int NDIM = 2;

        double[] center;
        int count;
        bool debug = true;
        int g;
        int g_next;
        int i;
        int i1;
        int i2;
        int i3;
        int ix1;
        int ix2;
        int j;
        int jp1;
        int k;
        int[] nodtri;
        double[] normal;
        int s;
        int sp1;
        int s_next = 0;
        double[] t = new double[NDIM * 3];
        int[] tnbr;
        int v;
        int v_inf;
        int v_next;
        int v_old;
        int v_save;
        //
        //  Compute the Delaunay triangulation.
        //
        nodtri = new int[3 * 2 * g_num];
        tnbr = new int[3 * 2 * g_num];

        Delauney.dtris2(g_num, 0, ref g_xy, ref v_num, ref nodtri, ref tnbr);
        //
        //  At this point, you could extend the NODTRI and TNBR data structures
        //  to account for the "vertices at infinity."
        //
        v_inf = Triangle.tri_augment(v_num, ref nodtri);

        switch (debug)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("  The generators that form each Delaunay triangle:");
                Console.WriteLine("  (Negative values are fictitious nodes at infinity.)");
                Console.WriteLine("");

                typeMethods.i4mat_transpose_print(3, v_num + v_inf, nodtri, "  Triangle nodes:");
                break;
        }

        //
        //  Negative entries in TNBR indicate a semi-infinite Voronoi side.
        //  However, DTRIS2 uses a peculiar numbering.  Renumber them.
        //
        i_num = 0;
        for (v = 0; v < v_num; v++)
        {
            for (i = 0; i < 3; i++)
            {
                switch (tnbr[i + v * 3])
                {
                    case < 0:
                        i_num += 1;
                        tnbr[i + v * 3] = -i_num;
                        break;
                }
            }
        }

        switch (debug)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("  Neighboring triangles of each Delaunay triangle:");
                Console.WriteLine("  Negative values indicate no finite neighbor.");
                Console.WriteLine("");

                typeMethods.i4mat_transpose_print(3, v_num, tnbr, "  Neighbor triangles:");
                break;
        }

        //
        //  Determine the degree of each cell.
        //
        for (g = 0; g < g_num; g++)
        {
            g_degree[g] = 0;
        }

        for (j = 0; j < v_num + v_inf; j++)
        {
            for (i = 0; i < 3; i++)
            {
                k = nodtri[i + j * 3];
                switch (k)
                {
                    case > 0:
                        g_degree[k - 1] += 1;
                        break;
                }
            }
        }

        switch (debug)
        {
            case true:
                typeMethods.i4vec_print(g_num, g_degree, "  Voronoi cell degrees:");
                break;
        }

        //
        //  Each (finite) Delaunay triangle contains a vertex of the Voronoi polygon,
        //  at the triangle's circumcenter.
        //
        for (j = 0; j < v_num; j++)
        {
            i1 = nodtri[0 + j * 3];
            i2 = nodtri[1 + j * 3];
            i3 = nodtri[2 + j * 3];

            t[0 + 0 * 2] = g_xy[0 + (i1 - 1) * 2];
            t[1 + 0 * 2] = g_xy[1 + (i1 - 1) * 2];
            t[0 + 1 * 2] = g_xy[0 + (i2 - 1) * 2];
            t[1 + 1 * 2] = g_xy[1 + (i2 - 1) * 2];
            t[0 + 2 * 2] = g_xy[0 + (i3 - 1) * 2];
            t[1 + 2 * 2] = g_xy[1 + (i3 - 1) * 2];

            center = typeMethods.triangle_circumcenter_2d(t);

            v_xy[0 + j * 2] = center[0];
            v_xy[1 + j * 2] = center[1];
        }

        switch (debug)
        {
            case true:
                typeMethods.r8mat_transpose_print(2, v_num, v_xy, "  The Voronoi vertices:");
                break;
        }

        //
        //  For each generator G:
        //    Determine if its region is infinite.
        //      Find a Delaunay triangle containing G.
        //      Seek another triangle containing the next node in that triangle.
        //
        count = 0;
        for (g = 0; g < g_num; g++)
        {
            g_start[g] = 0;
        }

        for (g = 1; g <= g_num; g++)
        {
            v_next = 0;
            for (v = 1; v <= v_num + v_inf; v++)
            {
                for (s = 1; s <= 3; s++)
                {
                    if (nodtri[s - 1 + (v - 1) * 3] == g)
                    {
                        v_next = v;
                        s_next = s;
                        break;
                    }
                }

                if (v_next != 0)
                {
                    break;
                }
            }

            v_save = v_next;

            for (;;)
            {
                s_next = typeMethods.i4_wrap(s_next + 1, 1, 3);
                g_next = nodtri[s_next - 1 + (v_next - 1) * 3];

                if (g_next == g)
                {
                    s_next = typeMethods.i4_wrap(s_next + 1, 1, 3);
                    g_next = nodtri[s_next - 1 + (v_next - 1) * 3];
                }

                v_old = v_next;
                v_next = 0;

                for (v = 1; v <= v_num + v_inf; v++)
                {
                    if (v == v_old)
                    {
                        continue;
                    }

                    for (s = 1; s <= 3; s++)
                    {
                        if (nodtri[s - 1 + (v - 1) * 3] == g)
                        {
                            sp1 = typeMethods.i4_wrap(s + 1, 1, 3);
                            if (nodtri[sp1 - 1 + (v - 1) * 3] == g_next)
                            {
                                v_next = v;
                                s_next = sp1;
                                break;
                            }

                            sp1 = typeMethods.i4_wrap(s + 2, 1, 3);
                            if (nodtri[sp1 - 1 + (v - 1) * 3] == g_next)
                            {
                                v_next = v;
                                s_next = sp1;
                                break;
                            }
                        }
                    }

                    if (v_next != 0)
                    {
                        break;
                    }
                }

                if (v_next == v_save)
                {
                    break;
                }

                if (v_next == 0)
                {
                    v_next = v_old;
                    break;
                }
            }

            //
            //  Now, starting in the current triangle, V_NEXT, cycle again,
            //  and copy the list of nodes into the array.
            //
            v_save = v_next;

            count += 1;
            g_start[g - 1] = count;
            g_face[count - 1] = v_next;

            for (;;)
            {
                s_next = typeMethods.i4_wrap(s_next + 1, 1, 3);
                g_next = nodtri[s_next - 1 + (v_next - 1) * 3];

                if (g_next == g)
                {
                    s_next = typeMethods.i4_wrap(s_next + 1, 1, 3);
                    g_next = nodtri[s_next - 1 + (v_next - 1) * 3];
                }

                v_old = v_next;
                v_next = 0;
                for (v = 1; v <= v_num + v_inf; v++)
                {
                    if (v == v_old)
                    {
                        continue;
                    }

                    for (s = 1; s <= 3; s++)
                    {
                        if (nodtri[s - 1 + (v - 1) * 3] == g)
                        {
                            sp1 = typeMethods.i4_wrap(s + 1, 1, 3);
                            if (nodtri[sp1 - 1 + (v - 1) * 3] == g_next)
                            {
                                v_next = v;
                                s_next = sp1;
                                break;
                            }

                            sp1 = typeMethods.i4_wrap(s + 2, 1, 3);
                            if (nodtri[sp1 - 1 + (v - 1) * 3] == g_next)
                            {
                                v_next = v;
                                s_next = sp1;
                                break;
                            }
                        }
                    }

                    if (v_next != 0)
                    {
                        break;
                    }
                }

                if (v_next == v_save)
                {
                    break;
                }

                if (v_next == 0)
                {
                    break;
                }

                count += 1;
                g_face[count - 1] = v_next;
            }
        }

        //
        //  Mark all the vertices at infinity with a negative sign,
        //  so that the data in G_FACE is easier to interpret.
        //
        for (i = 0; i < count; i++)
        {
            if (v_num < g_face[i])
            {
                g_face[i] = -g_face[i];
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  G_START: The index of the first Voronoi vertex");
        Console.WriteLine("  G_FACE: The Voronoi vertices");
        Console.WriteLine("");
        Console.WriteLine("   G  G_START  G_FACE");
        Console.WriteLine("");

        for (j = 0; j < g_num; j++)
        {
            k = g_start[j] - 1;
            Console.WriteLine("");
            Console.WriteLine("  "
                              + (j + 1).ToString().PadLeft(4) + "  "
                              + (k + 1).ToString().PadLeft(4) + "  "
                              + g_face[k].ToString().PadLeft(4) + "");
            for (i = 1; i < g_degree[j]; i++)
            {
                k += 1;
                Console.WriteLine("              "
                                  + g_face[k].ToString().PadLeft(4) + "");
            }
        }

        //
        //  For each (finite) Delaunay triangle, I
        //    For each side J,
        //
        for (i = 0; i < v_num; i++)
        {
            for (j = 0; j < 3; j++)
            {
                k = tnbr[j + i * 3];
                switch (k)
                {
                    //
                    //      If there is no neighboring triangle on that side,
                    //        extend a line from the circumcenter of I in the direction of the
                    //        outward normal to that side.  This is an infinite edge of
                    //        an infinite Voronoi polygon.
                    //
                    case < 0:
                        ix1 = nodtri[j + i * 3];
                        //      x1 = g_xy[0+(ix1-1)*2];
                        //      y1 = g_xy[1+(ix1-1)*2];
                        jp1 = typeMethods.i4_wrap(j + 1, 0, 2);

                        ix2 = nodtri[jp1 + i * 3];
                        //      x2 = g_xy[0+(ix2-1)*2];
                        //      y2 = g_xy[1+(ix2-1)*2];
                        //
                        //  Compute the direction I_XY(1:2,-K).
                        //
                        normal = Burkardt.LineNS.Geometry.line_exp_normal_2d(g_xy, g_xy,  p1Index:+ (ix1 - 1) * 2, p2Index:+ (ix2 - 1) * 2);

                        i_xy[0 + (-k - 1) * 2] = normal[0];
                        i_xy[1 + (-k - 1) * 2] = normal[1];
                        break;
                }
            }
        }
    }

}