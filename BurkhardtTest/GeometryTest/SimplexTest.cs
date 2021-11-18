using System;
using Burkardt.SimplexNS;
using Burkardt.Types;

namespace GeometryTest;

public static class SimplexTest
{
    public static void test1788()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST17888 tests SIMPLEX_LATTICE_LAYER_POINT_NEXT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int TEST_NUM = 4;

        int[] c;
        int i;
        int j;
        int layer;
        bool more;
        int n;
        int[] n_test = {1, 2, 3, 4};
        int test;
        int[] v;
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine("TEST1788");
        Console.WriteLine("  SIMPLEX_LATTICE_LAYER_POINT_NEXT returns the next");
        Console.WriteLine("  point in an N-dimensional simplex lattice layer defined by:");
        Console.WriteLine("");
        Console.WriteLine("    C(N+1) - 1 <= X[0]/C[0] + X(2)/C(2) + ... + X(N)/C(N) <= C(N+1).");

        for (test = 0; test < TEST_NUM; test++)
        {
            n = n_test[test];
            c = new int[n + 1];
            v = new int[n];

            for (i = 0; i < n + 1; i++)
            {
                c[i] = i + 2;
            }

            for (i = 0; i < n; i++)
            {
                v[i] = 0;
            }

            Console.WriteLine("");
            Console.WriteLine("  N = " + n + "");
            cout = "  C = ";
            for (i = 0; i < n; i++)
            {
                cout += "  " + c[i].ToString(CultureInfo.InvariantCulture).PadLeft(4);
            }

            Console.WriteLine(cout);
            Console.WriteLine("");

            for (layer = 0; layer <= 2; layer++)
            {
                Console.WriteLine("");
                Console.WriteLine("  Layer " + layer + "");
                Console.WriteLine("");

                c[n] = layer;
                more = false;
                i = 0;

                for (;;)
                {
                    Geometry.simplex_lattice_layer_point_next(n, c, ref v, ref more);
                    if (!more)
                    {
                        Console.WriteLine("  No more.");
                        break;
                    }

                    i += 1;
                    cout = "  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4);
                    for (j = 0; j < n; j++)
                    {
                        cout += "  " + v[j].ToString(CultureInfo.InvariantCulture).PadLeft(4);
                    }

                    Console.WriteLine(cout);
                }
            }
        }

    }

    public static void test1789()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST1789 tests SIMPLEX_LATTICE_POINT_NEXT.
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
        int TEST_NUM = 4;

        int[] c;
        int i;
        int j;
        bool more;
        int n;
        int[] n_test = {1, 2, 3, 4};
        int test;
        int[] v;
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine("TEST1789");
        Console.WriteLine("  SIMPLEX_LATTICE_POINT_NEXT returns the next lattice");
        Console.WriteLine("  point in an N-dimensional simplex defined by:");
        Console.WriteLine("");
        Console.WriteLine("    0 <= X(1)/C(1) + X(2)/C(2) + ... + X(N)/C(N) <= C(N+1).");

        for (test = 0; test < TEST_NUM; test++)
        {
            n = n_test[test];

            c = new int[n + 1];
            v = new int[n];

            for (i = 0; i < n + 1; i++)
            {
                c[i] = n + 1 - i;
            }

            for (i = 0; i < n; i++)
            {
                v[i] = 0;
            }

            more = false;

            Console.WriteLine("");
            Console.WriteLine("  N = " + n + "");
            cout = "  C = ";
            for (i = 0; i < n + 1; i++)
            {
                cout += "  " + c[i].ToString(CultureInfo.InvariantCulture).PadLeft(4);
            }

            Console.WriteLine(cout);
            Console.WriteLine("");

            i = 0;

            for (;;)
            {
                Geometry.simplex_lattice_point_next(n, c, v, ref more);
                if (!more)
                {
                    Console.WriteLine("  No more.");
                    break;
                }

                i += 1;
                cout = "  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4);
                for (j = 0; j < n; j++)
                {
                    cout += "  " + v[j].ToString(CultureInfo.InvariantCulture).PadLeft(4);
                }

                Console.WriteLine(cout);
            }
        }

    }
        
    public static void test1805 ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST1805 tests SIMPLEX_VOLUME_ND and TETRAHEDRON_VOLUME_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;

        double[] tetra = {
            0.000000,  0.942809, -0.333333,
            -0.816496, -0.816496, -0.333333,
            0.816496, -0.816496, -0.333333,
            0.000000,  0.000000,  1.000000 };
        double volume;

        Console.WriteLine("");
        Console.WriteLine("TEST1805");
        Console.WriteLine("  For an N-dimensional simplex,");
        Console.WriteLine("  SIMPLEX_VOLUME_ND computes the volume.");
        Console.WriteLine("  Here, we check the routine by comparing it");
        Console.WriteLine("  with TETRAHEDRON_VOLUME_3D.");

        typeMethods.r8mat_transpose_print ( DIM_NUM, DIM_NUM+1, tetra, "  Simplex vertices:" );

        volume = Burkardt.TetrahedronNS.Geometry.tetrahedron_volume_3d ( tetra );

        Console.WriteLine("");
        Console.WriteLine("  Volume computed by TETRAHEDRON_VOLUME_3D:");
        Console.WriteLine("  " + volume + "");

        volume = Geometry.simplex_volume_nd ( DIM_NUM, tetra );

        Console.WriteLine("");
        Console.WriteLine("  Volume computed by SIMPLEX_VOLUME_ND:");
        Console.WriteLine("  " + volume + "");

    }

        

}