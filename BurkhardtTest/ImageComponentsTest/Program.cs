using System;
using Burkardt.ImageNS;
using Burkardt.Table;
using Burkardt.Types;

namespace ImageComponentsTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for IMAGE_COMPONENTS_TEST.
        //
        //  Discussion:
        //
        //    IMAGE_COMPONENTS_TEST tests the IMAGE_COMPONENTS library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 March 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("IMAGE_COMPONENTS_TEST");
        Console.WriteLine("  Test the IMAGE_COMPONENTS library.");

        test01();
        test02();
        test03();

        Console.WriteLine("");
        Console.WriteLine("IMAGE_COMPONENTS_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests I4VEC_COMPONENTS on a simple case.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 March 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 28;

        int[] a =
        {
            0, 0, 1, 2, 4, 0, 0, 4, 0, 0,
            0, 8, 9, 9, 1, 2, 3, 0, 0, 5,
            0, 1, 6, 0, 0, 0, 4, 0
        };
        int[] c = new int[N];
        int j;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  I4VEC_COMPONENTS finds and labels connected");
        Console.WriteLine("  components in a 1D integer vector.");

        Console.WriteLine("");
        Console.WriteLine("  A:");
        Console.WriteLine("");
        string cout = "    ";
        for (j = 0; j < N; j++)
        {
            cout += a[j];
        }

        Console.WriteLine(cout);

        int component_num = ImageComponents.i4vec_components(N, a, ref c);

        Console.WriteLine("");
        Console.WriteLine("  Number of components = " + component_num + "");
        Console.WriteLine("");
        Console.WriteLine("  C:");
        Console.WriteLine("");
        cout = "    ";
        for (j = 0; j < N; j++)
        {
            cout += c[j];
        }

        Console.WriteLine(cout);
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests I4MAT_COMPONENTS on a simple case.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 March 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int M = 9;
        const int N = 17;

        int[] a =
        {
            0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 1, 0, 0, 0,
            0, 1, 1, 0, 1, 1, 1, 0, 0,
            0, 1, 1, 1, 1, 1, 1, 0, 0,
            0, 0, 1, 1, 1, 0, 0, 0, 0,
            0, 0, 1, 1, 1, 0, 0, 0, 0,
            0, 1, 1, 1, 0, 1, 0, 1, 0,
            0, 1, 1, 0, 0, 1, 0, 1, 0,
            0, 0, 1, 0, 0, 0, 0, 1, 0,
            0, 0, 0, 0, 1, 0, 1, 1, 0,
            0, 1, 0, 1, 1, 0, 1, 0, 0,
            0, 1, 1, 1, 1, 1, 0, 0, 0,
            0, 0, 1, 1, 0, 1, 0, 1, 0,
            0, 0, 1, 1, 0, 1, 0, 1, 0,
            0, 1, 1, 0, 1, 0, 1, 1, 0,
            0, 1, 0, 0, 1, 0, 1, 1, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0
        };
        int i;
        int j;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  I4MAT_COMPONENTS finds and labels connected");
        Console.WriteLine("  components in a 2D integer array.");

        Console.WriteLine("");
        Console.WriteLine("  A:");
        Console.WriteLine("");
        string cout = "";
        for (i = 0; i < M; i++)
        {
            cout = "    ";
            for (j = 0; j < N; j++)
            {
                cout += a[i + j * M];
            }

            Console.WriteLine(cout);
        }

        int[] c = new int[M * N];

        int component_num = ImageComponents.i4mat_components(M, N, a, ref c);

        Console.WriteLine("");
        Console.WriteLine("  Number of components = " + component_num + "");
        Console.WriteLine("");
        Console.WriteLine("  C:");
        Console.WriteLine("");
        cout = "";
        for (i = 0; i < M; i++)
        {
            cout = "    ";
            for (j = 0; j < N; j++)
            {
                cout += c[i + j * M];
            }

            Console.WriteLine(cout);
        }
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests I4BLOCK_COMPONENTS on a simple case.
        //
        //  Discussion:
        //
        //    This calculation is also done by a program called REGION.
        //    The two programs differ in the number of components discovered
        //    because REGION uses the full 3x3 block of pixels, resulting
        //    in 26 potential neighbors, whereas I4BLOCK_COMPONENTS uses only
        //    the north/south, east/west, up/down directions for 8 neighbors.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 March 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int L = 64;
        const int M = 64;
        const int N = 26;

        int[] a = new int[L * M * N];
        int[] c = new int[L * M * N];
        int i;
        int j;
        int j1;
        int k;

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  I4BLOCK_COMPONENTS finds and labels connected");
        Console.WriteLine("  components in a 3D integer block.");

        Console.WriteLine("");
        Console.WriteLine("  A is a 3D block of order " + L
                                                        + " * " + M
                                                        + " * " + N + "");

        for (k = 0; k < N; k++)
        {
            for (j = 0; j < M; j++)
            {
                for (i = 0; i < L; i++)
                {
                    a[i + j * L + k * L * M] = 0;
                }
            }
        }

        //
        //  Retrieve the indices of nonzero data in A by reading a file.
        //
        string filename = "indices.txt";

        TableHeader h = typeMethods.i4mat_header_read(filename);
        int m1 = h.m;
        int n1 = h.n;

        int[] indices = typeMethods.i4mat_data_read(filename, m1, n1);

        for (j1 = 0; j1 < n1; j1++)
        {
            i = indices[0 + j1 * 3] - 1;
            j = indices[1 + j1 * 3] - 1;
            k = indices[2 + j1 * 3] - 1;
            a[i + j * L + k * L * M] = 1;
        }

        Console.WriteLine("");
        Console.WriteLine("  Number of nonzero A values is " + n1 + "");
        //
        //  Determine the components.
        //
        int component_num = ImageComponents.i4block_components(L, M, N, a, ref c);

        int[] s = new int[component_num];

        for (i = 0; i < component_num; i++)
        {
            s[i] = 0;
        }

        Console.WriteLine("");
        Console.WriteLine("  Number of components = " + component_num + "");

        for (k = 0; k < N; k++)
        {
            for (j = 0; j < M; j++)
            {
                for (i = 0; i < L; i++)
                {
                    if (c[i + j * L + k * L * M] != 0)
                    {
                        s[c[i + j * L + k * L * M] - 1] += 1;
                    }
                }
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  Component  Size");
        Console.WriteLine("");
        int s_total = 0;
        for (i = 0; i < component_num; i++)
        {
            Console.WriteLine("  " + (i + 1).ToString().PadLeft(4)
                                   + "  " + s[i].ToString().PadLeft(8) + "");
            s_total += s[i];
        }

        Console.WriteLine("------  --------");
        Console.WriteLine(" Total  " + s_total.ToString().PadLeft(8) + "");
    }
}