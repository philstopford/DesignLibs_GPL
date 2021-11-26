using System;
using Burkardt.Storage;
using Burkardt.Types;
using Burkardt.Uniform;

namespace SparseTripletToCompressedColumnTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for ST_TO_ccs_TEST.
        //
        //  Discussion:
        //
        //    ST_TO_ccs_TEST tests st_to_ccs.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("ST_TO_ccs_TEST");
        Console.WriteLine("  Test st_to_ccs.");

        test01();
        test02();
        test03();
        test04();

        Console.WriteLine("");
        Console.WriteLine("ST_TO_ccs_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests st_to_ccs using a tiny matrix.
        //
        //  Discussion:
        //
        //    This test uses a trivial matrix whose full representation is:
        //
        //          2  3  0  0  0
        //          3  0  4  0  6
        //      A = 0 -1 -3  2  0
        //          0  0  1  0  0
        //          0  4  2  0  1
        //
        //    A (1-based) ST representation, reading in order by rows is:
        //
        //      I  J   A
        //     -- --  --
        //      1  1   2
        //      1  2   3
        //
        //      2  1   3
        //      2  3   4
        //      2  5   6
        //
        //      3  2  -1
        //      3  3  -3
        //      3  4   2
        //
        //      4  3   1
        //
        //      5  2   4
        //      5  3   2
        //      5  5   1
        //
        //    The CCS representation (which goes in order by columns) is
        //
        //      #   I  JC   A
        //     --  --  --  --
        //      1   1   1   2
        //      2   2       3
        //
        //      3   1   3   3
        //      4   3      -1
        //      5   5       4
        //
        //      6   2   6   4
        //      7   3      -3
        //      8   4       1
        //      9   5       2
        //
        //     10   3  10   2
        //
        //     11   2  11   6
        //     12   5       1
        //
        //     13   *  13
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int NST = 12;

        double[] ast =  {
                2.0, 3.0,
                3.0, 4.0, 6.0,
                -1.0, -3.0, 2.0,
                1.0,
                4.0, 2.0, 1.0
            }
            ;
        int[] ist =  {
                1, 1,
                2, 2, 2,
                3, 3, 3,
                4,
                5, 5, 5
            }
            ;
        int[] jst =  {
                1, 2,
                1, 3, 5,
                2, 3, 4,
                3,
                2, 3, 5
            }
            ;
        const int m = 5;
        const int n = 5;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Convert a sparse matrix from ST to CCS format.");
        Console.WriteLine("  ST:  sparse triplet,    I, J,  A.");
        Console.WriteLine("  CCS: compressed column, I, CC, A.");

        int i_min = typeMethods.i4vec_min(NST, ist);
        int i_max = typeMethods.i4vec_max(NST, ist);
        int j_min = typeMethods.i4vec_min(NST, jst);
        int j_max = typeMethods.i4vec_max(NST, jst);

        SparseTriplet.st_header_print(i_min, i_max, j_min, j_max, m, n, NST);
        //
        //  Decrement the 1-based data.
        //
        typeMethods.i4vec_dec(NST, ref ist);
        typeMethods.i4vec_dec(NST, ref jst);
        //
        //  Print the ST matrix.
        //
        SparseTriplet.st_print(m, n, NST, ist, jst, ast, "  The matrix in ST format:");
        //
        //  Get the CCS size.
        //
        int ncc = SparseTriplet.st_to_ccs_size(NST, ist, jst);

        Console.WriteLine("");
        Console.WriteLine("  Number of CCS values = " + ncc + "");
        //
        //  Create the CCS indices.
        //
        int[] icc = new int[ncc];
        int[] ccc = new int[n + 1];

        SparseTriplet.st_to_ccs_index(NST, ist, jst, ncc, n, ref icc, ref ccc);
        //
        //  Create the CCS values.
        //
        double[] acc = SparseTriplet.st_to_ccs_values(NST, ist, jst, ast, ncc, n, icc, ccc);
        //
        //  Print the CCS matrix.
        //
        CompressedColumnStorage.ccs_print(m, n, ncc, icc, ccc, acc, "  CCS Matrix:");
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests st_to_ccs on a matrix stored in a file.
        //
        //  Discussion:
        //
        //    We assume no prior knowledge about the matrix except the filename.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const string filename_st = "west_st.txt";
        int i_max = 0;
        int i_min = 0;
        int j_max = 0;
        int j_min = 0;
        int m = 0;
        int n = 0;
        int nst = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Convert a sparse matrix from ST to CCS format.");
        Console.WriteLine("  ST:  sparse triplet,    I, J,  A.");
        Console.WriteLine("  CCS: compressed column, I, CC, A.");
        Console.WriteLine("  This matrix is read from the file '" + filename_st + "'");
        //
        //  Get the size of the ST matrix.
        //
        SparseTriplet.st_header_read(filename_st, ref i_min, ref i_max, ref j_min, ref j_max, ref m, ref n, ref nst);

        SparseTriplet.st_header_print(i_min, i_max, j_min, j_max, m, n, nst);
        //
        //  Allocate space.
        //
        int[] ist = new int[nst];
        int[] jst = new int[nst];
        double[] ast = new double[nst];
        //
        //  Read the ST matrix.
        //
        SparseTriplet.st_data_read(filename_st, m, n, nst, ref ist, ref jst, ref ast);
        //
        //  Decrement the 1-based data.
        //
        typeMethods.i4vec_dec(nst, ref ist);
        typeMethods.i4vec_dec(nst, ref jst);
        //
        //  Print the ST matrix.
        //
        SparseTriplet.st_print(m, n, nst, ist, jst, ast, "  The matrix in ST format:");
        //
        //  Get the CCS size.
        //
        int ncc = SparseTriplet.st_to_ccs_size(nst, ist, jst);

        Console.WriteLine("");
        Console.WriteLine("  Number of CCS values = " + ncc + "");
        //
        //  Create the CCS indices.
        //
        int[] icc = new int[ncc];
        int[] ccc = new int[n + 1];

        SparseTriplet.st_to_ccs_index(nst, ist, jst, ncc, n, ref icc, ref ccc);
        //
        //  Create the CCS values.
        //
        double[] acc = SparseTriplet.st_to_ccs_values(nst, ist, jst, ast, ncc, n, icc, ccc);
        //
        //  Print the CCS matrix.
        //
        CompressedColumnStorage.ccs_print(m, n, ncc, icc, ccc, acc, "  CCS Matrix:");
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 creates a CCS sparse matrix file from an ST file.
        //
        //  Discussion:
        //
        //    We assume no prior knowledge about the matrix except the filename.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const string filename_acc = "west_acc.txt";
        const string filename_ccc = "west_ccc.txt";
        const string filename_icc = "west_icc.txt";
        const string filename_st = "west_st.txt";
        int i_max = 0;
        int i_min = 0;
        int j_max = 0;
        int j_min = 0;
        int m = 0;
        int n = 0;
        int nst = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  Convert a sparse matrix from ST to CCS format.");
        Console.WriteLine("  ST: sparse triplet,    I, J,  A.");
        Console.WriteLine("  CC: compressed column, I, CC, A.");
        Console.WriteLine("  The ST matrix is read from the file '" + filename_st + "'");
        Console.WriteLine("  and the CCS matrix is written to the files:");
        Console.WriteLine("    '" + filename_icc + "',");
        Console.WriteLine("    '" + filename_ccc + "', and");
        Console.WriteLine("    '" + filename_acc + "'.");
        //
        //  Get the size of the ST matrix.
        //
        SparseTriplet.st_header_read(filename_st, ref i_min, ref i_max, ref j_min, ref j_max, ref m, ref n, ref nst);

        SparseTriplet.st_header_print(i_min, i_max, j_min, j_max, m, n, nst);
        //
        //  Allocate space.
        //
        int[] ist = new int[nst];
        int[] jst = new int[nst];
        double[] ast = new double[nst];
        //
        //  Read the ST matrix.
        //
        SparseTriplet.st_data_read(filename_st, m, n, nst, ref ist, ref jst, ref ast);
        //
        //  Decrement the 1-based data.
        //
        typeMethods.i4vec_dec(nst, ref ist);
        typeMethods.i4vec_dec(nst, ref jst);
        //
        //  Get the CCS size.
        //
        int ncc = SparseTriplet.st_to_ccs_size(nst, ist, jst);

        Console.WriteLine("");
        Console.WriteLine("  Number of CCS values = " + ncc + "");
        //
        //  Create the CCS indices.
        //
        int[] icc = new int[ncc];
        int[] ccc = new int[n + 1];

        SparseTriplet.st_to_ccs_index(nst, ist, jst, ncc, n, ref icc, ref ccc);
        //
        //  Create the CCS values.
        //
        double[] acc = SparseTriplet.st_to_ccs_values(nst, ist, jst, ast, ncc, n, icc, ccc);
        //
        //  Write the CCS matrix.
        //
        typeMethods.i4vec_write(filename_icc, ncc, icc);
        typeMethods.i4vec_write(filename_ccc, n + 1, ccc);
        typeMethods.r8vec_write(filename_acc, ncc, acc);
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 works with a CCS sparse matrix with many repeated index pairs.
        //
        //  Discussion:
        //
        //    To complete this test, I want to compare AST * X and ACC * X.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  Convert a sparse matrix from ST to CCS format.");
        Console.WriteLine("  ST: sparse triplet,    I, J,  A.");
        Console.WriteLine("  CC: compressed column, I, CC, A.");
        Console.WriteLine("  The ST matrix is the Wathen finite element matrix.");
        Console.WriteLine("  It has many repeated index pairs.");
        Console.WriteLine("  To check, compare ACC*X - AST*X for a random X.");
        //
        //  Get the size of the ST matrix.
        //
        const int nx = 3;
        const int ny = 3;
        int nst = SparseTriplet.wathen_st_size(nx, ny);

        Console.WriteLine("");
        Console.WriteLine("  Number of ST values = " + nst + "");
        //
        //  Set the formal matrix size
        //
        const int m = 3 * nx * ny + 2 * nx + 2 * ny + 1;
        //
        //  Set a random vector.
        //
        int seed = 123456789;
        double[] x = UniformRNG.r8vec_uniform_01_new(m, ref seed);
        //
        //  Allocate space.
        //
        int[] ist = new int[nst];
        int[] jst = new int[nst];
        //
        //  Create the ST matrix.
        //
        seed = 123456789;
        double[] ast = SparseTriplet.wathen_st(nx, ny, nst, ref seed, ref ist, ref jst);

        int i_min = typeMethods.i4vec_min(nst, ist);
        int i_max = typeMethods.i4vec_max(nst, ist);
        int j_min = typeMethods.i4vec_min(nst, jst);
        int j_max = typeMethods.i4vec_max(nst, jst);

        SparseTriplet.st_header_print(i_min, i_max, j_min, j_max, m, m, nst);
        //
        //  Compute B1 = AST * X
        //
        double[] b1 = SparseTriplet.st_mv(m, m, nst, ist, jst, ast, x);
        //
        //  Get the CCS size.
        //
        int ncc = SparseTriplet.st_to_ccs_size(nst, ist, jst);

        Console.WriteLine("  Number of CCS values = " + ncc + "");
        //
        //  Create the CCS indices.
        //
        int[] icc = new int[ncc];
        int[] ccc = new int[m + 1];
        SparseTriplet.st_to_ccs_index(nst, ist, jst, ncc, m, ref icc, ref ccc);
        //
        //  Create the CCS values.
        //
        double[] acc = SparseTriplet.st_to_ccs_values(nst, ist, jst, ast, ncc, m, icc, ccc);
        //
        //  Compute B2 = ACC * X.
        //
        double[] b2 = CompressedColumnStorage.ccs_mv(m, m, ncc, icc, ccc, acc, x);
        //
        //  Compare B1 and B2.
        //
        double r = typeMethods.r8vec_diff_norm(m, b1, b2);
        Console.WriteLine("  || ACC*X - AST*X|| = " + r + "");
    }
}