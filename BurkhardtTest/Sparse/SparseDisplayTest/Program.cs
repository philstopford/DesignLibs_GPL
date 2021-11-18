using System;
using Burkardt.MatrixNS;
using Burkardt.SparsityNS;

namespace SparseDisplayTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SPARSE_DISPLAY_TEST.
        //
        //  Discussion:
        //
        //    SPARSE_DISPLAY_TEST tests the SPARSE_DISPLAY library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("SPARSE_DISPLAY_TEST");
        Console.WriteLine("  Test the SPARSE_DISPLAY library.");

        test01();
        test02();
        test03();

        Console.WriteLine("");
        Console.WriteLine("SPARSE_DISPLAY_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests SPY_GE for a general storage matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        string header = "wathen_ge";

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  SPY_GE generates a sparsity plot for a matrix stored");
        Console.WriteLine("  in general (GE) format.");

        int nx = 5;
        int ny = 5;
        int n = WathenMatrix.wathen_order(nx, ny);
        int seed = 123456789;
        double[] a = WathenMatrix.wathen_ge(nx, ny, n, ref seed);

        Sparsity.spy_ge(n, n, a, header);
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests SPY_FILE in cases where the indices are in a file.
        //
        //  Discussion:
        //
        //    The files used in this example actually use negative column indices
        //    because they were output by DEAL.II and intended to be passed directly
        //    into GNUPLOT without any clever commands.
        //
        //    So my own "SPY_FILE" is currently set up to deal exactly with such
        //    files, and hence, given sensible data will actually show a spy plot
        //    that is transposed - still readable and all, but wrong way round.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  SPY_FILE generates a sparsity plot for a matrix for");
        Console.WriteLine("  which there exists a file containing all the pairs");
        Console.WriteLine("  (I,J) corresponding to nonzero entries.");

        Sparsity.spy_file("before", "before_data.txt");
        Sparsity.spy_file("after", "after_data.txt");
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests SPY_GE for a general storage matrix.
        //
        //  Discussion:
        //
        //    It is not clear to me whether the plot being created is correctly
        //    oriented.  I might be seeing the transpose of the matrix.
        //    One way to check is to set up a moderate sized, highly asymmetric matrix.
        //    In this case, I will create a certain upper triangular matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        string header = "20x30";

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  SPY_GE generates a sparsity plot for a matrix stored");
        Console.WriteLine("  in general (GE) format.");
        Console.WriteLine("  Just to orient ourselves, generate an upper triangular matrix.");

        int m = 20;
        int n = 30;
        double[] a = new double[m * n];

        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (j < i || (i - j) % (i + 1) != 0)
                {
                    a[i + j * m] = 0.0;
                }
                else
                {
                    a[i + j * m] = 1.0;
                }
            }
        }

        Sparsity.spy_ge(m, n, a, header);
    }
}