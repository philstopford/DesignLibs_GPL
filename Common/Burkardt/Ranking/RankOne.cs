using System;
using System.Globalization;
using Burkardt.Types;

namespace Burkardt.RankingNS;

public static partial class Ranking
{
    public static void rank_one_print_test(int m, int n, double[] a, double[] u,
            double[] s, double[] v)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RANK_ONE_PRINT_TEST prints the sums of rank one matrices.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in the matrix A.
        //
        //    Input, double A[M*N], the matrix whose singular value
        //    decomposition we are investigating.
        //
        //    Input, double U[M*M], S[M*N], V[N*N], the factors
        //    that form the singular value decomposition of A.
        //
    {
        int i;
        int j;
        int k;
        int r;
        double[] svt;
        string title = "";
        double[] usvt;

        //a_norm = r8mat_norm_fro ( m, n, a );

        Console.WriteLine("");
        Console.WriteLine("RANK_ONE_PRINT_TEST:");
        Console.WriteLine("  Print the sums of R rank one matrices.");

        for (r = 0; r <= Math.Min(m, n); r++)
        {
            svt = new double[r * n];
            for (i = 0; i < r; i++)
            {
                for (j = 0; j < n; j++)
                {
                    svt[i + j * r] = 0.0;
                    for (k = 0; k < r; k++)
                    {
                        svt[i + j * r] += s[i + k * m] * v[j + k * n];
                    }
                }
            }

            usvt = new double[m * n];

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < n; j++)
                {
                    usvt[i + j * m] = 0.0;
                    for (k = 0; k < r; k++)
                    {
                        usvt[i + j * m] += u[i + k * m] * svt[k + j * r];
                    }
                }
            }

            Console.WriteLine(title + "  Rank R = " + r);
            typeMethods.r8mat_print(m, n, usvt, title);


        }

        typeMethods.r8mat_print(m, n, a, "  Original matrix A:");

    }

    public static void rank_one_test(int m, int n, double[] a, double[] u, double[] s,
            double[] v)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RANK_ONE_TEST compares A to the sum of rank one matrices.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in the matrix A.
        //
        //    Input, double A[M*N], the matrix whose singular value
        //    decomposition we are investigating.
        //
        //    Input, double U[M*M], S[M*N], V[N*N], the factors
        //    that form the singular value decomposition of A.
        //
    {
        double a_norm;
        double dif_norm;
        int i;
        int j;
        int k;
        int r;
        double[] svt;
        double[] usvt;

        a_norm = typeMethods.r8mat_norm_fro(m, n, a);

        Console.WriteLine("");
        Console.WriteLine("RANK_ONE_TEST:");
        Console.WriteLine("  Compare A to the sum of R rank one matrices.");
        Console.WriteLine("");
        Console.WriteLine("         R    Absolute      Relative");
        Console.WriteLine("              Error         Error");
        Console.WriteLine("");

        for (r = 0; r <= Math.Min(m, n); r++)
        {
            svt = new double[r * n];
            for (i = 0; i < r; i++)
            {
                for (j = 0; j < n; j++)
                {
                    svt[i + j * r] = 0.0;
                    for (k = 0; k < r; k++)
                    {
                        svt[i + j * r] += s[i + k * m] * v[j + k * n];
                    }
                }
            }

            usvt = new double[m * n];

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < n; j++)
                {
                    usvt[i + j * m] = 0.0;
                    for (k = 0; k < r; k++)
                    {
                        usvt[i + j * m] += u[i + k * m] * svt[k + j * r];
                    }
                }
            }

            dif_norm = typeMethods.r8mat_dif_fro(m, n, a, usvt);

            Console.WriteLine("  " + r.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + dif_norm.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + (dif_norm / a_norm).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

        }
    }
}