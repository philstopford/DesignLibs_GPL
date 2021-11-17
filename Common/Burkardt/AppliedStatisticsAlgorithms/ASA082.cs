using System;

namespace Burkardt.AppliedStatistics;

public static partial class Algorithms
{
    public static void detq(double[] a, int n, ref double d, ref int ifault )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    detq computes the determinant of an orthogonal matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 January 2020
        //
        //  Author:
        //
        //    Original FORTRAN77 version by J C Gower.
        //    This C++ version by John Burkardt
        //
        //  Reference:
        //
        //    J C Gower,
        //    Algorithm AS 82:
        //    The determinant of an orthogonal matrix,
        //    Applied Statistics,
        //    Volume 24, Number 1, 1975, page 150-153.
        //
        //  Input:
        //
        //    double A[N*N], the orthogonal matrix stored by rows or columns.
        //
        //    int N, the order of the matrix.
        //
        //  Output:
        //
        //    double *D, the determinant of A.
        //
        //    integer *IFAULT, 
        //    0, no error occurred.
        //    1, an error was detected.
        //
    {
        int i;
        int j;
        int k;
        int p;
        int q;
        int r;
        int s;
        double tol;
        double x;
        double y;

        d = 0.0;
        ifault = 0;
        tol = 0.0001;

        switch (n)
        {
            case <= 0:
                Console.WriteLine("");
                Console.WriteLine("detq - Fatal error!");
                Console.WriteLine("  n <= 0");
                ifault = 1;
                return;
        }

        double[] a2 = new double[n * n];

        for (k = 0; k < n * n; k++)
        {
            a2[k] = a[k];
        }

        d = 1.0;
        r = 0;

        for (k = 2; k <= n + 1; k++)
        {
            q = r;
            x = a2[r];
            y = x switch
            {
                < 0.0 => -1.0,
                _ => +1.0
            };

            d *= y;
            y = -1.0 / (x + y);
            x = Math.Abs(x) - 1.0;

            if (tol < Math.Abs(x))
            {
                switch (x)
                {
                    case > 0.0:
                        Console.WriteLine("");
                        Console.WriteLine("detq - Fatal error!");
                        Console.WriteLine("  x < 0.0");
                        Console.WriteLine("  x = " + x + "");
                        ifault = 1;
                        return;
                }

                if (k == n + 1)
                {
                    Console.WriteLine("");
                    Console.WriteLine("detq - Fatal error!");
                    Console.WriteLine("  k == n + 1");
                    ifault = 1;
                    return;
                }

                for (i = k; i <= n; i++)
                {
                    q += n;
                    x = a2[q] * y;
                    p = r;
                    s = q;
                    for (j = k; j <= n; j++)
                    {
                        p += 1;
                        s += 1;
                        a2[s] += x * a2[p];
                    }
                }
            }

            r = r + n + 1;
        }
    }
}