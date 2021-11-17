namespace Burkardt.Types;

public static partial class typeMethods
{

    public static void r8mat_scale(int m, int n, double s, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_SCALE multiplies an R8MAT by a scalar.
        //
        //  Discussion:
        //
        //    An R8MAT is an array of R8 values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 December 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double S, the scale factor.
        //
        //    Input/output, double A[M*N], the matrix to be scaled.
        //
    {
        int j;

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                a[i + j * m] *= s;
            }
        }
    }

    public static double[] r8mat_scale_01(int m, int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_SCALE_01 shifts and scales an R8MAT so columns have min 0 and max 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 October 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double X[M*N], the array to be rescaled.
        //
        //    Output, double R8MAT_SCALE_01[M*N], the rescaled array.
        //
    {
        int j;

        double[] xmax = r8mat_max_columns(m, n, x);
        double[] xmin = r8mat_min_columns(m, n, x);

        double[] xs = new double[m * n];

        for (j = 0; j < n; j++)
        {
            int i;
            switch (xmax[j] - xmin[j])
            {
                case > 0:
                {
                    for (i = 0; i < m; i++)
                    {
                        xs[i + j * m] = (x[i + j * m] - xmin[j]) / (xmax[j] - xmin[j]);
                    }

                    break;
                }
                default:
                {
                    for (i = 0; i < m; i++)
                    {
                        xs[i + j * m] = 0.5;
                    }

                    break;
                }
            }
        }

        return xs;
    }

    public static double[] r8mat_scale_ab(int m, int n, double[] x, double a, double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_SCALE_AB shifts and scales an R8MAT so columns have min A and max B.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 October 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double X[M*N], the array to be rescaled.
        //
        //    Input, double A, B, the new scale limits.
        //
        //    Output, double R8MAT_SCALE_AB[M*N], the rescaled array.
        //
    {
        int j;

        double[] xmax = r8mat_max_columns(m, n, x);
        double[] xmin = r8mat_min_columns(m, n, x);

        double[] xs = new double[m * n];

        for (j = 0; j < n; j++)
        {
            int i;
            switch (xmax[j] - xmin[j])
            {
                case > 0:
                {
                    for (i = 0; i < m; i++)
                    {
                        xs[i + j * m] = a + (b - a) * (x[i + j * m] - xmin[j]) / (xmax[j] - xmin[j]);
                    }

                    break;
                }
                default:
                {
                    for (i = 0; i < m; i++)
                    {
                        xs[i + j * m] = (a + b) / 2.0;
                    }

                    break;
                }
            }
        }

        return xs;
    }
        
}