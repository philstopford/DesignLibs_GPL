namespace Burkardt.Values;

public static class LinearSystem
{

    public static void linear_system_values(ref int n_data, ref int nrow, ref int ncol, ref int nsys,
            ref double[] a, ref double[] x, ref double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINEAR_SYSTEM_VALUES returns some linear systems.
        //
        //  Discussion:
        //
        //    Each call to this routine returns scalars NROW, NCOL and NSYS,
        //    which give the dimensions of the linear system
        //
        //      A(NROW,NCOL) * X(NCOL,NSYS) = B(NROW,NSYS)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 May 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, ref int NROW, int NCOL, the number of rows and columns of A.
        //
        //    Output, ref int NSYS, the number of systems.
        //
        //    Output, double **A[NROW*NCOL], the matrix.
        //
        //    Output, double **X[NCOL*NSYS], the solutions of the linear system.
        //
        //    Output, double **B[NROW*NSYS], the right hand sides.
        //
    {
        const int N_MAX = 4;

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        n_data += 1;

        if (N_MAX < n_data)
        {
            n_data = 0;
            nrow = 0;
            ncol = 0;
            nsys = 0;
            a = null;
            b = null;
            x = null;
        }
        else
        {
            switch (n_data)
            {
                case 1:
                    nrow = 3;
                    ncol = 3;
                    nsys = 2;

                    a = new double[nrow * ncol];
                    x = new double[ncol * nsys];
                    b = new double[nrow * nsys];

                    a[0] = 1.0;
                    a[1] = 0.0;
                    a[2] = 0.0;
                    a[0 + 1 * nrow] = 0.0;
                    a[1 + 1 * nrow] = 2.0;
                    a[2 + 1 * nrow] = 0.0;
                    a[0 + 2 * nrow] = 0.0;
                    a[1 + 2 * nrow] = 0.0;
                    a[2 + 2 * nrow] = 3.0;

                    x[0] = 1.0;
                    x[1] = 0.0;
                    x[2] = 0.0;
                    x[0 + 1 * ncol] = 1.0;
                    x[1 + 1 * ncol] = 1.0;
                    x[2 + 1 * ncol] = 1.0;

                    b[0] = 1.0;
                    b[1] = 0.0;
                    b[2] = 0.0;
                    b[0 + 1 * nrow] = 1.0;
                    b[1 + 1 * nrow] = 2.0;
                    b[2 + 1 * nrow] = 3.0;
                    break;
                case 2:
                    nrow = 3;
                    ncol = 3;
                    nsys = 2;

                    a = new double[nrow * ncol];
                    x = new double[ncol * nsys];
                    b = new double[nrow * nsys];

                    a[0] = 1.0;
                    a[1] = 2.0;
                    a[2] = 3.0;
                    a[0 + 1 * nrow] = 2.0;
                    a[1 + 1 * nrow] = 2.0;
                    a[2 + 1 * nrow] = 3.0;
                    a[0 + 2 * nrow] = 3.0;
                    a[1 + 2 * nrow] = 3.0;
                    a[2 + 2 * nrow] = 3.0;

                    x[0] = 1.0;
                    x[1] = 1.0;
                    x[2] = 1.0;
                    x[0 + 1 * ncol] = 1.0;
                    x[1 + 1 * ncol] = 2.0;
                    x[2 + 1 * ncol] = 3.0;

                    b[0] = 6.0;
                    b[1] = 7.0;
                    b[2] = 9.0;
                    b[0 + 1 * nrow] = 14.0;
                    b[1 + 1 * nrow] = 15.0;
                    b[2 + 1 * nrow] = 18.0;
                    break;
                case 3:
                    nrow = 5;
                    ncol = 5;
                    nsys = 2;

                    a = new double[nrow * ncol];
                    x = new double[ncol * nsys];
                    b = new double[nrow * nsys];

                    a[0] = 1.0;
                    a[1] = 2.0;
                    a[2] = 3.0;
                    a[3] = 4.0;
                    a[4] = 5.0;
                    a[0 + 1 * nrow] = 2.0;
                    a[1 + 1 * nrow] = 3.0;
                    a[2 + 1 * nrow] = 4.0;
                    a[3 + 1 * nrow] = 5.0;
                    a[4 + 1 * nrow] = 1.0;
                    a[0 + 2 * nrow] = 3.0;
                    a[1 + 2 * nrow] = 4.0;
                    a[2 + 2 * nrow] = 5.0;
                    a[3 + 2 * nrow] = 1.0;
                    a[4 + 2 * nrow] = 2.0;
                    a[0 + 3 * nrow] = 4.0;
                    a[1 + 3 * nrow] = 5.0;
                    a[2 + 3 * nrow] = 1.0;
                    a[3 + 3 * nrow] = 2.0;
                    a[4 + 3 * nrow] = 3.0;
                    a[0 + 4 * nrow] = 5.0;
                    a[1 + 4 * nrow] = 1.0;
                    a[2 + 4 * nrow] = 2.0;
                    a[3 + 4 * nrow] = 3.0;
                    a[4 + 4 * nrow] = 4.0;

                    x[0] = 0.066667;
                    x[1] = 0.066667;
                    x[2] = 0.066667;
                    x[3] = 0.066667;
                    x[4] = 0.066667;
                    x[0 + 1 * ncol] = 1.0;
                    x[1 + 1 * ncol] = 0.0;
                    x[2 + 1 * ncol] = 0.0;
                    x[3 + 1 * ncol] = 0.0;
                    x[4 + 1 * ncol] = 0.0;

                    b[0] = 1.0;
                    b[1] = 1.0;
                    b[2] = 1.0;
                    b[3] = 1.0;
                    b[4] = 1.0;
                    b[0 + 1 * nrow] = 1.0;
                    b[1 + 1 * nrow] = 2.0;
                    b[2 + 1 * nrow] = 3.0;
                    b[3 + 1 * nrow] = 4.0;
                    b[4 + 1 * nrow] = 5.0;
                    break;
                case 4:
                    nrow = 5;
                    ncol = 5;
                    nsys = 2;

                    a = new double[nrow * ncol];
                    x = new double[ncol * nsys];
                    b = new double[nrow * nsys];

                    a[0] = 1.4;
                    a[1] = 1.6;
                    a[2] = 3.8;
                    a[3] = 4.6;
                    a[4] = 2.6;
                    a[0 + 1 * nrow] = 2.1;
                    a[1 + 1 * nrow] = 1.5;
                    a[2 + 1 * nrow] = 8.0;
                    a[3 + 1 * nrow] = 8.2;
                    a[4 + 1 * nrow] = 2.9;
                    a[0 + 2 * nrow] = 2.1;
                    a[1 + 2 * nrow] = 1.1;
                    a[2 + 2 * nrow] = 9.6;
                    a[3 + 2 * nrow] = 8.4;
                    a[4 + 2 * nrow] = 0.1;
                    a[0 + 3 * nrow] = 7.4;
                    a[1 + 3 * nrow] = 0.7;
                    a[2 + 3 * nrow] = 5.4;
                    a[3 + 3 * nrow] = 0.4;
                    a[4 + 3 * nrow] = 9.6;
                    a[0 + 4 * nrow] = 9.6;
                    a[1 + 4 * nrow] = 5.0;
                    a[2 + 4 * nrow] = 8.8;
                    a[3 + 4 * nrow] = 8.0;
                    a[4 + 4 * nrow] = 7.7;

                    x[0] = -5.313077;
                    x[1] = 5.735670;
                    x[2] = -2.507606;
                    x[3] = -1.058741;
                    x[4] = 0.999381;
                    x[0 + 1 * ncol] = 31.601006;
                    x[1 + 1 * ncol] = -28.594793;
                    x[2 + 1 * ncol] = 13.389395;
                    x[3 + 1 * ncol] = 2.780322;
                    x[4 + 1 * ncol] = -3.008797;

                    b[0] = 1.1;
                    b[1] = 1.6;
                    b[2] = 4.7;
                    b[3] = 9.1;
                    b[4] = 0.1;
                    b[0 + 1 * nrow] = 4.0;
                    b[1 + 1 * nrow] = 9.3;
                    b[2 + 1 * nrow] = 8.4;
                    b[3 + 1 * nrow] = 0.4;
                    b[4 + 1 * nrow] = 4.1;
                    break;
            }
        }
    }

}