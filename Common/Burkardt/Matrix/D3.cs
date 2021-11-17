using System;
using Burkardt.Uniform;

namespace Burkardt.MatrixNS;

public static class D3
{
    public static double[] d3_mxv ( int n, double[] a, double[] x )

//****************************************************************************80
//
//  Purpose:
//
//    D3_MXV multiplies a D3 matrix times a vector.
//
//  Discussion:
//
//    The D3 storage format is used for a tridiagonal matrix.
//    The superdiagonal is stored in entries (1,2:N), the diagonal in
//    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
//    original matrix is "collapsed" vertically into the array.
//
//  Example:
//
//    Here is how a D3 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the linear system.
//
//    Input, double A[3*N], the D3 matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double D3_MXV[N], the product A * x.
//
    {
        double[] b;
        int i;

        b = new double[n];

        for ( i = 0; i < n; i++ )
        {
            b[i] =        a[1+i*3] * x[i];
        }
        for ( i = 0; i < n-1; i++ )
        {
            b[i] += a[0+(i+1)*3] * x[i+1];
        }
        for ( i = 1; i < n; i++ )
        {
            b[i] += a[2+(i-1)*3] * x[i-1];
        }

        return b;
    }

    public static double[] d3_np_fs ( int n, double[] a, double[] b )

//****************************************************************************80
//
//  Purpose:
//
//    D3_NP_FS factors and solves a D3 system.
//
//  Discussion:
//
//    The D3 storage format is used for a tridiagonal matrix.
//    The superdiagonal is stored in entries (1,2:N), the diagonal in
//    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
//    original matrix is "collapsed" vertically into the array.
//
//    This algorithm requires that each diagonal entry be nonzero.
//    It does not use pivoting, and so can fail on systems that
//    are actually nonsingular.
//
//  Example:
//
//    Here is how a D3 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the linear system.
//
//    Input/output, double A[3*N].
//    On input, the nonzero diagonals of the linear system.
//    On output, the data in these vectors has been overwritten
//    by factorization information.
//
//    Input, double B[N], the right hand side.
//
//    Output, double D3_NP_FS[N], the solution of the linear system.
//    This is NULL if there was an error because one of the diagonal
//    entries was zero.
//
    {
        int i;
        double[] x;
        double xmult;
//
//  Check.
//
        for ( i = 0; i < n; i++ )
        {
            switch (a[1+i*3])
            {
                case 0.0:
                    return null;
            }
        }
        x = new double[n];

        for ( i = 0; i < n; i++ )
        {
            x[i] = b[i];
        }

        for ( i = 1; i < n; i++ )
        {
            xmult = a[2+(i-1)*3] / a[1+(i-1)*3];
            a[1+i*3] -= xmult * a[0+i*3];
            x[i] -= xmult * x[i-1];
        }

        x[n-1] /= a[1+(n-1)*3];
        for ( i = n-2; 0 <= i; i-- )
        {
            x[i] = ( x[i] - a[0+(i+1)*3] * x[i+1] ) / a[1+i*3];
        }

        return x;
    }

    public static void d3_print ( int n, double[] a, string title )

//****************************************************************************80
//
//  Purpose:
//
//    D3_PRINT prints a D3 matrix.
//
//  Discussion:
//
//    The D3 storage format is used for a tridiagonal matrix.
//    The superdiagonal is stored in entries (1,2:N), the diagonal in
//    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
//    original matrix is "collapsed" vertically into the array.
//
//  Example:
//
//    Here is how a D3 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[3*N], the D3 matrix.
//
//    Input, string TITLE, a title to print.
//
    {
        Console.WriteLine("");
        Console.WriteLine(title + "");
        Console.WriteLine("");

        d3_print_some ( n, a, 1, 1, n, n );

    }

    public static void d3_print_some ( int n, double[] a, int ilo, int jlo, int ihi, int jhi )

//****************************************************************************80
//
//  Purpose:
//
//    D3_PRINT_SOME prints some of a D3 matrix.
//
//  Discussion:
//
//    The D3 storage format is used for a tridiagonal matrix.
//    The superdiagonal is stored in entries (1,2:N), the diagonal in
//    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
//    original matrix is "collapsed" vertically into the array.
//
//  Example:
//
//    Here is how a D3 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input, double A[3*N], the D3 matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column, to be printed.
//
    {
        int INCX = 5;

        int i;
        int i2hi;
        int i2lo;
        int inc;
        int j;
        int j2;
        int j2hi;
        int j2lo;
//
//  Print the columns of the matrix, in strips of 5.
//
        for ( j2lo = jlo; j2lo <= jhi; j2lo += INCX )
        {
            j2hi = j2lo + INCX - 1;
            j2hi = Math.Min ( j2hi, n );
            j2hi = Math.Min ( j2hi, jhi );

            inc = j2hi + 1 - j2lo;

            Console.WriteLine("");
            string cout = "  Col: ";
            for ( j = j2lo; j <= j2hi; j++ )
            {
                j2 = j + 1 - j2lo;
                cout += j.ToString().PadLeft(7) + "       ";
            }

            Console.WriteLine(cout);
            Console.WriteLine("  Row");
            Console.WriteLine("  ---");
//
//  Determine the range of the rows in this strip.
//
            i2lo = Math.Max ( ilo, 1 );
            i2lo = Math.Max ( i2lo, j2lo - 1 );

            i2hi = Math.Min ( ihi, n );
            i2hi = Math.Min ( i2hi, j2hi + 1 );

            for ( i = i2lo; i <= i2hi; i++ )
            {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
                cout = i.ToString().PadLeft(6) + "  ";

                for ( j2 = 1; j2 <= inc; j2++ )
                {
                    j = j2lo - 1 + j2;

                    if ( 1 < i-j || 1 < j-i )
                    {
                        cout += "              ";
                    }
                    else if ( j == i+1 )
                    {
                        cout += a[0+(j-1)*3].ToString().PadLeft(12) + "  ";
                    }
                    else if ( j == i )
                    {
                        cout += a[1+(j-1)*3].ToString().PadLeft(12) + "  ";
                    }
                    else if ( j == i-1 )
                    {
                        cout += a[2+(j-1)*3].ToString().PadLeft(12) + "  ";
                    }

                }
                Console.WriteLine(cout);
            }
        }
    }

    public static double[] d3_uniform ( int n, ref int seed )

//****************************************************************************80
//
//  Purpose:
//
//    D3_UNIFORM randomizes a D3 matrix.
//
//  Discussion:
//
//    The D3 storage format is used for a tridiagonal matrix.
//    The superdiagonal is stored in entries (1,2:N), the diagonal in
//    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
//    original matrix is "collapsed" vertically into the array.
//
//  Example:
//
//    Here is how a D3 matrix of order 5 would be stored:
//
//       *  A12 A23 A34 A45
//      A11 A22 A33 A44 A55
//      A21 A32 A43 A54  *
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the linear system.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double D3_UNIFORM[3*N], the D3 matrix.
//
    {
        double[] a;
        int i;
        double[] u;
        double[] v;
        double[] w;

        a = new double[3*n];

        u = UniformRNG.r8vec_uniform_new ( n-1, 0.0, 1.0, ref seed );
        v = UniformRNG.r8vec_uniform_new ( n,   0.0, 1.0, ref seed );
        w = UniformRNG.r8vec_uniform_new ( n-1, 0.0, 1.0, ref seed );

        a[0+0*3] = 0.0;
        for ( i = 1; i < n; i++ )
        {
            a[0+i*3] = u[i-1];
        }
        for ( i = 0; i < n; i++ )
        {
            a[1+i*3] = v[i];
        }
        for ( i = 0; i < n-1; i++ )
        {
            a[2+i*3] = w[i];
        }
        a[2+(n-1)*3] = 0.0;

        return a;
    }        
}