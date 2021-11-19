namespace Burkardt.SolveNS;

public static class Pentadiagonal
{
    public static double[] penta ( int n, double[] a1, double[] a2, double[] a3, double[] a4, 
            double[] a5, double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PENTA solves a pentadiagonal system of linear equations.
        //
        //  Discussion:
        //
        //    The matrix A is pentadiagonal.  It is entirely zero, except for
        //    the main diagaonal, and the two immediate sub- and super-diagonals.
        //
        //    The entries of Row I are stored as:
        //
        //      A(I,I-2) -> A1(I)
        //      A(I,I-1) -> A2(I)
        //      A(I,I)   -> A3(I)
        //      A(I,I+1) -> A4(I)
        //      A(I,I-2) -> A5(I)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Cheney, Kincaid,
        //    Numerical Mathematics and Computing,
        //    1985, pages 233-236.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A1[N], A2[N], A3[N], A4[N], A5[N], the nonzero
        //    elements of the matrix.  Note that the data in A2, A3 and A4
        //    is overwritten by this routine during the solution process.
        //
        //    Input, double B[N], the right hand side of the linear system.
        //
        //    Output, double PENTA[N], the solution of the linear system.
        //
    {
        int i;
        double xmult;

        double[] x = new double[n];

        for ( i = 1; i < n - 1; i++ )
        {
            xmult = a2[i] / a3[i-1];
            a3[i] -= xmult * a4[i-1];
            a4[i] -= xmult * a5[i-1];
            b[i] -= xmult * b[i-1];
            xmult = a1[i+1] / a3[i-1];
            a2[i+1] -= xmult * a4[i-1];
            a3[i+1] -= xmult * a5[i-1];
            b[i+1] -= xmult * b[i-1];
        }

        xmult = a2[n-1] / a3[n-2];
        a3[n-1] -= xmult * a4[n-2];
        x[n-1] = ( b[n-1] - xmult * b[n-2] ) / a3[n-1];
        x[n-2] = ( b[n-2] - a4[n-2] * x[n-1] ) / a3[n-2];
        for ( i = n - 3; 0 <= i; i-- )
        {
            x[i] = ( b[i] - a4[i] * x[i+1] - a5[i] * x[i+2] ) / a3[i];
        }

        return x;
    }
}