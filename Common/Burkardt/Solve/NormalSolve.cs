using Burkardt.ChebyshevNS;
using Burkardt.Table;
using Burkardt.Types;

namespace Burkardt.SolveNS;

public static class NormalSolve
{
    public static double[] normal_solve(int m, int n, double[] a, double[] b, ref int flag )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_SOLVE solves a linear system using the normal equations.
        //
        //  Discussion:
        //
        //    Given a presumably rectangular MxN system of equations A*x=b, this routine
        //    sets up the NxN system A'*A*x=A'b.  Assuming N <= M, and that A has full
        //    column rank, the system will be solvable, and the vector x that is returned
        //    will minimize the Euclidean norm of the residual.
        //
        //    One drawback to this approach is that the condition number of the linear
        //    system A'*A is effectively the square of the condition number of A, 
        //    meaning that there is a substantial loss of accuracy.
        //
        //    Thanks to David Doria for pointing out that this procedure was missing
        //    "delete[]" statements at the end, 08 April 2013.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 April 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    David Kahaner, Cleve Moler, Steven Nash,
        //    Numerical Methods and Software,
        //    Prentice Hall, 1989,
        //    ISBN: 0-13-627258-4,
        //    LC: TA345.K34.
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of A.
        //
        //    Input, int N, the number of columns of A.
        //    It must be the case that N <= M.
        //
        //    Input, double A[M*N], the matrix.
        //    The matrix must have full column rank.
        //
        //    Input, double B[M], the right hand side.
        //
        //    Output, int &FLAG,
        //    0, no error was detected.
        //    1, an error occurred.
        //
        //    Output, double NORMAL_SOLVE[N], the least squares solution.
        //
    {
        flag = 0;

        if (m < n)
        {
            flag = 1;
            return null;
        }

        double[] at = typeMethods.r8mat_transpose_new(m, n, a);

        double[] ata = typeMethods.r8mat_mm_new(n, m, n, at, a);

        double[] ata_c = typeMethods.r8mat_cholesky_factor(n, ata, ref flag);

        if (flag != 0)
        {
            return null;
        }

        double[] atb = typeMethods.r8mat_mv_new(n, m, at, b);

        double[] x = typeMethods.r8mat_cholesky_solve(n, ata_c, atb);

        return x;
    }
}