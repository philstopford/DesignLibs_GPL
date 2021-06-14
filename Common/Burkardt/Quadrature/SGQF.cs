using System;

namespace Burkardt
{
    public static class SGQF
    {
        public static void sgqf(int nt, double[] aj, ref double[] bj, double zemu, ref double[] t,
        ref double[] wts )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SGQF computes knots and weights of a Gauss Quadrature formula.
        //
        //  Discussion:
        //
        //    This routine computes all the knots and weights of a Gauss quadrature
        //    formula with simple knots from the Jacobi matrix and the zero-th
        //    moment of the weight function, using the Golub-Welsch technique.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 January 2010
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Sylvan Elhay, Jaroslav Kautsky,
        //    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
        //    Interpolatory Quadrature,
        //    ACM Transactions on Mathematical Software,
        //    Volume 13, Number 4, December 1987, pages 399-415.
        //
        //  Parameters:
        //
        //    Input, int NT, the number of knots.
        //
        //    Input, double AJ[NT], the diagonal of the Jacobi matrix.
        //
        //    Input/output, double BJ[NT], the subdiagonal of the Jacobi 
        //    matrix, in entries 1 through NT-1.  On output, BJ has been overwritten.
        //
        //    Input, double ZEMU, the zero-th moment of the weight function.
        //
        //    Output, double T[NT], the knots.
        //
        //    Output, double WTS[NT], the weights.
        //
        {
            int i;
            //
            //  Exit if the zero-th moment is not positive.
            //
            if (zemu <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("SGQF - Fatal error!");
                Console.WriteLine("  ZEMU <= 0.");
                return;
            }

            //
            //  Set up vectors for IMTQLX.
            //
            for (i = 0; i < nt; i++)
            {
                t[i] = aj[i];
            }

            wts[0] = Math.Sqrt(zemu);
            for (i = 1; i < nt; i++)
            {
                wts[i] = 0.0;
            }

            //
            //  Diagonalize the Jacobi matrix.
            //
            IMTQLX.imtqlx(nt, t, bj, ref wts);

            for (i = 0; i < nt; i++)
            {
                wts[i] = wts[i] * wts[i];
            }
        }
    }
}