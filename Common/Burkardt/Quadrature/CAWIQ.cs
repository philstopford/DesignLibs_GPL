using System;
using Burkardt.MatrixNS;
using Burkardt.Types;

namespace Burkardt.Quadrature;

public static class CAWIQ
{
    public static double[] cawiq(int nt, double[] t, int[] mlt, int nwts, ref int[] ndx, int key,
            int nst, ref double[] aj, ref double[] bj, ref int jdf, double zemu )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CAWIQ computes quadrature weights for a given set of knots.
        //
        //  Discussion:
        //
        //    This routine is given a set of distinct knots, T, their multiplicities MLT, 
        //    the Jacobi matrix associated with the polynomials orthogonal with respect 
        //    to the weight function W(X), and the zero-th moment of W(X).
        //
        //    It computes the weights of the quadrature formula
        //
        //      sum ( 1 <= J <= NT ) sum ( 0 <= I <= MLT(J) - 1 ) wts(j) d^i/dx^i f(t(j))
        //
        //    which is to approximate
        //
        //      integral ( a < x < b ) f(t) w(t) dt
        //
        //    The routine makes various checks, as indicated below, sets up
        //    various vectors and, if necessary, calls for the diagonalization
        //    of the Jacobi matrix that is associated with the polynomials
        //    orthogonal with respect to W(X) on the interval A, B. 
        //
        //    Then for each knot, the weights of which are required, it calls the 
        //    routine CWIQD which to compute the weights.
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
        //    Input, double T[NT], the knots.
        //
        //    Input, int MLT[NT], the multiplicity of the knots.
        //
        //    Input, int NWTS, the number of weights.
        //
        //    Input/output, int NDX[NT], associates with each distinct 
        //    knot T(J), an integer NDX(J) which is such that the weight to the I-th 
        //    derivative value of F at the J-th knot, is stored in
        //      WTS(abs(NDX(J))+I) for J = 1,2,...,NT, and I = 0,1,2,...,MLT(J)-1.
        //    The sign of NDX includes the following information:
        //    > 0, weights are wanted for this knot
        //    < 0, weights not wanted for this knot but it is included in the quadrature
        //    = 0. means ignore this knot completely.
        //
        //    Input, int KEY, indicates structure of WTS and NDX.
        //    KEY is an integer with absolute value between 1 and 4.
        //    The sign of KEY choosed the form of WTS:
        //    0 < KEY, WTS in standard form.
        //    0 > KEY, J]WTS(J) required.
        //    The absolute value has the following effect:
        //    1, set up pointers in NDX for all knots in T array (routine CAWIQ does 
        //    this).  the contents of NDX are not tested on input and weights are 
        //    packed sequentially in WTS as indicated above.
        //    2, set up pointers only for knots which have nonzero NDX on input.  All 
        //    knots which have a non-zero flag are allocated space in WTS.
        //    3, set up pointers only for knots which have NDX > 0 on input.  Space in 
        //    WTS allocated only for knots with NDX > 0.
        //    4, NDX assumed to be preset as pointer array on input.
        //
        //    Input, int NST, the dimension of the Jacobi matrix.  
        //    NST should be between (N+1)/2 and N.  The usual choice will be (N+1)/2.
        //
        //    Input/output, double AJ[NST], BJ[NST].
        //    If JDF = 0 then AJ contains the  diagonal of the Jacobi matrix and
        //    BJ(1:NST-1) contains the subdiagonal.
        //    If JDF = 1, AJ contains the eigenvalues of the Jacobi matrix and
        //    BJ contains the squares of the elements of the first row of U, the
        //    orthogonal matrix which diagonalized the Jacobi matrix as U*D*U'.
        //
        //    Input/output, int *JDF, indicates whether the Jacobi
        //    matrix needs to be diagonalized.
        //    0, diagonalization required;
        //    1, diagonalization not required.
        //
        //    Input, double ZEMU, the zero-th moment of the weight 
        //    function W(X).
        //
        //    Output, double CAWIQ[NWTS], the weights.
        //
    {
        int i;
        int j;
        int k;
        int l;
        int n = 0;
        double tmp;

        double prec = typeMethods.r8_epsilon();

        switch (nt)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("CAWIQ - Fatal error!");
                Console.WriteLine("  NT < 1.");
                return null;
            //
            //  Check for indistinct knots.
            //
            case > 1:
            {
                k = nt - 1;
                for (i = 1; i <= k; i++)
                {
                    tmp = t[i - 1];
                    l = i + 1;
                    for (j = l; j <= nt; j++)
                    {
                        if (!(Math.Abs(tmp - t[j - 1]) <= prec))
                        {
                            continue;
                        }

                        Console.WriteLine("");
                        Console.WriteLine("CAWIQ - Fatal error!");
                        Console.WriteLine("  Knots too close.");
                        return null;
                    }
                }

                break;
            }
        }

        //
        //  Check multiplicities,
        //  Set up various useful parameters and
        //  set up or check pointers to WTS array.
        //
        l = Math.Abs(key);

        switch (l)
        {
            case < 1:
            case > 4:
                Console.WriteLine("");
                Console.WriteLine("CAWIQ - Fatal error!");
                Console.WriteLine("  Magnitude of KEY not between 1 and 4.");
                return null;
        }

        k = 1;

        switch (l)
        {
            case 1:
            {
                for (i = 1; i <= nt; i++)
                {
                    ndx[i - 1] = k;
                    switch (mlt[i - 1])
                    {
                        case < 1:
                            Console.WriteLine("");
                            Console.WriteLine("CAWIQ - Fatal error!");
                            Console.WriteLine("  MLT(I) < 1.");
                            return null;
                        default:
                            k += mlt[i - 1];
                            break;
                    }
                }

                n = k - 1;
                break;
            }
            case 2:
            case 3:
            {
                n = 0;

                for (i = 1; i <= nt; i++)
                {
                    switch (ndx[i - 1])
                    {
                        case 0:
                            continue;
                    }

                    switch (mlt[i - 1])
                    {
                        case < 1:
                            Console.WriteLine("");
                            Console.WriteLine("CAWIQ - Fatal error!");
                            Console.WriteLine("  MLT(I) < 1.");
                            return null;
                    }

                    n += mlt[i - 1];

                    switch (ndx[i - 1])
                    {
                        case < 0 when l == 3:
                            continue;
                        default:
                            ndx[i - 1] = Math.Abs(k) * typeMethods.i4_sign(ndx[i - 1]);
                            k += mlt[i - 1];
                            break;
                    }
                }

                if (nwts + 1 < k)
                {
                    Console.WriteLine("");
                    Console.WriteLine("CAWIQ - Fatal error!");
                    Console.WriteLine("  NWTS + 1 < K.");
                    return null;
                }

                break;
            }
            case 4:
            {
                for (i = 1; i <= nt; i++)
                {
                    int ip = Math.Abs(ndx[i - 1]);

                    switch (ip)
                    {
                        case 0:
                            continue;
                    }

                    if (nwts < ip + mlt[i - 1])
                    {
                        Console.WriteLine("");
                        Console.WriteLine("CAWIQ - Fatal error!");
                        Console.WriteLine("  NWTS < IPM.");
                        return null;
                    }

                    if (i == nt)
                    {
                        break;
                    }

                    l = i + 1;
                    for (j = l; j <= nt; j++)
                    {
                        int jp = Math.Abs(ndx[j - 1]);
                        if (jp == 0)
                        {
                            continue;
                        }

                        if (jp <= ip + mlt[i - 1] && ip <= jp + mlt[j - 1])
                        {
                            break;
                        }
                    }
                }

                break;
            }
        }

        //
        //  Test some parameters.
        //
        if (nst < (n + 1) / 2)
        {
            Console.WriteLine("");
            Console.WriteLine("CAWIQ - Fatal error!");
            Console.WriteLine("  NST < ( N + 1 ) / 2.");
            return null;
        }

        switch (zemu)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("CAWIQ - Fatal error!");
                Console.WriteLine("  ZEMU <= 0.");
                return null;
        }

        double[] wts = new double[nwts];
        switch (n)
        {
            //
            //  Treat a quadrature formula with 1 simple knot first.
            //
            case <= 1:
            {
                for (i = 0; i < nt; i++)
                {
                    switch (ndx[i])
                    {
                        case > 0:
                            wts[Math.Abs(ndx[i]) - 1] = zemu;
                            return wts;
                    }
                }

                break;
            }
        }

        switch (jdf)
        {
            //
            //  Carry out diagonalization if not already done.
            //
            case 0:
            {
                //
                //  Set unit vector in work field to get back first row of Q.
                //
                double[] z = new double[nst];

                for (i = 0; i < nst; i++)
                {
                    z[i] = 0.0;
                }

                z[0] = 1.0;
                //
                //  Diagonalize the Jacobi matrix.
                //
                IMTQLX.imtqlx(nst, ref aj, ref bj, ref z);
                //
                //  Signal Jacobi matrix now diagonalized successfully.
                //
                jdf = 1;
                //
                //  Save squares of first row of U in subdiagonal array.
                //
                for (i = 0; i < nst; i++)
                {
                    bj[i] = z[i] * z[i];
                }

                break;
            }
        }

        //
        //  Find all the weights for each knot flagged.
        //
        for (i = 1; i <= nt; i++)
        {
            switch (ndx[i - 1])
            {
                case <= 0:
                    continue;
            }

            int m = mlt[i - 1];
            int mnm = Math.Max(n - m, 1);
            l = Math.Min(m, n - m + 1);
            //
            //  Set up K-hat matrix for CWIQD with knots according to their multiplicities.
            //
            double[] xk = new double[mnm];

            k = 1;
            for (j = 1; j <= nt; j++)
            {
                if (ndx[j - 1] == 0 || j == i)
                {
                    continue;
                }

                int jj;
                for (jj = 1; jj <= mlt[j - 1]; jj++)
                {
                    xk[k - 1] = t[j - 1];
                    k += 1;
                }
            }

            //
            //  Set up the right principal vector.
            //
            double[] r = new double[l];

            r[0] = 1.0 / zemu;
            for (j = 1; j < l; j++)
            {
                r[j] = 0.0;
            }

            //
            //  Pick up pointer for the location of the weights to be output.
            //
            k = ndx[i - 1];
            //
            //  Find all the weights for this knot.
            //
            double[] wtmp = CWIQD.cwiqd(m, mnm, l, t[i - 1], xk, nst, aj, bj, r);

            for (j = 0; j < m; j++)
            {
                wts[k - 1 + j] = wtmp[j];
            }

            switch (key)
            {
                case < 0:
                    continue;
            }

            //
            //  Divide by factorials for weights in standard form.
            //
            tmp = 1.0;
            for (j = 1; j < m - 1; j++)
            {
                double p = j;
                tmp *= p;
                wts[k - 1 + j] /= tmp;
            }
        }

        return wts;
    }
}