using System;
using Burkardt.Types;

namespace Burkardt.MinpackNS;

public static partial class Minpack
{
    public static void r1mpyq(int m, int n, ref double[] a, int lda, double[] v, double[] w, int aIndex = 0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    r1mpyq() multiplies an M by N matrix A by the Q factor.
        //
        //  Discussion:
        //
        //    Given an M by N matrix A, this function computes a*q where
        //    q is the product of 2*(n - 1) transformations
        //
        //      gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
        //
        //    and gv(i), gw(i) are givens rotations in the (i,n) plane which
        //    eliminate elements in the i-th and n-th planes, respectively.
        //
        //    Q itself is not given, rather the information to recover the
        //    GV and GW rotations is supplied.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 April 2010
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Jorge More, Burt Garbow, Ken Hillstrom.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Jorge More, Burton Garbow, Kenneth Hillstrom,
        //    User Guide for MINPACK-1,
        //    Technical Report ANL-80-74,
        //    Argonne National Laboratory, 1980.
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of A.
        //
        //    Input, int N, the number of columns of A.
        //
        //    Input/output, double A[M*N].  On input, the matrix to be postmultiplied 
        //    by the orthogonal matrix Q described above.  On output, the value of A*Q.
        //
        //    Input, int LDA, a positive value not less than M
        //    which specifies the leading dimension of the array A.
        //
        //    Input, double V[N].  V(I) must contain the information necessary to 
        //    recover the givens rotation GV(I) described above.
        //
        //    Input, double W[N], contains the information necessary to recover the
        //    Givens rotation gw(i) described above.
        //
    {
        double c;
        int i;
        int j;
        double s;
        double temp;
        //
        //  Apply the first set of Givens rotations to A.
        //
        for (j = n - 2; 0 <= j; j--)
        {
            switch (Math.Abs(v[j]))
            {
                case > 1.0:
                    c = 1.0 / v[j];
                    s = Math.Sqrt(1.0 - c * c);
                    break;
                default:
                    s = v[j];
                    c = Math.Sqrt(1.0 - s * s);
                    break;
            }

            for (i = 0; i < m; i++)
            {
                temp = c * a[aIndex + i + j * lda] - s * a[aIndex + i + (n - 1) * lda];
                a[aIndex + i + (n - 1) * lda] = s * a[aIndex + i + j * lda] + c * a[aIndex + i + (n - 1) * lda];
                a[aIndex + i + j * lda] = temp;
            }
        }

        //
        //  Apply the second set of Givens rotations to A.
        //
        for (j = 0; j < n - 1; j++)
        {
            switch (Math.Abs(w[j]))
            {
                case > 1.0:
                    c = 1.0 / w[j];
                    s = Math.Sqrt(1.0 - c * c);
                    break;
                default:
                    s = w[j];
                    c = Math.Sqrt(1.0 - s * s);
                    break;
            }

            for (i = 0; i < m; i++)
            {
                temp = c * a[aIndex + i + j * lda] + s * a[aIndex + i + (n - 1) * lda];
                a[aIndex + i + (n - 1) * lda] = -s * a[aIndex + i + j * lda] + c * a[aIndex + i + (n - 1) * lda];
                a[aIndex + i + j * lda] = temp;
            }
        }
    }

    public static bool r1updt(int m, int n, ref double[] s, int ls, double[] u, ref double[] v,
            double[] w, int sIndex = 0, int uIndex = 0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    r1upd() updates the Q factor after a rank one update of the matrix.
        //
        //  Discussion:
        //
        //    Given an M by N lower trapezoidal matrix S, an M-vector U,
        //    and an N-vector V, the problem is to determine an
        //    orthogonal matrix Q such that
        //
        //      (S + U*V') * Q
        //
        //    is again lower trapezoidal.
        //
        //    This function determines q as the product of 2*(n - 1) transformations
        //
        //      gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
        //
        //    where gv(i), gw(i) are givens rotations in the (i,n) plane
        //    which eliminate elements in the i-th and n-th planes,
        //    respectively. 
        //
        //    Q itself is not accumulated, rather the
        //    information to recover the gv, gw rotations is returned.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 April 2010
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Jorge More, Burt Garbow, Ken Hillstrom.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Jorge More, Burton Garbow, Kenneth Hillstrom,
        //    User Guide for MINPACK-1,
        //    Technical Report ANL-80-74,
        //    Argonne National Laboratory, 1980.
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of S.
        //
        //    Input, int N, the number of columns of S.  N must not exceed M.
        //
        //    Input/output, double S[LS].  On input, the lower trapezoidal matrix S 
        //    stored by columns.  On output, the lower trapezoidal matrix produced 
        //    as described above.
        //
        //    Input, int LS, a positive value not less than (N*(2*M-N+1))/2.
        //
        //    Input, double U[M], contains the vector U.
        //
        //    Input/output, double V[N].  On input, the vector v.  On output,
        //    information necessary to recover the givens rotation gv(i) described above.
        //
        //       w is an output array of length m. w(i) contains information
        //         necessary to recover the givens rotation gw(i) described
        //         above.
        //
        //       sing is a logical output variable. sing is set true if any
        //         of the diagonal elements of the output s are zero. otherwise
        //         sing is set false.
        //
    {
        double cotan;
        double cs;
        double giant;
        int i;
        int j;
        int jj;
        int l;
        int nm1;
        const double p25 = 0.25;
        const double p5 = 0.5;
        double sn;
        bool sing;
        double tan;
        double tau;
        double temp;
        //
        //  Because of the computation of the pointer JJ, this function was
        //  converted from FORTRAN77 to C++ in a conservative way.  All computations
        //  are the same, and only array indexing is adjusted.
        //
        //  GIANT is the largest magnitude.
        //
        giant = typeMethods.r8_huge();
        //
        //  Initialize the diagonal element pointer.
        //
        jj = n * (2 * m - n + 1) / 2 - (m - n);
        //
        //  Move the nontrivial part of the last column of S into W.
        //
        l = jj;
        for (i = n; i <= m; i++)
        {
            w[i - 1] = s[sIndex + (l - 1)];
            l += 1;
        }

        //
        //  Rotate the vector V into a multiple of the N-th unit vector
        //  in such a way that a spike is introduced into W.
        //
        nm1 = n - 1;

        for (j = n - 1; 1 <= j; j--)
        {
            jj -= (m - j + 1);
            w[j - 1] = 0.0;

            if (v[j - 1] != 0.0)
            {
                //
                //  Determine a Givens rotation which eliminates the J-th element of V.
                //
                if (Math.Abs(v[n - 1]) < Math.Abs(v[j - 1]))
                {
                    cotan = v[n - 1] / v[j - 1];
                    sn = p5 / Math.Sqrt(p25 + p25 * cotan * cotan);
                    cs = sn * cotan;
                    tau = (Math.Abs(cs) * giant) switch
                    {
                        > 1.0 => 1.0 / cs,
                        _ => 1.0
                    };
                }
                else
                {
                    tan = v[j - 1] / v[n - 1];
                    cs = p5 / Math.Sqrt(p25 + p25 * tan * tan);
                    sn = cs * tan;
                    tau = sn;
                }

                //
                //  Apply the transformation to V and store the information
                //  necessary to recover the Givens rotation.
                //
                v[n - 1] = sn * v[j - 1] + cs * v[n - 1];
                v[j - 1] = tau;
                //
                //  Apply the transformation to S and extend the spike in W.
                //
                l = jj;
                for (i = j; i <= m; i++)
                {
                    temp = cs * s[sIndex + (l - 1)] - sn * w[i - 1];
                    w[i - 1] = sn * s[sIndex + (l - 1)] + cs * w[i - 1];
                    s[sIndex + (l - 1)] = temp;
                    l += 1;
                }
            }
        }

        //
        //  Add the spike from the rank 1 update to W.
        //
        for (i = 1; i <= m; i++)
        {
            w[i - 1] += v[n - 1] * u[uIndex + (i - 1)];
        }

        //
        //  Eliminate the spike.
        //
        sing = false;

        for (j = 1; j <= nm1; j++)
        {
            //
            //  Determine a Givens rotation which eliminates the
            //  J-th element of the spike.
            //
            if (w[j - 1] != 0.0)
            {

                if (Math.Abs(s[sIndex + (jj - 1)]) < Math.Abs(w[j - 1]))
                {
                    cotan = s[sIndex + (jj - 1)] / w[j - 1];
                    sn = p5 / Math.Sqrt(p25 + p25 * cotan * cotan);
                    cs = sn * cotan;
                    tau = (Math.Abs(cs) * giant) switch
                    {
                        > 1.0 => 1.0 / cs,
                        _ => 1.0
                    };
                }
                else
                {
                    tan = w[j - 1] / s[sIndex + (jj - 1)];
                    cs = p5 / Math.Sqrt(p25 + p25 * tan * tan);
                    sn = cs * tan;
                    tau = sn;
                }

                //
                //  Apply the transformation to s and reduce the spike in w.
                //
                l = jj;

                for (i = j; i <= m; i++)
                {
                    temp = cs * s[sIndex + (l - 1)] + sn * w[i - 1];
                    w[i - 1] = -sn * s[sIndex + (l - 1)] + cs * w[i - 1];
                    s[sIndex + (l - 1)] = temp;
                    l += 1;
                }

                //
                //  Store the information necessary to recover the givens rotation.
                //
                w[j - 1] = tau;
            }

            sing = s[sIndex + (jj - 1)] switch
            {
                //
                //  Test for zero diagonal elements in the output s.
                //
                0.0 => true,
                _ => sing
            };

            jj = jj + (m - j) + 1;
        }

        //
        //  Move W back into the last column of the output S.
        //
        l = jj;
        for (i = n; i <= m; i++)
        {
            s[sIndex + (l - 1)] = w[i - 1];
            l += 1;
        }

        sing = s[sIndex + (jj - 1)] switch
        {
            0.0 => true,
            _ => sing
        };

        return sing;
    }
}