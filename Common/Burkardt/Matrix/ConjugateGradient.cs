using Burkardt.Types;

namespace Burkardt.MatrixNS;

public static class ConjugateGradient
{
    public static void cg_gb(int n, int ml, int mu, double[] a, double[] b, ref double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CG_GB uses the conjugate gradient method for a general banded (GB) matrix.
        //
        //  Discussion:
        //
        //    The linear system has the form A*x=b, where A is a positive-definite
        //    symmetric matrix.
        //
        //    The method is designed to reach the solution to the linear system
        //      A * x = b
        //    after N computational steps.  However, roundoff may introduce
        //    unacceptably large errors for some problems.  In such a case,
        //    calling the routine a second time, using the current solution estimate
        //    as the new starting guess, should result in improved results.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Frank Beckman,
        //    The Solution of Linear Equations by the Conjugate Gradient Method,
        //    in Mathematical Methods for Digital Computers,
        //    edited by John Ralston, Herbert Wilf,
        //    Wiley, 1967,
        //    ISBN: 0471706892,
        //    LC: QA76.5.R3.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int ML, MU, the lower and upper bandwidths.
        //
        //    Input, double A[(2*ML+MU+1)*N], the band matrix.
        //
        //    Input, double B[N], the right hand side vector.
        //
        //    Input/output, double X[N].
        //    On input, an estimate for the solution, which may be 0.
        //    On output, the approximate solution vector.  
        //
    {
        double alpha;
        double[] ap;
        double beta;
        int i;
        int it;
        double[] p;
        double pap;
        double pr;
        double[] r;
        double rap;
        //
        //  Initialize
        //    AP = A * x,
        //    R  = b - A * x,
        //    P  = b - A * x.
        //
        ap = MatbyVector.mv_gb(n, n, ml, mu, a, x);

        r = new double[n];
        for (i = 0; i < n; i++)
        {
            r[i] = b[i] - ap[i];
        }

        p = new double[n];
        for (i = 0; i < n; i++)
        {
            p[i] = b[i] - ap[i];
        }

        //
        //  Do the N steps of the conjugate gradient method.
        //
        for (it = 1; it <= n; it++)
        {
            //
            //  Compute the matrix*vector product AP = A*P.
            //
            ap = MatbyVector.mv_gb(n, n, ml, mu, a, p);
            //
            //  Compute the dot products
            //    PAP = P*AP,
            //    PR  = P*R
            //  Set
            //    ALPHA = PR / PAP.
            //
            pap = typeMethods.r8vec_dot_product(n, p, ap);
            pr = typeMethods.r8vec_dot_product(n, p, r);

            if (pap == 0.0)
            {
                break;
            }

            alpha = pr / pap;
            //
            //  Set
            //    X = X + ALPHA * P
            //    R = R - ALPHA * AP.
            //
            for (i = 0; i < n; i++)
            {
                x[i] += alpha * p[i];
            }

            for (i = 0; i < n; i++)
            {
                r[i] -= alpha * ap[i];
            }
            //
            //  Compute the vector dot product
            //    RAP = R*AP
            //  Set
            //    BETA = - RAP / PAP.
            //

            rap = typeMethods.r8vec_dot_product(n, r, ap);

            beta = -rap / pap;
            //
            //  Update the perturbation vector
            //    P = R + BETA * P.
            //
            for (i = 0; i < n; i++)
            {
                p[i] = r[i] + beta * p[i];
            }
        }

    }

    public static void cg_ge(int n, double[] a, double[] b, ref double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CG_GE uses the conjugate gradient method for a general storage (GE) matrix.
        //
        //  Discussion:
        //
        //    The linear system has the form A*x=b, where A is a positive-definite
        //    symmetric matrix, stored as a full storage matrix.
        //
        //    The method is designed to reach the solution to the linear system
        //      A * x = b
        //    after N computational steps.  However, roundoff may introduce
        //    unacceptably large errors for some problems.  In such a case,
        //    calling the routine a second time, using the current solution estimate
        //    as the new starting guess, should result in improved results.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Frank Beckman,
        //    The Solution of Linear Equations by the Conjugate Gradient Method,
        //    in Mathematical Methods for Digital Computers,
        //    edited by John Ralston, Herbert Wilf,
        //    Wiley, 1967,
        //    ISBN: 0471706892,
        //    LC: QA76.5.R3.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N*N], the matrix.
        //
        //    Input, double B[N], the right hand side vector.
        //
        //    Input/output, double X[N].
        //    On input, an estimate for the solution, which may be 0.
        //    On output,  the approximate solution vector.  
        //
    {
        double alpha;
        double[] ap;
        double beta;
        int i;
        int it;
        double[] p;
        double pap;
        double pr;
        double[] r;
        double rap;
        //
        //  Initialize
        //    AP = A * x,
        //    R  = b - A * x,
        //    P  = b - A * x.
        //
        ap = MatbyVector.mv_ge(n, n, a, x);

        r = new double[n];
        for (i = 0; i < n; i++)
        {
            r[i] = b[i] - ap[i];
        }

        p = new double[n];
        for (i = 0; i < n; i++)
        {
            p[i] = b[i] - ap[i];
        }

        //
        //  Do the N steps of the conjugate gradient method.
        //
        for (it = 1; it <= n; it++)
        {
            //
            //  Compute the matrix*vector product AP = A*P.
            //

            ap = MatbyVector.mv_ge(n, n, a, p);
            //
            //  Compute the dot products
            //    PAP = P*AP,
            //    PR  = P*R
            //  Set
            //    ALPHA = PR / PAP.
            //
            pap = typeMethods.r8vec_dot_product(n, p, ap);
            pr = typeMethods.r8vec_dot_product(n, p, r);

            if (pap == 0.0)
            {
                break;
            }

            alpha = pr / pap;
            //
            //  Set
            //    X = X + ALPHA * P
            //    R = R - ALPHA * AP.
            //
            for (i = 0; i < n; i++)
            {
                x[i] += alpha * p[i];
            }

            for (i = 0; i < n; i++)
            {
                r[i] -= alpha * ap[i];
            }

            //
            //  Compute the vector dot product
            //    RAP = R*AP
            //  Set
            //    BETA = - RAP / PAP.
            //
            rap = typeMethods.r8vec_dot_product(n, r, ap);

            beta = -rap / pap;
            //
            //  Update the perturbation vector
            //    P = R + BETA * P.
            //
            for (i = 0; i < n; i++)
            {
                p[i] = r[i] + beta * p[i];
            }
        }
    }

    public static void cg_st(int n, int nz_num, int[] row, int[] col, double[] a, double[] b,
            ref double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CG_ST uses the conjugate gradient method for a sparse triplet (ST) matrix.
        //
        //  Discussion:
        //
        //    The linear system has the form A*x=b, where A is a positive-definite
        //    symmetric matrix, stored as a full storage matrix.
        //
        //    The method is designed to reach the solution to the linear system
        //      A * x = b
        //    after N computational steps.  However, roundoff may introduce
        //    unacceptably large errors for some problems.  In such a case,
        //    calling the routine a second time, using the current solution estimate
        //    as the new starting guess, should result in improved results.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Frank Beckman,
        //    The Solution of Linear Equations by the Conjugate Gradient Method,
        //    in Mathematical Methods for Digital Computers,
        //    edited by John Ralston, Herbert Wilf,
        //    Wiley, 1967,
        //    ISBN: 0471706892,
        //    LC: QA76.5.R3.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzeros.
        //
        //    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and column 
        //    indices of the nonzero entries.
        //
        //    Input, double A[NZ_NUM], the nonzero entries.
        //
        //    Input, double B[N], the right hand side vector.
        //
        //    Input/output, double X[N].
        //    On input, an estimate for the solution, which may be 0.
        //    On output, the approximate solution vector.  
        //
    {
        double alpha;
        double[] ap;
        double beta;
        int i;
        int it;
        double[] p;
        double pap;
        double pr;
        double[] r;
        double rap;
        //
        //  Initialize
        //    AP = A * x,
        //    R  = b - A * x,
        //    P  = b - A * x.
        //
        ap = MatbyVector.mv_st(n, n, nz_num, row, col, a, x);

        r = new double[n];
        for (i = 0; i < n; i++)
        {
            r[i] = b[i] - ap[i];
        }

        p = new double[n];
        for (i = 0; i < n; i++)
        {
            p[i] = b[i] - ap[i];
        }

        //
        //  Do the N steps of the conjugate gradient method.
        //
        for (it = 1; it <= n; it++)
        {
            //
            //  Compute the matrix*vector product AP = A*P.
            //
            ap = MatbyVector.mv_st(n, n, nz_num, row, col, a, p);
            //
            //  Compute the dot products
            //    PAP = P*AP,
            //    PR  = P*R
            //  Set
            //    ALPHA = PR / PAP.
            //
            pap = typeMethods.r8vec_dot_product(n, p, ap);
            pr = typeMethods.r8vec_dot_product(n, p, r);

            if (pap == 0.0)
            {
                break;
            }

            alpha = pr / pap;
            //
            //  Set
            //    X = X + ALPHA * P
            //    R = R - ALPHA * AP.
            //
            for (i = 0; i < n; i++)
            {
                x[i] += alpha * p[i];
            }

            for (i = 0; i < n; i++)
            {
                r[i] -= alpha * ap[i];
            }

            //
            //  Compute the vector dot product
            //    RAP = R*AP
            //  Set
            //    BETA = - RAP / PAP.
            //
            rap = typeMethods.r8vec_dot_product(n, r, ap);

            beta = -rap / pap;
            //
            //  Update the perturbation vector
            //    P = R + BETA * P.
            //
            for (i = 0; i < n; i++)
            {
                p[i] = r[i] + beta * p[i];
            }
        }
    }
}