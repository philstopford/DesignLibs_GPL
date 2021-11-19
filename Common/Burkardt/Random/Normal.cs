using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.RandomNS;

public static class Normal
{
    public static double[] normal(int dim_num, int n, double[] r, double[] mu, ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL creates normally distributed points in DIM_NUM space.
        //
        //  Discussion:
        //
        //    The multivariate normal distribution for the DIM_NUM dimensional vector X
        //    has the form:
        //
        //      pdf(X) = (2*pi*det(V))^(-DIM_NUM/2) * exp(-0.5*(X-MU)'*inverse(V)*(X-MU))
        //
        //    where MU is the mean vector, and V is a positive definite symmetric
        //    matrix called the variance-covariance matrix.
        //
        //    This routine requires that the user supply the upper triangular
        //    Cholesky factor R, which has the property that
        //
        //      V = R' * R
        //
        //    This factorization always exists if V is actually symmetric and
        //    positive definite.  This factorization can be computed by the
        //    routine SPO_FA.
        //
        //    The user also supplies the mean vector MU.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 August 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the dimension of the space.
        //
        //    Input, int N, the number of points.
        //
        //    Input, double R[DIM_NUM*DIM_NUM], the upper triangular Cholesky factor
        //    of the variance-covariance matrix.
        //
        //    Input, double MU[DIM_NUM], the mean vector.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double NORMAL[DIM_NUM*N], the random points.
        //
    {
        int j;

        double[] v = new double[dim_num];
        double[] x = new double[dim_num * n];
        //
        //  Get a matrix V of normal data.
        //  Compute X = MU + R' * V.
        //  We actually carry out this computation in the equivalent form X' * R.
        //
        for (j = 0; j < n; j++)
        {
            typeMethods.r8vec_normal_01(dim_num, ref seed, ref v);

            int i;
            for (i = 0; i < dim_num; i++)
            {
                x[i + j * dim_num] = mu[i];
                int k;
                for (k = 0; k <= i; k++)
                {
                    x[i + j * dim_num] += v[k] * r[k + i * dim_num];
                }
            }
        }

        return x;
    }

    public static double[] normal_circular(int dim_num, int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_CIRCULAR creates circularly normal points in 2 space.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 August 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz and Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    US Department of Commerce, 1964, page 936.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the dimension of the space, which must be 2.
        //
        //    Input, int N, the number of points.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double NORMAL_CIRULAR[DIM_NUM*N], the random points.
        //
    {
        int j;

        double[] r = new double[n];
        double[] t = new double[n];
        double[] x = new double[dim_num * n];
        //
        //  The angle varies uniformly from 0 to 2 pi.
        //
        UniformRNG.r8vec_uniform_01(n, ref seed, ref t);

        for (j = 0; j < n; j++)
        {
            t[j] = 2.0 * Math.PI * t[j];
        }

        //
        //  The radius is normally distributed.
        //
        typeMethods.r8vec_normal_01(n, ref seed, ref r);

        for (j = 0; j < n; j++)
        {
            x[0 + j * dim_num] = r[j] * Math.Cos(t[j]);
            x[1 + j * dim_num] = r[j] * Math.Sin(t[j]);
        }

        return x;
    }

    public static double[] normal_multivariate(int m, int n, double[] r, double[] mu,
            ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_MULTIVARIATE samples a multivariate normal distribution.
        //
        //  Discussion:
        //
        //    The multivariate normal distribution for the M dimensional vector X
        //    has the form:
        //
        //      pdf(X) = (2*pi*det(V))^(-M/2) * exp(-0.5*(X-MU)'*inverse(V)*(X-MU))
        //
        //    where MU is the mean vector, and V is a positive definite symmetric
        //    matrix called the variance-covariance matrix.
        //
        //    This routine samples points associated with the M dimensional
        //    normal distribution with mean MU and covariance matrix V.
        //
        //    This routine requires that the user supply the upper triangular
        //    Cholesky factor R of V, which has the property that
        //
        //      V = R' * R
        //
        //    This factorization always exists if V is actually symmetric and
        //    positive definite.  This factorization can be computed by the
        //    routine R8PO_FA.
        //
        //    The user also supplies the mean vector MU.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 August 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Russell Cheng,
        //    Random Variate Generation,
        //    in Handbook of Simulation,
        //    edited by Jerry Banks,
        //    Wiley, 1998, pages 167-168.
        //
        //  Parameters:
        //
        //    Input, int M, the dimension of the space.
        //
        //    Input, int N, the number of points.
        //
        //    Input, double R[M*M], the upper triangular Cholesky factor
        //    of the variance-covariance matrix.
        //
        //    Input, double MU[M], the mean vector.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double NORMAL_MULTIVARIATE[DIM_NUM*N], corresponding
        //    points associated with the multivariate normal distribution.
        //
    {
        int j;

        double[] v = new double[m];
        double[] x = new double[m * n];
        //
        //  Compute X = MU + R' * V.
        //  We actually carry out this computation in the equivalent form MU + V' * R.
        //
        for (j = 0; j < n; j++)
        {
            typeMethods.r8vec_normal_01(m, ref seed, ref v);

            int i;
            for (i = 0; i < m; i++)
            {
                x[i + j * m] = mu[i];
                int k;
                for (k = 0; k <= i; k++)
                {
                    x[i + j * m] += v[k] * r[k + i * m];
                }
            }
        }
            
        return x;
    }

    public static double[] normal_simple(int dim_num, int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_SIMPLE creates normally distributed points in DIM_NUM space.
        //
        //  Discussion:
        //
        //    The multivariate normal distribution has the form:
        //
        //      f(x) = (2*pi*det(V))^(-DIM_NUM/2) * exp(-0.5*(x-mu)'*inverse(V)*(x-mu))
        //
        //    where mu is the mean vector, and V is a positive definite symmetric
        //    matrix called the variance-covariance matrix.
        //
        //    This routine implements the simplest version of a multivariate
        //    normal distribution.  The variance-covariance matrix is the identity,
        //    and the mean vector is entirely zero.  Thus, a sample on N points
        //    is simply DIM_NUM*N scalar values generated under the univariate
        //    normal distribution with zero mean and unit variance.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 August 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the dimension of the space.
        //
        //    Input, int N, the number of points.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double NORMAL_SIMPLE[DIM_NUM*N], the random points.
        //
    {
        double[] x = new double[dim_num * n];

        typeMethods.r8vec_normal_01(dim_num * n, ref seed, ref x);

        return x;
    }
}