﻿using System;
using System.Numerics;
using Burkardt.Types;

namespace Burkardt.PowerMethodNS;

public static class PowerMethod
{
    public static void power_method(int n, double[] a, double[] y, int it_max, double tol,
            ref double lambda, ref int it_num )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POWER_METHOD applies several steps of the power method.
        //
        //  Discussion:
        //
        //    For a given NxN matrix A and an N vector Y, the power method produces
        //    a series of estimates for LAMBDA, the largest eigenvalue, and Y,
        //    the eigenvector corresponding to LAMBDA.
        //
        //    The iteration repeats the following steps
        //
        //      AY     = A * Y
        //      LAMBDA = || AY ||
        //      Y      = AY / LAMBDA
        //
        //    If the matrix A has a single real eigenvalue of maximum modulus,
        //    then this iteration will generally produce a good estimate for that
        //    eigenvalue and its corresponding eigenvector.
        //
        //    If there are multiple distinct eigenvalues of the same modulus,
        //    perhaps two values of opposite sign, or complex eigenvalues, then
        //    the situation is more complicated.
        //
        //    Separate issues:
        //
        //    * when estimating the value of LAMBDA, we use the Rayleigh quotient,
        //    LAMBDA = ( y' * A * y ) / ( y' * y ).  Since we normalize Y, the
        //    bottom of the fraction is 1.  Using this estimate allows us to
        //    easily capture the sign of LAMDBA.  Using the eucldean norm
        //    instead, for instance, would always give a positive value.
        //
        //    * If the dominant eigenvalue is negative, then the iteration 
        //    as given will produce eigenvector iterates that alternate in sign.  
        //   
        //    * It is worth knowing whether the successive eigenvector estimates
        //    are tending to some value.  Since an eigenvector is really a direction,
        //    we need to normalize the vectors, and we need to somehow treat both
        //    a vector and its negative as holding the same information.  This
        //    means that the proper way of measuring the difference between two
        //    eigenvector estimates is to normalize them both, and then compute
        //    the cosine between them as y1'y2, followed by the sine, which is
        //    sqrt ( 1 - ( y1'y2)^2 ).  If this sine is small, the vectors y1 and y2
        //    are "close" in the sense of direction.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 July 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N*N], the matrix.
        //
        //    Input/output, double Y[N], the estimate for the eigenvector.
        //
        //    Input, int IT_MAX, the maximum number of iterations to take.
        //    1 <= IT_MAX.
        //
        //    Input, double TOL, an error tolerance.
        //
        //    Output, double *LAMBDA, the estimate for the eigenvalue.
        //
        //    Output, int *IT_NUM, the number of iterations taken.
        //
    {
        int i;

        double[] ay = new double[n];
        double[] y_old = new double[n];
        
        //
        //  Force Y to be a vector of unit norm.
        //
        double norm = typeMethods.r8vec_norm_l2(n, y);

        for (i = 0; i < n; i++)
        {
            y[i] /= norm;
        }

        int it;

        for (i = 0; i < n; i++)
        {
            y_old[i] = y[i];
        }

        typeMethods.r8mat_mv(n, n, a, y, ref ay);
        lambda = typeMethods.r8vec_dot(n, y, ay);
        norm = typeMethods.r8vec_norm_l2(n, ay);
        for (i = 0; i < n; i++)
        {
            y[i] = ay[i] / norm;
        }

        switch (lambda)
        {
            case < 0.0:
            {
                for (i = 0; i < n; i++)
                {
                    y[i] = -y[i];
                }

                break;
            }
        }

        for (it = 1; it <= it_max; it++)
        {
            double lambda_old = lambda;
            for (i = 0; i < n; i++)
            {
                y_old[i] = y[i];
            }

            typeMethods.r8mat_mv(n, n, a, y, ref ay);
            lambda = typeMethods.r8vec_dot(n, y, ay);
            norm = typeMethods.r8vec_norm_l2(n, ay);
            for (i = 0; i < n; i++)
            {
                y[i] = ay[i] / norm;
            }

            switch (lambda)
            {
                case < 0.0:
                {
                    for (i = 0; i < n; i++)
                    {
                        y[i] = -y[i];
                    }

                    break;
                }
            }

            double val_dif = Math.Abs(lambda - lambda_old);

            if (val_dif <= tol)
            {
                break;
            }
        }

        for (i = 0; i < n; i++)
        {
            y[i] = ay[i] / lambda;
        }

        it_num = it;
    }

    public static void power_method2(int n, double[] a, double[] x_init, int it_max,
            double tol, ref Complex lambda, Complex[] v, ref int it_num )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POWER_METHOD2 applies the power method for possibly complex eigenvalues.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 July 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Eric VanDeVelde,
        //    Concurrent Scientific Programming,
        //    Springer, 1994,
        //    ISBN: 0-387-94195-9,
        //    LC: QA76.58.V35.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N*N], the matrix.
        //
        //    Input, double X_INIT[N], the initial estimate for the eigenvector.
        //
        //    Input, int IT_MAX, the maximum number of iterations to take.
        //    1 <= IT_MAX.
        //
        //    Input, double TOL, an error tolerance.
        //
        //    Output, complex <double> *LAMBDA, the estimate for the eigenvalue.
        //
        //    Output, complex <double> V[N], the estimate for the eigenvector.
        //
        //    Output, int *IT_NUM, the number of iterations taken.
        //
    {
        int it;

        double[] x = new double[n];
        double[] y = new double[n];
        double[] z = new double[n];

        it_num = 0;
        //
        //  Compute data necessary to start the iteration.
        //
        typeMethods.r8vec_copy(n, x_init, ref x);

        double pi_xx = typeMethods.r8vec_dot(n, x, x);

        typeMethods.r8vec_divide(n, ref x, pi_xx);

        typeMethods.r8mat_mv(n, n, a, x, ref y);

        double pi_xy = typeMethods.r8vec_dot(n, x, y);
        double pi_yy = typeMethods.r8vec_dot(n, y, y);

        for (it = 1; it <= it_max; it++)
        {
            int i;
            if (pi_yy - pi_xy * pi_xy < tol * tol * pi_yy)
            {
                lambda = new Complex(pi_xy, 0.0);
                for (i = 0; i < n; i++)
                {
                    v[i] = new Complex(y[i], 0.0) / pi_yy;
                }

                return;
            }

            typeMethods.r8mat_mv(n, n, a, y, ref z);

            double pi_xz = typeMethods.r8vec_dot(n, x, z);
            double pi_yz = typeMethods.r8vec_dot(n, y, z);
            double pi_zz = typeMethods.r8vec_dot(n, z, z);

            double alpha = -(pi_yz - pi_xy * pi_xz) / (pi_yy - pi_xy * pi_xy);

            double beta = (pi_xy * pi_yz - pi_yy * pi_xz) / (pi_yy - pi_xy * pi_xy);

            double gamma = pi_zz + alpha * alpha * pi_yy + beta * beta
                           + 2.0 * (alpha * pi_yz + beta * pi_xz + alpha * beta * pi_xy);

            if (gamma < tol * tol * pi_zz && alpha * alpha < 4.0 * beta)
            {
                double lambda_real = -alpha / 2.0;
                double lambda_imag = Math.Sqrt(4.0 * beta - alpha * alpha) / 2.0;
                lambda = new Complex(lambda_real, lambda_imag);

                for (i = 0; i < n; i++)
                {
                    v[i] = (lambda * y[i] - z[i])
                           / Math.Sqrt(beta * pi_yy + alpha * pi_yz + pi_zz);
                }

                return;
            }

            typeMethods.r8vec_copy(n, y, ref x);
            typeMethods.r8vec_divide(n, ref x, Math.Sqrt(pi_yy));

            typeMethods.r8vec_copy(n, z, ref y);
            typeMethods.r8vec_divide(n, ref y, Math.Sqrt(pi_yy));

            pi_xy = pi_yz / pi_yy;
            pi_yy = pi_zz / pi_yy;

            it_num = it;
        }

        Console.WriteLine("");
        Console.WriteLine("POWER_METHOD2 - Fatal error!");
        Console.WriteLine("  Convergence was not reached.");
    }
}