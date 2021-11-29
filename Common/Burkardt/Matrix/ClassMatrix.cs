﻿using System;
using Burkardt.Types;
using Burkardt.Weight;

namespace Burkardt.MatrixNS;

public static partial class Matrix
{
    public static double class_matrix(int kind, int m, double alpha, double beta, ref double[] aj,
            ref double[] bj )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CLASS_MATRIX computes the Jacobi matrix for a quadrature rule.
        //
        //  Discussion:
        //
        //    This routine computes the diagonal AJ and sub-diagonal BJ
        //    elements of the order M tridiagonal symmetric Jacobi matrix
        //    associated with the polynomials orthogonal with respect to
        //    the weight function specified by KIND.
        //
        //    For weight functions 1-7, M elements are defined in BJ even
        //    though only M-1 are needed.  For weight function 8, BJ(M) is
        //    set to zero.
        //
        //    The zero-th moment of the weight function is returned in ZEMU.
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
        //    Input, int KIND, the rule.
        //    1, Legendre,             (a,b)       1.0
        //    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
        //    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
        //    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
        //    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
        //    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
        //    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
        //    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
        //
        //    Input, int M, the order of the Jacobi matrix.
        //
        //    Input, double ALPHA, the value of Alpha, if needed.
        //
        //    Input, double BETA, the value of Beta, if needed.
        //
        //    Output, double AJ[M], BJ[M], the diagonal and subdiagonal
        //    of the Jacobi matrix.
        //
        //    Output, double CLASS_MATRIX, the zero-th moment.
        //
    {
        double ab;
        double abi;
        double abj;
        int i;

        double zemu = 0;

        double temp = typeMethods.r8_epsilon();

        PARCHK.parchk(kind, 2 * m - 1, alpha, beta);

        const double temp2 = 0.5;

        if (500.0 * temp < Math.Abs(Math.Pow(Helpers.Gamma(temp2), 2) - Math.PI))
        {
            Console.WriteLine("");
            Console.WriteLine("CLASS_MATRIX - Fatal error!");
            Console.WriteLine("  Gamma function does not match machine parameters.");
            return 1;
        }

        switch (kind)
        {
            case 1:
            {
                ab = 0.0;

                zemu = 2.0 / (ab + 1.0);

                for (i = 0; i < m; i++)
                {
                    aj[i] = 0.0;
                }

                for (i = 1; i <= m; i++)
                {
                    abi = i + ab * (i % 2);
                    abj = 2 * i + ab;
                    bj[i - 1] = Math.Sqrt(abi * abi / (abj * abj - 1.0));
                }

                break;
            }
            case 2:
            {
                zemu = Math.PI;

                for (i = 0; i < m; i++)
                {
                    aj[i] = 0.0;
                }

                bj[0] = Math.Sqrt(0.5);
                for (i = 1; i < m; i++)
                {
                    bj[i] = 0.5;
                }

                break;
            }
            case 3:
            {
                ab = alpha * 2.0;
                zemu = Math.Pow(2.0, ab + 1.0) * Math.Pow(Helpers.Gamma(alpha + 1.0), 2)
                       / Helpers.Gamma(ab + 2.0);

                for (i = 0; i < m; i++)
                {
                    aj[i] = 0.0;
                }

                bj[0] = Math.Sqrt(1.0 / (2.0 * alpha + 3.0));
                for (i = 2; i <= m; i++)
                {
                    bj[i - 1] = Math.Sqrt(i * (i + ab) / (4.0 * Math.Pow(i + alpha, 2) - 1.0));
                }

                break;
            }
            case 4:
            {
                ab = alpha + beta;
                abi = 2.0 + ab;
                zemu = Math.Pow(2.0, ab + 1.0) * Helpers.Gamma(alpha + 1.0)
                                               * Helpers.Gamma(beta + 1.0) / Helpers.Gamma(abi);
                aj[0] = (beta - alpha) / abi;
                bj[0] = Math.Sqrt(4.0 * (1.0 + alpha) * (1.0 + beta)
                                  / ((abi + 1.0) * abi * abi));
                double a2b2 = beta * beta - alpha * alpha;

                for (i = 2; i <= m; i++)
                {
                    abi = 2.0 * i + ab;
                    aj[i - 1] = a2b2 / ((abi - 2.0) * abi);
                    abi *= abi;
                    bj[i - 1] = Math.Sqrt(4.0 * i * (i + alpha) * (i + beta) * (i + ab)
                                          / ((abi - 1.0) * abi));
                }

                break;
            }
            case 5:
            {
                zemu = Helpers.Gamma(alpha + 1.0);

                for (i = 1; i <= m; i++)
                {
                    aj[i - 1] = 2.0 * i - 1.0 + alpha;
                    bj[i - 1] = Math.Sqrt(i * (i + alpha));
                }

                break;
            }
            case 6:
            {
                zemu = Helpers.Gamma((alpha + 1.0) / 2.0);

                for (i = 0; i < m; i++)
                {
                    aj[i] = 0.0;
                }

                for (i = 1; i <= m; i++)
                {
                    bj[i - 1] = Math.Sqrt((i + alpha * (i % 2)) / 2.0);
                }

                break;
            }
            case 7:
            {
                ab = alpha;
                zemu = 2.0 / (ab + 1.0);

                for (i = 0; i < m; i++)
                {
                    aj[i] = 0.0;
                }

                for (i = 1; i <= m; i++)
                {
                    abi = i + ab * (i % 2);
                    abj = 2 * i + ab;
                    bj[i - 1] = Math.Sqrt(abi * abi / (abj * abj - 1.0));
                }

                break;
            }
            case 8:
            {
                ab = alpha + beta;
                zemu = Helpers.Gamma(alpha + 1.0) * Helpers.Gamma(-(ab + 1.0))
                       / Helpers.Gamma(-beta);
                double apone = alpha + 1.0;
                double aba = ab * apone;
                aj[0] = -apone / (ab + 2.0);
                bj[0] = -aj[0] * (beta + 1.0) / (ab + 2.0) / (ab + 3.0);
                double abti;
                for (i = 2; i <= m; i++)
                {
                    abti = ab + 2.0 * i;
                    aj[i - 1] = aba + 2.0 * (ab + i) * (i - 1);
                    aj[i - 1] = -aj[i - 1] / abti / (abti - 2.0);
                }

                for (i = 2; i <= m - 1; i++)
                {
                    abti = ab + 2.0 * i;
                    bj[i - 1] = i * (alpha + i) / (abti - 1.0) * (beta + i)
                        / (abti * abti) * (ab + i) / (abti + 1.0);
                }

                bj[m - 1] = 0.0;
                for (i = 0; i < m; i++)
                {
                    bj[i] = Math.Sqrt(bj[i]);
                }

                break;
            }
            case 9:
            {
                zemu = Math.PI / 2.0;

                for ( i = 0; i < m; i++ )
                {
                    aj[i] = 0.0;
                }

                for ( i = 0; i < m; i++ )
                {
                    bj[i] = 0.5;
                }

                break;
            }
        }

        return zemu;
    }
}