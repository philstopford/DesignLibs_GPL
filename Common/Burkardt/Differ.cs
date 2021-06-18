using System;
using Burkardt.Types;

namespace Burkardt
{
    public static class Differ
    {
        public static void differ_backward(double h, int o, int p, ref double[] c, ref double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIFFER_BACKWARD computes backward difference coefficients.
        //
        //  Discussion:
        //
        //    We determine coefficients C to approximate the derivative at X0
        //    of order O and precision P, using equally spaced backward
        //    differences, so that 
        //
        //      d^o f(x)/dx^o = sum ( 0 <= i <= o+p-1 ) c(i) f(x-ih) + O(h^(p))
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 November 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double H, the spacing.  0 < H.
        //
        //    Input, int O, the order of the derivative to be 
        //    approximated.  1 <= O.
        //
        //    Input, int P, the order of the error, as a power of H.
        //
        //    Output, double C[O+P], the coefficients.
        //
        //    Output, double X[O+P], the evaluation points.
        //
        {
            double[] b;
            int i;
            int info = 0;
            int job;
            int n;
            double t;

            n = o + p;

            for (i = 0; i < n; i++)
            {
                x[i] = (double) (i + 1 - n) * h;
            }

            b = new double[n];

            for (i = 0; i < n; i++)
            {
                b[i] = 0.0;
            }

            b[o] = 1.0;

            job = 0;
            typeMethods.r8vm_sl(n, x, b, job, c, ref info);

            if (info != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("DIFFER_BACKWARD - Fatal error!");
                Console.WriteLine("  Vandermonde linear system is singular.");
                return;
            }

            t = typeMethods.r8_factorial(o);
            for (i = 0; i < n; i++)
            {
                c[i] = c[i] * t;
            }
        }

        public static void differ_central(double h, int o, int p, ref double[] c, ref double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIFFER_CENTRAL computes central difference coefficients.
        //
        //  Discussion:
        //
        //    We determine coefficients C to approximate the derivative at X0
        //    of order O and precision P, using equally spaced central
        //    differences, so that 
        //
        //      d^o f(x)/dx^o = sum ( 0 <= i <= o+p-1 ) c(i) f(x+(2*i-o-p+1)*h/2) 
        //        + O(h^(p))
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 November 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double H, the spacing.  0 < H.
        //
        //    Input, int O, the order of the derivative to 
        //    be approximated.  1 <= O.
        //
        //    Input, int P, the order of the error, as a power of H.
        //
        //    Output, double C[O+P], the coefficients.
        //
        //    Output, double X[O+P], the evaluation points.
        //
        {
            double[] b;
            int i;
            int info = 0;
            int job;
            int n;
            double t;

            n = o + p;

            for (i = 0; i < n; i++)
            {
                x[i] = (double) (-n + 1 + 2 * i) * h / 2.0;
            }

            b = new double[n];
            for (i = 0; i < n; i++)
            {
                b[i] = 0.0;
            }

            b[o] = 1.0;

            job = 0;
            typeMethods.r8vm_sl(n, x, b, job, c, ref info);

            if (info != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("DIFFER_CENTRAL - Fatal error!");
                Console.WriteLine("  Vandermonde linear system is singular.");
                return;
            }

            t = typeMethods.r8_factorial(o);
            for (i = 0; i < n; i++)
            {
                c[i] = c[i] * t;
            }
        }

        public static void differ_forward(double h, int o, int p, ref double[] c, ref double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIFFER_FORWARD computes forward difference coefficients.
        //
        //  Discussion:
        //
        //    We determine coefficients C to approximate the derivative at X0
        //    of order O and precision P, using equally spaced forward
        //    differences, so that 
        //
        //      d^o f(x)/dx^o = sum ( 0 <= i <= o+p-1 ) c(i) f(x+ih) + O(h^(p))
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 November 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, real H, the spacing.  0 < H.
        //
        //    Input, integer O, the order of the derivative to be approximated.
        //    1 <= O.
        //
        //    Input, integer P, the order of the error, as a power of H.
        //
        //    Output, real C[O+P], the coefficients.
        //
        //    Output, real X[O+P], the evaluation points.
        //
        {
            double[] b;
            int i;
            int info = 0;
            int job;
            int n;
            double t;

            n = o + p;

            for (i = 0; i < n; i++)
            {
                x[i] = (double) (i) * h;
            }

            b = new double[n];
            for (i = 0; i < n; i++)
            {
                b[i] = 0.0;
            }

            b[o] = 1.0;

            job = 0;
            typeMethods.r8vm_sl(n, x, b, job, c, ref info);

            if (info != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("DIFFER_FORWARD - Fatal error!");
                Console.WriteLine("  Vandermonde linear system is singular.");
                return;
            }

            t = typeMethods.r8_factorial(o);
            for (i = 0; i < n; i++)
            {
                c[i] = c[i] * t;
            }
        }

        public static double[] differ_inverse(int n, double[] stencil)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DIFFER_INVERSE returns the inverse of the DIFFER matrix.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 November 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the matrix.
            //
            //    Input, double STENCIL[N], the values that define A.
            //
            //    Output, double DIFFER_INVERSE[N*N], the matrix.
            //
        {
            double[] a;
            int i;
            int indx;
            int j;
            int k;

            a = new double[n * n];

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < n; i++)
                {
                    if (j == 0)
                    {
                        a[i + j * n] = 1.0;
                    }
                    else
                    {
                        a[i + j * n] = 0.0;
                    }
                }
            }

            for (i = 0; i < n; i++)
            {
                indx = 0;

                for (k = 0; k < n; k++)
                {
                    if (k != i)
                    {
                        for (j = indx + 1; 0 <= j; j--)
                        {
                            a[i + j * n] = -stencil[k] * a[i + j * n] / (stencil[i] - stencil[k]);

                            if (0 < j)
                            {
                                a[i + j * n] = a[i + j * n] + a[i + (j - 1) * n] / (stencil[i] - stencil[k]);
                            }
                        }

                        indx = indx + 1;
                    }
                }
            }

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < n; i++)
                {
                    a[i + j * n] = a[i + j * n] / stencil[i];
                }
            }

            return a;
        }

        public static double[] differ_matrix(int n, double[] stencil)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DIFFER_MATRIX computes the stencil matrix from the stencil vector.
            //
            //  Discussion:
            //
            //    If N = 4, and STENCIL = ( -3, -2, -1, 1 ), then A will be
            //
            //    -3  -2  -1  1
            //     9   4   1  1
            //   -27  -8  -1  1
            //    81  16   1  1
            //
            //    This matrix is a generalized form of a Vandermonde matrix A2:
            //
            //     1   1   1  1
            //    -3  -2  -1  1
            //     9   4   1  1
            //   -27  -8  -1  1    
            //
            //    and if A * x = b, the A2 * x2 = b, where x2(i) = x(i) * stencil(i)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 November 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of stencil points.
            //
            //    Input, double STENCIL[N], the stencil vector.
            //    The entries in this vector must be distinct.
            //    No entry of STENCIL may be 0.
            //
            //    Output, double DIFFER_MATRIX[N*N], the stencil matrix.
            //
        {
            double[] a;
            int i;
            int j;

            a = new double[n * n];

            for (j = 0; j < n; j++)
            {
                a[0 + j * n] = stencil[j];
                for (i = 1; i < n; i++)
                {
                    a[i + j * n] = a[i - 1 + j * n] * stencil[j];
                }
            }

            return a;
        }

        public static double[] differ_solve(int n, double[] stencil, int order )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIFFER_SOLVE solves for finite difference coefficients.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 November 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of stencil points.
        //
        //    Input, double STENCIL[N], the stencil vector.
        //    The entries in this vector must be distinct.
        //    No entry of STENCIL may be 0.
        //
        //    Input, int ORDER, the order of the derivative to
        //    be approximated.  1 <= ORDER <= N.
        //
        //    Output, double DIFFER_SOLVE[N], the coefficients to be used
        //    to multiply U(STENCIL(I))-U(0), so that the sum forms an
        //    approximation to the derivative of order ORDER, with error 
        //    of order H^(N+1-ORDER).
        //
        {
            double[] a;
            double[] b;
            double[] c;
            int i;

            a = differ_matrix(n, stencil);

            b = new double[n];
            for (i = 0; i < n; i++)
            {
                b[i] = 0.0;
            }

            b[order - 1] = 1.0;
            //
            //  Solve A * C = B.
            //
            c = typeMethods.r8mat_fs_new(n, a, b);

            return c;
        }

        public static void differ_stencil(double x0, int o, int p, double[] x, ref double[] c )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIFFER_STENCIL computes finite difference coefficients.
        //
        //  Discussion:
        //
        //    We determine coefficients C to approximate the derivative at X0
        //    of order O and precision P, using finite differences, so that 
        //
        //      d^o f(x)/dx^o (x0) = sum ( 0 <= i <= o+p-1 ) c(i) f(x(i)) 
        //        + O(h^(p))
        //
        //    where H is the maximum spacing between X0 and any X(I).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 November 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X0, the point where the derivative is to 
        //    be approximated.
        //
        //    Input, int O, the order of the derivative to be 
        //    approximated.  1 <= O.
        //
        //    Input, int P, the order of the error, as a power of H.
        //
        //    Input, double X[O+P], the evaluation points.
        //
        //    Output, double C[O+P], the coefficients.
        //
        {
            double[] b;
            double[] dx;
            int i;
            int info = 0;
            int job;
            int n;
            double t;

            n = o + p;

            dx = new double[n];

            for (i = 0; i < n; i++)
            {
                dx[i] = x[i] - x0;
            }

            b = new double[n];

            for (i = 0; i < n; i++)
            {
                b[i] = 0.0;
            }

            b[o] = 1.0;

            job = 0;
            typeMethods.r8vm_sl(n, dx, b, job, c, ref info);

            if (info != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("DIFFER_STENCIL - Fatal error!");
                Console.WriteLine("  Vandermonde linear system is singular.");
                return;
            }

            t = typeMethods.r8_factorial(o);
            for (i = 0; i < n; i++)
            {
                c[i] = c[i] * t;
            }
        }
    }
}