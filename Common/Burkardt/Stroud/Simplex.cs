using System;
using Burkardt.Types;

namespace Burkardt.Stroud
{
    public static class Simplex
    {
        public static double simplex_nd(int setting, Func<int, int, double[], double> func, int n, ref double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIMPLEX_ND approximates an integral inside a simplex in ND.
            //
            //  Discussion:
            //
            //    An N+1 point second degree formula is used.
            //
            //    The integration region is the simplex bounded by the origin and a 
            //    convex combination of N points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    11 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Arthur Stroud,
            //    Approximate Calculation of Multiple Integrals,
            //    Prentice Hall, 1971,
            //    ISBN: 0130438936,
            //    LC: QA311.S85.
            //
            //  Parameters:
            //
            //    Input, Func < int, double[], double> func, the name of the 
            //    user supplied function which evaluates F(X).
            //
            //    Input, int N, the dimension of the space.
            //
            //    Input/output, double V[N*(N+1)].  On input, each of the
            //    N+1 columns of V contains the N coordinates of one of the
            //    "corners" of the simplex in entries 1 through N, with
            //    the last column being left free.
            //    On output, V has been overwritten in the process of
            //    computing the volume of the simplex.
            //
            //    Output, double SIMPLEX_ND, the approximate integral of the function.
            //
        {
            double c;
            int i;
            int j;
            double quad;
            double result;
            double s;
            double volume;
            double w;
            double[] x;

            x = new double[n];

            c = 1.0 / Math.Sqrt((double)(n + 2));
            w = 1.0 / (double)(n + 1);

            for (j = 0; j < n; j++)
            {
                s = 0.0;
                for (i = 0; i < n + 1; i++)
                {
                    s = s + v[i + j * n];
                }

                x[j] = w * (1.0 - c) * s;
            }

            quad = 0.0;

            for (j = 0; j < n + 1; j++)
            {
                for (i = 0; i < n; i++)
                {
                    x[i] = x[i] + c * v[i + j * n];
                }

                quad = quad + w * func(setting, n, x);
                for (i = 0; i < n; i++)
                {
                    x[i] = x[i] - c * v[i + j * n];
                }
            }

            volume = simplex_volume_nd(n, v);
            result = quad * volume;

            return result;
        }

        public static double simplex_unit_01_nd(int setting, Func<int, int, double[], double> func, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIMPLEX_UNIT_01_ND approximates an integral inside the unit simplex in ND.
            //
            //  Integration region:
            //
            //      0 <= X(1:N),
            //    and
            //      sum ( X(1:N) ) <= 1.
            //
            //  Discussion:
            //
            //    A 1 point formula of degree 1 is used.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    16 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Axel Grundmann, Michael Moeller,
            //    Invariant Integration Formulas for the N-Simplex by Combinatorial Methods,
            //    SIAM Journal on Numerical Analysis,
            //    Volume 15, Number 2, April 1978, pages 282-290.
            //
            //    Arthur Stroud,
            //    Approximate Calculation of Multiple Integrals,
            //    Prentice Hall, 1971,
            //    ISBN: 0130438936,
            //    LC: QA311.S85.
            //
            //  Parameters:
            //
            //    Input, Func < int, double[], double> func, the name of the 
            //    user supplied function to be integrated.
            //
            //    Input, int N, the dimension of the space.  
            //
            //    Output, double RESULT, the approximate integral of the function.
            //
        {
            double coef = 1.0;
            int i;
            double quad;
            double result;
            double volume;
            double[] x;

            x = new double[n];

            quad = 0.0;

            for (i = 0; i < n; i++)
            {
                x[i] = 1.0 / (double)(n);
            }

            quad = quad + coef * func(setting, n, x);

            volume = simplex_unit_volume_nd(n);
            result = quad * volume;

            return result;
        }

        public static double simplex_unit_03_nd(int setting, Func<int, int, double[], double> func, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIMPLEX_UNIT_03_ND approximates an integral inside the unit simplex in ND.
            //
            //  Integration region:
            //
            //      0 <= X(1:N),
            //    and
            //      sum ( X(1:N) ) <= 1.
            //
            //  Discussion:
            //
            //    An N+2 point formula of degree 3 is used.  This is Stroud TN:3-1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Axel Grundmann, Michael Moeller,
            //    Invariant Integration Formulas for the N-Simplex by Combinatorial Methods,
            //    SIAM Journal on Numerical Analysis,
            //    Volume 15, Number 2, April 1978, pages 282-290.
            //
            //    Arthur Stroud,
            //    Approximate Calculation of Multiple Integrals,
            //    Prentice Hall, 1971,
            //    ISBN: 0130438936,
            //    LC: QA311.S85.
            //
            //  Parameters:
            //
            //    Input, Func < int, double[], double> func, the name of the 
            //    user supplied function which is to be integrated.
            //
            //    Input, int N, the dimension of the space.  
            //
            //    Output, double SIMPLEX_UNIT_03_ND, the approximate integral of the function.
            //
        {
            double a;
            double b;
            double coef;
            int i;
            double quad;
            double result;
            double volume;
            double[] x;

            x = new double[n];

            quad = 0.0;

            for (i = 0; i < n; i++)
            {
                x[i] = 1.0 / (double)(n + 1);
            }

            coef = -0.25 * (double)((n + 1) * (n + 1)) / (double)(n + 2);
            quad = quad + coef * func(setting, n, x);

            a = 1.0 / (double)(n + 3);
            b = 3.0 / (double)(n + 3);

            for (i = 0; i < n; i++)
            {
                x[i] = a;
            }

            coef = 0.25 * (double)((n + 3) * (n + 3))
                   / (double)((n + 1) * (n + 2));
            quad = quad + coef * func(setting, n, x);

            for (i = 0; i < n; i++)
            {
                x[i] = b;
                quad = quad + coef * func(setting, n, x);
                x[i] = a;
            }

            volume = simplex_unit_volume_nd(n);
            result = quad * volume;

            return result;
        }

        public static double simplex_unit_05_nd(int setting, Func<int, int, double[], double> func, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIMPLEX_UNIT_05_ND approximates an integral inside the unit simplex in ND.
            //
            //  Integration region:
            //
            //      0 <= X(1:N),
            //    and
            //      sum ( X(1:N) ) <= 1.
            //
            //  Discussion:
            //
            //    An N^2 + 3 N + 3 point formula of degree 5 is used.  This is
            //    Stroud formula TN:5-1.
            //
            //    (For N = 2, the number of points is actually only 7, and
            //     for N = 3, the number of points is actually only 15.)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Arthur Stroud,
            //    A Fifth Degree Integration Formula for the N-Simplex,
            //    SIAM Journal on Numerical Analysis,
            //    Volume 6, Number 1, March 1969.
            //
            //  Parameters:
            //
            //    Input, Func < int, double[], double> func, the name of the 
            //    user supplied function is to be integrated.
            //
            //    Input, int N, the dimension of the space.  For this routine,
            //    it must be the case that 2 <= N <= 16.
            //
            //    Output, double SIMPLEX_UNIT_05_ND, the approximate integral of the function.
            //
        {
            double[] coef1 =
            {
                0.0,
                0.225,
                0.118518518519,
                0.0631521898883,
                0.235714285714,
                0.791575476992,
                1.85798728021,
                3.53666958042,
                5.90844340844,
                9.03765432098,
                12.9758241758,
                17.7645108738,
                23.4375030259,
                30.0224941950,
                37.5423613501,
                46.0161454949
            };
            double[] coef21 =
            {
                0.0,
                0.12593918054483,
                0.0719370837790,
                0.0470456145702,
                0.0333009774677,
                0.0248633014592,
                0.0192679696358,
                0.0153322153879,
                0.0124316229901,
                0.0102112988361,
                0.00845730697460,
                0.00703433430999,
                0.00585330520067,
                0.00485356735291,
                0.00399261092720,
                0.00323988713017
            };
            double[] coef22 =
            {
                0.0,
                0.13239415278851,
                0.0690682072263,
                0.0371530185868,
                -0.0719253160920,
                -0.264323879461,
                -0.537926779961,
                -0.886895605701,
                -1.30409181465,
                -1.78227048964,
                -2.31462336314,
                -2.89499045158,
                -3.51790849765,
                -4.17858310668,
                -4.87282884913,
                -5.59699944261
            };
            double[] coef31 =
            {
                0.0,
                0.0,
                0.0529100529100,
                0.0261368740713,
                0.0499020181331,
                0.0782233395867,
                0.109041040862,
                0.140874828568,
                0.172735353396,
                0.203992490408,
                0.234263814181,
                0.263332763315,
                0.291091849264,
                0.317504208212,
                0.342577872069,
                0.366348654344
            };
            double[] coef32 =
            {
                0.0,
                0.0,
                0.0,
                0.0254485903613,
                0.0165000982690,
                0.0115218303668,
                0.00850478779483,
                0.00655297510968,
                0.00522372456259,
                0.00428017828134,
                0.00358722367033,
                0.00306362964360,
                0.00265836687133,
                0.00233816221525,
                0.00208061510846,
                0.00187022027571
            };
            int i;
            int j;
            double quad;
            double r1;
            double r2;
            double result;
            double s1;
            double s2;
            double u1;
            double u2;
            double v1;
            double v2;
            double volume;
            double[] x;

            result = 0.0;

            if (n < 2 || 16 < n)
            {
                Console.WriteLine("");
                Console.WriteLine("SIMPLEX_UNIT_05_ND - Fatal error!");
                Console.WriteLine("  Input spatial dimension N out of range.");
                Console.WriteLine("  N = " + n + "");
                return (1);
            }

            x = new double[n];

            quad = 0.0;
            //
            //  S1
            //
            for (i = 0; i < n; i++)
            {
                x[i] = 1.0 / (double)(n + 1);
            }

            quad = quad + coef1[n - 1] * func(setting, n, x);
            //
            //  S21
            //
            r1 = ((double)(n + 4) - Math.Sqrt(15.0))
                 / (double)(n * n + 8 * n + 1);
            s1 = 1.0 - (double)(n) * r1;

            for (i = 0; i < n; i++)
            {
                x[i] = r1;
            }

            for (i = 0; i < n + 1; i++)
            {
                quad = quad + coef21[n - 1] * func(setting, n, x);

                if (0 < i)
                {
                    x[i - 1] = r1;
                }

                if (i < n)
                {
                    x[i] = s1;
                }
            }

            //
            //  S22
            //
            r2 = ((double)(n + 4) + Math.Sqrt(15.0))
                 / (double)(n * n + 8 * n + 1);
            s2 = 1.0 - (double)(n) * r2;

            for (i = 0; i < n; i++)
            {
                x[i] = r2;
            }

            for (i = 0; i < n + 1; i++)
            {
                quad = quad + coef22[n - 1] * func(setting, n, x);

                if (0 < i)
                {
                    x[i - 1] = r2;
                }

                if (i < n)
                {
                    x[i] = s2;
                }
            }

            //
            //  S31
            //
            u1 = ((double)(n + 7) + 2.0 * Math.Sqrt(15.0))
                 / (double)(n * n + 14 * n - 11);
            v1 = ((double)(4 * n - 2)
                  - (double)(n - 1) * Math.Sqrt(15.0))
                 / (double)(n * n + 14 * n - 11);

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    x[j] = u1;
                }

                x[i] = v1;

                for (j = i; j < n; j++)
                {
                    if (i < j - 1)
                    {
                        x[j - 1] = u1;
                    }

                    x[j] = v1;

                    quad = quad + coef31[n - 1] * func(setting, n, x);
                }
            }

            //
            //  S32
            //
            u2 = ((double)(n + 7) - 2.0 * Math.Sqrt(15.0))
                 / (double)(n * n + 14 * n - 11);
            v2 = ((double)(4 * n - 2)
                  + (double)(n - 1) * Math.Sqrt(15.0))
                 / (double)(n * n + 14 * n - 11);

            for (i = 0; i < n; i++)
            {

                for (j = 0; j < n; j++)
                {
                    x[j] = u2;
                }

                x[i] = v2;

                for (j = i; j < n; j++)
                {
                    if (i < j - 1)
                    {
                        x[j - 1] = u2;
                    }

                    x[j] = v2;

                    quad = quad + coef32[n - 1] * func(setting, n, x);
                }
            }

            volume = simplex_unit_volume_nd(n);
            result = quad * volume;

            return result;
        }

        public static double simplex_unit_05_2_nd(int setting, Func<int, int, double[], double> func, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIMPLEX_UNIT_05_2_ND approximates an integral in the unit simplex in ND.
            //
            //  Integration region:
            //
            //      0 <= X(1:N),
            //    and
            //      sum ( X(1:N) ) <= 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Axel Grundmann, Michael Moeller,
            //    Invariant Integration Formulas for the N-Simplex by Combinatorial Methods,
            //    SIAM Journal on Numerical Analysis,
            //    Volume 15, Number 2, April 1978, pages 282-290.
            //
            //  Parameters:
            //
            //    Input, Func < int, double[], double> func, the name of the 
            //    user supplied function to be integrated.
            //
            //    Input, int N, the dimension of the space.
            //
            //    Output, double SIMPLEX_UNIT_05_2_ND, the approximate integral of the function.
            //
        {
            double a;
            double b;
            double coef;
            int i;
            int j;
            double quad;
            double result;
            double volume;
            double[] x;

            x = new double[n];

            quad = 0.0;
            //
            //  Group 1
            //
            for (i = 0; i < n; i++)
            {
                x[i] = 1.0 / (double)(n + 1);
            }

            coef = (double)((int)Math.Pow(n + 1, 4))
                   / (double)(32 * (n + 2) * (n + 3));
            quad = quad + coef * func(setting, n, x);
            //
            //  Group 2
            //
            a = 1.0 / (double)(n + 3);
            b = 3.0 / (double)(n + 3);

            for (i = 0; i < n; i++)
            {
                x[i] = a;
            }

            coef = -(double)((int)Math.Pow(n + 3, 4))
                   / (double)(16 * (n + 1) * (n + 2) * (n + 4));
            quad = quad + coef * func(setting, n, x);

            for (i = 0; i < n; i++)
            {
                x[i] = b;
                quad = quad + coef * func(setting, n, x);
                x[i] = a;
            }

            //
            //  Group 3
            //
            a = 1.0 / (double)(n + 5);
            b = 5.0 / (double)(n + 5);

            for (i = 0; i < n; i++)
            {
                x[i] = a;
            }

            coef = (double)((int)Math.Pow(n + 5, 4))
                   / (double)(16 * (n + 1) * (n + 2) * (n + 3) * (n + 4));
            quad = quad + coef * func(setting, n, x);

            for (i = 0; i < n; i++)
            {
                x[i] = b;
                quad = quad + coef * func(setting, n, x);
                x[i] = a;
            }

            //
            //  Group 4
            //
            a = 1.0 / (double)(n + 5);
            b = 3.0 / (double)(n + 5);

            coef = (double)((int)Math.Pow(n + 5, 4))
                   / (double)(16 * (n + 1) * (n + 2) * (n + 3) * (n + 4));

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    x[j] = a;
                }

                x[i] = b;
                quad = quad + coef * func(setting, n, x);

                for (j = i + 1; j < n; j++)
                {
                    x[j] = b;
                    quad = quad + coef * func(setting, n, x);
                    x[j] = a;
                }
            }

            volume = simplex_unit_volume_nd(n);
            result = quad * volume;

            return result;
        }

        public static double simplex_unit_volume_nd(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIMPLEX_UNIT_VOLUME_ND returns the volume of the unit simplex in ND.
            //
            //  Integration region:
            //
            //      0 <= X(1:N),
            //    and
            //      sum ( X(1:N) ) <= 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 April 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the dimension of the space.
            //
            //    Output, double SIMPLEX_UNIT_VOLUME_ND, the volume of the
            //    unit simplex.
            //
        {
            double value;

            value = 1.0 / typeMethods.r8_factorial(n);

            return value;
        }

        public static double simplex_volume_nd(int n, double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIMPLEX_VOLUME_ND returns the volume of a simplex in ND.
            //
            //  Integration region:
            //
            //    The simplex bounded by the origin and a convex combination of N points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the dimension of the space.
            //
            //    Input, double V[N*(N+1)], the coordinates of the vertices.
            //
            //    Output, double SIMPLEX_VOLUME_ND, the volume of 
            //    the unit simplex.
            //
        {
            double det;
            int i;
            int info;
            int j;
            int[] pivot;
            double volume;
            double[] w;

            w = new double[n * n];

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    w[i + j * n] = v[i + (j + 1) * n] - v[i + 0 * n];
                }
            }

            pivot = new int[n];

            info = typeMethods.r8ge_fa(n, ref w, ref pivot);

            if (info != 0)
            {
                volume = 0.0;
            }
            else
            {
                det = typeMethods.r8ge_det(n, w, pivot);
                //
                //  Multiply by the volume of the unit simplex, which serves as a
                //  conversion factor between a parallelipiped and the simplex.
                //
                volume = Math.Abs(det) * simplex_unit_volume_nd(n);
            }

            return volume;
        }


    }
}