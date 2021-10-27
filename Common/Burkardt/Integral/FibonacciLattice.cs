using System;
using Burkardt.Types;

namespace Burkardt.IntegralNS
{
     public static class FibonacciLattice
     {
          public static double e_01_2d(int dim_num, double[] a, double[] b)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    E_01_2D is the exact integral of 2d test function #1.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license. 
               //
               //  Modified:
               //
               //    18 April 2003
               //
               //  Author:
               //
               //    John Burkardt
               //
               //  Reference:
               //
               //    Ian Sloan, Stephen Joe,
               //    Lattice Methods for Multiple Integration,
               //    Oxford, 1994,
               //    ISBN: 0198534728,
               //    LC: QA311.S56
               //
               //  Parameters:
               //
               //    Input, int DIM_NUM, the spatial dimension.
               //
               //    Input, double A[DIM_NUM], B[DIM_NUM], the integration limits.
               //
               //    Output, double E_01_2D, the integral of the function 
               //    over the limits.
               //
          {
               double value;

               value = 1.0;

               return value;
          }

          public static double f_01_2d(int dim_num, double[] x)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    F_01_2D is the 2D test function #1.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license. 
               //
               //  Modified:
               //
               //    19 November 2008
               //
               //  Author:
               //
               //    John Burkardt
               //
               //  Reference:
               //
               //    Ian Sloan, Stephen Joe,
               //    Lattice Methods for Multiple Integration,
               //    Oxford, 1994,
               //    ISBN: 0198534728,
               //    LC: QA311.S56
               //
               //  Parameters:
               //
               //    Input, int DIM_NUM, the spatial dimension.
               //
               //    Input, double X[DIM_NUM], the point where the function 
               //    is to be evaluated.
               //
               //    Output, double F_01_2D, the value of the function at X.
               //
          {
               double e = 2.718281828459045;
               double value;

               value = x[1] * Math.Exp(x[0] * x[1]) / (e - 2.0);

               return value;
          }

          public static double f2(double x)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    F2 evaluates a function of a scalar used in defining P2(Q).
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license. 
               //
               //  Modified:
               //
               //    19 November 2008
               //
               //  Author:
               //
               //    John Burkardt
               //
               //  Reference:
               //
               //    Ian Sloan, Stephen Joe,
               //    Lattice Methods for Multiple Integration,
               //    Oxford, 1994,
               //    ISBN: 0198534728,
               //    LC: QA311.S56
               //
               //  Parameters:
               //
               //    Input, double X, the value of the argument.
               //
               //    Output, double F2, the value of F2(X).
               //
          {
               double value;

               value = 1.0 + 2.0 * Math.PI * Math.PI * (x * x - x + 1.0 / 6.0);

               return value;
          }

          public static double f20_s(int dim_num, double[] x)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    F20_S evaluates a function of a vector used in defining P2(Q).
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license. 
               //
               //  Modified:
               //
               //    19 November 2008
               //
               //  Author:
               //
               //    John Burkardt
               //
               //  Reference:
               //
               //    Ian Sloan, Stephen Joe,
               //    Lattice Methods for Multiple Integration,
               //    Oxford, 1994,
               //    ISBN: 0198534728,
               //    LC: QA311.S56
               //
               //  Parameters:
               //
               //    Input, int DIM_NUM, the spatial dimension.
               //
               //    Input, double X[DIM_NUM], the value of the argument.
               //
               //    Output, double F20_S, the value of F20_S(X).
               //
          {
               int i;
               double value;

               value = 1.0;
               for (i = 0; i < dim_num; i++)
               {
                    value = value * (1.0 + (f2(x[i]) - 1.0));
               }

               value = value - 1.0;

               return value;
          }

          public static double fibonacci_lattice_b(int k, Func<int, double[], double> f)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    FIBONACCI_LATTICE_B applies an optimal Fibonacci lattice integration rule in 2D.
               //
               //  Discussion:
               //
               //    This routine may be applied to integrands which are not periodic.
               //
               //    When K is odd, this is the same as the symmetric Fibonacci lattice
               //    integration rule.  But when K is even, a correction is made to the
               //    corner weights which is expected to improve the results.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license. 
               //
               //  Modified:
               //
               //    19 November 2008
               //
               //  Author:
               //
               //    John Burkardt
               //
               //  Reference:
               //
               //    Ian Sloan, Stephen Joe,
               //    Lattice Methods for Multiple Integration,
               //    Oxford, 1994,
               //    ISBN: 0198534728,
               //    LC: QA311.S56
               //
               //  Parameters:
               //
               //    Input, int K, the index of the Fibonacci number to be used.
               //    K must be at least 3.
               //
               //    Input, double F ( int DIM_NUM, double X[] ), the name of the 
               //    user-supplied routine which evaluates the function.
               //
               //    Output, double FIBONACCI_LATTICE_B, the estimated integral.
               //
          {
               double delta;
               int dim;
               int dim_num = 2;
               int j;
               int m;
               int n;
               double quad;
               double quad1;
               double quad2;
               int rank;
               double[] w = new double[2 * 2];
               double[] x;
               int[] z;

               x = new double[dim_num];
               z = new int[dim_num];

               quad = 0.0;

               m = Helpers.fibonacci(k);
               n = Helpers.fibonacci(k - 1);
               //
               //  Get the corner weights.
               //
               if ((k % 2) == 1)
               {
                    w[0 + 0 * 2] = 1.0 / (double) (4 * m);
                    w[1 + 0 * 2] = 1.0 / (double) (4 * m);
                    w[0 + 1 * 2] = 1.0 / (double) (4 * m);
                    w[1 + 1 * 2] = 1.0 / (double) (4 * m);
               }
               else
               {
                    delta = 0.0;
                    for (j = 1; j <= m - 1; j++)
                    {
                         delta = delta + (double) (j * ((j * n) % m))
                              / (double) (m * m);
                    }

                    w[0 + 0 * 2] = 0.25 - delta / (double) (m);

                    delta = 0.0;
                    for (j = 1; j <= m - 1; j++)
                    {
                         delta = delta + (double) (j * (m - ((j * n) % m)))
                              / (double) (m * m);
                    }

                    w[0 + 1 * 2] = 0.25 - delta / (double) (m);

                    w[1 + 0 * 2] = w[0 + 1 * 2];
                    w[1 + 1 * 2] = w[0 + 0 * 2];
               }

               //
               //  Get all the corner values.
               //
               rank = 0;
               quad1 = 0.0;

               for (;;)
               {
                    BTuple.tuple_next(0, 1, dim_num, ref rank, ref z);

                    if (rank == 0)
                    {
                         break;
                    }

                    for (dim = 0; dim < dim_num; dim++)
                    {
                         x[dim] = (double) z[dim];
                    }

                    quad1 = quad1 + w[z[0] + z[1] * 2] * f(dim_num, x);
               }

               //
               //  Get the interior values.
               //
               z[0] = 1;
               z[1] = Helpers.fibonacci(k - 1);

               quad2 = 0.0;
               for (j = 1; j <= m - 1; j++)
               {
                    for (dim = 0; dim < dim_num; dim++)
                    {
                         x[dim] = ((double) (j * z[dim]) / (double) (m) % 1.0);
                    }

                    quad2 = quad2 + f(dim_num, x);
               }

               quad = quad1 + quad2 / (double) (m);


               return quad;
          }

          public static double fibonacci_lattice_q(int k, Func<int, double[], double> f)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    FIBONACCI_LATTICE_Q applies a Fibonacci lattice integration rule in 2D.
               //
               //  Discussion:
               //
               //    Because this is a standard lattice rule, it is really only suited
               //    for functions which are periodic, of period 1, in both X and Y.
               //
               //    The related routines FIBONACCI_LATTICE_S and FIBONACCI_LATTICE_B
               //    may be helpful in cases where the integrand does not satisfy this
               //    requirement.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license. 
               //
               //  Modified:
               //
               //    19 November 2008
               //
               //  Author:
               //
               //    John Burkardt
               //
               //  Reference:
               //
               //    Ian Sloan, Stephen Joe,
               //    Lattice Methods for Multiple Integration,
               //    Oxford, 1994,
               //    ISBN: 0198534728,
               //    LC: QA311.S56
               //
               //  Parameters:
               //
               //    Input, int K, the index of the Fibonacci number to be used.
               //    K must be at least 3.
               //
               //    Input, double F ( int DIM_NUM, double X[] ), the name of the 
               //    user-supplied routine which evaluates the function.
               //
               //    Output, double FIBONACCI_LATTICE_Q, the estimated integral.
               //
          {
               int dim;
               int dim_num = 2;
               int j;
               int m;
               double quad;
               double[] x;
               int[] z;

               x = new double[dim_num];
               z = new int[dim_num];

               quad = 0.0;

               m = Helpers.fibonacci(k);

               z[0] = 1;
               z[1] = Helpers.fibonacci(k - 1);

               for (j = 0; j <= m - 1; j++)
               {
                    for (dim = 0; dim < dim_num; dim++)
                    {
                         x[dim] = ((double) (j * z[dim]) / (double) (m) % 1.0);
                    }

                    quad = quad + f(dim_num, x);
               }

               quad = quad / (double) (m);

               return quad;
          }

          public static double[] fibonacci_lattice_q_nodes(int k)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    FIBONACCI_LATTICE_Q_NODES returns Fibonacci lattice nodes in 2D.
               //
               //  Discussion:
               //
               //    Because this is a standard lattice rule, it is really only suited
               //    for functions which are periodic, of period 1, in both X and Y.
               //
               //    The number of nodes returned is 
               //
               //      M = fibonacci ( k ).
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license. 
               //
               //  Modified:
               //
               //    19 November 2008
               //
               //  Author:
               //
               //    John Burkardt
               //
               //  Reference:
               //
               //    Ian Sloan, Stephen Joe,
               //    Lattice Methods for Multiple Integration,
               //    Oxford, 1994,
               //    ISBN: 0198534728,
               //    LC: QA311.S56
               //
               //  Parameters:
               //
               //    Input, int K, the index of the Fibonacci number to be used.
               //    K must be at least 3.
               //
               //    Output, double X[2*M], the nodes.
               //
          {
               int dim;
               int dim_num = 2;
               int j;
               int m;
               double[] x;
               int[] z;

               m = Helpers.fibonacci(k);

               x = new double[2 * m];
               z = new int[dim_num];

               z[0] = 1;
               z[1] = Helpers.fibonacci(k - 1);

               for (j = 0; j <= m - 1; j++)
               {
                    for (dim = 0; dim < dim_num; dim++)
                    {
                         x[dim + j * dim_num] = ((double) (j * z[dim])
                              / (double) (m) % 1.0);
                    }
               }

               return x;
          }

          public static double fibonacci_lattice_q1(int k, Func<int, double[], double> f)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    FIBONACCI_LATTICE_Q1 applies a Fibonacci lattice integration rule in 2D.
               //
               //  Discussion:
               //
               //    This is a modification of the algorithm in FIBONACCI_LATTICE_Q.
               //    It uses a nonlinear transformation on the integrand, which makes
               //    the lattice rule more suitable for nonperiodic integrands.
               //
               //    The transformation replaces the integration variable X by
               //
               //      PHI(X) = 3*X^2 - 2*X**3
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license. 
               //
               //  Modified:
               //
               //    21 November 2008
               //
               //  Author:
               //
               //    John Burkardt
               //
               //  Reference:
               //
               //    Ian Sloan, Stephen Joe,
               //    Lattice Methods for Multiple Integration,
               //    Oxford, 1994,
               //    ISBN: 0198534728,
               //    LC: QA311.S56
               //
               //  Parameters:
               //
               //    Input, int K, the index of the Fibonacci number to be used.
               //    K must be at least 3.
               //
               //    Input, double F ( int DIM_NUM, double X[] ), the name of the 
               //    user-supplied routine which evaluates the function.
               //
               //    Output, double FIBONACCI_LATTICE_Q1, the estimated integral.
               //
          {
               int dim_num = 2;
               double dphi;
               int i;
               int j;
               int m;
               double quad;
               double[] x;
               double[] y;
               int[] z;

               x = new double[dim_num];
               y = new double[dim_num];
               z = new int[dim_num];

               quad = 0.0;

               m = Helpers.fibonacci(k);

               z[0] = 1;
               z[1] = Helpers.fibonacci(k - 1);

               for (j = 0; j <= m - 1; j++)
               {
                    for (i = 0; i < dim_num; i++)
                    {
                         x[i] = ((double) (j * z[i]) / (double) (m) % 1.0);
                    }

                    dphi = 1.0;
                    for (i = 0; i < dim_num; i++)
                    {
                         y[i] = (3.0 - 2.0 * x[i]) * x[i] * x[i];
                         dphi = dphi * 6.0 * (1.0 - x[i]) * x[i];
                    }

                    quad = quad + f(dim_num, y) * dphi;
               }

               quad = quad / (double) (m);

               return quad;
          }

          public static double fibonacci_lattice_q2(int k, Func<int, double[], double> f)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    FIBONACCI_LATTICE_Q2 applies a Fibonacci lattice integration rule in 2D.
               //
               //  Discussion:
               //
               //    This is a modification of the algorithm in FIBONACCI_LATTICE_Q.
               //    It uses a nonlinear transformation on the integrand, which makes
               //    the lattice rule more suitable for nonperiodic integrands.
               //
               //    The transformation replaces the integration variable X by
               //
               //      PHI(X) = 3*X^3 *( 10 - 15 * X + 6 * X^2 )
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license. 
               //
               //  Modified:
               //
               //    21 November 2008
               //
               //  Author:
               //
               //    John Burkardt
               //
               //  Reference:
               //
               //    Ian Sloan, Stephen Joe,
               //    Lattice Methods for Multiple Integration,
               //    Oxford, 1994,
               //    ISBN: 0198534728,
               //    LC: QA311.S56
               //
               //  Parameters:
               //
               //    Input, int K, the index of the Fibonacci number to be used.
               //    K must be at least 3.
               //
               //    Input, double F ( int DIM_NUM, double X[] ), the name of the 
               //    user-supplied routine which evaluates the function.
               //
               //    Output, double FIBONACCI_LATTICE_Q2, the estimated integral.
               //
          {
               int dim_num = 2;
               double dphi;
               int i;
               int j;
               int m;
               double quad;
               double[] x;
               double[] y;
               int[] z;

               x = new double[dim_num];
               y = new double[dim_num];
               z = new int[dim_num];

               quad = 0.0;

               m = Helpers.fibonacci(k);

               z[0] = 1;
               z[1] = Helpers.fibonacci(k - 1);

               for (j = 0; j <= m - 1; j++)
               {
                    for (i = 0; i < dim_num; i++)
                    {
                         x[i] = ((double) (j * z[i]) / (double) (m) % 1.0);
                    }

                    dphi = 1.0;
                    for (i = 0; i < dim_num; i++)
                    {
                         y[i] = (10.0 - 15.0 * x[i] + 6.0 * Math.Pow(x[i], 2)) * Math.Pow(x[i], 3);
                         dphi = dphi * 30.0 * Math.Pow(1.0 - x[i], 2) * Math.Pow(x[i], 2);
                    }

                    quad = quad + f(dim_num, y) * dphi;
               }

               quad = quad / (double) (m);

               return quad;
          }

          public static double fibonacci_lattice_q3(int k, Func<int, double[], double> f)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    FIBONACCI_LATTICE_Q3 applies a Fibonacci lattice integration rule in 2D.
               //
               //  Discussion:
               //
               //    This is a modification of the algorithm in FIBONACCI_LATTICE_Q.
               //    It uses a nonlinear transformation on the integrand, which makes
               //    the lattice rule more suitable for nonperiodic integrands.
               //
               //    The transformation replaces the integration variable X by
               //
               //      PHI(X) = X - sin ( 2 * PI * X ) / ( 2 * PI )
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license. 
               //
               //  Modified:
               //
               //    21 November 2008
               //
               //  Author:
               //
               //    John Burkardt
               //
               //  Reference:
               //
               //    Ian Sloan, Stephen Joe,
               //    Lattice Methods for Multiple Integration,
               //    Oxford, 1994,
               //    ISBN: 0198534728,
               //    LC: QA311.S56
               //
               //  Parameters:
               //
               //    Input, int K, the index of the Fibonacci number to be used.
               //    K must be at least 3.
               //
               //    Input, double F ( int DIM_NUM, double X[] ), the name of the 
               //    user-supplied routine which evaluates the function.
               //
               //    Output, double FIBONACCI_LATTICE_Q3, the estimated integral.
               //
          {
               int dim_num = 2;
               double dphi;
               int i;
               int j;
               int m;
               double quad;
               double two_pi;
               double[] x;
               double[] y;
               int[] z;

               x = new double[dim_num];
               y = new double[dim_num];
               z = new int[dim_num];

               quad = 0.0;

               m = Helpers.fibonacci(k);

               z[0] = 1;
               z[1] = Helpers.fibonacci(k - 1);

               two_pi = 2.0 * Math.PI;

               for (j = 0; j <= m - 1; j++)
               {
                    for (i = 0; i < dim_num; i++)
                    {
                         x[i] = ((double) (j * z[i]) / (double) (m) % 1.0);
                    }

                    dphi = 1.0;
                    for (i = 0; i < dim_num; i++)
                    {
                         y[i] = x[i] - Math.Sin(two_pi * x[i]) / two_pi;
                         dphi = dphi * (1.0 - Math.Cos(two_pi * x[i]));
                    }

                    quad = quad + f(dim_num, y) * dphi;
               }

               quad = quad / (double) (m);

               return quad;
          }

          public static double fibonacci_lattice_t(int k, Func<int, double[], double> f)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    FIBONACCI_LATTICE_T applies a symmetric Fibonacci lattice integration rule in 2D.
               //
               //  Discussion:
               //
               //    This routine may be applied to integrands which are not periodic.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license. 
               //
               //  Modified:
               //
               //    21 November 2008
               //
               //  Author:
               //
               //    John Burkardt
               //
               //  Reference:
               //
               //    Ian Sloan, Stephen Joe,
               //    Lattice Methods for Multiple Integration,
               //    Oxford, 1994,
               //    ISBN: 0198534728,
               //    LC: QA311.S56
               //
               //  Parameters:
               //
               //    Input, int K, the index of the Fibonacci number to be used.
               //    K must be at least 3.
               //
               //    Input, double F ( int DIM_NUM, double X[] ), the name of the 
               //    user-supplied routine which evaluates the function.
               //
               //    Output, double FIBONACCI_LATTICE_T, the estimated integral.
               //
          {
               int dim_num = 2;
               int i;
               int j;
               int m;
               double quad;
               double quad1;
               double quad2;
               int rank;
               double w;
               double[] x;
               int[] z;

               x = new double[dim_num];
               z = new int[dim_num];

               quad = 0.0;

               m = Helpers.fibonacci(k);
               //
               //  Get all the corner values.
               //
               rank = 0;
               quad1 = 0.0;
               w = 1.0 / (double) (int)Math.Pow(2, dim_num);

               for (;;)
               {
                    BTuple.tuple_next(0, 1, dim_num, ref rank, ref z);

                    if (rank == 0)
                    {
                         break;
                    }

                    for (i = 0; i < dim_num; i++)
                    {
                         x[i] = (double) z[i];
                    }

                    quad1 = quad1 + w * f(dim_num, x);
               }

               //
               //  Get the interior values.
               //
               z[0] = 1;
               z[1] = Helpers.fibonacci(k - 1);

               quad2 = 0.0;
               for (j = 1; j <= m - 1; j++)
               {
                    for (i = 0; i < dim_num; i++)
                    {
                         x[i] = ((double) (j * z[i]) / (double) (m) % 1.0);
                    }

                    quad2 = quad2 + f(dim_num, x);
               }

               quad = (quad1 + quad2) / (double) (m);

               return quad;
          }

          public static int[] find_z20(int dim_num, int m)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    FIND_Z20 finds the appropriate Z vector to minimize P2(QS).
               //
               //  Discussion:
               //
               //    For the method of good lattice points, a number of points M, and
               //    a single generator vector Z is chosen.  The integrand is assumed
               //    to be periodic of period 1 in each argument, and is evaluated at
               //    each of the points X(I)(1:DIM_NUM) = I * Z(1:DIM_NUM) / M, 
               //    for I = 0 to M-1.  The integral is then approximated by the average
               //    of these values.
               //
               //    Assuming that DIM_NUM and M are known, and that the integrand is not
               //    known beforehand, the accuracy of the method depends entirely
               //    on the choice of Z.  One method of choosing Z is to search for
               //    the Z among all candidates which minimizes a particular quantity
               //    P_ALPHA(QS).
               //
               //    We only need to look at vectors Z of the form
               //    (1, L, L^2, ..., L^(DIM_NUM-1)),
               //    for L = 1 to M/2.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license. 
               //
               //  Modified:
               //
               //    23 November 2008
               //
               //  Author:
               //
               //    John Burkardt
               //
               //  Reference:
               //
               //    Ian Sloan, Stephen Joe,
               //    Lattice Methods for Multiple Integration,
               //    Oxford, 1994,
               //    ISBN: 0198534728,
               //    LC: QA311.S56
               //
               //  Parameters:
               //
               //    Input, int DIM_NUM, the spatial dimension.
               //
               //    Input, int M, the number of points to be used.
               //
               //    Output, int FIND_Z20[DIM_NUM], the optimal vector.
               //
          {
               int dim;
               int i;
               double q0;
               double q0_min;
               int value;
               int[] z;
               int[] z_min;

               z = new int[dim_num];
               z_min = new int[dim_num];

               q0_min = typeMethods.r8_huge();

               for (i = 1; i <= m / 2; i++)
               {
                    value = 1;
                    for (dim = 0; dim < dim_num; dim++)
                    {
                         z[dim] = value;
                         value = (value * i) % m;
                    }

                    //
                    //  Use this Z and the lattice integral method Q0 of order M,
                    //  to approximate the integral of P2.
                    //
                    q0 = Lattice.lattice(dim_num, m, z, f20_s);
                    //
                    //  If this result is the smallest so far, save the corresponding Z.
                    //
                    if (q0 < q0_min)
                    {
                         q0_min = q0;
                         for (dim = 0; dim < dim_num; dim++)
                         {
                              z_min[dim] = z[dim];
                         }
                    }
               }

               //
               //  Return the best Z.
               //
               return z_min;
          }
     }
}