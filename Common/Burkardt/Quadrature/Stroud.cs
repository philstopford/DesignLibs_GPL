using System;
using Burkardt.IntegralNS;
using Burkardt.Types;

namespace Burkardt.Quadrature
{
     public static class Stroud
     {

          public static void cn_leg_03_1(int n, int o, ref double[] x, ref double[] w)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    CN_LEG_03_1 implements the Stroud rule CN:3-1 for region CN_LEG.
               //
               //  Discussion:
               //
               //    The rule has order 
               //
               //      O = 2 * N.
               //
               //    The rule has precision P = 3.
               //
               //    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
               //
               //      w(x) = 1.
               //
               //    The necessary treatment of the final coordinate of points when
               //    N is odd seems to vary from what Stroud declares! 
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    03 March 2010
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
               //    Input, int N, the spatial dimension.
               //
               //    Input, int O, the order.
               //
               //    Output, double X[N*O], the abscissas.
               //
               //    Output, double W[O], the weights.
               //
          {
               double arg;
               int expon;
               int i;
               int j;
               double pi = 3.141592653589793;
               int r;
               double volume;

               expon = 0;
               volume = C1.c1_leg_monomial_integral(expon);
               volume = Math.Pow(volume, n);

               for (j = 0; j < o; j++)
               {
                    i = 0;
                    for (r = 1; r <= (n / 2); r++)
                    {
                         arg = (double)((2 * r - 1) * (j + 1)) * pi / (double)(n);

                         x[i + j * n] = Math.Sqrt(2.0) * Math.Cos(arg) / Math.Sqrt(3.0);
                         i = i + 1;
                         x[i + j * n] = Math.Sqrt(2.0) * Math.Sin(arg) / Math.Sqrt(3.0);
                         i = i + 1;
                    }

                    //
                    //  The following code does not correspond to what Stroud declares.
                    //
                    if (i < n)
                    {
                         if (n == 1)
                         {
                              x[i + j * n] = typeMethods.r8_mop(j + 1) / Math.Sqrt(3.0);
                         }
                         else
                         {
                              x[i + j * n] = Math.Sqrt(2.0) * typeMethods.r8_mop(j + 1) / Math.Sqrt(3.0);
                         }

                         i = i + 1;
                    }
               }

               for (j = 0; j < o; j++)
               {
                    w[j] = volume / (double)(o);
               }

               return;
          }

          public static int cn_leg_03_1_size(int n)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    CN_LEG_03_1_SIZE sizes the Stroud rule CN:3-1 for region CN_LEG.
               //
               //  Discussion:
               //
               //    The rule has order 
               //
               //      O = 2 * N.
               //
               //    The rule has precision P = 3.
               //
               //    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
               //
               //      w(x) = 1.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    03 March 2010
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
               //    Input, int N, the spatial dimension.
               //
               //    Output, int CN_LEG_03_1_SIZE, the order.
               //
          {
               int o;

               o = 2 * n;

               return o;
          }


          public static void cn_leg_05_1(int n, int option, int o, ref double[] x, ref double[] w)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    CN_LEG_05_1 implements the Stroud rule CN:5-1 for region CN_LEG.
               //
               //  Discussion:
               //
               //    The rule has order 
               //
               //      O = N^2 + N + 2.
               //
               //    The rule has precision P = 5.
               //
               //    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
               //
               //      w(x) = 1.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    03 March 2010
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
               //    Input, int N, the spatial dimension.
               //    N must be 4, 5, or 6.
               //
               //    Input, int OPTION, is only used in case N = 4 or 5.
               //    In that case, OPTION should be 1 or 2 to select the
               //    two available variants of the rule.
               //
               //    Input, int O, the order.
               //
               //    Output, double X[N*O], the abscissas.
               //
               //    Output, double W[O], the weights.
               //
          {
               double a = 0;
               double b = 0;
               double c = 0;
               double eta = 0;
               int expon = 0;
               double gamma = 0;
               int i = 0;
               int i1 = 0;
               int i2 = 0;
               int k = 0;
               double lambda = 0;
               double mu = 0;
               double volume = 0;
               double xsi = 0;

               if (n < 4 || 6 < n)
               {
                    Console.WriteLine("");
                    Console.WriteLine("CN_LEG_05_1 - Fatal error!");
                    Console.WriteLine("  The value of N must be 4, 5, or 6.");
                    return;
               }

               if (n == 4 || n == 5)
               {
                    if (option < 1 || 2 < option)
                    {
                         Console.WriteLine("");
                         Console.WriteLine("CN_LEG_05_1 - Fatal error!");
                         Console.WriteLine("  When N = 4 or 5, the value of OPTION must be 1 or 2.");
                         return;
                    }
               }

               expon = 0;
               volume = C1.c1_leg_monomial_integral(expon);
               volume = Math.Pow(volume, n);

               if (n == 4 && option == 1)
               {
                    eta = 0.778984505799815E+00;
                    lambda = 1.284565137874656E+00;
                    xsi = -0.713647298819253E+00;
                    mu = -0.715669761974162E+00;
                    gamma = 0.217089151000943E+00;
                    a = 0.206186096875899E-01 * volume;
                    b = 0.975705820221664E-02 * volume;
                    c = 0.733921929172573E-01 * volume;
               }
               else if (n == 4 && option == 2)
               {
                    eta = 0.546190755827425E+00;
                    lambda = 0.745069130115661E+00;
                    xsi = -0.413927294508700E+00;
                    mu = -0.343989637454535E+00;
                    gamma = 1.134017894600344E+00;
                    a = 0.853094758323323E-01 * volume;
                    b = 0.862099000096395E-01 * volume;
                    c = 0.116418206881849E-01 * volume;
               }
               else if (n == 5 && option == 1)
               {
                    eta = 0.522478547481276E+00;
                    lambda = 0.936135175985774E+00;
                    xsi = -0.246351362101519E+00;
                    mu = -0.496308106093758E+00;
                    gamma = 0.827180176822930E+00;
                    a = 0.631976901960153E-01 * volume;
                    b = 0.511464127430166E-01 * volume;
                    c = 0.181070246088902E-01 * volume;
               }
               else if (n == 5 && option == 2)
               {
                    eta = 0.798317301388741E+00;
                    lambda = 0.637344273885728E+00;
                    xsi = -0.455245909918377E+00;
                    mu = -1.063446229997311E+00;
                    gamma = 0.354482076665770E+00;
                    a = 0.116952384292206E-01 * volume;
                    b = 0.701731258612708E-01 * volume;
                    c = 0.137439132264426E-01 * volume;
               }
               else if (n == 6)
               {
                    eta = 0.660225291773525E+00;
                    lambda = 1.064581294844754E+00;
                    xsi = 0.000000000000000E+00;
                    mu = -0.660225291773525E+00;
                    gamma = 0.660225291773525E+00;
                    a = 0.182742214532872E-01 * volume;
                    b = 0.346020761245675E-01 * volume;
                    c = 0.182742214532872E-01 * volume;
               }

               k = -1;

               k = k + 1;
               for (i = 0; i < n; i++)
               {
                    x[i + k * n] = eta;
               }

               w[k] = a;

               k = k + 1;
               for (i = 0; i < n; i++)
               {
                    x[i + k * n] = -eta;
               }

               w[k] = a;

               for (i1 = 0; i1 < n; i1++)
               {
                    k = k + 1;
                    for (i = 0; i < n; i++)
                    {
                         x[i + k * n] = xsi;
                    }

                    x[i1 + k * n] = lambda;
                    w[k] = b;
               }

               for (i1 = 0; i1 < n; i1++)
               {
                    k = k + 1;
                    for (i = 0; i < n; i++)
                    {
                         x[i + k * n] = -xsi;
                    }

                    x[i1 + k * n] = -lambda;
                    w[k] = b;
               }

               for (i1 = 0; i1 < n - 1; i1++)
               {
                    for (i2 = i1 + 1; i2 < n; i2++)
                    {
                         k = k + 1;
                         for (i = 0; i < n; i++)
                         {
                              x[i + k * n] = gamma;
                         }

                         x[i1 + k * n] = mu;
                         x[i2 + k * n] = mu;
                         w[k] = c;
                    }
               }

               for (i1 = 0; i1 < n - 1; i1++)
               {
                    for (i2 = i1 + 1; i2 < n; i2++)
                    {
                         k = k + 1;
                         for (i = 0; i < n; i++)
                         {
                              x[i + k * n] = -gamma;
                         }

                         x[i1 + k * n] = -mu;
                         x[i2 + k * n] = -mu;
                         w[k] = c;
                    }
               }
          }

          public static int cn_leg_05_1_size(int n)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    CN_LEG_05_1_SIZE sizes the Stroud rule CN:5-1 for region CN_LEG.
               //
               //  Discussion:
               //
               //    The rule has order 
               //
               //      O = N^2 + N + 2.
               //
               //    The rule has precision P = 5.
               //
               //    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
               //
               //      w(x) = 1.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    03 March 2010
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
               //    Input, int N, the spatial dimension.
               //
               //    Output, int CN_LEG_05_1_SIZE, the order.
               //
          {
               int o;

               o = n * n + n + 2;

               return o;
          }

          public static void cn_leg_05_2(int n, int o, ref double[] x, ref double[] w)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    CN_LEG_05_2 implements the Stroud rule CN:5-2 for region CN_LEG.
               //
               //  Discussion:
               //
               //    The rule has order 
               //
               //      O = 2 N^2 + 1.
               //
               //    The rule has precision P = 5.
               //
               //    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
               //
               //      w(x) = 1.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    03 March 2010
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
               //    Input, int N, the spatial dimension.
               //    N must be at least 2.
               //
               //    Input, int O, the order.
               //
               //    Output, double X[N*O], the abscissas.
               //
               //    Output, double W[O], the weights.
               //
          {
               double b0;
               double b1;
               double b2;
               int expon;
               int i;
               int i1;
               int i2;
               int k;
               double r;
               double volume;

               if (n < 2)
               {
                    Console.WriteLine("");
                    Console.WriteLine("CN_LEG_05_2 - Fatal error!");
                    Console.WriteLine("  N must be at least 2.");
                    return;
               }

               expon = 0;
               volume = C1.c1_leg_monomial_integral(expon);
               volume = Math.Pow(volume, n);

               b0 = (double)(25 * n * n - 115 * n + 162) * volume / 162.0;
               b1 = (double)(70 - 25 * n) * volume / 162.0;
               b2 = 25.0 * volume / 324.0;

               r = Math.Sqrt(3.0 / 5.0);

               k = -1;

               k = k + 1;
               for (i = 0; i < n; i++)
               {
                    x[i + k * n] = 0.0;
               }

               w[k] = b0;

               for (i1 = 0; i1 < n; i1++)
               {
                    k = k + 1;
                    for (i = 0; i < n; i++)
                    {
                         x[i + k * n] = 0.0;
                    }

                    x[i1 + k * n] = +r;
                    w[k] = b1;

                    k = k + 1;
                    for (i = 0; i < n; i++)
                    {
                         x[i + k * n] = 0.0;
                    }

                    x[i1 + k * n] = -r;
                    w[k] = b1;
               }

               for (i1 = 0; i1 < n - 1; i1++)
               {
                    for (i2 = i1 + 1; i2 < n; i2++)
                    {
                         k = k + 1;
                         for (i = 0; i < n; i++)
                         {
                              x[i + k * n] = 0.0;
                         }

                         x[i1 + k * n] = +r;
                         x[i2 + k * n] = +r;
                         w[k] = b2;

                         k = k + 1;
                         for (i = 0; i < n; i++)
                         {
                              x[i + k * n] = 0.0;
                         }

                         x[i1 + k * n] = +r;
                         x[i2 + k * n] = -r;
                         w[k] = b2;

                         k = k + 1;
                         for (i = 0; i < n; i++)
                         {
                              x[i + k * n] = 0.0;
                         }

                         x[i1 + k * n] = -r;
                         x[i2 + k * n] = +r;
                         w[k] = b2;

                         k = k + 1;
                         for (i = 0; i < n; i++)
                         {
                              x[i + k * n] = 0.0;
                         }

                         x[i1 + k * n] = -r;
                         x[i2 + k * n] = -r;
                         w[k] = b2;
                    }
               }

               return;
          }

          public static int cn_leg_05_2_size(int n)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    CN_LEG_05_2_SIZE sizes the Stroud rule CN:5-2 for region CN_LEG.
               //
               //  Discussion:
               //
               //    The rule has order 
               //
               //      O = 2 N^2 + 1.
               //
               //    The rule has precision P = 5.
               //
               //    CN_LEG is the cube [-1,+1]^N with the Legendre weight function
               //
               //      w(x) = 1.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    03 March 2010
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
               //    Input, int N, the spatial dimension.
               //
               //    Output, int CN_LEG_05_2_SIZE, the order.
               //
          {
               int o;

               o = 2 * n * n + 1;

               return o;
          }

          public static void cn_geg_00_1(int n, double alpha, int o, ref double[] x, ref double[] w)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    CN_GEG_00_1 implements the midpoint rule for region CN_GEG.
               //
               //  Discussion:
               //
               //    The rule has order O = 1.
               //
               //    The rule has precision P = 0.
               //
               //    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
               //
               //      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
               //
               //    with -1.0 < alpha.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    30 January 2010
               //
               //  Author:
               //
               //    John Burkardt
               //
               //  Parameters:
               //
               //    Input, int N, the spatial dimension.
               //
               //    Input, double ALPHA, the parameter.
               //    -1.0 < ALPHA.
               //
               //    Input, int O, the order.
               //
               //    Output, double X[N*O], the abscissas.
               //
               //    Output, double W[O], the weights.
               //
          {
               int expon;
               int k;
               double volume;

               if (alpha <= -1.0)
               {
                    Console.WriteLine("");
                    Console.WriteLine("CN_GEG_00_1 - Fatal error!");
                    Console.WriteLine("  ALPHA <= -1.0");
                    return;
               }

               expon = 0;
               volume = C1.c1_geg_monomial_integral(alpha, expon);
               volume = Math.Pow(volume, n);

               typeMethods.r8vec_zero(n * o, ref x);

               k = -1;
               //
               //  1 point.
               //
               k = k + 1;
               w[k] = volume;

               return;
          }

          public static int cn_geg_00_1_size(int n, double alpha)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    CN_GEG_00_1_SIZE sizes the midpoint rule for region CN_GEG.
               //
               //  Discussion:
               //
               //    The rule has order O = 1.
               //
               //    The rule has precision P = 0.
               //
               //    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
               //
               //      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
               //
               //    with -1.0 < alpha.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    30 January 2010
               //
               //  Author:
               //
               //    John Burkardt
               //
               //  Parameters:
               //
               //    Input, int N, the spatial dimension.
               //
               //    Input, double ALPHA, the parameter.
               //    -1.0 < ALPHA.
               //
               //    Output, int CN_GEG_00_1_SIZE, the order.
               //
          {
               int o;

               if (alpha <= -1.0)
               {
                    Console.WriteLine("");
                    Console.WriteLine("CN_GEG_00_1_SIZE - Fatal error!");
                    Console.WriteLine("  ALPHA <= -1.0");
                    return (1);
               }

               o = 1;

               return o;
          }

          public static void cn_geg_01_1(int n, double alpha, int o, ref double[] x, ref double[] w)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    CN_GEG_01_1 implements a precision 1 rule for region CN_GEG.
               //
               //  Discussion:
               //
               //    The rule has order O = 1.
               //
               //    The rule has precision P = 1.
               //
               //    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
               //
               //      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
               //
               //    with -1.0 < alpha.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    30 January 2010
               //
               //  Author:
               //
               //    John Burkardt
               //
               //  Parameters:
               //
               //    Input, int N, the spatial dimension.
               //
               //    Input, double ALPHA, the parameter.
               //    -1.0 < ALPHA.
               //
               //    Input, int O, the order.
               //
               //    Output, double X[N*O], the abscissas.
               //
               //    Output, double W[[O], the weights.
               //
          {
               int expon;
               int i;
               int k;
               double value1;
               double value2;
               double volume;

               if (alpha <= -1.0)
               {
                    Console.WriteLine("");
                    Console.WriteLine("CN_GEG_01_1 - Fatal error!");
                    Console.WriteLine("  ALPHA <= -1.0");
                    return;
               }

               expon = 0;
               value1 = C1.c1_geg_monomial_integral(alpha, expon);
               volume = Math.Pow(value1, n);

               expon = 1;
               value2 = C1.c1_geg_monomial_integral(alpha, expon);

               typeMethods.r8vec_zero(n * o, ref x);

               k = -1;
               //
               //  1 point.
               //
               k = k + 1;
               for (i = 0; i < n; i++)
               {
                    x[i + k * n] = value2 / value1;
               }

               w[k] = volume;
          }

          public static int cn_geg_01_1_size(int n, double alpha)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    CN_GEG_01_1_SIZE sizes a precision 1 rule for region CN_GEG.
               //
               //  Discussion:
               //
               //    The rule has order O = 1.
               //
               //    The rule has precision P = 1.
               //
               //    CN_GEG is the cube [-1,+1]^N with the Gegenbauer weight function
               //
               //      w(alpha;x) = product ( 1 <= i <= n ) (1-x(i)^2)^alpha.
               //
               //    with -1.0 < alpha.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    30 January 2010
               //
               //  Author:
               //
               //    John Burkardt
               //
               //  Parameters:
               //
               //    Input, int N, the spatial dimension.
               //
               //    Input, double ALPHA, the parameter.
               //    -1.0 < ALPHA.
               //
               //    Output, int CN_GEG_01_1_SIZE, the order.
               //
          {
               int o;

               if (alpha <= -1.0)
               {
                    Console.WriteLine("");
                    Console.WriteLine("CN_GEG_01_1_SIZE - Fatal error!");
                    Console.WriteLine("  ALPHA <= -1.0");
                    return (1);
               }

               o = 1;

               return o;
          }

          public static void en_her_01_1(int n, int o, ref double[] x, ref double[] w)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_HER_01_1 implements the Stroud rule 1.1 for region EN_HER.
               //
               //  Discussion:
               //
               //    The rule has order O = 1.
               //
               //    The rule has precision P = 1.
               //
               //    EN_HER is the entire N-dimensional space with weight function
               //
               //      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    23 January 2010
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
               //    Input, int N, the spatial dimension.
               //
               //    Input, int O, the order.
               //
               //    Output, double X[N*O], the abscissas.
               //
               //    Output, double W[O], the weights.
               //
          {
               int k;
               double pi = 3.141592653589793E+00;
               double volume;
               double volume_1d;

               volume_1d = Math.Sqrt(pi);
               volume = Math.Pow(volume_1d, n);

               typeMethods.r8vec_zero(n * o, ref x);

               k = 0;
               //
               //  1 point.
               //
               w[k] = volume;
          }

          public static int en_her_01_1_size(int n)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_HER_01_1_SIZE sizes the Stroud rule 1.1 for region EN_HER.
               //
               //  Discussion:
               //
               //    The rule has order O = 1.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    23 January 2010
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
               //    Input, int N, the spatial dimension.
               //
               //    Output, int EN_HER_01_1_SIZE, the order.
               //
          {
               int o;

               o = 1;

               return o;
          }
          
          public static void en_her_03_1 ( int n, int o, ref double[] x, ref double[] w )

          //****************************************************************************80
          //
          //  Purpose:
          //
          //    EN_HER_03_1 implements the Stroud rule 3.1 for region EN_HER.
          //
          //  Discussion:
          //
          //    The rule has order O = 2 * N.
          //
          //    The rule has precision P = 3.
          //
          //    EN_HER is the entire N-dimensional space with weight function
          //
          //      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
          //
          //  Licensing:
          //
          //    This code is distributed under the GNU LGPL license.
          //
          //  Modified:
          //
          //    23 January 2010
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
          //    Input, int N, the spatial dimension.
          //
          //    Input, int O, the order.
          //
          //    Output, double X[N*O], the abscissas.
          //
          //    Output, double W[O], the weights.
          //
          {
               double a;
               int i;
               int k;
               double pi = 3.141592653589793;
               double r;
               double volume;
               double volume_1d;

               volume_1d = Math.Sqrt(pi);
               volume = Math.Pow(volume_1d, n);

               a = volume / (double) (o);
               r = Math.Sqrt((double) (n) / 2.0);

               typeMethods.r8vec_zero(n * o, ref x);

               k = -1;
               //
               //  2 * N points.
               //
               for (i = 0; i < n; i++)
               {
                    k = k + 1;
                    x[i + k * n] = -r;
                    w[k] = a;
                    k = k + 1;
                    x[i + k * n] = +r;
                    w[k] = a;
               }
          }

          public static int en_her_03_1_size(int n)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_HER_03_1_SIZE sizes the Stroud rule 3.1 for region EN_HER.
               //
               //  Discussion:
               //
               //    The rule has order O = 2 * N.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    23 January 2010
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
               //    Input, int N, the spatial dimension.
               //
               //    Output, int EN_HER_03_1_SIZE, the order.
               //
          {
               int o;

               o = 2 * n;

               return o;
          }

          public static void en_her_05_1(int n, int option, int o, ref double[] x, ref double[] w)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_HER_05_1 implements the Stroud rule 5.1 for region EN_HER.
               //
               //  Discussion:
               //
               //    The rule has order O = N^2 + N + 2.
               //
               //    The rule has precision P = 5.
               //
               //    EN_HER is the entire N-dimensional space with weight function
               //
               //      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
               //
               //    For N = 3, 5 and 6, there are two versions of the rule, chosen by setting 
               //    the OPTION variable to 1 or 2.
               //
               //    Versions of this rule are only available for N = 2 through 7.
               //
               //    There is a typographical error in the reference.
               //    For the second version of the rule for N = 2, the line
               //      gamma =    0.313300683022281E+00
               //    should read
               //      gamma =    0.312200683022281E+00
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    03 March 2010
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
               //    Input, int N, the spatial dimension.
               //    2 <= N <= 7.
               //
               //    Input, int OPTION, selects option 1 or 2.
               //
               //    Input, int O, the order.
               //
               //    Output, double X[N*O], the abscissas.
               //
               //    Output, double W[O], the weights.
               //
          {
               double a = 0;
               double b = 0;
               double c = 0;
               double eta = 0;
               double gamma = 0;
               int i = 0;
               int i1 = 0;
               int j = 0;
               int k = 0;
               double lambda = 0;
               double mu = 0;
               double pi = 3.141592653589793E+00;
               double volume = 0;
               double volume_1d = 0;
               double xsi = 0;

               if (n < 2 || 7 < n)
               {
                    Console.WriteLine("");
                    Console.WriteLine("EN_HER_05_1 - Fatal error!");
                    Console.WriteLine("  2 <= N <= 7 required.");
                    return;
               }

               if (option < 1 || 2 < option)
               {
                    Console.WriteLine("");
                    Console.WriteLine("EN_HER_05_1 - Fatal error!");
                    Console.WriteLine("  1 <= OPTION <= 2 required.");
                    return;
               }

               if (option == 2)
               {
                    if (n != 3 && n != 5 && n != 6)
                    {
                         Console.WriteLine("");
                         Console.WriteLine("EN_HER_05_1 - Fatal error!");
                         Console.WriteLine("  OPTION = 2 requires N = 3, 5 or 6.");
                         return;
                    }
               }

               volume_1d = Math.Sqrt(pi);
               volume = Math.Pow(volume_1d, n);

               if (n == 2)
               {
                    eta = 0.446103183094540E+00;
                    lambda = 0.136602540378444E+01;
                    xsi = -0.366025403784439E+00;
                    mu = 0.198167882945871E+01;
                    gamma = 0.000000000000000E+00;
                    a = 0.328774019778636E+00 * volume;
                    b = 0.833333333333333E-01 * volume;
                    c = 0.455931355469736E-02 * volume;
               }
               else if (n == 3 && option == 1)
               {
                    eta = 0.476731294622796E+00;
                    lambda = 0.935429018879534E+00;
                    xsi = -0.731237647787132E+00;
                    mu = 0.433155309477649E+00;
                    gamma = 0.266922328697744E+01;
                    a = 0.242000000000000E+00 * volume;
                    b = 0.810000000000000E-01 * volume;
                    c = 0.500000000000000E-02 * volume;
               }
               //
               //  The value of gamma that follows corrects an error in the reference.
               //
               else if (n == 3 && option == 2)
               {
                    eta = 0.476731294622796E+00;
                    lambda = 0.128679320334269E+01;
                    xsi = -0.379873463323979E+00;
                    mu = -0.192386729447751E+01;
                    gamma = 0.312200683022281E+00;
                    a = 0.242000000000000E+00 * volume;
                    b = 0.810000000000000E-01 * volume;
                    c = 0.500000000000000E-02 * volume;
               }
               else if (n == 4)
               {
                    eta = 0.523945658287507E+00;
                    lambda = 0.119433782552719E+01;
                    xsi = -0.398112608509063E+00;
                    mu = -0.318569372920112E+00;
                    gamma = 0.185675837424096E+01;
                    a = 0.155502116982037E+00 * volume;
                    b = 0.777510584910183E-01 * volume;
                    c = 0.558227484231506E-02 * volume;
               }
               else if (n == 5 && option == 1)
               {
                    eta = 0.214972564378798E+01;
                    lambda = 0.464252986016289E+01;
                    xsi = -0.623201054093728E+00;
                    mu = -0.447108700673434E+00;
                    gamma = 0.812171426076311E+00;
                    a = 0.487749259189752E-03 * volume;
                    b = 0.487749259189752E-03 * volume;
                    c = 0.497073504444862E-01 * volume;
               }
               else if (n == 5 && option == 2)
               {
                    eta = 0.615369528365158E+00;
                    lambda = 0.132894698387445E+01;
                    xsi = -0.178394363877324E+00;
                    mu = -0.745963266507289E+00;
                    gamma = 0.135503972310817E+01;
                    a = 0.726415024414905E-01 * volume;
                    b = 0.726415024414905E-01 * volume;
                    c = 0.641509853510569E-02 * volume;
               }
               else if (n == 6 && option == 1)
               {
                    eta = 0.100000000000000E+01;
                    lambda = 0.141421356237309E+01;
                    xsi = 0.000000000000000E+00;
                    mu = -0.100000000000000E+01;
                    gamma = 0.100000000000000E+01;
                    a = 0.781250000000000E-02 * volume;
                    b = 0.625000000000000E-01 * volume;
                    c = 0.781250000000000E-02 * volume;
               }
               else if (n == 6 && option == 2)
               {
                    eta = 0.100000000000000E+01;
                    lambda = 0.942809041582063E+00;
                    xsi = -0.471404520791032E+00;
                    mu = -0.166666666666667E+01;
                    gamma = 0.333333333333333E+00;
                    a = 0.781250000000000E-02 * volume;
                    b = 0.625000000000000E-01 * volume;
                    c = 0.781250000000000E-02 * volume;
               }
               else if (n == 7)
               {
                    eta = 0.000000000000000E+00;
                    lambda = 0.959724318748357E+00;
                    xsi = -0.772326488820521E+00;
                    mu = -0.141214270131942E+01;
                    gamma = 0.319908106249452E+00;
                    a = 0.111111111111111E+00 * volume;
                    b = 0.138888888888889E-01 * volume;
                    c = 0.138888888888889E-01 * volume;
               }

               typeMethods.r8vec_zero(n * o, ref x);

               k = -1;
               //
               //  2 points.
               //
               k = k + 1;
               for (i1 = 0; i1 < n; i1++)
               {
                    x[i1 + k * n] = -eta;
               }

               w[k] = a;
               k = k + 1;
               for (i1 = 0; i1 < n; i1++)
               {
                    x[i1 + k * n] = +eta;
               }

               w[k] = a;
               //
               //  2 * N points.
               //
               for (i = 0; i < n; i++)
               {
                    k = k + 1;
                    for (i1 = 0; i1 < n; i1++)
                    {
                         x[i1 + k * n] = -xsi;
                    }

                    x[i + k * n] = -lambda;
                    w[k] = b;
                    k = k + 1;
                    for (i1 = 0; i1 < n; i1++)
                    {
                         x[i1 + k * n] = +xsi;
                    }

                    x[i + k * n] = +lambda;
                    w[k] = b;
               }

               //
               //  2 * ( N * ( N - 1 ) / 2 ) points.
               //
               for (i = 0; i < n - 1; i++)
               {
                    for (j = i + 1; j < n; j++)
                    {
                         k = k + 1;
                         for (i1 = 0; i1 < n; i1++)
                         {
                              x[i1 + k * n] = -gamma;
                         }

                         x[i + k * n] = -mu;
                         x[j + k * n] = -mu;
                         w[k] = c;
                         k = k + 1;
                         for (i1 = 0; i1 < n; i1++)
                         {
                              x[i1 + k * n] = +gamma;
                         }

                         x[i + k * n] = +mu;
                         x[j + k * n] = +mu;
                         w[k] = c;
                    }
               }
          }

          public static int en_her_05_1_size(int n)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_HER_05_1_SIZE sizes the Stroud rule 5.1 for region EN_HER.
               //
               //  Discussion:
               //
               //    The rule has order O = N^2 + N + 2.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    03 March 2010
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
               //    Input, int N, the spatial dimension.
               //
               //    Output, int EN_HER_05_1_SIZE, the order.
               //
          {
               int o;

               o = n * n + n + 2;

               return o;
          }

          public static void en_her_05_2(int n, int o, ref double[] x, ref double[] w)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_HER_05_2 implements the Stroud rule 5.2 for region EN_HER.
               //
               //  Discussion:
               //
               //    The rule has order O = 2 * N^2 + 1.
               //
               //    The rule has precision P = 5.
               //
               //    EN_HER is the entire N-dimensional space with weight function
               //
               //      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 ) 
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    24 January 2010
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
               //    Input, int N, the spatial dimension.
               //
               //    Input, int O, the order.
               //
               //    Output, double X[N*O], the abscissas.
               //
               //    Output, double W[O], the weights.
               //
          {
               double a;
               double b;
               double c;
               int i;
               int j;
               int k;
               double pi = 3.141592653589793;
               double r;
               double s;
               double volume;
               double volume_1d;

               volume_1d = Math.Sqrt(pi);
               volume = Math.Pow(volume_1d, n);

               a = 2.0E+00 * volume / (double) (n + 2);
               b = (double) (4 - n) * volume / 2.0
                                             / (double) ((n + 2) * (n + 2));
               c = volume / (double) ((n + 2) * (n + 2));

               r = Math.Sqrt((double) (n + 2) / 2.0);
               s = Math.Sqrt((double) (n + 2) / 4.0);

               typeMethods.r8vec_zero(n * o, ref x);

               k = -1;
               //
               //  1 point.
               //
               k = k + 1;
               w[k] = a;
               //
               //  2 * N points.
               //
               for (i = 0; i < n; i++)
               {
                    k = k + 1;
                    x[i + k * n] = -r;
                    w[k] = b;
                    k = k + 1;
                    x[i + k * n] = +r;
                    w[k] = b;
               }

               //
               //  4 * ( N * ( N - 1 ) / 2 ) points.
               //
               for (i = 0; i < n - 1; i++)
               {
                    for (j = i + 1; j < n; j++)
                    {
                         k = k + 1;
                         x[i + k * n] = -s;
                         x[j + k * n] = -s;
                         w[k] = c;
                         k = k + 1;
                         x[i + k * n] = -s;
                         x[j + k * n] = +s;
                         w[k] = c;
                         k = k + 1;
                         x[i + k * n] = +s;
                         x[j + k * n] = -s;
                         w[k] = c;
                         k = k + 1;
                         x[i + k * n] = +s;
                         x[j + k * n] = +s;
                         w[k] = c;
                    }
               }
          }

          public static int en_her_05_2_size(int n)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_HER_05_2_SIZE sizes the Stroud rule 5.2 for region EN_HER.
               //
               //  Discussion:
               //
               //    The rule has order O = 2 * N^2 + 1.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    23 January 2010
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
               //    Input, int N, the spatial dimension.
               //
               //    Output, int EN_HER_01_1_SIZE, the order.
               //
          {
               int o;

               o = 2 * n * n + 1;

               return o;
          }


          public static void en_r2_01_1(int n, int o, ref double[] x, ref double[] w)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_R2_01_1 implements the Stroud rule 1.1 for region EN_R2.
               //
               //  Discussion:
               //
               //    The rule has order O = 1.
               //
               //    The rule has precision P = 1.
               //
               //    EN_R2 is the entire N-dimensional space with weight function
               //
               //      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    23 January 2010
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
               //    Input, int N, the spatial dimension.
               //
               //    Input, int O, the order.
               //
               //    Output, double X[N*O], the abscissas.
               //
               //    Output, double W[O], the weights.
               //
          {
               int k;
               double pi = 3.141592653589793E+00;
               double volume;

               volume = Math.Sqrt(Math.Pow(pi, n));

               typeMethods.r8vec_zero(n * o, ref x);

               k = 0;
               //
               //  1 point.
               //
               w[k] = volume;

               return;
          }

          public static int en_r2_01_1_size(int n)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_R2_01_1_SIZE sizes the Stroud rule 1.1 for region EN_R2.
               //
               //  Discussion:
               //
               //    The rule has order O = 1.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    23 January 2010
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
               //    Input, int N, the spatial dimension.
               //
               //    Output, int EN_R2_01_1_SIZE, the order.
               //
          {
               int o;

               o = 1;

               return o;
          }


          public static void en_r2_03_1(int n, int o, ref double[] x, ref double[] w)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_R2_03_1 implements the Stroud rule 3.1 for region EN_R2.
               //
               //  Discussion:
               //
               //    The rule has order O = 2 * N.
               //
               //    The rule has precision P = 3.
               //
               //    EN_R2 is the entire N-dimensional space with weight function
               //
               //      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    23 January 2010
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
               //    Input, int N, the spatial dimension.
               //
               //    Input, int O, the order.
               //
               //    Output, double X[N*O], the abscissas.
               //
               //    Output, double W[O], the weights.
               //
          {
               double a;
               int i;
               int k;
               double pi = 3.141592653589793;
               double r;
               double volume;

               volume = Math.Sqrt(Math.Pow(pi, n));

               a = volume / (double)(o);
               r = Math.Sqrt((double)(n) / 2.0);

               typeMethods.r8vec_zero(n * o, ref x);

               k = -1;
               //
               //  2 * N points.
               //
               for (i = 0; i < n; i++)
               {
                    k = k + 1;
                    x[i + k * n] = -r;
                    w[k] = a;
                    k = k + 1;
                    x[i + k * n] = +r;
                    w[k] = a;
               }

               return;
          }

          public static int en_r2_03_1_size(int n)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_R2_03_1_SIZE sizes the Stroud rule 3.1 for region EN_R2.
               //
               //  Discussion:
               //
               //    The rule has order O = 2 * N.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    23 January 2010
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
               //    Input, int N, the spatial dimension.
               //
               //    Output, int EN_R2_03_1_SIZE, the order.
               //
          {
               int o;

               o = 2 * n;

               return o;
          }

          public static void en_r2_03_2(int n, int o, ref double[] x, ref double[] w)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_R2_03_2 implements the Stroud rule 3.2 for region EN_R2.
               //
               //  Discussion:
               //
               //    The rule has order O = 2^N.
               //
               //    The rule has precision P = 3.
               //
               //    EN_R2 is the entire N-dimensional space with weight function
               //
               //      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    24 January 2010
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
               //    Input, int N, the spatial dimension.
               //
               //    Input, int O, the order.
               //
               //    Output, double X[N*O], the abscissas.
               //
               //    Output, double W[O], the weights.
               //
          {
               double a;
               int i;
               int i1;
               int k;
               bool more;
               double pi = 3.141592653589793E+00;
               double r;
               double volume;

               volume = Math.Sqrt(Math.Pow(pi, n));

               a = volume / (double)(o);
               r = Math.Sqrt(0.5);

               typeMethods.r8vec_zero(n * o, ref x);

               k = -1;
               //
               //  2^N points.
               //
               k = k + 1;
               for (i1 = 0; i1 < n; i1++)
               {
                    x[i1 + k * n] = -r;
               }

               w[k] = a;
               more = true;

               while (more)
               {
                    more = false;
                    for (i = n - 1; 0 <= i; i--)
                    {
                         if (x[i + k * n] < 0.0)
                         {
                              k = k + 1;
                              for (i1 = 0; i1 < i; i1++)
                              {
                                   x[i1 + k * n] = x[i1 + (k - 1) * n];
                              }

                              x[i + k * n] = +r;
                              for (i1 = i + 1; i1 < n; i1++)
                              {
                                   x[i1 + k * n] = -r;
                              }

                              w[k] = a;
                              more = true;
                              break;
                         }
                    }
               }

               return;
          }

          public static int en_r2_03_2_size(int n)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_R2_03_2_SIZE sizes the Stroud rule 3.2 for region EN_R2.
               //
               //  Discussion:
               //
               //    The rule has order O = 2^N.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    23 January 2010
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
               //    Input, int N, the spatial dimension.
               //
               //    Output, int EN_R2_03_2_SIZE, the order.
               //
          {
               int o;

               o = (int)Math.Pow(2, n);

               return o;
          }


          public static void en_r2_05_1(int n, int option, int o, ref double[] x, ref double[] w)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_R2_05_1 implements the Stroud rule 5.1 for region EN_R2.
               //
               //  Discussion:
               //
               //    The rule has order O = N^2 + N + 2.
               //
               //    The rule has precision P = 5.
               //
               //    EN_R2 is the entire N-dimensional space with weight function
               //
               //      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
               //
               //    For N = 3, 5 and 6, there are two versions of the rule, chosen by setting 
               //    the OPTION variable to 1 or 2.
               //
               //    Versions of this rule are only available for N = 2 through 7.
               //
               //    There is a typographical error in the reference.
               //    For the second version of the rule for N = 2, the line
               //      gamma =    0.313300683022281E+00
               //    should read
               //      gamma =    0.312200683022281E+00
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    24 January 2010
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
               //    Input, int N, the spatial dimension.
               //    2 <= N <= 7.
               //
               //    Input, int OPTION, selects option 1 or 2.
               //
               //    Input, int O, the order.
               //
               //    Output, double X[N*O], the abscissas.
               //
               //    Output, double W[O], the weights.
               //
          {
               double a = 0;
               double b = 0;
               double c = 0;
               double eta = 0;
               double gamma = 0;
               int i = 0;
               int i1 = 0;
               int j = 0;
               int k = 0;
               double lambda = 0;
               double mu = 0;
               double pi = 3.141592653589793E+00;
               double volume = 0;
               double xsi = 0;

               if (n < 2 || 7 < n)
               {
                    Console.WriteLine("");
                    Console.WriteLine("EN_R2_05_1 - Fatal error!");
                    Console.WriteLine("  2 <= N <= 7 required.");
                    return;
               }

               if (option < 1 || 2 < option)
               {
                    Console.WriteLine("");
                    Console.WriteLine("EN_R2_05_1 - Fatal error!");
                    Console.WriteLine("  1 <= OPTION <= 2 required.");
                    return;
               }

               if (option == 2)
               {
                    if (n != 3 && n != 5 && n != 6)
                    {
                         Console.WriteLine("");
                         Console.WriteLine("EN_R2_05_1 - Fatal error!");
                         Console.WriteLine("  OPTION = 2 requires N = 3, 5 or 6.");
                         return;
                    }
               }

               volume = Math.Sqrt(Math.Pow(pi, n));

               if (n == 2)
               {
                    eta = 0.446103183094540E+00;
                    lambda = 0.136602540378444E+01;
                    xsi = -0.366025403784439E+00;
                    mu = 0.198167882945871E+01;
                    gamma = 0.000000000000000E+00;
                    a = 0.328774019778636E+00 * volume;
                    b = 0.833333333333333E-01 * volume;
                    c = 0.455931355469736E-02 * volume;
               }
               else if (n == 3 && option == 1)
               {
                    eta = 0.476731294622796E+00;
                    lambda = 0.935429018879534E+00;
                    xsi = -0.731237647787132E+00;
                    mu = 0.433155309477649E+00;
                    gamma = 0.266922328697744E+01;
                    a = 0.242000000000000E+00 * volume;
                    b = 0.810000000000000E-01 * volume;
                    c = 0.500000000000000E-02 * volume;
               }
               //
               //  The value of gamma that follows corrects an error in the reference.
               //
               else if (n == 3 && option == 2)
               {
                    eta = 0.476731294622796E+00;
                    lambda = 0.128679320334269E+01;
                    xsi = -0.379873463323979E+00;
                    mu = -0.192386729447751E+01;
                    gamma = 0.312200683022281E+00;
                    a = 0.242000000000000E+00 * volume;
                    b = 0.810000000000000E-01 * volume;
                    c = 0.500000000000000E-02 * volume;
               }
               else if (n == 4)
               {
                    eta = 0.523945658287507E+00;
                    lambda = 0.119433782552719E+01;
                    xsi = -0.398112608509063E+00;
                    mu = -0.318569372920112E+00;
                    gamma = 0.185675837424096E+01;
                    a = 0.155502116982037E+00 * volume;
                    b = 0.777510584910183E-01 * volume;
                    c = 0.558227484231506E-02 * volume;
               }
               else if (n == 5 && option == 1)
               {
                    eta = 0.214972564378798E+01;
                    lambda = 0.464252986016289E+01;
                    xsi = -0.623201054093728E+00;
                    mu = -0.447108700673434E+00;
                    gamma = 0.812171426076311E+00;
                    a = 0.487749259189752E-03 * volume;
                    b = 0.487749259189752E-03 * volume;
                    c = 0.497073504444862E-01 * volume;
               }
               else if (n == 5 && option == 2)
               {
                    eta = 0.615369528365158E+00;
                    lambda = 0.132894698387445E+01;
                    xsi = -0.178394363877324E+00;
                    mu = -0.745963266507289E+00;
                    gamma = 0.135503972310817E+01;
                    a = 0.726415024414905E-01 * volume;
                    b = 0.726415024414905E-01 * volume;
                    c = 0.641509853510569E-02 * volume;
               }
               else if (n == 6 && option == 1)
               {
                    eta = 0.100000000000000E+01;
                    lambda = 0.141421356237309E+01;
                    xsi = 0.000000000000000E+00;
                    mu = -0.100000000000000E+01;
                    gamma = 0.100000000000000E+01;
                    a = 0.781250000000000E-02 * volume;
                    b = 0.625000000000000E-01 * volume;
                    c = 0.781250000000000E-02 * volume;
               }
               else if (n == 6 && option == 2)
               {
                    eta = 0.100000000000000E+01;
                    lambda = 0.942809041582063E+00;
                    xsi = -0.471404520791032E+00;
                    mu = -0.166666666666667E+01;
                    gamma = 0.333333333333333E+00;
                    a = 0.781250000000000E-02 * volume;
                    b = 0.625000000000000E-01 * volume;
                    c = 0.781250000000000E-02 * volume;
               }
               else if (n == 7)
               {
                    eta = 0.000000000000000E+00;
                    lambda = 0.959724318748357E+00;
                    xsi = -0.772326488820521E+00;
                    mu = -0.141214270131942E+01;
                    gamma = 0.319908106249452E+00;
                    a = 0.111111111111111E+00 * volume;
                    b = 0.138888888888889E-01 * volume;
                    c = 0.138888888888889E-01 * volume;
               }

               typeMethods.r8vec_zero(n * o, ref x);

               k = -1;
               //
               //  2 points.
               //
               k = k + 1;
               for (i1 = 0; i1 < n; i1++)
               {
                    x[i1 + k * n] = -eta;
               }

               w[k] = a;
               k = k + 1;
               for (i1 = 0; i1 < n; i1++)
               {
                    x[i1 + k * n] = +eta;
               }

               w[k] = a;
               //
               //  2 * N points.
               //
               for (i = 0; i < n; i++)
               {
                    k = k + 1;
                    for (i1 = 0; i1 < n; i1++)
                    {
                         x[i1 + k * n] = -xsi;
                    }

                    x[i + k * n] = -lambda;
                    w[k] = b;
                    k = k + 1;
                    for (i1 = 0; i1 < n; i1++)
                    {
                         x[i1 + k * n] = +xsi;
                    }

                    x[i + k * n] = +lambda;
                    w[k] = b;
               }

               //
               //  2 * ( N * ( N - 1 ) / 2 ) points.
               //
               for (i = 0; i < n - 1; i++)
               {
                    for (j = i + 1; j < n; j++)
                    {
                         k = k + 1;
                         for (i1 = 0; i1 < n; i1++)
                         {
                              x[i1 + k * n] = -gamma;
                         }

                         x[i + k * n] = -mu;
                         x[j + k * n] = -mu;
                         w[k] = c;
                         k = k + 1;
                         for (i1 = 0; i1 < n; i1++)
                         {
                              x[i1 + k * n] = +gamma;
                         }

                         x[i + k * n] = +mu;
                         x[j + k * n] = +mu;
                         w[k] = c;
                    }
               }

               return;
          }

          public static int en_r2_05_1_size(int n)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_R2_05_1_SIZE sizes the Stroud rule 5.1 for region EN_R2.
               //
               //  Discussion:
               //
               //    The rule has order O = N^2 + N + 2.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    23 January 2010
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
               //    Input, int N, the spatial dimension.
               //
               //    Output, int EN_R2_05_1_SIZE, the order.
               //
          {
               int o;

               o = n * n + n + 2;

               return o;
          }

          public static void en_r2_05_2(int n, int o, ref double[] x, ref double[] w)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_R2_05_2 implements the Stroud rule 5.2 for region EN_R2.
               //
               //  Discussion:
               //
               //    The rule has order O = 2 * N^2 + 1.
               //
               //    The rule has precision P = 5.
               //
               //    EN_R2 is the entire N-dimensional space with weight function
               //
               //      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    24 January 2010
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
               //    Input, int N, the spatial dimension.
               //
               //    Input, int O, the order.
               //
               //    Output, double X[N*O], the abscissas.
               //
               //    Output, double W[O], the weights.
               //
          {
               double a;
               double b;
               double c;
               int i;
               int j;
               int k;
               double pi = 3.141592653589793E+00;
               double r;
               double s;
               double volume;

               volume = Math.Sqrt(Math.Pow(pi, n));

               a = 2.0E+00 * volume / (double)(n + 2);
               b = (double)(4 - n) * volume / 2.0E+00
                                            / (double)((n + 2) * (n + 2));
               c = volume / (double)((n + 2) * (n + 2));

               r = Math.Sqrt((double)(n + 2) / 2.0E+00);
               s = Math.Sqrt((double)(n + 2) / 4.0E+00);

               typeMethods.r8vec_zero(n * o, ref x);

               k = -1;
               //
               //  1 point.
               //
               k = k + 1;
               w[k] = a;
               //
               //  2 * N points.
               //
               for (i = 0; i < n; i++)
               {
                    k = k + 1;
                    x[i + k * n] = -r;
                    w[k] = b;
                    k = k + 1;
                    x[i + k * n] = +r;
                    w[k] = b;
               }

               //
               //  4 * ( N * ( N - 1 ) / 2 ) points.
               //
               for (i = 0; i < n - 1; i++)
               {
                    for (j = i + 1; j < n; j++)
                    {
                         k = k + 1;
                         x[i + k * n] = -s;
                         x[j + k * n] = -s;
                         w[k] = c;
                         k = k + 1;
                         x[i + k * n] = -s;
                         x[j + k * n] = +s;
                         w[k] = c;
                         k = k + 1;
                         x[i + k * n] = +s;
                         x[j + k * n] = -s;
                         w[k] = c;
                         k = k + 1;
                         x[i + k * n] = +s;
                         x[j + k * n] = +s;
                         w[k] = c;
                    }
               }

               return;
          }

          public static int en_r2_05_2_size(int n)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_R2_05_2_SIZE sizes the Stroud rule 5.2 for region EN_R2.
               //
               //  Discussion:
               //
               //    The rule has order O = 2 * N^2 + 1.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    23 January 2010
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
               //    Input, int N, the spatial dimension.
               //
               //    Output, int EN_R2_01_1_SIZE, the order.
               //
          {
               int o;

               o = 2 * n * n + 1;

               return o;
          }

          public static void en_r2_05_3(int n, int o, ref double[] x, ref double[] w)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_R2_05_3 implements the Stroud rule 5.3 for region EN_R2.
               //
               //  Discussion:
               //
               //    The rule has order O = 2^N + 2 * N.
               //
               //    The rule has precision P = 5.
               //
               //    EN_R2 is the entire N-dimensional space with weight function
               //
               //      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
               //
               //    The rule requires 3 <= N.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    22 January 2010
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
               //    Input, int N, the spatial dimension.
               //    3 <= N.
               //
               //    Input, int O, the order.
               //
               //    Output, double X[N*O], the abscissas.
               //
               //    Output, double W[O], the weights.
               //
          {
               double a;
               double b;
               int i;
               int i1;
               int k;
               bool more;
               double pi = 3.141592653589793E+00;
               double r;
               double s;
               double volume;

               if (n < 3)
               {
                    Console.WriteLine("");
                    Console.WriteLine("EN_R2_05_3 - Fatal error!");
                    Console.WriteLine("  3 <= N is required.");
                    return;
               }

               volume = Math.Sqrt(Math.Pow(pi, n));

               a = 4.0E+00 * volume / (double)((n + 2) * (n + 2));
               b = (double)((n - 2) * (n - 2)) * volume / (double)((int)Math.Pow(2, n))
                                                        / (double)((n + 2) * (n + 2));
               r = Math.Sqrt((double)(n + 2) / 4.0E+00);
               s = Math.Sqrt((double)(n + 2) / 2.0E+00 / (double)(n - 2));

               typeMethods.r8vec_zero(n * o, ref x);

               k = -1;
               //
               //  2 * N points.
               //
               for (i = 0; i < n; i++)
               {
                    k = k + 1;
                    x[i + k * n] = -r;
                    w[k] = a;
                    k = k + 1;
                    x[i + k * n] = +r;
                    w[k] = a;
               }

               //
               //  2^N points.
               //
               k = k + 1;
               for (i1 = 0; i1 < n; i1++)
               {
                    x[i1 + k * n] = -s;
               }

               w[k] = b;
               more = true;
               while (more)
               {
                    more = false;
                    for (i = n - 1; 0 <= i; i--)
                    {
                         if (x[i + k * n] < 0.0E+00)
                         {
                              k = k + 1;
                              for (i1 = 0; i1 < n; i1++)
                              {
                                   x[i1 + k * n] = x[i1 + (k - 1) * n];
                              }

                              x[i + k * n] = +s;
                              for (i1 = i + 1; i1 < n; i1++)
                              {
                                   x[i1 + k * n] = -s;
                              }

                              w[k] = b;
                              more = true;
                              break;
                         }
                    }
               }

               return;
          }

          public static int en_r2_05_3_size(int n)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_R2_05_3_SIZE sizes the Stroud rule 5.3 for region EN_R2.
               //
               //  Discussion:
               //
               //    The rule has order O = 2^N + 2 * N.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    23 January 2010
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
               //    Input, int N, the spatial dimension.
               //
               //    Output, int EN_R2_05_3_SIZE, the order.
               //
          {
               int o;

               o = (int)Math.Pow(2, n) + 2 * n;

               return o;
          }

          public static void en_r2_05_4(int n, int o, ref double[] x, ref double[] w)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_R2_05_4 implements the Stroud rule 5.4 for region EN_R2.
               //
               //  Discussion:
               //
               //    The rule has order O = 2^(N+1) - 1.
               //
               //    The rule has precision P = 5.
               //
               //    EN_R2 is the entire N-dimensional space with weight function
               //
               //      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    22 January 2010
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
               //    Input, int N, the spatial dimension.
               //
               //    Input, int O, the order.
               //
               //    Output, double X[N*O], the abscissas.
               //
               //    Output, double W[O], the weights.
               //
          {
               double b;
               int i;
               int i1;
               int j;
               int k;
               bool more;
               double pi = 3.141592653589793E+00;
               double r;
               double s;
               double volume;

               volume = Math.Sqrt(Math.Pow(pi, n));

               s = Math.Sqrt(0.5E+00);

               typeMethods.r8vec_zero(n * o, ref x);

               k = -1;
               //
               //  2^N + 2^(N-1) + 2^(N-2) + ... + 1 = 2^(N+1)-1 points.
               //  but do the last point separately.
               //
               for (i = 0; i < n; i++)
               {
                    r = Math.Sqrt((double)(i + 3) / 2.0E+00);
                    b = Math.Pow(2.0E+00, i + 1 - n) * volume / (double)(i + 2)
                                                              / (double)(i + 3);

                    k = k + 1;
                    x[i + k * n] = -r;
                    for (i1 = i + 1; i1 < n; i1++)
                    {
                         x[i1 + k * n] = -s;
                    }

                    w[k] = b;
                    more = true;
                    while (more)
                    {
                         more = false;
                         for (j = n - 1; 0 <= j; j--)
                         {
                              if (x[j + k * n] < 0.0E+00)
                              {
                                   k = k + 1;
                                   for (i1 = 0; i1 < n; i1++)
                                   {
                                        x[i1 + k * n] = x[i1 + (k - 1) * n];
                                   }

                                   x[j + k * n] = Math.Abs(x[j + k * n]);
                                   for (i1 = j + 1; i1 < n; i1++)
                                   {
                                        x[i1 + k * n] = -Math.Abs(x[i1 + k * n]);
                                   }

                                   w[k] = b;
                                   more = true;
                                   break;
                              }
                         }
                    }
               }

               //
               //  Last point.
               //
               k = k + 1;
               w[k] = 2.0E+00 * volume / (double)(n + 2);

               return;
          }

          public static int en_r2_05_4_size(int n)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_R2_05_4_SIZE sizes the Stroud rule 5.4 for region EN_R2.
               //
               //  Discussion:
               //
               //    The rule has order O = 2^(N+1) - 1.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    23 January 2010
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
               //    Input, int N, the spatial dimension.
               //
               //    Output, int EN_R2_05_4_SIZE, the order.
               //
          {
               int o;

               o = (int)Math.Pow(2, n + 1) - 1;

               return o;
          }

          public static void en_r2_05_5(int n, int o, ref double[] x, ref double[] w)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_R2_05_5 implements the Stroud rule 5.5 for region EN_R2.
               //
               //  Discussion:
               //
               //    The rule has order O = N * 2^N + 1.
               //
               //    The rule has precision P = 5.
               //
               //    EN_R2 is the entire N-dimensional space with weight function
               //
               //      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
               //
               //    There is a second version of this rule however it results in
               //    complex abscissas, and so it has been disabled.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    21 January 2010
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
               //    Input, int N, the spatial dimension.
               //
               //    Input, int O, the order.
               //
               //    Output, double X[N*O], the abscissas.
               //
               //    Output, double W[O], the weights.
               //
          {
               double a = 0;
               double b = 0;
               int i = 0;
               int i1 = 0;
               int j = 0;
               int k = 0;
               bool more = false;
               double n_r8 = 0;
               int option = 0;
               double pi = 3.141592653589793E+00;
               double r = 0;
               double s = 0;
               double volume = 0;

               volume = Math.Sqrt(Math.Pow(pi, n));

               n_r8 = (double)(n);

               a = 2.0E+00 * volume / (n_r8 + 2.0E+00);
               b = volume / (n_r8 + 2.0E+00) / Math.Pow(2.0, n);

               option = 1;

               if (option == 1)
               {
                    r = Math.Sqrt((n_r8 + 2.0E+00
                                        + (n_r8 - 1.0E+00) * Math.Sqrt(2.0E+00 * (n_r8 + 2.0E+00)))
                                  / 2.0E+00 / n_r8);
                    s = Math.Sqrt((n_r8 + 2.0E+00
                                   - Math.Sqrt(2.0E+00 * (n_r8 + 2.0E+00)))
                                  / 2.0E+00 / n_r8);
               }
               else if (option == 2)
               {
                    r = Math.Sqrt((n_r8 + 2.0E+00
                                   - (n_r8 - 1.0E+00) * Math.Sqrt(2.0E+00 * (n_r8 + 2.0E+00)))
                                  / 2.0E+00 / n_r8);
                    s = Math.Sqrt((n_r8 + 2.0E+00
                                        + Math.Sqrt(2.0E+00 * (n_r8 + 2.0E+00)))
                                  / 2.0E+00 / n_r8);
               }

               typeMethods.r8vec_zero(n * o, ref x);

               k = -1;
               //
               //  1 point.
               //
               k = k + 1;
               w[k] = a;
               //
               //  N * 2^N points:
               //  N choices for location of R, 2^N choices of sign pattern.
               //
               for (i = 0; i < n; i++)
               {
                    k = k + 1;
                    for (i1 = 0; i1 < n; i1++)
                    {
                         x[i1 + k * n] = -s;
                    }

                    x[i + k * n] = -r;
                    w[k] = b;

                    more = true;

                    while (more)
                    {
                         more = false;
                         for (j = n - 1; 0 <= j; j--)
                         {
                              if (x[j + k * n] < 0.0E+00)
                              {
                                   k = k + 1;
                                   for (i1 = 0; i1 < n; i1++)
                                   {
                                        x[i1 + k * n] = x[i1 + (k - 1) * n];
                                   }

                                   x[j + k * n] = Math.Abs(x[j + k * n]);
                                   for (i1 = j + 1; i1 < n; i1++)
                                   {
                                        x[i1 + k * n] = -Math.Abs(x[i1 + k * n]);
                                   }

                                   w[k] = b;
                                   more = true;
                                   break;
                              }
                         }
                    }
               }

               return;
          }

          public static int en_r2_05_5_size(int n)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_R2_05_5_SIZE sizes the Stroud rule 5.5 for region EN_R2.
               //
               //  Discussion:
               //
               //    The rule has order O = N * 2^N + 1.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    23 January 2010
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
               //    Input, int N, the spatial dimension.
               //
               //    Output, int EN_R2_05_5_SIZE, the order.
               //
          {
               int o;

               o = n * (int)Math.Pow(2, n) + 1;

               return o;
          }

          public static void en_r2_05_6(int n, int o, ref double[] x, ref double[] w)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_R2_05_6 implements the Stroud rule 5.6 for region EN_R2.
               //
               //  Discussion:
               //
               //    The rule has order O = ( N + 1 ) * 2^N.
               //
               //    The rule has precision P = 5.
               //
               //    EN_R2 is the entire N-dimensional space with weight function
               //
               //      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
               //
               //    The rule requires 5 <= N.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    24 January 2010
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
               //    Input, int N, the spatial dimension.
               //    5 <= N.
               //
               //    Input, int O, the order.
               //
               //    Output, double X[N*O], the abscissas.
               //
               //    Output, double W[O], the weights.
               //
          {
               double a;
               int i;
               int i1;
               int j;
               int k;
               bool more;
               double n_r8;
               double pi = 3.141592653589793E+00;
               double r;
               double s;
               double t;
               double volume;

               if (n < 5)
               {
                    Console.WriteLine("");
                    Console.WriteLine("EN_R2_05_6 - Fatal error!");
                    Console.WriteLine("  5 <= N is required.");
                    return;
               }

               volume = Math.Sqrt(Math.Pow(pi, n));

               n_r8 = (double)(n);

               a = volume / Math.Pow(2.0, n) / (n_r8 + 1.0E+00);

               r = Math.Sqrt((n_r8 - Math.Sqrt(2.0E+00)
                              + (n_r8 - 1.0E+00) * Math.Sqrt(2.0E+00 * (n_r8 + 1.0E+00)))
                             / 2.0E+00 / n_r8);
               s = Math.Sqrt((n_r8 - Math.Sqrt(2.0E+00)
                                   - Math.Sqrt(2.0E+00 * (n_r8 + 1.0E+00)))
                             / 2.0E+00 / n_r8);
               t = Math.Sqrt((1.0E+00 + Math.Sqrt(2.0E+00)) / 2.0E+00);

               typeMethods.r8vec_zero(n * o, ref x);

               k = -1;
               //
               //  N * 2^N points.
               //
               for (i = 0; i < n; i++)
               {
                    k = k + 1;
                    for (i1 = 0; i1 < n; i1++)
                    {
                         x[i1 + k * n] = -s;
                    }

                    x[i + k * n] = -r;
                    w[k] = a;

                    more = true;

                    while (more)
                    {
                         more = false;
                         for (j = n - 1; 0 <= j; j--)
                         {
                              if (x[j + k * n] < 0.0E+00)
                              {
                                   k = k + 1;
                                   for (i1 = 0; i1 < n; i1++)
                                   {
                                        x[i1 + k * n] = x[i1 + (k - 1) * n];
                                   }

                                   x[j + k * n] = Math.Abs(x[j + k * n]);
                                   for (i1 = j + 1; i1 < n; i1++)
                                   {
                                        x[i1 + k * n] = -Math.Abs(x[i1 + k * n]);
                                   }

                                   w[k] = a;
                                   more = true;
                                   break;
                              }
                         }
                    }
               }

               //
               //  2^N points.
               //
               k = k + 1;
               for (i1 = 0; i1 < n; i1++)
               {
                    x[i1 + k * n] = -t;
               }

               w[k] = a;
               more = true;
               while (more)
               {
                    more = false;
                    for (j = n - 1; 0 <= j; j--)
                    {
                         if (x[j + k * n] < 0.0E+00)
                         {
                              k = k + 1;
                              for (i1 = 0; i1 < n; i1++)
                              {
                                   x[i1 + k * n] = x[i1 + (k - 1) * n];
                              }

                              x[j + k * n] = Math.Abs(x[j + k * n]);
                              for (i1 = j + 1; i1 < n; i1++)
                              {
                                   x[i1 + k * n] = -Math.Abs(x[i1 + k * n]);
                              }

                              w[k] = a;
                              more = true;
                              break;
                         }
                    }
               }

               return;
          }

          public static int en_r2_05_6_size(int n)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_R2_05_6_SIZE sizes the Stroud rule 5.6 for region EN_R2.
               //
               //  Discussion:
               //
               //    The rule has order O = ( N + 1 ) * 2^N.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    23 January 2010
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
               //    Input, int N, the spatial dimension.
               //
               //    Output, int EN_R2_05_6_SIZE, the order.
               //
          {
               int o;

               o = (n + 1) * (int)Math.Pow(2, n);

               return o;
          }

          public static void en_r2_07_1(int n, int option, int o, ref double[] x, ref double[] w)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_R2_07_1 implements the Stroud rule 7.1 for region EN_R2.
               //
               //  Discussion:
               //
               //    The rule has order O = 2^N + 2 * N^2 + 1.
               //
               //    The rule has precision P = 7.
               //
               //    EN_R2 is the entire N-dimensional space with weight function
               //
               //      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
               //
               //    There are two versions of the rule, chosen by setting the
               //    OPTION variable to 1 or 2.  
               //
               //    Option 1 is only valid for N = 3, 4, 6 or 7.
               //    Option 2 is only valid for N = 3 or 4.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    23 January 2010
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
               //    Input, int N, the spatial dimension.
               //    N = 3, 4, 6 or 7.
               //
               //    Input, int OPTION, chooses rule option 1 or 2.
               //
               //    Input, int O, the order.
               //
               //    Output, double X[N*O], the abscissas.
               //
               //    Output, double W[O], the weights.
               //
          {
               double a = 0;
               double b = 0;
               double c = 0;
               double d = 0;
               int i = 0;
               int i1 = 0;
               int j = 0;
               int k = 0;
               bool more = false;
               double n_r8 = 0;
               double pi = 3.141592653589793E+00;
               double r = 0;
               double s = 0;
               double t = 0;
               double volume = 0;

               if (option < 1 || 2 < option)
               {
                    Console.WriteLine("");
                    Console.WriteLine("EN_R2_07_1 - Fatal error!");
                    Console.WriteLine("  1 <= OPTION <= 2 required.");
                    return;
               }

               if (option == 1)
               {
                    if (n != 3 && n != 4 && n != 6 && n != 7)
                    {
                         Console.WriteLine("");
                         Console.WriteLine("EN_R2_07_1 - Fatal error!");
                         Console.WriteLine("  OPTION 1 requires N =  3, 4, 6 or 7.");
                         return;
                    }
               }

               if (option == 2)
               {
                    if (n != 3 && n != 4)
                    {
                         Console.WriteLine("");
                         Console.WriteLine("EN_R2_07_1 - Fatal error!");
                         Console.WriteLine("  OPTION 2 requires N =  3 or 4.");
                         return;
                    }
               }

               volume = Math.Sqrt(Math.Pow(pi, n));

               n_r8 = (double)(n);

               if (option == 1)
               {
                    r = Math.Sqrt((3.0E+00 * (8.0E+00 - n_r8) - (n_r8 - 2.0E+00)
                         * Math.Sqrt(3.0E+00 * (8.0E+00 - n_r8))) / 2.0E+00 / (5.0E+00 - n_r8));
                    s = Math.Sqrt((3.0E+00 * n_r8 - 2.0E+00
                              * Math.Sqrt(3.0E+00 * (8.0E+00 - n_r8))) / 2.0E+00
                                                                       / (3.0E+00 * n_r8 - 8.0E+00));
                    t = Math.Sqrt((6.0E+00 + Math.Sqrt(3.0E+00 * (8.0E+00 - n_r8))) / 2.0E+00);
               }
               else if (option == 2)
               {
                    r = Math.Sqrt((3.0E+00 * (8.0E+00 - n_r8) + (n_r8 - 2.0E+00)
                         * Math.Sqrt(3.0E+00 * (8.0E+00 - n_r8))) / 2.0E+00 / (5.0E+00 - n_r8));
                    s = Math.Sqrt((3.0E+00 * n_r8 + 2.0E+00
                              * Math.Sqrt(3.0E+00 * (8.0E+00 - n_r8))) / 2.0E+00
                                                                       / (3.0E+00 * n_r8 - 8.0E+00));
                    t = Math.Sqrt((6.0E+00 - Math.Sqrt(3.0E+00 * (8.0E+00 - n_r8))) / 2.0E+00);
               }

               b = (8.0E+00 - n_r8) * volume / 8.0E+00 / Math.Pow(r, 6);
               c = volume / Math.Pow(2.0E+00, n + 3) / Math.Pow(s, 6);
               d = volume / 16.0E+00 / Math.Pow(t, 6);
               a = volume - 2.0E+00 * n_r8 * b - Math.Pow(2.0E+00, n) * c - 2.0E+00 * n_r8
                    * (n_r8 - 1.0E+00) * d;

               typeMethods.r8vec_zero(n * o, ref x);

               k = -1;
               //
               //  1 point.
               //
               k = k + 1;
               w[k] = a;
               //
               //  2 * N points.
               //
               for (i = 0; i < n; i++)
               {
                    k = k + 1;
                    x[i + k * n] = -r;
                    w[k] = b;
                    k = k + 1;
                    x[i + k * n] = +r;
                    w[k] = b;
               }

               //
               //  2^N points.
               //
               k = k + 1;
               for (i1 = 0; i1 < n; i1++)
               {
                    x[i1 + k * n] = -s;
               }

               w[k] = c;
               more = true;
               while (more)
               {
                    more = false;
                    for (i = n - 1; 0 <= i; i--)
                    {
                         if (x[i + k * n] < 0.0E+00)
                         {
                              k = k + 1;
                              for (i1 = 0; i1 < n; i1++)
                              {
                                   x[i1 + k * n] = x[i1 + (k - 1) * n];
                              }

                              x[i + k * n] = Math.Abs(x[i + k * n]);
                              for (i1 = i + 1; i1 < n; i1++)
                              {
                                   x[i1 + k * n] = -Math.Abs(x[i1 + k * n]);
                              }

                              w[k] = c;
                              more = true;
                              break;
                         }
                    }
               }

               //
               //  2 * ( N * ( N - 1 ) / 2 ) points.
               //
               for (i = 0; i < n - 1; i++)
               {
                    for (j = i + 1; j < n; j++)
                    {
                         k = k + 1;
                         x[i + k * n] = -t;
                         x[j + k * n] = -t;
                         w[k] = d;
                         k = k + 1;
                         x[i + k * n] = -t;
                         x[j + k * n] = +t;
                         w[k] = d;
                         k = k + 1;
                         x[i + k * n] = +t;
                         x[j + k * n] = -t;
                         w[k] = d;
                         k = k + 1;
                         x[i + k * n] = +t;
                         x[j + k * n] = +t;
                         w[k] = d;
                    }
               }

               return;
          }

          public static int en_r2_07_1_size(int n)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_R2_07_1_SIZE sizes the Stroud rule 7.1 for region EN_R2.
               //
               //  Discussion:
               //
               //    The rule has order O = 2^N + 2 * N^2 + 1.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    23 January 2010
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
               //    Input, int N, the spatial dimension.
               //
               //    Output, int EN_R2_07_1_SIZE, the order.
               //
          {
               int o;

               o = (int)Math.Pow(2, n) + 2 * n * n + 1;

               return o;
          }

          public static void en_r2_07_2(int n, int o, ref double[] x, ref double[] w)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_R2_07_2 implements the Stroud rule 7.2 for region EN_R2.
               //
               //  Discussion:
               //
               //    The rule has order O = 2^(N+1) + 4 * N^2.
               //
               //    The rule has precision P = 7.
               //
               //    EN_R2 is the entire N-dimensional space with weight function
               //
               //      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
               //
               //    The rule requires 3 <= N.
               //
               //    The reference has a typographical error in the description of this rule.
               //    The formula:
               //
               //      (t,t,t,...,t)FS
               //
               //    should read
               //
               //      (t,t,0,...,0)FS.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    23 January 2010
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
               //    Input, int N, the spatial dimension.
               //    3 <= N.
               //
               //    Input, int O, the order.
               //
               //    Output, double X[N*O], the abscissas.
               //
               //    Output, double W[O], the weights.
               //
          {
               double a1;
               double a2;
               double b;
               double c;
               double d;
               int i;
               int i1;
               int j;
               int k;
               bool more;
               double n_r8;
               double pi = 3.141592653589793E+00;
               double r;
               double rho1;
               double rho2;
               double s;
               double t;
               double volume;

               if (n < 3)
               {
                    Console.WriteLine("");
                    Console.WriteLine("EN_R2_07_2 - Fatal error!");
                    Console.WriteLine("  3 <= N is required.");
                    return;
               }

               volume = Math.Sqrt(Math.Pow(pi, n));

               n_r8 = (double)(n);

               rho1 = Math.Sqrt((n_r8 + 2.0E+00 - Math.Sqrt(2.0E+00 * (n_r8 + 2.0E+00)))
                                / 2.0E+00);
               rho2 = Math.Sqrt((n_r8 + 2.0E+00 + Math.Sqrt(2.0E+00 * (n_r8 + 2.0E+00)))
                                / 2.0E+00);
               a1 = (n_r8 + 2.0E+00 + Math.Sqrt(2.0E+00 * (n_r8 + 2.0E+00))) / 2.0E+00
                                                                             / (n_r8 + 2.0E+00);
               a2 = (n_r8 + 2.0E+00 - Math.Sqrt(2.0E+00 * (n_r8 + 2.0E+00))) / 2.0E+00
                                                                             / (n_r8 + 2.0E+00);

               r = 1.0E+00;
               s = Math.Sqrt(1.0E+00 / n_r8);
               t = Math.Sqrt(0.5E+00);
               b = (8.0E+00 - n_r8) * volume / n_r8 / (n_r8 + 2.0E+00) / (n_r8 + 4.0E+00);
               c = Math.Pow(n_r8, 3) * volume / Math.Pow(2.0E+00, n) / n_r8 / (n_r8 + 2.0E+00)
                   / (n_r8 + 4.0E+00);
               d = 4.0E+00 * volume / n_r8 / (n_r8 + 2.0E+00) / (n_r8 + 4.0E+00);

               typeMethods.r8vec_zero(n * o, ref x);

               k = -1;
               //
               //  2 * 2 * N points.
               //
               for (i = 0; i < n; i++)
               {
                    k = k + 1;
                    x[i + k * n] = -rho1 * r;
                    w[k] = a1 * b;
                    k = k + 1;
                    x[i + k * n] = -rho2 * r;
                    w[k] = a2 * b;
                    k = k + 1;
                    x[i + k * n] = +rho1 * r;
                    w[k] = a1 * b;
                    k = k + 1;
                    x[i + k * n] = +rho2 * r;
                    w[k] = a2 * b;
               }

               //
               //  2 * 2^N points.
               //
               k = k + 1;
               for (i1 = 0; i1 < n; i1++)
               {
                    x[i1 + k * n] = -rho1 * s;
               }

               w[k] = a1 * c;
               k = k + 1;
               for (i1 = 0; i1 < n; i1++)
               {
                    x[i1 + k * n] = -rho2 * s;
               }

               w[k] = a2 * c;
               more = true;
               while (more)
               {
                    more = false;
                    for (i = n - 1; 0 <= i; i--)
                    {
                         if (x[i + k * n] < 0.0E+00)
                         {
                              k = k + 1;
                              for (i1 = 0; i1 < n; i1++)
                              {
                                   x[i1 + k * n] = x[i1 + (k - 2) * n];
                                   ;
                              }

                              x[i + k * n] = Math.Abs(x[i + k * n]);
                              for (i1 = i + 1; i1 < n; i1++)
                              {
                                   x[i1 + k * n] = -Math.Abs(x[i1 + k * n]);
                              }

                              w[k] = a1 * c;
                              k = k + 1;
                              for (i1 = 0; i1 < n; i1++)
                              {
                                   x[i1 + k * n] = x[i1 + (k - 2) * n];
                                   ;
                              }

                              x[i + k * n] = Math.Abs(x[i + k * n]);
                              for (i1 = i + 1; i1 < n; i1++)
                              {
                                   x[i1 + k * n] = -Math.Abs(x[i1 + k * n]);
                              }

                              w[k] = a2 * c;
                              more = true;
                              break;
                         }
                    }
               }

               //
               //  2 * 4 * ( N * ( N - 1 ) / 2 ) points.
               //
               for (i = 0; i < n - 1; i++)
               {
                    for (j = i + 1; j < n; j++)
                    {
                         k = k + 1;
                         x[i + k * n] = -rho1 * t;
                         x[j + k * n] = -rho1 * t;
                         w[k] = a1 * d;
                         k = k + 1;
                         x[i + k * n] = -rho1 * t;
                         x[j + k * n] = +rho1 * t;
                         w[k] = a1 * d;
                         k = k + 1;
                         x[i + k * n] = +rho1 * t;
                         x[j + k * n] = -rho1 * t;
                         w[k] = a1 * d;
                         k = k + 1;
                         x[i + k * n] = +rho1 * t;
                         x[j + k * n] = +rho1 * t;
                         w[k] = a1 * d;
                         k = k + 1;
                         x[i + k * n] = -rho2 * t;
                         x[j + k * n] = -rho2 * t;
                         w[k] = a2 * d;
                         k = k + 1;
                         x[i + k * n] = -rho2 * t;
                         x[j + k * n] = +rho2 * t;
                         w[k] = a2 * d;
                         k = k + 1;
                         x[i + k * n] = +rho2 * t;
                         x[j + k * n] = -rho2 * t;
                         w[k] = a2 * d;
                         k = k + 1;
                         x[i + k * n] = +rho2 * t;
                         x[j + k * n] = +rho2 * t;
                         w[k] = a2 * d;
                    }
               }

               return;
          }

          public static int en_r2_07_2_size(int n)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_R2_07_2_SIZE sizes the Stroud rule 7.2 for region EN_R2.
               //
               //  Discussion:
               //
               //    The rule has order O = 2^(N+1) + 4 * N^2.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    23 January 2010
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
               //    Input, int N, the spatial dimension.
               //
               //    Output, int EN_R2_07_2_SIZE, the order.
               //
          {
               int o;

               o = (int)Math.Pow(2, n + 1) + 4 * n * n;

               return o;
          }

          public static void en_r2_07_3(int n, int option, int o, ref double[] x, ref double[] w)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_R2_07_3 implements the Stroud rule 7.3 for region EN_R2.
               //
               //  Discussion:
               //
               //    The rule has order O = ( 4 * N^3 + 8 * N + 3 ) / 3.
               //
               //    The rule has precision P = 7.
               //
               //    EN_R2 is the entire N-dimensional space with weight function
               //
               //      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
               //
               //    There are two versions of each rule, chosen by setting the
               //    OPTION variable to 1 or 2.
               //
               //    The rule as tabulated by Stenger is available for N = 2 through 20.
               //    This function accepts N = 3 through 6.
               //
               //     N    O
               //    __  ___
               //     3   45
               //     4   97
               //     5  181
               //     6  305
               //
               //    The reference has a typographical error for N = 5, OPTION 1, B4:
               //
               //      -(1)0.736330882774831
               //
               //    should read
               //
               //      (-1)0.736330882774831
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    23 January 2010
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
               //    Input, int N, the spatial dimension.
               //    3 <= N <= 6.
               //
               //    Input, int OPTION, chooses rule option 1 or 2.
               //
               //    Input, int O, the order.
               //
               //    Output, double X[N*O], the abscissas.
               //
               //    Output, double W[O], the weights.
               //
          {
               double b0 = 0;
               double b1 = 0;
               double b2 = 0;
               double b3 = 0;
               double b4 = 0;
               double b5 = 0;
               int i = 0;
               int j = 0;
               int k = 0;
               int l = 0;
               double u = 0;
               double v = 0;

               if (n < 3 || 6 < n)
               {
                    Console.WriteLine("");
                    Console.WriteLine("EN_R2_07_3 - Fatal error!");
                    Console.WriteLine("  3 <= N <= 6 required.");
                    return;
               }

               if (option < 1 || 2 < option)
               {
                    Console.WriteLine("");
                    Console.WriteLine("EN_R2_07_3 - Fatal error!");
                    Console.WriteLine("  1 <= OPTION <= 2 required.");
                    return;
               }

               if (n == 3 && option == 1)
               {
                    u = 0.524647623275290E+00;
                    v = 0.165068012388578E+01;
                    b0 = -0.166705761599566E+02;
                    b1 = 0.100296981655678E+02;
                    b2 = 0.161699246687754E+00;
                    b3 = -0.604719151221535E+01;
                    b4 = 0.234381399489666E-01;
                    b5 = 0.417194501880647E+01;
               }
               else if (n == 3 && option == 2)
               {
                    u = 0.165068012388578E+01;
                    v = 0.524647623275290E+00;
                    b0 = 0.166705761599566E+02;
                    b1 = 0.178903161957074E+00;
                    b2 = -0.665808190965810E+01;
                    b3 = 0.148361823143070E-01;
                    b4 = 0.229669852539758E+01;
                    b5 = 0.430097881732984E-02;
               }
               else if (n == 4 && option == 1)
               {
                    u = 0.524647623275290E+00;
                    v = 0.165068012388578E+01;
                    b0 = -0.167539329651562E+03;
                    b1 = 0.687922329603575E+02;
                    b2 = 0.203518409659014E+00;
                    b3 = -0.255075279116885E+02;
                    b4 = 0.415430214106084E-01;
                    b5 = 0.739458001434961E+01;
               }
               else if (n == 4 && option == 2)
               {
                    u = 0.165068012388578E+01;
                    v = 0.524647623275290E+00;
                    b0 = 0.688432856406677E+02;
                    b1 = 0.294997847268286E+00;
                    b2 = -0.199427272118378E+02;
                    b3 = 0.110498755408511E-01;
                    b4 = 0.407079214570997E+01;
                    b5 = 0.762328646743931E-02;
               }
               else if (n == 5 && option == 1)
               {
                    u = 0.524647623275290E+00;
                    v = 0.165068012388578E+01;
                    b0 = -0.826940846964452E+03;
                    b1 = 0.264779097660331E+03;
                    b2 = 0.213460812375320E+00;
                    b3 = -0.714240197186780E+02;
                    b4 = 0.736330882774831E-01;
                    b5 = 0.131065518222629E+02;
               }
               else if (n == 5 && option == 2)
               {
                    u = 0.165068012388578E+01;
                    v = 0.524647623275290E+00;
                    b0 = 0.220502344940121E+03;
                    b1 = 0.537746975313769E+00;
                    b2 = -0.497781460739792E+02;
                    b3 = -0.743845245712926E-02;
                    b4 = 0.721529121489956E+01;
                    b5 = 0.135119234557687E-01;
               }
               else if (n == 6 && option == 1)
               {
                    u = 0.524647623275290E+00;
                    v = 0.165068012388578E+01;
                    b0 = -0.309679578630802E+04;
                    b1 = 0.815423321880237E+03;
                    b2 = 0.117326937169073E+00;
                    b3 = -0.173057295296448E+03;
                    b4 = 0.130511250871491E+00;
                    b5 = 0.232307582494626E+02;
               }
               else if (n == 6 && option == 2)
               {
                    u = 0.165068012388578E+01;
                    v = 0.524647623275290E+00;
                    b0 = 0.616293651884027E+03;
                    b1 = 0.107529736766179E+01;
                    b2 = -0.113807008098269E+03;
                    b3 = -0.610828352270520E-01;
                    b4 = 0.127887706992535E+02;
                    b5 = 0.239492607623178E-01;
               }

               typeMethods.r8vec_zero(n * o, ref x);

               k = -1;
               //
               //  1 point.
               //
               k = k + 1;
               w[k] = b0;
               //
               //  2 * N points.
               //
               for (i = 0; i < n; i++)
               {
                    k = k + 1;
                    x[i + k * n] = -u;
                    w[k] = b1;
                    k = k + 1;
                    x[i + k * n] = +u;
                    w[k] = b1;
               }

               //
               //  2 * N points.
               //
               for (i = 0; i < n; i++)
               {
                    k = k + 1;
                    x[i + k * n] = -v;
                    w[k] = b2;
                    k = k + 1;
                    x[i + k * n] = +v;
                    w[k] = b2;
               }

               //
               //  4 * ( N * ( N - 1 ) / 2 ) points.
               //
               for (i = 0; i < n - 1; i++)
               {
                    for (j = i + 1; j < n; j++)
                    {
                         k = k + 1;
                         x[i + k * n] = -u;
                         x[j + k * n] = -u;
                         w[k] = b3;
                         k = k + 1;
                         x[i + k * n] = -u;
                         x[j + k * n] = +u;
                         w[k] = b3;
                         k = k + 1;
                         x[i + k * n] = +u;
                         x[j + k * n] = -u;
                         w[k] = b3;
                         k = k + 1;
                         x[i + k * n] = +u;
                         x[j + k * n] = +u;
                         w[k] = b3;
                    }
               }

               //
               //  4 * ( N * ( N - 1 ) / 2 ) points.
               //
               for (i = 0; i < n - 1; i++)
               {
                    for (j = i + 1; j < n; j++)
                    {
                         k = k + 1;
                         x[i + k * n] = -v;
                         x[j + k * n] = -v;
                         w[k] = b4;
                         k = k + 1;
                         x[i + k * n] = -v;
                         x[j + k * n] = +v;
                         w[k] = b4;
                         k = k + 1;
                         x[i + k * n] = +v;
                         x[j + k * n] = -v;
                         w[k] = b4;
                         k = k + 1;
                         x[i + k * n] = +v;
                         x[j + k * n] = +v;
                         w[k] = b4;
                    }
               }

               //
               //  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
               //
               for (i = 0; i < n - 2; i++)
               {
                    for (j = i + 1; j < n - 1; j++)
                    {
                         for (l = j + 1; l < n; l++)
                         {
                              k = k + 1;
                              x[i + k * n] = -u;
                              x[j + k * n] = -u;
                              x[l + k * n] = -u;
                              w[k] = b5;
                              k = k + 1;
                              x[i + k * n] = -u;
                              x[j + k * n] = -u;
                              x[l + k * n] = +u;
                              w[k] = b5;
                              k = k + 1;
                              x[i + k * n] = -u;
                              x[j + k * n] = +u;
                              x[l + k * n] = -u;
                              w[k] = b5;
                              k = k + 1;
                              x[i + k * n] = -u;
                              x[j + k * n] = +u;
                              x[l + k * n] = +u;
                              w[k] = b5;
                              k = k + 1;
                              x[i + k * n] = +u;
                              x[j + k * n] = -u;
                              x[l + k * n] = -u;
                              w[k] = b5;
                              k = k + 1;
                              x[i + k * n] = +u;
                              x[j + k * n] = -u;
                              x[l + k * n] = +u;
                              w[k] = b5;
                              k = k + 1;
                              x[i + k * n] = +u;
                              x[j + k * n] = +u;
                              x[l + k * n] = -u;
                              w[k] = b5;
                              k = k + 1;
                              x[i + k * n] = +u;
                              x[j + k * n] = +u;
                              x[l + k * n] = +u;
                              w[k] = b5;
                         }
                    }
               }

               return;
          }

          public static int en_r2_07_3_size(int n)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_R2_07_3_SIZE sizes the Stroud rule 7.3 for region EN_R2.
               //
               //  Discussion:
               //
               //    The rule has order O = ( 4 * N^3 + 8 * N + 3 ) / 3.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    23 January 2010
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
               //    Input, int N, the spatial dimension.
               //
               //    Output, int EN_R2_07_3_SIZE, the order.
               //
          {
               int o;

               o = (4 * n * n * n + 8 * n + 3) / 3;

               return o;
          }

          public static void en_r2_09_1(int n, int option, int o, ref double[] x, ref double[] w)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_R2_09_1 implements the Stroud rule 9.1 for region EN_R2.
               //
               //  Discussion:
               //
               //    The rule has order O = ( 2 * N^4 - 4 * N^3 + 22 * N^2 - 8 * N + 3 ) / 3.
               //
               //    The rule has precision P = 9.
               //
               //    EN_R2 is the entire N-dimensional space with weight function
               //
               //      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
               //
               //    There are two versions of each rule, chosen by setting the 
               //    OPTION variable to 1 or 2.
               //
               //    The rule as tabulated by Stenger is available for N = 2 through 20.
               //    This function accepts N = 3 through 6.
               //
               //     N    O
               //    __  ___
               //     3   77
               //     4  193
               //     5  421
               //     6  825
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    24 January 2010
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
               //    Input, int N, the spatial dimension.
               //    3 <= N <= 6.
               //
               //    Input, int OPTION, chooses rule option 1 or 2.
               //
               //    Input, int O, the order.
               //
               //    Output, double X[N*O], the abscissas.
               //
               //    Output, double W[O], the weights.
               //
          {
               double b0 = 0;
               double b1 = 0;
               double b2 = 0;
               double b3 = 0;
               double b4 = 0;
               double b5 = 0;
               double b6 = 0;
               double b7 = 0;
               double b8 = 0;
               int i = 0;
               int j = 0;
               int k = 0;
               int l = 0;
               int m = 0;
               double u = 0;
               double v = 0;

               if (n < 3 || 6 < n)
               {
                    Console.WriteLine("");
                    Console.WriteLine("EN_R2_09_1 - Fatal error!");
                    Console.WriteLine("  3 <= N <= 6 required.");
                    return;
               }

               if (option < 1 || 2 < option)
               {
                    Console.WriteLine("");
                    Console.WriteLine("EN_R2_09_1 - Fatal error!");
                    Console.WriteLine("  1 <= OPTION <= 2 required.");
                    return;
               }

               if (n == 3)
               {
                    u = 0.202018287045609E+01;
                    v = 0.958572464613819E+00;
                    b0 = 0.676448734429924E+00;
                    b1 = 0.511989106291551E-02;
                    b2 = 0.448595723493744E+00;
                    b3 = 0.235223454595606E-03;
                    b4 = 0.915390713080005E-01;
                    b5 = 0.139208199920793E-01;
                    b6 = 0.235223454595606E-03;
                    b7 = 0.915390713080008E-01;
                    b8 = 0.000000000000000E+00;
               }
               else if (n == 4 && option == 1)
               {
                    u = 0.202018287045609E+01;
                    v = 0.958572464613819E+00;
                    b0 = -0.860452945007048E+00;
                    b1 = -0.405511998533795E-01;
                    b2 = 0.107026475449715E+01;
                    b3 = 0.138974239307092E-03;
                    b4 = -0.162248779448181E+00;
                    b5 = 0.246740110027234E-01;
                    b6 = 0.138974239307094E-03;
                    b7 = 0.162248779448181E+00;
                    b8 = 0.138974239307094E-03;
               }
               else if (n == 4 && option == 2)
               {
                    u = 0.958572464613819E+00;
                    v = 0.202018287045609E+01;
                    b0 = 0.265029088766810E-02;
                    b1 = 0.637601342635332E+00;
                    b2 = -0.394394059389228E-01;
                    b3 = 0.540829264827264E-01;
                    b4 = -0.416922717921281E-03;
                    b5 = 0.246740110027234E-01;
                    b6 = 0.540829264827270E-01;
                    b7 = 0.416922717921281E-03;
                    b8 = 0.540829264827269E-01;
               }
               else if (n == 5 && option == 1)
               {
                    u = 0.202018287045609E+01;
                    v = 0.958572464613819E+00;
                    b0 = -0.827347006200826E+01;
                    b1 = -0.160820174530905E+00;
                    b2 = 0.353499863758467E+01;
                    b3 = 0.738976276909564E-03;
                    b4 = -0.862735421812943E+00;
                    b5 = 0.437335458190621E-01;
                    b6 = -0.246325425636523E-03;
                    b7 = 0.287578473937648E+00;
                    b8 = 0.246325425636523E-03;
               }
               else if (n == 5 && option == 2)
               {
                    u = 0.958572464613819E+00;
                    v = 0.202018287045609E+01;
                    b0 = -0.624416791055272E+00;
                    b1 = 0.467494915583104E+00;
                    b2 = -0.152937760910536E+00;
                    b3 = 0.287578473937646E+00;
                    b4 = -0.221692883072871E-02;
                    b5 = 0.437335458190621E-01;
                    b6 = -0.958594913125490E-01;
                    b7 = 0.738976276909568E-03;
                    b8 = 0.958594913125492E-01;
               }
               else if (n == 6 && option == 1)
               {
                    u = 0.202018287045609E+01;
                    v = 0.958572464613819E+00;
                    b0 = -0.361840434143098E+02;
                    b1 = -0.447936529138517E+00;
                    b2 = 0.112077863004144E+02;
                    b3 = 0.392940404320855E-02;
                    b4 = -0.254859786784158E+01;
                    b5 = 0.775156917007496E-01;
                    b6 = -0.130980134773619E-02;
                    b7 = 0.509719573568315E+00;
                    b8 = 0.436600449245395E-03;
               }
               else if (n == 6 && option == 2)
               {
                    u = 0.958572464613819E+00;
                    v = 0.202018287045609E+01;
                    b0 = 0.448873836333650E+01;
                    b1 = -0.238473566140736E+01;
                    b2 = -0.413008493198885E+00;
                    b3 = 0.152915872070494E+01;
                    b4 = -0.654900673868093E-02;
                    b5 = 0.775156917007496E-01;
                    b6 = -0.509719573568314E+00;
                    b7 = 0.130980134773618E-02;
                    b8 = 0.169906524522772E+00;
               }

               typeMethods.r8vec_zero(n * o, ref x);

               k = -1;
               //
               //  1 point.
               //
               k = k + 1;
               w[k] = b0;
               //
               //  2 * N points.
               //
               for (i = 0; i < n; i++)
               {
                    k = k + 1;
                    x[i + k * n] = -u;
                    w[k] = b1;
                    k = k + 1;
                    x[i + k * n] = +u;
                    w[k] = b1;
               }

               //
               //  2 * N points.
               //
               for (i = 0; i < n; i++)
               {
                    k = k + 1;
                    x[i + k * n] = -v;
                    w[k] = b2;
                    k = k + 1;
                    x[i + k * n] = +v;
                    w[k] = b2;
               }

               //
               //  4 * ( N * ( N - 1 ) / 2 ) points.
               //
               for (i = 0; i < n - 1; i++)
               {
                    for (j = i + 1; j < n; j++)
                    {
                         k = k + 1;
                         x[i + k * n] = -u;
                         x[j + k * n] = -u;
                         w[k] = b3;
                         k = k + 1;
                         x[i + k * n] = -u;
                         x[j + k * n] = +u;
                         w[k] = b3;
                         k = k + 1;
                         x[i + k * n] = +u;
                         x[j + k * n] = -u;
                         w[k] = b3;
                         k = k + 1;
                         x[i + k * n] = +u;
                         x[j + k * n] = +u;
                         w[k] = b3;
                    }
               }

               //
               //  4 * ( N * ( N - 1 ) / 2 ) points.
               //
               for (i = 0; i < n - 1; i++)
               {
                    for (j = i + 1; j < n; j++)
                    {
                         k = k + 1;
                         x[i + k * n] = -v;
                         x[j + k * n] = -v;
                         w[k] = b4;
                         k = k + 1;
                         x[i + k * n] = -v;
                         x[j + k * n] = +v;
                         w[k] = b4;
                         k = k + 1;
                         x[i + k * n] = +v;
                         x[j + k * n] = -v;
                         w[k] = b4;
                         k = k + 1;
                         x[i + k * n] = +v;
                         x[j + k * n] = +v;
                         w[k] = b4;
                    }
               }

               //
               //  4 * ( N * ( N - 1 ) ) points.
               //
               for (i = 0; i < n - 1; i++)
               {
                    for (j = i + 1; j < n; j++)
                    {
                         k = k + 1;
                         x[i + k * n] = -u;
                         x[j + k * n] = -v;
                         w[k] = b5;
                         k = k + 1;
                         x[i + k * n] = -u;
                         x[j + k * n] = +v;
                         w[k] = b5;
                         k = k + 1;
                         x[i + k * n] = +u;
                         x[j + k * n] = -v;
                         w[k] = b5;
                         k = k + 1;
                         x[i + k * n] = +u;
                         x[j + k * n] = +v;
                         w[k] = b5;
                         k = k + 1;
                         x[i + k * n] = -v;
                         x[j + k * n] = -u;
                         w[k] = b5;
                         k = k + 1;
                         x[i + k * n] = -v;
                         x[j + k * n] = +u;
                         w[k] = b5;
                         k = k + 1;
                         x[i + k * n] = +v;
                         x[j + k * n] = -u;
                         w[k] = b5;
                         k = k + 1;
                         x[i + k * n] = +v;
                         x[j + k * n] = +u;
                         w[k] = b5;
                    }
               }

               //
               //  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
               //
               for (i = 0; i < n - 2; i++)
               {
                    for (j = i + 1; j < n - 1; j++)
                    {
                         for (l = j + 1; l < n; l++)
                         {
                              k = k + 1;
                              x[i + k * n] = -u;
                              x[j + k * n] = -u;
                              x[l + k * n] = -u;
                              w[k] = b6;
                              k = k + 1;
                              x[i + k * n] = -u;
                              x[j + k * n] = -u;
                              x[l + k * n] = +u;
                              w[k] = b6;
                              k = k + 1;
                              x[i + k * n] = -u;
                              x[j + k * n] = +u;
                              x[l + k * n] = -u;
                              w[k] = b6;
                              k = k + 1;
                              x[i + k * n] = -u;
                              x[j + k * n] = +u;
                              x[l + k * n] = +u;
                              w[k] = b6;
                              k = k + 1;
                              x[i + k * n] = +u;
                              x[j + k * n] = -u;
                              x[l + k * n] = -u;
                              w[k] = b6;
                              k = k + 1;
                              x[i + k * n] = +u;
                              x[j + k * n] = -u;
                              x[l + k * n] = +u;
                              w[k] = b6;
                              k = k + 1;
                              x[i + k * n] = +u;
                              x[j + k * n] = +u;
                              x[l + k * n] = -u;
                              w[k] = b6;
                              k = k + 1;
                              x[i + k * n] = +u;
                              x[j + k * n] = +u;
                              x[l + k * n] = +u;
                              w[k] = b6;
                         }
                    }
               }

               //
               //  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
               //
               for (i = 0; i < n - 2; i++)
               {
                    for (j = i + 1; j < n - 1; j++)
                    {
                         for (l = j + 1; l < n; l++)
                         {
                              k = k + 1;
                              x[i + k * n] = -v;
                              x[j + k * n] = -v;
                              x[l + k * n] = -v;
                              w[k] = b7;
                              k = k + 1;
                              x[i + k * n] = -v;
                              x[j + k * n] = -v;
                              x[l + k * n] = +v;
                              w[k] = b7;
                              k = k + 1;
                              x[i + k * n] = -v;
                              x[j + k * n] = +v;
                              x[l + k * n] = -v;
                              w[k] = b7;
                              k = k + 1;
                              x[i + k * n] = -v;
                              x[j + k * n] = +v;
                              x[l + k * n] = +v;
                              w[k] = b7;
                              k = k + 1;
                              x[i + k * n] = +v;
                              x[j + k * n] = -v;
                              x[l + k * n] = -v;
                              w[k] = b7;
                              k = k + 1;
                              x[i + k * n] = +v;
                              x[j + k * n] = -v;
                              x[l + k * n] = +v;
                              w[k] = b7;
                              k = k + 1;
                              x[i + k * n] = +v;
                              x[j + k * n] = +v;
                              x[l + k * n] = -v;
                              w[k] = b7;
                              k = k + 1;
                              x[i + k * n] = +v;
                              x[j + k * n] = +v;
                              x[l + k * n] = +v;
                              w[k] = b7;
                         }
                    }
               }

               //
               //  16 * ( N * ( N - 1 ) * ( N - 2 ) * ( N - 3 ) / 24 ) points.
               //
               for (i = 0; i < n - 3; i++)
               {
                    for (j = i + 1; j < n - 2; j++)
                    {
                         for (l = j + 1; l < n - 1; l++)
                         {
                              for (m = l + 1; m < n; m++)
                              {
                                   k = k + 1;
                                   x[i + k * n] = -u;
                                   x[j + k * n] = -u;
                                   x[l + k * n] = -u;
                                   x[m + k * n] = -u;
                                   w[k] = b8;
                                   k = k + 1;
                                   x[i + k * n] = -u;
                                   x[j + k * n] = -u;
                                   x[l + k * n] = -u;
                                   x[m + k * n] = +u;
                                   w[k] = b8;
                                   k = k + 1;
                                   x[i + k * n] = -u;
                                   x[j + k * n] = -u;
                                   x[l + k * n] = +u;
                                   x[m + k * n] = -u;
                                   w[k] = b8;
                                   k = k + 1;
                                   x[i + k * n] = -u;
                                   x[j + k * n] = -u;
                                   x[l + k * n] = +u;
                                   x[m + k * n] = +u;
                                   w[k] = b8;
                                   k = k + 1;
                                   x[i + k * n] = -u;
                                   x[j + k * n] = +u;
                                   x[l + k * n] = -u;
                                   x[m + k * n] = -u;
                                   w[k] = b8;
                                   k = k + 1;
                                   x[i + k * n] = -u;
                                   x[j + k * n] = +u;
                                   x[l + k * n] = -u;
                                   x[m + k * n] = +u;
                                   w[k] = b8;
                                   k = k + 1;
                                   x[i + k * n] = -u;
                                   x[j + k * n] = +u;
                                   x[l + k * n] = +u;
                                   x[m + k * n] = -u;
                                   w[k] = b8;
                                   k = k + 1;
                                   x[i + k * n] = -u;
                                   x[j + k * n] = +u;
                                   x[l + k * n] = +u;
                                   x[m + k * n] = +u;
                                   w[k] = b8;
                                   k = k + 1;
                                   x[i + k * n] = +u;
                                   x[j + k * n] = -u;
                                   x[l + k * n] = -u;
                                   x[m + k * n] = -u;
                                   w[k] = b8;
                                   k = k + 1;
                                   x[i + k * n] = +u;
                                   x[j + k * n] = -u;
                                   x[l + k * n] = -u;
                                   x[m + k * n] = +u;
                                   w[k] = b8;
                                   k = k + 1;
                                   x[i + k * n] = +u;
                                   x[j + k * n] = -u;
                                   x[l + k * n] = +u;
                                   x[m + k * n] = -u;
                                   w[k] = b8;
                                   k = k + 1;
                                   x[i + k * n] = +u;
                                   x[j + k * n] = -u;
                                   x[l + k * n] = +u;
                                   x[m + k * n] = +u;
                                   w[k] = b8;
                                   k = k + 1;
                                   x[i + k * n] = +u;
                                   x[j + k * n] = +u;
                                   x[l + k * n] = -u;
                                   x[m + k * n] = -u;
                                   w[k] = b8;
                                   k = k + 1;
                                   x[i + k * n] = +u;
                                   x[j + k * n] = +u;
                                   x[l + k * n] = -u;
                                   x[m + k * n] = +u;
                                   w[k] = b8;
                                   k = k + 1;
                                   x[i + k * n] = +u;
                                   x[j + k * n] = +u;
                                   x[l + k * n] = +u;
                                   x[m + k * n] = -u;
                                   w[k] = b8;
                                   k = k + 1;
                                   x[i + k * n] = +u;
                                   x[j + k * n] = +u;
                                   x[l + k * n] = +u;
                                   x[m + k * n] = +u;
                                   w[k] = b8;
                              }
                         }
                    }
               }

               return;
          }

          public static int en_r2_09_1_size(int n)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_R2_09_1_SIZE sizes the Stroud rule 9.1 for region EN_R2.
               //
               //  Discussion:
               //
               //    The rule has order O = ( 2 * N^4 - 4 * N^3 + 22 * N^2 - 8 * N + 3 ) / 3.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    23 January 2010
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
               //    Input, int N, the spatial dimension.
               //
               //    Output, int EN_R2_09_1_SIZE, the order.
               //
          {
               int o;

               o = (2 * (int)Math.Pow(n, 4)
                    - 4 * (int)Math.Pow(n, 3)
                    + 22 * (int)Math.Pow(n, 2)
                    - 8 * n
                    + 3) / 3;

               return o;
          }

          public static void en_r2_11_1(int n, int option, int o, ref double[] x, ref double[] w)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_R2_11_1 implements the Stroud rule 11.1 for region EN_R2.
               //
               //  Discussion:
               //
               //    The rule has order 
               //
               //      O = ( 4 * N^5 - 20 * N^4 + 140 * N^3 - 130 * N^2 + 96 * N + 15 ) / 15.
               //
               //    The rule has precision P = 11.
               //
               //    EN_R2 is the entire N-dimensional space with weight function
               //
               //      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 ) 
               //
               //    There are two versions of each rule, chosen by setting the
               //    OPTION variable to 1 or 2.
               //
               //    The rule as tabulated by Stenger is available for N = 2 through 20.
               //    This function accepts N = 3 through 5.
               //
               //     N    O
               //    __  ___
               //     3  151
               //     4  417
               //     5  983
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    24 January 2010
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
               //    Input, int N, the spatial dimension.
               //    3 <= N <= 5.
               //
               //    Input, int OPTION, chooses rule option 1 or 2.
               //
               //    Input, int O, the order.
               //
               //    Output, double X[N*O], the abscissas.
               //
               //    Output, double W[O], the weights.
               //
          {
               double b0 = 0;
               double b1 = 0;
               double b2 = 0;
               double b3 = 0;
               double b4 = 0;
               double b5 = 0;
               double b6 = 0;
               double b7 = 0;
               double b8 = 0;
               double b9 = 0;
               double b10 = 0;
               double b11 = 0;
               double b12 = 0;
               double b13 = 0;
               double b14 = 0;
               double b15 = 0;
               int i = 0;
               int i1 = 0;
               int i2 = 0;
               int i3 = 0;
               int i4 = 0;
               int i5 = 0;
               int j = 0;
               int k = 0;
               int l = 0;
               int m = 0;
               double u = 0;
               double v = 0;
               double w2 = 0;

               if (n < 3 || 5 < n)
               {
                    Console.WriteLine("");
                    Console.WriteLine("EN_R2_11_1 - Fatal error!");
                    Console.WriteLine("  3 <= N <= 5 required.");
                    return;
               }

               if (option < 1 || 2 < option)
               {
                    Console.WriteLine("");
                    Console.WriteLine("EN_R2_11_1 - Fatal error!");
                    Console.WriteLine("  1 <= OPTION <= 2 required.");
                    return;
               }

               if (n == 3 && option == 1)
               {
                    u = 0.235060497367449E+01;
                    v = 0.436077411927617E+00;
                    w2 = 0.133584907401370E+01;
                    b0 = -0.881591029957858E+01;
                    b1 = -0.751996143360650E-01;
                    b2 = 0.621743189471515E+01;
                    b3 = 0.241426451456494E+00;
                    b4 = -0.120709739276065E-02;
                    b5 = -0.427751221210138E+01;
                    b6 = 0.550169924840163E-01;
                    b7 = 0.237084999634707E-01;
                    b8 = -0.169791992887741E-02;
                    b9 = -0.252266276123350E-04;
                    b10 = 0.326777873717691E+01;
                    b11 = 0.968469949206802E-02;
                    b12 = 0.789754514877422E-03;
                    b13 = 0.000000000000000E+00;
                    b14 = 0.000000000000000E+00;
                    b15 = 0.000000000000000E+00;
               }
               else if (n == 3 && option == 2)
               {
                    u = 0.235060497367449E+01;
                    v = 0.133584907401370E+01;
                    w2 = 0.436077411927617E+00;
                    b0 = -0.141214037032900E+02;
                    b1 = -0.803730274707282E-01;
                    b2 = 0.235546545595906E+00;
                    b3 = 0.888123191556611E+01;
                    b4 = 0.142467131155533E-03;
                    b5 = 0.582993124006494E-01;
                    b6 = -0.561099173155661E+01;
                    b7 = -0.204028691521686E-02;
                    b8 = 0.252880089932256E-01;
                    b9 = -0.814378678627283E-04;
                    b10 = 0.804353953375146E-02;
                    b11 = 0.393451849690453E+01;
                    b12 = 0.171183493169724E-03;
                    b13 = 0.000000000000000E+00;
                    b14 = 0.000000000000000E+00;
                    b15 = 0.000000000000000E+00;
               }
               else if (n == 4 && option == 1)
               {
                    u = 0.235060497367449E+01;
                    v = 0.436077411927617E+00;
                    w2 = 0.133584907401370E+01;
                    b0 = 0.241502736147339E+03;
                    b1 = -0.196095938531478E+00;
                    b2 = -0.128675737999280E+03;
                    b3 = 0.307568784278696E+00;
                    b4 = -0.480908422319460E-02;
                    b5 = 0.698087019367085E+02;
                    b6 = 0.631837143743771E-01;
                    b7 = 0.392226151971179E-01;
                    b8 = -0.300948471646799E-02;
                    b9 = -0.650235306755170E-04;
                    b10 = -0.386951974646715E+02;
                    b11 = 0.171656829095787E-01;
                    b12 = 0.139980343116450E-02;
                    b13 = 0.101552487093372E-04;
                    b14 = 0.222435922356439E+02;
                    b15 = 0.000000000000000E+00;
               }
               else if (n == 4 && option == 2)
               {
                    u = 0.235060497367449E+01;
                    v = 0.133584907401370E+01;
                    w2 = 0.436077411927617E+00;
                    b0 = -0.151944464736584E+03;
                    b1 = -0.223498438689039E+00;
                    b2 = 0.243574919068010E+00;
                    b3 = 0.634373877008693E+02;
                    b4 = -0.782065187814018E-04;
                    b5 = 0.911833754536616E-01;
                    b6 = -0.238927288245914E+02;
                    b7 = -0.422314408318853E-02;
                    b8 = 0.448218289217760E-01;
                    b9 = -0.138053374667391E-03;
                    b10 = 0.607473265800655E-02;
                    b11 = 0.697375246129742E+01;
                    b12 = 0.303414841680135E-03;
                    b13 = -0.314574391771792E-05;
                    b14 = 0.409103498175100E-02;
                    b15 = 0.000000000000000E+00;
               }
               else if (n == 5 && option == 1)
               {
                    u = 0.235060497367449E+01;
                    v = 0.436077411927617E+00;
                    w2 = 0.133584907401370E+01;
                    b0 = 0.255885269311763E+04;
                    b1 = -0.439598677491526E+00;
                    b2 = -0.106541406144610E+04;
                    b3 = 0.453540909054264E+00;
                    b4 = -0.132100905623778E-01;
                    b5 = 0.418606568954203E+03;
                    b6 = 0.511394563043680E-01;
                    b7 = 0.645581013845604E-01;
                    b8 = -0.533417277494500E-02;
                    b9 = -0.137981626254496E-03;
                    b10 = -0.147436933189884E+03;
                    b11 = 0.304253807765057E-01;
                    b12 = 0.248108698207828E-02;
                    b13 = 0.113652094546015E-04;
                    b14 = 0.394257407160391E+02;
                    b15 = 0.331725011358320E-05;
               }
               else if (n == 5 && option == 2)
               {
                    u = 0.235060497367449E+01;
                    v = 0.133584907401370E+01;
                    w2 = 0.436077411927617E+00;
                    b0 = -0.761305347548192E+03;
                    b1 = -0.536360805019297E+00;
                    b2 = 0.110669832078736E+00;
                    b3 = 0.246421088923968E+03;
                    b4 = -0.773649327968607E-03;
                    b5 = 0.169088641205970E+00;
                    b6 = -0.670700680243651E+02;
                    b7 = -0.856090560229205E-02;
                    b8 = 0.794446232770302E-01;
                    b9 = -0.220272863263544E-03;
                    b10 = -0.373515812228225E-02;
                    b11 = 0.123606544052884E+02;
                    b12 = 0.537788804557843E-03;
                    b13 = -0.122101861480881E-04;
                    b14 = 0.725117070759373E-02;
                    b15 = 0.331725011358320E-05;
               }

               typeMethods.r8vec_zero(n * o, ref x);

               k = -1;
               //
               //  1 point.
               //
               k = k + 1;
               w[k] = b0;
               //
               //  2 * N points.
               //
               for (i = 0; i < n; i++)
               {
                    k = k + 1;
                    x[i + k * n] = -u;
                    w[k] = b1;
                    k = k + 1;
                    x[i + k * n] = +u;
                    w[k] = b1;
               }

               //
               //  2 * N points.
               //
               for (i = 0; i < n; i++)
               {
                    k = k + 1;
                    x[i + k * n] = -v;
                    w[k] = b2;
                    k = k + 1;
                    x[i + k * n] = +v;
                    w[k] = b2;
               }

               //
               //  2 * N points.
               //
               for (i = 0; i < n; i++)
               {
                    k = k + 1;
                    x[i + k * n] = -w2;
                    w[k] = b3;
                    k = k + 1;
                    x[i + k * n] = +w2;
                    w[k] = b3;
               }

               //
               //  4 * ( N * ( N - 1 ) / 2 ) points.
               //
               for (i = 0; i < n - 1; i++)
               {
                    for (j = i + 1; j < n; j++)
                    {
                         k = k + 1;
                         x[i + k * n] = -u;
                         x[j + k * n] = -u;
                         w[k] = b4;
                         k = k + 1;
                         x[i + k * n] = -u;
                         x[j + k * n] = +u;
                         w[k] = b4;
                         k = k + 1;
                         x[i + k * n] = +u;
                         x[j + k * n] = -u;
                         w[k] = b4;
                         k = k + 1;
                         x[i + k * n] = +u;
                         x[j + k * n] = +u;
                         w[k] = b4;
                    }
               }

               //
               //  4 * ( N * ( N - 1 ) / 2 ) points.
               //
               for (i = 0; i < n - 1; i++)
               {
                    for (j = i + 1; j < n; j++)
                    {
                         k = k + 1;
                         x[i + k * n] = -v;
                         x[j + k * n] = -v;
                         w[k] = b5;
                         k = k + 1;
                         x[i + k * n] = -v;
                         x[j + k * n] = +v;
                         w[k] = b5;
                         k = k + 1;
                         x[i + k * n] = +v;
                         x[j + k * n] = -v;
                         w[k] = b5;
                         k = k + 1;
                         x[i + k * n] = +v;
                         x[j + k * n] = +v;
                         w[k] = b5;
                    }
               }

               //
               //  4 * ( N * ( N - 1 ) / 2 ) points.
               //
               for (i = 0; i < n - 1; i++)
               {
                    for (j = i + 1; j < n; j++)
                    {
                         k = k + 1;
                         x[i + k * n] = -w2;
                         x[j + k * n] = -w2;
                         w[k] = b6;
                         k = k + 1;
                         x[i + k * n] = -w2;
                         x[j + k * n] = +w2;
                         w[k] = b6;
                         k = k + 1;
                         x[i + k * n] = +w2;
                         x[j + k * n] = -w2;
                         w[k] = b6;
                         k = k + 1;
                         x[i + k * n] = +w2;
                         x[j + k * n] = +w2;
                         w[k] = b6;
                    }
               }

               //
               //  4 * ( N * ( N - 1 ) ) points.
               //
               for (i = 0; i < n - 1; i++)
               {
                    for (j = i + 1; j < n; j++)
                    {
                         k = k + 1;
                         x[i + k * n] = -u;
                         x[j + k * n] = -v;
                         w[k] = b7;
                         k = k + 1;
                         x[i + k * n] = -u;
                         x[j + k * n] = +v;
                         w[k] = b7;
                         k = k + 1;
                         x[i + k * n] = +u;
                         x[j + k * n] = -v;
                         w[k] = b7;
                         k = k + 1;
                         x[i + k * n] = +u;
                         x[j + k * n] = +v;
                         w[k] = b7;
                         k = k + 1;
                         x[i + k * n] = -v;
                         x[j + k * n] = -u;
                         w[k] = b7;
                         k = k + 1;
                         x[i + k * n] = -v;
                         x[j + k * n] = +u;
                         w[k] = b7;
                         k = k + 1;
                         x[i + k * n] = +v;
                         x[j + k * n] = -u;
                         w[k] = b7;
                         k = k + 1;
                         x[i + k * n] = +v;
                         x[j + k * n] = +u;
                         w[k] = b7;
                    }
               }

               //
               //  4 * ( N * ( N - 1 ) ) points.
               //
               for (i = 0; i < n - 1; i++)
               {
                    for (j = i + 1; j < n; j++)
                    {
                         k = k + 1;
                         x[i + k * n] = -u;
                         x[j + k * n] = -w2;
                         w[k] = b8;
                         k = k + 1;
                         x[i + k * n] = -u;
                         x[j + k * n] = +w2;
                         w[k] = b8;
                         k = k + 1;
                         x[i + k * n] = +u;
                         x[j + k * n] = -w2;
                         w[k] = b8;
                         k = k + 1;
                         x[i + k * n] = +u;
                         x[j + k * n] = +w2;
                         w[k] = b8;
                         k = k + 1;
                         x[i + k * n] = -w2;
                         x[j + k * n] = -u;
                         w[k] = b8;
                         k = k + 1;
                         x[i + k * n] = -w2;
                         x[j + k * n] = +u;
                         w[k] = b8;
                         k = k + 1;
                         x[i + k * n] = +w2;
                         x[j + k * n] = -u;
                         w[k] = b8;
                         k = k + 1;
                         x[i + k * n] = +w2;
                         x[j + k * n] = +u;
                         w[k] = b8;
                    }
               }

               //
               //  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
               //
               for (i = 0; i < n - 2; i++)
               {
                    for (j = i + 1; j < n - 1; j++)
                    {
                         for (l = j + 1; l < n; l++)
                         {
                              k = k + 1;
                              x[i + k * n] = -u;
                              x[j + k * n] = -u;
                              x[l + k * n] = -u;
                              w[k] = b9;
                              k = k + 1;
                              x[i + k * n] = -u;
                              x[j + k * n] = -u;
                              x[l + k * n] = +u;
                              w[k] = b9;
                              k = k + 1;
                              x[i + k * n] = -u;
                              x[j + k * n] = +u;
                              x[l + k * n] = -u;
                              w[k] = b9;
                              k = k + 1;
                              x[i + k * n] = -u;
                              x[j + k * n] = +u;
                              x[l + k * n] = +u;
                              w[k] = b9;
                              k = k + 1;
                              x[i + k * n] = +u;
                              x[j + k * n] = -u;
                              x[l + k * n] = -u;
                              w[k] = b9;
                              k = k + 1;
                              x[i + k * n] = +u;
                              x[j + k * n] = -u;
                              x[l + k * n] = +u;
                              w[k] = b9;
                              k = k + 1;
                              x[i + k * n] = +u;
                              x[j + k * n] = +u;
                              x[l + k * n] = -u;
                              w[k] = b9;
                              k = k + 1;
                              x[i + k * n] = +u;
                              x[j + k * n] = +u;
                              x[l + k * n] = +u;
                              w[k] = b9;
                         }
                    }
               }

               //
               //  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
               //
               for (i = 0; i < n - 2; i++)
               {
                    for (j = i + 1; j < n - 1; j++)
                    {
                         for (l = j + 1; l < n; l++)
                         {
                              k = k + 1;
                              x[i + k * n] = -v;
                              x[j + k * n] = -v;
                              x[l + k * n] = -v;
                              w[k] = b10;
                              k = k + 1;
                              x[i + k * n] = -v;
                              x[j + k * n] = -v;
                              x[l + k * n] = +v;
                              w[k] = b10;
                              k = k + 1;
                              x[i + k * n] = -v;
                              x[j + k * n] = +v;
                              x[l + k * n] = -v;
                              w[k] = b10;
                              k = k + 1;
                              x[i + k * n] = -v;
                              x[j + k * n] = +v;
                              x[l + k * n] = +v;
                              w[k] = b10;
                              k = k + 1;
                              x[i + k * n] = +v;
                              x[j + k * n] = -v;
                              x[l + k * n] = -v;
                              w[k] = b10;
                              k = k + 1;
                              x[i + k * n] = +v;
                              x[j + k * n] = -v;
                              x[l + k * n] = +v;
                              w[k] = b10;
                              k = k + 1;
                              x[i + k * n] = +v;
                              x[j + k * n] = +v;
                              x[l + k * n] = -v;
                              w[k] = b10;
                              k = k + 1;
                              x[i + k * n] = +v;
                              x[j + k * n] = +v;
                              x[l + k * n] = +v;
                              w[k] = b10;
                         }
                    }
               }

               //
               //  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
               //
               for (i = 0; i < n - 2; i++)
               {
                    for (j = i + 1; j < n - 1; j++)
                    {
                         for (l = j + 1; l < n; l++)
                         {
                              k = k + 1;
                              x[i + k * n] = -w2;
                              x[j + k * n] = -w2;
                              x[l + k * n] = -w2;
                              w[k] = b11;
                              k = k + 1;
                              x[i + k * n] = -w2;
                              x[j + k * n] = -w2;
                              x[l + k * n] = +w2;
                              w[k] = b11;
                              k = k + 1;
                              x[i + k * n] = -w2;
                              x[j + k * n] = +w2;
                              x[l + k * n] = -w2;
                              w[k] = b11;
                              k = k + 1;
                              x[i + k * n] = -w2;
                              x[j + k * n] = +w2;
                              x[l + k * n] = +w2;
                              w[k] = b11;
                              k = k + 1;
                              x[i + k * n] = +w2;
                              x[j + k * n] = -w2;
                              x[l + k * n] = -w2;
                              w[k] = b11;
                              k = k + 1;
                              x[i + k * n] = +w2;
                              x[j + k * n] = -w2;
                              x[l + k * n] = +w2;
                              w[k] = b11;
                              k = k + 1;
                              x[i + k * n] = +w2;
                              x[j + k * n] = +w2;
                              x[l + k * n] = -w2;
                              w[k] = b11;
                              k = k + 1;
                              x[i + k * n] = +w2;
                              x[j + k * n] = +w2;
                              x[l + k * n] = +w2;
                              w[k] = b11;
                         }
                    }
               }

               //
               //  8 * ( N * ( N - 1 ) * ( N - 2 ) / 2 ) points.
               //
               for (i = 0; i < n - 2; i++)
               {
                    for (j = i + 1; j < n - 1; j++)
                    {
                         for (l = j + 1; l < n; l++)
                         {
                              k = k + 1;
                              x[i + k * n] = -u;
                              x[j + k * n] = -u;
                              x[l + k * n] = -v;
                              w[k] = b12;
                              k = k + 1;
                              x[i + k * n] = -u;
                              x[j + k * n] = -u;
                              x[l + k * n] = +v;
                              w[k] = b12;
                              k = k + 1;
                              x[i + k * n] = -u;
                              x[j + k * n] = +u;
                              x[l + k * n] = -v;
                              w[k] = b12;
                              k = k + 1;
                              x[i + k * n] = -u;
                              x[j + k * n] = +u;
                              x[l + k * n] = +v;
                              w[k] = b12;
                              k = k + 1;
                              x[i + k * n] = +u;
                              x[j + k * n] = -u;
                              x[l + k * n] = -v;
                              w[k] = b12;
                              k = k + 1;
                              x[i + k * n] = +u;
                              x[j + k * n] = -u;
                              x[l + k * n] = +v;
                              w[k] = b12;
                              k = k + 1;
                              x[i + k * n] = +u;
                              x[j + k * n] = +u;
                              x[l + k * n] = -v;
                              w[k] = b12;
                              k = k + 1;
                              x[i + k * n] = +u;
                              x[j + k * n] = +u;
                              x[l + k * n] = +v;
                              w[k] = b12;
                              k = k + 1;
                              x[i + k * n] = -u;
                              x[j + k * n] = -v;
                              x[l + k * n] = -u;
                              w[k] = b12;
                              k = k + 1;
                              x[i + k * n] = -u;
                              x[j + k * n] = -v;
                              x[l + k * n] = +u;
                              w[k] = b12;
                              k = k + 1;
                              x[i + k * n] = -u;
                              x[j + k * n] = +v;
                              x[l + k * n] = -u;
                              w[k] = b12;
                              k = k + 1;
                              x[i + k * n] = -u;
                              x[j + k * n] = +v;
                              x[l + k * n] = +u;
                              w[k] = b12;
                              k = k + 1;
                              x[i + k * n] = +u;
                              x[j + k * n] = -v;
                              x[l + k * n] = -u;
                              w[k] = b12;
                              k = k + 1;
                              x[i + k * n] = +u;
                              x[j + k * n] = -v;
                              x[l + k * n] = +u;
                              w[k] = b12;
                              k = k + 1;
                              x[i + k * n] = +u;
                              x[j + k * n] = +v;
                              x[l + k * n] = -u;
                              w[k] = b12;
                              k = k + 1;
                              x[i + k * n] = +u;
                              x[j + k * n] = +v;
                              x[l + k * n] = +u;
                              w[k] = b12;
                              k = k + 1;
                              x[i + k * n] = -v;
                              x[j + k * n] = -u;
                              x[l + k * n] = -u;
                              w[k] = b12;
                              k = k + 1;
                              x[i + k * n] = -v;
                              x[j + k * n] = -u;
                              x[l + k * n] = +u;
                              w[k] = b12;
                              k = k + 1;
                              x[i + k * n] = -v;
                              x[j + k * n] = +u;
                              x[l + k * n] = -u;
                              w[k] = b12;
                              k = k + 1;
                              x[i + k * n] = -v;
                              x[j + k * n] = +u;
                              x[l + k * n] = +u;
                              w[k] = b12;
                              k = k + 1;
                              x[i + k * n] = +v;
                              x[j + k * n] = -u;
                              x[l + k * n] = -u;
                              w[k] = b12;
                              k = k + 1;
                              x[i + k * n] = +v;
                              x[j + k * n] = -u;
                              x[l + k * n] = +u;
                              w[k] = b12;
                              k = k + 1;
                              x[i + k * n] = +v;
                              x[j + k * n] = +u;
                              x[l + k * n] = -u;
                              w[k] = b12;
                              k = k + 1;
                              x[i + k * n] = +v;
                              x[j + k * n] = +u;
                              x[l + k * n] = +u;
                              w[k] = b12;
                         }
                    }
               }

               //
               //  16 * ( N * ( N - 1 ) * ( N - 2 ) * ( N - 3 ) / 24 ) points.
               //
               for (i = 0; i < n - 3; i++)
               {
                    for (j = i + 1; j < n - 2; j++)
                    {
                         for (l = j + 1; l < n - 1; l++)
                         {
                              for (m = l + 1; m < n; m++)
                              {
                                   k = k + 1;
                                   x[i + k * n] = -u;
                                   x[j + k * n] = -u;
                                   x[l + k * n] = -u;
                                   x[m + k * n] = -u;
                                   w[k] = b13;
                                   k = k + 1;
                                   x[i + k * n] = -u;
                                   x[j + k * n] = -u;
                                   x[l + k * n] = -u;
                                   x[m + k * n] = +u;
                                   w[k] = b13;
                                   k = k + 1;
                                   x[i + k * n] = -u;
                                   x[j + k * n] = -u;
                                   x[l + k * n] = +u;
                                   x[m + k * n] = -u;
                                   w[k] = b13;
                                   k = k + 1;
                                   x[i + k * n] = -u;
                                   x[j + k * n] = -u;
                                   x[l + k * n] = +u;
                                   x[m + k * n] = +u;
                                   w[k] = b13;
                                   k = k + 1;
                                   x[i + k * n] = -u;
                                   x[j + k * n] = +u;
                                   x[l + k * n] = -u;
                                   x[m + k * n] = -u;
                                   w[k] = b13;
                                   k = k + 1;
                                   x[i + k * n] = -u;
                                   x[j + k * n] = +u;
                                   x[l + k * n] = -u;
                                   x[m + k * n] = +u;
                                   w[k] = b13;
                                   k = k + 1;
                                   x[i + k * n] = -u;
                                   x[j + k * n] = +u;
                                   x[l + k * n] = +u;
                                   x[m + k * n] = -u;
                                   w[k] = b13;
                                   k = k + 1;
                                   x[i + k * n] = -u;
                                   x[j + k * n] = +u;
                                   x[l + k * n] = +u;
                                   x[m + k * n] = +u;
                                   w[k] = b13;
                                   k = k + 1;
                                   x[i + k * n] = +u;
                                   x[j + k * n] = -u;
                                   x[l + k * n] = -u;
                                   x[m + k * n] = -u;
                                   w[k] = b13;
                                   k = k + 1;
                                   x[i + k * n] = +u;
                                   x[j + k * n] = -u;
                                   x[l + k * n] = -u;
                                   x[m + k * n] = +u;
                                   w[k] = b13;
                                   k = k + 1;
                                   x[i + k * n] = +u;
                                   x[j + k * n] = -u;
                                   x[l + k * n] = +u;
                                   x[m + k * n] = -u;
                                   w[k] = b13;
                                   k = k + 1;
                                   x[i + k * n] = +u;
                                   x[j + k * n] = -u;
                                   x[l + k * n] = +u;
                                   x[m + k * n] = +u;
                                   w[k] = b13;
                                   k = k + 1;
                                   x[i + k * n] = +u;
                                   x[j + k * n] = +u;
                                   x[l + k * n] = -u;
                                   x[m + k * n] = -u;
                                   w[k] = b13;
                                   k = k + 1;
                                   x[i + k * n] = +u;
                                   x[j + k * n] = +u;
                                   x[l + k * n] = -u;
                                   x[m + k * n] = +u;
                                   w[k] = b13;
                                   k = k + 1;
                                   x[i + k * n] = +u;
                                   x[j + k * n] = +u;
                                   x[l + k * n] = +u;
                                   x[m + k * n] = -u;
                                   w[k] = b13;
                                   k = k + 1;
                                   x[i + k * n] = +u;
                                   x[j + k * n] = +u;
                                   x[l + k * n] = +u;
                                   x[m + k * n] = +u;
                                   w[k] = b13;
                              }
                         }
                    }
               }

               //
               //  16 * ( N * ( N - 1 ) * ( N - 2 ) * ( N - 3 ) / 24 ) points.
               //
               for (i = 0; i < n - 3; i++)
               {
                    for (j = i + 1; j < n - 2; j++)
                    {
                         for (l = j + 1; l < n - 1; l++)
                         {
                              for (m = l + 1; m < n; m++)
                              {
                                   k = k + 1;
                                   x[i + k * n] = -v;
                                   x[j + k * n] = -v;
                                   x[l + k * n] = -v;
                                   x[m + k * n] = -v;
                                   w[k] = b14;
                                   k = k + 1;
                                   x[i + k * n] = -v;
                                   x[j + k * n] = -v;
                                   x[l + k * n] = -v;
                                   x[m + k * n] = +v;
                                   w[k] = b14;
                                   k = k + 1;
                                   x[i + k * n] = -v;
                                   x[j + k * n] = -v;
                                   x[l + k * n] = +v;
                                   x[m + k * n] = -v;
                                   w[k] = b14;
                                   k = k + 1;
                                   x[i + k * n] = -v;
                                   x[j + k * n] = -v;
                                   x[l + k * n] = +v;
                                   x[m + k * n] = +v;
                                   w[k] = b14;
                                   k = k + 1;
                                   x[i + k * n] = -v;
                                   x[j + k * n] = +v;
                                   x[l + k * n] = -v;
                                   x[m + k * n] = -v;
                                   w[k] = b14;
                                   k = k + 1;
                                   x[i + k * n] = -v;
                                   x[j + k * n] = +v;
                                   x[l + k * n] = -v;
                                   x[m + k * n] = +v;
                                   w[k] = b14;
                                   k = k + 1;
                                   x[i + k * n] = -v;
                                   x[j + k * n] = +v;
                                   x[l + k * n] = +v;
                                   x[m + k * n] = -v;
                                   w[k] = b14;
                                   k = k + 1;
                                   x[i + k * n] = -v;
                                   x[j + k * n] = +v;
                                   x[l + k * n] = +v;
                                   x[m + k * n] = +v;
                                   w[k] = b14;
                                   k = k + 1;
                                   x[i + k * n] = +v;
                                   x[j + k * n] = -v;
                                   x[l + k * n] = -v;
                                   x[m + k * n] = -v;
                                   w[k] = b14;
                                   k = k + 1;
                                   x[i + k * n] = +v;
                                   x[j + k * n] = -v;
                                   x[l + k * n] = -v;
                                   x[m + k * n] = +v;
                                   w[k] = b14;
                                   k = k + 1;
                                   x[i + k * n] = +v;
                                   x[j + k * n] = -v;
                                   x[l + k * n] = +v;
                                   x[m + k * n] = -v;
                                   w[k] = b14;
                                   k = k + 1;
                                   x[i + k * n] = +v;
                                   x[j + k * n] = -v;
                                   x[l + k * n] = +v;
                                   x[m + k * n] = +v;
                                   w[k] = b14;
                                   k = k + 1;
                                   x[i + k * n] = +v;
                                   x[j + k * n] = +v;
                                   x[l + k * n] = -v;
                                   x[m + k * n] = -v;
                                   w[k] = b14;
                                   k = k + 1;
                                   x[i + k * n] = +v;
                                   x[j + k * n] = +v;
                                   x[l + k * n] = -v;
                                   x[m + k * n] = +v;
                                   w[k] = b14;
                                   k = k + 1;
                                   x[i + k * n] = +v;
                                   x[j + k * n] = +v;
                                   x[l + k * n] = +v;
                                   x[m + k * n] = -v;
                                   w[k] = b14;
                                   k = k + 1;
                                   x[i + k * n] = +v;
                                   x[j + k * n] = +v;
                                   x[l + k * n] = +v;
                                   x[m + k * n] = +v;
                                   w[k] = b14;
                              }
                         }
                    }
               }

               //
               //  All quintuples UUUUU with 32 sign combinations.
               //
               for (i1 = 0; i1 < n - 4; i1++)
               {
                    for (i2 = i1 + 1; i2 < n - 3; i2++)
                    {
                         for (i3 = i2 + 1; i3 < n - 2; i3++)
                         {
                              for (i4 = i3 + 1; i4 < n - 1; i4++)
                              {
                                   for (i5 = i4 + 1; i5 < n; i5++)
                                   {
                                        k = k + 1;
                                        x[i1 + k * n] = -u;
                                        x[i2 + k * n] = -u;
                                        x[i3 + k * n] = -u;
                                        x[i4 + k * n] = -u;
                                        x[i5 + k * n] = -u;
                                        w[k] = b15;
                                        k = k + 1;
                                        x[i1 + k * n] = -u;
                                        x[i2 + k * n] = -u;
                                        x[i3 + k * n] = -u;
                                        x[i4 + k * n] = -u;
                                        x[i5 + k * n] = +u;
                                        w[k] = b15;
                                        k = k + 1;
                                        x[i1 + k * n] = -u;
                                        x[i2 + k * n] = -u;
                                        x[i3 + k * n] = -u;
                                        x[i4 + k * n] = +u;
                                        x[i5 + k * n] = -u;
                                        w[k] = b15;
                                        k = k + 1;
                                        x[i1 + k * n] = -u;
                                        x[i2 + k * n] = -u;
                                        x[i3 + k * n] = -u;
                                        x[i4 + k * n] = +u;
                                        x[i5 + k * n] = +u;
                                        w[k] = b15;
                                        k = k + 1;
                                        x[i1 + k * n] = -u;
                                        x[i2 + k * n] = -u;
                                        x[i3 + k * n] = +u;
                                        x[i4 + k * n] = -u;
                                        x[i5 + k * n] = -u;
                                        w[k] = b15;
                                        k = k + 1;
                                        x[i1 + k * n] = -u;
                                        x[i2 + k * n] = -u;
                                        x[i3 + k * n] = +u;
                                        x[i4 + k * n] = -u;
                                        x[i5 + k * n] = +u;
                                        w[k] = b15;
                                        k = k + 1;
                                        x[i1 + k * n] = -u;
                                        x[i2 + k * n] = -u;
                                        x[i3 + k * n] = +u;
                                        x[i4 + k * n] = +u;
                                        x[i5 + k * n] = -u;
                                        w[k] = b15;
                                        k = k + 1;
                                        x[i1 + k * n] = -u;
                                        x[i2 + k * n] = -u;
                                        x[i3 + k * n] = +u;
                                        x[i4 + k * n] = +u;
                                        x[i5 + k * n] = +u;
                                        w[k] = b15;
                                        k = k + 1;
                                        x[i1 + k * n] = -u;
                                        x[i2 + k * n] = +u;
                                        x[i3 + k * n] = -u;
                                        x[i4 + k * n] = -u;
                                        x[i5 + k * n] = -u;
                                        w[k] = b15;
                                        k = k + 1;
                                        x[i1 + k * n] = -u;
                                        x[i2 + k * n] = +u;
                                        x[i3 + k * n] = -u;
                                        x[i4 + k * n] = -u;
                                        x[i5 + k * n] = +u;
                                        w[k] = b15;
                                        k = k + 1;
                                        x[i1 + k * n] = -u;
                                        x[i2 + k * n] = +u;
                                        x[i3 + k * n] = -u;
                                        x[i4 + k * n] = +u;
                                        x[i5 + k * n] = -u;
                                        w[k] = b15;
                                        k = k + 1;
                                        x[i1 + k * n] = -u;
                                        x[i2 + k * n] = +u;
                                        x[i3 + k * n] = -u;
                                        x[i4 + k * n] = +u;
                                        x[i5 + k * n] = +u;
                                        w[k] = b15;
                                        k = k + 1;
                                        x[i1 + k * n] = -u;
                                        x[i2 + k * n] = +u;
                                        x[i3 + k * n] = +u;
                                        x[i4 + k * n] = -u;
                                        x[i5 + k * n] = -u;
                                        w[k] = b15;
                                        k = k + 1;
                                        x[i1 + k * n] = -u;
                                        x[i2 + k * n] = +u;
                                        x[i3 + k * n] = +u;
                                        x[i4 + k * n] = -u;
                                        x[i5 + k * n] = +u;
                                        w[k] = b15;
                                        k = k + 1;
                                        x[i1 + k * n] = -u;
                                        x[i2 + k * n] = +u;
                                        x[i3 + k * n] = +u;
                                        x[i4 + k * n] = +u;
                                        x[i5 + k * n] = -u;
                                        w[k] = b15;
                                        k = k + 1;
                                        x[i1 + k * n] = -u;
                                        x[i2 + k * n] = +u;
                                        x[i3 + k * n] = +u;
                                        x[i4 + k * n] = +u;
                                        x[i5 + k * n] = +u;
                                        w[k] = b15;
                                        k = k + 1;
                                        x[i1 + k * n] = +u;
                                        x[i2 + k * n] = -u;
                                        x[i3 + k * n] = -u;
                                        x[i4 + k * n] = -u;
                                        x[i5 + k * n] = -u;
                                        w[k] = b15;
                                        k = k + 1;
                                        x[i1 + k * n] = +u;
                                        x[i2 + k * n] = -u;
                                        x[i3 + k * n] = -u;
                                        x[i4 + k * n] = -u;
                                        x[i5 + k * n] = +u;
                                        w[k] = b15;
                                        k = k + 1;
                                        x[i1 + k * n] = +u;
                                        x[i2 + k * n] = -u;
                                        x[i3 + k * n] = -u;
                                        x[i4 + k * n] = +u;
                                        x[i5 + k * n] = -u;
                                        w[k] = b15;
                                        k = k + 1;
                                        x[i1 + k * n] = +u;
                                        x[i2 + k * n] = -u;
                                        x[i3 + k * n] = -u;
                                        x[i4 + k * n] = +u;
                                        x[i5 + k * n] = +u;
                                        w[k] = b15;
                                        k = k + 1;
                                        x[i1 + k * n] = +u;
                                        x[i2 + k * n] = -u;
                                        x[i3 + k * n] = +u;
                                        x[i4 + k * n] = -u;
                                        x[i5 + k * n] = -u;
                                        w[k] = b15;
                                        k = k + 1;
                                        x[i1 + k * n] = +u;
                                        x[i2 + k * n] = -u;
                                        x[i3 + k * n] = +u;
                                        x[i4 + k * n] = -u;
                                        x[i5 + k * n] = +u;
                                        w[k] = b15;
                                        k = k + 1;
                                        x[i1 + k * n] = +u;
                                        x[i2 + k * n] = -u;
                                        x[i3 + k * n] = +u;
                                        x[i4 + k * n] = +u;
                                        x[i5 + k * n] = -u;
                                        w[k] = b15;
                                        k = k + 1;
                                        x[i1 + k * n] = +u;
                                        x[i2 + k * n] = -u;
                                        x[i3 + k * n] = +u;
                                        x[i4 + k * n] = +u;
                                        x[i5 + k * n] = +u;
                                        w[k] = b15;
                                        k = k + 1;
                                        x[i1 + k * n] = +u;
                                        x[i2 + k * n] = +u;
                                        x[i3 + k * n] = -u;
                                        x[i4 + k * n] = -u;
                                        x[i5 + k * n] = -u;
                                        w[k] = b15;
                                        k = k + 1;
                                        x[i1 + k * n] = +u;
                                        x[i2 + k * n] = +u;
                                        x[i3 + k * n] = -u;
                                        x[i4 + k * n] = -u;
                                        x[i5 + k * n] = +u;
                                        w[k] = b15;
                                        k = k + 1;
                                        x[i1 + k * n] = +u;
                                        x[i2 + k * n] = +u;
                                        x[i3 + k * n] = -u;
                                        x[i4 + k * n] = +u;
                                        x[i5 + k * n] = -u;
                                        w[k] = b15;
                                        k = k + 1;
                                        x[i1 + k * n] = +u;
                                        x[i2 + k * n] = +u;
                                        x[i3 + k * n] = -u;
                                        x[i4 + k * n] = +u;
                                        x[i5 + k * n] = +u;
                                        w[k] = b15;
                                        k = k + 1;
                                        x[i1 + k * n] = +u;
                                        x[i2 + k * n] = +u;
                                        x[i3 + k * n] = +u;
                                        x[i4 + k * n] = -u;
                                        x[i5 + k * n] = -u;
                                        w[k] = b15;
                                        k = k + 1;
                                        x[i1 + k * n] = +u;
                                        x[i2 + k * n] = +u;
                                        x[i3 + k * n] = +u;
                                        x[i4 + k * n] = -u;
                                        x[i5 + k * n] = +u;
                                        w[k] = b15;
                                        k = k + 1;
                                        x[i1 + k * n] = +u;
                                        x[i2 + k * n] = +u;
                                        x[i3 + k * n] = +u;
                                        x[i4 + k * n] = +u;
                                        x[i5 + k * n] = -u;
                                        w[k] = b15;
                                        k = k + 1;
                                        x[i1 + k * n] = +u;
                                        x[i2 + k * n] = +u;
                                        x[i3 + k * n] = +u;
                                        x[i4 + k * n] = +u;
                                        x[i5 + k * n] = +u;
                                        w[k] = b15;
                                   }
                              }
                         }
                    }
               }

               return;
          }

          public static int en_r2_11_1_size(int n)

               //****************************************************************************80
               //
               //  Purpose:
               //
               //    EN_R2_11_1_SIZE sizes the Stroud rule 11.1 for region EN_R2.
               //
               //  Discussion:
               //
               //    The rule has order 
               //    O = ( 4 * N^5 - 20 * N^4 + 140 * N^3 - 130 * N^2 + 96 * N + 15 ) / 15.
               //
               //  Licensing:
               //
               //    This code is distributed under the GNU LGPL license.
               //
               //  Modified:
               //
               //    23 January 2010
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
               //    Input, int N, the spatial dimension.
               //
               //    Output, int EN_R2_11_1_SIZE, the order.
               //
          {
               int o;

               o = (4 * (int)Math.Pow(n, 5)
                    - 20 * (int)Math.Pow(n, 4)
                    + 140 * (int)Math.Pow(n, 3)
                    - 130 * (int)Math.Pow(n, 2)
                    + 96 * n
                    + 15) / 15;

               return o;
          }

     }
}