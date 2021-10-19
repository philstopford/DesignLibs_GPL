using System;
using Burkardt.Types;

namespace Burkardt
{
       public static class AlpertRule
       {
              public static int a_log(int rule)

                     //****************************************************************************80
                     //
                     //  Purpose:
                     //
                     //    A_LOG returns the value of A for an Alpert rule for log singular functions.
                     //
                     //  Licensing:
                     //
                     //    This code is distributed under the GNU LGPL license.
                     //
                     //  Modified:
                     //
                     //    05 December 2015
                     //
                     //  Author:
                     //
                     //    John Burkardt
                     //
                     //  Parameters:
                     //
                     //    Input, int RULE, the index of the rule, between 1 and 10.
                     //
                     //    Output, int A_LOG, the value of A.
                     //
              {
                     int a;
                     int[] a_vec =
                     {
                            1, 2, 2, 3, 3,
                            5, 6, 7, 9, 10
                     };

                     if (rule < 1 || 10 < rule)
                     {
                            Console.WriteLine("");
                            Console.WriteLine("A_LOG - Fatal error!");
                            Console.WriteLine("  Input value of RULE is not between 1 and 10.");
                            return (-1);
                     }

                     a = a_vec[rule - 1];

                     return a;
              }

              public static int a_power(int rule)

                     //****************************************************************************80
                     //
                     //  Purpose:
                     //
                     //    A_POWER returns A for an Alpert rule for power singular functions.
                     //
                     //  Licensing:
                     //
                     //    This code is distributed under the GNU LGPL license.
                     //
                     //  Modified:
                     //
                     //    05 December 2015
                     //
                     //  Author:
                     //
                     //    John Burkardt
                     //
                     //  Parameters:
                     //
                     //    Input, int RULE, the index of the rule, between 1 and 12.
                     //
                     //    Output, int A_POWER, the value of A.
                     //
              {
                     int a;
                     int[] a_vec =
                     {
                            1, 2, 2, 2, 2,
                            3, 4, 5, 6, 8,
                            9, 10
                     };

                     if (rule < 1 || 12 < rule)
                     {
                            Console.WriteLine("");
                            Console.WriteLine("A_POWER - Fatal error!");
                            Console.WriteLine("  Input value of RULE is not between 1 and 12.");
                            return (-1);
                     }

                     a = a_vec[rule - 1];

                     return a;
              }

              public static int a_regular(int rule)

                     //****************************************************************************80
                     //
                     //  Purpose:
                     //
                     //    A_REGULAR returns the value of A for an Alpert rule for regular functions.
                     //
                     //  Licensing:
                     //
                     //    This code is distributed under the GNU LGPL license.
                     //
                     //  Modified:
                     //
                     //    05 December 2015
                     //
                     //  Author:
                     //
                     //    John Burkardt
                     //
                     //  Parameters:
                     //
                     //    Input, int RULE, the index of the rule, between 1 and 12.
                     //
                     //    Output, int A_REGULAR, the value of A.
                     //
              {
                     int a;
                     int[] a_vec =
                     {
                            1, 2, 2, 3, 3,
                            4, 5, 7, 9, 10,
                            12, 14
                     };

                     if (rule < 1 || 12 < rule)
                     {
                            Console.WriteLine("");
                            Console.WriteLine("A_REGULAR - Fatal error!");
                            Console.WriteLine("  Input value of RULE is not between 1 and 12.");
                            return (-1);
                     }

                     a = a_vec[rule - 1];

                     return a;
              }

              public static double integral_log()

                     //****************************************************************************80
                     //
                     //  Purpose:
                     //
                     //    INTEGRAL_LOG evaluates the test integral with logarithmic singularity.
                     //
                     //  Licensing:
                     //
                     //    This code is distributed under the GNU LGPL license.
                     //
                     //  Modified:
                     //
                     //    05 December 2015
                     //
                     //  Author:
                     //
                     //    John Burkardt
                     //
                     //  Parameters:
                     //
                     //    Output, double INTEGRAL_LOG, the integral of the test integrand 
                     //    from 0 to 1.
                     //
              {
                     double value;

                     value = -0.012771107587415899716E+00;

                     return value;
              }

              public static double integral_power()

                     //****************************************************************************80
                     //
                     //  Purpose:
                     //
                     //    INTEGRAL_POWER evaluates the test integral with power singularity.
                     //
                     //  Licensing:
                     //
                     //    This code is distributed under the GNU LGPL license.
                     //
                     //  Modified:
                     //
                     //    05 December 2015
                     //
                     //  Author:
                     //
                     //    John Burkardt
                     //
                     //  Parameters:
                     //
                     //    Output, double INTEGRAL POWER, the integral of the test integrand 
                     //    from 0 to 1.
                     //
              {
                     double value;

                     value = 0.079321002746971411182E+00;

                     return value;
              }

              public static double integral_regular()

                     //****************************************************************************80
                     //
                     //  Purpose:
                     //
                     //    INTEGRAL_REGULAR evaluates the regular test integral.
                     //
                     //  Licensing:
                     //
                     //    This code is distributed under the GNU LGPL license.
                     //
                     //  Modified:
                     //
                     //    05 December 2015
                     //
                     //  Author:
                     //
                     //    John Burkardt
                     //
                     //  Parameters:
                     //
                     //    Output, double INTEGRAL_REGULAR, the integral of the test integrand 
                     //    from 0 to 1.
                     //
              {
                     double value;

                     value = (-Math.Sin(0.3E+00) + Math.Sin(200.0E+00)
                                                 + Math.Sin(200.3E+00)) / 200.0E+00;

                     return value;
              }

              public static double[] integrand_log(int n, double[] x)

                     //****************************************************************************80
                     //
                     //  Purpose:
                     //
                     //    INTEGRAND_LOG evaluates the test integrand with logarithmic singularity.
                     //
                     //  Licensing:
                     //
                     //    This code is distributed under the GNU LGPL license.
                     //
                     //  Modified:
                     //
                     //    05 December 2015
                     //
                     //  Author:
                     //
                     //    John Burkardt
                     //
                     //  Parameters:
                     //
                     //    Input, int N, the number of evaluation points.
                     //
                     //    Input, double X[N], the evaluation points.
                     //
                     //    Output, double INTEGRAND_LOG[N], the integrand at the evaluation points.
                     //
              {
                     double[] f;
                     int i;

                     f = new double[n];

                     for (i = 0; i < n; i++)
                     {
                            f[i] = Math.Cos(200.0E+00 * x[i]) * Math.Log(x[i])
                                   + Math.Cos(200.0E+00 * x[i] + 0.3E+00);
                     }

                     return f;
              }

              public static double[] integrand_power(int n, double[] x)

                     //****************************************************************************80
                     //
                     //  Purpose:
                     //
                     //    INTEGRAND_POWER evaluates the test integrand with power singularity.
                     //
                     //  Licensing:
                     //
                     //    This code is distributed under the GNU LGPL license.
                     //
                     //  Modified:
                     //
                     //    05 December 2015
                     //
                     //  Author:
                     //
                     //    John Burkardt
                     //
                     //  Parameters:
                     //
                     //    Input, int N, the number of evaluation points.
                     //
                     //    Input, double X[N], the evaluation points.
                     //
                     //    Output, double INTEGRAND_POWER[N], the integrand at the evaluation points.
                     //
              {
                     double[] f;
                     int i;

                     f = new double[n];

                     for (i = 0; i < n; i++)
                     {
                            f[i] = Math.Cos(200.0E+00 * x[i]) * Math.Pow(x[i], -0.5E+00)
                                   + Math.Cos(200.0E+00 * x[i] + 0.3E+00);
                     }

                     return f;
              }

              public static double[] integrand_regular(int n, double[] x)

                     //****************************************************************************80
                     //
                     //  Purpose:
                     //
                     //    INTEGRAND_REGULAR evaluates the regular test integrand.
                     //
                     //  Licensing:
                     //
                     //    This code is distributed under the GNU LGPL license.
                     //
                     //  Modified:
                     //
                     //    05 December 2015
                     //
                     //  Author:
                     //
                     //    John Burkardt
                     //
                     //  Parameters:
                     //
                     //    Input, int N, the number of evaluation points.
                     //
                     //    Input, double X[N], the evaluation points.
                     //
                     //    Output, double INTEGRAND_REGULAR[N], the integrand at the evaluation points.
                     //
              {
                     double[] f;
                     int i;

                     f = new double[n];

                     for (i = 0; i < n; i++)
                     {
                            f[i] = Math.Cos(200.0E+00 * x[i]) + Math.Cos(200.0E+00 * x[i] + 0.3E+00);
                     }

                     return f;
              }

              public static int j_log(int rule)

                     //****************************************************************************80
                     //
                     //  Purpose:
                     //
                     //    J_LOG returns the value of J for an Alpert rule for log singular functions.
                     //
                     //  Licensing:
                     //
                     //    This code is distributed under the GNU LGPL license.
                     //
                     //  Modified:
                     //
                     //    05 December 2015
                     //
                     //  Author:
                     //
                     //    John Burkardt
                     //
                     //  Parameters:
                     //
                     //    Input, int RULE, the index of the rule, between 1 and 10.
                     //
                     //    Output, int J_LOG, the value of J.
                     //
              {
                     int j;
                     int[] j_vec =
                     {
                            1, 2, 3, 4, 5,
                            7, 10, 11, 14, 15
                     };

                     if (rule < 1 || 10 < rule)
                     {
                            Console.WriteLine("");
                            Console.WriteLine("J_LOG - Fatal error!");
                            Console.WriteLine("  Input value of RULE is not between 1 and 10.");
                            return (-1);
                     }

                     j = j_vec[rule - 1];

                     return j;
              }

              public static int j_power(int rule)

                     //****************************************************************************80
                     //
                     //  Purpose:
                     //
                     //    J_POWER returns J for an Alpert rule for power singular functions.
                     //
                     //  Licensing:
                     //
                     //    This code is distributed under the GNU LGPL license.
                     //
                     //  Modified:
                     //
                     //    05 December 2015
                     //
                     //  Author:
                     //
                     //    John Burkardt
                     //
                     //  Parameters:
                     //
                     //    Input, int RULE, the index of the rule, between 1 and 12.
                     //
                     //    Output, int J_POWER, the value of J.
                     //
              {
                     int j;
                     int[] j_vec =
                     {
                            1, 2, 2, 3, 3,
                            4, 6, 8, 10, 12,
                            14, 16
                     };

                     if (rule < 1 || 12 < rule)
                     {
                            Console.WriteLine("");
                            Console.WriteLine("J_POWER - Fatal error!");
                            Console.WriteLine("  Input value of RULE is not between 1 and 12.");
                            return (-1);
                     }

                     j = j_vec[rule - 1];

                     return j;
              }

              public static int j_regular(int rule)

                     //****************************************************************************80
                     //
                     //  Purpose:
                     //
                     //    J_REGULAR returns the value of J for an Alpert rule for regular functions.
                     //
                     //  Licensing:
                     //
                     //    This code is distributed under the GNU LGPL license.
                     //
                     //  Modified:
                     //
                     //    05 December 2015
                     //
                     //  Author:
                     //
                     //    John Burkardt
                     //
                     //  Parameters:
                     //
                     //    Input, int RULE, the index of the rule, between 1 and 12.
                     //
                     //    Output, int J, the value of J.
                     //
              {
                     int j;
                     int[] j_vec =
                     {
                            1, 2, 2, 3, 3,
                            4, 6, 8, 10, 12,
                            14, 16
                     };

                     if (rule < 1 || 12 < rule)
                     {
                            Console.WriteLine("");
                            Console.WriteLine("J_REGULAR - Fatal error!");
                            Console.WriteLine("  Input value of RULE is not between 1 and 12.");
                            return (-1);
                     }

                     j = j_vec[rule - 1];

                     return j;
              }

              public static int num_log()

                     //****************************************************************************80
                     //
                     //  Purpose:
                     //
                     //    NUM_LOG returns the number of Alpert rules for log singular functions.
                     //
                     //  Licensing:
                     //
                     //    This code is distributed under the GNU LGPL license.
                     //
                     //  Modified:
                     //
                     //    05 December 2015
                     //
                     //  Author:
                     //
                     //    John Burkardt
                     //
                     //  Parameters:
                     //
                     //    Output, int NUM_LOG, the number of rules.
                     //
              {
                     int num;

                     num = 10;

                     return num;
              }

              public static int num_power()

                     //****************************************************************************80
                     //
                     //  Purpose:
                     //
                     //    NUM_POWER returns the number of Alpert rules for power singular functions.
                     //
                     //  Licensing:
                     //
                     //    This code is distributed under the GNU LGPL license.
                     //
                     //  Modified:
                     //
                     //    05 December 2015
                     //
                     //  Author:
                     //
                     //    John Burkardt
                     //
                     //  Parameters:
                     //
                     //    Output, int NUM_POWER, the number of rules.
                     //
              {
                     int num;

                     num = 12;

                     return num;
              }

              public static int num_regular()

                     //****************************************************************************80
                     //
                     //  Purpose:
                     //
                     //    NUM_REGULAR returns the number of Alpert rules for regular functions.
                     //
                     //  Licensing:
                     //
                     //    This code is distributed under the GNU LGPL license.
                     //
                     //  Modified:
                     //
                     //    05 December 2015
                     //
                     //  Author:
                     //
                     //    John Burkardt
                     //
                     //  Parameters:
                     //
                     //    Output, int NUM_REGULAR, the number of rules.
                     //
              {
                     int num;

                     num = 12;

                     return num;
              }

              public static int order_log(int rule)

                     //****************************************************************************80
                     //
                     //  Purpose:
                     //
                     //    ORDER_LOG returns the order of an Alpert rule for log singular functions.
                     //
                     //  Licensing:
                     //
                     //    This code is distributed under the GNU LGPL license.
                     //
                     //  Modified:
                     //
                     //    05 December 2015
                     //
                     //  Author:
                     //
                     //    John Burkardt
                     //
                     //  Parameters:
                     //
                     //    Input, int RULE, the index of the rule, between 1 and 10.
                     //
                     //    Output, int ORDER_LOG, the order of the rule.
                     //
              {
                     int order;
                     int[] order_vec =
                     {
                            2, 3, 4, 5, 6,
                            8, 10, 12, 14, 16
                     };

                     if (rule < 1 || 10 < rule)
                     {
                            Console.WriteLine("");
                            Console.WriteLine("ORDER_LOG - Fatal error!");
                            Console.WriteLine("  Input value of RULE is not between 1 and 10.");
                            return (-1);
                     }

                     order = order_vec[rule - 1];

                     return order;
              }

              public static double order_power(int rule)

                     //****************************************************************************80
                     //
                     //  Purpose:
                     //
                     //    ORDER_POWER returns the order of an Alpert rule for power singular functions.
                     //
                     //  Licensing:
                     //
                     //    This code is distributed under the GNU LGPL license.
                     //
                     //  Modified:
                     //
                     //    05 December 2015
                     //
                     //  Author:
                     //
                     //    John Burkardt
                     //
                     //  Parameters:
                     //
                     //    Input, int RULE, the index of the rule, between 1 and 10.
                     //
                     //    Output, double ORDER_POWER, the order of the rule.
                     //
              {
                     double order;
                     double[] order_vec =
                     {
                            1.5, 2.0, 2.5, 3.0, 3.5,
                            4.0, 6.0, 8.0, 10.0, 12.0,
                            14.0, 16.0
                     };

                     if (rule < 1 || 12 < rule)
                     {
                            Console.WriteLine("");
                            Console.WriteLine("ORDER_POWER - Fatal error!");
                            Console.WriteLine("  Input value of RULE is not between 1 and 12.");
                            return (-1);
                     }

                     order = order_vec[rule - 1];

                     return order;
              }

              public static int order_regular(int rule)

                     //****************************************************************************80
                     //
                     //  Purpose:
                     //
                     //    ORDER_REGULAR returns the order of an Alpert rule for regular functions.
                     //
                     //  Licensing:
                     //
                     //    This code is distributed under the GNU LGPL license.
                     //
                     //  Modified:
                     //
                     //    05 December 2015
                     //
                     //  Author:
                     //
                     //    John Burkardt
                     //
                     //  Parameters:
                     //
                     //    Input, int RULE, the index of the rule, between 1 and 12.
                     //
                     //    Output, int ORDER_REGULAR, the order of the rule.
                     //
              {
                     int order;
                     int[] order_vec =
                     {
                            3, 4, 5, 6, 7,
                            8, 12, 16, 20, 24,
                            28, 32
                     };

                     if (rule < 1 || 12 < rule)
                     {
                            Console.WriteLine("");
                            Console.WriteLine("ORDER_REGULAR - Fatal error!");
                            Console.WriteLine("  Input value of RULE is not between 1 and 12.");
                            return (-1);
                     }

                     order = order_vec[rule - 1];

                     return order;
              }

              public static void rule_log(int rule, int j, ref double[] x, ref double[] w)

                     //****************************************************************************80
                     //
                     //  Purpose:
                     //
                     //    RULE_LOG returns an Alpert rule for log singular functions.
                     //
                     //  Licensing:
                     //
                     //    This code is distributed under the GNU LGPL license.
                     //
                     //  Modified:
                     //
                     //    05 December 2015
                     //
                     //  Author:
                     //
                     //    John Burkardt
                     //
                     //  Reference:
                     //
                     //    Bradley Alpert,
                     //    Hybrid Gauss-Trapezoidal Quadrature Rules,
                     //    SIAM Journal on Scientific Computing,
                     //    Volume 20, Number 5, pages 1551-1584, 1999.
                     //
                     //  Parameters:
                     //
                     //    Input, int RULE, the index of the rule, between 1 and 10.
                     //
                     //    Input, int J, the number of points in the rule.
                     //
                     //    Output, double X[J], W[J], the points and weights for the rule.
                     //
              {
                     double[] x01 =
                     {
                            1.591549430918953E-01
                     };
                     double[] w01 =
                     {
                            5.0E-01
                     };
                     double[] x02 =
                     {
                            1.150395811972836E-01,
                            9.365464527949632E-01
                     };
                     double[] w02 =
                     {
                            3.913373788753340E-01,
                            1.108662621124666E+00
                     };
                     double[] x03 =
                     {
                            2.379647284118974E-02,
                            2.935370741501914E-01,
                            1.023715124251890E+00
                     };
                     double[] w03 =
                     {
                            8.795942675593887E-02,
                            4.989017152913699E-01,
                            9.131388579526912E-01
                     };
                     double[] x04 =
                     {
                            2.339013027203800E-02,
                            2.854764931311984E-01,
                            1.005403327220700E+00,
                            1.994970303994294E+00
                     };
                     double[] w04 =
                     {
                            8.609736556158105E-02,
                            4.847019685417959E-01,
                            9.152988869123725E-01,
                            1.013901778984250E+00
                     };
                     double[] x05 =
                     {
                            4.004884194926570E-03,
                            7.745655373336686E-02,
                            3.972849993523248E-01,
                            1.075673352915104E+00,
                            2.003796927111872E+00
                     };
                     double[] w05 =
                     {
                            1.671879691147102E-02,
                            1.636958371447360E-01,
                            4.981856569770637E-01,
                            8.372266245578912E-01,
                            9.841730844088381E-01
                     };
                     double[] x06 =
                     {
                            6.531815708567918E-03,
                            9.086744584657729E-02,
                            3.967966533375878E-01,
                            1.027856640525646E+00,
                            1.945288592909266E+00,
                            2.980147933889640E+00,
                            3.998861349951123E+00
                     };
                     double[] w06 =
                     {
                            2.462194198995203E-02,
                            1.701315866854178E-01,
                            4.609256358650077E-01,
                            7.947291148621894E-01,
                            1.008710414337933E+00,
                            1.036093649726216E+00,
                            1.004787656533285E+00
                     };
                     double[] x07 =
                     {
                            1.175089381227308E-03,
                            1.877034129831289E-02,
                            9.686468391426860E-02,
                            3.004818668002884E-01,
                            6.901331557173356E-01,
                            1.293695738083659E+00,
                            2.090187729798780E+00,
                            3.016719313149212E+00,
                            4.001369747872486E+00,
                            5.000025661793423E+00
                     };
                     double[] w07 =
                     {
                            4.560746882084207E-03,
                            3.810606322384757E-02,
                            1.293864997289512E-01,
                            2.884360381408835E-01,
                            4.958111914344961E-01,
                            7.077154600594529E-01,
                            8.741924365285083E-01,
                            9.661361986515218E-01,
                            9.957887866078700E-01,
                            9.998665787423845E-01
                     };
                     double[] x08 =
                     {
                            1.674223682668368E-03,
                            2.441110095009738E-02,
                            1.153851297429517E-01,
                            3.345898490480388E-01,
                            7.329740531807683E-01,
                            1.332305048525433E+00,
                            2.114358752325948E+00,
                            3.026084549655318E+00,
                            4.003166301292590E+00,
                            5.000141170055870E+00,
                            6.000001002441859E+00
                     };
                     double[] w08 =
                     {
                            6.364190780720557E-03,
                            4.723964143287529E-02,
                            1.450891158385963E-01,
                            3.021659470785897E-01,
                            4.984270739715340E-01,
                            6.971213795176096E-01,
                            8.577295622757315E-01,
                            9.544136554351155E-01,
                            9.919938052776484E-01,
                            9.994621875822987E-01,
                            9.999934408092805E-01
                     };
                     double[] x09 =
                     {
                            9.305182368545380E-04,
                            1.373832458434617E-02,
                            6.630752760779359E-02,
                            1.979971397622003E-01,
                            4.504313503816532E-01,
                            8.571888631101634E-01,
                            1.434505229617112E+00,
                            2.175177834137754E+00,
                            3.047955068386372E+00,
                            4.004974906813428E+00,
                            4.998525901820967E+00,
                            5.999523015116678E+00,
                            6.999963617883990E+00,
                            7.999999488130134E+00
                     };
                     double[] w09 =
                     {
                            3.545060644780164E-03,
                            2.681514031576498E-02,
                            8.504092035093420E-02,
                            1.854526216643691E-01,
                            3.251724374883192E-01,
                            4.911553747260108E-01,
                            6.622933417369036E-01,
                            8.137254578840510E-01,
                            9.235595514944174E-01,
                            9.821609923744658E-01,
                            1.000047394596121E+00,
                            1.000909336693954E+00,
                            1.000119534283784E+00,
                            1.000002835746089E+00
                     };
                     double[] x10 =
                     {
                            8.371529832014113E-04,
                            1.239382725542637E-02,
                            6.009290785739468E-02,
                            1.805991249601928E-01,
                            4.142832599028031E-01,
                            7.964747731112430E-01,
                            1.348993882467059E+00,
                            2.073471660264395E+00,
                            2.947904939031494E+00,
                            3.928129252248612E+00,
                            4.957203086563112E+00,
                            5.986360113977494E+00,
                            6.997957704791519E+00,
                            7.999888757524622E+00,
                            8.999998754306120E+00
                     };
                     double[] w10 =
                     {
                            3.190919086626234E-03,
                            2.423621380426338E-02,
                            7.740135521653088E-02,
                            1.704889420286369E-01,
                            3.029123478511309E-01,
                            4.652220834914617E-01,
                            6.401489637096768E-01,
                            8.051212946181061E-01,
                            9.362411945698647E-01,
                            1.014359775369075E+00,
                            1.035167721053657E+00,
                            1.020308624984610E+00,
                            1.004798397441514E+00,
                            1.000395017352309E+00,
                            1.000007149422537E+00
                     };

                     if (rule < 1 || 10 < rule)
                     {
                            Console.WriteLine("");
                            Console.WriteLine("RULE_LOG - Fatal error!");
                            Console.WriteLine("  Input value of RULE is not between 1 and 10.");
                            return;
                     }

                     if (rule == 1)
                     {
                            typeMethods.r8vec_copy(j, x01, ref x);
                            typeMethods.r8vec_copy(j, w01, ref w);
                     }
                     else if (rule == 2)
                     {
                            typeMethods.r8vec_copy(j, x02, ref x);
                            typeMethods.r8vec_copy(j, w02, ref w);
                     }
                     else if (rule == 3)
                     {
                            typeMethods.r8vec_copy(j, x03, ref x);
                            typeMethods.r8vec_copy(j, w03, ref w);
                     }
                     else if (rule == 4)
                     {
                            typeMethods.r8vec_copy(j, x04, ref x);
                            typeMethods.r8vec_copy(j, w04, ref w);
                     }
                     else if (rule == 5)
                     {
                            typeMethods.r8vec_copy(j, x05, ref x);
                            typeMethods.r8vec_copy(j, w05, ref w);
                     }
                     else if (rule == 6)
                     {
                            typeMethods.r8vec_copy(j, x06, ref x);
                            typeMethods.r8vec_copy(j, w06, ref w);
                     }
                     else if (rule == 7)
                     {
                            typeMethods.r8vec_copy(j, x07, ref x);
                            typeMethods.r8vec_copy(j, w07, ref w);
                     }
                     else if (rule == 8)
                     {
                            typeMethods.r8vec_copy(j, x08, ref x);
                            typeMethods.r8vec_copy(j, w08, ref w);
                     }
                     else if (rule == 9)
                     {
                            typeMethods.r8vec_copy(j, x09, ref x);
                            typeMethods.r8vec_copy(j, w09, ref w);
                     }
                     else if (rule == 10)
                     {
                            typeMethods.r8vec_copy(j, x10, ref x);
                            typeMethods.r8vec_copy(j, w10, ref w);
                     }

                     return;
              }

              public static void rule_power(int rule, int j, ref double[] x, ref double[] w)

                     //****************************************************************************80
                     //
                     //  Purpose:
                     //
                     //    RULE_POWER returns an Alpert rule for power singular functions.
                     //
                     //  Licensing:
                     //
                     //    This code is distributed under the GNU LGPL license.
                     //
                     //  Modified:
                     //
                     //    05 December 2015
                     //
                     //  Author:
                     //
                     //    John Burkardt
                     //
                     //  Reference:
                     //
                     //    Bradley Alpert,
                     //    Hybrid Gauss-Trapezoidal Quadrature Rules,
                     //    SIAM Journal on Scientific Computing,
                     //    Volume 20, Number 5, pages 1551-1584, 1999.
                     //
                     //  Parameters:
                     //
                     //    Input, int RULE, the index of the rule, between 1 and 12.
                     //
                     //    Input, int J, the number of points in the rule.
                     //
                     //    Output, double X[J], W[J], the points and weights for the rule.
                     //
              {
                     double[] x01 =
                     {
                            1.172258571393266E-01
                     };
                     double[] w01 =
                     {
                            5.000000000000000E-01
                     };
                     double[] x02 =
                     {
                            9.252112715421378E-02,
                            1.000000000000000E-00
                     };
                     double[] w02 =
                     {
                            4.198079625266162E-01,
                            1.080192037473384E+00
                     };
                     double[] x03 =
                     {
                            6.023873796408450E-02,
                            8.780704050676215E-01
                     };
                     double[] w03 =
                     {
                            2.858439990420468E-01,
                            1.214156000957953E+00
                     };
                     double[] x04 =
                     {
                            7.262978413470474E-03,
                            2.246325512521893E-01,
                            1.000000000000000E+00
                     };
                     double[] w04 =
                     {
                            3.907638767531813E-02,
                            4.873484056646474E-01,
                            9.735752066600344E-01
                     };
                     double[] x05 =
                     {
                            1.282368909458828E-02,
                            2.694286346792474E-01,
                            1.018414523786358E+00
                     };
                     double[] w05 =
                     {
                            6.363996663105925E-02,
                            5.077434578043636E-01,
                            9.286165755645772E-01
                     };
                     double[] x06 =
                     {
                            1.189242434021285E-02,
                            2.578220434738662E-01,
                            1.007750064585281E+00,
                            2.000000000000000E+00
                     };
                     double[] w06 =
                     {
                            5.927215035616424E-02,
                            4.955981740306228E-01,
                            9.427131290628058E-01,
                            1.002416546550407E+00
                     };
                     double[] x07 =
                     {
                            3.317925942699451E-03,
                            8.283019705296352E-02,
                            4.136094925726231E-01,
                            1.088744373688402E+00,
                            2.006482101852379E+00,
                            3.000000000000000E+00
                     };
                     double[] w07 =
                     {
                            1.681780929883469E-02,
                            1.755244404544475E-01,
                            5.039350503858001E-01,
                            8.266241339680867E-01,
                            9.773065848981277E-01,
                            9.997919809947032E-01
                     };
                     double[] x08 =
                     {
                            1.214130606523435E-03,
                            3.223952700027058E-02,
                            1.790935383649920E-01,
                            5.437663805244631E-01,
                            1.176116628396759E+00,
                            2.031848210716014E+00,
                            3.001961225690812E+00,
                            4.000000000000000E+00
                     };
                     double[] w08 =
                     {
                            6.199844884297793E-03,
                            7.106286791720044E-02,
                            2.408930104410471E-01,
                            4.975929263668960E-01,
                            7.592446540441226E-01,
                            9.322446399614420E-01,
                            9.928171438160095E-01,
                            9.999449125689846E-01
                     };
                     double[] x09 =
                     {
                            1.745862989163252E-04,
                            8.613670540457314E-03,
                            6.733385088703690E-02,
                            2.514488774733840E-01,
                            6.341845573737690E-01,
                            1.248404055083152E+00,
                            2.065688031953401E+00,
                            3.009199358662542E+00,
                            4.000416269690208E+00,
                            5.000000000000000E+00
                     };
                     double[] w09 =
                     {
                            1.016950985948944E-03,
                            2.294670686517670E-02,
                            1.076657968022888E-01,
                            2.734577662465576E-01,
                            4.978815591924992E-01,
                            7.256208919565360E-01,
                            8.952638690320078E-01,
                            9.778157465381624E-01,
                            9.983390781399277E-01,
                            9.999916342408948E-01
                     };
                     double[] x10 =
                     {
                            5.710218427206990E-04,
                            1.540424351115548E-02,
                            8.834248407196555E-02,
                            2.824462054509770E-01,
                            6.574869892305580E-01,
                            1.246541060977993E+00,
                            2.039218495130811E+00,
                            2.979333487049800E+00,
                            3.985772595393049E+00,
                            4.997240804311428E+00,
                            5.999868793951190E+00,
                            7.000000000000000E+00
                     };
                     double[] w10 =
                     {
                            2.921018926912141E-03,
                            3.431130611256885E-02,
                            1.224669495638615E-01,
                            2.761108242022520E-01,
                            4.797809643010337E-01,
                            6.966555677271379E-01,
                            8.790077941972658E-01,
                            9.868622449294327E-01,
                            1.015142389688201E+00,
                            1.006209712632210E+00,
                            1.000528829922287E+00,
                            1.000002397796838E+00
                     };
                     double[] x11 =
                     {
                            3.419821460249725E-04,
                            9.296593430187960E-03,
                            5.406214771755252E-02,
                            1.763945096508648E-01,
                            4.218486605653738E-01,
                            8.274022895884040E-01,
                            1.410287585637014E+00,
                            2.160997505238153E+00,
                            3.043504749358223E+00,
                            4.005692579069439E+00,
                            4.999732707905968E+00,
                            5.999875191971098E+00,
                            6.999994560568667E+00,
                            8.000000000000000E+00
                     };
                     double[] w11 =
                     {
                            1.750957243202047E-03,
                            2.080726584287380E-02,
                            7.586830616433430E-02,
                            1.766020526671851E-01,
                            3.206624362072232E-01,
                            4.934405290553812E-01,
                            6.707497030698472E-01,
                            8.244959025366557E-01,
                            9.314646742162802E-01,
                            9.845768443163154E-01,
                            9.992852769154770E-01,
                            1.000273112957723E+00,
                            1.000022857402321E+00,
                            1.000000081405180E+00
                     };
                     double[] x12 =
                     {
                            2.158438988280793E-04,
                            5.898432743709196E-03,
                            3.462795956896131E-02,
                            1.145586495070213E-01,
                            2.790344218856415E-01,
                            5.600113798653321E-01,
                            9.814091242883119E-01,
                            1.553594853974655E+00,
                            2.270179114036658E+00,
                            3.108234601715371E+00,
                            4.032930893996553E+00,
                            5.006803270228157E+00,
                            6.000815466735179E+00,
                            7.000045035079542E+00,
                            8.000000738923901E+00,
                            9.000000000000000E+00
                     };
                     double[] w12 =
                     {
                            1.105804873501181E-03,
                            1.324499944707956E-02,
                            4.899842307592144E-02,
                            1.165326192868815E-01,
                            2.178586693194957E-01,
                            3.481766016945031E-01,
                            4.964027915911545E-01,
                            6.469026189623831E-01,
                            7.823688971783889E-01,
                            8.877772445893361E-01,
                            9.551665077035583E-01,
                            9.876285579741800E-01,
                            9.979929183863017E-01,
                            9.998470620634641E-01,
                            9.999962891645340E-01,
                            9.999999946893169E-01
                     };

                     if (rule < 1 || 12 < rule)
                     {
                            Console.WriteLine("");
                            Console.WriteLine("RULE_POWER - Fatal error!");
                            Console.WriteLine("  Input value of RULE is not between 1 and 12.");
                            return;
                     }

                     if (rule == 1)
                     {
                            typeMethods.r8vec_copy(j, x01, ref x);
                            typeMethods.r8vec_copy(j, w01, ref w);
                     }
                     else if (rule == 2)
                     {
                            typeMethods.r8vec_copy(j, x02, ref x);
                            typeMethods.r8vec_copy(j, w02, ref w);
                     }
                     else if (rule == 3)
                     {
                            typeMethods.r8vec_copy(j, x03, ref x);
                            typeMethods.r8vec_copy(j, w03, ref w);
                     }
                     else if (rule == 4)
                     {
                            typeMethods.r8vec_copy(j, x04, ref x);
                            typeMethods.r8vec_copy(j, w04, ref w);
                     }
                     else if (rule == 5)
                     {
                            typeMethods.r8vec_copy(j, x05, ref x);
                            typeMethods.r8vec_copy(j, w05, ref w);
                     }
                     else if (rule == 6)
                     {
                            typeMethods.r8vec_copy(j, x06, ref x);
                            typeMethods.r8vec_copy(j, w06, ref w);
                     }
                     else if (rule == 7)
                     {
                            typeMethods.r8vec_copy(j, x07, ref x);
                            typeMethods.r8vec_copy(j, w07, ref w);
                     }
                     else if (rule == 8)
                     {
                            typeMethods.r8vec_copy(j, x08, ref x);
                            typeMethods.r8vec_copy(j, w08, ref w);
                     }
                     else if (rule == 9)
                     {
                            typeMethods.r8vec_copy(j, x09, ref x);
                            typeMethods.r8vec_copy(j, w09, ref w);
                     }
                     else if (rule == 10)
                     {
                            typeMethods.r8vec_copy(j, x10, ref x);
                            typeMethods.r8vec_copy(j, w10, ref w);
                     }
                     else if (rule == 11)
                     {
                            typeMethods.r8vec_copy(j, x11, ref x);
                            typeMethods.r8vec_copy(j, w11, ref w);
                     }
                     else if (rule == 12)
                     {
                            typeMethods.r8vec_copy(j, x12, ref x);
                            typeMethods.r8vec_copy(j, w12, ref w);
                     }

                     return;
              }

              public static void rule_regular(int rule, int j, ref double[] x, ref double[] w)

                     //****************************************************************************80
                     //
                     //  Purpose:
                     //
                     //    RULE_REGULAR returns an Alpert rule for regular functions.
                     //
                     //  Licensing:
                     //
                     //    This code is distributed under the GNU LGPL license.
                     //
                     //  Modified:
                     //
                     //    05 December 2015
                     //
                     //  Author:
                     //
                     //    John Burkardt
                     //
                     //  Reference:
                     //
                     //    Bradley Alpert,
                     //    Hybrid Gauss-Trapezoidal Quadrature Rules,
                     //    SIAM Journal on Scientific Computing,
                     //    Volume 20, Number 5, pages 1551-1584, 1999.
                     //
                     //  Parameters:
                     //
                     //    Input, int RULE, the index of the rule, between 1 and 12.
                     //
                     //    Input, int J, the number of points in the rule.
                     //
                     //    Output, double X[J], W[J], the points and weights for the rule.
                     //
              {
                     double[] x01 =
                     {
                            1.666666666666667E-01
                     };
                     double[] w01 =
                     {
                            5.000000000000000E-01
                     };
                     double[] x02 =
                     {
                            2.000000000000000E-01,
                            1.000000000000000E+00
                     };
                     double[] w02 =
                     {
                            5.208333333333333E-01,
                            9.791666666666667E-01
                     };
                     double[] x03 =
                     {
                            2.245784979812614E-01,
                            1.013719374359164E+00
                     };
                     double[] w03 =
                     {
                            5.540781643606372E-01,
                            9.459218356393628E-01
                     };
                     double[] x04 =
                     {
                            2.250991042610971E-01,
                            1.014269060987992E+00,
                            2.000000000000000E+00
                     };
                     double[] w04 =
                     {
                            5.549724327164180E-01,
                            9.451317411845473E-01,
                            9.998958260990347E-01
                     };
                     double[] x05 =
                     {
                            2.180540672543505E-01,
                            1.001181873031216E+00,
                            1.997580526418033E+00
                     };
                     double[] w05 =
                     {
                            5.408088967208193E-01,
                            9.516615045823566E-01,
                            1.007529598696824E+00
                     };
                     double[] x06 =
                     {
                            2.087647422032129E-01,
                            9.786087373714483E-01,
                            1.989541386579751E+00,
                            3.000000000000000E+00
                     };
                     double[] w06 =
                     {
                            5.207988277246498E-01,
                            9.535038018555888E-01,
                            1.024871626402471E+00,
                            1.000825744017291E+00
                     };
                     double[] x07 =
                     {
                            7.023955461621939E-02,
                            4.312297857227970E-01,
                            1.117752734518115E+00,
                            2.017343724572518E+00,
                            3.000837842847590E+00,
                            4.000000000000000E+00
                     };
                     double[] w07 =
                     {
                            1.922315977843698E-01,
                            5.348399530514687E-01,
                            8.170209442488760E-01,
                            9.592111521445966E-01,
                            9.967143408044999E-01,
                            9.999820119661890E-01
                     };
                     double[] x08 =
                     {
                            9.919337841451028E-02,
                            5.076592669645529E-01,
                            1.184972925827278E+00,
                            2.047493467134072E+00,
                            3.007168911869310E+00,
                            4.000474996776184E+00,
                            5.000007879022339E+00,
                            6.000000000000000E+00
                     };
                     double[] w08 =
                     {
                            2.528198928766921E-01,
                            5.550158230159486E-01,
                            7.852321453615224E-01,
                            9.245915673876714E-01,
                            9.839350200445296E-01,
                            9.984463448413151E-01,
                            9.999592378464547E-01,
                            9.999999686258662E-01
                     };
                     double[] x09 =
                     {
                            9.209200446233291E-02,
                            4.752021947758861E-01,
                            1.124687945844539E+00,
                            1.977387385642367E+00,
                            2.953848957822108E+00,
                            3.976136786048776E+00,
                            4.994354281979877E+00,
                            5.999469539335291E+00,
                            6.999986704874333E+00,
                            8.000000000000000E+00
                     };
                     double[] w09 =
                     {
                            2.351836144643984E-01,
                            5.248820509085946E-01,
                            7.634026409869887E-01,
                            9.284711336658351E-01,
                            1.010969886587741E+00,
                            1.024959725311073E+00,
                            1.010517534639652E+00,
                            1.001551595797932E+00,
                            1.000061681794188E+00,
                            1.000000135843597E+00
                     };
                     double[] x10 =
                     {
                            6.001064731474805E-02,
                            3.149685016229433E-01,
                            7.664508240518316E-01,
                            1.396685781342510E+00,
                            2.175195903206602E+00,
                            3.062320575880355E+00,
                            4.016440988792476E+00,
                            5.002872064275734E+00,
                            6.000285453310164E+00,
                            7.000012964962529E+00,
                            8.000000175554469E+00,
                            9.000000000000000E+00
                     };
                     double[] w10 =
                     {
                            1.538932104518340E-01,
                            3.551058128559424E-01,
                            5.449200036280007E-01,
                            7.104078497715549E-01,
                            8.398780940253654E-01,
                            9.272767950890611E-01,
                            9.750605697371132E-01,
                            9.942629650823470E-01,
                            9.992421778421898E-01,
                            9.999534370786161E-01,
                            9.999990854912925E-01,
                            9.999999989466828E-01
                     };
                     double[] x11 =
                     {
                            6.234360533194102E-02,
                            3.250286721702614E-01,
                            7.837350794282182E-01,
                            1.415673112616924E+00,
                            2.189894250061313E+00,
                            3.070053877483040E+00,
                            4.018613756218047E+00,
                            5.002705902035397E+00,
                            5.999929741810400E+00,
                            6.999904720846024E+00,
                            7.999986894843540E+00,
                            8.999999373380393E+00,
                            9.999999992002911E+00,
                            1.100000000000000E+01
                     };
                     double[] w11 =
                     {
                            1.595975279734157E-01,
                            3.637046028193864E-01,
                            5.498753177297441E-01,
                            7.087986792086956E-01,
                            8.335172275501195E-01,
                            9.204446510608518E-01,
                            9.710881776552090E-01,
                            9.933296578555239E-01,
                            9.994759087910050E-01,
                            1.000133030254421E+00,
                            1.000032915011460E+00,
                            1.000002261653775E+00,
                            1.000000042393520E+00,
                            1.000000000042872E+00
                     };
                     double[] x12 =
                     {
                            5.899550614325259E-02,
                            3.082757062227814E-01,
                            7.463707253079130E-01,
                            1.355993726494664E+00,
                            2.112943217346336E+00,
                            2.987241496545946E+00,
                            3.944798920961176E+00,
                            4.950269202842798E+00,
                            5.972123043117706E+00,
                            6.989783558137742E+00,
                            7.997673019512965E+00,
                            8.999694932747039E+00,
                            9.999979225211805E+00,
                            1.099999938266130E+01,
                            1.199999999462073E+01,
                            1.300000000000000E+01
                     };
                     double[] w12 =
                     {
                            1.511076023874179E-01,
                            3.459395921169090E-01,
                            5.273502805146873E-01,
                            6.878444094543021E-01,
                            8.210319140034114E-01,
                            9.218382875515803E-01,
                            9.873027487553060E-01,
                            1.018251913441155E+00,
                            1.021933430349293E+00,
                            1.012567983413513E+00,
                            1.004052289554521E+00,
                            1.000713413344501E+00,
                            1.000063618302950E+00,
                            1.000002486385216E+00,
                            1.000000030404477E+00,
                            1.000000000020760E+00
                     };

                     if (rule < 1 || 12 < rule)
                     {
                            Console.WriteLine("");
                            Console.WriteLine("RULE_REGULAR - Fatal error!");
                            Console.WriteLine("  Input value of RULE is not between 1 and 12.");
                            return;
                     }

                     if (rule == 1)
                     {
                            typeMethods.r8vec_copy(j, x01, ref x);
                            typeMethods.r8vec_copy(j, w01, ref w);
                     }
                     else if (rule == 2)
                     {
                            typeMethods.r8vec_copy(j, x02, ref x);
                            typeMethods.r8vec_copy(j, w02, ref w);
                     }
                     else if (rule == 3)
                     {
                            typeMethods.r8vec_copy(j, x03, ref x);
                            typeMethods.r8vec_copy(j, w03, ref w);
                     }
                     else if (rule == 4)
                     {
                            typeMethods.r8vec_copy(j, x04, ref x);
                            typeMethods.r8vec_copy(j, w04, ref w);
                     }
                     else if (rule == 5)
                     {
                            typeMethods.r8vec_copy(j, x05, ref x);
                            typeMethods.r8vec_copy(j, w05, ref w);
                     }
                     else if (rule == 6)
                     {
                            typeMethods.r8vec_copy(j, x06, ref x);
                            typeMethods.r8vec_copy(j, w06, ref w);
                     }
                     else if (rule == 7)
                     {
                            typeMethods.r8vec_copy(j, x07, ref x);
                            typeMethods.r8vec_copy(j, w07, ref w);
                     }
                     else if (rule == 8)
                     {
                            typeMethods.r8vec_copy(j, x08, ref x);
                            typeMethods.r8vec_copy(j, w08, ref w);
                     }
                     else if (rule == 9)
                     {
                            typeMethods.r8vec_copy(j, x09, ref x);
                            typeMethods.r8vec_copy(j, w09, ref w);
                     }
                     else if (rule == 10)
                     {
                            typeMethods.r8vec_copy(j, x10, ref x);
                            typeMethods.r8vec_copy(j, w10, ref w);
                     }
                     else if (rule == 11)
                     {
                            typeMethods.r8vec_copy(j, x11, ref x);
                            typeMethods.r8vec_copy(j, w11, ref w);
                     }
                     else if (rule == 12)
                     {
                            typeMethods.r8vec_copy(j, x12, ref x);
                            typeMethods.r8vec_copy(j, w12, ref w);
                     }

              }

       }
}