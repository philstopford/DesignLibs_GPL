using System;
using Burkardt.Quadrature;

namespace Burkardt.Stroud;

public static class Tetrahedron
{
       public static double tetra_07(int setting, Func<int, double, double, double, double> func, double[] x,
                     double[] y, double[] z)

              //****************************************************************************80
              //
              //  Purpose:
              //
              //    TETRA_07 approximates an integral inside a tetrahedron in 3D.
              //
              //  Integration region:
              //
              //    Points inside a tetrahedron whose four corners are given.
              //
              //  Discussion:
              //
              //    A 64 point 7-th degree conical product Gauss formula is used,
              //    Stroud number T3:7-1.
              //
              //  Licensing:
              //
              //    This code is distributed under the GNU LGPL license. 
              //
              //  Modified:
              //
              //    15 March 2008
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
              //    Arthur Stroud, Don Secrest,
              //    Gaussian Quadrature Formulas,
              //    Prentice Hall, 1966, pages 42-43,
              //    LC: QA299.4G3S7
              //
              //  Parameters:
              //
              //    Input, Func< double, double, double, double > func, the name of the 
              //    user supplied function to be integrated.
              //
              //    Input, double X[4], Y[4], Z[4], the coordinates of 
              //    the vertices.
              //
              //    Output, double TETRAQ_07, the approximate integral of the function.
              //
       {
              double a;
              double b;
              double c;
              double d;
              int i;
              int j;
              int k;
              int order = 4;
              double quad;
              double result;
              double t;
              double u;
              double v;
              double volume;
              double w;
              double[] weight1 = new double[4];
              double[] weight2 =
              {
                     0.1355069134, 0.2034645680, 0.1298475476, 0.0311809709
              };
              double[] weight3 =
              {
                     0.1108884156, 0.1434587898, 0.0686338872, 0.0103522407
              };
              double[] xtab1 = new double[4];
              double[] xtab2 =
              {
                     0.0571041961, 0.2768430136, 0.5835904324, 0.8602401357
              };
              double[] xtab3 =
              {
                     0.0485005495, 0.2386007376, 0.5170472951, 0.7958514179
              };
              double xval;
              double yval;
              double zval;
              //
              //  Get the Gauss-Legendre weights and abscissas for [-1,1].
              //
              LegendreQuadrature.legendre_set(order, ref xtab1, ref weight1);
              //
              //  Adjust the rule for the interval [0,1].
              //
              a = -1.0;
              b = +1.0;

              c = 0.0;
              d = 1.0;

              QuadratureRule.rule_adjust(a, b, c, d, order, ref xtab1, ref weight1);
              //
              //  Carry out the quadrature.
              //
              quad = 0.0;

              for (i = 0; i < order; i++)
              {
                     for (j = 0; j < order; j++)
                     {
                            for (k = 0; k < order; k++)
                            {
                                   //
                                   //  Compute the barycentric coordinates of the point in the unit triangle.
                                   //
                                   t = xtab3[k];
                                   u = xtab2[j] * (1.0 - xtab3[k]);
                                   v = xtab1[i] * (1.0 - xtab2[j]) * (1.0 - xtab3[k]);
                                   w = 1.0 - t - u - v;
                                   //
                                   //  Compute the corresponding point in the triangle.
                                   //
                                   xval = t * x[0] + u * x[1] + v * x[2] + w * x[3];
                                   yval = t * y[0] + u * y[1] + v * y[2] + w * y[3];
                                   zval = t * z[0] + u * z[1] + v * z[2] + w * z[3];

                                   quad += 6.0 * weight1[i] * weight2[j] * weight3[k]
                                           * func(setting, xval, yval, zval);
                            }
                     }
              }

              volume = tetra_volume(x, y, z);
              result = quad * volume;

              return result;
       }

       public static double tetra_sum(int setting, Func<int, double, double, double, double> func, double[] x,
                     double[] y, double[] z, int order, double[] xtab, double[] ytab,
                     double[] ztab, double[] weight)

              //****************************************************************************80
              //
              //  Purpose:
              //
              //    TETRA_SUM carries out a quadrature rule in a tetrahedron in 3D.
              //
              //  Integration region:
              //
              //    A tetrahedron whose vertices are specified.
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
              //  Parameters:
              //
              //    Input, Func< double, double, double, double > func, the name of the
              //    user supplied function which is to be integrated.
              //
              //    Input, double X[4], Y[4], Z[4], the vertices.
              //
              //    Input, int ORDER, the order of the rule.
              //
              //    Input, double XTAB[ORDER], YTAB[ORDER], ZTAB[ORDER], the
              //    abscissas.
              //
              //    Input, double WEIGHT[ORDER], the weights.
              //
              //    Output, double TETRA_SUM, the approximate integral of the function.
              //
       {
              int i;
              double quad;
              double result;
              double volume;
              double xval;
              double yval;
              double zval;

              quad = 0.0;
              for (i = 0; i < order; i++)
              {
                     xval = xtab[i] * x[0]
                            + ytab[i] * x[1]
                            + ztab[i] * x[2]
                            + (1.0 - xtab[i] - ytab[i] - ztab[i]) * x[3];

                     yval = xtab[i] * y[0]
                            + ytab[i] * y[1]
                            + ztab[i] * y[2]
                            + (1.0 - xtab[i] - ytab[i] - ztab[i]) * y[3];

                     zval = xtab[i] * z[0]
                            + ytab[i] * z[1]
                            + ztab[i] * z[2]
                            + (1.0 - xtab[i] - ytab[i] - ztab[i]) * z[3];

                     quad += weight[i] * func(setting, xval, yval, zval);
              }

              volume = tetra_volume(x, y, z);
              result = quad * volume;

              return result;
       }

       public static double tetra_tproduct(int setting, Func<int, double, double, double, double> func,
                     int order, double[] x, double[] y, double[] z)

              //****************************************************************************80
              //
              //  Purpose:
              //
              //    TETRA_TPRODUCT approximates an integral in a tetrahedron in 3D.
              //
              //  Discussion:
              //
              //    Integration is carried out over the points inside an arbitrary
              //    tetrahedron whose four vertices are given.
              //
              //    An ORDER**3 point (2*ORDER-1)-th degree triangular product
              //    Gauss-Legendre rule is used.
              //
              //    With ORDER = 8, this routine is equivalent to the routine TETR15
              //    in the reference, page 367.
              //
              //    Thanks to Joerg Behrens, jbehren@gwdg.de, for numerous suggestions
              //    and corrections.
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
              //    Approximate Calculation of Multiple Integrals,
              //    Prentice Hall, 1971,
              //    ISBN: 0130438936,
              //    LC: QA311.S85.
              //
              //  Parameters:
              //
              //    Input, Func< double, double, double, double > func, the name of the 
              //    user supplied function to be integrated.
              //
              //    Input, int ORDER, the order of the basic quadrature rules.
              //    ORDER should be between 1 and 9.
              //
              //    Input, double X[4], Y[4], Z[4], the vertices
              //    of the tetrahedron.
              //
              //    Output, double TETRA_TPRODUCT, the approximate integral of the function.
              //
       {
              double a;
              double b;
              double c;
              double d;
              int i;
              int j;
              int k;
              double quad;
              double result;
              double volume;
              double[] weight0;
              double[] weight1;
              double[] weight2;
              double[] xtab0;
              double[] xtab1;
              double[] xtab2;
              double xval;
              double yval;
              double zval;

              switch (order)
              {
                     case < 1:
                     case > 9:
                            Console.WriteLine("");
                            Console.WriteLine("TETRA_TPRODUCT - Fatal error!");
                            Console.WriteLine("  The quadrature rule orders must be between 1 and 9.");
                            Console.WriteLine("  The input value was ORDER = " + order + "");
                            return 1;
              }

              //
              //  Get the Gauss-Legendre ORDER point rules on [-1,1] for integrating
              //    F(X),
              //    X * F(X),
              //    X * X * F(X).
              //
              xtab0 = new double[order];
              xtab1 = new double[order];
              xtab2 = new double[order];
              weight0 = new double[order];
              weight1 = new double[order];
              weight2 = new double[order];

              LegendreQuadrature.legendre_set(order, ref xtab0, ref weight0);
              LegendreQuadrature.legendre_set_x1(order, ref xtab1, ref weight1);
              LegendreQuadrature.legendre_set_x2(order, ref xtab2, ref weight2);
              //
              //  Adjust the rules from [-1,1] to [0,1].
              //
              a = -1.0;
              b = +1.0;
              c = 0.0;
              d = 1.0;

              QuadratureRule.rule_adjust(a, b, c, d, order, ref xtab0, ref weight0);

              QuadratureRule.rule_adjust(a, b, c, d, order, ref xtab1, ref weight1);

              QuadratureRule.rule_adjust(a, b, c, d, order, ref xtab2, ref weight2);
              //
              //  For rules with a weight function that is not 1, the weight vectors
              //  require further adjustment.
              //
              for (i = 0; i < order; i++)
              {
                     weight1[i] /= 2.0;
              }

              for (i = 0; i < order; i++)
              {
                     weight2[i] /= 4.0;
              }

              //
              //  Carry out the quadrature.
              //
              quad = 0.0;

              for (k = 0; k < order; k++)
              {
                     for (j = 0; j < order; j++)
                     {
                            for (i = 0; i < order; i++)
                            {
                                   xval = x[0] + (((x[3] - x[2]) * xtab0[i]
                                                   + (x[2] - x[1])) * xtab1[j]
                                                  + (x[1] - x[0])) * xtab2[k];

                                   yval = y[0] + (((y[3] - y[2]) * xtab0[i]
                                                   + (y[2] - y[1])) * xtab1[j]
                                                  + (y[1] - y[0])) * xtab2[k];

                                   zval = z[0] + (((z[3] - z[2]) * xtab0[i]
                                                   + (z[2] - z[1])) * xtab1[j]
                                                  + (z[1] - z[0])) * xtab2[k];

                                   quad += 6.0 * weight0[i] * weight1[j] * weight2[k]
                                           * func(setting, xval, yval, zval);
                            }
                     }
              }

              //
              //  Compute the volume of the tetrahedron.
              //
              volume = tetra_volume(x, y, z);
              result = quad * volume;

              return result;
       }

       public static void tetra_unit_set(int rule, int order, ref double[] xtab, ref double[] ytab,
                     ref double[] ztab, ref double[] weight)

              //****************************************************************************80
              //
              //  Purpose:
              //
              //    TETRA_UNIT_SET sets quadrature weights and abscissas in the unit tetrahedron.
              //
              //  Integration region:
              //
              //      0 <= X,
              //    and
              //      0 <= Y,
              //    and
              //      0 <= Z, 
              //    and
              //      X + Y + Z <= 1.
              //
              //  Licensing:
              //
              //    This code is distributed under the GNU LGPL license. 
              //
              //  Modified:
              //
              //    19 March 2008
              //
              //  Author:
              //
              //    John Burkardt
              //
              //  Reference:
              //
              //    Hermann Engels,
              //    Numerical Quadrature and Cubature,
              //    Academic Press, 1980,
              //    ISBN: 012238850X,
              //    LC: QA299.3E5.
              //
              //    Patrick Keast,
              //    Moderate Degree Tetrahedral Quadrature Formulas,
              //    Computer Methods in Applied Mechanics and Engineering,
              //    Volume 55, Number 3, May 1986, pages 339-348.
              //
              //    Olgierd Zienkiewicz,
              //    The Finite Element Method,
              //    Sixth Edition,
              //    Butterworth-Heinemann, 2005,
              //    ISBN: 0750663200,
              //    LC: TA640.2.Z54
              //
              //  Parameters:
              //
              //    Input, int RULE, the index of the rule.
              //     1, order 1, precision 0, Newton Cotes formula #0, Zienkiewicz #1.
              //     2, order 4, precision 1, Newton Cotes formula #1.
              //     3, order 4, precision 2, Zienkiewicz #2.
              //     4, order 10, precision 2, Newton Cotes formula #2
              //     5, order 5, precision 3, Zienkiewicz #3.
              //     6, order 8, precision 3, Newton Cotes formula #3.
              //     7, order 35, precision 4, Newton Cotes formula #4.
              //     8, order 11, precision 4, a Keast rule.
              //
              //    Input, int ORDER, the order of the rule.
              //
              //    Output, double XTAB[ORDER], YTAB[ORDER], ZTAB[ORDER],
              //    the abscissas.
              //
              //    Output, double WEIGHT[ORDER], the weights.
              //
       {
              double a;
              double b;
              double c;
              double d;
              double e;
              double f;
              double g;
              double h;
              double z;
              switch (rule)
              {
                     //
                     //  Newton Cotes #0.
                     //
                     case 1:
                            xtab[0] = 0.25;
                            ytab[0] = 0.25;
                            ztab[0] = 0.25;
                            weight[0] = 1.0;
                            break;
                     //
                     //  Newton Cotes #1.
                     //
                     case 2:
                            a = 1.0;
                            b = 1.0 / 4.0;
                            z = 0.0;

                            xtab[0] = z;
                            xtab[1] = a;
                            xtab[2] = z;
                            xtab[3] = z;

                            ytab[0] = z;
                            ytab[1] = z;
                            ytab[2] = a;
                            ytab[3] = z;

                            ztab[0] = z;
                            ztab[1] = z;
                            ztab[2] = z;
                            ztab[3] = a;

                            weight[0] = b;
                            weight[1] = b;
                            weight[2] = b;
                            weight[3] = b;
                            break;
                     //
                     //  Zienkiewicz #2.
                     //
                     case 3:
                            a = 0.5854101966249685;
                            b = 0.1381966011250105;
                            c = 0.25;

                            xtab[0] = a;
                            xtab[1] = b;
                            xtab[2] = b;
                            xtab[3] = b;

                            ytab[0] = b;
                            ytab[1] = a;
                            ytab[2] = b;
                            ytab[3] = b;

                            ztab[0] = b;
                            ztab[1] = b;
                            ztab[2] = a;
                            ztab[3] = b;

                            weight[0] = c;
                            weight[1] = c;
                            weight[2] = c;
                            weight[3] = c;
                            break;
                     //
                     //  Newton Cotes #2.
                     //
                     case 4:
                            a = 1.0;
                            b = 0.5;
                            c = -1.0 / 20.0;
                            d = 4.0 / 20.0;
                            z = 0.0;

                            xtab[0] = z;
                            xtab[1] = a;
                            xtab[2] = z;
                            xtab[3] = z;
                            xtab[4] = b;
                            xtab[5] = z;
                            xtab[6] = z;
                            xtab[7] = b;
                            xtab[8] = b;
                            xtab[9] = z;

                            ytab[0] = z;
                            ytab[1] = z;
                            ytab[2] = a;
                            ytab[3] = z;
                            ytab[4] = z;
                            ytab[5] = b;
                            ytab[6] = z;
                            ytab[7] = b;
                            ytab[8] = z;
                            ytab[9] = b;

                            ztab[0] = z;
                            ztab[1] = z;
                            ztab[2] = z;
                            ztab[3] = a;
                            ztab[4] = z;
                            ztab[5] = z;
                            ztab[6] = b;
                            ztab[7] = z;
                            ztab[8] = b;
                            ztab[9] = b;

                            weight[0] = c;
                            weight[1] = c;
                            weight[2] = c;
                            weight[3] = c;
                            weight[4] = d;
                            weight[5] = d;
                            weight[6] = d;
                            weight[7] = d;
                            weight[8] = d;
                            weight[9] = d;
                            break;
                     //
                     //  Zienkiewicz #3.
                     //
                     case 5:
                            a = 1.0 / 6.0;
                            b = 0.25;
                            c = 0.5;
                            d = -0.8;
                            e = 0.45;

                            xtab[0] = b;
                            xtab[1] = c;
                            xtab[2] = a;
                            xtab[3] = a;
                            xtab[4] = a;

                            ytab[0] = b;
                            ytab[1] = a;
                            ytab[2] = c;
                            ytab[3] = a;
                            ytab[4] = a;

                            ztab[0] = b;
                            ztab[1] = a;
                            ztab[2] = a;
                            ztab[3] = c;
                            ztab[4] = a;

                            weight[0] = d;
                            weight[1] = e;
                            weight[2] = e;
                            weight[3] = e;
                            weight[4] = e;
                            break;
                     //
                     //  Newton Cotes #3.
                     //  (This is actually formally a 20 point rule, but with 12 zero coefficients.)
                     //
                     case 6:
                            a = 1.0;
                            b = 1.0 / 40.0;
                            c = 1.0 / 3.0;
                            d = 9.0 / 40.0;
                            z = 0.0;

                            xtab[0] = z;
                            xtab[1] = a;
                            xtab[2] = z;
                            xtab[3] = z;
                            xtab[4] = c;
                            xtab[5] = c;
                            xtab[6] = z;
                            xtab[7] = c;

                            ytab[0] = z;
                            ytab[1] = z;
                            ytab[2] = a;
                            ytab[3] = z;
                            ytab[4] = c;
                            ytab[5] = z;
                            ytab[6] = c;
                            ytab[7] = c;

                            ztab[0] = z;
                            ztab[1] = z;
                            ztab[2] = z;
                            ztab[3] = a;
                            ztab[4] = z;
                            ztab[5] = c;
                            ztab[6] = c;
                            ztab[7] = c;

                            weight[0] = b;
                            weight[1] = b;
                            weight[2] = b;
                            weight[3] = b;
                            weight[4] = d;
                            weight[5] = d;
                            weight[6] = d;
                            weight[7] = d;
                            break;
                     //
                     //  Newton Cotes #4.
                     //
                     case 7:
                            a = 0.25;
                            b = 0.50;
                            c = 0.75;
                            d = 1.00;
                            e = -5.0 / 420.0;
                            f = -12.0 / 420.0;
                            g = 16.0 / 420.0;
                            h = 128.0 / 420.0;
                            z = 0.0;

                            xtab[0] = z;
                            xtab[1] = d;
                            xtab[2] = z;
                            xtab[3] = z;
                            xtab[4] = a;
                            xtab[5] = z;
                            xtab[6] = z;
                            xtab[7] = c;
                            xtab[8] = c;
                            xtab[9] = c;
                            xtab[10] = z;
                            xtab[11] = a;
                            xtab[12] = z;
                            xtab[13] = z;
                            xtab[14] = a;
                            xtab[15] = z;
                            xtab[16] = b;
                            xtab[17] = z;
                            xtab[18] = z;
                            xtab[19] = b;
                            xtab[20] = b;
                            xtab[21] = z;
                            xtab[22] = a;
                            xtab[23] = b;
                            xtab[24] = a;
                            xtab[25] = a;
                            xtab[26] = b;
                            xtab[27] = z;
                            xtab[28] = b;
                            xtab[29] = z;
                            xtab[30] = a;
                            xtab[31] = a;
                            xtab[32] = z;
                            xtab[33] = a;
                            xtab[34] = a;

                            ytab[0] = z;
                            ytab[1] = z;
                            ytab[2] = d;
                            ytab[3] = z;
                            ytab[4] = z;
                            ytab[5] = a;
                            ytab[6] = z;
                            ytab[7] = z;
                            ytab[8] = a;
                            ytab[9] = z;
                            ytab[10] = c;
                            ytab[11] = c;
                            ytab[12] = c;
                            ytab[13] = z;
                            ytab[14] = z;
                            ytab[15] = a;
                            ytab[16] = z;
                            ytab[17] = b;
                            ytab[18] = z;
                            ytab[19] = b;
                            ytab[20] = z;
                            ytab[21] = b;
                            ytab[22] = a;
                            ytab[23] = a;
                            ytab[24] = b;
                            ytab[25] = z;
                            ytab[26] = z;
                            ytab[27] = a;
                            ytab[28] = a;
                            ytab[29] = b;
                            ytab[30] = b;
                            ytab[31] = z;
                            ytab[32] = a;
                            ytab[33] = a;
                            ytab[34] = a;

                            ztab[0] = z;
                            ztab[1] = z;
                            ztab[2] = z;
                            ztab[3] = d;
                            ztab[4] = z;
                            ztab[5] = z;
                            ztab[6] = a;
                            ztab[7] = z;
                            ztab[8] = z;
                            ztab[9] = a;
                            ztab[10] = z;
                            ztab[11] = z;
                            ztab[12] = a;
                            ztab[13] = c;
                            ztab[14] = c;
                            ztab[15] = c;
                            ztab[16] = z;
                            ztab[17] = z;
                            ztab[18] = b;
                            ztab[19] = z;
                            ztab[20] = b;
                            ztab[21] = b;
                            ztab[22] = z;
                            ztab[23] = z;
                            ztab[24] = z;
                            ztab[25] = a;
                            ztab[26] = a;
                            ztab[27] = a;
                            ztab[28] = a;
                            ztab[29] = a;
                            ztab[30] = a;
                            ztab[31] = b;
                            ztab[32] = b;
                            ztab[33] = b;
                            ztab[34] = a;

                            weight[0] = e;
                            weight[1] = e;
                            weight[2] = e;
                            weight[3] = e;
                            weight[4] = g;
                            weight[5] = g;
                            weight[6] = g;
                            weight[7] = g;
                            weight[8] = g;
                            weight[9] = g;
                            weight[10] = g;
                            weight[11] = g;
                            weight[12] = g;
                            weight[13] = g;
                            weight[14] = g;
                            weight[15] = g;
                            weight[16] = f;
                            weight[17] = f;
                            weight[18] = f;
                            weight[19] = f;
                            weight[20] = f;
                            weight[21] = f;
                            weight[22] = g;
                            weight[23] = g;
                            weight[24] = g;
                            weight[25] = g;
                            weight[26] = g;
                            weight[27] = g;
                            weight[28] = g;
                            weight[29] = g;
                            weight[30] = g;
                            weight[31] = g;
                            weight[32] = g;
                            weight[33] = g;
                            weight[34] = h;
                            break;
                     //
                     //  Keast Rule of order 11
                     //
                     case 8:
                            a = 0.25;
                            b = 11.0 / 14.0;
                            c = 1.0 / 14.0;
                            d = 0.25 * (1.0 + Math.Sqrt(5.0 / 14.0));
                            e = 0.25 * (1.0 - Math.Sqrt(5.0 / 14.0));
                            f = -74.0 / 5625.0;
                            g = 343.0 / 45000.0;
                            h = 56.0 / 2250.0;

                            xtab[0] = a;
                            xtab[1] = b;
                            xtab[2] = c;
                            xtab[3] = c;
                            xtab[4] = c;
                            xtab[5] = d;
                            xtab[6] = d;
                            xtab[7] = d;
                            xtab[8] = e;
                            xtab[9] = e;
                            xtab[10] = e;

                            ytab[0] = a;
                            ytab[1] = c;
                            ytab[2] = b;
                            ytab[3] = c;
                            ytab[4] = c;
                            ytab[5] = d;
                            ytab[6] = e;
                            ytab[7] = e;
                            ytab[8] = d;
                            ytab[9] = d;
                            ytab[10] = e;

                            ztab[0] = a;
                            ztab[1] = c;
                            ztab[2] = c;
                            ztab[3] = b;
                            ztab[4] = c;
                            ztab[5] = e;
                            ztab[6] = d;
                            ztab[7] = e;
                            ztab[8] = d;
                            ztab[9] = e;
                            ztab[10] = d;

                            weight[0] = f;
                            weight[1] = g;
                            weight[2] = g;
                            weight[3] = g;
                            weight[4] = g;
                            weight[5] = h;
                            weight[6] = h;
                            weight[7] = h;
                            weight[8] = h;
                            weight[9] = h;
                            weight[10] = h;
                            break;
                     default:
                            Console.WriteLine("");
                            Console.WriteLine("TETRA_UNIT_SET - Fatal error!");
                            Console.WriteLine("  Illegal value of RULE = " + rule + "");
                            break;
              }
       }

       public static int tetra_unit_size(int rule)

              //****************************************************************************80
              //
              //  Purpose:
              //
              //    TETRA_UNIT_SIZE sizes quadrature weights and abscissas in the unit tetrahedron.
              //
              //  Integration region:
              //
              //      0 <= X,
              //    and
              //      0 <= Y,
              //    and
              //      0 <= Z, 
              //    and
              //      X + Y + Z <= 1.
              //
              //  Licensing:
              //
              //    This code is distributed under the GNU LGPL license. 
              //
              //  Modified:
              //
              //    19 March 2008
              //
              //  Author:
              //
              //    John Burkardt
              //
              //  Reference:
              //
              //    Hermann Engels,
              //    Numerical Quadrature and Cubature,
              //    Academic Press, 1980,
              //    ISBN: 012238850X,
              //    LC: QA299.3E5.
              //
              //    Patrick Keast,
              //    Moderate Degree Tetrahedral Quadrature Formulas,
              //    Computer Methods in Applied Mechanics and Engineering,
              //    Volume 55, Number 3, May 1986, pages 339-348.
              //
              //    Olgierd Zienkiewicz,
              //    The Finite Element Method,
              //    Sixth Edition,
              //    Butterworth-Heinemann, 2005,
              //    ISBN: 0750663200,
              //    LC: TA640.2.Z54
              //
              //  Parameters:
              //
              //    Input, int RULE, the index of the rule.
              //     1, order 1, precision 0, Newton Cotes formula #0, Zienkiewicz #1.
              //     2, order 4, precision 1, Newton Cotes formula #1.
              //     3, order 4, precision 2, Zienkiewicz #2.
              //     4, order 10, precision 2, Newton Cotes formula #2
              //     5, order 5, precision 3, Zienkiewicz #3.
              //     6, order 8, precision 3, Newton Cotes formula #3.
              //     7, order 35, precision 4, Newton Cotes formula #4.
              //     8, order 11, precision 4, a Keast rule.
              //
              //    Output, int TETRA_UNIT_SET, the order of the rule.
              //
       {
              int order;
              switch (rule)
              {
                     //
                     //  Newton Cotes #0.
                     //
                     case 1:
                            order = 1;
                            break;
                     //
                     //  Newton Cotes #1.
                     //
                     case 2:
                     //
                     //  Zienkiewicz #2.
                     //
                     case 3:
                            order = 4;
                            break;
                     //
                     //  Newton Cotes #2.
                     //
                     case 4:
                            order = 10;
                            break;
                     //
                     //  Zienkiewicz #3.
                     //
                     case 5:
                            order = 5;
                            break;
                     //
                     //  Newton Cotes #3.
                     //  (This is actually formally a 20 point rule, but with 12 zero coefficients//)
                     //
                     case 6:
                            order = 8;
                            break;
                     //
                     //  Newton Cotes #4.
                     //
                     case 7:
                            order = 35;
                            break;
                     //
                     //  Keast Rule of order 11
                     //
                     case 8:
                            order = 11;
                            break;
                     default:
                            Console.WriteLine("");
                            Console.WriteLine("TETRA_UNIT_SIZE - Fatal error!");
                            Console.WriteLine("  Illegal value of RULE = " + rule + "");
                            return 1;
              }

              return order;
       }

       public static double tetra_unit_sum(int setting, Func<int, double, double, double, double> func,
                     int order, double[] xtab, double[] ytab, double[] ztab, double[] weight)

              //****************************************************************************80
              //
              //  Purpose:
              //
              //    TETRA_UNIT_SUM carries out a quadrature rule in the unit tetrahedron in 3D.
              //
              //  Integration region:
              //
              //      0 <= X,
              //    and
              //      0 <= Y,
              //    and
              //      0 <= Z, 
              //    and
              //      X + Y + Z <= 1.
              //
              //  Licensing:
              //
              //    This code is distributed under the GNU LGPL license. 
              //
              //  Modified:
              //
              //    15 March 2008
              //
              //  Author:
              //
              //    John Burkardt
              //
              //  Parameters:
              //
              //    Input, Func< double, double, double, double > func, the name of the 
              //    user supplied function to be integrated.
              //
              //    Input, int ORDER, the order of the rule.
              //
              //    Input, double XTAB[ORDER], YTAB[ORDER], ZTAB[ORDER], the
              //    abscissas.
              //
              //    Input, double WEIGHT[ORDER], the weights.
              //
              //    Output, double RESULT, the approximate integral of the function.
              //
       {
              int i;
              double quad;
              double result;
              double volume;

              quad = 0.0;

              for (i = 0; i < order; i++)
              {
                     quad += weight[i] * func(setting, xtab[i], ytab[i], ztab[i]);
              }

              volume = tetra_unit_volume();
              result = quad * volume;

              return result;
       }

       public static double tetra_unit_volume()

              //****************************************************************************80
              //
              //  Purpose:
              //
              //    TETRA_UNIT_VOLUME returns the volume of the unit tetrahedron.
              //
              //  Discussion:
              //
              //    The integration region is:
              //
              //      0 <= X,
              //      0 <= Y,
              //      0 <= Z, 
              //      X + Y + Z <= 1.
              //
              //  Licensing:
              //
              //    This code is distributed under the GNU LGPL license. 
              //
              //  Modified:
              //
              //    12 March 2008
              //
              //  Author:
              //
              //    John Burkardt
              //
              //  Parameters:
              //
              //    Output, double TETRA_UNIT_VOLUME, the volume.
              //
       {
              double volume;

              volume = 1.0 / 6.0;

              return volume;
       }

       public static double tetra_volume(double[] x, double[] y, double[] z)

              //****************************************************************************80
              //
              //  Purpose:
              //
              //    TETRA_VOLUME computes the volume of a tetrahedron.
              //
              //  Integration region:
              //
              //    Points inside a tetrahedron whose four vertices are given.
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
              //  Parameters:
              //
              //    Input, double X[4], Y[4], Z[4], the vertices.
              //
              //    Output, double TETRA_VOLUME, the volume of the tetrahedron.
              //
       {
              double volume;

              volume = Parallelipiped.parallelipiped_volume_3d(x, y, z);

              volume *= tetra_unit_volume();

              return volume;
       }

}