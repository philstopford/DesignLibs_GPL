using System;

namespace Burkardt.TriangleNS
{
    public class Triangle
    {
        public static void triangle_unit_set(int rule, double[] xtab, double[] ytab,
        double[] weight )

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    TRIANGLE_UNIT_SET sets a quadrature rule in a unit triangle.
        //
        //  Integration region:
        //
        //    Points (X,Y) such that
        //
        //      0 <= X,
        //      0 <= Y, and
        //      X + Y <= 1.
        //
        //  Graph:
        //
        //      ^
        //    1 | *
        //      | ..
        //    Y | . .
        //      | .  .
        //    0 | *---*
        //      +------->
        //        0 X 1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  References:
        //
        //    H R Schwarz,
        //    Methode der Finiten Elemente,
        //    Teubner Studienbuecher, 1980.
        //
        //    Strang and Fix,
        //    An Analysis of the Finite Element Method,
        //    Prentice Hall, 1973, page 184.
        //
        //    Arthur H Stroud,
        //    Approximate Calculation of Multiple Integrals,
        //    Prentice Hall, 1971.
        //
        //    O C Zienkiewicz,
        //    The Finite Element Method,
        //    McGraw Hill, Third Edition, 1977, page 201.
        //
        //  Parameters:
        //
        //    Input, int RULE, the index of the rule.
        //
        //     1, NORDER =  1, precision 1, Zienkiewicz #1.
        //     2, NORDER =  3, precision 1, the "vertex rule".
        //     3, NORDER =  3, precision 2, Strang and Fix formula #1.
        //     4, NORDER =  3, precision 2, Strang and Fix formula #2, Zienkiewicz #2.
        //     5, NORDER =  4, precision 3, Strang and Fix formula #3, Zienkiewicz #3.
        //     6, NORDER =  6, precision 3, Strang and Fix formula #4.
        //     7, NORDER =  6, precision 3, Stroud formula T2:3-1.
        //     8, NORDER =  6, precision 4, Strang and Fix formula #5.
        //     9, NORDER =  7, precision 4, Strang and Fix formula #6.
        //    10, NORDER =  7, precision 5, Strang and Fix formula #7,
        //        Stroud formula T2:5-1, Zienkiewicz #4, Schwarz Table 2.2.
        //    11, NORDER =  9, precision 6, Strang and Fix formula #8.
        //    12, NORDER = 12, precision 6, Strang and Fix formula #9.
        //    13, ORDER = 13, precision 7, Strang and Fix formula #10.
        //        Note that there is a typographical error in Strang and Fix
        //        which lists the value of the XSI(3) component of the
        //        last generator point as 0.4869... when it should be 0.04869...
        //    14, ORDER =  7, precision ?.
        //    15, ORDER = 16, precision 7, conical product Gauss, Stroud formula T2:7-1.
        //    16, ORDER = 64, precision 15, triangular product Gauss rule.
        //    17, ORDER = 19, precision 8, from CUBTRI, ACM TOMS #584.
        //    18, ORDER = 19, precision 9, from TRIEX, Lyness and Jespersen.
        //    19, ORDER = 28, precision 11, from TRIEX, Lyness and Jespersen.
        //    20, ORDER = 37, precision 13, from ACM TOMS #706.
        //
        //    Output, double XTAB[NORDER], YTAB[NORDER], the abscissas.
        //
        //    Output, double WEIGHT[NORDER], the weights of the rule.
        //
        {
            double a;
            double b;
            double c;
            double d;
            double e;
            double f;
            double g;
            int i;
            int j;
            int k;
            int order2;
            double p;
            double q;
            double r;
            double s;
            double t;
            double u;
            double v;
            double w;
            double w1;
            double w2;
            double w3;
            double w4;
            double w5;
            double w6;
            double w7;
            double w8;
            double w9;
            double[] weight1 = new double[8];
            double[] weight2 = new double[8];
            double wx;
            double x;
            double[] xtab1 = new double[8];
            double[] xtab2 = new double[8];
            double y;
            double z;
            //
            //  1 point, precision 1.
            //
            if (rule == 1)
            {
                xtab[0] = 1.0 / 3.0;
                ytab[0] = 1.0 / 3.0;
                weight[0] = 1.0;
            }
            //
            //  3 points, precision 1, the "vertex rule".
            //
            else if (rule == 2)
            {
                xtab[0] = 1.0;
                xtab[1] = 0.0;
                xtab[2] = 0.0;

                ytab[0] = 0.0;
                ytab[1] = 1.0;
                ytab[2] = 0.0;

                weight[0] = 1.0 / 3.0;
                weight[1] = 1.0 / 3.0;
                weight[2] = 1.0 / 3.0;
            }
            //
            //  3 points, precision 2, Strang and Fix formula #1.
            //
            else if (rule == 3)
            {
                xtab[0] = 4.0 / 6.0;
                xtab[1] = 1.0 / 6.0;
                xtab[2] = 1.0 / 6.0;

                ytab[0] = 1.0 / 6.0;
                ytab[1] = 4.0 / 6.0;
                ytab[2] = 1.0 / 6.0;

                weight[0] = 1.0 / 3.0;
                weight[1] = 1.0 / 3.0;
                weight[2] = 1.0 / 3.0;
            }
            //
            //  3 points, precision 2, Strang and Fix formula #2.
            //
            else if (rule == 4)
            {
                xtab[0] = 0.0;
                xtab[1] = 1.0 / 2.0;
                xtab[2] = 1.0 / 2.0;

                ytab[0] = 1.0 / 2.0;
                ytab[1] = 0.0;
                ytab[2] = 1.0 / 2.0;

                weight[0] = 1.0 / 3.0;
                weight[1] = 1.0 / 3.0;
                weight[2] = 1.0 / 3.0;
            }
            //
            //  4 points, precision 3, Strang and Fix formula #3.
            //
            else if (rule == 5)
            {
                xtab[0] = 10.0 / 30.0;
                xtab[1] = 18.0 / 30.0;
                xtab[2] = 6.0 / 30.0;
                xtab[3] = 6.0 / 30.0;

                ytab[0] = 10.0 / 30.0;
                ytab[1] = 6.0 / 30.0;
                ytab[2] = 18.0 / 30.0;
                ytab[3] = 6.0 / 30.0;

                weight[0] = -27.0 / 48.0;
                weight[1] = 25.0 / 48.0;
                weight[2] = 25.0 / 48.0;
                weight[3] = 25.0 / 48.0;
            }
            //
            //  6 points, precision 3, Strang and Fix formula #4.
            //
            else if (rule == 6)
            {
                xtab[0] = 0.659027622374092;
                xtab[1] = 0.659027622374092;
                xtab[2] = 0.231933368553031;
                xtab[3] = 0.231933368553031;
                xtab[4] = 0.109039009072877;
                xtab[5] = 0.109039009072877;

                ytab[0] = 0.231933368553031;
                ytab[1] = 0.109039009072877;
                ytab[2] = 0.659027622374092;
                ytab[3] = 0.109039009072877;
                ytab[4] = 0.659027622374092;
                ytab[5] = 0.231933368553031;

                weight[0] = 1.0 / 6.0;
                weight[1] = 1.0 / 6.0;
                weight[2] = 1.0 / 6.0;
                weight[3] = 1.0 / 6.0;
                weight[4] = 1.0 / 6.0;
                weight[5] = 1.0 / 6.0;
            }
            //
            //  6 points, precision 3, Stroud T2:3-1.
            //
            else if (rule == 7)
            {
                xtab[0] = 0.0;
                xtab[1] = 0.5;
                xtab[2] = 0.5;
                xtab[3] = 2.0 / 3.0;
                xtab[4] = 1.0 / 6.0;
                xtab[5] = 1.0 / 6.0;

                ytab[0] = 0.5;
                ytab[1] = 0.0;
                ytab[2] = 0.5;
                ytab[3] = 1.0 / 6.0;
                ytab[4] = 2.0 / 3.0;
                ytab[5] = 1.0 / 6.0;

                weight[0] = 1.0 / 30.0;
                weight[1] = 1.0 / 30.0;
                weight[2] = 1.0 / 30.0;
                weight[3] = 3.0 / 10.0;
                weight[4] = 3.0 / 10.0;
                weight[5] = 3.0 / 10.0;
            }
            //
            //  6 points, precision 4, Strang and Fix, formula #5.
            //
            else if (rule == 8)
            {
                xtab[0] = 0.816847572980459;
                xtab[1] = 0.091576213509771;
                xtab[2] = 0.091576213509771;
                xtab[3] = 0.108103018168070;
                xtab[4] = 0.445948490915965;
                xtab[5] = 0.445948490915965;

                ytab[0] = 0.091576213509771;
                ytab[1] = 0.816847572980459;
                ytab[2] = 0.091576213509771;
                ytab[3] = 0.445948490915965;
                ytab[4] = 0.108103018168070;
                ytab[5] = 0.445948490915965;

                weight[0] = 0.109951743655322;
                weight[1] = 0.109951743655322;
                weight[2] = 0.109951743655322;
                weight[3] = 0.223381589678011;
                weight[4] = 0.223381589678011;
                weight[5] = 0.223381589678011;
            }
            //
            //  7 points, precision 4, Strang and Fix formula #6.
            //
            else if (rule == 9)
            {
                xtab[0] = 1.0 / 3.0;
                xtab[1] = 0.736712498968435;
                xtab[2] = 0.736712498968435;
                xtab[3] = 0.237932366472434;
                xtab[4] = 0.237932366472434;
                xtab[5] = 0.025355134551932;
                xtab[6] = 0.025355134551932;

                ytab[0] = 1.0 / 3.0;
                ytab[1] = 0.237932366472434;
                ytab[2] = 0.025355134551932;
                ytab[3] = 0.736712498968435;
                ytab[4] = 0.025355134551932;
                ytab[5] = 0.736712498968435;
                ytab[6] = 0.237932366472434;

                weight[0] = 0.375;
                weight[1] = 0.1041666666666667;
                weight[2] = 0.1041666666666667;
                weight[3] = 0.1041666666666667;
                weight[4] = 0.1041666666666667;
                weight[5] = 0.1041666666666667;
                weight[6] = 0.1041666666666667;
            }
            //
            //  7 points, precision 5, Strang and Fix formula #7, Stroud T2:5-1
            //
            else if (rule == 10)
            {
                xtab[0] = 1.0 / 3.0;
                xtab[1] = (9.0 + 2.0 * Math.Sqrt(15.0)) / 21.0;
                xtab[2] = (6.0 - Math.Sqrt(15.0)) / 21.0;
                xtab[3] = (6.0 - Math.Sqrt(15.0)) / 21.0;
                xtab[4] = (9.0 - 2.0 * Math.Sqrt(15.0)) / 21.0;
                xtab[5] = (6.0 + Math.Sqrt(15.0)) / 21.0;
                xtab[6] = (6.0 + Math.Sqrt(15.0)) / 21.0;

                ytab[0] = 1.0 / 3.0;
                ytab[1] = (6.0 - Math.Sqrt(15.0)) / 21.0;
                ytab[2] = (9.0 + 2.0 * Math.Sqrt(15.0)) / 21.0;
                ytab[3] = (6.0 - Math.Sqrt(15.0)) / 21.0;
                ytab[4] = (6.0 + Math.Sqrt(15.0)) / 21.0;
                ytab[5] = (9.0 - 2.0 * Math.Sqrt(15.0)) / 21.0;
                ytab[6] = (6.0 + Math.Sqrt(15.0)) / 21.0;

                weight[0] = 0.225;
                weight[1] = (155.0 - Math.Sqrt(15.0)) / 1200.0;
                weight[2] = (155.0 - Math.Sqrt(15.0)) / 1200.0;
                weight[3] = (155.0 - Math.Sqrt(15.0)) / 1200.0;
                weight[4] = (155.0 + Math.Sqrt(15.0)) / 1200.0;
                weight[5] = (155.0 + Math.Sqrt(15.0)) / 1200.0;
                weight[6] = (155.0 + Math.Sqrt(15.0)) / 1200.0;
            }
            //
            //  9 points, precision 6, Strang and Fix formula #8.
            //
            else if (rule == 11)
            {
                xtab[0] = 0.124949503233232;
                xtab[1] = 0.437525248383384;
                xtab[2] = 0.437525248383384;
                xtab[3] = 0.797112651860071;
                xtab[4] = 0.797112651860071;
                xtab[5] = 0.165409927389841;
                xtab[6] = 0.165409927389841;
                xtab[7] = 0.037477420750088;
                xtab[8] = 0.037477420750088;

                ytab[0] = 0.437525248383384;
                ytab[1] = 0.124949503233232;
                ytab[2] = 0.437525248383384;
                ytab[3] = 0.165409927389841;
                ytab[4] = 0.037477420750088;
                ytab[5] = 0.797112651860071;
                ytab[6] = 0.037477420750088;
                ytab[7] = 0.797112651860071;
                ytab[8] = 0.165409927389841;

                weight[0] = 0.205950504760887;
                weight[1] = 0.205950504760887;
                weight[2] = 0.205950504760887;
                weight[3] = 0.063691414286223;
                weight[4] = 0.063691414286223;
                weight[5] = 0.063691414286223;
                weight[6] = 0.063691414286223;
                weight[7] = 0.063691414286223;
                weight[8] = 0.063691414286223;
            }
            //
            //  12 points, precision 6, Strang and Fix, formula #9.
            //
            else if (rule == 12)
            {
                xtab[0] = 0.873821971016996;
                xtab[1] = 0.063089014491502;
                xtab[2] = 0.063089014491502;
                xtab[3] = 0.249286745170910;
                xtab[4] = 0.501426509658179;
                xtab[5] = 0.249286745170910;
                xtab[6] = 0.636502499121399;
                xtab[7] = 0.636502499121399;
                xtab[8] = 0.310352451033785;
                xtab[9] = 0.310352451033785;
                xtab[10] = 0.053145049844816;
                xtab[11] = 0.053145049844816;

                ytab[0] = 0.063089014491502;
                ytab[1] = 0.873821971016996;
                ytab[2] = 0.063089014491502;
                ytab[3] = 0.501426509658179;
                ytab[4] = 0.249286745170910;
                ytab[5] = 0.249286745170910;
                ytab[6] = 0.310352451033785;
                ytab[7] = 0.053145049844816;
                ytab[8] = 0.636502499121399;
                ytab[9] = 0.053145049844816;
                ytab[10] = 0.636502499121399;
                ytab[11] = 0.310352451033785;

                weight[0] = 0.050844906370207;
                weight[1] = 0.050844906370207;
                weight[2] = 0.050844906370207;
                weight[3] = 0.116786275726379;
                weight[4] = 0.116786275726379;
                weight[5] = 0.116786275726379;
                weight[6] = 0.082851075618374;
                weight[7] = 0.082851075618374;
                weight[8] = 0.082851075618374;
                weight[9] = 0.082851075618374;
                weight[10] = 0.082851075618374;
                weight[11] = 0.082851075618374;
            }
            //
            //  13 points, precision 7, Strang and Fix, formula #10.
            //
            else if (rule == 13)
            {
                xtab[0] = 0.479308067841923;
                xtab[1] = 0.260345966079038;
                xtab[2] = 0.260345966079038;
                xtab[3] = 0.869739794195568;
                xtab[4] = 0.065130102902216;
                xtab[5] = 0.065130102902216;
                xtab[6] = 0.638444188569809;
                xtab[7] = 0.638444188569809;
                xtab[8] = 0.312865496004875;
                xtab[9] = 0.312865496004875;
                xtab[10] = 0.048690315425316;
                xtab[11] = 0.048690315425316;
                xtab[12] = 1.0 / 3.0;

                ytab[0] = 0.260345966079038;
                ytab[1] = 0.479308067841923;
                ytab[2] = 0.260345966079038;
                ytab[3] = 0.065130102902216;
                ytab[4] = 0.869739794195568;
                ytab[5] = 0.065130102902216;
                ytab[6] = 0.312865496004875;
                ytab[7] = 0.048690315425316;
                ytab[8] = 0.638444188569809;
                ytab[9] = 0.048690315425316;
                ytab[10] = 0.638444188569809;
                ytab[11] = 0.312865496004875;
                ytab[12] = 1.0 / 3.0;

                weight[0] = 0.175615257433204;
                weight[1] = 0.175615257433204;
                weight[2] = 0.175615257433204;
                weight[3] = 0.053347235608839;
                weight[4] = 0.053347235608839;
                weight[5] = 0.053347235608839;
                weight[6] = 0.077113760890257;
                weight[7] = 0.077113760890257;
                weight[8] = 0.077113760890257;
                weight[9] = 0.077113760890257;
                weight[10] = 0.077113760890257;
                weight[11] = 0.077113760890257;
                weight[12] = -0.149570044467670;
            }
            //
            //  7 points, precision ?.
            //
            else if (rule == 14)
            {
                a = 1.0 / 3.0;
                b = 1.0;
                c = 0.5;
                z = 0.0;

                u = 27.0 / 60.0;
                v = 3.0 / 60.0;
                w = 8.0 / 60.0;

                xtab[0] = a;
                xtab[1] = b;
                xtab[2] = z;
                xtab[3] = z;
                xtab[4] = z;
                xtab[5] = c;
                xtab[6] = c;

                ytab[0] = a;
                ytab[1] = z;
                ytab[2] = b;
                ytab[3] = z;
                ytab[4] = c;
                ytab[5] = z;
                ytab[6] = c;

                weight[0] = u;
                weight[1] = v;
                weight[2] = v;
                weight[3] = v;
                weight[4] = w;
                weight[5] = w;
                weight[6] = w;
                weight[7] = w;
            }
            //
            //  16 points.
            //
            else if (rule == 15)
            {
                //
                //  Legendre rule of order 4.
                //
                order2 = 4;

                xtab1[0] = -0.861136311594052575223946488893;
                xtab1[1] = -0.339981043584856264802665759103;
                xtab1[2] = 0.339981043584856264802665759103;
                xtab1[3] = 0.861136311594052575223946488893;

                weight1[0] = 0.347854845137453857373063949222;
                weight1[1] = 0.652145154862546142626936050778;
                weight1[2] = 0.652145154862546142626936050778;
                weight1[3] = 0.347854845137453857373063949222;

                for (i = 0; i < order2; i++)
                {
                    xtab1[i] = 0.5 * (xtab1[i] + 1.0);
                }

                weight2[0] = 0.1355069134;
                weight2[1] = 0.2034645680;
                weight2[2] = 0.1298475476;
                weight2[3] = 0.0311809709;

                xtab2[0] = 0.0571041961;
                xtab2[1] = 0.2768430136;
                xtab2[2] = 0.5835904324;
                xtab2[3] = 0.8602401357;

                k = 0;
                for (i = 0; i < order2; i++)
                {
                    for (j = 0; j < order2; j++)
                    {
                        xtab[k] = xtab2[j];
                        ytab[k] = xtab1[i] * (1.0 - xtab2[j]);
                        weight[k] = weight1[i] * weight2[j];
                        k = k + 1;
                    }
                }
            }
            //
            //  64 points, precision 15.
            //
            else if (rule == 16)
            {
                //
                //  Legendre rule of order 8.
                //
                order2 = 8;

                xtab1[0] = -0.960289856497536231683560868569;
                xtab1[1] = -0.796666477413626739591553936476;
                xtab1[2] = -0.525532409916328985817739049189;
                xtab1[3] = -0.183434642495649804939476142360;
                xtab1[4] = 0.183434642495649804939476142360;
                xtab1[5] = 0.525532409916328985817739049189;
                xtab1[6] = 0.796666477413626739591553936476;
                xtab1[7] = 0.960289856497536231683560868569;

                weight1[0] = 0.101228536290376259152531354310;
                weight1[1] = 0.222381034453374470544355994426;
                weight1[2] = 0.313706645877887287337962201987;
                weight1[3] = 0.362683783378361982965150449277;
                weight1[4] = 0.362683783378361982965150449277;
                weight1[5] = 0.313706645877887287337962201987;
                weight1[6] = 0.222381034453374470544355994426;
                weight1[7] = 0.101228536290376259152531354310;

                weight2[0] = 0.00329519144;
                weight2[1] = 0.01784290266;
                weight2[2] = 0.04543931950;
                weight2[3] = 0.07919959949;
                weight2[4] = 0.10604735944;
                weight2[5] = 0.11250579947;
                weight2[6] = 0.09111902364;
                weight2[7] = 0.04455080436;

                xtab2[0] = 0.04463395529;
                xtab2[1] = 0.14436625704;
                xtab2[2] = 0.28682475714;
                xtab2[3] = 0.45481331520;
                xtab2[4] = 0.62806783542;
                xtab2[5] = 0.78569152060;
                xtab2[6] = 0.90867639210;
                xtab2[7] = 0.98222008485;

                k = 0;
                for (j = 0; j < order2; j++)
                {
                    for (i = 0; i < order2; i++)
                    {
                        xtab[k] = 1.0 - xtab2[j];
                        ytab[k] = 0.5 * (1.0 + xtab1[i]) * xtab2[j];
                        weight[k] = weight1[i] * weight2[j];
                        k = k + 1;
                    }
                }
            }
            //
            //  19 points, precision 8.
            //
            else if (rule == 17)
            {
                a = 1.0 / 3.0;
                b = (9.0 + 2.0 * Math.Sqrt(15.0)) / 21.0;
                c = (6.0 - Math.Sqrt(15.0)) / 21.0;
                d = (9.0 - 2.0 * Math.Sqrt(15.0)) / 21.0;
                e = (6.0 + Math.Sqrt(15.0)) / 21.0;
                f = (40.0 - 10.0 * Math.Sqrt(15.0)
                     + 10.0 * Math.Sqrt(7.0) + 2.0 * Math.Sqrt(105.0)) / 90.0;
                g = (25.0 + 5.0 * Math.Sqrt(15.0)
                     - 5.0 * Math.Sqrt(7.0) - Math.Sqrt(105.0)) / 90.0;
                p = (40.0 + 10.0 * Math.Sqrt(15.0)
                          + 10.0 * Math.Sqrt(7.0) - 2.0 * Math.Sqrt(105.0)) / 90.0;
                q = (25.0 - 5.0 * Math.Sqrt(15.0)
                          - 5.0 * Math.Sqrt(7.0) + Math.Sqrt(105.0)) / 90.0;
                r = (40.0 + 10.0 * Math.Sqrt(7.0)) / 90.0;
                s = (25.0 + 5.0 * Math.Sqrt(15.0) - 5.0 * Math.Sqrt(7.0)
                                             - Math.Sqrt(105.0)) / 90.0;
                t = (25.0 - 5.0 * Math.Sqrt(15.0) - 5.0 * Math.Sqrt(7.0)
                     + Math.Sqrt(105.0)) / 90.0;

                w1 = (7137.0 - 1800.0 * Math.Sqrt(7.0)) / 62720.0;
                w2 = -9301697.0 / 4695040.0 - 13517313.0 * Math.Sqrt(15.0)
                    / 23475200.0 + 764885.0 * Math.Sqrt(7.0) / 939008.0
                                 + 198763.0 * Math.Sqrt(105.0) / 939008.0;
                w2 = w2 / 3.0;
                w3 = -9301697.0 / 4695040.0 + 13517313.0 * Math.Sqrt(15.0)
                                            / 23475200.0
                                            + 764885.0 * Math.Sqrt(7.0) / 939008.0
                     - 198763.0 * Math.Sqrt(105.0) / 939008.0;
                w3 = w3 / 3.0;
                w4 = (102791225.0 - 23876225.0 * Math.Sqrt(15.0)
                                  - 34500875.0 * Math.Sqrt(7.0)
                      + 9914825.0 * Math.Sqrt(105.0)) / 59157504.0;
                w4 = w4 / 3.0;
                w5 = (102791225.0 + 23876225.0 * Math.Sqrt(15.0)
                      - 34500875.0 * Math.Sqrt(7.0)
                      - 9914825 * Math.Sqrt(105.0)) / 59157504.0;
                w5 = w5 / 3.0;
                w6 = (11075.0 - 3500.0 * Math.Sqrt(7.0)) / 8064.0;
                w6 = w6 / 6.0;

                xtab[0] = a;
                xtab[1] = b;
                xtab[2] = c;
                xtab[3] = c;
                xtab[4] = d;
                xtab[5] = e;
                xtab[6] = e;
                xtab[7] = f;
                xtab[8] = g;
                xtab[9] = g;
                xtab[10] = p;
                xtab[11] = q;
                xtab[12] = q;
                xtab[13] = r;
                xtab[14] = r;
                xtab[15] = s;
                xtab[16] = s;
                xtab[17] = t;
                xtab[18] = t;

                ytab[0] = a;
                ytab[1] = c;
                ytab[2] = b;
                ytab[3] = c;
                ytab[4] = e;
                ytab[5] = d;
                ytab[6] = e;
                ytab[7] = g;
                ytab[8] = f;
                ytab[9] = g;
                ytab[10] = q;
                ytab[11] = p;
                ytab[12] = q;
                ytab[13] = s;
                ytab[14] = t;
                ytab[15] = r;
                ytab[16] = t;
                ytab[17] = r;
                ytab[18] = s;

                weight[0] = w1;
                weight[1] = w2;
                weight[2] = w2;
                weight[3] = w2;
                weight[4] = w3;
                weight[5] = w3;
                weight[6] = w3;
                weight[7] = w4;
                weight[8] = w4;
                weight[9] = w4;
                weight[10] = w5;
                weight[11] = w5;
                weight[12] = w5;
                weight[13] = w6;
                weight[14] = w6;
                weight[15] = w6;
                weight[16] = w6;
                weight[17] = w6;
                weight[18] = w6;
            }
            //
            //  19 points, precision 9.
            //
            else if (rule == 18)
            {
                a = 1.0 / 3.0;
                b = 0.02063496160252593;
                c = 0.4896825191987370;
                d = 0.1258208170141290;
                e = 0.4370895914929355;
                f = 0.6235929287619356;
                g = 0.1882035356190322;
                r = 0.9105409732110941;
                s = 0.04472951339445297;
                t = 0.7411985987844980;
                u = 0.03683841205473626;
                v = 0.22196288916076574;

                w1 = 0.09713579628279610;
                w2 = 0.03133470022713983;
                w3 = 0.07782754100477543;
                w4 = 0.07964773892720910;
                w5 = 0.02557767565869810;
                w6 = 0.04328353937728940;

                xtab[0] = a;
                xtab[1] = b;
                xtab[2] = c;
                xtab[3] = c;
                xtab[4] = d;
                xtab[5] = e;
                xtab[6] = e;
                xtab[7] = f;
                xtab[8] = g;
                xtab[9] = g;
                xtab[10] = r;
                xtab[11] = s;
                xtab[12] = s;
                xtab[13] = t;
                xtab[14] = t;
                xtab[15] = u;
                xtab[16] = u;
                xtab[17] = v;
                xtab[18] = v;

                ytab[0] = a;
                ytab[1] = c;
                ytab[2] = b;
                ytab[3] = c;
                ytab[4] = e;
                ytab[5] = d;
                ytab[6] = e;
                ytab[7] = g;
                ytab[8] = f;
                ytab[9] = g;
                ytab[10] = s;
                ytab[11] = r;
                ytab[12] = s;
                ytab[13] = u;
                ytab[14] = v;
                ytab[15] = t;
                ytab[16] = v;
                ytab[17] = t;
                ytab[18] = u;

                weight[0] = w1;
                weight[1] = w2;
                weight[2] = w2;
                weight[3] = w2;
                weight[4] = w3;
                weight[5] = w3;
                weight[6] = w3;
                weight[7] = w4;
                weight[8] = w4;
                weight[9] = w4;
                weight[10] = w5;
                weight[11] = w5;
                weight[12] = w5;
                weight[13] = w6;
                weight[14] = w6;
                weight[15] = w6;
                weight[16] = w6;
                weight[17] = w6;
                weight[18] = w6;
            }
            //
            //  28 points, precision 11.
            //
            else if (rule == 19)
            {
                a = 1.0 / 3.0;
                b = 0.9480217181434233;
                c = 0.02598914092828833;
                d = 0.8114249947041546;
                e = 0.09428750264792270;
                f = 0.01072644996557060;
                g = 0.4946367750172147;
                p = 0.5853132347709715;
                q = 0.2073433826145142;
                r = 0.1221843885990187;
                s = 0.4389078057004907;
                t = 0.6779376548825902;
                u = 0.04484167758913055;
                v = 0.27722066752827925;
                w = 0.8588702812826364;
                x = 0.0;
                y = 0.1411297187173636;

                w1 = 0.08797730116222190;
                w2 = 0.008744311553736190;
                w3 = 0.03808157199393533;
                w4 = 0.01885544805613125;
                w5 = 0.07215969754474100;
                w6 = 0.06932913870553720;
                w7 = 0.04105631542928860;
                w8 = 0.007362383783300573;

                xtab[0] = a;
                xtab[1] = b;
                xtab[2] = c;
                xtab[3] = c;
                xtab[4] = d;
                xtab[5] = e;
                xtab[6] = e;
                xtab[7] = f;
                xtab[8] = g;
                xtab[9] = g;
                xtab[10] = p;
                xtab[11] = q;
                xtab[12] = q;
                xtab[13] = r;
                xtab[14] = s;
                xtab[15] = s;
                xtab[16] = t;
                xtab[17] = t;
                xtab[18] = u;
                xtab[19] = u;
                xtab[20] = v;
                xtab[21] = v;
                xtab[22] = w;
                xtab[23] = w;
                xtab[24] = x;
                xtab[25] = x;
                xtab[26] = y;
                xtab[27] = y;

                ytab[0] = a;
                ytab[1] = c;
                ytab[2] = b;
                ytab[3] = c;
                ytab[4] = e;
                ytab[5] = d;
                ytab[6] = e;
                ytab[7] = g;
                ytab[8] = f;
                ytab[9] = g;
                ytab[10] = q;
                ytab[11] = p;
                ytab[12] = q;
                ytab[13] = s;
                ytab[14] = r;
                ytab[15] = s;
                ytab[16] = u;
                ytab[17] = v;
                ytab[18] = t;
                ytab[19] = v;
                ytab[20] = t;
                ytab[21] = u;
                ytab[22] = x;
                ytab[23] = y;
                ytab[24] = w;
                ytab[25] = y;
                ytab[26] = w;
                ytab[27] = x;

                weight[0] = w1;
                weight[1] = w2;
                weight[2] = w2;
                weight[3] = w2;
                weight[4] = w3;
                weight[5] = w3;
                weight[6] = w3;
                weight[7] = w4;
                weight[8] = w4;
                weight[9] = w4;
                weight[10] = w5;
                weight[11] = w5;
                weight[12] = w5;
                weight[13] = w6;
                weight[14] = w6;
                weight[15] = w6;
                weight[16] = w7;
                weight[17] = w7;
                weight[18] = w7;
                weight[19] = w7;
                weight[20] = w7;
                weight[21] = w7;
                weight[22] = w8;
                weight[23] = w8;
                weight[24] = w8;
                weight[25] = w8;
                weight[26] = w8;
                weight[27] = w8;
            }
            //
            //  37 points, precision 13.
            //
            else if (rule == 20)
            {
                a = 1.0 / 3.0;
                b = 0.950275662924105565450352089520;
                c = 0.024862168537947217274823955239;
                d = 0.171614914923835347556304795551;
                e = 0.414192542538082326221847602214;
                f = 0.539412243677190440263092985511;
                g = 0.230293878161404779868453507244;

                w1 = 0.051739766065744133555179145422;
                w2 = 0.008007799555564801597804123460;
                w3 = 0.046868898981821644823226732071;
                w4 = 0.046590940183976487960361770070;
                w5 = 0.031016943313796381407646220131;
                w6 = 0.010791612736631273623178240136;
                w7 = 0.032195534242431618819414482205;
                w8 = 0.015445834210701583817692900053;
                w9 = 0.017822989923178661888748319485;
                wx = 0.037038683681384627918546472190;

                xtab[0] = a;
                xtab[1] = b;
                xtab[2] = c;
                xtab[3] = c;
                xtab[4] = d;
                xtab[5] = e;
                xtab[6] = e;
                xtab[7] = f;
                xtab[8] = g;
                xtab[9] = g;

                ytab[0] = a;
                ytab[1] = c;
                ytab[2] = b;
                ytab[3] = c;
                ytab[4] = e;
                ytab[5] = d;
                ytab[6] = e;
                ytab[7] = g;
                ytab[8] = f;
                ytab[9] = g;

                weight[0] = w1;
                weight[1] = w2;
                weight[2] = w2;
                weight[3] = w2;
                weight[4] = w3;
                weight[5] = w3;
                weight[6] = w3;
                weight[7] = w4;
                weight[8] = w4;
                weight[9] = w4;
                weight[10] = w5;
                weight[11] = w5;
                weight[12] = w5;
                weight[13] = w6;
                weight[14] = w6;
                weight[15] = w6;
                weight[16] = w7;
                weight[17] = w7;
                weight[18] = w7;
                weight[19] = w8;
                weight[20] = w8;
                weight[21] = w8;
                weight[22] = w8;
                weight[23] = w8;
                weight[24] = w8;
                weight[25] = w9;
                weight[26] = w9;
                weight[27] = w9;
                weight[28] = w9;
                weight[29] = w9;
                weight[30] = w9;
                weight[31] = wx;
                weight[32] = wx;
                weight[33] = wx;
                weight[34] = wx;
                weight[35] = wx;
                weight[36] = wx;

                a = 0.772160036676532561750285570113;
                b = 0.113919981661733719124857214943;

                xtab[10] = a;
                ytab[10] = b;

                xtab[11] = b;
                ytab[11] = a;

                xtab[12] = b;
                ytab[12] = b;

                a = 0.009085399949835353883572964740;
                b = 0.495457300025082323058213517632;

                xtab[13] = a;
                ytab[13] = b;

                xtab[14] = b;
                ytab[14] = a;

                xtab[15] = b;
                ytab[15] = b;

                a = 0.062277290305886993497083640527;
                b = 0.468861354847056503251458179727;

                xtab[16] = a;
                ytab[16] = b;

                xtab[17] = b;
                ytab[17] = a;

                xtab[18] = b;
                ytab[18] = b;

                a = 0.022076289653624405142446876931;
                b = 0.851306504174348550389457672223;
                c = 1.0 - a - b;

                xtab[19] = a;
                ytab[19] = b;

                xtab[20] = a;
                ytab[20] = c;

                xtab[21] = b;
                ytab[21] = a;

                xtab[22] = b;
                ytab[22] = c;

                xtab[23] = c;
                ytab[23] = a;

                xtab[24] = c;
                ytab[24] = b;

                a = 0.018620522802520968955913511549;
                b = 0.689441970728591295496647976487;
                c = 1.0 - a - b;

                xtab[25] = a;
                ytab[25] = b;

                xtab[26] = a;
                ytab[26] = c;

                xtab[27] = b;
                ytab[27] = a;

                xtab[28] = b;
                ytab[28] = c;

                xtab[29] = c;
                ytab[29] = a;

                xtab[30] = c;
                ytab[30] = b;

                a = 0.096506481292159228736516560903;
                b = 0.635867859433872768286976979827;
                c = 1.0 - a - b;

                xtab[31] = a;
                ytab[31] = b;

                xtab[32] = a;
                ytab[32] = c;

                xtab[33] = b;
                ytab[33] = a;

                xtab[34] = b;
                ytab[34] = c;

                xtab[35] = c;
                ytab[35] = a;

                xtab[36] = c;
                ytab[36] = b;
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("TRIANGLE_UNIT_SET - Fatal error!");
                Console.WriteLine("  Illegal value of RULE = " + rule + "");
            }
        }

        public static int triangle_unit_size(int rule)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGLE_UNIT_SIZE returns the "size" of a unit triangle quadrature rule.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    06 September 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Gilbert Strang, George Fix,
            //    An Analysis of the Finite Element Method,
            //    Prentice Hall, 1973,
            //    TA335.S77.
            //
            //    Olgierd Zienkiewicz,
            //    The Finite Element Method,
            //    McGraw Hill, Third Edition, 1977, page 202.
            //
            //  Parameters:
            //
            //    Input, int RULE, the index of the rule.
            //
            //     1, ORDER =  1, precision 1, Zienkiewicz #1.
            //     2, ORDER =  3, precision 1, the "vertex rule".
            //     3, ORDER =  3, precision 2, Strang and Fix formula #1.
            //     4, ORDER =  3, precision 2, Strang and Fix formula #2, Zienkiewicz #2.
            //     5, ORDER =  4, precision 3, Strang and Fix formula #3, Zienkiewicz #3.
            //     6, ORDER =  6, precision 3, Strang and Fix formula #4.
            //     7, ORDER =  6, precision 3, Stroud formula T2:3-1.
            //     8, ORDER =  6, precision 4, Strang and Fix formula #5.
            //     9, ORDER =  7, precision 4, Strang and Fix formula #6.
            //    10, ORDER =  7, precision 5, Strang and Fix formula #7,
            //        Stroud formula T2:5-1, Zienkiewicz #4, Schwarz Table 2.2.
            //    11, ORDER =  9, precision 6, Strang and Fix formula #8.
            //    12, ORDER = 12, precision 6, Strang and Fix formula #9.
            //    13, ORDER = 13, precision 7, Strang and Fix formula #10.
            //    14, ORDER =  7, precision ?.
            //    15, ORDER = 16, precision 7, conical product Gauss, Stroud formula T2:7-1.
            //    16, ORDER = 64, precision 15, triangular product Gauss rule.
            //    17, ORDER = 19, precision 8, from CUBTRI, ACM TOMS #584.
            //    18, ORDER = 19, precision 9, from TRIEX, Lyness and Jespersen.
            //    19, ORDER = 28, precision 11, from TRIEX, Lyness and Jespersen.
            //    20, ORDER = 37, precision 13, from ACM TOMS #706.
            //
            //    Output, int TRIANGLE_UNIT_SIZE, the order of the rule.
            //
        {
            int value;

            if (rule == 1)
            {
                value = 1;
            }
            else if (rule == 2)
            {
                value = 3;
            }
            else if (rule == 3)
            {
                value = 3;
            }
            else if (rule == 4)
            {
                value = 3;
            }
            else if (rule == 5)
            {
                value = 4;
            }
            else if (rule == 6)
            {
                value = 6;
            }
            else if (rule == 7)
            {
                value = 6;
            }
            else if (rule == 8)
            {
                value = 6;
            }
            else if (rule == 9)
            {
                value = 7;
            }
            else if (rule == 10)
            {
                value = 7;
            }
            else if (rule == 11)
            {
                value = 9;
            }
            else if (rule == 12)
            {
                value = 12;
            }
            else if (rule == 13)
            {
                value = 13;
            }
            else if (rule == 14)
            {
                value = 7;
            }
            else if (rule == 15)
            {
                value = 16;
            }
            else if (rule == 16)
            {
                value = 64;
            }
            else if (rule == 17)
            {
                value = 19;
            }
            else if (rule == 18)
            {
                value = 19;
            }
            else if (rule == 19)
            {
                value = 28;
            }
            else if (rule == 20)
            {
                value = 37;
            }
            else
            {
                value = -1;
            }

            return value;
        }

    }
}