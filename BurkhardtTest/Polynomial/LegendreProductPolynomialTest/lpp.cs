using System;
using System.Globalization;
using Burkardt.Composition;
using Burkardt.PolynomialNS;
using Burkardt.Uniform;

namespace LegendreProductPolynomialTest;

public static class lppTest
{
    public static void lpp_to_polynomial_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LPP_TO_POLYNOMIAL_TEST tests LPP_TO_POLYNOMIAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int m = 2;
        int o = 0;
        int rank;

        Console.WriteLine("");
        Console.WriteLine("LPP_TO_POLYNOMIAL_TEST:");
        Console.WriteLine("  LPP_TO_POLYNOMIAL is given a Legendre product polynomial");
        Console.WriteLine("  and determines its polynomial representation.");

        Console.WriteLine("");
        Console.WriteLine("  Using spatial dimension M = " + m + "");

        for (rank = 1; rank <= 11; rank++)
        {
            int[] l = Comp.comp_unrank_grlex(m, rank);

            int o_max = 1;
            int i;
            for (i = 0; i < m; i++)
            {
                o_max = o_max * (l[i] + 2) / 2;
            }

            double[] c = new double[o_max];
            int[] e = new int[o_max];

            Legendre.lpp_to_polynomial(m, l, o_max, ref o, ref c, ref e);

            string label = "  LPP #" + rank
                                     + " = L(" + l[0]
                                     + ",X)*L(" + l[1]
                                     + ",Y) = ";

            Console.WriteLine("");
            Polynomial.polynomial_print(m, o, c, e, label);
        }
    }

    public static void lpp_value_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LPP_VALUE_TEST tests LPP_VALUE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        const int m = 3;
        const int n = 1;
        int o = 0;
        int rank;

        Console.WriteLine("");
        Console.WriteLine("LPP_VALUE_TEST:");
        Console.WriteLine("  LPP_VALUE evaluates a Legendre product polynomial.");

        const double xlo = -1.0;
        const double xhi = +1.0;
        int seed = 123456789;
        double[] x = UniformRNG.r8vec_uniform_ab_new(m, xlo, xhi, ref seed);

        Console.WriteLine("");
        string cout = "  Evaluate at X = ";
        for (i = 0; i < m; i++)
        {
            cout += "  " + x[i + 0 * m];
        }

        Console.WriteLine(cout);
        Console.WriteLine("");
        Console.WriteLine("  Rank  I1  I2  I3:  L(I1,X1)*L(I2,X2)*L(I3,X3)    P(X1,X2,X3)");
        Console.WriteLine("");

        for (rank = 1; rank <= 20; rank++)
        {
            int[] l = Comp.comp_unrank_grlex(m, rank);
            //
            //  Evaluate the LPP directly.
            //
            double[] v1 = Legendre.lpp_value(m, n, l, x);
            //
            //  Convert the LPP to a polynomial.
            //
            int o_max = 1;
            for (i = 0; i < m; i++)
            {
                o_max = o_max * (l[i] + 2) / 2;
            }

            double[] c = new double[o_max];
            int[] e = new int[o_max];

            Legendre.lpp_to_polynomial(m, l, o_max, ref o, ref c, ref e);
            //
            //  Evaluate the polynomial.
            //
            double[] v2 = Polynomial.polynomial_value(m, o, c, e, n, x);
            //
            //  Compare results.
            //
            Console.WriteLine(rank.ToString().PadLeft(6) + "  "
                                                         + l[0].ToString().PadLeft(2) + "  "
                                                         + l[1].ToString().PadLeft(2) + "  "
                                                         + l[2].ToString().PadLeft(2) + "  "
                                                         + v1[0].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "  "
                                                         + v2[0].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }
}