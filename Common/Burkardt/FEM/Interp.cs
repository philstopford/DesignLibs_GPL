using System;
using System.Globalization;
using Burkardt.OrderNS;
using Burkardt.PolynomialNS;

namespace Burkardt.FEM;

using Polynomial = Polynomial;
public static class Interp
{
    public static void interp(string code, int element_order, double r, double s,
            double[] ubase, ref double u, ref double dudr, ref double duds)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    INTERP interpolates a quantity in an element from basis node values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 March 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string CODE, identifies the element.
        //    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 
        //    'T3', 'T6' and 'T10'.
        //
        //    Input, int ELEMENT_ORDER, order of the element.
        //
        //    Input, double R, S, the reference coordinates of a point.
        //
        //    Input, double UBASE[ELEMENT_ORDER], the value of the quantity 
        //    at the basis nodes.
        //
        //    Output, ref double U, *DUDR, *DUDS, the interpolated value of the
        //    quantity and its derivatives at the point (R,S).
        //
    {
        int i;

        double[] dtdr = new double[element_order];
        double[] dtds = new double[element_order];
        double[] t = new double[element_order];

        Shape.shape(code, r, s, ref t, ref dtdr, ref dtds);

        u = 0.0;
        for (i = 0; i < element_order; i++)
        {
            u += ubase[i] * t[i];
        }

        dudr = 0.0;
        for (i = 0; i < element_order; i++)
        {
            dudr += ubase[i] * dtdr[i];
        }

        duds = 0.0;
        for (i = 0; i < element_order; i++)
        {
            duds += ubase[i] * dtds[i];
        }
    }

    public static void interp_test(string code)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INTERP_TEST tests the interpolation property of an element.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string CODE, identifies the element.
        //    Legal values include "Q4", "Q8", "Q9", "Q12", "Q16", "QL", 
        //    "T3", "T4", "T6" and "T10".
        //
    {
        double area = 0;
        double dudr = 0;
        double duds = 0;
        int i;
        double r = 0;
        double s = 0;
        const int test_num = 5;
        double u = 0;

        switch (code)
        {
            case "T4":
                Console.WriteLine("");
                Console.WriteLine("INTERP_TEST - Warning!");
                Console.WriteLine("  Skipping test for element \"T4\".");
                return;
        }

        Console.WriteLine("");
        Console.WriteLine("INTERP_TEST for element \"" + code + "\".");

        int element_order = Order.order_code(code);

        Console.WriteLine("  Number of nodes = " + element_order + "");

        double[] node_r = new double [element_order];
        double[] node_s = new double [element_order];
        double[] node_u = new double [element_order];
        int[] rexp = new int [element_order];
        int[] sexp = new int[element_order];
        //
        //  Get the coordinates of the reference nodes.
        //
        NodeReference.node_reference(code, ref node_r, ref node_s, ref area);
        //
        //  Get the monomial exponents for which the element is exact.
        //
        Polynomial.poly(code, ref rexp, ref sexp);

        int seed = 123456789;

        for (i = 0; i < element_order; i++)
        {
            Console.WriteLine("");
            Console.WriteLine("  Interpolate R^" + rexp[i] + " * S^" + sexp[i] + "");
            Console.WriteLine("");
            //
            //  Evaluate R**REXP(I) * S**SEXP(I) at the nodes.  This is our data.
            //
            int node;
            for (node = 0; node < element_order; node++)
            {
                r = node_r[node];
                s = node_s[node];
                double r_factor = rexp[i] switch
                {
                    0 => 1.0,
                    _ => Math.Pow(r, rexp[i])
                };

                double s_factor = sexp[i] switch
                {
                    0 => 1.0,
                    _ => Math.Pow(s, sexp[i])
                };

                node_u[node] = r_factor * s_factor;
                Console.WriteLine("  (R,S,U):     "
                                  + "  " + r.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                  + "  " + s.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                  + "  " + node_u[node].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            }

            //
            //  Now pick random points in the element, and compute the interpolated
            //  value of R**REXP(*) * S**SEXP(I) there.  Mathematically, these
            //  values should be exact.
            //
            int test;
            for (test = 1; test <= test_num; test++)
            {
                Reference.reference_sample(code, ref seed, ref r, ref s);

                Console.WriteLine("");
                Console.WriteLine("  (R,S):"
                                  + "  " + r.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                  + "  " + s.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");

                double u_exact = Math.Pow(r, rexp[i]) * Math.Pow(s, sexp[i]);

                double dudr_exact = rexp[i]
                                    * Math.Pow(r, rexp[i] - 1) * Math.Pow(s, sexp[i]);

                double duds_exact = Math.Pow(r, rexp[i]) * sexp[i]
                                                         * Math.Pow(s, sexp[i] - 1);

                interp(code, element_order, r, s, node_u, ref u, ref dudr, ref duds);

                Console.WriteLine("  (U ,U* ,Error): "
                                  + "  " + u_exact.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                  + "  " + u.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                  + "  " + Math.Abs(u_exact - u).ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");

                Console.WriteLine("  (Ur,Ur*,Error): "
                                  + "  " + dudr_exact.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                  + "  " + dudr.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                  + "  " + Math.Abs(dudr_exact - dudr).ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");

                Console.WriteLine("  (Us,Us*,Error): "
                                  + "  " + duds_exact.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                  + "  " + duds.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                  + "  " + Math.Abs(duds_exact - duds).ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            }
        }
    }

}