using System;
using Burkardt.PolynomialNS;

namespace Burkardt.FEM
{
    using Polynomial = Burkardt.PolynomialNS.Polynomial;
    public class Interp
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
            double[] dtdr;
            double[] dtds;
            int i;
            double[] t;

            dtdr = new double[element_order];
            dtds = new double[element_order];
            t = new double[element_order];

            Shape.shape(code, r, s, ref t, ref dtdr, ref dtds);

            u = 0.0;
            for (i = 0; i < element_order; i++)
            {
                u = u + ubase[i] * t[i];
            }

            dudr = 0.0;
            for (i = 0; i < element_order; i++)
            {
                dudr = dudr + ubase[i] * dtdr[i];
            }

            duds = 0.0;
            for (i = 0; i < element_order; i++)
            {
                duds = duds + ubase[i] * dtds[i];
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
            double dudr_exact;
            double duds = 0;
            double duds_exact;
            int element_order;
            int i;
            int node;
            double[] node_r;
            double[] node_s;
            double[] node_u;
            double r = 0;
            double r_factor;
            int[] rexp;
            double s = 0;
            double s_factor;
            int[] sexp;
            int seed;
            int test;
            int test_num = 5;
            double u = 0;
            double u_exact;

            if (code == "T4")
            {
                Console.WriteLine("");
                Console.WriteLine("INTERP_TEST - Warning!");
                Console.WriteLine("  Skipping test for element \"T4\".");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("INTERP_TEST for element \"" + code + "\".");

            element_order = Order.order_code(code);

            Console.WriteLine("  Number of nodes = " + element_order + "");

            node_r = new double [element_order];
            node_s = new double [element_order];
            node_u = new double [element_order];
            rexp = new int [element_order];
            sexp = new int[element_order];
            //
            //  Get the coordinates of the reference nodes.
            //
            NodeReference.node_reference(code, ref node_r, ref node_s, ref area);
            //
            //  Get the monomial exponents for which the element is exact.
            //
            Polynomial.poly(code, ref rexp, ref sexp);

            seed = 123456789;

            for (i = 0; i < element_order; i++)
            {
                Console.WriteLine("");
                Console.WriteLine("  Interpolate R^" + rexp[i] + " * S^" + sexp[i] + "");
                Console.WriteLine("");
                //
                //  Evaluate R**REXP(I) * S**SEXP(I) at the nodes.  This is our data.
                //
                for (node = 0; node < element_order; node++)
                {
                    r = node_r[node];
                    s = node_s[node];
                    if (rexp[i] == 0)
                    {
                        r_factor = 1.0;
                    }
                    else
                    {
                        r_factor = Math.Pow(r, rexp[i]);
                    }

                    if (sexp[i] == 0)
                    {
                        s_factor = 1.0;
                    }
                    else
                    {
                        s_factor = Math.Pow(s, sexp[i]);
                    }

                    node_u[node] = r_factor * s_factor;
                    Console.WriteLine("  (R,S,U):     "
                                      + "  " + r.ToString().PadLeft(12)
                                      + "  " + s.ToString().PadLeft(12)
                                      + "  " + node_u[node].ToString().PadLeft(12) + "");
                }

                //
                //  Now pick random points in the element, and compute the interpolated
                //  value of R**REXP(*) * S**SEXP(I) there.  Mathematically, these
                //  values should be exact.
                //
                for (test = 1; test <= test_num; test++)
                {
                    Reference.reference_sample(code, ref seed, ref r, ref s);

                    Console.WriteLine("");
                    Console.WriteLine("  (R,S):"
                                      + "  " + r.ToString().PadLeft(12)
                                      + "  " + s.ToString().PadLeft(12) + "");

                    u_exact = Math.Pow(r, rexp[i]) * Math.Pow(s, sexp[i]);

                    dudr_exact = (double) (rexp[i])
                                 * Math.Pow(r, rexp[i] - 1) * Math.Pow(s, sexp[i]);

                    duds_exact = Math.Pow(r, rexp[i]) * (double) (sexp[i])
                                                      * Math.Pow(s, sexp[i] - 1);

                    interp(code, element_order, r, s, node_u, ref u, ref dudr, ref duds);

                    Console.WriteLine("  (U ,U* ,Error): "
                                      + "  " + u_exact.ToString().PadLeft(12)
                                      + "  " + u.ToString().PadLeft(12)
                                      + "  " + Math.Abs(u_exact - u).ToString().PadLeft(12) + "");

                    Console.WriteLine("  (Ur,Ur*,Error): "
                                      + "  " + dudr_exact.ToString().PadLeft(12)
                                      + "  " + dudr.ToString().PadLeft(12)
                                      + "  " + Math.Abs(dudr_exact - dudr).ToString().PadLeft(12) + "");

                    Console.WriteLine("  (Us,Us*,Error): "
                                      + "  " + duds_exact.ToString().PadLeft(12)
                                      + "  " + duds.ToString().PadLeft(12)
                                      + "  " + Math.Abs(duds_exact - duds).ToString().PadLeft(12) + "");
                }
            }
        }

    }
}