using System;
using Burkardt.PolynomialNS;

namespace LegendreProductPolynomialTest
{
    public static class polyTest
    {
        public static void polynomial_compress_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYNOMIAL_COMPRESS_TEST tests POLYNOMIAL_COMPRESS.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    27 October 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] c =  {
                7.0, -5.0, 5.0, 9.0, 11.0, 3.0, 6.0, 0.0, -13.0, 1.0E-20
            }
            ;
            double[] c2 = new double[10];
            int m = 3;
            int[] e =  {
                1, 2, 2, 4, 5, 5, 5, 12, 33, 35
            }
            ;
            int[] e2 = new int[10];
            int o = 10;
            int o2 = 0;
            string title;

            Console.WriteLine("");
            Console.WriteLine("POLYNOMIAL_COMPRESS_TEST");
            Console.WriteLine("  POLYNOMIAL_COMPRESS compresses a polynomial.");

            Console.WriteLine("");
            title = "  Uncompressed P(X) = ";
            Polynomial.polynomial_print(m, o, c, e, title);

            Polynomial.polynomial_compress(o, c, e, ref o2, ref c2, ref e2);

            Console.WriteLine("");
            title = "  Compressed P(X) = ";
            Polynomial.polynomial_print(m, o2, c2, e2, title);
        }

        public static void polynomial_print_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYNOMIAL_PRINT_TEST tests POLYNOMIAL_PRINT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    27 October 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] c =  {
                7.0, -5.0, 9.0, 11.0, 0.0, -13.0
            }
            ;
            int m = 3;
            int[] e =  {
                1, 2, 4, 5, 12, 33
            }
            ;
            int o = 6;
            string title = "  P1(X) =";

            Console.WriteLine("");
            Console.WriteLine("POLYNOMIAL_PRINT_TEST");
            Console.WriteLine("  POLYNOMIAL_PRINT prints a polynomial.");

            Console.WriteLine("");
            Polynomial.polynomial_print(m, o, c, e, title);
        }

        public static void polynomial_sort_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYNOMIAL_SORT_TEST tests POLYNOMIAL_SORT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    27 October 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] c =  {
                0.0, 9.0, -5.0, -13.0, 7.0, 11.0
            }
            ;
            int m = 3;
            int[] e =  {
                12, 4, 2, 33, 1, 5
            }
            ;
            int o = 6;
            string title;

            Console.WriteLine("");
            Console.WriteLine("POLYNOMIAL_SORT_TEST");
            Console.WriteLine("  POLYNOMIAL_SORT sorts a polynomial by exponent index.");

            Console.WriteLine("");
            title = "  Unsorted polynomial:";
            Polynomial.polynomial_print(m, o, c, e, title);

            Polynomial.polynomial_sort(o, ref c, ref e);

            Console.WriteLine("");
            title = "  Sorted polynomial:";
            Polynomial.polynomial_print(m, o, c, e, title);
        }

        public static void polynomial_value_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYNOMIAL_VALUE_TEST tests POLYNOMIAL_VALUE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    27 October 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] c =  {
                7.0, -5.0, 9.0, 11.0, 0.0, -13.0
            }
            ;
            int m = 3;
            int[] e =  {
                1, 2, 4, 5, 12, 33
            }
            ;
            int j;
            int nx = 2;
            int o = 6;
            double[] p;
            string title = "  P(X) =";
            double[] x =  {
                1.0, 2.0, 3.0,
                -2.0, 4.0, 1.0
            }
            ;

            Console.WriteLine("");
            Console.WriteLine("POLYNOMIAL_VALUE_TEST");
            Console.WriteLine("  POLYNOMIAL_VALUE evaluates a polynomial.");

            Console.WriteLine("");
            Polynomial.polynomial_print(m, o, c, e, title);

            p = Polynomial.polynomial_value(m, o, c, e, nx, x);

            Console.WriteLine("");
            for (j = 0; j < nx; j++)
            {
                Console.WriteLine("  P(" + x[0 + j * m]
                    + "," + x[1 + j * m]
                    + "," + x[2 + j * m]
                    + ") = " + p[j] + "");
            }
        }
    }
}