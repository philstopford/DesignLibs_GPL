using System;
using TestValues;

namespace TestValuesTest
{
    public static class FrobeniusTest
    {
        public static void frobenius_number_data_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FROBENIUS_NUMBER_DATA_VALUES_TEST tests FROBENIUS_NUMBER_DATA_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    27 January 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int[] c;
            int f = 0;
            int i;
            int n_data;
            int order = 0;
            Console.WriteLine("");
            Console.WriteLine("FROBENIUS_NUMBER_DATA_VALUES_TEST:");
            Console.WriteLine("  FROBENIUS_NUMBER_DATA_VALUES returns the corresponding");
            Console.WriteLine("  coin denominations.");
            n_data = 0;
            for (;;)
            {
                Frobenius.frobenius_number_order_values(ref n_data, ref order);
                if (n_data == 0)
                {
                    break;
                }

                c = new int[order];
                Frobenius.frobenius_number_data_values(ref n_data, order, ref c, ref f);
                Console.WriteLine("");
                Console.WriteLine("  Order = " + order + "");
                string cout = "";
                for (i = 0; i < order; i++)
                {
                    cout += "  " + c[i].ToString().PadLeft(8);
                }

                Console.WriteLine(cout);
                Console.WriteLine("  Frobenius number = " + f + "");
            }
        }

        public static void frobenius_number_order_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FROBENIUS_NUMBER_ORDER_VALUES tests FROBENIUS_NUMBER_ORDER_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    27 January 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n_data;
            int order = 0;
            Console.WriteLine("");
            Console.WriteLine("FROBENIUS_NUMBER_ORDER_VALUES_TEST:");
            Console.WriteLine("  FROBENIUS_NUMBER_ORDER_VALUES returns the order for");
            Console.WriteLine("  a Frobenius problem;");
            Console.WriteLine("");
            Console.WriteLine("   Problem   ORDER");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Frobenius.frobenius_number_order_values(ref n_data, ref order);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("");
                Console.WriteLine("  " + n_data.ToString().PadLeft(4)
                                       + "  " + order.ToString().PadLeft(4) + "");
            }
        }

        public static void frobenius_number_order2_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FROBENIUS_NUMBER_ORDER2_VALUES_TEST tests FROBENIUS_NUMBER_ORDER2_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 January 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int c1 = 0;
            int c2 = 0;
            int f = 0;
            int n_data;
            Console.WriteLine("");
            Console.WriteLine("FROBENIUS_NUMBER_ORDER2_VALUES_TEST:");
            Console.WriteLine("  FROBENIUS_NUMBER_ORDER2_VALUES returns values of ");
            Console.WriteLine("  the Frobenius number of order 2.");
            Console.WriteLine("");
            Console.WriteLine("         C1        C2          F(C1,C2)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Frobenius.frobenius_number_order2_values(ref n_data, ref c1, ref c2, ref f);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + c1.ToString().PadLeft(8)
                                       + "  " + c2.ToString().PadLeft(8)
                                       + "  " + f.ToString().PadLeft(8) + "");
            }
        }

    }
}