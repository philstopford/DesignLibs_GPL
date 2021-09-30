using System;
using Burkardt.Values;

namespace TestValuesTest
{
    public class WignerTest
    {
        public static void nine_j_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NINE_J_VALUES_TEST demonstrates NINE_J_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fx = 0;
            double j1 = 0;
            double j2 = 0;
            double j3 = 0;
            double j4 = 0;
            double j5 = 0;
            double j6 = 0;
            double j7 = 0;
            double j8 = 0;
            double j9 = 0;
            int n_data;
            Console.WriteLine("");
            Console.WriteLine("NINE_J_VALUES_TEST:");
            Console.WriteLine("  NINE_J_VALUES returns values of");
            Console.WriteLine("  the Wigner 9J coefficient.");
            Console.WriteLine("");
            Console.WriteLine("      J1      J2      J3      J4      J5      J6"
                              + "      J7      J8      J9        NINE_J");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Wigner.nine_j_values(ref n_data, ref j1, ref j2, ref j3, ref j4, ref j5, ref j6, ref j7, ref j8, ref j9,
                    ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + j1.ToString().PadLeft(6)
                                       + "  " + j2.ToString().PadLeft(6)
                                       + "  " + j3.ToString().PadLeft(6)
                                       + "  " + j4.ToString().PadLeft(6)
                                       + "  " + j5.ToString().PadLeft(6)
                                       + "  " + j6.ToString().PadLeft(6)
                                       + "  " + j7.ToString().PadLeft(6)
                                       + "  " + j8.ToString().PadLeft(6)
                                       + "  " + j9.ToString().PadLeft(6)
                                       + "  " + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void six_j_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIX_J_VALUES_TEST tests SIX_J_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fx = 0;
            double j1 = 0;
            double j2 = 0;
            double j3 = 0;
            double j4 = 0;
            double j5 = 0;
            double j6 = 0;
            int n_data;
            Console.WriteLine("");
            Console.WriteLine("SIX_J_VALUES_TEST:");
            Console.WriteLine("  SIX_J_VALUES returns values of ");
            Console.WriteLine("  the Wigner 6J coefficient.");
            Console.WriteLine("");
            Console.WriteLine(
                "      J1      J2      J3      J4      J5      J6        SIX_J");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Wigner.six_j_values(ref n_data, ref j1, ref j2, ref j3, ref j4, ref j5, ref j6, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + j1.ToString().PadLeft(6)
                                       + "  " + j2.ToString().PadLeft(6)
                                       + "  " + j3.ToString().PadLeft(6)
                                       + "  " + j4.ToString().PadLeft(6)
                                       + "  " + j5.ToString().PadLeft(6)
                                       + "  " + j6.ToString().PadLeft(6)
                                       + "  " + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void three_j_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    THREE_J_VALUES_TEST tests THREE_J_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fx = 0;
            double j1 = 0;
            double j2 = 0;
            double j3 = 0;
            double m1 = 0;
            double m2 = 0;
            double m3 = 0;
            int n_data;
            Console.WriteLine("");
            Console.WriteLine("THREE_J_VALUES_TEST:");
            Console.WriteLine("  THREE_J_VALUES returns values of");
            Console.WriteLine("  the Wigner 3J coefficient.");
            Console.WriteLine("");
            Console.WriteLine("      J1      J2      J3      M1      M2      M3        THREE_J");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Wigner.three_j_values(ref n_data, ref j1, ref j2, ref j3, ref m1, ref m2, ref m3, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + j1.ToString().PadLeft(6)
                                       + "  " + j2.ToString().PadLeft(6)
                                       + "  " + j3.ToString().PadLeft(6)
                                       + "  " + m1.ToString().PadLeft(6)
                                       + "  " + m2.ToString().PadLeft(6)
                                       + "  " + m3.ToString().PadLeft(6)
                                       + "  " + fx.ToString("0.################").PadLeft(24) + "");
            }
        }
    }
}