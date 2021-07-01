using System;
using System.Numerics;
using TestValues;

namespace TestValuesTest
{
    public class AiryTest
    {
        public static void airy_ai_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    AIRY_AI_VALUES_TEST tests AIRY_AI_VALUES.
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
            double ai = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("AIRY_AI_VALUES_TEST:");
            Console.WriteLine("  AIRY_AI_VALUES stores values of ");
            Console.WriteLine("  the Airy functions Ai(X).");
            Console.WriteLine("");
            Console.WriteLine("                X                     Ai(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Airy.airy_ai_values(ref n_data, ref x, ref ai);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + ai.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void airy_ai_int_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    AIRY_AI_INT_VALUES_TEST tests AIRY_AI_INT_VALUES.
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
            double fx = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("AIRY_AI_INT_VALUES_TEST:");
            Console.WriteLine("  AIRY_AI_INT_VALUES stores values of ");
            Console.WriteLine("  the integral of the Airy Ai function.");
            Console.WriteLine("");
            Console.WriteLine("                X                     FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Airy.airy_ai_int_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void airy_ai_prime_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    AIRY_AI_PRIME_VALUES_TEST tests AIRY_AI_PRIME_VALUES.
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
            double aip = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("AIRY_AI_PRIME_VALUES_TEST:");
            Console.WriteLine("  AIRY_AI_PRIME_VALUES stores values of ");
            Console.WriteLine("  the derivative of the Airy function Ai'(X).");
            Console.WriteLine("");
            Console.WriteLine("                X                    Ai'");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Airy.airy_ai_prime_values(ref n_data, ref x, ref aip);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                     + x.ToString("0.################").PadLeft(24) + "  "
                     + aip.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void airy_bi_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    AIRY_BI_VALUES_TEST tests AIRY_BI_VALUES.
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
            double bi = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("AIRY_BI_VALUES_TEST:");
            Console.WriteLine("  AIRY_BI_VALUES stores values of ");
            Console.WriteLine("  the Airy function Bi.");
            Console.WriteLine("");
            Console.WriteLine("                X                     Bi");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Airy.airy_bi_values(ref n_data, ref x, ref bi);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + bi.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void airy_bi_int_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    AIRY_BI_INT_VALUES_TEST tests AIRY_BI_INT_VALUES.
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
            double fx = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("AIRY_BI_INT_VALUES_TEST:");
            Console.WriteLine("  AIRY_BI_INT_VALUES stores values of ");
            Console.WriteLine("  the integral of the Airy Bi function.");
            Console.WriteLine("");
            Console.WriteLine("                X                     FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Airy.airy_bi_int_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void airy_bi_prime_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    AIRY_BI_PRIME_VALUES_TEST tests AIRY_BI_PRIME_VALUES.
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
            double bip = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("AIRY_BI_PRIME_VALUES_TEST:");
            Console.WriteLine("  AIRY_BI_PRIME_VALUES stores values of ");
            Console.WriteLine("  the derivative of Airy function Bi'(X).");
            Console.WriteLine("");
            Console.WriteLine("                X                     Bi'");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Airy.airy_bi_prime_values(ref n_data, ref x, ref bip);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                     + x.ToString("0.################").PadLeft(24) + "  "
                     + bip.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void airy_cai_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    AIRY_CAI_VALUES_TEST tests AIRY_CAI_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    13 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Complex cai = new Complex();
            int n_data;
            Complex x = new Complex();
            Console.WriteLine("");
            Console.WriteLine("AIRY_CAI_VALUES_TEST:");
            Console.WriteLine("  AIRY_CAI_VALUES stores values of ");
            Console.WriteLine("  the Airy functions Ai(X) for complex argument.");
            Console.WriteLine("");
            Console.WriteLine("                X                     Ai");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Airy.airy_cai_values(ref n_data, ref x, ref cai);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                     + (x.Real).ToString("0.######").PadLeft(14) + "  "
                     + (x.Imaginary).ToString("0.######").PadLeft(14) + "  "
                     + (cai.Real).ToString("0.################").PadLeft(24) + "  "
                     + (cai.Imaginary).ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void airy_cbi_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    AIRY_CBI_VALUES_TEST tests AIRY_CBI_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    13 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Complex cbi = new Complex();
            int n_data;
            Complex x = new Complex();
            Console.WriteLine("");
            Console.WriteLine("AIRY_CBI_VALUES_TEST:");
            Console.WriteLine("  AIRY_CBI_VALUES stores values of ");
            Console.WriteLine("  the Airy functions Bi(X) for complex argument.");
            Console.WriteLine("");
            Console.WriteLine("                X                     Bi");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Airy.airy_cbi_values(ref n_data, ref x, ref cbi);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + (x.Real).ToString("0.######").PadLeft(14) + "  "
                                  + (x.Imaginary).ToString("0.######").PadLeft(14) + "  "
                                  + (cbi.Real).ToString("0.################").PadLeft(24) + "  "
                                  + (cbi.Imaginary).ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void airy_gi_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    AIRY_GI_VALUES_TEST tests AIRY_GI_VALUES.
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
            double bip = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("AIRY_GI_VALUES_TEST:");
            Console.WriteLine("  AIRY_GI_VALUES stores values of ");
            Console.WriteLine("  the modified Airy function Gi(X).");
            Console.WriteLine("");
            Console.WriteLine("                X                     Gi");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Airy.airy_gi_values(ref n_data, ref x, ref bip);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + bip.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void airy_hi_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    AIRY_HI_VALUES_TEST tests AIRY_HI_VALUES.
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
            double bip = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("AIRY_HI_VALUES_TEST:");
            Console.WriteLine("  AIRY_HI_VALUES stores values of ");
            Console.WriteLine("  the modified Airy function Hi(X).");
            Console.WriteLine("");
            Console.WriteLine("                X                     Hi");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Airy.airy_hi_values(ref n_data, ref x, ref bip);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + bip.ToString("0.################").PadLeft(24) + "");
            }
        }

    }
}