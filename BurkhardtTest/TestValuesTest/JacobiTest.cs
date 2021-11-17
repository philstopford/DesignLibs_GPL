    using System;
    using Burkardt.Values;

    namespace TestValuesTest; 

    public class JacobiTest
    {
        public static void jacobi_cn_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    jacobi_cn_values_test tests jacobi_cn_values().
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 November 2020
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a = 0;
            double fx = 0;
            double k = 0;
            double m = 0;
            int n_data;
            double u = 0;
            Console.WriteLine("");
            Console.WriteLine("jacobi_cn_values_test:");
            Console.WriteLine("  jacobi_cn_values() returns values of ");
            Console.WriteLine("  the Jacobi elliptic CN function.");
            Console.WriteLine("");
            Console.WriteLine("      U         M       CN(U,M)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Jacobi.jacobi_cn_values(ref n_data, ref u, ref a, ref k, ref m, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + u.ToString().PadLeft(10)
                                       + "  " + m.ToString().PadLeft(10)
                                       + "  " + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void jacobi_dn_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    jacobi_dn_values_test tests jacobi_dn_values().
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 November 2020
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a = 0;
            double fx = 0;
            double k = 0;
            double m = 0;
            int n_data;
            double u = 0;
            Console.WriteLine("");
            Console.WriteLine("jacobi_dn_values_test:");
            Console.WriteLine("  jacobi_dn_values() returns values of ");
            Console.WriteLine("  the Jacobi elliptic DN function.");
            Console.WriteLine("");
            Console.WriteLine("      U         M       DN(U,M)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Jacobi.jacobi_dn_values(ref n_data, ref u, ref a, ref k, ref m, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + u.ToString().PadLeft(10)
                                       + "  " + m.ToString().PadLeft(10)
                                       + "  " + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void jacobi_poly_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    JACOBI_POLY_VALUES_TEST tests JACOBI_POLY_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 April 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a = 0;
            double b = 0;
            double fx = 0;
            int n = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("JACOBI_POLY_VALUES_TEST:");
            Console.WriteLine("  JACOBI_POLY_VALUES returns values of");
            Console.WriteLine("  the Jacobi polynomial.");
            Console.WriteLine("");
            Console.WriteLine("       N         A         B      X       J(N,A,B)(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Jacobi.jacobi_poly_values(ref n_data, ref n, ref a, ref b, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + n.ToString().PadLeft(6) + "  "
                                  + a.ToString().PadLeft(8) + "  "
                                  + b.ToString().PadLeft(8) + "  "
                                  + x.ToString().PadLeft(10) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void jacobi_sn_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    jacobi_sn_values_test tests jacobi_sn_values().
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 November 2020
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a = 0;
            double fx = 0;
            double k = 0;
            double m = 0;
            int n_data;
            double u = 0;
            Console.WriteLine("");
            Console.WriteLine("jacobi_sn_values_test:");
            Console.WriteLine("  jacobi_sn_values() returns values of ");
            Console.WriteLine("  the Jacobi elliptic SN function.");
            Console.WriteLine("");
            Console.WriteLine("      U         M       SN(U,M)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Jacobi.jacobi_sn_values(ref n_data, ref u, ref a, ref k, ref m, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + u.ToString().PadLeft(10)
                                       + "  " + m.ToString().PadLeft(10)
                                       + "  " + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

    }