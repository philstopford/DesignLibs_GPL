using System;
using Burkardt;
using Burkardt.Sequence;
using Burkardt.Types;

namespace HammersleyTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for HAMMERSLEY_TEST.
            //
            //  Discussion:
            //
            //    HAMMERSLEY_TEST tests the HAMMERSLEY library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 August 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("HAMMERSLEY_TEST:");
            Console.WriteLine("  Test the HAMMERSLEY library.");

            hammersley_test();
            hammersley_inverse_test();
            hammersley_sequence_test();

            Console.WriteLine("");
            Console.WriteLine("HAMMERSLEY_TEST:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void hammersley_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HAMMERSLEY_TEST tests HAMMERSLEY.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 August 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            int j;
            int m;
            int n;
            double[] r;

            Console.WriteLine("");
            Console.WriteLine("HAMMERSLEY_TEST");
            Console.WriteLine("  HAMMERSLEY returns the I-th element of an M-dimensional");
            Console.WriteLine("  Hammersley sequence.");
            Console.WriteLine("");
            Console.WriteLine("    I          HAMMERSLEY(I)");

            n = 16;

            for (m = 1; m <= 3; m++)
            {
                Console.WriteLine("");
                Console.WriteLine("  Use M = " + m + "");
                Console.WriteLine("      N = " + n + "");
                Console.WriteLine("");
                for (i = 0; i <= 10; i++)
                {
                    string cout = "  " + i.ToString().PadLeft(3);
                    r = Hammersley.hammersley(i, m, n);
                    for (j = 0; j < m; j++)
                    {
                        cout += "  " + r[j].ToString().PadLeft(14);;
                    }

                    Console.WriteLine(cout);
                }
            }
        }

        static void hammersley_inverse_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HAMMERSLEY_INVERSE_TEST tests HAMMERSLEY_INVERSE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 August 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            int i2;
            int j;
            int m;
            int n;
            double[] r;

            Console.WriteLine("");
            Console.WriteLine("HAMMERSLEY_INVERSE_TEST");
            Console.WriteLine("  HAMMERSLEY_INVERSE inverts an element of a Hammersley sequence.");
            Console.WriteLine("");
            Console.WriteLine("    I        R=HAMMERSLEY(I,3)  HAMMERSLEY_INVERSE(R,3)");
            Console.WriteLine("");

            m = 3;
            n = 16;

            for (i = 0; i <= 10; i++)
            {
                string cout = "  " + i.ToString().PadLeft(3);
                r = Hammersley.hammersley(i, m, n);
                for (j = 0; j < m; j++)
                {
                    cout += "  " + r[j].ToString().PadLeft(14);;
                }

                i2 = Hammersley.hammersley_inverse(r, m, n);
                Console.WriteLine(cout + "  " + i2.ToString().PadLeft(3) + "");
            }
        }

        static void hammersley_sequence_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HAMMERSLEY_SEQUENCE_TEST tests HAMMERSLEY_SEQUENCE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 August 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int m;
            int n;
            double[] r;

            Console.WriteLine("");
            Console.WriteLine("HAMMERSLEY_SEQUENCE_TEST");
            Console.WriteLine("  HAMMERSLEY_SEQUENCE returns the elements I1 through I2");
            Console.WriteLine("  of an M-dimensional Hammersley sequence.");

            n = 16;

            for (m = 1; m <= 3; m++)
            {
                Console.WriteLine("");
                Console.WriteLine("  HAMMERSLEY_SEQUENCE(0,10," + m + ",N,R):");
                r = Hammersley.hammersley_sequence(0, 10, m, n);
                typeMethods.r8mat_print(m, 11, r, "  R:");
            }

            m = 3;
            n = 16;
            Console.WriteLine("");
            string cout = "  HAMMERSLEY_SEQUENCE(10,0," + m + ",N,R):";
            r = Hammersley.hammersley_sequence(10, 0, m, n);
            typeMethods.r8mat_print(m, 11, r, cout + "  R:");
        }
    }
}