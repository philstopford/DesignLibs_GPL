using System;
using System.Globalization;
using Burkardt.Sequence;
using Burkardt.Types;

namespace HammersleyTest;

internal static class Program
{
    private static void Main()
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

    private static void hammersley_test()

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
        int m;

        Console.WriteLine("");
        Console.WriteLine("HAMMERSLEY_TEST");
        Console.WriteLine("  HAMMERSLEY returns the I-th element of an M-dimensional");
        Console.WriteLine("  Hammersley sequence.");
        Console.WriteLine("");
        Console.WriteLine("    I          HAMMERSLEY(I)");

        int n = 16;

        for (m = 1; m <= 3; m++)
        {
            Console.WriteLine("");
            Console.WriteLine("  Use M = " + m + "");
            Console.WriteLine("      N = " + n + "");
            Console.WriteLine("");
            int i;
            for (i = 0; i <= 10; i++)
            {
                string cout = "  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(3);
                double[] r = Hammersley.hammersley(i, m, n);
                int j;
                for (j = 0; j < m; j++)
                {
                    cout += "  " + r[j].ToString(CultureInfo.InvariantCulture).PadLeft(14);
                }

                Console.WriteLine(cout);
            }
        }
    }

    private static void hammersley_inverse_test()

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

        Console.WriteLine("");
        Console.WriteLine("HAMMERSLEY_INVERSE_TEST");
        Console.WriteLine("  HAMMERSLEY_INVERSE inverts an element of a Hammersley sequence.");
        Console.WriteLine("");
        Console.WriteLine("    I        R=HAMMERSLEY(I,3)  HAMMERSLEY_INVERSE(R,3)");
        Console.WriteLine("");

        int m = 3;
        int n = 16;

        for (i = 0; i <= 10; i++)
        {
            string cout = "  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(3);
            double[] r = Hammersley.hammersley(i, m, n);
            int j;
            for (j = 0; j < m; j++)
            {
                cout += "  " + r[j].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            }

            int i2 = Hammersley.hammersley_inverse(r, m, n);
            Console.WriteLine(cout + "  " + i2.ToString(CultureInfo.InvariantCulture).PadLeft(3) + "");
        }
    }

    private static void hammersley_sequence_test()

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
        double[] r;

        Console.WriteLine("");
        Console.WriteLine("HAMMERSLEY_SEQUENCE_TEST");
        Console.WriteLine("  HAMMERSLEY_SEQUENCE returns the elements I1 through I2");
        Console.WriteLine("  of an M-dimensional Hammersley sequence.");

        int n = 16;

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