using System;
using Burkardt;
using Burkardt.Sequence;
using Burkardt.Types;

namespace HaltonTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for HALTON_TEST.
        //
        //  Discussion:
        //
        //    HALTON_TEST tests HALTON.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("HALTON_TEST:");
        Console.WriteLine("  Test the HALTON library");

        halton_test();
        halton_sequence_test();
        halton_inverse_test();
        halton_base_test();

        Console.WriteLine("");
        Console.WriteLine("HALTON_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void halton_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HALTON_TEST tests HALTON.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        int m;
        double[] r;

        Console.WriteLine("");
        Console.WriteLine("HALTON_TEST");
        Console.WriteLine("  HALTON returns the I-th element of an M-dimensional");
        Console.WriteLine("  Halton sequence.");
        Console.WriteLine("");
        Console.WriteLine("    I          HALTON(I)");

        for (m = 1; m <= 3; m++)
        {
            Console.WriteLine("");
            Console.WriteLine("  Use M = " + m + "");
            Console.WriteLine("");
            for (i = 0; i <= 10; i++)
            {
                r = Halton.halton(i, m);
                string cout = "  " + i.ToString().PadLeft(3);
                for (j = 0; j < m; j++)
                {
                    cout += "  " + r[j].ToString().PadLeft(14);
                }

                Console.WriteLine(cout);
                    
            }
        }
    }

    private static void halton_base_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HALTON_BASE_TEST tests HALTON_BASE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] b1 =  {
                2, 3, 5
            }
            ;
        int[] b2 =  {
                3, 10, 2
            }
            ;
        int i;
        int j;
        int m;
        double[] r;

        Console.WriteLine("");
        Console.WriteLine("HALTON_BASE_TEST");
        Console.WriteLine("  HALTON_BASE returns the I-th element of an M-dimensional");
        Console.WriteLine("  Halton sequence, using user-specified bases.");

        m = 3;
        Console.WriteLine("");
        Console.WriteLine("  M = " + m + "");
        string cout = "  B:";
        for (j = 0; j < m; j++)
        {
            cout += "  " + b1[j].ToString().PadLeft(14);
        }

        Console.WriteLine(cout);
        Console.WriteLine("");
        for (i = 0; i <= 10; i++)
        {
            r = Halton.halton_base(i, m, b1);
            cout = "  " + i.ToString().PadLeft(3);
            for (j = 0; j < m; j++)
            {
                cout += "  " + r[j].ToString().PadLeft(14);
            }

            Console.WriteLine(cout);
                
        }

        m = 3;
        Console.WriteLine("");
        Console.WriteLine("  M = " + m + "");
        cout = "  B:";
        for (j = 0; j < m; j++)
        {
            cout += "  " + b2[j].ToString().PadLeft(14);
        }

        Console.WriteLine(cout);
        Console.WriteLine("");
        for (i = 0; i <= 10; i++)
        {
            r = Halton.halton_base(i, m, b2);
            cout = "  " + i.ToString().PadLeft(3);
            for (j = 0; j < m; j++)
            {
                cout += "  " + r[j].ToString().PadLeft(14);
            }

            Console.WriteLine(cout);
                
        }
    }

    private static void halton_inverse_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HALTON_INVERSE_TEST tests HALTON_INVERSE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 August 2016
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
        double[] r;

        Console.WriteLine("");
        Console.WriteLine("HALTON_INVERSE_TEST");
        Console.WriteLine("  HALTON_INVERSE inverts an element of a Halton sequence.");
        Console.WriteLine("");
        Console.WriteLine("    I        R=HALTON(I,3)  HALTON_INVERSE(R,3)");
        Console.WriteLine("");

        m = 3;

        for (i = 0; i <= 10; i++)
        {
            r = Halton.halton(i, m);
            i2 = Halton.halton_inverse(r, m);
            string cout = "  " + i.ToString().PadLeft(3);
            for (j = 0; j < m; j++)
            {
                cout += "  " + r[j].ToString().PadLeft(14);
            }

            Console.WriteLine(cout + "  " + i2.ToString().PadLeft(3) + "");
                
        }
    }

    private static void halton_sequence_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HALTON_SEQUENCE_TEST tests HALTON_SEQUENCE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int m;
        double[] r;

        Console.WriteLine("");
        Console.WriteLine("HALTON_SEQUENCE_TEST");
        Console.WriteLine("  HALTON_SEQUENCE returns the elements I1 through I2");
        Console.WriteLine("  of an M-dimensional Halton sequence.");

        for (m = 1; m <= 3; m++)
        {
            Console.WriteLine("");
            Console.WriteLine("  HALTON_SEQUENCE(0,10," + m + ",R):");
            r = Halton.halton_sequence(0, 10, m);
            typeMethods.r8mat_print(m, 11, r, "  R:");
                
        }

        m = 3;
        Console.WriteLine("");
        Console.WriteLine("  HALTON_SEQUENCE(10,0," + m + ",R):");
        r = Halton.halton_sequence(10, 0, m);
        typeMethods.r8mat_print(m, 11, r, "  R:");
    }
}