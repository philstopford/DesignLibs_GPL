using System;
using Burkardt.Sequence;
using Burkardt.Types;

namespace VanderCorputSeqTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VAN_DER_CORPUT_TEST tests the VAN_DER_CORPUT library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("VAN_DER_CORPUT_TEST");
        Console.WriteLine("  Test the VAN_DER_CORPUT library.");

        vdc_test();
        vdc_inverse_test();
        vdc_sequence_test();
        vdc_base_test();

        Console.WriteLine("");
        Console.WriteLine("VAN_DER_CORPUT_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");

    }

    private static void vdc_base_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VDC_BASE_TEST tests VDC_BASE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        double r2;
        double r3;
        double r5;

        Console.WriteLine("");
        Console.WriteLine("VDC_BASE_TEST");
        Console.WriteLine("  VDC_BASE returns the I-th element of a van der Corput");
        Console.WriteLine("  sequence in base B.");
        Console.WriteLine("");
        Console.WriteLine("    I          VDC_BASE(I,2)   VDC_BASE(I,3)   VDC_BASE(I,5)");
        Console.WriteLine("");

        for (i = -10; i <= 10; i++)
        {
            r2 = VanDerCorput.vdc_base(i, 2);
            r3 = VanDerCorput.vdc_base(i, 3);
            r5 = VanDerCorput.vdc_base(i, 5);
            Console.WriteLine("  " + i.ToString().PadLeft(3)
                                   + "        " + r2.ToString().PadLeft(14)
                                   + "  " + r3.ToString().PadLeft(14)
                                   + "  " + r5.ToString().PadLeft(14) + "");
        }

    }

    private static void vdc_inverse_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VDC_INVERSE_TEST tests VDC_INVERSE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int i2;
        double r;

        Console.WriteLine("");
        Console.WriteLine("VDC_INVERSE_TEST");
        Console.WriteLine("  VDC_INVERSE inverts an element of a van der Corput sequence.");
        Console.WriteLine("");
        Console.WriteLine("    I        R=VDC(I)  VDC_INVERSE(R)");
        Console.WriteLine("");

        for (i = -10; i <= 10; i++)
        {
            r = VanDerCorput.vdc(i);
            i2 = VanDerCorput.vdc_inverse(r);
            Console.WriteLine("  " + i.ToString().PadLeft(3)
                                   + "  " + r.ToString().PadLeft(14)
                                   + "  " + i2.ToString().PadLeft(3) + "");
        }

    }

    private static void vdc_sequence_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VDC_SEQUENCE_TEST tests VDC_SEQUENCE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i1;
        int i2;
        int n;
        double[] r;

        Console.WriteLine("");
        Console.WriteLine("VDC_SEQUENCE_TEST");
        Console.WriteLine("  VDC_SEQUENCE returns elements I1 through I2 of");
        Console.WriteLine("  a van der Corput sequence.");

        i1 = 7;
        i2 = 7;
        n = Math.Abs(i2 - i1) + 1;
        r = VanDerCorput.vdc_sequence(i1, i2);
        Console.WriteLine("");
        typeMethods.r8vec_transpose_print(n, r, "  R=VDC_SEQUENCE(  7,  7):");

        i1 = 0;
        i2 = 8;
        n = Math.Abs(i2 - i1) + 1;
        r = VanDerCorput.vdc_sequence(i1, i2);
        Console.WriteLine("");
        typeMethods.r8vec_transpose_print(n, r, "  R=VDC_SEQUENCE(  0,  8):");

        i1 = 8;
        i2 = 0;
        n = Math.Abs(i2 - i1) + 1;
        r = VanDerCorput.vdc_sequence(i1, i2);
        Console.WriteLine("");
        typeMethods.r8vec_transpose_print(n, r, "  R=VDC_SEQUENCE(  8,  0):");

        i1 = -3;
        i2 = +5;
        n = Math.Abs(i2 - i1) + 1;
        r = VanDerCorput.vdc_sequence(i1, i2);
        Console.WriteLine("");
        typeMethods.r8vec_transpose_print(n, r, "  R=VDC_SEQUENCE( -3,  5):");

        i1 = 100;
        i2 = 105;
        n = Math.Abs(i2 - i1) + 1;
        r = VanDerCorput.vdc_sequence(i1, i2);
        Console.WriteLine("");
        typeMethods.r8vec_transpose_print(n, r, "  R=VDC_SEQUENCE(100,105):");

    }

    private static void vdc_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VDC_TEST tests VDC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        double r;

        Console.WriteLine("");
        Console.WriteLine("VDC_TEST");
        Console.WriteLine("  VDC returns the I-th element of a van der Corput sequence.");
        Console.WriteLine("");
        Console.WriteLine("    I          VDC(I)");
        Console.WriteLine("");

        for (i = -10; i <= 10; i++)
        {
            r = VanDerCorput.vdc(i);
            Console.WriteLine("  " + i.ToString().PadLeft(3)
                                   + "  " + r.ToString().PadLeft(14) + "");
        }

    }
}