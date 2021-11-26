﻿using System;
using System.Globalization;
using Burkardt.MatrixNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace MatrixConditionTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for CONDITION_TEST.
        //
        //  Discussion:
        //
        //    CONDITION_TEST tests the CONDITION library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("CONDITION_TEST");

        Console.WriteLine("  Test the CONDITION library.");
        Console.WriteLine("  This test also requires the R8LIB library.");

        condition_linpack_test();
        condition_sample1_test();
        condition_hager_test();
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("CONDITION_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void condition_linpack_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CONDITION_LINPACK_TEST tests CONDITION_LINPACK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("CONDITION_LINPACK_TEST");
        Console.WriteLine("  CONDITION_LINPACK estimates the L1 condition number");
        Console.WriteLine("  For a matrix in general storage.");
        Console.WriteLine("");
        Console.WriteLine("  Matrix               Order   Condition         Linpack");
        Console.WriteLine("");
        //
        //  Combinatorial matrix.
        //
        string name = "Combinatorial";
        int n = 4;
        double alpha = 2.0;
        const double beta = 3.0;
        double[] a = Matrix.combin(alpha, beta, n);
        double[] a_inverse = Matrix.combin_inverse(alpha, beta, n);
        double a_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a);
        double a_inverse_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a_inverse);
        double cond_l1 = a_norm_l1 * a_inverse_norm_l1;
        double cond = Matrix.condition_linpack(n, a);
        Console.WriteLine("  " + name.ToString(CultureInfo.InvariantCulture).PadLeft(20)
                               + "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                               + "  " + cond_l1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + cond.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        //
        //  CONEX1
        //
        name = "CONEX1";
        n = 4;
        alpha = 100.0;
        a = Matrix.conex1(alpha);
        a_inverse = Matrix.conex1_inverse(alpha);
        a_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a);
        a_inverse_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a_inverse);
        cond_l1 = a_norm_l1 * a_inverse_norm_l1;
        cond = Matrix.condition_linpack(n, a);
        Console.WriteLine("  " + name.ToString(CultureInfo.InvariantCulture).PadLeft(20)
                               + "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                               + "  " + cond_l1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + cond.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        //
        //  CONEX2
        //
        name = "CONEX2";
        n = 3;
        alpha = 100.0;
        a = Matrix.conex2(alpha);
        a_inverse = Matrix.conex2_inverse(alpha);
        a_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a);
        a_inverse_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a_inverse);
        cond_l1 = a_norm_l1 * a_inverse_norm_l1;
        cond = Matrix.condition_linpack(n, a);
        Console.WriteLine("  " + name.ToString(CultureInfo.InvariantCulture).PadLeft(20)
                               + "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                               + "  " + cond_l1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + cond.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        //
        //  CONEX3
        //
        name = "CONEX3";
        n = 5;
        a = Matrix.conex3(n);
        a_inverse = Matrix.conex3_inverse(n);
        a_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a);
        a_inverse_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a_inverse);
        cond_l1 = a_norm_l1 * a_inverse_norm_l1;
        cond = Matrix.condition_linpack(n, a);
        Console.WriteLine("  " + name.ToString(CultureInfo.InvariantCulture).PadLeft(20)
                               + "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                               + "  " + cond_l1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + cond.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        //
        //  CONEX4
        //
        name = "CONEX4";
        n = 4;
        a = Matrix.conex4();
        a_inverse = Matrix.conex4_inverse();
        a_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a);
        a_inverse_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a_inverse);
        cond_l1 = a_norm_l1 * a_inverse_norm_l1;
        cond = Matrix.condition_linpack(n, a);
        Console.WriteLine("  " + name.ToString(CultureInfo.InvariantCulture).PadLeft(20)
                               + "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                               + "  " + cond_l1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + cond.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        //
        //  KAHAN
        //
        name = "KAHAN";
        n = 4;
        alpha = 0.25;
        a = Matrix.kahan(alpha, n, n);
        a_inverse = Matrix.kahan_inverse(alpha, n);
        a_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a);
        a_inverse_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a_inverse);
        cond_l1 = a_norm_l1 * a_inverse_norm_l1;
        cond = Matrix.condition_linpack(n, a);
        Console.WriteLine("  " + name.ToString(CultureInfo.InvariantCulture).PadLeft(20)
                               + "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                               + "  " + cond_l1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + cond.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        //
        //  Random
        //
        int seed = 123456789;

        for (i = 1; i <= 5; i++)
        {
            name = "RANDOM";
            n = 4;
            a = UniformRNG.r8mat_uniform_01_new(n, n, ref seed);
            double[] a_lu = typeMethods.r8mat_copy_new(n, n, a);
            int[] pivot = new int[n];
            typeMethods.r8ge_fa(n, ref a_lu, ref pivot);
            a_inverse = typeMethods.r8ge_inverse(n, a_lu, pivot);
            a_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a);
            a_inverse_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a_inverse);
            cond_l1 = a_norm_l1 * a_inverse_norm_l1;
            cond = Matrix.condition_linpack(n, a);
            Console.WriteLine("  " + name.ToString(CultureInfo.InvariantCulture).PadLeft(20)
                                   + "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + cond_l1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + cond.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void condition_sample1_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CONDITION_SAMPLE1_TEST tests CONDITION_SAMPLE1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 August 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double cond;
        int i;
        int j;
        int m;
        int[] m_test =  {
                10, 1000, 100000
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("CONDITION_SAMPLE1_TEST");
        Console.WriteLine("  CONDITION_SAMPLE1 estimates the L1 condition number using sampling");
        Console.WriteLine("  for a matrix in general storage.");
        Console.WriteLine("");
        Console.WriteLine("  Matrix                 Samples Order   Condition        Estimate");
        //
        //  Combinatorial matrix.
        //
        string name = "Combinatorial";
        int n = 4;
        double alpha = 2.0;
        const double beta = 3.0;
        double[] a = Matrix.combin(alpha, beta, n);
        double[] a_inverse = Matrix.combin_inverse(alpha, beta, n);
        double a_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a);
        double a_inverse_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a_inverse);
        double cond_l1 = a_norm_l1 * a_inverse_norm_l1;
        Console.WriteLine("");
        for (i = 0; i < 3; i++)
        {
            m = m_test[i];
            cond = Matrix.condition_sample1(n, a, m);
            Console.WriteLine("  " + name.ToString(CultureInfo.InvariantCulture).PadLeft(20)
                                   + "  " + m.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + cond_l1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + cond.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        //
        //  CONEX1
        //
        name = "CONEX1";
        n = 4;
        alpha = 100.0;
        a = Matrix.conex1(alpha);
        a_inverse = Matrix.conex1_inverse(alpha);
        a_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a);
        a_inverse_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a_inverse);
        cond_l1 = a_norm_l1 * a_inverse_norm_l1;
        Console.WriteLine("");
        for (i = 0; i < 3; i++)
        {
            m = m_test[i];
            cond = Matrix.condition_sample1(n, a, m);
            Console.WriteLine("  " + name.ToString(CultureInfo.InvariantCulture).PadLeft(20)
                                   + "  " + m.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + cond_l1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + cond.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        //
        //  CONEX2
        //
        name = "CONEX2";
        n = 3;
        alpha = 100.0;
        a = Matrix.conex2(alpha);
        a_inverse = Matrix.conex2_inverse(alpha);
        a_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a);
        a_inverse_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a_inverse);
        cond_l1 = a_norm_l1 * a_inverse_norm_l1;
        Console.WriteLine("");
        for (i = 0; i < 3; i++)
        {
            m = m_test[i];
            cond = Matrix.condition_sample1(n, a, m);
            Console.WriteLine("  " + name.ToString(CultureInfo.InvariantCulture).PadLeft(20)
                                   + "  " + m.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + cond_l1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + cond.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        //
        //  CONEX3
        //
        name = "CONEX3";
        n = 5;
        a = Matrix.conex3(n);
        a_inverse = Matrix.conex3_inverse(n);
        a_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a);
        a_inverse_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a_inverse);
        cond_l1 = a_norm_l1 * a_inverse_norm_l1;
        Console.WriteLine("");
        for (i = 0; i < 3; i++)
        {
            m = m_test[i];
            cond = Matrix.condition_sample1(n, a, m);
            Console.WriteLine("  " + name.ToString(CultureInfo.InvariantCulture).PadLeft(20)
                                   + "  " + m.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + cond_l1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + cond.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        //
        //  CONEX4
        //
        name = "CONEX4";
        n = 4;
        a = Matrix.conex4();
        a_inverse = Matrix.conex4_inverse();
        a_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a);
        a_inverse_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a_inverse);
        cond_l1 = a_norm_l1 * a_inverse_norm_l1;
        Console.WriteLine("");
        for (i = 0; i < 3; i++)
        {
            m = m_test[i];
            cond = Matrix.condition_sample1(n, a, m);
            Console.WriteLine("  " + name.ToString(CultureInfo.InvariantCulture).PadLeft(20)
                                   + "  " + m.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + cond_l1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + cond.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        //
        //  KAHAN
        //
        name = "KAHAN";
        n = 4;
        alpha = 0.25;
        a = Matrix.kahan(alpha, n, n);
        a_inverse = Matrix.kahan_inverse(alpha, n);
        a_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a);
        a_inverse_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a_inverse);
        cond_l1 = a_norm_l1 * a_inverse_norm_l1;
        Console.WriteLine("");
        for (i = 0; i < 3; i++)
        {
            m = m_test[i];
            cond = Matrix.condition_sample1(n, a, m);
            Console.WriteLine("  " + name.ToString(CultureInfo.InvariantCulture).PadLeft(20)
                                   + "  " + m.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + cond_l1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + cond.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        //
        //  Random
        //
        int seed = 123456789;

        for (j = 1; j <= 5; j++)
        {
            name = "RANDOM";
            n = 4;
            a = UniformRNG.r8mat_uniform_01_new(n, n, ref seed);
            double[] a_lu = typeMethods.r8mat_copy_new(n, n, a);
            int[] pivot = new int[n];
            typeMethods.r8ge_fa(n, ref a_lu, ref pivot);
            a_inverse = typeMethods.r8ge_inverse(n, a_lu, pivot);
            a_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a);
            a_inverse_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a_inverse);
            cond_l1 = a_norm_l1 * a_inverse_norm_l1;
            Console.WriteLine("");
            for (i = 0; i < 3; i++)
            {
                m = m_test[i];
                cond = Matrix.condition_sample1(n, a, m);
                Console.WriteLine("  " + name.ToString(CultureInfo.InvariantCulture).PadLeft(20)
                                       + "  " + m.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                       + "  " + cond_l1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + cond.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }
        }
    }

    private static void condition_hager_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CONDITION_HAGER_TEST tests CONDITION_HAGER.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("CONDITION_HAGER_TEST");
        Console.WriteLine("  CONDITION_HAGER estimates the L1 condition number");
        Console.WriteLine("  for a matrix in general storage.");
        Console.WriteLine("");
        Console.WriteLine("  Matrix               Order   Condition         Hager");
        Console.WriteLine("");
        //
        //  Combinatorial matrix.
        //
        string name = "Combinatorial";
        int n = 4;
        double alpha = 2.0;
        const double beta = 3.0;
        double[] a = Matrix.combin(alpha, beta, n);
        double[] a_inverse = Matrix.combin_inverse(alpha, beta, n);
        double a_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a);
        double a_inverse_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a_inverse);
        double cond_l1 = a_norm_l1 * a_inverse_norm_l1;
        double cond = Matrix.condition_hager(n, a);
        Console.WriteLine("  " + name.ToString(CultureInfo.InvariantCulture).PadLeft(20)
                               + "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                               + "  " + cond_l1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + cond.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        //
        //  CONEX1
        //
        name = "CONEX1";
        n = 4;
        alpha = 100.0;
        a = Matrix.conex1(alpha);
        a_inverse = Matrix.conex1_inverse(alpha);
        a_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a);
        a_inverse_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a_inverse);
        cond_l1 = a_norm_l1 * a_inverse_norm_l1;
        cond = Matrix.condition_hager(n, a);
        Console.WriteLine("  " + name.ToString(CultureInfo.InvariantCulture).PadLeft(20)
                               + "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                               + "  " + cond_l1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + cond.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        //
        //  CONEX2
        //
        name = "CONEX2";
        n = 3;
        alpha = 100.0;
        a = Matrix.conex2(alpha);
        a_inverse = Matrix.conex2_inverse(alpha);
        a_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a);
        a_inverse_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a_inverse);
        cond_l1 = a_norm_l1 * a_inverse_norm_l1;
        cond = Matrix.condition_hager(n, a);
        Console.WriteLine("  " + name.ToString(CultureInfo.InvariantCulture).PadLeft(20)
                               + "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                               + "  " + cond_l1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + cond.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        //
        //  CONEX3
        //
        name = "CONEX3";
        n = 5;
        a = Matrix.conex3(n);
        a_inverse = Matrix.conex3_inverse(n);
        a_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a);
        a_inverse_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a_inverse);
        cond_l1 = a_norm_l1 * a_inverse_norm_l1;
        cond = Matrix.condition_hager(n, a);
        Console.WriteLine("  " + name.ToString(CultureInfo.InvariantCulture).PadLeft(20)
                               + "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                               + "  " + cond_l1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + cond.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        //
        //  CONEX4
        //
        name = "CONEX4";
        n = 4;
        a = Matrix.conex4();
        a_inverse = Matrix.conex4_inverse();
        a_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a);
        a_inverse_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a_inverse);
        cond_l1 = a_norm_l1 * a_inverse_norm_l1;
        cond = Matrix.condition_hager(n, a);
        Console.WriteLine("  " + name.ToString(CultureInfo.InvariantCulture).PadLeft(20)
                               + "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                               + "  " + cond_l1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + cond.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        //
        //  KAHAN
        //
        name = "KAHAN";
        n = 4;
        alpha = 0.25;
        a = Matrix.kahan(alpha, n, n);
        a_inverse = Matrix.kahan_inverse(alpha, n);
        a_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a);
        a_inverse_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a_inverse);
        cond_l1 = a_norm_l1 * a_inverse_norm_l1;
        cond = Matrix.condition_hager(n, a);
        Console.WriteLine("  " + name.ToString(CultureInfo.InvariantCulture).PadLeft(20)
                               + "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                               + "  " + cond_l1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + cond.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        //
        //  Random
        //
        int seed = 123456789;

        for (i = 1; i <= 5; i++)
        {
            name = "RANDOM";
            n = 4;
            a = UniformRNG.r8mat_uniform_01_new(n, n, ref seed);
            double[] a_lu = typeMethods.r8mat_copy_new(n, n, a);
            int[] pivot = new int[n];
            typeMethods.r8ge_fa(n, ref a_lu, ref pivot);
            a_inverse = typeMethods.r8ge_inverse(n, a_lu, pivot);
            a_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a);
            a_inverse_norm_l1 = typeMethods.r8mat_norm_l1(n, n, a_inverse);
            cond_l1 = a_norm_l1 * a_inverse_norm_l1;
            cond = Matrix.condition_hager(n, a);
            Console.WriteLine("  " + name.ToString(CultureInfo.InvariantCulture).PadLeft(20)
                                   + "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + cond_l1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + cond.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }
}