﻿using System;
using Burkardt.Probability;

namespace Burkardt.ProbabilityTest
{
    partial class Program
    {
        static void tfn_test()

//****************************************************************************80
//
//  Purpose:
//
//    TFN_TEST tests TFN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 April 2016
//
//  Author:
//
//    John Burkardt
//
        {
            double a = 0;
            double h = 0;
            double t = 0;

            Console.WriteLine("");
            Console.WriteLine("TFN_TEST");
            Console.WriteLine("  TFN evaluates Owen's T function;");
            Console.WriteLine("");
            Console.WriteLine("      H             A           T(H,A)          Exact");
            Console.WriteLine("");

            int n_data = 0;

            for (;;)
            {
                Owen.owen_values(ref n_data, ref h, ref a, ref t);

                if (n_data <= 0)
                {
                    break;
                }

                double t2 = Owen.tfn(h, a);

                Console.WriteLine("  "
                                  + h.ToString().PadLeft(14) + "  "
                                  + a.ToString().PadLeft(14) + "  "
                                  + t2.ToString().PadLeft(14) + "  "
                                  + t.ToString().PadLeft(14) + "");
            }
        }
        
    }
}