﻿using System;
using System.Numerics;
using Burkardt.TOMSNS;
using TestValues;

namespace TOMS243Test
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TOMS243_TEST tests TOMS243.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 January 2019
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Complex fx1 = new Complex(0, 0);
            Complex fx2;
            int n_data;
            Complex x = new Complex(0, 0);

            Console.WriteLine("");
            Console.WriteLine("TOMS243_TEST:");
            Console.WriteLine("  TOMS243 computes the natural logarithm of a complex value.");
            Console.WriteLine("");
            Console.WriteLine("               X                               FX exact");
            Console.WriteLine("                                               FX computed");
            Console.WriteLine("");

            n_data = 0;

            while (true)
            {
                Cmplex.c8_log_values(ref n_data, ref x, ref fx1);

                if (n_data == 0)
                {
                    break;
                }

                fx2 = TOMS.toms243(x);

                Console.WriteLine("  ( "
                                  + (x.Real.ToString("0.####").PadLeft(8)) + ","
                                  + (x.Imaginary.ToString("0.####").PadLeft(8)) + ")  ( "
                                  + (fx1.Real.ToString("0.############").PadLeft(18)) + ","
                                  + (fx1.Imaginary.ToString("0.############").PadLeft(18)) + ")");
                Console.WriteLine("                        ( "
                                  + (fx2.Real.ToString("0.############").PadLeft(18)) + ","
                                  + (fx2.Imaginary.ToString("0.############").PadLeft(18)) + ")");
            }

            Console.WriteLine("");
            Console.WriteLine("TOMS243_TEST:");
            Console.WriteLine("  Normal end of execution:");
            Console.WriteLine("");
        }
    }
}