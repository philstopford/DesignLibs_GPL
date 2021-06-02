using System;
using Burkardt;

namespace Burkardt.BacktrackBinaryTest
{
    class Program
    {
        static void Main(string[] args)
//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for BACKTRACK_BINARY_RC_TEST.
//
//  Discussion:
//
//    BACKTRACK_BINARY_RC_TEST tests BACKTRACK_BINARY_RC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 January 2014
//
//  Author:
//
//    John Burkardt
//
        {
            Console.WriteLine("");
            Console.WriteLine("BACKTRACK_BINARY_RC_TEST:");
            Console.WriteLine("  C++ version");
            Console.WriteLine("  Test the BACKTRACK_BINARY_RC library.");

            test01();
            test02();

            Console.WriteLine("");
            Console.WriteLine("BACKTRACK_BINARY_RC_TEST:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01()
//****************************************************************************80
//
//  Purpose:
//
//    TEST01 seeks a selection of binary powers that have a given sum.
//
//  Discussion:
//
//    We consider the binary powers 1, 2, 4, ... 2^(n-1).
//
//    We wish to select some of these powers, so that the sum is equal
//    to a given target value.  We are actually simply seeking the binary
//    representation of an integer.
//
//    A partial solution is acceptable if it is less than the target value.
//
//    We list the powers in descending order, so that the bactracking
//    procedure makes the most significant choices first, thus quickly
//    eliminating many unsuitable choices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 January 2014
//
//  Author:
//
//    John Burkardt
//
        {
            int call_num;
            int[] choice = new int[8];
            int factor;
            int i;
            int n = 8;
            int n2;
            bool reject = false;
            int result;
            int target;
            int[] targets =  {
                73, 299, -3
            }
            ;
            int test;
            int test_num = 3;

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  Use BACKBIN_RC to find the binary expansion of");
            Console.WriteLine("  an integer between 0 and 255.");
            Console.WriteLine("  The choices are 0/1 for the 8 digits.");

            for (test = 0; test < test_num; test++)
            {
                target = targets[test];
                Console.WriteLine("");
                Console.WriteLine("  TARGET = " + target + "");
                call_num = 0;
                n2 = -1;

                for (;;)
                {
                    BacktrackBinary.backbin_rc(n, reject, ref n2, ref choice);
                    call_num = call_num + 1;

                    if (n2 == -1)
                    {
                        Console.WriteLine("  Termination without solution.");
                        break;
                    }

//
//  Evaluate the integer determined by the choices.
//
                    factor = 1;
                    for (i = n; n2 < i; i--)
                    {
                        factor = factor * 2;
                    }

                    result = 0;
                    for (i = 0; i < n2; i++)
                    {
                        result = result * 2 + choice[i];
                    }

                    result = result * factor;
//
//  If the integer is too big, then we reject it, and
//  all the related integers formed by making additional choices.
//
                    reject = (target < result);
//
//  If we hit the target, then in this case, we can exit because
//  the solution is unique.
//
                    if (result == target)
                    {
                        break;
                    }

                }

                Console.WriteLine("  Number of calls = " + call_num + "");
                Console.WriteLine("  Binary search space = " + Math.Pow(2, n) + "");
                string cout = "  ";
                for (i = 0; i < n; i++)
                {
                    cout += choice[i].ToString().PadLeft(2);
                }

                Console.WriteLine(cout);
            }
       }

        static void test02()
//****************************************************************************80
//
//  Purpose:
//
//    TEST02 seeks a subset of a set of numbers which add to a given sum.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 January 2014
//
//  Author:
//
//    John Burkardt
//
        {
            int call_num;
            int[] choice = new int[8];
            int i;
            int n = 8;
            int n2;
            bool reject = false;
            int result;
            int target = 53;
            int[] w =  {
                15, 22, 14, 26, 32, 9, 16, 8
            }
            ;

            Console.WriteLine("");
            Console.WriteLine("TEST02");
            Console.WriteLine("  Use BACKBIN_RC to seek subsets of a set W");
            Console.WriteLine("  that sum to a given target value.");
            Console.WriteLine("  The choices are 0/1 to select each element of W.");

            Console.WriteLine("");
            Console.WriteLine("  TARGET = " + target + "");
            Console.WriteLine("");
            call_num = 0;
            n2 = -1;

            for (;;)
            {
                BacktrackBinary.backbin_rc(n, reject, ref n2, ref choice);
                call_num = call_num + 1;

                if (n2 == -1)
                {
                    break;
                }

//
//  Evaluate the partial sum.
//
                result = 0;
                for (i = 0; i < n2; i++)
                {
                    result = result + choice[i] * w[i];
                }

//
//  If the sum is too big, then we reject it, and
//  all the related sums formed by making additional choices.
//
                reject = (target < result);
//
//  If we hit the target, print out the information.
//
                if (result == target && n2 == n)
                {
                    string cout = "  ";
                    for (i = 0; i < n; i++)
                    {
                        cout += choice[i].ToString().PadLeft(2);
                    }

                    Console.WriteLine(cout);
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  Number of calls = " + call_num + "");
            Console.WriteLine("  Binary search space = " + Math.Pow(2, n) + "");
        }
    }
}