using System;
using Burkardt.AppliedStatistics;

namespace Burkardt.ASA299Test
{
    class Program
    {
        static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for ASA299_TEST.
        //
        //  Discussion:
        //
        //    ASA299_TEST tests the ASA299 library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            Console.WriteLine("");
            Console.WriteLine("ASA299_TEST");
            Console.WriteLine("  Test the ASA299 library.");

            test01 ( );

            Console.WriteLine("");
            Console.WriteLine("ASA299_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
            
        }
        
        static void test01 ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests SIMPLEX_LATTICE_POINT_NEXT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            int N = 4;

            int i;
            int j;
            bool more;
            int t = 4;
            int[] x = new int[N];

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  SIMPLEX_LATTICE_POINT_NEXT generates lattice points");
            Console.WriteLine("  in the simplex");
            Console.WriteLine("    0 <= X");
            Console.WriteLine("    sum ( X(1:N) ) <= T");
            Console.WriteLine("  Here N = " + N + "");
            Console.WriteLine("  and T =  " + t + "");
            Console.WriteLine("");
            Console.WriteLine("     Index        X(1)      X(2)      X(3)      X(4)");
            Console.WriteLine("");

            more = false;

            i = 0;

            for ( ; ; )
            {
                Algorithms.simplex_lattice_point_next ( N, t, ref more, ref x );

                i = i + 1;

                string cout = "  " + i.ToString().PadLeft(8);
                cout += "  ";
                for ( j = 0; j < N; j++ )
                {
                    cout += x[j].ToString().PadLeft(8);
                }
                Console.WriteLine(cout);

                if ( !more )
                {
                    break;
                }
            }
        }
    }
}