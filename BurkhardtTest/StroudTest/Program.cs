using System;

namespace StroudTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for STROUD_TEST.
        //
        //  Discussion:
        //
        //    STROUD_TEST tests.tests the STROUD library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 February 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("STROUD_TEST");
        Console.WriteLine("  Test the STROUD library.");
            
        tests.test01 (  );
        tests.test02 (  );
        tests.test03 (  );
        tests.test04 (  );
        tests.test045 (  );
        tests.test05 (  );
        tests.test052 (  );
        tests.test054 (  );
        tests.test07 (  );
        tests.test08 (  );
        tests.test085 (  );
        tests.test09 (  );

        tests2.test10 (  );
        tests2.test11 (  );
        tests2.test12 (  );
        tests2.test13 (  );
        tests2.test14 (  );
        tests2.test15 (  );
        tests2.test16 (  );
        tests2.test163 (  );
        tests2.test165 (  );
        tests2.test167 (  );
        tests2.test17 (  );
        tests2.test18 (  );
        tests2.test19 (  );

        tests3.test20 (  );
        tests3.test205 (  );
        tests3.test207 (  );
        tests3.test2075 (  );
        tests3.test208 (  );
        tests3.test21 (  );
        tests3.test215 (  );
        tests3.test22 (  );
        tests3.test23 (  );
        tests3.test24 (  );
        tests3.test25 (  );
        tests3.test255 (  );
        tests3.test26 (  );
        tests3.test27 (  );
        tests3.test28 (  );
        tests3.test29 (  );

        tests4.test30 (  );
        tests4.test31 (  );
        tests4.test32 (  );
        tests4.test322 (  );
        tests4.test324 (  );
        tests4.test326 (  );
        tests4.test33 (  );
        tests4.test335 (  );
        tests4.test34 (  );
        tests4.test345 (  );
        tests4.test35 (  );
        tests4.test36 (  );
        tests4.test37 (  );
        tests4.test38 (  );
        tests4.test39 (  );

        tests5.test40 (  );
        tests5.test41 (  );
        tests5.test42 (  );
        tests5.test425 (  );
        tests5.test43 (  );
        tests5.test44 (  );
        tests5.test45 (  );
        tests5.test46 (  );
        tests5.test47 (  );
        tests5.test48 (  );
        tests5.test49 (  );

        Console.WriteLine("");
        Console.WriteLine("STROUD_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}