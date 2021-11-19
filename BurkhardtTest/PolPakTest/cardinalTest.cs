using System;
using Burkardt.Function;
using Burkardt.Types;

namespace PolPakTest;

public static class cardinalTest
{
    public static void cardinal_cos_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CARDINAL_COS_TEST tests CARDINAL_COS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] c;
        int i;
        int j;
        int m = 11;
        const double r8_pi = 3.141592653589793;
        double[] t;

        Console.WriteLine("");
        Console.WriteLine("CARDINAL_COS_TEST");
        Console.WriteLine("  CARDINAL_COS evaluates cardinal cosine functions.");
        Console.WriteLine("  Ci(Tj) = Delta(i,j), where Tj = cos(pi*i/(n+1)).");
        Console.WriteLine("  A simple check of all pairs should form the identity matrix.");

        Console.WriteLine("");
        Console.WriteLine("  The CARDINAL_COS test matrix:");
        Console.WriteLine("");

        t = typeMethods.r8vec_linspace_new(m + 2, 0.0, r8_pi);

        for (j = 0; j <= m + 1; j++)
        {
            c = Cardinal.cardinal_cos(j, m, m + 2, t);
            string cout = "";
            for (i = 0; i <= m + 1; i++)
            {
                cout += "  " + c[i].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
        }

    }

    public static void cardinal_sin_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CARDINAL_SIN_TEST tests CARDINAL_SIN.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        int m = 11;
        const double r8_pi = 3.141592653589793;
        double[] s;
        double[] t;

        Console.WriteLine("");
        Console.WriteLine("CARDINAL_SIN_TEST");
        Console.WriteLine("  CARDINAL_SIN evaluates cardinal sine functions.");
        Console.WriteLine("  Si(Tj) = Delta(i,j), where Tj = cos(pi*i/(n+1)).");
        Console.WriteLine("  A simple check of all pairs should form the identity matrix.");

        t = typeMethods.r8vec_linspace_new(m + 2, 0.0, r8_pi);

        Console.WriteLine("");
        Console.WriteLine("  The CARDINAL_SIN test matrix:");
        Console.WriteLine("");
        for (j = 0; j <= m + 1; j++)
        {
            s = Cardinal.cardinal_sin(j, m, m + 2, t);
            string cout = "";
            for (i = 0; i <= m + 1; i++)
            {
                cout += "  " + s[i].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
        }

    }

}