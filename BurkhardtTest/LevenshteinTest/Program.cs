using System;
using Burkardt.MatrixNS;

namespace LevenshteinTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //   MAIN is the main program for LEVENSHTEIN_TEST.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 March 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("LEVENSHTEIN_TEST");
        Console.WriteLine("  Test the LEVENSHTEIN library.");

        levenshtein_distance_test();
        levenshtein_matrix_test();
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("LEVENSHTEIN_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void levenshtein_distance_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //   Levenshtein.levenshtein_distance_TEST tests Levenshtein.levenshtein_distance.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 March 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const string s1 = "water";
        const string s2 = "kitten";
        const string s3 = "saturday";
        const string s4 = "pheromones";
        const string t1 = "wine";
        const string t2 = "sitting";
        const string t3 = "sunday";
        const string t4 = "photographer";

        Console.WriteLine("");
        Console.WriteLine("Levenshtein.levenshtein_distance_TEST:");
        Console.WriteLine("  Levenshtein.levenshtein_distance computes the Levenshtein distance");
        Console.WriteLine("  between two strings.");

        int m = s1.Length;
        int n = t1.Length;
        int d1 = Levenshtein.levenshtein_distance(m, s1, n, t1);
        int d2 = 3;
        Console.WriteLine("");
        Console.WriteLine("  S = '" + s1 + "'");
        Console.WriteLine("  T = '" + t1 + "'");
        Console.WriteLine("  Computed distance = " + d1 + "");
        Console.WriteLine("  Correct distance  = " + d2 + "");

        m = s2.Length;
        n = t2.Length;
        d1 = Levenshtein.levenshtein_distance(m, s2, n, t2);
        d2 = 3;
        Console.WriteLine("");
        Console.WriteLine("  S = '" + s2 + "'");
        Console.WriteLine("  T = '" + t2 + "'");
        Console.WriteLine("  Computed distance = " + d1 + "");
        Console.WriteLine("  Correct distance  = " + d2 + "");

        m = s3.Length;
        n = t3.Length;
        d1 = Levenshtein.levenshtein_distance(m, s3, n, t3);
        d2 = 3;
        Console.WriteLine("");
        Console.WriteLine("  S = '" + s3 + "'");
        Console.WriteLine("  T = '" + t3 + "'");
        Console.WriteLine("  Computed distance = " + d1 + "");
        Console.WriteLine("  Correct distance  = " + d2 + "");

        m = s4.Length;
        n = t4.Length;
        d1 = Levenshtein.levenshtein_distance(m, s4, n, t4);
        d2 = 8;
        Console.WriteLine("");
        Console.WriteLine("  S = '" + s4 + "'");
        Console.WriteLine("  T = '" + t4 + "'");
        Console.WriteLine("  Computed distance = " + d1 + "");
        Console.WriteLine("  Correct distance  = " + d2 + "");
    }

    private static void levenshtein_matrix_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        // -  Levenshtein.levenshtein_matrix_TEST tests Levenshtein.levenshtein_matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 March 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        const string s1 = "water";
        const string s2 = "kitten";
        const string s3 = "saturday";
        const string s4 = "pheromones";
        const string t1 = "wine";
        const string t2 = "sitting";
        const string t3 = "sunday";
        const string t4 = "photographer";
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine("Levenshtein.levenshtein_matrix_TEST:");
        Console.WriteLine("  Levenshtein.levenshtein_matrix computes the Levenshtein matrix");
        Console.WriteLine("  associated with the computation of the Levenshtein");
        Console.WriteLine("  distance between two strings.");

        int m = s1.Length;
        int n = t1.Length;
        int[] d = Levenshtein.levenshtein_matrix(m, s1, n, t1);
        Console.WriteLine("");
        Console.WriteLine("  S = '" + s1 + "'");
        Console.WriteLine("  T = '" + t1 + "'");
        for (i = 0; i <= m; i++)
        {
            for (j = 0; j <= n; j++)
            {
                cout += " " + d[i + j * (m + 1)].ToString().PadLeft(2);
            }

            Console.WriteLine(cout);
            cout = "";
        }

        m = s2.Length;
        n = t2.Length;
        d = Levenshtein.levenshtein_matrix(m, s2, n, t2);
        Console.WriteLine("");
        Console.WriteLine("  S = '" + s2 + "'");
        Console.WriteLine("  T = '" + t2 + "'");
        for (i = 0; i <= m; i++)
        {
            for (j = 0; j <= n; j++)
            {
                cout += " " + d[i + j * (m + 1)].ToString().PadLeft(2);
            }

            Console.WriteLine(cout);
            cout = "";
        }

        m = s3.Length;
        n = t3.Length;
        d = Levenshtein.levenshtein_matrix(m, s3, n, t3);
        Console.WriteLine("");
        Console.WriteLine("  S = '" + s3 + "'");
        Console.WriteLine("  T = '" + t3 + "'");
        for (i = 0; i <= m; i++)
        {
            for (j = 0; j <= n; j++)
            {
                cout += " " + d[i + j * (m + 1)].ToString().PadLeft(2);
            }

            Console.WriteLine(cout);
            cout = "";
        }

        m = s4.Length;
        n = t4.Length;
        d = Levenshtein.levenshtein_matrix(m, s4, n, t4);
        Console.WriteLine("");
        Console.WriteLine("  S = '" + s4 + "'");
        Console.WriteLine("  T = '" + t4 + "'");
        for (i = 0; i <= m; i++)
        {
            for (j = 0; j <= n; j++)
            {
                cout += " " + d[i + j * (m + 1)].ToString().PadLeft(2);
            }

            Console.WriteLine(cout);
            cout = "";
        }
    }
}