using System;

namespace Burkardt
{
    public static class Levenshtein
    {
        public static int levenshtein_distance(int m, string s, int n, string t)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //   LEVENSHTEIN_DISTANCE computes the Levenshtein distance between strings.
            //
            //  Discussion:
            //
            //    Let S and T be source and target strings.  Consider the task of
            //    converting S to T in the minimal number of steps, involving
            //    * Insertion: insert a new character
            //    * Deletion: delete a character
            //    * Substitution: swap one character for another.
            //    The Levenshtein distance from S to T is the minimal number of such
            //    steps required to transform S into T.
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
            //  Parameters:
            //
            //    Input, int M, the length of string S.
            //
            //    Input, string S, the first string.
            //
            //    Input, int N, the length of string T.
            //
            //    Input, string T, the second string.
            //
            //    Output, int LEVENSHTEIN_DISTANCE, the Levenshtein distance between the
            //    two strings.
            //
        {
            int[] d;
            int distance;
            int i;
            int j;
            int substitution_cost;

            d = new int[(m + 1) * (n + 1)];

            d[0 + 0 * (m + 1)] = 0;
            //
            //  Source prefixes can be transformed into empty string by
            //  dropping all characters,
            //
            for (i = 1; i <= m; i++)
            {
                d[i + 0 * (m + 1)] = i;
            }

            //
            //  Target prefixes can be reached from empty source prefix
            //  by inserting every character.
            //
            for (j = 1; j <= n; j++)
            {
                d[0 + j * (m + 1)] = j;
            }

            for (j = 1; j <= n; j++)
            {
                for (i = 1; i <= m; i++)
                {
                    if (s[i - 1] == t[j - 1])
                    {
                        substitution_cost = 0;
                    }
                    else
                    {
                        substitution_cost = 1;
                    }

                    d[i + j * (m + 1)] = Math.Min(d[i - 1 + j * (m + 1)] + 1,
                        Math.Min(d[i + (j - 1) * (m + 1)] + 1,
                            d[i - 1 + (j - 1) * (m + 1)] + substitution_cost));
                }
            }

            distance = d[m + n * (m + 1)];

            return distance;
        }

        public static int[] levenshtein_matrix(int m, string s, int n, string t)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //   LEVENSHTEIN_MATRIX computes the Levenshtein distance matrix between strings.
            //
            //  Discussion:
            //
            //    Let S and T be source and target strings.  Consider the task of
            //    converting S to T in the minimal number of steps, involving
            //    * Insertion: insert a new character
            //    * Deletion: delete a character
            //    * Substitution: swap one character for another.
            //    The Levenshtein distance from S to T is the minimal number of such
            //    steps required to transform S into T.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 March 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //  Parameters:
            //
            //    Input, int M, the length of string S.
            //
            //    Input, string S, the first string.
            //
            //    Input, int N, the length of string T.
            //
            //    Input, string T, the second string.
            //
            //    Output, int LEVENSHTEIN_MATRIX[(M+1)*(N+1)], the Levenshtein matrix.
            //
        {
            int[] d;
            int i;
            int j;
            int substitution_cost;

            d = new int[(m + 1) * (n + 1)];

            d[0 + 0 * (m + 1)] = 0;
            //
            //  Source prefixes can be transformed into empty string by
            //  dropping all characters,
            //
            for (i = 1; i <= m; i++)
            {
                d[i + 0 * (m + 1)] = i;
            }

            //
            //  Target prefixes can be reached from empty source prefix
            //  by inserting every character.
            //
            for (j = 1; j <= n; j++)
            {
                d[0 + j * (m + 1)] = j;
            }

            for (j = 1; j <= n; j++)
            {
                for (i = 1; i <= m; i++)
                {
                    if (s[i - 1] == t[j - 1])
                    {
                        substitution_cost = 0;
                    }
                    else
                    {
                        substitution_cost = 1;
                    }

                    d[i + j * (m + 1)] = Math.Min(d[i - 1 + j * (m + 1)] + 1,
                        Math.Min(d[i + (j - 1) * (m + 1)] + 1,
                            d[i - 1 + (j - 1) * (m + 1)] + substitution_cost));
                }
            }

            return d;
        }
    }
}