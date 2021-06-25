using System;
using Burkardt.Uniform;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        
        public static void perm_check0 ( int n, int[] p )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM_CHECK0 checks a 0-based permutation.
            //
            //  Discussion:
            //
            //    The routine verifies that each of the integers from 0 to
            //    to N-1 occurs among the N entries of the permutation.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 October 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries.
            //
            //    Input, int P[N], the array to check.
            //
        {
            int ierror;
            int location;
            int value;

            for ( value = 0; value < n; value++ )
            {
                ierror = 1;

                for ( location = 0; location < n; location++ )
                {
                    if ( p[location] == value )
                    {
                        ierror = 0;
                        break;
                    }
                }

                if ( ierror != 0 )
                {
                    Console.WriteLine("");
                    Console.WriteLine("PERM_CHECK0 - Fatal error!");
                    Console.WriteLine("  Permutation is missing value " + value + "");
                    return;
                }

            }
        }
        
        public static int[] perm_uniform(int n, int base_, ref int seed)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM_UNIFORM selects a random permutation of N objects.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 October 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Albert Nijenhuis, Herbert Wilf,
            //    Combinatorial Algorithms,
            //    Academic Press, 1978, second edition,
            //    ISBN 0-12-519260-6.
            //
            //  Parameters:
            //
            //    Input, int N, the number of objects to be permuted.
            //
            //    Input, int BASE, is 0 for a 0-based permutation and 1 for 
            //    a 1-based permutation.
            //
            //    Input/output, int *SEED, a seed for the random number generator.
            //
            //    Output, int PERM_UNIFORM[N], a permutation of (BASE, BASE+1, ..., BASE+N-1).
            //
        {
            int[] p = new int[n];

            for (int i = 0; i < n; i++)
            {
                p[i] = i + base_;
            }

            for (int i = 0; i < n; i++)
            {
                int j = UniformRNG.i4_uniform(i, n - 1, ref seed);
                int k = p[i];
                p[i] = p[j];
                p[j] = k;
            }

            return p;
        }

        public static int[] perm_uniform_new(int n, ref int seed)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM_UNIFORM_NEW selects a random permutation of N objects.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 August 2012
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Albert Nijenhuis, Herbert Wilf,
            //    Combinatorial Algorithms for Computers and Calculators,
            //    Second Edition,
            //    Academic Press, 1978,
            //    ISBN: 0-12-519260-6,
            //    LC: QA164.N54.
            //
            //  Parameters:
            //
            //    Input, int N, the number of objects to be permuted.
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, int PERM_UNIFORM_NEW[N], a permutation of (1, 2, ..., N).
            //
        {
            int[] p = new int[n];

            for (int i = 0; i < n; i++)
            {
                p[i] = i + 1;
            }

            for (int i = 0; i < n; i++)
            {
                int j = UniformRNG.i4_uniform_ab(i, n - 1, ref seed);
                int k = p[i];
                p[i] = p[j];
                p[j] = k;
            }

            return p;
        }

        public static bool perm_check(int n, int[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM_CHECK checks that a vector represents a permutation.
            //
            //  Discussion:
            //
            //    The routine verifies that each of the integers from 1
            //    to N occurs among the N entries of the permutation.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 January 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries.
            //
            //    Input, int P[N], the array to check.
            //
            //    Output, bool PERM_CHECK, is TRUE if the permutation is OK.
            //
        {
            bool found;
            int i;
            int seek;

            for (seek = 1; seek <= n; seek++)
            {
                found = false;

                for (i = 0; i < n; i++)
                {
                    if (p[i] == seek)
                    {
                        found = true;
                        break;
                    }
                }

                if (!found)
                {
                    return false;
                }

            }

            return true;
        }

        public static void _perm_check(int n, int[] p)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM_CHECK checks that a vector represents a permutation.
            //
            //  Discussion:
            //
            //    The routine verifies that each of the integers from 1
            //    to N occurs among the N entries of the permutation.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    06 May 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries.
            //
            //    Input, int P[N], the permutation, in standard index form.
            //
        {
            bool error;
            int ifind;
            int iseek;

            for (iseek = 1; iseek <= n; iseek++)
            {
                error = true;

                for (ifind = 1; ifind <= n; ifind++)
                {
                    if (p[ifind - 1] == iseek)
                    {
                        error = false;
                        break;
                    }
                }

                if (error)
                {
                    Console.WriteLine();
                    Console.WriteLine("PERM_CHECK - Fatal error!");
                    Console.WriteLine("  The input permutation is not legal.");
                }
            }
        }

        public static void perm_print(int n, int[] p, string title)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM_PRINT prints a permutation.
            //
            //  Discussion:
            //
            //    The permutation is assumed to be zero-based.
            //
            //  Example:
            //
            //    Input:
            //
            //      P = 6 1 2 0 4 2 5
            //
            //    Printed output:
            //
            //      "This is the permutation:"
            //
            //      0 1 2 3 4 5 6
            //      6 1 2 0 4 2 5
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 May 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of objects permuted.
            //
            //    Input, int P[N], the permutation, in standard index form.
            //
            //    Input, string TITLE, a title.
            //    If no title is supplied, then only the permutation is printed.
            //
        {
            int inc = 20;

            if (s_len_trim(title) != 0)
            {
                Console.WriteLine();
                Console.WriteLine(title);

                for (int ilo = 0; ilo < n; ilo = ilo + inc)
                {
                    int ihi = ilo + inc;
                    if (n < ihi)
                    {
                        ihi = n;
                    }

                    Console.WriteLine();
                    string cout = "  ";
                    for (int i = ilo; i < ihi; i++)
                    {
                        cout += i.ToString().PadLeft(4);
                    }

                    Console.WriteLine(cout);

                    Console.WriteLine();
                    cout = "  ";
                    for (int i = ilo; i < ihi; i++)
                    {
                        cout += p[i].ToString().PadLeft(4);
                    }

                    Console.WriteLine(cout);
                    Console.WriteLine();
                }
            }
            else
            {
                for (int ilo = 0; ilo < n; ilo = ilo + inc)
                {
                    int ihi = ilo + inc;
                    if (n < ihi)
                    {
                        ihi = n;
                    }

                    string cout = "  ";
                    for (int i = ilo; i < ihi; i++)
                    {
                        cout += p[i].ToString().PadLeft(4);
                    }

                    Console.WriteLine(cout);
                    Console.WriteLine();
                }
            }
        }

        public static bool perm0_check(int n, int[] p)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM0_CHECK checks a permutation of ( 0, ..., N-1 ).
            //
            //  Discussion:
            //
            //    The routine verifies that each of the integers from 0 to
            //    to N-1 occurs among the N entries of the permutation.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 May 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries.
            //
            //    Input, int P[N], the array to check.
            //
            //    Output, bool PERM0_CHECK, is 
            //    TRUE if P is a legal permutation of 0,...,N-1.
            //    FALSE if P is not a legal permuation of 0,...,N-1.
            //
        {
            bool check;
            int location;
            int value;

            check = true;

            for (value = 0; value < n; value++)
            {
                check = false;

                for (location = 0; location < n; location++)
                {
                    if (p[location] == value)
                    {
                        check = true;
                        break;
                    }
                }

                if (!check)
                {
                    Console.WriteLine("");
                    Console.WriteLine("PERM0_CHECK - Fatal error!");
                    Console.WriteLine("  Permutation is missing value " + value + "");
                    break;
                }

            }

            return check;
        }

        public static int[] perm0_uniform_new(int n, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM0_UNIFORM_NEW selects a random permutation of 0,...,N-1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 May 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Albert Nijenhuis, Herbert Wilf,
            //    Combinatorial Algorithms,
            //    Academic Press, 1978, second edition,
            //    ISBN 0-12-519260-6.
            //
            //  Parameters:
            //
            //    Input, int N, the number of objects to be permuted.
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, int PERM0_UNIFORM_NEW[N], a permutation of
            //    (0, 1, ..., N-1).
            //
        {
            int i;
            int j;
            int k;
            int[] p;

            p = new int[n];

            for (i = 0; i < n; i++)
            {
                p[i] = i;
            }

            for (i = 0; i < n; i++)
            {
                j = UniformRNG.i4_uniform_ab(i, n - 1, ref seed);
                k = p[i];
                p[i] = p[j];
                p[j] = k;
            }

            return p;
        }

        public static bool perm1_check(int n, int[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM1_CHECK checks a permutation of (1, ..., N ).
            //
            //  Discussion:
            //
            //    The routine verifies that each of the integers from 0 to
            //    to N-1 occurs among the N entries of the permutation.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 May 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries.
            //
            //    Input, int P[N], the array to check.
            //
            //    Output, bool PERM1_CHECK, is 
            //    TRUE if P is a legal permutation of 1,...,N.
            //    FALSE if P is not a legal permuation of 1,...,N.
            //
        {
            bool check;
            int location;
            int value;

            check = true;

            for (value = 1; value <= n; value++)
            {
                check = false;

                for (location = 0; location < n; location++)
                {
                    if (p[location] == value)
                    {
                        check = true;
                        break;
                    }
                }

                if (!check)
                {
                    Console.WriteLine("");
                    Console.WriteLine("PERM1_CHECK - Fatal error!");
                    Console.WriteLine("  Permutation is missing value " + value + "");
                    break;
                }

            }

            return check;
        }

        public static int[] perm1_uniform_new(int n, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM1_UNIFORM_NEW selects a random permutation of 1,...,N.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 May 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Albert Nijenhuis, Herbert Wilf,
            //    Combinatorial Algorithms,
            //    Academic Press, 1978, second edition,
            //    ISBN 0-12-519260-6.
            //
            //  Parameters:
            //
            //    Input, int N, the number of objects to be permuted.
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, int PERM1_UNIFORM_NEW[N], a permutation of
            //    (1, ..., N).
            //
        {
            int i;
            int j;
            int k;
            int[] p;

            p = new int[n];

            for (i = 0; i < n; i++)
            {
                p[i] = i + 1;
            }

            for (i = 0; i < n; i++)
            {
                j = UniformRNG.i4_uniform_ab(i, n - 1, ref seed);
                k = p[i];
                p[i] = p[j];
                p[j] = k;
            }

            return p;
        }


    }
}