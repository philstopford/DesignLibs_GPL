﻿using System;
using System.Globalization;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{

    public static void perm_check0(int n, int[] p)

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
        int value;

        for (value = 0; value < n; value++)
        {
            int ierror = 1;

            int location;
            for (location = 0; location < n; location++)
            {
                if (p[location] != value)
                {
                    continue;
                }

                ierror = 0;
                break;
            }

            if (ierror == 0)
            {
                continue;
            }

            Console.WriteLine("");
            Console.WriteLine("PERM_CHECK0 - Fatal error!");
            Console.WriteLine("  Permutation is missing value " + value + "");
            return;

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
            (p[i], p[j]) = (p[j], p[i]);
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
            (p[i], p[j]) = (p[j], p[i]);
        }

        return p;
    }

    public static bool perm_check(int n, int[] p, int base_ = 1)

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
        int seek;

        for (seek = base_; seek < base_ + n; seek++)
        {
            bool found = false;

            int i;
            for (i = 0; i < n; i++)
            {
                if (p[i] != seek)
                {
                    continue;
                }

                found = true;
                break;
            }

            switch (found)
            {
                case false:
                    Console.WriteLine("");
                    Console.WriteLine("PERM_CHECK - Fatal error!");
                    Console.WriteLine("  Did not find " + found + "");
                    return false;
            }

        }

        return true;
    }
        
    public static int[] perm_inverse3 ( int n, int[] perm )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM_INVERSE3 produces the inverse of a given permutation.
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
        //    Input, int N, the number of items permuted.
        //
        //    Input, int PERM[N], a permutation.
        //
        //    Output, int PERM_INVERSE3[N], the inverse permutation.
        //
    {
        int i;

        int[] perm_inv = new int[n];

        for ( i = 0; i < n; i++ )
        {
            perm_inv[(perm[i] + perm_inv.Length) % perm_inv.Length] = i;
        }

        return perm_inv;
    }
        
    public static void perm_inverse3 ( int n, int[] perm, ref int[] perm_inv )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM_INVERSE3 produces the inverse of a given permutation.
        //
        //  Discussion:
        //
        //    This function assumes the permutation is 0-based.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 September 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of items permuted.
        //
        //    Input, int PERM[N], a permutation.
        //
        //    Output, int PERM_INV[N], the inverse permutation.
        //
    {
        int i;

        for ( i = 0; i < n; i++ )
        {
            perm_inv[perm[i]] = i;
        }
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
        int iseek;

        for (iseek = 1; iseek <= n; iseek++)
        {
            bool error = true;

            int ifind;
            for (ifind = 1; ifind <= n; ifind++)
            {
                if (p[ifind - 1] != iseek)
                {
                    continue;
                }

                error = false;
                break;
            }

            switch (error)
            {
                case true:
                    Console.WriteLine();
                    Console.WriteLine("PERM_CHECK - Fatal error!");
                    Console.WriteLine("  The input permutation is not legal.");
                    break;
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
        const int inc = 20;

        if (s_len_trim(title) != 0)
        {
            Console.WriteLine();
            Console.WriteLine(title);

            for (int ilo = 0; ilo < n; ilo += inc)
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
                    cout += i.ToString(CultureInfo.InvariantCulture).PadLeft(4);
                }

                Console.WriteLine(cout);

                Console.WriteLine();
                cout = "  ";
                for (int i = ilo; i < ihi; i++)
                {
                    cout += p[i].ToString(CultureInfo.InvariantCulture).PadLeft(4);
                }

                Console.WriteLine(cout);
                Console.WriteLine();
            }
        }
        else
        {
            for (int ilo = 0; ilo < n; ilo += inc)
            {
                int ihi = ilo + inc;
                if (n < ihi)
                {
                    ihi = n;
                }

                string cout = "  ";
                for (int i = ilo; i < ihi; i++)
                {
                    cout += p[i].ToString(CultureInfo.InvariantCulture).PadLeft(4);
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
        int value;

        bool check = true;

        for (value = 0; value < n; value++)
        {
            check = false;

            int location;
            for (location = 0; location < n; location++)
            {
                if (p[location] != value)
                {
                    continue;
                }

                check = true;
                break;
            }

            if (check)
            {
                continue;
            }

            Console.WriteLine("");
            Console.WriteLine("PERM0_CHECK - Fatal error!");
            Console.WriteLine("  Permutation is missing value " + value + "");
            break;

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

        int[] p = new int[n];

        for (i = 0; i < n; i++)
        {
            p[i] = i;
        }

        for (i = 0; i < n; i++)
        {
            int j = UniformRNG.i4_uniform_ab(i, n - 1, ref seed);
            (p[i], p[j]) = (p[j], p[i]);
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
        int value;

        bool check = true;

        for (value = 1; value <= n; value++)
        {
            check = false;

            int location;
            for (location = 0; location < n; location++)
            {
                if (p[location] != value)
                {
                    continue;
                }

                check = true;
                break;
            }

            if (check)
            {
                continue;
            }

            Console.WriteLine("");
            Console.WriteLine("PERM1_CHECK - Fatal error!");
            Console.WriteLine("  Permutation is missing value " + value + "");
            break;

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

        int[] p = new int[n];

        for (i = 0; i < n; i++)
        {
            p[i] = i + 1;
        }

        for (i = 0; i < n; i++)
        {
            int j = UniformRNG.i4_uniform_ab(i, n - 1, ref seed);
            (p[i], p[j]) = (p[j], p[i]);
        }

        return p;
    }

    public static bool perm_check2(int n, int[] p, int base_)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM_CHECK2 checks that a vector represents a permutation.
        //
        //  Discussion:
        //
        //    The routine verifies that each of the integers from BASE to
        //    to BASE+N-1 occurs among the N entries of the permutation.
        //
        //    Set the input quantity BASE to 0, if P is a 0-based permutation,
        //    or to 1 if P is a 1-based permutation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 June 2009
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
        //    Input, int BASE, the index base.
        //
        //    Output, bool PERM_CHECK2, is TRUE if the permutation is OK.
        //
    {
        int seek;

        for (seek = base_; seek < base_ + n; seek++)
        {
            bool found = false;

            int i;
            for (i = 0; i < n; i++)
            {
                if (p[i] != seek)
                {
                    continue;
                }

                found = true;
                break;
            }

            switch (found)
            {
                case false:
                    return false;
            }

        }

        return true;
    }

    public static void perm_inverse(int n, ref int[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM_INVERSE inverts a permutation "in place".
        //
        //  Discussion:
        //
        //    This algorithm assumes that the entries in the permutation vector are
        //    strictly positive.  In particular, the value 0 must not occur.
        //
        //    When necessary, this function shifts the data temporarily so that
        //    this requirement is satisfied.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of objects being permuted.
        //
        //    Input/output, int P[N], the permutation, in standard index form.
        //    On output, P describes the inverse permutation
        //
    {
        int i;
        int i1;
        int i2;

        switch (n)
        {
            case <= 0:
                Console.WriteLine("");
                Console.WriteLine("PERM_INVERSE - Fatal error!");
                Console.WriteLine("  Input value of N = " +n +"");
                return;
        }

        //
        //  Find the least value, and shift data so it begins at 1.
        //
        int p_min = i4vec_min(n, p);
        const int base_ = 1;

        for (i = 0; i < n; i++)
        {
            p[i] = p[i] - p_min + base_;
        }

        //
        //  Now we can safely check the permutation.
        //
        if (!perm_check2(n, p, base_))
        {
            Console.WriteLine("");
            Console.WriteLine("PERM_INVERSE - Fatal error!");
            Console.WriteLine("  PERM_CHECK rejects this permutation.");
            return;
        }
        //
        //  Now we can invert the permutation.
        //

        for (i = 1; i <= n; i++)
        {
            i1 = p[i - 1];

            while (i < i1)
            {
                i2 = p[i1 - 1];
                p[i1 - 1] = -i2;
                i1 = i2;
            }

            int is_ = -i4_sign(p[i - 1]);
            p[i - 1] = i4_sign( is_ ) *Math.Abs(p[i - 1]);
        }

        for (i = 1; i <= n; i++)
        {
            i1 = -p[i - 1];

            switch (i1)
            {
                case >= 0:
                {
                    int i0 = i;

                    for (;;)
                    {
                        i2 = p[i1 - 1];
                        p[i1 - 1] = i0;

                        if (i2 < 0)
                        {
                            break;
                        }

                        i0 = i1;
                        i1 = i2;
                    }

                    break;
                }
            }
        }

        //
        //  Now we can restore the permutation.
        //
        for (i = 0; i < n; i++)
        {
            p[i] = p[i] + p_min - base_;
        }
    }

}