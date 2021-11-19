using System;
using System.Globalization;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.SubsetNS;

public static class Equiv
{
    public static void equiv_next(int n, ref int npart, ref int[] jarray, ref int[] iarray, ref bool more)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EQUIV_NEXT computes the partitions of a set one at a time.
        //
        //  Discussion:
        //
        //    A partition of a set assigns each element to exactly one subset.
        //
        //    The number of partitions of a set of size N is the Bell number B(N).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 May 2003
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
        //    Input, int N, number of elements in the set to be partitioned.
        //
        //    Output, int &NPART, number of subsets in the partition.
        //
        //    Output, int JARRAY[N].  JARRAY[I] is the number of elements
        //    in the I-th subset of the partition.
        //
        //    Output, int IARRAY[N].  IARRAY(I) is the class to which
        //    element I belongs.
        //
        //    Input/output, bool &MORE.  Set MORE = FALSE before first call.
        //    It is reset and held at TRUE as long as
        //    the partition returned is not the last one.
        //    When MORE is returned FALSE, all the partitions
        //    have been computed and returned.
        //
    {
        switch (more)
        {
            case false:
            {
                npart = 1;
                int i;
                for (i = 0; i < n; i++)
                {
                    iarray[i] = 1;
                }

                jarray[0] = n;
                break;
            }
            default:
            {
                int m = n;

                while (jarray[iarray[m - 1] - 1] == 1)
                {
                    iarray[m - 1] = 1;
                    m -= 1;
                }

                int l = iarray[m - 1];
                npart = npart + m - n;
                jarray[0] = jarray[0] + n - m;

                if (l == npart)
                {
                    npart += 1;
                    jarray[npart - 1] = 0;
                }

                iarray[m - 1] = l + 1;
                jarray[l - 1] -= 1;
                jarray[l] += 1;
                break;
            }
        }

        more = npart != n;
    }

    public static void equiv_next2(ref bool done, ref int[] iarray, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EQUIV_NEXT2 computes, one at a time, the partitions of a set.
        //
        //  Discussion:
        //
        //    A partition of a set assigns each element to exactly one subset.
        //
        //    The number of partitions of a set of size N is the Bell number B(N).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 June 2003
        //
        //  Author:
        //
        //    John Burkardt.
        //
        //  Parameters:
        //
        //    Input/output, bool &DONE.  Before the very first call, the
        //    user should set DONE to TRUE, which prompts the program
        //    to initialize its data, and return the first partition.
        //    Thereafter, the user should call again, for the next
        //    partition, and so on, until the routine returns with DONE
        //    equal to TRUE, at which point there are no more partitions
        //    to compute.
        //
        //    Input/output, int IARRAY[N], contains the information
        //    defining the current partition.  The user should not alter
        //    IARRAY between calls.  Except for the very first
        //    call, the routine uses the previous output value of IARRAY to compute
        //    the next value.
        //    The entries of IARRAY are the partition subset to which each
        //    element of the original set belongs.  If there are NPART distinct
        //    parts of the partition, then each entry of IARRAY will be a
        //    number between 1 and NPART.  Every number from 1 to NPART will
        //    occur somewhere in the list.  If the entries of IARRAY are
        //    examined in order, then each time a new partition subset occurs,
        //    it will be the next unused integer.
        //    For instance, for N = 4, the program will describe the set
        //    where each element is in a separate subset as 1, 2, 3, 4,
        //    even though such a partition might also be described as
        //    4, 3, 2, 1 or even 1, 5, 8, 19.
        //
        //    Input, int N, the number of elements in the set.
        //
    {
        int i;

        switch (done)
        {
            case true:
            {
                done = false;
                for (i = 0; i < n; i++)
                {
                    iarray[i] = 1;
                }

                break;
            }
            default:
            {
                //
                //  Find the last element J that can be increased by 1.
                //  This is the element that is not equal to its maximum possible value,
                //  which is the maximum value of all preceding elements +1.
                //
                int jmax = iarray[0];
                int imax = 1;

                int j;
                for (j = 2; j <= n; j++)
                {
                    if (jmax < iarray[j - 1])
                    {
                        jmax = iarray[j - 1];
                    }
                    else
                    {
                        imax = j;
                    }
                }

                switch (imax)
                {
                    //
                    //  If no element can be increased by 1, we are done.
                    //
                    case 1:
                        done = true;
                        return;
                }

                //
                //  Increase the value of the IMAX-th element by 1, set its successors to 1.
                //
                done = false;
                iarray[imax - 1] += 1;
                for (i = imax; i < n; i++)
                {
                    iarray[i] = 1;
                }

                break;
            }
        }
    }

    public static void equiv_print(int n, int[] iarray, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EQUIV_PRINT prints a partition of a set.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 July 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, number of elements in set to be partitioned.
        //
        //    Input, int IARRAY[N], defines the partition or set of equivalence
        //    classes.  Element I belongs to subset IARRAY[I].
        //
        //    Input, string TITLE, a title.
        //
    {
        int s;

        switch (title.Length)
        {
            case > 0:
                Console.WriteLine("");
                Console.WriteLine(title + "");
                break;
        }

        Console.WriteLine("  Set    Size  Elements");

        int[] karray = new int[n];

        int s_min = typeMethods.i4vec_min(n, iarray);
        int s_max = typeMethods.i4vec_max(n, iarray);

        for (s = s_min; s <= s_max; s++)
        {
            int k = 0;

            int j;
            for (j = 0; j < n; j++)
            {
                if (iarray[j] != s)
                {
                    continue;
                }

                karray[k] = j + 1;
                k += 1;
            }

            switch (k)
            {
                case > 0:
                {
                    string cout = "  "
                                  + s.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                                  + k.ToString(CultureInfo.InvariantCulture).PadLeft(4) + " :: ";
                    int kk;
                    for (kk = 0; kk < k; kk++)
                    {
                        cout += karray[kk].ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  ";
                    }

                    Console.WriteLine(cout);
                    break;
                }
            }
        }
    }

    public static void equiv_print2(int n, int[] s, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EQUIV_PRINT2 prints a partition of a set.
        //
        //  Discussion:
        //
        //    The partition is printed using the parenthesis format.
        //
        //    For example, here are the partitions of a set of 4 elements:
        //
        //      (1,2,3,4)
        //      (1,2,3)(4)
        //      (1,2,4)(3)
        //      (1,2)(3,4)
        //      (1,2)(3)(4)
        //      (1,3,4)(2)
        //      (1,3)(2,4)
        //      (1,3)(2)(4)
        //      (1,4)(2,3)
        //      (1)(2,3,4)
        //      (1)(2,3)(4)
        //      (1,4)(2)(3)
        //      (1)(2,4)(3)
        //      (1)(2)(3,4)
        //      (1)(2)(3)(4)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, number of elements in the set.
        //
        //    Input, int S[N], defines the partition.  
        //    Element I belongs to subset S[I].
        //
        //    Input, string TITLE, a title.
        //
    {
        int j;

        switch (title.Length)
        {
            case > 0:
                Console.WriteLine("");
                Console.WriteLine(title + "");
                break;
        }

        int s_min = typeMethods.i4vec_min(n, s);
        int s_max = typeMethods.i4vec_max(n, s);

        string cout = "  ";
        for (j = s_min; j <= s_max; j++)
        {
            cout += "(";
            int size = 0;
            int i;
            for (i = 0; i < n; i++)
            {
                if (s[i] != j)
                {
                    continue;
                }

                switch (size)
                {
                    case > 0:
                        cout += ",";
                        break;
                }

                cout += i;
                size += 1;
            }

            cout += ")";
        }

        Console.WriteLine(cout);
    }

    public static void equiv_random(int n, ref int seed, ref int npart, ref int[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EQUIV_RANDOM selects a random partition of a set.
        //
        //  Discussion:
        //
        //    The user does not control the number of parts in the partition.
        //
        //    The equivalence classes are numbered in no particular order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 May 2015
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
        //    Input, int N, the number of elements in the set to be partitioned.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, int &NPART, the number of classes or parts in the 
        //    partition.  NPART will be between 1 and N.
        //
        //    Output, int A[N], indicates the class to which each element
        //    is assigned.
        //
    {
        int i;
        int k;
        int l;

        double[] b = new double[n];

        b[0] = 1.0;

        for (l = 1; l <= n - 1; l++)
        {
            double sum1 = 1.0 / l;
            for (k = 1; k <= l - 1; k++)
            {
                sum1 = (sum1 + b[k - 1]) / (l - k);
            }

            b[l] = (sum1 + b[l - 1]) / (l + 1);
        }

        int m = n;
        npart = 0;

        for (;;)
        {
            double z = UniformRNG.r8_uniform_01(ref seed);
            z = m * b[m - 1] * z;
            k = 0;
            npart += 1;

            while (0.0 <= z)
            {
                a[m - 1] = npart;
                m -= 1;

                if (m == 0)
                {
                    break;
                }

                z -= b[m - 1];
                k += 1;
                z *= k;
            }

            if (m == 0)
            {
                break;
            }
        }

        //
        //  Randomly permute the assignments.
        //
        for (i = 0; i < n - 1; i++)
        {
            int j = UniformRNG.i4_uniform_ab(i, n - 1, ref seed);

            k = a[i];
            a[i] = a[j];
            a[j] = k;
        }
    }
}