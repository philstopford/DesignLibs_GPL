using System;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void binary_vector_next(int n, ref int[] bvec)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BINARY_VECTOR_NEXT generates the next binary vector.
        //
        //  Discussion:
        //
        //    A binary vector is a vector whose entries are 0 or 1.
        //
        //    The user inputs an initial zero vector to start.  The program returns
        //    the "next" vector.
        //
        //    The vectors are produced in the order:
        //
        //    ( 0, 0, 0, ..., 0 )
        //    ( 1, 0, 0, ..., 0 )
        //    ( 0, 1, 0, ..., 0 )
        //    ( 1, 1, 0, ..., 0 )
        //    ( 0, 0, 1, ..., 0 )
        //    ( 1, 0, 1, ..., 0 )
        //               ...
        //    ( 1, 1, 1, ..., 1)
        //
        //    and the "next" vector after (1,1,...,1) is (0,0,...,0).  That is,
        //    we allow wrap around.
        //
        //  Example:
        //
        //    N = 3
        //
        //    Input      Output
        //    -----      ------
        //    0 0 0  =>  1 0 0
        //    1 0 0  =>  0 1 0
        //    0 1 0  =>  1 1 0
        //    1 1 0  =>  0 0 1
        //    0 0 1  =>  1 0 1
        //    1 0 1  =>  0 1 1
        //    0 1 1  =>  1 1 1
        //    1 1 1  =>  0 0 0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 September 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the dimension of the vectors.
        //
        //    Input/output, int BVEC[N], on output, the successor
        //    to the input vector.
        //
    {
        int i;

        for (i = 0; i < n; i++)
        {
            if (bvec[i] == 1)
            {
                bvec[i] = 0;
            }
            else
            {
                bvec[i] = 1;
                break;
            }
        }
    }

    public class VecNextData
    {
        public int kount { get; set; }
        public int last { get; set; }

        public VecNextData()
        {
            kount = 0;
            last = 0;
        }
    }

    public static void vec_next(ref VecNextData data, int n, int ibase, ref int[] iarray, ref bool more)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VEC_NEXT generates all N-vectors of integers modulo a given base.
        //
        //  Discussion:
        //
        //    The items are produced one at a time.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2013
        //
        //  Parameters:
        //
        //    Input, int N, the size of the vectors to be used.
        //
        //    Input, int IBASE, the base to be used.  IBASE = 2 will
        //    give vectors of 0's and 1's, for instance.
        //
        //    Input/output, int IARRAY[N].  On each return,
        //    IARRAY will contain entries in the range 0 to IBASE-1.
        //
        //    Input/output, ref int MORE.  Set this variable 0 before
        //    the first call.  Normally, MORE will be returned 1 but
        //    once all the vectors have been generated, MORE will be
        //    reset 0 and you should stop calling the program.
        //
    {
        int i;

        switch (more)
        {
            case false:
            {
                data.kount = 1;
                data.last = (int) Math.Pow(ibase, n);
                more = true;
                for (i = 0; i < n; i++)
                {
                    iarray[i] = 0;
                }

                break;
            }
            default:
            {
                data.kount += 1;

                if (data.kount == data.last)
                {
                    more = false;
                }

                iarray[n - 1] += 1;

                for (i = 1; i <= n; i++)
                {
                    int nn = n - i;

                    if (iarray[nn] < ibase)
                    {
                        return;
                    }

                    iarray[nn] = 0;

                    if (nn != 0)
                    {
                        iarray[nn - 1] += 1;
                    }
                }

                break;
            }
        }
    }

    public static void vec_random(int n, int base_, ref int seed, int[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VEC_RANDOM selects a random N-vector of integers modulo a given base.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 March 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the size of the vector to be generated.
        //
        //    Input, int BASE, the base to be used.
        //
        //    Input/output, ref int SEED, a random number seed.
        //
        //    Output, int A[N], a list of N random values between
        //    0 and BASE-1.
        //
    {
        int i;

        for (i = 0; i < n; i++)
        {
            a[i] = UniformRNG.i4_uniform_ab(0, base_ - 1, ref seed);
        }

    }

    public static int[] vec_random_new(int n, int base_, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VEC_RANDOM_NEW selects a random N-vector of integers modulo a given base.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 March 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the size of the vector to be generated.
        //
        //    Input, int BASE, the base to be used.
        //
        //    Input/output, ref int SEED, a random number seed.
        //
        //    Output, int VEC_RANDOM_NEW[N], a list of N random values between
        //    0 and BASE-1.
        //
    {
        int i;

        int[] a = new int[n];

        for (i = 0; i < n; i++)
        {
            a[i] = UniformRNG.i4_uniform_ab(0, base_ - 1, ref seed);
        }

        return a;
    }

    public static void vec_lex_next(int dim_num, int base_, ref int[] a, ref bool more)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VEC_LEX_NEXT generates vectors in lex order.
        //
        //  Discussion:
        //
        //    The vectors are produced in lexical order, starting with
        //    (0,0,...,0),
        //    (0,0,...,1), 
        //    ...
        //    (BASE-1,BASE-1,...,BASE-1).
        //
        //  Example:
        //
        //    DIM_NUM = 2,
        //    BASE = 3
        //
        //    0   0
        //    0   1
        //    0   2
        //    1   0
        //    1   1
        //    1   2
        //    2   0
        //    2   1
        //    2   2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the size of the vectors to be used.
        //
        //    Input, int BASE, the base to be used.  BASE = 2 will
        //    give vectors of 0's and 1's, for instance.
        //
        //    Output, int A[DIM_NUM], the next vector.  
        //
        //    Input/output, bool &MORE.  Set this variable false before
        //    the first call.  On return, MORE is TRUE if another vector has
        //    been computed.  If MORE is returned FALSE, ignore the output 
        //    vector and stop calling the routine.
        //
    {
        int i;

        switch (more)
        {
            case false:
            {
                for (i = 0; i < dim_num; i++)
                {
                    a[i] = 0;
                }

                more = true;
                break;
            }
            default:
            {
                for (i = dim_num - 1; 0 <= i; i--)
                {
                    a[i] += 1;

                    if (a[i] < base_)
                    {
                        return;
                    }

                    a[i] = 0;
                }

                more = false;
                break;
            }
        }
    }

    public static void vec_colex_next(int dim_num, int base_, ref int[] a, ref bool more)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VEC_COLEX_NEXT generates vectors in colex order.
        //
        //  Discussion:
        //
        //    The vectors are produced in colexical order, starting with
        //    (0,0,...,0),
        //    (1,0,...,0), 
        //    ...
        //    (BASE-1,BASE-1,...,BASE-1).
        //
        //  Example:
        //
        //    DIM_NUM = 2,
        //    BASE = 3
        //
        //    0   0
        //    1   0
        //    2   0
        //    0   1
        //    1   1
        //    2   1
        //    0   2
        //    1   2
        //    2   2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int BASE, the base to be used.  BASE = 2 will
        //    give vectors of 0's and 1's, for instance.
        //
        //    Output, int A[DIM_NUM], the next vector.
        //
        //    Input/output, bool &MORE.  Set this variable false before
        //    the first call.  On return, MORE is TRUE if another vector has
        //    been computed.  If MORE is returned FALSE, ignore the output 
        //    vector and stop calling the routine.
        //
    {
        int i;

        switch (more)
        {
            case false:
            {
                for (i = 0; i < dim_num; i++)
                {
                    a[i] = 0;
                }

                more = true;
                break;
            }
            default:
            {
                for (i = 0; i < dim_num; i++)
                {
                    a[i] += 1;

                    if (a[i] < base_)
                    {
                        return;
                    }

                    a[i] = 0;
                }

                more = false;
                break;
            }
        }
    }

    public static void vec_colex_next2(int dim_num, int[] base_, ref int[] a, ref bool more)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VEC_COLEX_NEXT2 generates vectors in colex order.
        //
        //  Discussion:
        //
        //    The vectors are produced in colexical order, starting with
        //
        //    (0,        0,        ...,0),
        //    (1,        0,        ...,0),
        //     ...
        //    (BASE(1)-1,0,        ...,0)
        //
        //    (0,        1,        ...,0)
        //    (1,        1,        ...,0)
        //    ...
        //    (BASE(1)-1,1,        ...,0)
        //
        //    (0,        2,        ...,0)
        //    (1,        2,        ...,0)
        //    ...
        //    (BASE(1)-1,BASE(2)-1,...,BASE(DIM_NUM)-1).
        //
        //  Example:
        //
        //    DIM_NUM = 2,
        //    BASE = { 3, 3 }
        //
        //    0   0
        //    1   0
        //    2   0
        //    0   1
        //    1   1
        //    2   1
        //    0   2
        //    1   2
        //    2   2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int BASE[DIM_NUM], the bases to be used in each dimension.
        //    In dimension I, entries will range from 0 to BASE[I]-1.
        //
        //    Output, int A[DIM_NUM], the next vector.
        //
        //    Input/output, bool &MORE.  Set this variable false before
        //    the first call.  On return, MORE is TRUE if another vector has
        //    been computed.  If MORE is returned FALSE, ignore the output 
        //    vector and stop calling the routine.
        //
    {
        int i;

        switch (more)
        {
            case false:
            {
                for (i = 0; i < dim_num; i++)
                {
                    a[i] = 0;
                }

                more = true;
                break;
            }
            default:
            {
                for (i = 0; i < dim_num; i++)
                {
                    a[i] += 1;

                    if (a[i] < base_[i])
                    {
                        return;
                    }

                    a[i] = 0;
                }

                more = false;
                break;
            }
        }
    }

    public static void vec_colex_next3(int dim_num, int[] base_, ref int[] a, ref bool more)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VEC_COLEX_NEXT3 generates vectors in colex order.
        //
        //  Discussion:
        //
        //    The vectors are produced in colexical order, starting with
        //
        //    (1,        1,        ...,1),
        //    (2,        1,        ...,1),
        //     ...
        //    (BASE(1),  1,        ...,1)
        //
        //    (1,        2,        ...,1)
        //    (2,        2,        ...,1)
        //    ...
        //    (BASE(1),  2,        ...,1)
        //
        //    (1,        3,        ...,1)
        //    (2,        3,        ...,1)
        //    ...
        //    (BASE(1),  BASE(2), ...,BASE(DIM_NUM)).
        //
        //  Example:
        //
        //    DIM_NUM = 2,
        //    BASE = { 3, 3 }
        //
        //    1   1
        //    2   1
        //    3   1
        //    1   2
        //    2   2
        //    3   2
        //    1   3
        //    2   3
        //    3   3
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 August 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int BASE[DIM_NUM], the bases to be used in each dimension.
        //    In dimension I, entries will range from 1 to BASE[I].
        //
        //    Output, int A[DIM_NUM], the next vector.
        //
        //    Input/output, bool &MORE.  Set this variable false before
        //    the first call.  On return, MORE is TRUE if another vector has
        //    been computed.  If MORE is returned FALSE, ignore the output 
        //    vector and stop calling the routine.
        //
    {
        int i;

        switch (more)
        {
            case false:
            {
                for (i = 0; i < dim_num; i++)
                {
                    a[i] = 1;
                }

                more = true;
                break;
            }
            default:
            {
                for (i = 0; i < dim_num; i++)
                {
                    a[i] += 1;

                    if (a[i] <= base_[i])
                    {
                        return;
                    }

                    a[i] = 1;
                }

                more = false;
                break;
            }
        }
    }


    public class VecGrayData
    {
        public int[] active;
        public int[] dir;
    }

    public static void vec_next_gray(ref VecGrayData data, int n, int[] base_, ref int[] a, ref bool done,
            ref int change)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VEC_NEXT_GRAY computes the elements of a product space.
        //
        //  Discussion:
        //
        //    The elements are produced one at a time.
        //
        //    This routine handles the case where the number of degrees of freedom may
        //    differ from one component to the next.
        //
        //    A method similar to the Gray code is used, so that successive
        //    elements returned by this routine differ by only a single element.
        //
        //    The routine uses internal static memory.
        //
        //  Example:
        //
        //    N = 2, BASE = ( 2, 3 ), DONE = TRUE
        //
        //     A    DONE  CHANGE
        //    ---  -----  ------
        //    0 0  FALSE    1
        //    0 1  FALSE    2
        //    0 2  FALSE    2
        //    1 2  FALSE    1
        //    1 1  FALSE    2
        //    1 0  FALSE    2
        //    1 0   TRUE   -1  
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 October 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Dennis Stanton, Dennis White,
        //    Constructive Combinatorics,
        //    Springer, 1986,
        //    ISBN: 0387963472.
        //
        //  Parameters:
        //
        //    Input, int N, the number of components.
        //
        //    Input, int BASE[N], contains the number of degrees of
        //    freedom of each component.  The output values of A will
        //    satisfy 0 <= A[I] < BASE[I].
        //
        //    Input/output, int A[N].  On the first call, the input value
        //    of A doesn't matter.  Thereafter, it should be the same as
        //    its output value from the previous call.  On output, if DONE
        //    is FALSE, then A contains the next element of the space.
        //
        //    Input/output, bool *DONE.  On the first call, the user must
        //    set DONE to TRUE.  This signals the program to initialize data.
        //    On every return, if DONE is FALSE, the program has computed
        //    another entry, which is contained in A.  If DONE is TRUE,
        //    then there are no more entries, and the program should not be
        //    called for any more.
        //
        //    Output, int *CHANGE, is set to the index of the element whose
        //    value was changed.  On return from the first call, CHANGE
        //    is 1, even though all the elements have been "changed".  On
        //    return with DONE equal to TRUE, CHANGE is -1.
        //
    {
        int i;
        switch (done)
        {
            //
            //  The user is calling for the first time.
            //
            case true:
            {
                done = false;
                for (i = 0; i < n; i++)
                {
                    a[i] = 0;
                }

                data.active = new int[n];
                data.dir = new int[n];

                for (i = 0; i < n; i++)
                {
                    data.dir[i] = 1;
                }

                for (i = 0; i < n; i++)
                {
                    data.active[i] = 1;
                }

                for (i = 0; i < n; i++)
                {
                    switch (base_[i])
                    {
                        case < 1:
                            Console.WriteLine("");
                            Console.WriteLine("VEC_NEXT_GRAY - Warning");
                            Console.WriteLine("  For index I = " + i + "");
                            Console.WriteLine("  the nonpositive value of BASE[I] = " + base_[i] + "");
                            Console.WriteLine("  which was reset to 1!");
                            base_[i] = 1;
                            data.active[i] = 0;
                            break;
                        case 1:
                            data.active[i] = 0;
                            break;
                    }
                }

                change = 0;
                return;
            }
        }

        //
        //  Find the maximum active index.
        //
        change = -1;

        for (i = 0; i < n; i++)
        {
            if (data.active[i] != 0)
            {
                change = i;
            }
        }

        switch (change)
        {
            //
            //  If there are NO active indices, we have generated all vectors.
            //
            case -1:
                data.dir = null;
                data.active = null;
                done = true;
                return;
        }

        //
        //  Increment the element with maximum active index.
        //
        a[change] += data.dir[change];
        //
        //  If we attained a minimum or maximum value, reverse the direction
        //  vector, and deactivate the index.
        //
        if (a[change] == 0 || a[change] == base_[change] - 1)
        {
            data.dir[change] = -data.dir[change];
            data.active[change] = 0;
        }

        //
        //  Activate all subsequent indices.
        //
        for (i = change + 1; i < n; i++)
        {
            data.active[i] = base_[i] switch
            {
                > 1 => 1,
                _ => data.active[i]
            };
        }
    }

    public static void vec_gray_next(int n, int[] base_, int[] a, ref bool done, int[] active,
            int[] dir, ref int change)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VEC_GRAY_NEXT computes the next vector in Gray code.
        //
        //  Discussion:
        //
        //    The elements are produced one at a time.
        //
        //    This routine handles the case where the number of degrees of freedom may
        //    differ from one component to the next.
        //
        //    A method similar to the Gray code is used, so that successive
        //    elements returned by this routine differ by only a single element.
        //
        //    A previous version of this routine used internal static memory.
        //
        //  Example:
        //
        //    N = 2, BASE = ( 2, 3 ), DONE = TRUE
        //
        //     A    DONE  CHANGE
        //    ---  -----  ------
        //    0 0  FALSE    0 (1)  <-- C++ routine returns 0-based CHANGE.
        //    0 1  FALSE    1 (2)
        //    0 2  FALSE    1 (2)
        //    1 2  FALSE    0 (1)
        //    1 1  FALSE    1 (2)
        //    1 0  FALSE    1 (2)
        //    1 0   TRUE   -1  
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Dennis Stanton, Dennis White,
        //    Constructive Combinatorics,
        //    Springer, 1986,
        //    ISBN: 0387963472,
        //    LC: QA164.S79.
        //
        //  Parameters:
        //
        //    Input, int N, the number of components.
        //
        //    Input, int BASE[N], contains the number of degrees of
        //    freedom of each component.  The output values of A will
        //    satisfy 0 <= A[I] < BASE[I].
        //
        //    Input/output, int A[N].  On the first call, the input value
        //    of A doesn't matter.  Thereafter, it should be the same as
        //    its output value from the previous call.  On output, if DONE
        //    is FALSE, then A contains the next element of the space.
        //
        //    Input/output, bool &DONE.  On the first call, the user must
        //    set DONE to TRUE.  This signals the program to initialize data.
        //    On every return, if DONE is FALSE, the program has computed
        //    another entry, which is contained in A.  If DONE is TRUE,
        //    then there are no more entries, and the program should not be
        //    called for any more.
        //
        //    Input/output, int ACTIVE[N], DIR[N], two work arrays.  The user must
        //    allocate memory for them before the first call.  Their values should
        //    not be altered by the user.
        //
        //    Output, int &CHANGE, is set to the index of the element whose
        //    value was changed.  On return from the first call, CHANGE
        //    is 0, even though all the elements have been "changed".  On
        //    return with DONE equal to TRUE, CHANGE is -1.  (Note that CHANGE
        //    is a vector index.  In this C++ version, it is zero-based.)
        //
    {
        int i;
        switch (done)
        {
            //
            //  The user is calling for the first time.
            //
            case true:
            {
                done = false;
                for (i = 0; i < n; i++)
                {
                    a[i] = 0;
                }

                for (i = 0; i < n; i++)
                {
                    dir[i] = 1;
                }

                for (i = 0; i < n; i++)
                {
                    active[i] = 1;
                }

                for (i = 0; i < n; i++)
                {
                    switch (base_[i])
                    {
                        case < 1:
                            Console.WriteLine("");
                            Console.WriteLine("VEC_GRAY_NEXT - Warning");
                            Console.WriteLine("  For index I = " + i + "");
                            Console.WriteLine("  the nonpositive value of BASE[I] = " + base_[i] + "");
                            Console.WriteLine("  which was reset to 1!");
                            base_[i] = 1;
                            active[i] = 0;
                            break;
                        case 1:
                            active[i] = 0;
                            break;
                    }
                }

                change = 0;
                return;
            }
        }

        //
        //  Find the maximum active index.
        //
        change = -1;

        for (i = 0; i < n; i++)
        {
            if (active[i] != 0)
            {
                change = i;
            }
        }

        switch (change)
        {
            //
            //  If there are NO active indices, we have generated all vectors.
            //
            case -1:
                done = true;
                return;
        }

        //
        //  Increment the element with maximum active index.
        //
        a[change] += dir[change];
        //
        //  If we attained a minimum or maximum value, reverse the direction
        //  vector, and deactivate the index.
        //
        if (a[change] == 0 || a[change] == base_[change] - 1)
        {
            dir[change] = -dir[change];
            active[change] = 0;
        }

        //
        //  Activate all subsequent indices.
        //
        for (i = change + 1; i < n; i++)
        {
            active[i] = base_[i] switch
            {
                > 1 => 1,
                _ => active[i]
            };
        }
    }

    public static void vec_random(int n, int base_, ref int seed, ref int[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VEC_RANDOM selects a random N-vector of integers modulo a given base.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the size of the vector to be generated.
        //
        //    Input, int BASE, the base to be used.
        //
        //    Input/output, int &SEED, a random number seed.
        //
        //    Output, int A[N], a list of N random values between
        //    0 and BASE-1.
        //
    {
        int i;

        for (i = 0; i < n; i++)
        {
            a[i] = UniformRNG.i4_uniform_ab(0, base_ - 1, ref seed);
        }
    }

    public static void vector_constrained_next(int n, int[] x_min, int[] x_max, ref int[] x,
            ref int constraint, ref bool more)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VECTOR_CONSTRAINED_NEXT returns the "next" constrained vector.
        //
        //  Discussion:
        //
        //    We consider all vectors of dimension N whose components
        //    satisfy X_MIN(1:N) <= X(1:N) <= X_MAX(1:N).
        //
        //    We are only interested in the subset of these vectors which
        //    satisfy the following constraint:
        //
        //      sum ( 1 <= I <= N ) ( ( X(I) - 1 ) / X_MAX(I) ) <= 1
        //
        //    We can carry out this check using integer arithmetic if we
        //    multiply through by P = product ( X_MAX(1:N) ):
        //
        //      sum ( 1 <= I <= N ) ( ( X(I) - 1 ) * ( P / X_MAX(I) ) ) <= P.
        //
        //    This routine returns, one at a time, and in right-to-left
        //    lexicographic order, exactly those vectors which satisfy
        //    the constraint.
        //
        //  Example:
        //
        //    N = 3
        //    X_MIN:   2   2   1
        //    X_MAX:   4   5   3
        //
        //    P = 60
        //
        //    #  X(1)  X(2)  X(3)  CONSTRAINT
        //
        //    1    2     2     1       27
        //    2    3     2     1       42
        //    3    4     2     1       57
        //    4    2     3     1       39
        //    5    3     3     1       54
        //    6    2     4     1       51
        //    7    2     2     2       47
        //    8    2     3     2       59
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of components in the vector.
        //
        //    Input, int X_MIN[N], X_MAX[N], the minimum and maximum
        //    values allowed in each component.
        //
        //    Input/output, int X[N].  On first call (with MORE = FALSE),
        //    the input value of X is not important.  On subsequent calls, the
        //    input value of X should be the output value from the previous call.
        //    On output, (with MORE = TRUE), the value of X will be the "next"
        //    vector in the reverse lexicographical list of vectors that satisfy
        //    the condition.  However, on output with MORE = FALSE, the vector
        //    X is meaningless, because there are no more vectors in the list.
        //
        //    Output, int &CONSTRAINT, the constraint value for X.  Valid vectors X
        //    will have a CONSTRAINT value between product(X_MIN(1:N)) (automatically)
        //    and product(X_MAX(1:N)) (because we skip over vectors with a
        //    constraint larger than this value).
        //
        //    Input/output, bool &MORE.  On input, if the user has set MORE
        //    FALSE, the user is requesting the initiation of a new sequence
        //    of values.  If MORE is TRUE, then the user is requesting "more"
        //    values in the current sequence.  On output, if MORE is TRUE,
        //    then another value was found and returned in X, but if MORE is
        //    FALSE, then there are no more values in the sequence, and X is
        //    NOT the next value.
        //
    {
        int j;

        int x_prod = 1;
        for (j = 0; j < n; j++)
        {
            x_prod *= x_max[j];
        }

        switch (more)
        {
            case false:
            {
                for (j = 0; j < n; j++)
                {
                    x[j] = x_min[j];
                }

                constraint = 0;
                for (j = 0; j < n; j++)
                {
                    constraint += (x[j] - 1) * (x_prod / x_max[j]);
                }

                more = x_prod >= constraint;

                break;
            }
            default:
            {
                int i = 0;

                for (;;)
                {
                    if (x[i] < x_max[i])
                    {
                        x[i] += 1;

                        constraint = 0;
                        for (j = 0; j < n; j++)
                        {
                            constraint += (x[j] - 1) * (x_prod / x_max[j]);
                        }

                        if (constraint <= x_prod)
                        {
                            break;
                        }
                    }

                    x[i] = x_min[i];

                    i += 1;

                    if (n > i)
                    {
                        continue;
                    }

                    more = false;
                    break;
                }

                break;
            }
        }
    }

    public static void vector_constrained_next2(int n, int[] x_min, int[] x_max, ref int[] x,
            ref int constraint, ref bool more)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VECTOR_CONSTRAINED_NEXT2 returns the "next" constrained vector.
        //
        //  Discussion:
        //
        //    We consider all vectors of dimension N whose components
        //    satisfy X_MIN(1:N) <= X(1:N) <= X_MAX(1:N).
        //
        //    We are only interested in the subset of these vectors which
        //    satisfy the following constraint:
        //
        //      sum ( 1 <= I <= N ) ( X(I) / X_MAX(I) ) <= 1
        //
        //    We can carry out this check using integer arithmetic if we
        //    multiply through by P = product ( X_MAX(1:N) ):
        //
        //      sum ( 1 <= I <= N ) ( X(I)  * ( P / X_MAX(I) ) ) <= P.
        //
        //    This routine returns, one at a time, and in right-to-left
        //    lexicographic order, exactly those vectors which satisfy
        //    the constraint.
        //
        //  Example:
        //
        //    N = 3
        //    X_MIN:   1   1   1
        //    X_MAX:   5   6   4
        //
        //    P = 120
        //
        //    #  X(1)  X(2)  X(3)  CONSTRAINT
        //
        //    1    1     1     1       74
        //    2    2     1     1       98
        //    3    1     2     1       94
        //    4    2     2     1      119
        //    5    1     3     1      114
        //    6    1     1     2      104
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of components in the vector.
        //
        //    Input, int X_MIN[N], X_MAX[N], the minimum and maximum
        //    values allowed in each component.
        //
        //    Input/output, int X[N].  On first call (with MORE = FALSE),
        //    the input value of X is not important.  On subsequent calls, the
        //    input value of X should be the output value from the previous call.
        //    On output, (with MORE = TRUE), the value of X will be the "next"
        //    vector in the reverse lexicographical list of vectors that satisfy
        //    the condition.  However, on output with MORE = FALSE, the vector
        //    X is meaningless, because there are no more vectors in the list.
        //
        //    Output, int &CONSTRAINT, the constraint value for X.  Valid vectors X
        //    will have a CONSTRAINT value between product(X_MIN(1:N)) (automatically)
        //    and product(X_MAX(1:N)) (because we skip over vectors with a
        //    constraint larger than this value).
        //
        //    Input/output, bool &MORE.  On input, if the user has set MORE
        //    FALSE, the user is requesting the initiation of a new sequence
        //    of values.  If MORE is TRUE, then the user is requesting "more"
        //    values in the current sequence.  On output, if MORE is TRUE,
        //    then another value was found and returned in X, but if MORE is
        //    FALSE, then there are no more values in the sequence, and X is
        //    NOT the next value.
        //
    {
        int j;

        int x_prod = 1;
        for (j = 0; j < n; j++)
        {
            x_prod *= x_max[j];
        }

        switch (more)
        {
            case false:
            {
                for (j = 0; j < n; j++)
                {
                    x[j] = x_min[j];
                }

                constraint = 0;
                for (j = 0; j < n; j++)
                {
                    constraint += x[j] * (x_prod / x_max[j]);
                }

                more = x_prod >= constraint;

                break;
            }
            default:
            {
                int i = 0;

                for (;;)
                {
                    if (x[i] < x_max[i])
                    {
                        x[i] += 1;

                        constraint = 0;
                        for (j = 0; j < n; j++)
                        {
                            constraint += x[j] * (x_prod / x_max[j]);
                        }

                        if (constraint <= x_prod)
                        {
                            break;
                        }
                    }

                    x[i] = x_min[i];

                    i += 1;

                    if (n > i)
                    {
                        continue;
                    }

                    more = false;
                    break;
                }

                break;
            }
        }

    }

    public static void vector_constrained_next3(int n, int[] x_min, int[] x_max, ref int[] x,
            ref double constraint, ref bool more)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VECTOR_CONSTRAINED_NEXT3 returns the "next" constrained vector.
        //
        //  Discussion:
        //
        //    This routine addresses the same problem as VECTOR_CONSTRAINED_NEXT2,
        //    and differs only in that real arithmetic is used, rather than
        //    integer arithmetic.  Integer arithmetic allows us to do an exact
        //    calculation, but we run into overflow problems in simple cases
        //    where N is 10 and the X_MAX entries are of order 10, for instance.
        //
        //    We consider all vectors of dimension N whose components
        //    satisfy X_MIN(1:N) <= X(1:N) <= X_MAX(1:N).
        //
        //    We are only interested in the subset of these vectors which
        //    satisfy the following constraint:
        //
        //      sum ( 1 <= I <= N ) ( X(I) / X_MAX(I) ) <= 1
        //
        //  Example:
        //
        //    N = 3
        //    X_MIN:   1   1   1
        //    X_MAX:   5   6   4
        //
        //    P = 120
        //
        //    #  X(1)  X(2)  X(3)  CONSTRAINT
        //
        //    1    1     1     1       0.62
        //    2    2     1     1       0.82
        //    3    1     2     1       0.78
        //    4    2     2     1       0.98
        //    5    1     3     1       0.95
        //    6    1     1     2       0.87
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 April 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of components in the vector.
        //
        //    Input, int X_MIN[N], X_MAX[N], the minimum and maximum
        //    values allowed in each component.
        //
        //    Input/output, int X[N].  On first call (with MORE = FALSE),
        //    the input value of X is not important.  On subsequent calls, the
        //    input value of X should be the output value from the previous call.
        //    On output, (with MORE = TRUE), the value of X will be the "next"
        //    vector in the reverse lexicographical list of vectors that satisfy
        //    the condition.  However, on output with MORE = FALSE, the vector
        //    X is meaningless, because there are no more vectors in the list.
        //
        //    Output, double &CONSTRAINT, the constraint value for X.  Valid vectors 
        //    X will have a CONSTRAINT value between 
        //      product(X_MIN(1:N)) / product(X_MAX(1:N))
        //    and 1.0.
        //
        //    Input/output, bool &MORE.  On input, if the user has set MORE
        //    FALSE, the user is requesting the initiation of a new sequence
        //    of values.  If MORE is TRUE, then the user is requesting "more"
        //    values in the current sequence.  On output, if MORE is TRUE,
        //    then another value was found and returned in X, but if MORE is
        //    FALSE, then there are no more values in the sequence, and X is
        //    NOT the next value.
        //
    {
        int j;

        switch (more)
        {
            case false:
            {
                for (j = 0; j < n; j++)
                {
                    x[j] = x_min[j];
                }

                constraint = 0.0;
                for (j = 0; j < n; j++)
                {
                    constraint += x[j] / (double) x_max[j];
                }

                more = constraint switch
                {
                    > 1.0 => false,
                    _ => true
                };

                break;
            }
            default:
            {
                int i = 0;

                for (;;)
                {
                    if (x[i] < x_max[i])
                    {
                        x[i] += 1;

                        constraint = 0;
                        for (j = 0; j < n; j++)
                        {
                            constraint += x[j] / (double) x_max[j];
                        }

                        if (constraint <= 1.0)
                        {
                            break;
                        }
                    }

                    x[i] = x_min[i];

                    i += 1;

                    if (n > i)
                    {
                        continue;
                    }

                    more = false;
                    break;
                }

                break;
            }
        }

    }

    public static void vector_constrained_next4(int n, double[] alpha, int[] x_min,
            int[] x_max, ref int[] x, double q, ref bool more)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VECTOR_CONSTRAINED_NEXT4 returns the "next" constrained vector.
        //
        //  Discussion:
        //
        //    This routine is similar to VECTOR_CONSTRAINED2 and VECTOR_CONSTRAINED3.
        //
        //    We consider all vectors X of dimension N whose components
        //    satisfy X_MIN(1:N) <= X(1:N) <= X_MAX(1:N).
        //
        //    We are only interested in the subset of these vectors which
        //    satisfy the following constraint:
        //
        //      sum ( 1 <= I <= N ) ( ALPHA(I) * X(I) ) <= Q
        //
        //  Example:
        //
        //    N = 3
        //    ALPHA    4.0  3.0  5.0
        //    Q       20.0
        //    X_MIN:   1   1   1
        //    X_MAX:   5   6   4
        //
        //    P = 120
        //
        //    #  X(1)  X(2)  X(3)      Total
        //
        //    1    1     1     1       12.0
        //    2    2     1     1       20.0
        //    3    1     2     1       15.0
        //    4    2     2     1       19.0
        //    5    1     3     1       18.0
        //    6    1     1     2       17.0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 May 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of components in the vector.
        //
        //    Input, double ALPHA[N], the coefficient vector.
        //
        //    Input, int X_MIN[N], X_MAX[N], the minimum and maximum
        //    values allowed in each component.
        //
        //    Input/output, integer X[N].  On first call (with MORE = FALSE),
        //    the input value of X is not important.  On subsequent calls, the
        //    input value of X should be the output value from the previous call.
        //    On output, (with MORE = TRUE), the value of X will be the "next"
        //    vector in the reverse lexicographical list of vectors that satisfy
        //    the condition.  However, on output with MORE = FALSE, the vector
        //    X is meaningless, because there are no more vectors in the list.
        //
        //    Input, double Q, the limit on the sum.
        //
        //    Input/output, bool *MORE.  On input, if the user has set MORE
        //    FALSE, the user is requesting the initiation of a new sequence
        //    of values.  If MORE is TRUE, then the user is requesting "more"
        //    values in the current sequence.  On output, if MORE is TRUE,
        //    then another value was found and returned in X, but if MORE is
        //    FALSE, then there are no more values in the sequence, and X is
        //    NOT the next value.
        //
    {
        int j;
        double total;

        switch (more)
        {
            case false:
            {
                for (j = 0; j < n; j++)
                {
                    x[j] = x_min[j];
                }

                total = 0.0;
                for (j = 0; j < n; j++)
                {
                    total += alpha[j] * x[j];
                }

                more = !(q < total);

                break;
            }
            default:
            {
                int i = 0;

                for (;;)
                {
                    if (x[i] < x_max[i])
                    {
                        x[i] += 1;

                        total = 0;
                        for (j = 0; j < n; j++)
                        {
                            total += alpha[j] * x[j];
                        }

                        if (total <= q)
                        {
                            break;
                        }
                    }

                    x[i] = x_min[i];

                    i += 1;

                    if (n > i)
                    {
                        continue;
                    }

                    more = false;
                    break;
                }

                break;
            }
        }
    }

    public static void vector_constrained_next5(int n, ref int[] x, int sum_min, int sum_max,
            ref int base_, ref bool more)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VECTOR_CONSTRAINED_NEXT5 returns the "next" constrained vector.
        //
        //  Discussion:
        //
        //    We consider all positive integer vectors of dimension N whose 
        //    components satisfy SUM_MIN <= X(1:N) <= SUM_MAX.
        //
        //    This routine returns, one at a time, and in right-to-left
        //    lexicographic order, exactly those vectors which satisfy
        //    the constraint.
        //
        //  Example:
        //
        //    N = 3
        //    SUM_MIN = 5
        //    SUM_MAX = 6
        //
        //    #  X(1)  X(2)  X(3)     SUM
        //
        //    1    3     1     1        5
        //    2    2     2     1        5
        //    3    2     1     2        5
        //    4    1     3     1        5
        //    5    1     2     2        5
        //    6    1     1     3        5
        //
        //    7    4     1     1        6
        //    8    3     2     1        6
        //    9    3     1     2        6
        //   10    2     3     1        6
        //   11    2     2     2        6
        //   12    2     1     3        6
        //   13    1     4     1        6
        //   14    1     3     2        6
        //   15    1     2     3        6
        //   16    1     1     4        6
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of components in the vector.
        //
        //    Input, int SUM_MIN, SUM_MAX, the minimum and maximum sums..
        //
        //    Input/output, int X(N).  On first call (with MORE = FALSE), 
        //    the input value of X is not important.  On subsequent calls, the
        //    input value of X should be the output value from the previous call.
        //    On output, (with MORE = TRUE), the value of X will be the "next"
        //    vector in the reverse lexicographical list of vectors that satisfy
        //    the condition.  However, on output with MORE = FALSE, the vector
        //    X is meaningless, because there are no more vectors in the list.
        //
        //    Input/output, int &BASE, a value that is controlled by this function.
        //    The user must declare this variable, and pass the output value from
        //    one call as the input value to the next call, but should not
        //    change the value.
        //
        //    Input/output, logical &MORE.  On input, if the user has set MORE
        //    FALSE, the user is requesting the initiation of a new sequence
        //    of values.  If MORE is TRUE, then the user is requesting "more"
        //    values in the current sequence.  On output, if MORE is TRUE,
        //    then another value was found and returned in X, but if MORE is
        //    FALSE, then there are no more values in the sequence, and X is
        //    NOT the next value.
        //
    {
        int i;
        switch (more)
        {
            //
            //  Initialization.
            //
            case false when sum_max < n:
                more = false;
                return;
            case false when sum_max < sum_min:
                more = false;
                return;
            case false:
            {
                more = true;

                base_ = sum_min;
                if (base_ < n)
                {
                    base_ = n;
                }

                x[0] = base_ - n + 1;
                for (i = 1; i < n; i++)
                {
                    x[i] = 1;
                }

                break;
            }
            //
            default:
            {
                //
                //  Search from the right, seeking an index I < N for which 1 < X(I).
                //
                for (i = n - 2; 0 <= i; i--)
                {
                    switch (x[i])
                    {
                        //
                        //  If you find such an I, decrease X(I) by 1, and add that to X(I+1).
                        //
                        case > 1:
                        {
                            x[i] -= 1;
                            x[i + 1] += 1;
                            //
                            //  Now grab all the "excess" 1's from the entries to the right of X(I+1).
                            //
                            int j;
                            for (j = i + 2; j < n; j++)
                            {
                                switch (x[j])
                                {
                                    case > 1:
                                        x[i + 1] = x[i + 1] + x[j] - 1;
                                        x[j] = 1;
                                        break;
                                }
                            }

                            return;
                        }
                    }
                }

                //
                //  The current vector is (1,1,1,...BASE-N+1).
                //  If BASE < SUM_MAX, then increase BASE by 1, and start the new series.
                //
                if (base_ < sum_max)
                {
                    base_ += 1;
                    x[0] = base_ - n + 1;
                    for (i = 1; i < n; i++)
                    {
                        x[i] = 1;
                    }

                    return;
                }

                //
                //  We returned the last legal vector on the previouis call.
                //  The calculation is done.
                //
                more = false;
                break;
            }
        }

    }

    public static void vector_constrained_next6(int n, double[] alpha, int[] x_min,
            int[] x_max, ref int[] x, double q_min, double q_max, ref bool more)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VECTOR_CONSTRAINED_NEXT6 returns the "next" constrained vector.
        //
        //  Discussion:
        //
        //    This routine is similar to VECTOR_CONSTRAINED_NEXT2,
        //    VECTOR_CONSTRAINED_NEXT3, and VECTOR_CONSTRAINED_NEXT4.
        //
        //    We consider all vectors X of dimension N whose components
        //    satisfy X_MIN(1:N) <= X(1:N) <= X_MAX(1:N).
        //
        //    We are only interested in the subset of these vectors which
        //    satisfy the following constraint:
        //
        //      Q_MIN <= sum ( 1 <= I <= N ) ALPHA(I) * X(I) <= Q_MAX
        //
        //    This routine returns, one at a time, and in right-to-left
        //    lexicographic order, exactly those vectors which satisfy
        //    the constraint.
        //
        //  Example:
        //
        //    N = 3
        //    ALPHA    4.0  3.0  5.0
        //    Q_MIN   16.0
        //    Q_MAX   20.0
        //    X_MIN:   1   1   1
        //    X_MAX:   5   6   4
        //
        //    #  X(1)  X(2)  X(3)     Total
        //
        //    1    2     1     1       20.0
        //    2    2     2     1       19.0
        //    3    1     3     1       18.0
        //    4    1     1     2       17.0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 March 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of components in the vector.
        //
        //    Input, double ALPHA[N], the coefficient vector.
        //
        //    Input, int X_MIN[N], X_MAX[N], the minimum and maximum
        //    values allowed in each component.
        //
        //    Input/output, int X[N].  On first call (with MORE = FALSE),
        //    the input value of X is not important.  On subsequent calls, the
        //    input value of X should be the output value from the previous call.
        //    On output, (with MORE = TRUE), the value of X will be the "next"
        //    vector in the reverse lexicographical list of vectors that satisfy
        //    the condition.  However, on output with MORE = FALSE, the vector
        //    X is meaningless, because there are no more vectors in the list.
        //
        //    Input, double Q_MIN, Q_MAX, the lower and upper
        //    limits on the sum.
        //
        //    Input/output, bool &MORE.  On input, if the user has set MORE
        //    FALSE, the user is requesting the initiation of a new sequence
        //    of values.  If MORE is TRUE, then the user is requesting "more"
        //    values in the current sequence.  On output, if MORE is TRUE,
        //    then another value was found and returned in X, but if MORE is
        //    FALSE, then there are no more values in the sequence, and X is
        //    NOT the next value.
        //
    {
        int i;
        double total;

        switch (more)
        {
            case false:
            {
                more = true;
                for (i = 0; i < n; i++)
                {
                    x[i] = x_min[i];
                }

                total = 0.0;
                for (i = 0; i < n; i++)
                {
                    total += alpha[i] * x[i];
                }

                if (q_min <= total && total <= q_max)
                {
                    return;
                }

                break;
            }
        }

        for (;;)
        {
            int j = n - 1;

            for (;;)
            {
                if (x[j] < x_max[j])
                {
                    break;
                }

                switch (j)
                {
                    case <= 0:
                        more = false;
                        return;
                    default:
                        j -= 1;
                        break;
                }
            }

            x[j] += 1;
            for (i = j + 1; i < n; i++)
            {
                x[i] = x_min[i];
            }

            total = 0.0;
            for (i = 0; i < n; i++)
            {
                total += alpha[i] * x[i];
            }

            if (q_min <= total && total <= q_max)
            {
                break;
            }
        }

    }

    public static void vector_constrained_next7(int n, double[] level_weight, int[] x_max,
            ref int[] x, double q_min, double q_max, ref bool more)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VECTOR_CONSTRAINED_NEXT7 returns the "next" constrained vector.
        //
        //  Discussion:
        //
        //    We consider vectors X of dimension N satisfying:
        //
        //      0 <= X(1:N) <= X_MAX(1:N).
        //
        //    We are only interested in the subset of these vectors which
        //    satisfy the following constraint:
        //
        //      Q_MIN < sum ( 1 <= I <= N ) LEVEL_WEIGHT(I) * X(I) <= Q_MAX
        //
        //    This routine returns, one at a time, and in right-to-left
        //    lexicographic order, exactly those vectors which satisfy
        //    the constraint.
        //
        //  Example:
        //
        //    N = 3
        //    LEVEL_WEIGHT    4.0  3.0  5.0
        //    Q_MIN   16.0
        //    Q_MAX   20.0
        //    X_MAX:   5   6   4
        //
        //    #  X(1)  X(2)  X(3)     Total
        //
        //    1    2     1     1       20.0
        //    2    2     2     1       19.0
        //    3    1     3     1       18.0
        //    4    1     1     2       17.0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of components in the vector.
        //
        //    Input, double LEVEL_WEIGHT[N], the coefficient vector.
        //
        //    Input, int X_MAX[N], the maximum values allowed in each component.
        //
        //    Input/output, int X[N].  On first call (with MORE = FALSE),
        //    the input value of X is not important.  On subsequent calls, the
        //    input value of X should be the output value from the previous call.
        //    On output, (with MORE = TRUE), the value of X will be the "next"
        //    vector in the reverse lexicographical list of vectors that satisfy
        //    the condition.  However, on output with MORE = FALSE, the vector
        //    X is meaningless, because there are no more vectors in the list.
        //
        //    Input, double Q_MIN, Q_MAX, the lower and upper
        //    limits on the sum.
        //
        //    Input/output, bool &MORE.  On input, if the user has set MORE
        //    FALSE, the user is requesting the initiation of a new sequence
        //    of values.  If MORE is TRUE, then the user is requesting "more"
        //    values in the current sequence.  On output, if MORE is TRUE,
        //    then another value was found and returned in X, but if MORE is
        //    FALSE, then there are no more values in the sequence, and X is
        //    NOT the next value.
        //
    {
        int i;
        double total;

        switch (more)
        {
            case false:
            {
                more = true;
                for (i = 0; i < n; i++)
                {
                    x[i] = 0;
                }

                total = 0.0;
                for (i = 0; i < n; i++)
                {
                    total += level_weight[i] * x[i];
                }

                if (q_min < total && total <= q_max)
                {
                    return;
                }

                break;
            }
        }

        for (;;)
        {
            int j = n - 1;

            for (;;)
            {
                if (x[j] < x_max[j])
                {
                    break;
                }

                switch (j)
                {
                    case <= 0:
                        more = false;
                        return;
                    default:
                        j -= 1;
                        break;
                }
            }

            x[j] += 1;
            for (i = j + 1; i < n; i++)
            {
                x[i] = 0;
            }

            total = 0.0;
            for (i = 0; i < n; i++)
            {
                total += level_weight[i] * x[i];
            }

            if (q_min < total && total <= q_max)
            {
                break;
            }
        }

    }

    public static void vector_next(int n, int[] x_min, int[] x_max, ref int[] x, ref bool more)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VECTOR_NEXT returns the "next" integer vector between two ranges.
        //
        //  Discussion:
        //
        //    We consider all integer vectors of dimension N satisfying:
        //
        //      X_MIN(1:N) <= X(1:N) <= X_MAX(1:N).
        //
        //    This routine returns, one at a time, and in right-to-left
        //    lexicographic order, all these vectors.
        //
        //  Example:
        //
        //    N = 3
        //    X_MIN:   2   2   0
        //    X_MAX:   4   3   1
        // 
        //    #  X(1)  X(2)  X(3)
        //
        //    1    2     2     0
        //    2    3     2     0
        //    3    4     2     0
        //    4    2     3     0
        //    5    3     3     0
        //    6    4     3     0
        //    7    2     2     1
        //    8    3     2     1
        //    9    4     2     1
        //   10    2     3     1
        //   11    3     3     1
        //   12    4     3     1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 July 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of components in the vector.
        //
        //    Input, int X_MIN[N], X_MAX[N], the minimum and maximum
        //    values allowed in each component.
        //
        //    Input/output, int X[N].  On first call, with 
        //    MORE = FALSE, the input value of X is not important.  On subsequent calls,
        //    the input value of X should be the output value from the previous call.
        //    On output, with MORE = TRUE, the value of X will be the "next"
        //    vector in the reverse lexicographical list of vectors.  However, on 
        //    output with MORE = FALSE, the vector X is meaningless, because there 
        //    are no more vectors in the list.
        //
        //    Input/output, bool &MORE.  On input, if the user has set MORE
        //    FALSE, the user is requesting the initiation of a new sequence
        //    of values.  If MORE is TRUE, then the user is requesting "more"
        //    values in the current sequence.  On output, if MORE is TRUE,
        //    then another value was found and returned in X, but if MORE is
        //    FALSE, then there are no more values in the sequence, and X is
        //    NOT the next value.
        //
    {
        int i;

        switch (more)
        {
            case false:
            {
                for (i = 0; i < n; i++)
                {
                    x[i] = x_min[i];
                }

                more = true;
                break;
            }
            default:
            {
                i = 0;

                for (;;)
                {
                    if (x[i] < x_max[i])
                    {
                        x[i] += 1;
                        break;
                    }

                    x[i] = x_min[i];

                    i += 1;

                    if (n > i)
                    {
                        continue;
                    }

                    more = false;
                    break;
                }

                break;
            }
        }
    }
}