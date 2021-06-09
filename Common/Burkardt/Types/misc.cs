using System;
using Burkardt.Uniform;
using entropyRNG;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static int get_seed()
        {
            return RNG.nextint(1, Int32.MaxValue);
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


        public static double[] monomial_value(int m, int n, int[] e, double[] x)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MONOMIAL_VALUE evaluates a monomial.
            //
            //  Discussion:
            //
            //    This routine evaluates a monomial of the form
            //
            //      product ( 1 <= i <= m ) x(i)^e(i)
            //
            //    where the exponents are nonnegative integers.  Note that
            //    if the combination 0^0 is encountered, it should be treated
            //    as 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 May 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the spatial dimension.
            //
            //    Input, int N, the number of points at which the
            //    monomial is to be evaluated.
            //
            //    Input, int E[M], the exponents.
            //
            //    Input, double X[M*N], the point coordinates.
            //
            //    Output, double MONOMIAL_VALUE[N], the value of the monomial.
            //
        {
            double[] v = new double[n];

            for (int j = 0; j < n; j++)
            {
                v[j] = 1.0;
            }

            for (int i = 0; i < m; i++)
            {
                if (0 != e[i])
                {
                    for (int j = 0; j < n; j++)
                    {
                        v[j] = v[j] * Math.Pow(x[i + j * m], e[i]);
                    }
                }
            }

            return v;
        }

        public static int inits(double[] dos, int nos, double eta)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    INITS initializes a Chebyshev series.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    15 September 2011
            //
            //  Author:
            //
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Roger Broucke,
            //    Algorithm 446:
            //    Ten Subroutines for the Manipulation of Chebyshev Series,
            //    Communications of the ACM,
            //    Volume 16, Number 4, April 1973, pages 254-256.
            //
            //  Parameters:
            //
            //    Input, double DOS[NOS], the Chebyshev coefficients.
            //
            //    Input, int NOS, the number of coefficients.
            //
            //    Input, double ETA, the desired accuracy.
            //
            //    Output, int INITS, the number of terms of the series needed
            //    to ensure the requested accuracy.
            //
        {
            double err;
            int i;
            int value;

            if (nos < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("INITS - Fatal error!");
                Console.WriteLine("  Number of coefficients < 1.");
                return (1);
            }

            err = 0.0;

            for (i = nos - 1; 0 <= i; i--)
            {
                err = err + Math.Abs(dos[i]);
                if (eta < err)
                {
                    value = i + 1;
                    return value;
                }
            }

            value = i;
            Console.WriteLine("");
            Console.WriteLine("INITS - Warning!");
            Console.WriteLine("  ETA may be too small.");

            return value;
        }

        public static double csevl(double x, double[] a, int n)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CSEVL evaluates a Chebyshev series.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    15 September 2011
            //
            //  Author:
            //
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Roger Broucke,
            //    Algorithm 446:
            //    Ten Subroutines for the Manipulation of Chebyshev Series,
            //    Communications of the ACM,
            //    Volume 16, Number 4, April 1973, pages 254-256.
            //
            //  Parameters:
            //
            //    Input, double X, the evaluation point.
            //
            //    Input, double CS[N], the Chebyshev coefficients.
            //
            //    Input, int N, the number of Chebyshev coefficients.
            //
            //    Output, double CSEVL, the Chebyshev series evaluated at X.
            //
        {
            double b0;
            double b1;
            double b2 = 0;
            int i;
            double twox;
            double value;

            if (n < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("CSEVL - Fatal error!");
                Console.WriteLine("  Number of terms <= 0.");
                return (1);
            }

            if (1000 < n)
            {
                Console.WriteLine("");
                Console.WriteLine("CSEVL - Fatal error!");
                Console.WriteLine("  Number of terms greater than 1000.");
                return (1);
            }

            if (x < -1.1 || 1.1 < x)
            {
                Console.WriteLine("");
                Console.WriteLine("CSEVL - Fatal error!");
                Console.WriteLine("  X outside (-1,+1).");
                return (1);
            }

            twox = 2.0 * x;
            b1 = 0.0;
            b0 = 0.0;

            for (i = n - 1; 0 <= i; i--)
            {
                b2 = b1;
                b1 = b0;
                b0 = twox * b1 - b2 + a[i];
            }

            value = 0.5 * (b0 - b2);

            return value;
        }
        
        public static int diaedg(double x0, double y0, double x1, double y1, double x2, double y2,
                double x3, double y3)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DIAEDG chooses a diagonal edge.
            //
            //  Discussion:
            //
            //    The routine determines whether 0--2 or 1--3 is the diagonal edge
            //    that should be chosen, based on the circumcircle criterion, where
            //    (X0,Y0), (X1,Y1), (X2,Y2), (X3,Y3) are the vertices of a simple
            //    quadrilateral in counterclockwise order.
            //
            //  Modified:
            //
            //    28 August 2003
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Barry Joe.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Barry Joe,
            //    GEOMPACK - a software package for the generation of meshes
            //    using geometric algorithms,
            //    Advances in Engineering Software,
            //    Volume 13, pages 325-331, 1991.
            //
            //  Parameters:
            //
            //    Input, double X0, Y0, X1, Y1, X2, Y2, X3, Y3, the coordinates of the
            //    vertices of a quadrilateral, given in counter clockwise order.
            //
            //    Output, int DIAEDG, chooses a diagonal:
            //    +1, if diagonal edge 02 is chosen;
            //    -1, if diagonal edge 13 is chosen;
            //     0, if the four vertices are cocircular.
            //
        {
            double ca;
            double cb;
            double dx10;
            double dx12;
            double dx30;
            double dx32;
            double dy10;
            double dy12;
            double dy30;
            double dy32;
            double s;
            double tol;
            double tola;
            double tolb;
            int value;

            tol = 100.0 * double.Epsilon;

            dx10 = x1 - x0;
            dy10 = y1 - y0;
            dx12 = x1 - x2;
            dy12 = y1 - y2;
            dx30 = x3 - x0;
            dy30 = y3 - y0;
            dx32 = x3 - x2;
            dy32 = y3 - y2;

            tola = tol * Math.Max(Math.Abs(dx10),
                Math.Max(Math.Abs(dy10),
                    Math.Max(Math.Abs(dx30), Math.Abs(dy30))));

            tolb = tol * Math.Max(Math.Abs(dx12),
                Math.Max(Math.Abs(dy12),
                    Math.Max(Math.Abs(dx32), Math.Abs(dy32))));

            ca = dx10 * dx30 + dy10 * dy30;
            cb = dx12 * dx32 + dy12 * dy32;

            if (tola < ca && tolb < cb)
            {
                value = -1;
            }
            else if (ca < -tola && cb < -tolb)
            {
                value = 1;
            }
            else
            {
                tola = Math.Max(tola, tolb);
                s = (dx10 * dy30 - dx30 * dy10) * cb
                    + (dx32 * dy12 - dx12 * dy32) * ca;

                if (tola < s)
                {
                    value = -1;
                }
                else if (s < -tola)
                {
                    value = 1;
                }
                else
                {
                    value = 0;
                }

            }

            return value;
        }


        public static int swapec(int i, ref int top, ref int btri, ref int bedg, int point_num,
                double[] point_xy, int tri_num, ref int[] tri_vert, ref int[] tri_nabe,
                int[] stack)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SWAPEC swaps diagonal edges until all triangles are Delaunay.
            //
            //  Discussion:
            //
            //    The routine swaps diagonal edges in a 2D triangulation, based on
            //    the empty circumcircle criterion, until all triangles are Delaunay,
            //    given that I is the index of the new vertex added to the triangulation.
            //
            //  Modified:
            //
            //    03 September 2003
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Barry Joe.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Barry Joe,
            //    GEOMPACK - a software package for the generation of meshes
            //    using geometric algorithms,
            //    Advances in Engineering Software,
            //    Volume 13, pages 325-331, 1991.
            //
            //  Parameters:
            //
            //    Input, int I, the index of the new vertex.
            //
            //    Input/output, int *TOP, the index of the top of the stack.
            //    On output, TOP is zero.
            //
            //    Input/output, int *BTRI, *BEDG; on input, if positive, are the
            //    triangle and edge indices of a boundary edge whose updated indices
            //    must be recorded.  On output, these may be updated because of swaps.
            //
            //    Input, int POINT_NUM, the number of points.
            //
            //    Input, double POINT_XY[POINT_NUM*2], the coordinates of the points.
            //
            //    Input, int TRI_NUM, the number of triangles.
            //
            //    Input/output, int TRI_VERT[TRI_NUM*3], the triangle incidence list.
            //    May be updated on output because of swaps.
            //
            //    Input/output, int TRI_NABE[TRI_NUM*3], the triangle neighbor list;
            //    negative values are used for links of the counter-clockwise linked
            //    list of boundary edges;  May be updated on output because of swaps.
            //
            //      LINK = -(3*I + J-1) where I, J = triangle, edge index.
            //
            //    Workspace, int STACK[MAXST]; on input, entries 1 through TOP
            //    contain the indices of initial triangles (involving vertex I)
            //    put in stack; the edges opposite I should be in interior;  entries
            //    TOP+1 through MAXST are used as a stack.
            //
            //    Output, int SWAPEC, is set to 8 for abnormal return.
            //
        {
            int a;
            int b;
            int c;
            int e;
            int ee;
            int em1;
            int ep1;
            int f;
            int fm1;
            int fp1;
            int l;
            int r;
            int s;
            int swap;
            int t;
            int tt;
            int u;
            double x;
            double y;
            //
            //  Determine whether triangles in stack are Delaunay, and swap
            //  diagonal edge of convex quadrilateral if not.
            //
            x = point_xy[2 * (i - 1) + 0];
            y = point_xy[2 * (i - 1) + 1];

            for (;;)
            {
                if (top <= 0)
                {
                    break;
                }

                t = stack[(top) - 1];
                top = top - 1;

                if (tri_vert[3 * (t - 1) + 0] == i)
                {
                    e = 2;
                    b = tri_vert[3 * (t - 1) + 2];
                }
                else if (tri_vert[3 * (t - 1) + 1] == i)
                {
                    e = 3;
                    b = tri_vert[3 * (t - 1) + 0];
                }
                else
                {
                    e = 1;
                    b = tri_vert[3 * (t - 1) + 1];
                }

                a = tri_vert[3 * (t - 1) + e - 1];
                u = tri_nabe[3 * (t - 1) + e - 1];

                if (tri_nabe[3 * (u - 1) + 0] == t)
                {
                    f = 1;
                    c = tri_vert[3 * (u - 1) + 2];
                }
                else if (tri_nabe[3 * (u - 1) + 1] == t)
                {
                    f = 2;
                    c = tri_vert[3 * (u - 1) + 0];
                }
                else
                {
                    f = 3;
                    c = tri_vert[3 * (u - 1) + 1];
                }

                swap = diaedg(x, y,
                    point_xy[2 * (a - 1) + 0], point_xy[2 * (a - 1) + 1],
                    point_xy[2 * (c - 1) + 0], point_xy[2 * (c - 1) + 1],
                    point_xy[2 * (b - 1) + 0], point_xy[2 * (b - 1) + 1]);

                if (swap == 1)
                {
                    em1 = i4_wrap(e - 1, 1, 3);
                    ep1 = i4_wrap(e + 1, 1, 3);
                    fm1 = i4_wrap(f - 1, 1, 3);
                    fp1 = i4_wrap(f + 1, 1, 3);

                    tri_vert[3 * (t - 1) + ep1 - 1] = c;
                    tri_vert[3 * (u - 1) + fp1 - 1] = i;
                    r = tri_nabe[3 * (t - 1) + ep1 - 1];
                    s = tri_nabe[3 * (u - 1) + fp1 - 1];
                    tri_nabe[3 * (t - 1) + ep1 - 1] = u;
                    tri_nabe[3 * (u - 1) + fp1 - 1] = t;
                    tri_nabe[3 * (t - 1) + e - 1] = s;
                    tri_nabe[3 * (u - 1) + f - 1] = r;

                    if (0 < tri_nabe[3 * (u - 1) + fm1 - 1])
                    {
                        top = top + 1;
                        stack[(top) - 1] = u;
                    }

                    if (0 < s)
                    {
                        if (tri_nabe[3 * (s - 1) + 0] == u)
                        {
                            tri_nabe[3 * (s - 1) + 0] = t;
                        }
                        else if (tri_nabe[3 * (s - 1) + 1] == u)
                        {
                            tri_nabe[3 * (s - 1) + 1] = t;
                        }
                        else
                        {
                            tri_nabe[3 * (s - 1) + 2] = t;
                        }

                        top = top + 1;

                        if (point_num < top)
                        {
                            return 8;
                        }

                        stack[(top) - 1] = t;
                    }
                    else
                    {
                        if (u == btri && fp1 == bedg)
                        {
                            btri = t;
                            bedg = e;
                        }

                        l = -(3 * t + e - 1);
                        tt = t;
                        ee = em1;

                        while (0 < tri_nabe[3 * (tt - 1) + ee - 1])
                        {
                            tt = tri_nabe[3 * (tt - 1) + ee - 1];

                            if (tri_vert[3 * (tt - 1) + 0] == a)
                            {
                                ee = 3;
                            }
                            else if (tri_vert[3 * (tt - 1) + 1] == a)
                            {
                                ee = 1;
                            }
                            else
                            {
                                ee = 2;
                            }

                        }

                        tri_nabe[3 * (tt - 1) + ee - 1] = l;

                    }

                    if (0 < r)
                    {
                        if (tri_nabe[3 * (r - 1) + 0] == t)
                        {
                            tri_nabe[3 * (r - 1) + 0] = u;
                        }
                        else if (tri_nabe[3 * (r - 1) + 1] == t)
                        {
                            tri_nabe[3 * (r - 1) + 1] = u;
                        }
                        else
                        {
                            tri_nabe[3 * (r - 1) + 2] = u;
                        }
                    }
                    else
                    {
                        if (t == btri && ep1 == bedg)
                        {
                            btri = u;
                            bedg = f;
                        }

                        l = -(3 * u + f - 1);
                        tt = u;
                        ee = fm1;

                        while (0 < tri_nabe[3 * (tt - 1) + ee - 1])
                        {
                            tt = tri_nabe[3 * (tt - 1) + ee - 1];

                            if (tri_vert[3 * (tt - 1) + 0] == b)
                            {
                                ee = 3;
                            }
                            else if (tri_vert[3 * (tt - 1) + 1] == b)
                            {
                                ee = 1;
                            }
                            else
                            {
                                ee = 2;
                            }
                        }

                        tri_nabe[3 * (tt - 1) + ee - 1] = l;
                    }
                }
            }

            return 0;
        }
    }
}