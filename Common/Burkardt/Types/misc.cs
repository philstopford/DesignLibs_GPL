using System;
using Burkardt.Uniform;
using entropyRNG;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void lvec_print(int n, bool[] a, string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LVEC_PRINT prints a logical vector.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 April 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of components of the vector.
            //
            //    Input, bool A[N], the vector to be printed.
            //
            //    Input, string TITLE, a title.
            //
        {
            int i;

            Console.WriteLine("");
            Console.WriteLine(title + "");
            Console.WriteLine("");
            for (i = 0; i < n; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(8)
                                       + "  " + a[i].ToString().PadLeft(1) + "");
            }
        }

        public static void lvec_next(int n, ref bool[] lvec)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LVEC_NEXT generates the next logical vector.
            //
            //  Discussion:
            //
            //    In the following discussion, we will let '0' stand for FALSE and
            //    '1' for TRUE.
            //
            //    The vectors have the order
            //
            //      (0,0,...,0),
            //      (0,0,...,1), 
            //      ...
            //      (1,1,...,1)
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
            //    0 0 0  =>  0 0 1
            //    0 0 1  =>  0 1 0
            //    0 1 0  =>  0 1 1
            //    0 1 1  =>  1 0 0
            //    1 0 0  =>  1 0 1
            //    1 0 1  =>  1 1 0
            //    1 1 0  =>  1 1 1
            //    1 1 1  =>  0 0 0
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 May 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the dimension of the vectors.
            //
            //    Input/output, bool LVEC[N], on output, the successor to the
            //    input vector.  
            //
        {
            int i;

            for (i = n - 1; 0 <= i; i--)
            {
                if (!lvec[i])
                {
                    lvec[i] = true;
                    return;
                }

                lvec[i] = false;
            }

            return;
        }

        public static int get_seed()
        {
            return RNG.nextint(1, Int32.MaxValue);
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


        public static int swapec(int i, ref int top, ref int btri, ref int bedg, int node_num,
            double[] node_xy, int triangle_num, ref int[] triangle_node,
        ref int[] triangle_neighbor, int[] stack )

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
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
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
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
        //
        //    Input, int TRIANGLE_NUM, the number of triangles.
        //
        //    Input/output, int TRIANGLE_NODE[3*TRIANGLE_NUM], the triangle incidence
        //    list.  May be updated on output because of swaps.
        //
        //    Input/output, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the triangle neighbor
        //    list; negative values are used for links of the counter-clockwise linked
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
            x = node_xy[2 * (i - 1) + 0];
            y = node_xy[2 * (i - 1) + 1];

            for (;;)
            {
                if (top <= 0)
                {
                    break;
                }

                t = stack[(top) - 1];
                top = top - 1;

                if (triangle_node[3 * (t - 1) + 0] == i)
                {
                    e = 2;
                    b = triangle_node[3 * (t - 1) + 2];
                }
                else if (triangle_node[3 * (t - 1) + 1] == i)
                {
                    e = 3;
                    b = triangle_node[3 * (t - 1) + 0];
                }
                else
                {
                    e = 1;
                    b = triangle_node[3 * (t - 1) + 1];
                }

                a = triangle_node[3 * (t - 1) + e - 1];
                u = triangle_neighbor[3 * (t - 1) + e - 1];

                if (triangle_neighbor[3 * (u - 1) + 0] == t)
                {
                    f = 1;
                    c = triangle_node[3 * (u - 1) + 2];
                }
                else if (triangle_neighbor[3 * (u - 1) + 1] == t)
                {
                    f = 2;
                    c = triangle_node[3 * (u - 1) + 0];
                }
                else
                {
                    f = 3;
                    c = triangle_node[3 * (u - 1) + 1];
                }

                swap = diaedg(x, y,
                    node_xy[2 * (a - 1) + 0], node_xy[2 * (a - 1) + 1],
                    node_xy[2 * (c - 1) + 0], node_xy[2 * (c - 1) + 1],
                    node_xy[2 * (b - 1) + 0], node_xy[2 * (b - 1) + 1]);

                if (swap == 1)
                {
                    em1 = i4_wrap(e - 1, 1, 3);
                    ep1 = i4_wrap(e + 1, 1, 3);
                    fm1 = i4_wrap(f - 1, 1, 3);
                    fp1 = i4_wrap(f + 1, 1, 3);

                    triangle_node[3 * (t - 1) + ep1 - 1] = c;
                    triangle_node[3 * (u - 1) + fp1 - 1] = i;
                    r = triangle_neighbor[3 * (t - 1) + ep1 - 1];
                    s = triangle_neighbor[3 * (u - 1) + fp1 - 1];
                    triangle_neighbor[3 * (t - 1) + ep1 - 1] = u;
                    triangle_neighbor[3 * (u - 1) + fp1 - 1] = t;
                    triangle_neighbor[3 * (t - 1) + e - 1] = s;
                    triangle_neighbor[3 * (u - 1) + f - 1] = r;

                    if (0 < triangle_neighbor[3 * (u - 1) + fm1 - 1])
                    {
                        top = top + 1;
                        stack[(top) - 1] = u;
                    }

                    if (0 < s)
                    {
                        if (triangle_neighbor[3 * (s - 1) + 0] == u)
                        {
                            triangle_neighbor[3 * (s - 1) + 0] = t;
                        }
                        else if (triangle_neighbor[3 * (s - 1) + 1] == u)
                        {
                            triangle_neighbor[3 * (s - 1) + 1] = t;
                        }
                        else
                        {
                            triangle_neighbor[3 * (s - 1) + 2] = t;
                        }

                        top = top + 1;

                        if (node_num < top)
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

                        while (0 < triangle_neighbor[3 * (tt - 1) + ee - 1])
                        {
                            tt = triangle_neighbor[3 * (tt - 1) + ee - 1];

                            if (triangle_node[3 * (tt - 1) + 0] == a)
                            {
                                ee = 3;
                            }
                            else if (triangle_node[3 * (tt - 1) + 1] == a)
                            {
                                ee = 1;
                            }
                            else
                            {
                                ee = 2;
                            }

                        }

                        triangle_neighbor[3 * (tt - 1) + ee - 1] = l;

                    }

                    if (0 < r)
                    {
                        if (triangle_neighbor[3 * (r - 1) + 0] == t)
                        {
                            triangle_neighbor[3 * (r - 1) + 0] = u;
                        }
                        else if (triangle_neighbor[3 * (r - 1) + 1] == t)
                        {
                            triangle_neighbor[3 * (r - 1) + 1] = u;
                        }
                        else
                        {
                            triangle_neighbor[3 * (r - 1) + 2] = u;
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

                        while (0 < triangle_neighbor[3 * (tt - 1) + ee - 1])
                        {
                            tt = triangle_neighbor[3 * (tt - 1) + ee - 1];

                            if (triangle_node[3 * (tt - 1) + 0] == b)
                            {
                                ee = 3;
                            }
                            else if (triangle_node[3 * (tt - 1) + 1] == b)
                            {
                                ee = 1;
                            }
                            else
                            {
                                ee = 2;
                            }

                        }

                        triangle_neighbor[3 * (tt - 1) + ee - 1] = l;

                    }

                }

            }

            return 0;
        }
    }
}