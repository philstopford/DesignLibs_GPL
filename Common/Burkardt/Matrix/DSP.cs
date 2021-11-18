using System;
using System.Globalization;
using Burkardt.FEM;
using Burkardt.Quadrature;
using Burkardt.Types;

namespace Burkardt.MatrixNS;

public static class DSP
{
    public static void dirichlet_apply_dsp ( int node_num, double[] node_xy, int[] node_condition,
            int nz_num, int[] ia, int[] ja, ref double[] a, ref double[] f, Func<int, double[], double[]> dirichlet_condition )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIRICHLET_APPLY_DSP accounts for Dirichlet boundary conditions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the coordinates of nodes.
        //
        //    Input, int NODE_CONDITION[NODE_NUM], reports the condition
        //    used to set the unknown associated with the node.
        //    0, unknown.
        //    1, finite element equation.
        //    2, Dirichlet condition;
        //    3, Neumann condition.
        //
        //    Input, int NZ_NUM, the number of nonzero entries.
        //
        //    Input, int IA[NZ_NUM], JA[NZ_NUM], the row and column
        //    indices of the nonzero entries.
        //
        //    Input/output, double A[NZ_NUM], the nonzero entries of the matrix.
        //    On output, adjusted to account for Dirichlet boundary conditions.
        //
        //    Input/output, double F[NODE_NUM], the right hand side.
        //    On output, adjusted to account for Dirichlet boundary conditions.
        //
    {
        int column;
        int DIRICHLET = 2;
        int node;
        double[] node_bc;
        int nz;
            
        node_bc = dirichlet_condition(node_num, node_xy);
        //
        //  Consider every matrix entry, NZ.
        //
        //  If the row I corresponds to a boundary node, then
        //  zero out all off diagonal matrix entries, set the diagonal to 1,
        //  and the right hand side to the Dirichlet boundary condition value.
        //
        for (nz = 0; nz < nz_num; nz++)
        {
            node = ia[nz];

            if (node_condition[node - 1] == DIRICHLET)
            {
                column = ja[nz];

                if (column == node)
                {
                    a[nz] = 1.0;
                    f[node - 1] = node_bc[node - 1];
                }
                else
                {
                    a[nz] = 0.0;
                }
            }
        }
    }

    public static void assemble_poisson_dsp ( int node_num, double[] node_xy,
            int element_num, int[] element_node, int quad_num, int nz_num, int[] ia,
            int[] ja, ref double[] a, ref double[] f, Func<int, double[], double[]> rhs, Func<int, double[], double[]> h_coef, Func<int, double[], double[]> k_coef )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ASSEMBLE_POISSON_DSP assembles the system for the Poisson equation.
        //
        //  Discussion:
        //
        //    The matrix is sparse, and stored in the DSP or "sparse triple" format.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the
        //    coordinates of nodes.
        //
        //    Input, int ELEMENT_NUM, the number of triangles.
        //
        //    Input, int ELEMENT_NODE[3*ELEMENT_NUM];
        //    ELEMENT_NODE(I,J) is the global index of local node I in triangle J.
        //
        //    Input, int QUAD_NUM, the number of quadrature points used in assembly.
        //
        //    Input, int NZ_NUM, the number of nonzero entries.
        //
        //    Input, int IA[NZ_NUM], JA[NZ_NUM], the row and column
        //    indices of the nonzero entries.
        //
        //    Output, double A[NZ_NUM], the nonzero entries of the matrix.
        //
        //    Output, double F(NODE_NUM), the right hand side.
        //
        //  Local parameters:
        //
        //    Local, double BI, DBIDX, DBIDY, the value of some basis function
        //    and its first derivatives at a quadrature point.
        //
        //    Local, double BJ, DBJDX, DBJDY, the value of another basis
        //    function and its first derivatives at a quadrature point.
        //
    {
        double area;
        int basis;
        double bi = 0;
        double bj = 0;
        double dbidx = 0;
        double dbidy = 0;
        double dbjdx = 0;
        double dbjdy = 0;
        int element;
        int i;
        int j;
        int k;
        int node;
        int nz;
        double[] p = new double[2];
        double[] phys_h;
        double[] phys_k;
        double[] phys_rhs;
        double[] phys_xy;
        int quad;
        double[] quad_w;
        double[] quad_xy;
        double[] t3 = new double[2 * 3];
        int test;
        double[] w;

        phys_h = new double[quad_num];
        phys_k = new double[quad_num];
        phys_rhs = new double[quad_num];
        phys_xy = new double[2 * quad_num];

        quad_w = new double[quad_num];
        quad_xy = new double[2 * quad_num];

        w = new double[quad_num];
        //
        //  Initialize the arrays to zero.
        //
        for (node = 0; node < node_num; node++)
        {
            f[node] = 0.0;
        }

        for (nz = 0; nz < nz_num; nz++)
        {
            a[nz] = 0.0;
        }

        //
        //  Get the quadrature weights and nodes.
        //
        QuadratureRule.quad_rule(quad_num, ref quad_w, ref quad_xy);
        //
        //  Add up all quantities associated with the ELEMENT-th element.
        //
        for (element = 0; element < element_num; element++)
        {
            //
            //  Make a copy of the element.
            //
            for (j = 0; j < 3; j++)
            {
                for (i = 0; i < 2; i++)
                {
                    t3[i + j * 2] = node_xy[i + (element_node[j + element * 3] - 1) * 2];
                }
            }

            //
            //  Map the quadrature points QUAD_XY to points XY in the physical element.
            //
            Reference.reference_to_physical_t3(t3, quad_num, quad_xy, ref phys_xy);

            area = Math.Abs(typeMethods.triangle_area_2d(t3));

            for (quad = 0; quad < quad_num; quad++)
            {
                w[quad] = quad_w[quad] * area;
            }

            phys_rhs = rhs(quad_num, phys_xy);
            phys_h = h_coef(quad_num, phys_xy);
            phys_k = k_coef(quad_num, phys_xy);
            //
            //  Consider the QUAD-th quadrature point.
            //
            for (quad = 0; quad < quad_num; quad++)
            {
                p[0] = phys_xy[0 + quad * 2];
                p[1] = phys_xy[1 + quad * 2];
                //
                //  Consider the TEST-th test function.
                //
                //  We generate an integral for every node associated with an unknown.
                //  But if a node is associated with a boundary condition, we do nothing.
                //
                for (test = 1; test <= 3; test++)
                {
                    i = element_node[test - 1 + element * 3];

                    Basis11.basis_one_t3(t3, test, p, ref bi, ref dbidx, ref dbidy);

                    f[i - 1] += w[quad] * phys_rhs[quad] * bi;
                    //
                    //  Consider the BASIS-th basis function, which is used to form the
                    //  value of the solution function.
                    //
                    for (basis = 1; basis <= 3; basis++)
                    {
                        j = element_node[basis - 1 + element * 3];

                        Basis11.basis_one_t3(t3, basis, p, ref bj, ref dbjdx, ref dbjdy);

                        k = dsp_ij_to_k(nz_num, ia, ja, i, j);

                        a[k - 1] += w[quad] * (
                            phys_h[quad] * (dbidx * dbjdx + dbidy * dbjdy)
                            + phys_k[quad] * bj * bi);
                    }
                }
            }
        }
    }

    public static int dsp_ij_to_k(int nz_num, int[] row, int[] col, int i, int j)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DSP_IJ_TO_K seeks the compressed index of the (I,J) entry of A.
        //
        //  Discussion:
        //
        //    If A(I,J) is nonzero, then its value is stored in location K.
        //
        //    This routine searches the DSP storage structure for the index K
        //    corresponding to (I,J), returning -1 if no such entry was found.
        //
        //    This routine assumes that the data structure has been sorted,
        //    so that the entries of ROW are ascending sorted, and that the
        //    entries of COL are ascending sorted, within the group of entries
        //    that have a common value of ROW.
        //
        //    The DSP storage format stores the row, column and value of each nonzero
        //    entry of a sparse matrix.
        //
        //    The DSP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP
        //    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NZ_NUM, the number of nonzero elements in
        //    the matrix.
        //
        //    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and
        //    column indices of the nonzero elements.
        //
        //    Input, int I, J, the row and column indices of the
        //    matrix entry.
        //
        //    Output, int DSP_IJ_TO_K, the DSP index of the (I,J) entry.
        //
    {
        int hi;
        int k;
        int lo;
        int md;

        lo = 1;
        hi = nz_num;

        for (;;)
        {
            if (hi < lo)
            {
                k = -1;
                break;
            }

            md = (lo + hi) / 2;

            if (row[md - 1] < i || row[md - 1] == i && col[md - 1] < j)
            {
                lo = md + 1;
            }
            else if (i < row[md - 1] || row[md - 1] == i && j < col[md - 1])
            {
                hi = md - 1;
            }
            else
            {
                k = md;
                break;
            }
        }

        return k;
    }

    public static void dsp_print_some(int m, int n, int nz_num, int[] row, int[] col,
            double[] a, int ilo, int jlo, int ihi, int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DSP_PRINT_SOME prints some of a DSP matrix.
        //
        //  Discussion:
        //
        //    This version of DSP_PRINT_SOME has been specifically modified to allow,
        //    and correctly handle, the case in which a single matrix location
        //    A(I,J) is referenced more than once by the sparse matrix structure.
        //    In such cases, the routine prints out the sum of all the values.
        //
        //    The DSP storage format stores the row, column and value of each nonzero
        //    entry of a sparse matrix.
        //
        //    It is possible that a pair of indices (I,J) may occur more than
        //    once.  Presumably, in this case, the intent is that the actual value
        //    of A(I,J) is the sum of all such entries.  This is not a good thing
        //    to do, but I seem to have come across this in MATLAB.
        //
        //    The DSP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP
        //    (nonsymmetric case), by MATLAB, and by SPARSEKIT ("COO" format).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero elements in the matrix.
        //
        //    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and column indices
        //    of the nonzero elements.
        //
        //    Input, double A[NZ_NUM], the nonzero elements of the matrix.
        //
        //    Input, int ILO, JLO, IHI, JHI, the first row and
        //    column, and the last row and column to be printed.
        //
        //    Input, string TITLE, a title to print.
        //
    {
        int INCX = 5;

        double[] aij = new double[INCX];
        int i;
        int i2hi;
        int i2lo;
        int inc;
        int j;
        int j2;
        int j2hi;
        int j2lo;
        int k;
        bool nonzero;

        Console.WriteLine("");
        Console.WriteLine(title + "");
        //
        //  Print the columns of the matrix, in strips of 5.
        //
        for (j2lo = jlo; j2lo <= jhi; j2lo += INCX)
        {
            j2hi = j2lo + INCX - 1;
            j2hi = Math.Min(j2hi, n);
            j2hi = Math.Min(j2hi, jhi);

            inc = j2hi + 1 - j2lo;

            Console.WriteLine("");

            string cout = "  Col:  ";
            for (j = j2lo; j <= j2hi; j++)
            {
                cout += j.ToString(CultureInfo.InvariantCulture).PadLeft(7) + "       ";
            }

            Console.WriteLine(cout);
            Console.WriteLine("  Row");
            Console.WriteLine("  ---");
            //
            //  Determine the range of the rows in this strip.
            //
            i2lo = Math.Max(ilo, 1);
            i2hi = Math.Min(ihi, m);

            for (i = i2lo; i <= i2hi; i++)
            {
                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                nonzero = false;
                for (j2 = 0; j2 < INCX; j2++)
                {
                    aij[j2] = 0.0;
                }

                for (k = 1; k <= nz_num; k++)
                {
                    if (i == row[k - 1] && j2lo <= col[k - 1] && col[k - 1] <= j2hi)
                    {
                        j2 = col[k - 1] - j2lo;

                        switch (a[k - 1])
                        {
                            case 0.0:
                                continue;
                            default:
                                nonzero = true;
                                aij[j2] += a[k - 1];
                                break;
                        }
                    }
                }

                switch (nonzero)
                {
                    case true:
                    {
                        cout = i.ToString(CultureInfo.InvariantCulture).PadLeft(6);
                        for (j2 = 0; j2 < inc; j2++)
                        {
                            cout += aij[j2].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  ";
                        }

                        Console.WriteLine(cout);
                        break;
                    }
                }
            }
        }
    }
}