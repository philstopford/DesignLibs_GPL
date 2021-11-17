namespace Burkardt.Grid;

public static class Transfer
{
    public static void ctof(int nc, double[] uc, int nf, ref double[] uf, int ucIndex = 0, int ufIndex = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CTOF transfers data from a coarse to a finer grid.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 December 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    William Hager,
        //    Applied Numerical Linear Algebra,
        //    Prentice-Hall, 1988,
        //    ISBN13: 978-0130412942,
        //    LC: QA184.H33.
        //
        //  Parameters:
        //
        //    Input, int NC, the number of coarse nodes.
        //
        //    Input, double UC[NC], the coarse correction data.
        //
        //    Input, int NF, the number of fine nodes.
        //
        //    Input/output, double UF[NF], on input, the fine grid data.
        //    On output, the data has been updated with prolonged coarse 
        //    correction data.
        //
    {
        int ic;
        int iff;

        for (ic = 0; ic < nc; ic++)
        {
            iff = 2 * ic;
            uf[ufIndex + iff] += uc[ucIndex + ic];
        }

        for (ic = 0; ic < nc - 1; ic++)
        {
            iff = 2 * ic + 1;
            uf[ufIndex + iff] += 0.5 * (uc[ucIndex + ic] + uc[ucIndex + ic + 1]);
        }

    }


    public static void ftoc(int nf, double[] uf, double[] rf, int nc, ref double[] uc,
            ref double[] rc, int ufIndex = 0, int rfIndex = 0, int ucIndex = 0, int rcIndex = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FTOC transfers data from a fine grid to a coarser grid.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 December 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    William Hager,
        //    Applied Numerical Linear Algebra,
        //    Prentice-Hall, 1988,
        //    ISBN13: 978-0130412942,
        //    LC: QA184.H33.
        //
        //  Parameters:
        //
        //    Input, int NF, the number of fine nodes.
        //
        //    Input, double UF[NF], the fine data.
        //
        //    Input, double RF[NF], the right hand side for the fine grid.
        //
        //    Input, int NC, the number of coarse nodes.
        //
        //    Output, double UC[NC], the coarse grid data, set to zero.
        //
        //    Output, double RC[NC], the right hand side for the coarse grid.
        //
    {
        int ic;
        int iff;

        for (ic = 0; ic < nc; ic++)
        {
            uc[ucIndex + ic] = 0.0;
        }

        rc[0] = 0.0;
        for (ic = 1; ic < nc - 1; ic++)
        {
            iff = 2 * ic;
            rc[rcIndex + ic] = 4.0 * (rf[rfIndex + iff] + uf[ufIndex + iff - 1] - 2.0 * uf[ufIndex + iff] + uf[ufIndex + iff + 1]);
        }

        rc[rcIndex + nc - 1] = 0.0;

    }
}