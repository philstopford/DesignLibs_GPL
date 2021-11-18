using Burkardt.Types;

namespace Burkardt.Square;

public static partial class MinimalRule
{
    public static double[] smr08()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SMR08 returns the SMR rule of degree 8.
        //
        //  Discussion:
        //
        //    R.COOLS SAYS THERE IS A RULE WITH 15 PTS.
        //
        //    DEGREE: 8
        //    POINTS CARDINALITY: 16
        //    NORM INF MOMS. RESIDUAL: 5.27356e-16
        //    SUM NEGATIVE WEIGHTS: 0.00000e+00,
        // 
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 February 2018
        //
        //  Author:
        //
        //    Original MATLAB version by Mattia Festa, Alvise Sommariva,
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Mattia Festa, Alvise Sommariva,
        //    Computing almost minimal formulas on the square,
        //    Journal of Computational and Applied Mathematics,
        //    Volume 17, Number 236, November 2012, pages 4296-4302.
        //
        //  Parameters:
        //
        //    Output, double *SMR08[3*16], the requested rule.
        //
    {
        const int degree = 8;
        double[] xw =
        {
            9.200233979246099e-01, -1.061811254185043e-01, 2.066851970823297e-01,
            8.930854460206108e-01, 7.538456868810103e-01, 1.637084884729354e-01,
            6.675369061590991e-01, -6.602769863240182e-01, 3.448725288612837e-01,
            9.736841855611369e-01, -9.418906397956112e-01, 4.154990098574050e-02,
            4.634059670545806e-01, 9.713146222867728e-01, 1.207701070516287e-01,
            5.257251055722159e-01, 3.284384783969826e-01, 4.761612149298055e-01,
            -3.171906088585705e-01, 8.431651719683217e-01, 8.714564872094989e-02,
            -6.189113622341841e-02, 6.853227004655966e-01, 3.755078795520529e-01,
            6.720825424065331e-02, -2.654137657769545e-01, 5.939493422195995e-01,
            1.957478656950412e-01, -9.531380985058940e-01, 1.576823081145121e-01,
            -7.280591295595631e-01, 9.245623167165428e-01, 1.261748042046440e-01,
            -5.538775451118830e-01, 1.849451143123379e-01, 4.968422862065046e-01,
            -4.448446986979669e-01, -7.417922043114903e-01, 3.611024723493892e-01,
            -9.421081322947572e-01, 5.652760298758271e-01, 1.413572491732567e-01,
            -8.900106969566213e-01, -3.626326971963237e-01, 2.363188753269229e-01,
            -8.841558907294489e-01, -9.418511077281306e-01, 7.017169674844463e-02
        };

        int order = square_minimal_rule_order(degree);
        double[] xw_copy = typeMethods.r8mat_copy_new(3, order, xw);

        return xw_copy;
    }
}