using System;

namespace StroudTest
{
    public static class functions
    {
        public static double function_1d(int function_1d_index, double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FUNCTION_1D evaluates the current 1D function.
            //
            //  Discussion:
            //
            //    This routine assumes that the global variable FUNCTION_1D_INDEX has been
            //    set, and is accessible.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 April 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, the value of the variable.
            //
            //    Output, double FUNCTION_1D, the value of the function.
            //
        {
            double value;

            if (function_1d_index == 1)
            {
                value = 1.0;
            }
            else if (function_1d_index == 2)
            {
                value = x;
            }
            else if (function_1d_index == 3)
            {
                value = x * x;
            }
            else if (function_1d_index == 4)
            {
                value = x * x * x * x;
            }
            else if (function_1d_index == 5)
            {
                value = x * x * x * x;
            }
            else if (function_1d_index == 6)
            {
                value = x * x * x * x * x;
            }
            else if (function_1d_index == 7)
            {
                value = x * x * x * x * x * x;
            }
            else if (function_1d_index == 8)
            {
                value = Math.Abs(x);
            }
            else if (function_1d_index == 9)
            {
                value = Math.Sin(x);
            }
            else if (function_1d_index == 10)
            {
                value = Math.Exp(x);
            }
            else if (function_1d_index == 11)
            {
                value = 1.0 / (1.0 + Math.Abs(x));
            }
            else if (function_1d_index == 12)
            {
                value = Math.Sqrt(Math.Abs(x));
            }
            else
            {
                value = 0.0;
            }

            return value;
        }

        public static void function_1d_name(int function_1d_index, ref string name)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FUNCTION_1D_NAME returns the name of the current 1D function.
            //
            //  Discussion:
            //
            //    This routine assumes that the global variable FUNCTION_1D_INDEX has been
            //    set, and is accessible.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 April 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, char NAME[8], the name of the current 1D function.
            //
        {
            if (function_1d_index == 1)
            {
                name = "      1";
            }
            else if (function_1d_index == 2)
            {
                name = "      X";
            }
            else if (function_1d_index == 3)
            {
                name = "    X^2";
            }
            else if (function_1d_index == 4)
            {
                name = "    X^3";
            }
            else if (function_1d_index == 5)
            {
                name = "    X^4";
            }
            else if (function_1d_index == 6)
            {
                name = "    X^5";
            }
            else if (function_1d_index == 7)
            {
                name = "    X^6";
            }
            else if (function_1d_index == 8)
            {
                name = "      R";
            }
            else if (function_1d_index == 9)
            {
                name = " SIN(X)";
            }
            else if (function_1d_index == 10)
            {
                name = " EXP(X)";
            }
            else if (function_1d_index == 11)
            {
                name = "1/(1+R)";
            }
            else if (function_1d_index == 12)
            {
                name = "SQRT(R)";
            }
            else
            {
                name = "???????";
            }

        }

        public static int function_1d_num()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FUNCTION_1D_NUM returns the number of 1D functions.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 April 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, int FUNCTION_1D_NUM, the number of 1D functions.
            //
        {
            int value = 12;

            return value;
        }

        public static double function_2d(int function_2d_index, double x, double y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FUNCTION_2D evaluates the current 2D function.
            //
            //  Discussion:
            //
            //    This routine assumes that the global variable FUNCTION_2D_INDEX has been
            //    set, and is accessible.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 April 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, double Y, the value of the variables.
            //
            //    Output, double FUNCTION_2D, the value of the function.
            //
        {
            double value;

            if (function_2d_index == 1)
            {
                value = 1.0;
            }
            else if (function_2d_index == 2)
            {
                value = x;
            }
            else if (function_2d_index == 3)
            {
                value = x * x;
            }
            else if (function_2d_index == 4)
            {
                value = x * x * x;
            }
            else if (function_2d_index == 5)
            {
                value = x * x * x * x;
            }
            else if (function_2d_index == 6)
            {
                value = x * x * x * x * x;
            }
            else if (function_2d_index == 7)
            {
                value = x * x * x * x * x * x;
            }
            else if (function_2d_index == 8)
            {
                value = Math.Sqrt(x * x + y * y);
            }
            else if (function_2d_index == 9)
            {
                value = Math.Sin(x);
            }
            else if (function_2d_index == 10)
            {
                value = Math.Exp(x);
            }
            else if (function_2d_index == 11)
            {
                value = 1.0 / (1.0 + Math.Sqrt(x * x + y * y));
            }
            else if (function_2d_index == 12)
            {
                value = Math.Sqrt(Math.Sqrt(x * x + y * y));
            }
            else
            {
                value = 0.0;
            }

            return value;
        }

        public static void function_2d_name(int function_2d_index, ref string name)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FUNCTION_2D_NAME returns the name of the current 2D function.
            //
            //  Discussion:
            //
            //    This routine assumes that the global variable FUNCTION_2D_INDEX has been
            //    set, and is accessible.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 April 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, char NAME[8], the name of the current 2D function.
            //
        {
            if (function_2d_index == 1)
            {
                name = "      1";
            }
            else if (function_2d_index == 2)
            {
                name = "      X";
            }
            else if (function_2d_index == 3)
            {
                name = "    X^2";
            }
            else if (function_2d_index == 4)
            {
                name = "    X^3";
            }
            else if (function_2d_index == 5)
            {
                name = "    X^4";
            }
            else if (function_2d_index == 6)
            {
                name = "    X^5";
            }
            else if (function_2d_index == 7)
            {
                name = "    X^6";
            }
            else if (function_2d_index == 8)
            {
                name = "      R";
            }
            else if (function_2d_index == 9)
            {
                name = " SIN(X)";
            }
            else if (function_2d_index == 10)
            {
                name = " EXP(X)";
            }
            else if (function_2d_index == 11)
            {
                name = "1/(1+R)";
            }
            else if (function_2d_index == 12)
            {
                name = "SQRT(R)";
            }
            else
            {
                name = "???????";
            }

        }

        public static int function_2d_num()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FUNCTION_2D_NUM returns the number of 2D functions.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 April 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, int FUNCTION_2D_NUM, the number of 2D functions.
            //
        {
            int value = 12;

            return value;
        }

        public static double function_3d(int function_3d_index, double x, double y, double z)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FUNCTION_3D evaluates a function F(X,Y,Z) of 3 variables.
            //
            ///  Discussion:
            //
            //    This routine assumes that the global variable FUNCTION_3D_INDEX has been
            //    set, and is accessible.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    31 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, Y, Z, the value of the variables.
            //
            //    Output, double FUNC_3D, the value of the function.
            //
        {
            double value;

            if (function_3d_index == 1)
            {
                value = 1.0;
            }
            else if (function_3d_index == 2)
            {
                value = x;
            }
            else if (function_3d_index == 3)
            {
                value = y;
            }
            else if (function_3d_index == 4)
            {
                value = z;
            }
            else if (function_3d_index == 5)
            {
                value = x * x;
            }
            else if (function_3d_index == 6)
            {
                value = x * y;
            }
            else if (function_3d_index == 7)
            {
                value = x * z;
            }
            else if (function_3d_index == 8)
            {
                value = y * y;
            }
            else if (function_3d_index == 9)
            {
                value = y * z;
            }
            else if (function_3d_index == 10)
            {
                value = z * z;
            }
            else if (function_3d_index == 11)
            {
                value = x * x * x;
            }
            else if (function_3d_index == 12)
            {
                value = x * y * z;
            }
            else if (function_3d_index == 13)
            {
                value = z * z * z;
            }
            else if (function_3d_index == 14)
            {
                value = x * x * x * x;
            }
            else if (function_3d_index == 15)
            {
                value = x * x * z * z;
            }
            else if (function_3d_index == 16)
            {
                value = z * z * z * z;
            }
            else if (function_3d_index == 17)
            {
                value = x * x * x * x * x;
            }
            else if (function_3d_index == 18)
            {
                value = Math.Pow(x, 6);
            }
            else if (function_3d_index == 19)
            {
                value = Math.Sqrt(x * x + y * y + z * z);
            }
            else if (function_3d_index == 20)
            {
                value = Math.Sin(x);
            }
            else if (function_3d_index == 21)
            {
                value = Math.Exp(x);
            }
            else if (function_3d_index == 22)
            {
                value = 1.0 / Math.Sqrt(1.0 + x * x + y * y + z * z);
            }
            else if (function_3d_index == 23)
            {
                value = Math.Sqrt(Math.Sqrt(x * x + y * y + z * z));
            }
            else
            {
                value = 0.0;
            }

            return value;
        }

        public static void function_3d_name(int function_3d_index, ref string name)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FUNCTION_3D_NAME returns the name of the current 3D function.
            //
            //  Discussion:
            //
            //    This routine assumes that the global variable FUNCTION_3D_INDEX has been
            //    set, and is accessible.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    31 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, char NAME[8], the name of the current 3D function.
            //
        {
            if (function_3d_index == 1)
            {
                name = "      1";
            }
            else if (function_3d_index == 2)
            {
                name = "      X";
            }
            else if (function_3d_index == 3)
            {
                name = "      Y";
            }
            else if (function_3d_index == 4)
            {
                name = "      Z";
            }
            else if (function_3d_index == 5)
            {
                name = "    X*X";
            }
            else if (function_3d_index == 6)
            {
                name = "    X*Y";
            }
            else if (function_3d_index == 7)
            {
                name = "    X*Z";
            }
            else if (function_3d_index == 8)
            {
                name = "    Y*Y";
            }
            else if (function_3d_index == 9)
            {
                name = "    Y*Z";
            }
            else if (function_3d_index == 10)
            {
                name = "    Z*Z";
            }
            else if (function_3d_index == 11)
            {
                name = "    X^3";
            }
            else if (function_3d_index == 12)
            {
                name = "  X*Y*Z";
            }
            else if (function_3d_index == 13)
            {
                name = "  Z*Z*Z";
            }
            else if (function_3d_index == 14)
            {
                name = "    X^4";
            }
            else if (function_3d_index == 15)
            {
                name = "X^2 Z^2";
            }
            else if (function_3d_index == 16)
            {
                name = "    Z^4";
            }
            else if (function_3d_index == 17)
            {
                name = "    X^5";
            }
            else if (function_3d_index == 18)
            {
                name = "    X^6";
            }
            else if (function_3d_index == 19)
            {
                name = "      R";
            }
            else if (function_3d_index == 20)
            {
                name = " SIN(X)";
            }
            else if (function_3d_index == 21)
            {
                name = " EXP(X)";
            }
            else if (function_3d_index == 22)
            {
                name = "1/(1+R)";
            }
            else if (function_3d_index == 23)
            {
                name = "SQRT(R)";
            }
            else
            {
                name = "???????";
            }

        }

        public static int function_3d_num()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FUNCTION_3D_NUM returns the number of 3D functions.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 April 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, int FUNCTION_3D_NUM, the number of 3D functions.
            //
        {
            int value = 23;

            return value;
        }

        public static double function_nd(int function_nd_index, int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FUNCTION_ND evaluates the current ND function.
            //
            //  Discussion:
            //
            //    This routine assumes that the global variable FUNCTION_ND_INDEX has been
            //    set, and is accessible.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of variables.
            //
            //    Input, double X[N], the value of the variables.
            //
            //    Output, double FUNCTION_ND, the value of the function.
            //
        {
            int i;
            double temp;
            double value;

            if (function_nd_index == 1)
            {
                value = 1.0;
            }
            else if (function_nd_index == 2)
            {
                value = x[0];
            }
            else if (function_nd_index == 3)
            {
                value = Math.Pow(x[0], 2);
            }
            else if (function_nd_index == 4)
            {
                value = Math.Pow(x[0], 3);
            }
            else if (function_nd_index == 5)
            {
                value = Math.Pow(x[0], 4);
            }
            else if (function_nd_index == 6)
            {
                value = Math.Pow(x[0], 5);
            }
            else if (function_nd_index == 7)
            {
                value = Math.Pow(x[0], 6);
            }
            else if (function_nd_index == 8)
            {
                temp = 0.0;
                for (i = 0; i < n; i++)
                {
                    temp = temp + x[i] * x[i];
                }

                value = Math.Sqrt(temp);
            }
            else if (function_nd_index == 9)
            {
                value = Math.Sin(x[0]);
            }
            else if (function_nd_index == 10)
            {
                value = Math.Exp(x[0]);
            }
            else if (function_nd_index == 11)
            {
                temp = 0.0;
                for (i = 0; i < n; i++)
                {
                    temp = temp + x[i] * x[i];
                }

                value = 1.0 / (1.0 + Math.Sqrt(temp));
            }
            else if (function_nd_index == 12)
            {
                temp = 0.0;
                for (i = 0; i < n; i++)
                {
                    temp = temp + x[i] * x[i];
                }

                value = Math.Sqrt(Math.Sqrt(temp));
            }
            else
            {
                value = 0.0;
            }

            return value;
        }

        public static void function_nd_name(int function_nd_index, ref string name)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FUNCTION_ND_NAME returns the name of the current ND function.
            //
            //  Discussion:
            //
            //    This routine assumes that the global variable FUNCTION_ND_INDEX has been
            //    set, and is accessible.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 March 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, char NAME[8], the name of the current ND function.
            //
        {
            if (function_nd_index == 1)
            {
                name = "      1";
            }
            else if (function_nd_index == 2)
            {
                name = "      X";
            }
            else if (function_nd_index == 3)
            {
                name = "    X^2";
            }
            else if (function_nd_index == 4)
            {
                name = "    X^3";
            }
            else if (function_nd_index == 5)
            {
                name = "    X^4";
            }
            else if (function_nd_index == 6)
            {
                name = "    X^5";
            }
            else if (function_nd_index == 7)
            {
                name = "    X^6";
            }
            else if (function_nd_index == 8)
            {
                name = "      R";
            }
            else if (function_nd_index == 9)
            {
                name = " SIN(X)";
            }
            else if (function_nd_index == 10)
            {
                name = " EXP(X)";
            }
            else if (function_nd_index == 11)
            {
                name = "1/(1+R)";
            }
            else if (function_nd_index == 12)
            {
                name = "SQRT(R)";
            }
            else
            {
                name = "???????";
            }

        }

        public static int function_nd_num()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FUNCTION_ND_NUM returns the number of ND functions.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 April 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, int FUNCTION_ND_NUM, the number of ND functions.
            //
        {
            int value = 12;

            return value;
        }

        public static double f_1_2d(double x, double y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    F_1_2D evaluates the function 1 in 2D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 May 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, Y, the arguments.
            //
            //    Output, double F_1_2D, the value of the function.
            //
        {
            double value;

            value = 1.0;

            return value;
        }

        public static double f_x_2d(double x, double y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    F_X_2D evaluates the function X in 2D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 May 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, Y, the arguments.
            //
            //    Output, double F_X_2D, the value of the function.
            //
        {
            double value;

            value = x;

            return value;
        }

        public static double f_r_2d(double x, double y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    F_R_2D evaluates the function Math.Sqrt ( X**2 + Y**2 ) in 2D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 May 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, Y, the arguments.
            //
            //    Output, double F_R_2D, the value of the function.
            //
        {
            double value;

            value = Math.Sqrt(x * x + y * y);

            return value;
        }


    }
}