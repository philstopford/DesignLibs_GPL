/*
 * Author: Patrik Lundin, patrik@lundin.info
 * Web: http://www.lundin.info
 * 
 * Source code released under the Microsoft Public License (Ms-PL) 
 * http://www.microsoft.com/en-us/openness/licenses.aspx#MPL
*/
using System;

namespace info.lundin.math
{
    public static class MathExtra
    {
        /// <summary>
        /// Calculates the factorial.
        /// </summary>
        /// <param name="val">the value to calcualte the factorial of</param>
        /// <returns>the factorial</returns>
        public static double
        Fac(double val)
        {

            if (!IsInteger(val))
            {
                return Double.NaN;
            }
            else if (val < 0)
            {
                return Double.NaN;
            }
            else if (val <= 1)
            {
                return 1;
            }

            return (val * Fac(val - 1));
        }

        /// <summary>
        /// Calculates the semi factorial.
        /// </summary>
        /// <param name="val">the value to calcualte the semi factorial of</param>
        /// <returns>the semi factorial</returns>
        public static double
        Sfac(double val)
        {
            if (!IsInteger(val))
            {
                return Double.NaN;
            }
            else if (val < 0)
            {
                return Double.NaN;
            }
            else if (val <= 1)
            {
                return 1;
            }

            return (val * Sfac(val - 2));
        }

        /// <summary>
        /// Returns the decimal part of the value
        /// </summary>
        /// <param name="val">the value to calculate the fpart for</param>
        /// <returns>the decimal part of the value</returns>
        public static double
        Fpart(double val)
        {
            if (val >= 0)
            {
                return (val - Math.Floor(val));
            }
            else
            {
                return (val - Math.Ceiling(val));
            }
        }


        /// <summary>
        /// Checks to see if the double value a can be considered to be a mathematical integer.</summary>
        /// <param name="a">the double value to check</param>
        /// <returns>true if the double value is an integer, false otherwise.</returns>
        private static bool
        IsInteger(double a)
        {
            return ((a - (int)a) == 0.0);
        }

        /// <summary>
        /// Checks to see if the int value a can be considered to be even. </summary>
        /// <param name="a">the int value to check</param>
        /// <returns>true if the int value is even, false otherwise.</returns>
        private static bool
        IsEven(int a)
        {
            return (IsInteger(a / 2));
        }
    }
}
