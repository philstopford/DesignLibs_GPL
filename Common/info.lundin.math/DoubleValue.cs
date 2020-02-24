/*
 * Author: Patrik Lundin, patrik@lundin.info
 * Web: http://www.lundin.info
 * 
 * Source code released under the Microsoft Public License (Ms-PL) 
 * http://www.microsoft.com/en-us/openness/licenses.aspx#MPL
*/
using System;
using System.Globalization;

namespace info.lundin.math
{
    /// <summary>
    /// Encapsulates a double value for parser to use as a value when evaluating expressions
    /// </summary>
    public class DoubleValue : Value
    {
        /// <summary>
        /// Creates an instace of type ValueType.Constant
        /// </summary>
        public DoubleValue()
        {
            Type = ValueType.Constant;
        }

        /// <summary>
        /// The double value
        /// </summary>
        public double Value { get; set; }

        /// <summary>
        /// Returns the value as a string
        /// </summary>
        /// <param name="format">format provider to use</param>
        /// <returns>value as string</returns>
        public override string ToString(IFormatProvider format)
        {
            return Value.ToString(format);
        }

        /// <summary>
        /// Returns the value as a string
        /// </summary>
        /// <returns>value as string</returns>
        public override string ToString()
        {
            return ToString(CultureInfo.InvariantCulture);
        }

        /// <summary>
        /// Returns the double value
        /// </summary>
        /// <param name="format">unused parameter</param>
        /// <returns>value as double</returns>
        public override double ToDouble(IFormatProvider format)
        {
            return Value;
        }

        /// <summary>
        /// Returns the double value
        /// </summary>
        /// <returns>value as double</returns>
        public override double ToDouble()
        {
            return ToDouble(CultureInfo.InvariantCulture);
        }

        /// <summary>
        /// Sets the value 
        /// </summary>
        /// <param name="value">value to set as string</param>
        public override void SetValue(string value)
        {
            SetValue(value, CultureInfo.InvariantCulture);
        }

        /// <summary>
        /// Sets the value
        /// </summary>
        /// <param name="value">value to set as double</param>
        public override void SetValue(double value)
        {
            SetValue(value, CultureInfo.InvariantCulture);
        }

        /// <summary>
        /// Sets the value
        /// </summary>
        /// <param name="value">value to set as string</param>
        /// <param name="format">format provider to use when parsing the string to a double</param>
        public override void SetValue(string value, IFormatProvider format)
        {
            Value = Double.Parse(value, format);
        }

        /// <summary>
        /// Sets the value
        /// </summary>
        /// <param name="value">value as double</param>
        /// <param name="format">unused parameter</param>
        public override void SetValue(double value, IFormatProvider format)
        {
            Value = value;
        }
    }
}
