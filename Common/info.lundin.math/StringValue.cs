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
    /// Encapsulates a string value (nested expression) in the parser.
    /// </summary>
    public class StringValue : Value
    {
        /// <summary>
        /// Constructs a value of type ValueType.String
        /// </summary>
        public StringValue()
        {
            Type = ValueType.String;
        }

        /// <summary>
        /// String value
        /// </summary>
        public string Value { get; set; }

        /// <summary>
        /// Returns the string value
        /// </summary>
        /// <param name="format">unused parameter</param>
        /// <returns>the string value</returns>
        public override string ToString(IFormatProvider format)
        {
            return Value;
        }

        /// <summary>
        /// Returns the string value
        /// </summary>
        /// <returns>the string value</returns>
        public override string ToString()
        {
            return ToString(CultureInfo.InvariantCulture);
        }

        /// <summary>
        /// Returns the parsed double value if the string can
        /// be parsed to a double
        /// </summary>
        /// <param name="format">format to use when parsing string</param>
        /// <returns>value as double</returns>
        public override double ToDouble(IFormatProvider format)
        {
            return Double.Parse(Value, format);
        }

        /// <summary>
        /// Returns the parsed double value if the string can
        /// be parsed to a double
        /// </summary>
        /// <returns>value as double</returns>
        public override double ToDouble()
        {
            return ToDouble(CultureInfo.InvariantCulture);
        }

        /// <summary>
        /// Sets the value
        /// </summary>
        /// <param name="value">value to set</param>
        public override void SetValue(string value)
        {
            SetValue(value, CultureInfo.InvariantCulture);
        }

        /// <summary>
        /// Sets the value
        /// </summary>
        /// <param name="value">value to set</param>
        public override void SetValue(double value)
        {
            SetValue(value, CultureInfo.InvariantCulture);
        }

        /// <summary>
        /// Sets the value
        /// </summary>
        /// <param name="value">value to set</param>
        /// <param name="format">unused parameter</param>
        public override void SetValue(string value, IFormatProvider format)
        {
            Value = value;
        }

        /// <summary>
        /// Sets the value
        /// </summary>
        /// <param name="value">value to set</param>
        /// <param name="format">format provider to use when converting the value to string</param>
        public override void SetValue(double value, IFormatProvider format)
        {
            Value = value.ToString(format);
        }
    }
}
