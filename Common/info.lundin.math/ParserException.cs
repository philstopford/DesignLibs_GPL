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
    /// <summary>
    /// Exception class for parser related exceptions
    /// </summary>
    public class ParserException : System.Exception
    {
        public ParserException()
            : base()
        {

        }

        public ParserException(string message)
            : base(message)
        {

        }

        public ParserException(string message, Exception inner)
            : base(message, inner)
        {

        }
    }
}
