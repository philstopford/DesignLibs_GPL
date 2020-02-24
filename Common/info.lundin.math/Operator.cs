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
    /// Class Operator, represents an Operator by holding information about it's symbol
    /// the number of arguments it takes, the operator precedence and the delegate to use
    /// for evaluating.
    /// </summary>
    [Serializable]
    public class Operator
    {
        private String op = ""; // the string operator 
        private int args = 0; // the number of arguments this operator takes
        private int prec = System.Int32.MaxValue; // the precedence this operator has
        private Evaluator evaluator;

        /// <summary>
        /// Delegate for evaluating this operator
        /// </summary>
        /// <param name="parser">reference to parser to use for evaluating arguments</param>
        /// <param name="arg1">first argument</param>
        /// <param name="arg2">second argument</param>
        /// <returns>result after operator evaluation</returns>
        public delegate double Evaluator(ExpressionParser parser, Node arg1, Node arg2);

        /// <summary>
        /// Creates an Operator
        /// </summary>
        /// <param name="_operator">string operator symbol</param>
        /// <param name="arguments">number of arguments</param>
        /// <param name="precedence">precedence relative to other operators</param>
        /// <param name="eval">delegate to use for evaluating operator</param>
        public Operator(String _operator, int arguments, int precedence, Evaluator eval)
        {
            this.op = _operator;
            this.args = arguments;
            this.prec = precedence;
            this.evaluator = eval;
        }

        /// <summary>
        /// Provides access to the delegate for evaluating the operator
        /// </summary>
        public Evaluator Eval { get { return evaluator; } }

        /// <summary>
        /// Returns the precedence for this Operator.
        /// </summary>
        public int Precedence
        {
            get { return (this.prec); }
        }

        /// <summary>
        /// Returns the string symbol of this Operator.
        /// </summary>
        public string Symbol
        {
            get { return (this.op); }
        }

        /// <summary>
        /// Returns the number of arguments this Operator can take.
        /// </summary>
        public int Arguments
        {
            get { return (this.args); }
        }

    } // End class Operator
}
