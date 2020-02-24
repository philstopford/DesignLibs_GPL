/*
 * Author: Patrik Lundin, patrik@lundin.info
 * Web: http://www.lundin.info
 * 
 * Source code released under the Microsoft Public License (Ms-PL) 
 * http://www.microsoft.com/en-us/openness/licenses.aspx#MPL
*/
namespace info.lundin.math
{

    using System;
    using System.Collections.Generic;
    using System.Globalization;

    /// <summary>
    ///	Class ExpressionParser, this class evaluates a mathematical expression given 
    ///	as a string to a double value.	
    /// </summary>
    public class ExpressionParser
    {
        private ValuesDictionary htbl;
        private ExpressionDictionary cachedExpressions;
        private TreeParser treeParser;

        private IDictionary<string, Operator> ops;
        private IDictionary<string, double> spconst;

        private CultureInfo culture;

        private int maxoplength;

        /// <summary>
		/// Default constructor, creates an ExpressionParser object
		/// </summary>
		public ExpressionParser()
        {
            // Cache for values and expressions
            htbl = new ValuesDictionary();
            cachedExpressions = new ExpressionDictionary();

            // Add all valid operators.
            DefaultOperators operators = new DefaultOperators();
            ops = new Dictionary<string, Operator>();

            foreach (Operator op in operators.Operators)
            {
                string symbol = op.Symbol;
                if (symbol.Length > maxoplength)
                {
                    maxoplength = symbol.Length;
                }

                ops.Add(symbol, op);
            }

            // Constants
            spconst = new Dictionary<string, double>();
            spconst.Add("euler", Math.E);
            spconst.Add("pi", Math.PI);
            spconst.Add("nan", double.NaN);
            spconst.Add("infinity", double.PositiveInfinity);
            spconst.Add("true", 1D);
            spconst.Add("false", 0D);

            treeParser = new TreeParser(ops, spconst);
            treeParser.ImplicitMultiplication = true;
            treeParser.RequireParentheses = true;

            Culture = CultureInfo.InvariantCulture;
        }

        /// <summary>
        /// Enables or disables the requirement to put all function arguments within parantheses.
        /// 
        /// Default is true making it required for all function arguments to be enclosed within () for example sin(x+5)
        /// failure to do so will generate an exception stating that parantheses are required.
        /// 
        /// Setting this property to false disables this requirement making expressions such as for example sinx or sin5 allowed.
        /// </summary>
        public bool RequireParentheses
        {
            get
            {
                return treeParser.RequireParentheses;
            }
            set
            {
                treeParser.RequireParentheses = value;
            }
        }

        /// <summary>
        /// Enables or disables the support for implicit multiplication.
        /// 
        /// Default is true making implicit multiplication allowed, for example 2x or 3sin(2x)
        /// setting this property to false disables support for implicit multiplication making the * operator required
        /// for example 2*x or 3*sin(2*x)
        /// 
        /// If implicit multiplication is disabled (false) parsing an expression that does not explicitly use the * operator may
        /// throw syntax errors with various error messages.
        /// </summary>
        public bool ImplicitMultiplication
        {
            get
            {
                return treeParser.ImplicitMultiplication;
            }
            set
            {
                treeParser.ImplicitMultiplication = value;
            }
        }

        /// <summary>
        /// Provides access to the ValueDictionary for adding values to the parser
        /// for example parser.Values.Add("x", 5)
        /// </summary>
        public ValuesDictionary Values
        {
            get
            {
                return htbl;
            }
        }

        /// <summary>
        /// Provides access to the ExpressionDictionary and the parsed expressions
        /// </summary>
        public ExpressionDictionary Expressions
        {
            get { return cachedExpressions; }
        }

        /// <summary>
        /// Gets or sets the culture to use when parsing strings that contain double values.
        /// </summary>
        public CultureInfo Culture
        {
            get
            {
                return culture;
            }

            set
            {
                if (value == null)
                    throw new ArgumentNullException("Culture cannot be null");

                // Cannot allow division operator as decimal separator (fa, fa-IR)
                if (value.NumberFormat.NumberDecimalSeparator == "/")
                {
                    throw new ArgumentOutOfRangeException(
                        String.Format("Unsupported decimal separator / culture {0}", value.Name));
                }

                // Cannot allow same separators
                if (value.NumberFormat.CurrencyGroupSeparator
                    == value.NumberFormat.NumberDecimalSeparator)
                {
                    throw new ArgumentOutOfRangeException(
                        String.Format("Same decimal and group separator is unsupported. Culture {0}", value.Name));
                }

                culture = value;
                treeParser.Culture = culture;
            }
        }

        /// <summary>
        /// Evaluates the mathematical infix expression using the values that have been
        /// added to the parser.
        /// </summary>
        /// <remarks></remarks>
        /// <param name="exp">the infix string expression to parse and evaluate.</param>
        /// <returns>resulting double value after evaluating</returns>
        public double
        Parse(string exp)
        {
            double ans = 0d;
            Expression expression;

            if (String.IsNullOrWhiteSpace(exp))
                throw new ParserException("First argument to method Parse is null or empty");

            try
            {
                if (cachedExpressions.TryGetValue(exp, out expression))
                {
                    ans = EvalExpression(expression);
                }
                else
                {
                    expression = treeParser.Parse(exp);
                    ans = EvalExpression(expression);
                    cachedExpressions.Add(exp, expression);
                }

                return ans;
            }
            catch (Exception e)
            {
                throw new ParserException(e.Message);
            }
        }

        /// <summary>
        /// Traverses and evaluates the datastructure created by the parse method.
        /// </summary>
        /// <param name="expression">expression to evaluate</param>
        /// <returns>resulting double value</returns>
        public double
        EvalExpression(Expression expression)
        {
            return EvalTree(expression.ExpressionTree);
        }

        /// <summary>
        /// Traverses and evaluates the datastructure
        /// </summary>
        /// <param name="tree">tree as a structure of Node(s)</param>
        /// <returns>resulting double value</returns>
        internal double
        EvalTree(Node tree)
        {
            string tmp;

            if (tree.Type == NodeType.Value)
            {
                return (tree.Value.ToDouble());
            }
            else if (tree.Type == NodeType.Variable)
            {
                // get value associated with variable
                Value value = GetValue(tree);

                if (value.Type == ValueType.String
                    && cachedExpressions.ContainsKey(value.ToString(culture))) // cached expression
                {
                    return EvalExpression(cachedExpressions[value.ToString(culture)]);
                }
                else if (value.Type == ValueType.Constant) // constant value
                {
                    return value.ToDouble(culture);
                }
                else
                {
                    // apparently a nested expression, parse and cache
                    htbl[tree.Variable] = new StringValue() { Value = value.ToString(culture) };

                    tmp = value.ToString(culture);
                    Expression expression = treeParser.Parse(tmp);

                    cachedExpressions.Add(tmp, expression);

                    return EvalExpression(expression);
                }
            }

            return tree.Operator.Eval(this, tree.FirstArgument, tree.SecondArgument);
        }

        /// <summary>
        /// Retrieves the value associated with the variable in the node
        /// </summary>
        /// <param name="node">node of type NodeType.Variable</param>
        /// <returns>value associated with variable</returns>
        private Value GetValue(Node node)
        {
            Value value;

            if (node.Value != null
                && node.Value.Type != ValueType.Invalid)
            {
                // Value already in node
                value = node.Value;
            }
            else
            {
                if (!htbl.TryGetValue(node.Variable, out value))
                {
                    throw new ParserException("No value associated with " + node.Variable);
                }

                // Save value reference in the node for next time
                node.Value = value;
            }

            return value;
        }


    } // End class ExpressionParser

} // End namespace info.lundin.math