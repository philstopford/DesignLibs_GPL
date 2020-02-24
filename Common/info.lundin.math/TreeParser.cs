/*
 * Author: Patrik Lundin, patrik@lundin.info
 * Web: http://www.lundin.info
 * 
 * Source code released under the Microsoft Public License (Ms-PL) 
 * http://www.microsoft.com/en-us/openness/licenses.aspx#MPL
*/
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;

namespace info.lundin.math
{

    public class TreeParser
    {
        private IDictionary<string, Operator> operators;
        private IDictionary<string, double> constants;
        private CultureInfo culture;

        private int maxoplength;

        bool bRequireParentheses;
        bool bImplicitMultiplication;

        /// <summary>
        /// Constructs a TreeParser instance
        /// </summary>
        /// <param name="operators">operators to use when creating tree representation of an expression</param>
        /// <param name="constants">constants to use when creating tree representation of an expression</param>
        public TreeParser(IDictionary<string, Operator> operators, IDictionary<string, double> constants)
        {
            this.operators = operators;
            this.constants = constants;

            foreach (var pair in operators)
            {
                string symbol = pair.Value.Symbol;
                if (symbol.Length > maxoplength)
                {
                    maxoplength = symbol.Length;
                }
            }

            culture = CultureInfo.InvariantCulture;
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
                return bRequireParentheses;
            }
            set
            {
                bRequireParentheses = value;
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
                return bImplicitMultiplication;
            }
            set
            {
                bImplicitMultiplication = value;
            }
        }

        /// <summary>
        /// Gets or sets the culture to use
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
            }
        }

        /// <summary>
        /// Checks the String expression to see if the syntax is valid.
        /// this method doesn't return anything, instead it throws an Exception
        /// if the syntax is invalid.
        ///	
        /// Examples of invalid syntax can be non matching paranthesis, non valid symbols appearing
        /// or a variable or operator name is invalid in the expression.
        /// </summary>
        /// <param name="exp">the string expression to check, infix notation.</param>
        /// 
        /// <remarks>This validates some syntax errors such as unbalanced paranthesis and invalid characters.
        /// Exceptions can also be thrown at evaluation time.</remarks>
        private void
        SyntaxCheck(string exp)
        {
            int i = 0, oplen = 0;
            string op = null;
            string nop = null;

            // Check if all paranthesis match in expression
            if (!matchParant(exp))
            {
                throw new ParserException("Unbalanced parenthesis");
            }

            int l = exp.Length;

            // Go through expression and validate syntax
            while (i < l)
            {
                try
                {

                    if ((op = getOp(exp, i)) != null)
                    {
                        // Found operator at position, check syntax for operators

                        oplen = op.Length;
                        i += oplen;

                        // If it's a function and we are missing parentheses around arguments and bRequireParantheses is true it is an error.
                        // Note that this only checks opening paranthesis, we checked whole expression for balanced paranthesis
                        // earlier but not for each individual function.
                        if (bRequireParentheses && !isTwoArgOp(op) && op != "!" && exp[i] != '(')
                        {
                            throw new ParserException("Paranthesis required for arguments -> " + exp.Substring(i - oplen));
                        }

                        // If we have an operator immediately following a function and it's not unary + or - then it's an error
                        nop = getOp(exp, i);
                        if (nop != null && isTwoArgOp(nop) && !(nop.Equals("+") || nop.Equals("-")))
                        {
                            throw new ParserException("Syntax error near -> " + exp.Substring(i - oplen));
                        }
                    }
                    else if (!isAlpha(exp[i]) && !isConstant(exp[i]) && !isAllowedSym(exp[i]))
                    {
                        // This cannot be a valid character, throw exception
                        throw new ParserException("Syntax error near -> " + exp.Substring(i));
                    }
                    else
                    {
                        // Count forward
                        i++;
                    }
                }
                catch (IndexOutOfRangeException)
                {
                    // Might happen when we are checking syntax, just move forward
                    i++;
                }
            }

            return;
        }



        /// <summary>
        /// Inserts the multiplication operator where needed.
        /// This method adds limited support for implicit multiplication in expressions.
        /// </summary>
        /// <remarks>	
        /// Implicit multiplication is supported in these type cases:
        ///
        /// case: variable jp one-arg-op , xcos(x)
        /// case: const jp variable or one-arg-op, 2x, 2tan(x)
        /// case: "const jp ( expr )" , 2(3+x)
        /// case: ( expr ) jp variable or one-arg-op , (2-x)x , (2-x)sin(x)
        /// case: var jp  ( expr ) , x(x+1) , x(1-sin(x))
        ///
        /// Note that this also puts extra limitations on variable names, they cannot
        /// contain digits within them or at the beginning, only at the end.
        /// 
        /// If the ImplicitMultiplication property is set to false this method just returns
        /// the same expression sent to it without modification, effectively disabling support
        /// for implicit multiplication. This hoewever does not change the variable naming requirements.
        /// </remarks>
        /// <param name="exp">the infix string expression to process</param>
        /// <returns>the processed infix expression</returns>
        private string
        ParseImplicit(string exp)
        {
            // Return expression unchanged if ImplicitMultiplication property is false.
            if (!bImplicitMultiplication) return exp;

            int i = 0, p = 0;
            String tmp = null;
            StringBuilder str = new StringBuilder(exp);

            int l = exp.Length;

            while (i < l)
            {
                try
                {
                    if ((tmp = getOp(exp, i)) != null && !isTwoArgOp(tmp) && (i - 1 >= 0 && isAlpha(exp[i - 1])))
                    {
                        // case: variable jp one-arg-op , xcos(x)
                        str.Insert(i + p, "*");
                        p++;
                    }
                    else if (isAlpha(exp[i]) && (i - 1 >= 0 && isConstant(exp[i - 1])))
                    {
                        // case: const jp variable or one-arg-op, 2x, 2tan(x)
                        str.Insert(i + p, "*");
                        p++;
                    }
                    else if (exp[i] == '(' && (i - 1 >= 0 && isConstant(exp[i - 1])))
                    {
                        // case: "const jp ( expr )" , 2(3+x)
                        str.Insert(i + p, "*");
                        p++;
                    }
                    else if (isAlpha(exp[i]) && (i - 1 >= 0 && exp[i - 1] == ')'))
                    {
                        // case: ( expr ) jp variable or one-arg-op , (2-x)x , (2-x)sin(x)
                        str.Insert(i + p, "*");
                        p++;
                    }
                    else if (exp[i] == '(' && (i - 1 >= 0 && exp[i - 1] == ')'))
                    {
                        // case: ( expr ) jp  ( expr ) , (2-x)(x+1) , sin(x)(2-x) 
                        str.Insert(i + p, "*");
                        p++;
                    }
                    else if (exp[i] == '(' && (i - 1 >= 0 && isAlpha(exp[i - 1])) && backTrack(exp.Substring(0, i)) == null)
                    {
                        // case: var jp  ( expr ) , x(x+1) , x(1-sin(x))
                        str.Insert(i + p, "*");
                        p++;
                    }
                }
                catch
                {
                    // checking indexes properly so try/catch should not be needed but lets keep it for now
                }

                if (tmp != null)
                {
                    i += tmp.Length;
                }
                else
                {
                    i++;
                }

                tmp = null;
            }

            return str.ToString();
        }


        /// <summary>
        /// Adds support for "scientific notation" by replacing the E operator with *10^
        /// </summary>
        /// <remarks>
        /// For example the value 1E-3 would be changed to 1*10^-3 which the parser will treat
        /// as a normal expression.
        /// </remarks>
        /// <param name="exp">the infix string expression to process</param>
        /// <returns>the processed infix expression</returns>
        private string
        ParseE(string exp)
        {

            int i, p, len;

            StringBuilder newstr = new StringBuilder(exp);

            i = p = 0;
            len = exp.Length;

            while (i < len)
            {
                try
                {
                    if (exp[i] == 'e' && (i - 1 >= 0 && Char.IsDigit(exp[i - 1])))
                    {
                        if ((i + 1 < len && Char.IsDigit(exp[i + 1])) || (i + 2 < len && ((exp[i + 1] == '-' || exp[i + 1] == '+') && Char.IsDigit(exp[i + 2]))))
                        {
                            // replace the 'e'
                            newstr[i + p] = '*';
                            // insert the rest
                            newstr.Insert(i + p + 1, "10^");
                            p = p + 3; // buffer growed by 3 chars
                        }
                    }
                }
                catch
                {
                    // checking indexes properly so try/catch should not be needed but lets keep it for now
                }

                i++;
            }

            return newstr.ToString();
        }

        /// <summary>
        /// Parses an infix String expression and creates a parse tree of Node's.
        /// </summary>
        /// <remarks>
        /// This is the heart of the parser, it takes a normal infix expression and creates
        /// a tree datastructure we can easily recurse when evaluating.
        /// </remarks>
        /// <param name="exp">the infix string expression to process</param>
        /// <returns>A tree datastructure of Node objects representing the expression</returns>
        private Node
        ParseInfix(string exp)
        {
            int i, ma, len;
            string farg, sarg, fop;
            Node tree = null;

            farg = sarg = fop = "";
            ma = i = 0;

            len = exp.Length;

            if (len == 0)
            {
                throw new ParserException("Wrong number of arguments to operator");
            }
            else if (exp[0] == '(' && ((ma = match(exp, 0)) == (len - 1)))
            {
                return (ParseInfix(exp.Substring(1, ma - 1)));
            }
            else if (isVariable(exp))
            {
                // If built in constant put in value otherwise the variable
                if (constants.ContainsKey(exp)) return new Node(constants[exp]);
                else return new Node(exp);
            }
            else if (isConstant(exp))
            {
                try
                {
                    return (new Node(Double.Parse(exp, NumberStyles.Any, culture)));
                }
                catch (FormatException)
                {
                    throw new ParserException("Syntax error-> " + exp + " (not using regional decimal separator?)");
                }
            }

            while (i < len)
            {
                if ((fop = getOp(exp, i)) == null)
                {
                    farg = arg(null, exp, i);
                    fop = getOp(exp, i + farg.Length);

                    if (fop == null) throw new Exception("Missing operator");

                    if (isTwoArgOp(fop))
                    {
                        sarg = arg(fop, exp, i + farg.Length + fop.Length);
                        if (sarg.Equals("")) throw new Exception("Wrong number of arguments to operator " + fop);
                        tree = new Node(operators[fop], ParseInfix(farg), ParseInfix(sarg));
                        i += farg.Length + fop.Length + sarg.Length;
                    }
                    else
                    {
                        if (farg.Equals("")) throw new Exception("Wrong number of arguments to operator " + fop);
                        tree = new Node(operators[fop], ParseInfix(farg));
                        i += farg.Length + fop.Length;
                    }
                }
                else
                {
                    if (isTwoArgOp(fop))
                    {
                        farg = arg(fop, exp, i + fop.Length);
                        if (farg.Equals("")) throw new Exception("Wrong number of arguments to operator " + fop);
                        if (tree == null)
                        {
                            if (fop.Equals("+") || fop.Equals("-"))
                            {
                                tree = new Node(0D);
                            }
                            else throw new Exception("Wrong number of arguments to operator " + fop);
                        }
                        tree = new Node(operators[fop], tree, ParseInfix(farg));
                        i += farg.Length + fop.Length;
                    }
                    else
                    {
                        farg = arg(fop, exp, i + fop.Length);
                        if (farg.Equals("")) throw new Exception("Wrong number of arguments to operator " + fop);
                        tree = new Node(operators[fop], ParseInfix(farg));
                        i += farg.Length + fop.Length;
                    }
                }
            }

            return tree;
        }

        public Expression Parse(string expression)
        {
            string tmp = skipSpaces(expression.ToLower());

            SyntaxCheck(tmp);

            Node tree = ParseInfix(ParseImplicit(ParseE(tmp)));

            return new Expression(tree);
        }

        /// <summary>Matches all paranthesis and returns true if they all match or false if they do not.</summary>
        /// <param name="exp">expression to check, infix notation</param>
        /// <returns>true if ok false otherwise</returns>
        private bool
        matchParant(String exp)
        {
            int count = 0;
            int i = 0;

            int l = exp.Length;

            for (i = 0; i < l; i++)
            {
                if (exp[i] == '(')
                {
                    count++;
                }
                else if (exp[i] == ')')
                {
                    count--;
                }
            }

            return (count == 0);
        }

        /// <summary>Checks if the character is alphabetic.</summary>
        /// <param name="ch">Character to check</param>
        /// <returns>true or false</returns>
        private bool
        isAlpha(char ch)
        {
            return ((ch >= 'a' && ch <= 'z') || (ch >= 'A' && ch <= 'Z'));
        }


        /// <summary>Checks if the string can be considered to be a valid variable name.</summary>
        /// <param name="str">The String to check</param>
        /// <returns>true or false</returns>
        private bool
        isVariable(String str)
        {
            int i = 0;
            int len = str.Length;

            if (isConstant(str)) return false;

            for (i = 0; i < len; i++)
            {
                if (getOp(str, i) != null || isAllowedSym(str[i]))
                {
                    return false;
                }
            }

            return true;
        }


        /// <summary>Checks if the character is a digit</summary>
        /// <param name="ch">Character to check</param>
        /// <returns>true or false</returns>
        private bool
        isConstant(char ch)
        {
            return (Char.IsDigit(ch));
        }



        /// <summary>Checks to se if a string is numeric</summary>
        /// <param name="exp">String to check</param>
        /// <returns>true if the string was numeric, false otherwise</returns>
        private bool
        isConstant(String exp)
        {
            double val = 0D;
            bool ok = Double.TryParse(exp, NumberStyles.Any, culture, out val);
            return (ok && !Double.IsNaN(val));
        }

        /// <summary>
        /// Checks to see if the string is the name of a acceptable operator.
        /// </summary>
        /// <param name="str">The string to check</param>
        /// <returns>true if it is an acceptable operator, false otherwise.</returns>
        private bool
        isOperator(String str)
        {
            return (operators.ContainsKey(str));
        }



        /// <summary>
        /// Checks to see if the operator name represented by str takes two arguments.
        /// </summary>
        /// <param name="str">The string to check</param>
        /// <returns>true if the operator takes two arguments, false otherwise.</returns>
        private bool
        isTwoArgOp(String str)
        {
            if (str == null) return false;
            Object o = operators[str];
            if (o == null) return false;
            return (((Operator)o).Arguments == 2);
        }

        /// <summary>
        /// Checks to see if the character is a valid symbol for this parser.
        /// </summary>
        /// <param name="s">the character to check</param>
        /// <returns>true if the char is valid, false otherwise.</returns>
        private bool
        isAllowedSym(char s)
        {
            return (
                s == ',' || s == '.' || s == ')' || s == '(' || s == '>' || s == '<' || s == '&' || s == '=' || s == '|'
                || culture.NumberFormat.CurrencyGroupSeparator.ToCharArray().Contains(s)
                || culture.NumberFormat.CurrencyDecimalSeparator.ToCharArray().Contains(s));
        }

        /// <summary>
        /// Parses out spaces from a string
        /// </summary>
        /// <param name="str">The string to process</param>
        /// <returns>A copy of the string stripped of all spaces</returns>
        private String
        skipSpaces(String str)
        {
            int i = 0;
            int len = str.Length;
            StringBuilder nstr = new StringBuilder(len);

            while (i < len)
            {
                if (!Char.IsWhiteSpace(str[i]))
                {
                    nstr.Append(str[i]);
                }
                i++;
            }

            return nstr.ToString();
        }

        /// <summary>
        /// Matches an opening left paranthesis.
        /// </summary>
        /// <param name="exp">the string to search in</param>
        /// <param name="index">the index of the opening left paranthesis</param>
        /// <returns>the index of the matching closing right paranthesis</returns>
        private int
        match(String exp, int index)
        {
            int len = exp.Length;
            int i = index;
            int count = 0;

            while (i < len)
            {
                if (exp[i] == '(')
                {
                    count++;
                }
                else if (exp[i] == ')')
                {
                    count--;
                }

                if (count == 0) return i;

                i++;
            }

            return index;
        }

        /// <summary>
        /// Parses out an operator from an infix string expression.
        /// </summary>
        /// <param name="exp">the infix string expression to look in</param>
        /// <param name="index">the index to start searching from</param>
        /// <returns>the operator if any or null.</returns>
        private String
        getOp(String exp, int index)
        {
            String tmp;
            int i = 0;
            int len = exp.Length;

            for (i = 0; i < maxoplength; i++)
            {
                if (index >= 0 && (index + maxoplength - i) <= len)
                {
                    tmp = exp.Substring(index, maxoplength - i);
                    if (isOperator(tmp))
                    {
                        return (tmp);
                    }
                }
            }

            return null;
        }

        /// <summary>
        /// Parses the infix expression for arguments to the specified operator.
        /// </summary>
        /// <param name="_operator">the operator we are interested in</param>
        /// <param name="exp">the infix string expression</param>
        /// <param name="index">the index to start the search from</param>
        /// <returns>the argument to the operator</returns>
        private String
        arg(String _operator, String exp, int index)
        {
            int ma, i, prec = -1;
            int len = exp.Length;
            String op = null;

            StringBuilder str = new StringBuilder();

            i = index;
            ma = 0;

            if (_operator == null)
            {
                prec = -1;
            }
            else
            {
                prec = (operators[_operator]).Precedence;
            }

            while (i < len)
            {

                if (exp[i] == '(')
                {
                    ma = match(exp, i);
                    str.Append(exp.Substring(i, ma + 1 - i));
                    i = ma + 1;
                }
                else if ((op = getOp(exp, i)) != null)
                {
                    // (_operator != null && _operator.Equals("&&") && op.Equals("||") ) || 
                    if (str.Length != 0 && !isTwoArgOp(backTrack(str.ToString())) && (operators[op]).Precedence >= prec)
                    {
                        return str.ToString();
                    }
                    str.Append(op);
                    i += op.Length;
                }
                else
                {
                    str.Append(exp[i]);
                    i++;
                }
            }

            return str.ToString();
        }



        /// <summary>
        /// Returns an operator at the end of the String str if present.
        /// </summary>
        /// <remarks>
        /// Used when parsing for arguments, the purpose is to recognize
        /// expressions like for example 10^-1
        /// </remarks>
        /// <param name="str">part of infix string expression to search</param>
        /// <returns>the operator if found or null otherwise</returns>
        private String
        backTrack(String str)
        {
            int i = 0;
            int len = str.Length;
            String op = null;

            try
            {
                for (i = 0; i <= maxoplength; i++)
                {
                    if ((op = getOp(str, (len - 1 - maxoplength + i))) != null
                    && (len - maxoplength - 1 + i + op.Length) == len)
                    {
                        return op;
                    }
                }
            }
            catch { }

            return null;
        }

    } // End class TreeParser

} // End namespace info.lundin.math

