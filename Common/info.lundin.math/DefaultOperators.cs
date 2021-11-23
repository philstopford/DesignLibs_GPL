/*
 * Author: Patrik Lundin, patrik@lundin.info
 * Web: http://www.lundin.info
 * 
 * Source code released under the Microsoft Public License (Ms-PL) 
 * http://www.microsoft.com/en-us/openness/licenses.aspx#MPL
*/
using System;
using System.Collections.Generic;

namespace info.lundin.math
{
    /// <summary>
    /// Provides operators for evaluation
    /// </summary>
    public class DefaultOperators
    {
        private IList<Operator> operators;

        public DefaultOperators()
        {
            operators = new List<Operator>
            {
                    new Operator( "^",  2, 3, (p, a, b) => { return Math.Pow(p.EvalTree(a), p.EvalTree(b));}),
                    new Operator( "+",  2, 6, (p, a, b) => { return p.EvalTree(a) + p.EvalTree(b);}),
                    new Operator( "-",  2, 6, (p, a, b) => { return p.EvalTree(a) - p.EvalTree(b);} ),
                    new Operator( "/",  2, 4, (p, a, b) => { return p.EvalTree(a) / p.EvalTree(b);} ),
                    new Operator( "*",  2, 4, (p, a, b) => { return p.EvalTree(a) * p.EvalTree(b);} ),
                    new Operator( "cos",    1, 2, (p, a, b) => { return Math.Cos(p.EvalTree(a));} ),
                    new Operator( "sin",    1, 2, (p, a, b) => { return Math.Sin(p.EvalTree(a));} ),
                    new Operator( "exp",    1, 2, (p, a, b) => { return Math.Exp(p.EvalTree(a));} ),
                    new Operator( "ln", 1, 2, (p, a, b) => { return Math.Log(p.EvalTree(a));} ),
                    new Operator( "tan",    1, 2, (p, a, b) => { return Math.Tan(p.EvalTree(a));} ),
                    new Operator( "acos",   1, 2, (p, a, b) => { return Math.Acos(p.EvalTree(a));} ),
                    new Operator( "asin",   1, 2, (p, a, b) => { return Math.Asin(p.EvalTree(a));} ),
                    new Operator( "atan",   1, 2, (p, a, b) => { return Math.Atan(p.EvalTree(a));} ),
                    new Operator( "cosh",   1, 2, (p, a, b) => { return Math.Cosh(p.EvalTree(a));} ),
                    new Operator( "sinh",   1, 2, (p, a, b) => { return Math.Sinh(p.EvalTree(a));} ),
                    new Operator( "tanh",   1, 2, (p, a, b) => { return Math.Tanh(p.EvalTree(a));} ),
                    new Operator( "sqrt",   1, 2, (p, a, b) => { return Math.Sqrt(p.EvalTree(a));} ),
                    new Operator( "cotan",1, 2, (p, a, b) => { return 1 / Math.Tan(p.EvalTree(a));} ),
                    new Operator( "fpart",1, 2, (p, a, b) => { return MathExtra.Fpart(p.EvalTree(a));}),
                    new Operator( "acotan",1, 2, (p, a, b) => { return Math.PI / 2 - Math.Atan(p.EvalTree(a));} ),
                    new Operator( "round", 1, 2, (p, a, b) => { return Math.Round(p.EvalTree(a));} ),
                    new Operator( "ceil",  1, 2, (p, a, b) => { return Math.Ceiling(p.EvalTree(a));} ),
                    new Operator( "floor",1, 2, (p, a, b) => { return Math.Floor(p.EvalTree(a));} ),
                    new Operator( "fac",    1, 2, (p, a, b) => { return MathExtra.Fac(p.EvalTree(a));}),
                    new Operator( "sfac",   1, 2, (p, a, b) => { return MathExtra.Sfac(p.EvalTree(a));}),
                    new Operator( "abs",    1, 2, (p, a, b) => { return Math.Abs(p.EvalTree(a));} ),
                    new Operator( "log",    1, 2, (p, a, b) => { return Math.Log10(p.EvalTree(a));}),
                    new Operator( "%",  2, 4, (p, a, b) => { return p.EvalTree(a) % p.EvalTree(b);} ),
                    new Operator( ">",  2, 7, (p, a, b) => { return p.EvalTree(a) > p.EvalTree(b) ? 1 : 0;} ),
                    new Operator( "<",  2, 7, (p, a, b) => { return p.EvalTree(a) < p.EvalTree(b) ? 1 : 0;} ),
                    new Operator( "&&", 2, 8, (p, a, b) => { return p.EvalTree(a) == 1 && p.EvalTree(b) == 1 ? 1 : 0;} ),
                    new Operator( "==", 2, 7, (p, a, b) => { return p.EvalTree(a) == p.EvalTree(b) ? 1 : 0;} ),
                    new Operator( "!=", 2, 7, (p, a, b) => { return p.EvalTree(a) != p.EvalTree(b) ? 1 : 0;} ),
                    new Operator( "||", 2, 9, (p, a, b) => { return p.EvalTree(a) == 1 || p.EvalTree(b) == 1 ? 1 : 0;} ),
                    new Operator( "!",  1, 1, (p, a, b) => { return !(p.EvalTree(a) == 1) ? 1 : 0;} ),
                    new Operator( ">=", 2, 7, (p, a, b) => { return p.EvalTree(a) >= p.EvalTree(b) ? 1 : 0;} ),
                    new Operator( "<=" ,    2, 7, (p, a, b) => { return p.EvalTree(a) <= p.EvalTree(b) ? 1 : 0;} )
            };
        }

        /// <summary>
        /// List of operators
        /// </summary>
        public IList<Operator> Operators
        {
            get
            {
                return operators;
            }
        }
    }
}
