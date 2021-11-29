/*
 * Author: Patrik Lundin, patrik@lundin.info
 * Web: http://www.lundin.info
 * 
 * Source code released under the Microsoft Public License (Ms-PL) 
 * http://www.microsoft.com/en-us/openness/licenses.aspx#MPL
*/

using System;
using System.Collections.Generic;

namespace info.lundin.math;

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
            new( "^",  2, 3, (p, a, b) => Math.Pow(p.EvalTree(a), p.EvalTree(b))),
            new( "+",  2, 6, (p, a, b) => p.EvalTree(a) + p.EvalTree(b)),
            new( "-",  2, 6, (p, a, b) => p.EvalTree(a) - p.EvalTree(b)),
            new( "/",  2, 4, (p, a, b) => p.EvalTree(a) / p.EvalTree(b)),
            new( "*",  2, 4, (p, a, b) => p.EvalTree(a) * p.EvalTree(b)),
            new( "cos",    1, 2, (p, a, b) => Math.Cos(p.EvalTree(a))),
            new( "sin",    1, 2, (p, a, b) => Math.Sin(p.EvalTree(a))),
            new( "exp",    1, 2, (p, a, b) => Math.Exp(p.EvalTree(a))),
            new( "ln", 1, 2, (p, a, b) => Math.Log(p.EvalTree(a))),
            new( "tan",    1, 2, (p, a, b) => Math.Tan(p.EvalTree(a))),
            new( "acos",   1, 2, (p, a, b) => Math.Acos(p.EvalTree(a))),
            new( "asin",   1, 2, (p, a, b) => Math.Asin(p.EvalTree(a))),
            new( "atan",   1, 2, (p, a, b) => Math.Atan(p.EvalTree(a))),
            new( "cosh",   1, 2, (p, a, b) => Math.Cosh(p.EvalTree(a))),
            new( "sinh",   1, 2, (p, a, b) => Math.Sinh(p.EvalTree(a))),
            new( "tanh",   1, 2, (p, a, b) => Math.Tanh(p.EvalTree(a))),
            new( "sqrt",   1, 2, (p, a, b) => Math.Sqrt(p.EvalTree(a))),
            new( "cotan",1, 2, (p, a, b) => 1 / Math.Tan(p.EvalTree(a))),
            new( "fpart",1, 2, (p, a, b) => MathExtra.Fpart(p.EvalTree(a))),
            new( "acotan",1, 2, (p, a, b) => Math.PI / 2 - Math.Atan(p.EvalTree(a))),
            new( "round", 1, 2, (p, a, b) => Math.Round(p.EvalTree(a))),
            new( "ceil",  1, 2, (p, a, b) => Math.Ceiling(p.EvalTree(a))),
            new( "floor",1, 2, (p, a, b) => Math.Floor(p.EvalTree(a))),
            new( "fac",    1, 2, (p, a, b) => MathExtra.Fac(p.EvalTree(a))),
            new( "sfac",   1, 2, (p, a, b) => MathExtra.Sfac(p.EvalTree(a))),
            new( "abs",    1, 2, (p, a, b) => Math.Abs(p.EvalTree(a))),
            new( "log",    1, 2, (p, a, b) => Math.Log10(p.EvalTree(a))),
            new( "%",  2, 4, (p, a, b) => p.EvalTree(a) % p.EvalTree(b)),
            new( ">",  2, 7, (p, a, b) => p.EvalTree(a) > p.EvalTree(b) ? 1 : 0),
            new( "<",  2, 7, (p, a, b) => p.EvalTree(a) < p.EvalTree(b) ? 1 : 0),
            new( "&&", 2, 8, (p, a, b) => Math.Abs(p.EvalTree(a) - 1) <= double.Epsilon && Math.Abs(p.EvalTree(b) - 1) <= double.Epsilon ? 1 : 0),
            new( "==", 2, 7, (p, a, b) => Math.Abs(p.EvalTree(a) - p.EvalTree(b)) <= double.Epsilon ? 1 : 0),
            new( "!=", 2, 7, (p, a, b) => Math.Abs(p.EvalTree(a) - p.EvalTree(b)) > double.Epsilon ? 1 : 0),
            new( "||", 2, 9, (p, a, b) => Math.Abs(p.EvalTree(a) - 1) <= double.Epsilon || Math.Abs(p.EvalTree(b) - 1) <= double.Epsilon ? 1 : 0),
            new( "!",  1, 1, (p, a, b) => !(Math.Abs(p.EvalTree(a) - 1) <= double.Epsilon) ? 1 : 0),
            new( ">=", 2, 7, (p, a, b) => p.EvalTree(a) >= p.EvalTree(b) ? 1 : 0),
            new( "<=" ,    2, 7, (p, a, b) => p.EvalTree(a) <= p.EvalTree(b) ? 1 : 0)
        };
    }

    /// <summary>
    /// List of operators
    /// </summary>
    public IList<Operator> Operators => operators;
}