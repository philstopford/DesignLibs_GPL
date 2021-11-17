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
            new( "^",  2, 3, (p, a, b) => { return Math.Pow(p.EvalTree(a), p.EvalTree(b));}),
            new( "+",  2, 6, (p, a, b) => { return p.EvalTree(a) + p.EvalTree(b);}),
            new( "-",  2, 6, (p, a, b) => { return p.EvalTree(a) - p.EvalTree(b);} ),
            new( "/",  2, 4, (p, a, b) => { return p.EvalTree(a) / p.EvalTree(b);} ),
            new( "*",  2, 4, (p, a, b) => { return p.EvalTree(a) * p.EvalTree(b);} ),
            new( "cos",    1, 2, (p, a, b) => { return Math.Cos(p.EvalTree(a));} ),
            new( "sin",    1, 2, (p, a, b) => { return Math.Sin(p.EvalTree(a));} ),
            new( "exp",    1, 2, (p, a, b) => { return Math.Exp(p.EvalTree(a));} ),
            new( "ln", 1, 2, (p, a, b) => { return Math.Log(p.EvalTree(a));} ),
            new( "tan",    1, 2, (p, a, b) => { return Math.Tan(p.EvalTree(a));} ),
            new( "acos",   1, 2, (p, a, b) => { return Math.Acos(p.EvalTree(a));} ),
            new( "asin",   1, 2, (p, a, b) => { return Math.Asin(p.EvalTree(a));} ),
            new( "atan",   1, 2, (p, a, b) => { return Math.Atan(p.EvalTree(a));} ),
            new( "cosh",   1, 2, (p, a, b) => { return Math.Cosh(p.EvalTree(a));} ),
            new( "sinh",   1, 2, (p, a, b) => { return Math.Sinh(p.EvalTree(a));} ),
            new( "tanh",   1, 2, (p, a, b) => { return Math.Tanh(p.EvalTree(a));} ),
            new( "sqrt",   1, 2, (p, a, b) => { return Math.Sqrt(p.EvalTree(a));} ),
            new( "cotan",1, 2, (p, a, b) => { return 1 / Math.Tan(p.EvalTree(a));} ),
            new( "fpart",1, 2, (p, a, b) => { return MathExtra.Fpart(p.EvalTree(a));}),
            new( "acotan",1, 2, (p, a, b) => { return Math.PI / 2 - Math.Atan(p.EvalTree(a));} ),
            new( "round", 1, 2, (p, a, b) => { return Math.Round(p.EvalTree(a));} ),
            new( "ceil",  1, 2, (p, a, b) => { return Math.Ceiling(p.EvalTree(a));} ),
            new( "floor",1, 2, (p, a, b) => { return Math.Floor(p.EvalTree(a));} ),
            new( "fac",    1, 2, (p, a, b) => { return MathExtra.Fac(p.EvalTree(a));}),
            new( "sfac",   1, 2, (p, a, b) => { return MathExtra.Sfac(p.EvalTree(a));}),
            new( "abs",    1, 2, (p, a, b) => { return Math.Abs(p.EvalTree(a));} ),
            new( "log",    1, 2, (p, a, b) => { return Math.Log10(p.EvalTree(a));}),
            new( "%",  2, 4, (p, a, b) => { return p.EvalTree(a) % p.EvalTree(b);} ),
            new( ">",  2, 7, (p, a, b) => { return p.EvalTree(a) > p.EvalTree(b) ? 1 : 0;} ),
            new( "<",  2, 7, (p, a, b) => { return p.EvalTree(a) < p.EvalTree(b) ? 1 : 0;} ),
            new( "&&", 2, 8, (p, a, b) => { return p.EvalTree(a) == 1 && p.EvalTree(b) == 1 ? 1 : 0;} ),
            new( "==", 2, 7, (p, a, b) => { return p.EvalTree(a) == p.EvalTree(b) ? 1 : 0;} ),
            new( "!=", 2, 7, (p, a, b) => { return p.EvalTree(a) != p.EvalTree(b) ? 1 : 0;} ),
            new( "||", 2, 9, (p, a, b) => { return p.EvalTree(a) == 1 || p.EvalTree(b) == 1 ? 1 : 0;} ),
            new( "!",  1, 1, (p, a, b) => { return !(p.EvalTree(a) == 1) ? 1 : 0;} ),
            new( ">=", 2, 7, (p, a, b) => { return p.EvalTree(a) >= p.EvalTree(b) ? 1 : 0;} ),
            new( "<=" ,    2, 7, (p, a, b) => { return p.EvalTree(a) <= p.EvalTree(b) ? 1 : 0;} )
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