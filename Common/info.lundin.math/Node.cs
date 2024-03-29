﻿/*
 * Author: Patrik Lundin, patrik@lundin.info
 * Web: http://www.lundin.info
 * 
 * Source code released under the Microsoft Public License (Ms-PL) 
 * http://www.microsoft.com/en-us/openness/licenses.aspx#MPL
*/

using System;

namespace info.lundin.math;

/// <summary>
/// Types of Node(s)
/// </summary>
public enum NodeType
{
    Undefined,
    Expression,
    Variable,
    Value
}

/// <summary>
/// Class Node, represents a Node in a tree data structure representation of a mathematical expression.
/// </summary>
[Serializable]
public class Node
{
    private Operator op;

    private Node arg1;
    private Node arg2;

    private int args;

    private NodeType type;

    /// <summary>
    /// Backing variable for property. Should absolutely not be serialized
    /// since that will cause stale values to be persisted.
    /// </summary>
    [NonSerialized]
    private Value value;

    private string variable;

    /// <summary>
    /// Creates a Node containing the specified Operator and arguments.
    /// This will automatically mark this Node as a TYPE_EXPRESSION
    /// </summary>
    /// <param name="_operator">the string representing an operator</param>
    /// <param name="_arg1">the first argument to the specified operator</param>
    /// <param name="_arg2">the second argument to the specified operator</param>
    internal Node(Operator op, Node arg1, Node arg2)
    {
        this.arg1 = arg1;
        this.arg2 = arg2;
        this.op = op;
        args = 2;

        type = NodeType.Expression;
    }

    /// <summary>
    /// Creates a Node containing the specified Operator and argument.
    /// This will automatically mark this Node as a TYPE_EXPRESSION
    /// </summary>
    /// <param name="_operator">the string representing an operator</param>
    /// <param name="_arg1">the argument to the specified operator</param>
    internal Node(Operator op, Node arg1)
    {
        this.arg1 = arg1;
        this.op = op;
        args = 1;

        type = NodeType.Expression;
    }

    /// <summary>
    /// Creates a Node containing the specified variable.
    /// This will automatically mark this Node as a TYPE_VARIABLE
    /// </summary>
    /// <param name="variable">the string representing a variable</param>
    internal Node(string variable)
    {
        this.variable = variable;
        type = NodeType.Variable;
    }

    /// <summary>
    /// Creates a Node containing the specified value.
    /// This will automatically mark this Node as a TYPE_CONSTANT
    /// </summary>
    /// <param name="value">the value for this Node</param>
    internal Node(double value)
    {
        this.value = new DoubleValue { Value = value };
        type = NodeType.Value;
    }

    /// <summary>
    /// Returns the String operator of this Node 
    /// </summary>
    internal Operator Operator => op;

    /// <summary>
    /// Gets or sets the value of this Node 
    /// </summary>
    internal Value Value
    {
        get => value;
        set => this.value = value;
    }

    /// <summary>
    /// Returns the String variable of this Node 
    /// </summary>
    internal string Variable => variable;

    /// <summary>
    /// Returns the number of arguments this Node has
    /// </summary>
    internal int Arguments => args;

    /// <summary>
    /// Returns the node type
    /// </summary>
    internal NodeType Type => type;

    /// <summary>
    /// Returns the first argument of this Node
    /// </summary>
    internal Node FirstArgument => arg1;

    /// <summary>
    /// Returns the second argument of this Node
    /// </summary>
    internal Node SecondArgument => arg2;
} // End class Node