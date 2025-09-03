using NUnit.Framework;
using info.lundin.math;
using System;
using System.Globalization;
using System.IO;

namespace UnitTests;

[TestFixture]
public class InfoLundinMathTests
{
    private ExpressionParser parser;

    [SetUp]
    public void Setup()
    {
        parser = new ExpressionParser();
    }

    [Test]
    public void Parser_Constructor_ShouldInitializeDefaultState()
    {
        var parser = new ExpressionParser();

        Assert.That(parser.RequireParentheses, Is.True);
        Assert.That(parser.ImplicitMultiplication, Is.True);
        Assert.That(parser.Culture, Is.EqualTo(CultureInfo.InvariantCulture));
        Assert.That(parser.Values, Is.Not.Null);
        Assert.That(parser.Expressions, Is.Not.Null);
    }

    [Test]
    public void Parse_SimpleArithmetic_ShouldEvaluateCorrectly()
    {
        Assert.That(parser.Parse("2 + 3"), Is.EqualTo(5.0).Within(1e-10));
        Assert.That(parser.Parse("10 - 4"), Is.EqualTo(6.0).Within(1e-10));
        Assert.That(parser.Parse("3 * 7"), Is.EqualTo(21.0).Within(1e-10));
        Assert.That(parser.Parse("15 / 3"), Is.EqualTo(5.0).Within(1e-10));
    }

    [Test]
    public void Parse_MathematicalFunctions_ShouldEvaluateCorrectly()
    {
        Assert.That(parser.Parse("sin(0)"), Is.EqualTo(0.0).Within(1e-10));
        Assert.That(parser.Parse("cos(0)"), Is.EqualTo(1.0).Within(1e-10));
        Assert.That(parser.Parse("sqrt(16)"), Is.EqualTo(4.0).Within(1e-10));
        Assert.That(parser.Parse("log(1)"), Is.EqualTo(0.0).Within(1e-10));
        Assert.That(parser.Parse("abs(-5)"), Is.EqualTo(5.0).Within(1e-10));
    }

    [Test]
    public void Parse_BuiltInConstants_ShouldEvaluateCorrectly()
    {
        Assert.That(parser.Parse("pi"), Is.EqualTo(Math.PI).Within(1e-10));
        Assert.That(parser.Parse("euler"), Is.EqualTo(Math.E).Within(1e-10));
        Assert.That(parser.Parse("true"), Is.EqualTo(1.0).Within(1e-10));
        Assert.That(parser.Parse("false"), Is.EqualTo(0.0).Within(1e-10));
        Assert.That(double.IsNaN(parser.Parse("nan")), Is.True);
        Assert.That(double.IsPositiveInfinity(parser.Parse("infinity")), Is.True);
    }

    [Test]
    public void Parse_Variables_ShouldUseAssignedValues()
    {
        parser.Values.Add("x", 5.0);
        parser.Values.Add("y", 3.0);

        Assert.That(parser.Parse("x + y"), Is.EqualTo(8.0).Within(1e-10));
        Assert.That(parser.Parse("x * y"), Is.EqualTo(15.0).Within(1e-10));
        Assert.That(parser.Parse("2 * x + y"), Is.EqualTo(13.0).Within(1e-10));
    }

    [Test]
    public void Parse_ComplexExpressions_ShouldEvaluateCorrectly()
    {
        parser.Values.Add("x", 2.0);

        Assert.That(parser.Parse("2 * x^2 + 3 * x + 1"), Is.EqualTo(15.0).Within(1e-10));
        Assert.That(parser.Parse("sin(pi/2)"), Is.EqualTo(1.0).Within(1e-10));
        Assert.That(parser.Parse("(2 + 3) * (4 - 1)"), Is.EqualTo(15.0).Within(1e-10));
    }

    [Test]
    public void Parse_ImplicitMultiplication_ShouldWork()
    {
        parser.Values.Add("x", 3.0);

        Assert.That(parser.Parse("2x"), Is.EqualTo(6.0).Within(1e-10));
        Assert.That(parser.Parse("3sin(0)"), Is.EqualTo(0.0).Within(1e-10));
        Assert.That(parser.Parse("2pi"), Is.EqualTo(2 * Math.PI).Within(1e-10));
    }

    [Test]
    public void Parse_WithImplicitMultiplicationDisabled_ShouldRequireExplicitOperator()
    {
        parser.ImplicitMultiplication = false;
        parser.Values.Add("x", 3.0);

        Assert.That(parser.Parse("2*x"), Is.EqualTo(6.0).Within(1e-10));
        Assert.Throws<ParserException>(() => parser.Parse("2x"));
    }

    [Test]
    public void Parse_WithParenthesesRequired_ShouldEnforceFunction()
    {
        parser.RequireParentheses = true;

        Assert.That(parser.Parse("sin(0)"), Is.EqualTo(0.0).Within(1e-10));
        Assert.Throws<ParserException>(() => parser.Parse("sin0"));
    }

    [Test]
    public void Parse_WithParenthesesNotRequired_ShouldAllowFunctionsWithoutParens()
    {
        parser.RequireParentheses = false;

        Assert.That(parser.Parse("sin0"), Is.EqualTo(0.0).Within(1e-10));
        Assert.That(parser.Parse("sin(0)"), Is.EqualTo(0.0).Within(1e-10));
    }

    [Test]
    public void Parse_NullOrEmptyExpression_ShouldThrowException()
    {
        Assert.Throws<ParserException>(() => parser.Parse(null));
        Assert.Throws<ParserException>(() => parser.Parse(""));
        Assert.Throws<ParserException>(() => parser.Parse("   "));
    }

    [Test]
    public void Parse_InvalidExpression_ShouldThrowException()
    {
        // Some expressions that should be invalid
        Assert.Throws<ParserException>(() => parser.Parse("2 +"));
        Assert.Throws<ParserException>(() => parser.Parse("* 3"));
        Assert.Throws<ParserException>(() => parser.Parse("((2 + 3)"));
        Assert.Throws<ParserException>(() => parser.Parse("2 3"));  // No operator between numbers
    }

    [Test]
    public void Culture_SetInvalidCulture_ShouldThrowException()
    {
        Assert.Throws<ArgumentNullException>(() => parser.Culture = null);

        // Test culture with division as decimal separator
        var persianCulture = new CultureInfo("fa-IR");
        if (persianCulture.NumberFormat.NumberDecimalSeparator == "/")
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => parser.Culture = persianCulture);
        }
    }

    [Test]
    public void Culture_SetValidCulture_ShouldWork()
    {
        var usCulture = new CultureInfo("en-US");
        parser.Culture = usCulture;

        Assert.That(parser.Culture, Is.EqualTo(usCulture));
    }

    [Test]
    public void Expression_Caching_ShouldReuseComputedExpressions()
    {
        // Parse the same expression multiple times - should cache and reuse
        var result1 = parser.Parse("2 + 3 * 4");
        var result2 = parser.Parse("2 + 3 * 4");

        Assert.That(result1, Is.EqualTo(result2));
        Assert.That(parser.Expressions.ContainsKey("2 + 3 * 4"), Is.True);
    }

    [Test]
    public void Expression_CachingWithVariables_ShouldReuseComputedExpressions()
    {
        parser.Values.Add("x", 5.0);
        var result = parser.Parse("2 + 3 * x");

        // Just verify the expression was processed correctly
        Assert.That(result, Is.EqualTo(17.0).Within(1e-10));

        // Check that the expression was cached
        Assert.That(parser.Expressions.ContainsKey("2 + 3 * x"), Is.True);

        // Parse the same expression again - should be cached
        var result2 = parser.Parse("2 + 3 * x");
        Assert.That(result2, Is.EqualTo(17.0).Within(1e-10));
    }

    [Test]
    public void EvalExpression_WithValidExpression_ShouldEvaluateCorrectly()
    {
        parser.Values.Add("x", 5.0);
        var result = parser.Parse("2 * x + 1");

        // Verify the expression was cached and can be retrieved
        Assert.That(parser.Expressions.ContainsKey("2 * x + 1"), Is.True);
        var expression = parser.Expressions["2 * x + 1"];

        var evalResult = parser.EvalExpression(expression);

        Assert.That(evalResult, Is.EqualTo(11.0).Within(1e-10));
    }

    [Test]
    public void Parse_DivisionByZero_ShouldReturnInfinity()
    {
        var result = parser.Parse("1 / 0");
        Assert.That(double.IsPositiveInfinity(result), Is.True);
    }

    [Test]
    public void Parse_NestedFunctions_ShouldEvaluateCorrectly()
    {
        Assert.That(parser.Parse("sin(cos(0))"), Is.EqualTo(Math.Sin(1.0)).Within(1e-10));
        Assert.That(parser.Parse("sqrt(abs(-16))"), Is.EqualTo(4.0).Within(1e-10));
    }

    [Test]
    public void Parse_PowerOperations_ShouldEvaluateCorrectly()
    {
        Assert.That(parser.Parse("2^3"), Is.EqualTo(8.0).Within(1e-10));
        Assert.That(parser.Parse("2^(1/2)"), Is.EqualTo(Math.Sqrt(2)).Within(1e-10));
        Assert.That(parser.Parse("(-2)^2"), Is.EqualTo(4.0).Within(1e-10));
    }

    [Test]
    public void Values_Dictionary_ShouldManageVariablesCorrectly()
    {
        parser.Values.Add("a", 10.0);
        parser.Values.Add("b", 20.0);

        Assert.That(parser.Values.ContainsKey("a"), Is.True);
        Assert.That(parser.Values.ContainsKey("b"), Is.True);
        Assert.That(parser.Values.ContainsKey("c"), Is.False);

        parser.Values.Remove("a");
        Assert.That(parser.Values.ContainsKey("a"), Is.False);
    }
}