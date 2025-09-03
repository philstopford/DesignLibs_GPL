using EmailNS;
using NUnit.Framework;

namespace UnitTests;

/// <summary>
/// Tests for the Email library functionality.
/// Note: These tests focus on input validation and error handling rather than 
/// actual email sending to avoid dependency on external SMTP servers.
/// </summary>
[TestFixture]
public class EmailTests
{
    [Test]
    public void Send_WithEmptyHost_ShouldReturnSilently()
    {
        // Test that empty host parameter causes silent return
        Assert.DoesNotThrow(() => Email.Send("", "587", true, "Test", "Content", "test@example.com", "password"));
    }

    [Test]
    public void Send_WithEmptyPort_ShouldReturnSilently()
    {
        // Test that empty port parameter causes silent return
        Assert.DoesNotThrow(() => Email.Send("smtp.example.com", "", true, "Test", "Content", "test@example.com", "password"));
    }

    [Test]
    public void Send_WithEmptyAddress_ShouldReturnSilently()
    {
        // Test that empty address parameter causes silent return
        Assert.DoesNotThrow(() => Email.Send("smtp.example.com", "587", true, "Test", "Content", "", "password"));
    }

    [Test]
    public void Send_WithEmptyPassword_ShouldReturnSilently()
    {
        // Test that empty password parameter causes silent return
        Assert.DoesNotThrow(() => Email.Send("smtp.example.com", "587", true, "Test", "Content", "test@example.com", ""));
    }

    [Test]
    public void Send_WithInvalidHost_ShouldThrowException()
    {
        // Test that invalid host causes exception
        var exception = Assert.Throws<Exception>(() =>
            Email.Send("invalid.nonexistent.host", "587", true, "Test", "Content", "test@example.com", "password"));

        Assert.That(exception.Message, Is.EqualTo("Email problem"));
    }

    [Test]
    public void Send_WithInvalidPort_ShouldThrowException()
    {
        // Test that invalid port format causes exception
        var exception = Assert.Throws<Exception>(() =>
            Email.Send("smtp.example.com", "invalid_port", true, "Test", "Content", "test@example.com", "password"));

        Assert.That(exception.Message, Is.EqualTo("Email problem"));
    }

    [Test]
    public void Send_WithAllEmptyParameters_ShouldReturnSilently()
    {
        // Test that all empty parameters cause silent return
        Assert.DoesNotThrow(() => Email.Send("", "", true, "", "", "", ""));
    }

    [Test]
    public void Send_WithValidParametersButUnreachableServer_ShouldThrowException()
    {
        // Test with syntactically valid but unreachable server
        var exception = Assert.Throws<Exception>(() =>
            Email.Send("127.0.0.1", "9999", false, "Test Subject", "Test Content", "test@example.com", "testpass"));

        Assert.That(exception.Message, Is.EqualTo("Email problem"));
    }
}