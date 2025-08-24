using MimeKit;
using SmtpClient = MailKit.Net.Smtp.SmtpClient;

namespace EmailNS;

/// <summary>
/// Provides simple email sending functionality using SMTP.
/// </summary>
public static class Email
{
    /// <summary>
    /// Sends an email message via SMTP.
    /// </summary>
    /// <param name="host">SMTP server hostname or IP address</param>
    /// <param name="port">SMTP server port number as string</param>
    /// <param name="ssl">Whether to use SSL/TLS encryption</param>
    /// <param name="subject">Email subject line</param>
    /// <param name="messageContent">Plain text message body</param>
    /// <param name="address">Email address for both sender and recipient authentication</param>
    /// <param name="password">Password for SMTP authentication</param>
    /// <exception cref="Exception">Thrown when email sending fails with message "Email problem"</exception>
    public static void Send(string host, string port, bool ssl, string subject, string messageContent, string address, string password)
    {
        if (host == "" || port == "" || address == "" || password == "")
        {
            return;
        }

        try
        {
            SmtpClient client = new() {ServerCertificateValidationCallback = (s, c, h, e) => true};
            client.Connect(host, Convert.ToInt32(port), ssl);
            MimeMessage message = new()
            {
                Subject = subject, Body = new TextPart("plain") {Text = messageContent}
            };
            message.From.Add(new MailboxAddress(address, address));
            message.To.Add(new MailboxAddress(address, address));
            client.Authenticate(address, password);
            client.Send(message);
            client.Disconnect(true);
        }
        catch
        {
            throw new Exception( "Email problem");
        }
    }
}