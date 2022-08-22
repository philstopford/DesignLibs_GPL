using MimeKit;
using SmtpClient = MailKit.Net.Smtp.SmtpClient;

namespace EmailNS;

public static class Email
{
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