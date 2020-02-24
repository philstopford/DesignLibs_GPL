using Eto.Forms;

namespace Error
{
    public static class ErrorReporter
    {
        public static void showMessage_OK(string stringToDisplay, string caption)
        {
            Application.Instance.Invoke(() =>
            {
                MessageBox.Show(stringToDisplay, caption, MessageBoxButtons.OK);
            });
        }
    }
}
