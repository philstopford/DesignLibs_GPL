using System;
using System.Windows;

namespace DynamicAero2.Styles;

partial class ContextMenu
{
    public ContextMenu()
    {
        InitializeComponent();

        AddPrivateControlStyle(typeof(Application).Assembly.GetType("System.Windows.Documents.TextEditorContextMenu+EditorContextMenu"), typeof(System.Windows.Controls.ContextMenu));
        AddPrivateControlStyle(typeof(Application).Assembly.GetType("System.Windows.Documents.TextEditorContextMenu+EditorMenuItem"), typeof(System.Windows.Controls.MenuItem));
    }

    private void AddPrivateControlStyle(Type targetType, Type baseType)
    {
        Style baseStyle = this[baseType] as Style;
        Style style = new(targetType, baseStyle);
        Add(targetType, style);
    }
}