using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Reflection;

namespace MiscUtil;

/// <summary>
/// Class allowing a user to choose a console application to run, after
/// examining an assembly for all classes containing a static Main method, either
/// parameterless or with a string array parameter.
/// </summary>
public class ApplicationChooser
{
    private const string Keys = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";

    /// <summary>
    /// Displays entry points and prompts the user to choose one.
    /// </summary>
    /// <param name="type">Type within the assembly containing the applications. This type is
    /// not included in the list of entry points to run.</param>
    /// <param name="args">Arguments to pass in for methods which have a single string[] parameter.</param>
    public static void Run(Type type, string[] args)
    {
        Assembly assembly = type.Assembly;

        List<MethodBase> entryPoints = [];
        foreach (Type candidate in assembly.GetTypes())
        {
            if (candidate == type)
            {
                continue;
            }
            MethodBase entryPoint = GetEntryPoint(candidate);
            if (entryPoint != null)
            {
                entryPoints.Add(entryPoint);
            }
        }

        entryPoints.Sort((x, y) =>
        {
            if (x.DeclaringType == null)
            {
                return 0;
            }

            return y.DeclaringType != null ? string.Compare(x.DeclaringType.Name, y.DeclaringType.Name, StringComparison.Ordinal) : 0;
        });

        switch (entryPoints.Count)
        {
            case 0:
                Console.WriteLine("No entry points found. Press return to exit.");
                Console.ReadLine();
                return;
        }

        for (int i = 0; i < entryPoints.Count; i++)
        {
            Console.WriteLine("{0}: {1}", Keys[i],
                GetEntryPointName(entryPoints[i]));
        }
        Console.WriteLine();
        Console.Write("Entry point to run? ");
        Console.Out.Flush();
        char key = Console.ReadKey().KeyChar;
        Console.WriteLine();

        switch (key)
        {
            // "Enter" means "Oops, let's just quit"
            case '\r':
                return;
        }

        int entry = Keys.IndexOf(char.ToUpper(key));
        if (entry == -1 || entry >= entryPoints.Count)
        {
            Console.WriteLine("Invalid choice");
        }
        else
        {
            try
            {
                MethodBase main = entryPoints[entry];
                main.Invoke(null, main.GetParameters().Length == 0 ? null : [args]);
            }
            catch (Exception e)
            {
                Console.WriteLine("Exception: {0}", e);
            }
        }
        Console.WriteLine();
        Console.WriteLine("Press return to exit.");
        Console.ReadLine();
    }

    private static object GetEntryPointName(MethodBase methodBase)
    {
        Type type = methodBase.DeclaringType;

        object[] descriptions = type.GetCustomAttributes(typeof(DescriptionAttribute), false);
        return descriptions.Length == 0 ?
            type.Name :
            string.Format("{0} [{1}]", type.Name, ((DescriptionAttribute)descriptions[0]).Description);
    }

    /// <summary>
    /// Returns the entry point for a method, or null if no entry points can be used.
    /// An entry point taking string[] is preferred to one with no parameters.
    /// </summary>
    internal static MethodBase GetEntryPoint(Type type)
    {
        if (type.IsGenericTypeDefinition || type.IsGenericType)
        {
            return null;
        }

        BindingFlags anyStatic = BindingFlags.DeclaredOnly | BindingFlags.Static | BindingFlags.Public | BindingFlags.NonPublic;
        // Can't use GetMethod directly as then we can't ignore generic methods :(

        MethodInfo[] methods = type.GetMethods(anyStatic);

        MethodInfo parameterless = null;
        MethodInfo stringArrayParameter = null;

        foreach (MethodInfo method in methods)
        {
            if (method.Name != "Main")
            {
                continue;
            }
            if (method.IsGenericMethod || method.IsGenericMethodDefinition)
            {
                continue;
            }
            ParameterInfo[] parameters = method.GetParameters();
            switch (parameters.Length)
            {
                case 0:
                    parameterless = method;
                    break;
                case 1 when !parameters[0].IsOut && !parameters[0].IsOptional && parameters[0].ParameterType == typeof(string[]):
                    stringArrayParameter = method;
                    break;
            }
        }

        // Prefer the version with parameters, return null if neither have been found
        return stringArrayParameter ?? parameterless;
    }
}