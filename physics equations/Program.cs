using System;
using System.Reflection;


Console.WriteLine("Enter function name (a, b, or c):");
string functionName = Console.ReadLine();

// Get the method info using the function name
MethodInfo method = typeof(Program).GetMethod(functionName, BindingFlags.Static | BindingFlags.NonPublic);

// Check if the method exists
if (method != null)
{
    // Invoke the method
    method.Invoke(null, null);
}
else
{
    Console.WriteLine("Invalid function name");
}