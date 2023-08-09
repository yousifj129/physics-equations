using physics_equations;
using System;
using System.Reflection;


while (true)
{

    Console.WriteLine("Enter a command:");
    string input = Console.ReadLine();
    if(input == "c")
    {
        Console.Clear();
        continue;
    }
    // Split the input into an array of strings
    string[] inputArray = input.Split(' ');

    // Extract the function name from the first element
    string functionName = inputArray[0];

    // Create an instance of the PhysicsEquations class
    var physicsEquations = new PhysicsEquations();

    // Get the method info using the function name
    MethodInfo method = typeof(PhysicsEquations).GetMethod(functionName, BindingFlags.Instance | BindingFlags.Public | BindingFlags.NonPublic);
    bool commandFound = false;
    // Check if the method exists and has the correct number of parameters

    if (method != null && method.GetParameters().Length <= inputArray.Length - 1)
    {
        // Convert the input values to doubles using TryParse
        double[] parsedValues = inputArray.Skip(1).Select(s =>
        {
            double value;
            if (double.TryParse(s, out value))
            {
                return value;
            }

            // Check if the input is a valid constant/variable
            FieldInfo constantOrVariable = typeof(PhysicsEquations).GetField(s, BindingFlags.Public | BindingFlags.Static);
            if (constantOrVariable != null)
            {
                // Retrieve the constant/variable value
                double constantOrVariableValue = (double)constantOrVariable.GetValue(null);
                return constantOrVariableValue;
            }

            // Invalid parameter entered
            return double.NaN;
        }).ToArray();

        // Check if any parsing errors occurred
        if (parsedValues.Contains(double.NaN))
        {
            Console.WriteLine("Invalid parameter(s) entered");
        }
        else
        {
            // Invoke the method on the instance with the modified parameters
            object result = method.Invoke(physicsEquations, parsedValues.Cast<object>().ToArray());

            // Check if the method has a return value
            if (result != null)
            {
                // Cast the returned value to the appropriate type
                double returnValue = (double)result;

                // Use the returned value as needed
                Console.WriteLine("Result: " + returnValue);
                commandFound = true;
            }
        }
    }
    else
    {
        // Check for a variable or constant with the given name
        FieldInfo variable = typeof(PhysicsEquations).GetField(functionName, BindingFlags.Instance | BindingFlags.Public | BindingFlags.NonPublic | BindingFlags.Static);
        if (variable != null)
        {
            // Print the value of the variable or constant
            object value = variable.GetValue(physicsEquations);
            Console.WriteLine("Variable or constant value: " + value);
            commandFound = true;
        }
    }
    if (!commandFound)
    {
        Console.WriteLine("Invalid command");
    }
}
