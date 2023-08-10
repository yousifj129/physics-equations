using physics_equations;
using System;
using System.Numerics;
using System.Reflection;
using static helperFunctions.HelperFunctions;


while (true)
{
    Console.WriteLine("Enter a command:");
    string input = Console.ReadLine();
    if(input == "clear")
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
        ParameterInfo[] parameters = method.GetParameters();
        object[] parsedValues = new object[parameters.Length];

        // Loop through the parameters and parse the corresponding input values
        for (int i = 0; i < parameters.Length; i++)
        {
            ParameterInfo parameter = parameters[i];
            Type parameterType = parameter.ParameterType;

            // Check if the input value can be parsed to the parameter type
            if (i + 1 < inputArray.Length)
            {
                string inputValue = inputArray[i + 1];
                object parsedValue = null;

                physicsEquations.random = new Random().NextDouble();
                
                // Check if the input is a valid constant/variable
                FieldInfo constantOrVariable = typeof(PhysicsEquations).GetField(inputValue, BindingFlags.Instance | BindingFlags.Public | BindingFlags.NonPublic | BindingFlags.Static);
                if (constantOrVariable != null)
                {
                    // Retrieve the constant/variable value
                    parsedValue = constantOrVariable.GetValue(physicsEquations);
                }
                else if (TryParseValue(inputValue, parameterType, out parsedValue))
                {
                    // Input value can be parsed to the parameter type
                }
                else
                {
                    // Invalid parameter entered
                    Console.WriteLine("Invalid parameter(s) entered");
                    break;
                }

                parsedValues[i] = parsedValue;
            }
            else
            {
                // Insufficient number of parameters
                Console.WriteLine("Insufficient number of parameters");
                break;
            }
        }
        // Check if all parameters were successfully parsed
        if (parsedValues.All(value => value != null))
        {
            // Invoke the method on the instance with the modified parameters
            object result = method.Invoke(physicsEquations, parsedValues);

            // Check if the method has a return value
            if (result != null)
            {
                Type returnType = method.ReturnType;

                if (returnType == typeof(double))
                {
                    // Cast the returned value to double
                    double returnValue = (double)result;
                    Console.WriteLine("Result: " + returnValue);
                }
                else if (returnType == typeof(Vector3))
                {
                    // Cast the returned value to Vector3
                    Vector3 returnVector = (Vector3)result;
                    Console.WriteLine("Result: " + returnVector);
                }
                else
                {
                    // Handle unsupported return type
                    Console.WriteLine("Unsupported return type");
                }
            }
            commandFound = true;
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
        try
        {
            if(method != null)
            {
                ParameterInfo[] p = method.GetParameters();
                string parametersNeeded = "\n";
                foreach (ParameterInfo pi in p)
                {
                    parametersNeeded += pi.Name + ": " + pi.ParameterType.ToString() + ", ";
                }
                Console.WriteLine("Invalid command, " + $"variables needed for the function {method.Name}:" + parametersNeeded + $"\n{p.Length} parameters is needed");
            }
           
        }
        catch
        {
            Console.WriteLine("Invalid command");

        }

    }
}
