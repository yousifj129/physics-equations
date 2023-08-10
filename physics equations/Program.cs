using physics_equations;
using System;
using System.Numerics;
using System.Reflection;
using static helperFunctions.HelperFunctions;

PhysicsEquations p = new PhysicsEquations();

while (true)
{
    Console.WriteLine("Enter a command:");
    string input = Console.ReadLine();
    if (input == "clear")
    {
        Console.Clear();
        break;
    }
    Console.WriteLine("results : " + p.userInput(input));
}

/*double a = 0; 
double b = 0;
double.TryParse( PhysicsEquations.userInput("GravitationalPotentialEnergy 55 10 9.81"), out a);
double.TryParse(PhysicsEquations.userInput("KineticEnergy 55 20"), out b);
double distance = 0;
double.TryParse(PhysicsEquations.userInput("distanceSquaredVector 3,10,3 3,0,3"), out distance);

Console.WriteLine(a + b); Console.WriteLine(distance*0.1);
*/