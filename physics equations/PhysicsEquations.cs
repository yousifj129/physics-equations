using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;

using static helperFunctions.HelperFunctions;
using System.Reflection;

namespace physics_equations
{
    public class PhysicsEquations
    {
        #region constants
        // Gravitational constant (units: m^3 kg^-1 s^-2)
        public const double GravitationalConstant = 6.67430e-11;

        public const double AccelerationConstant = 1.2e-10; // MOND acceleration constant in m/s^2

        // Planck constant (units: J s)
        public const double PlanckConstant = 6.62607015e-34;

        // Speed of light in a vacuum (units: m/s)
        public const double SpeedOfLight = 299792458;

        // Elementary charge (units: C)
        public const double ElementaryCharge = 1.602176634e-19;

        // Avogadro's number (units: mol^-1)
        public const double AvogadroNumber = 6.02214076e23;

        // Boltzmann constant (units: J K^-1)
        public const double BoltzmannConstant = 1.380649e-23;

        // Electron mass (units: kg)
        public const double ElectronMass = 9.10938356e-31;

        // Proton mass (units: kg)
        public const double ProtonMass = 1.67262192369e-27;

        // Speed of sound in air at 20°C (units: m/s)
        public const double SpeedOfSoundInAir = 343;

        // Universal gas constant (units: J K^-1 mol^-1)
        public const double UniversalGasConstant = 8.314462618;

        // Stefan-Boltzmann constant (units: W m^-2 K^-4)
        public const double StefanBoltzmannConstant = 5.670374419e-8;
        // Electron charge (units: C)
        public const double ElectronCharge = -1.602176634e-19;

        // Proton charge (units: C)
        public const double ProtonCharge = 1.602176634e-19;

        // Neutron charge (units: C)
        public const double NeutronCharge = 0;

        // Electron volt (units: J)
        public const double ElectronVolt = 1.602176634e-19;

        // Vacuum permeability (units: N A^-2)
        public const double VacuumPermeability = 4 * Math.PI * 1e-7;

        // Vacuum permittivity (units: F m^-1)
        public const double VacuumPermittivity = 8.8541878128e-12;

        // Planck constant over 2π (units: J s)
        public const double ReducedPlanckConstant = PlanckConstant / (2 * Math.PI);

        // Fine-structure constant
        public const double FineStructureConstant = 7.2973525693e-3;

        // Rydberg constant (units: m^-1)
        public const double RydbergConstant = 10973731.568508;

        // Bohr radius (units: m)
        public const double BohrRadius = 5.291772109038e-11;

        // Electron rest mass (units: kg)
        public const double ElectronRestMass = 9.10938356e-31;

        // Proton rest mass (units: kg)
        public const double ProtonRestMass = 1.67262192369e-27;

        // Neutron rest mass (units: kg)
        public const double NeutronRestMass = 1.67492749804e-27;

        // Speed of light in a vacuum (units: km/s)
        public const double SpeedOfLightInKmPerSec = SpeedOfLight / 1000;

        // Planck temperature (units: K)
        public const double PlanckTemperature = 1.416808e32;

        // Planck mass (units: kg)
        public const double PlanckMass = 2.176434e-8;

        // Planck length (units: m)
        public const double PlanckLength = 1.616255e-35;

        // Planck time (units: s)
        public const double PlanckTime = 5.39116e-44;
        
        public const double CoulombConstant = 8.9875517923e9; // Coulomb's constant in N·m^2/C^2
        
        public double random = new Random().NextDouble();

        #endregion

        #region basic math functions
        public double add(double a, double b)
        {
            return a + b;
        }
        public double divide(double a, double b)
        {
            return a / b;
        }
        public double div(double a, double b)
        {
            return a / b;
        }
        public double multiply(double a, double b)
        {
            return a * b;
        }
        public double mult(double a, double b)
        {
            return a * b;
        }
        public double subtract(double a, double b)
        {
            return a - b;
        }
        public double sub(double a, double b)
        {
            return a - b;
        }
        public double Random(double minValue, double maxValue)
        {
            Random random = new Random();
            return random.NextDouble() * (maxValue - minValue) + minValue;
        }

        #endregion

        #region uncertainty princples 
        public double CalculateUncertainty(double deltaP)
        {
            double uncertainty = ReducedPlanckConstant / (2 * deltaP);
            return uncertainty;
        }
        public double CalculateTimeUncertainty(double deltaE)
        {
            double uncertainty = ReducedPlanckConstant / (2 * deltaE);
            return uncertainty;
        }

        #endregion

        #region probability
        public double CalculateProbability(double estimateNumber, double NumberOfTries)
        {
            double probability = estimateNumber / NumberOfTries;
            return probability;
        }
        public double ComplementaryProbability(double probability)
        {
            double complementaryProbability = 1 - probability;
            return complementaryProbability;
        }
        public double JointProbability(double probabilityA, double probabilityB)
        {
            double jointProbability = probabilityA * probabilityB;
            return jointProbability;
        }
        public double ProbabilityDensity(double wavefunction)
        {
            double probabilityDensity = Math.Pow(wavefunction, 2);
            return probabilityDensity;
        }
        public double TransitionProbability(double initialWavefunction, double finalWavefunction)
        {
            double transitionProbability = Math.Pow(Math.Abs(initialWavefunction * finalWavefunction), 2);
            return transitionProbability;
        }
        public double FermiDiracDistribution(double energy, double temperature, double mu)
        {
            double k = 8.617333262145e-5; // Boltzmann constant
            double fermiDiracDistribution = 1 / (Math.Exp((energy - mu) / (k * temperature)) + 1);
            return fermiDiracDistribution;
        }
        #endregion

        #region wave functions

        public double BoxWaveFunction(double amplitude, int n, double boxLength, double position)
        {
            double wavefunction = amplitude * Math.Sin(n * Math.PI * position / boxLength);
            return wavefunction;
        }
        public  double GaussianWaveFunction(double amplitude, double sigma, double x, double x0)
        {
            double wavefunction = amplitude * Math.Exp(-(Math.Pow(x - x0, 2)) / (2 * Math.Pow(sigma, 2)));
            return wavefunction;
        }
        public double HarmonicOscillatorWaveFunction(double amplitude, int n, double x, double x0, double k)
        {
            double wavefunction = amplitude * Math.Pow(k / (Math.PI * n), 0.25) * Math.Exp(-(Math.Pow(x - x0, 2) * k) / (2 * n));
            return wavefunction;
        }
        public double CoulombWaveFunction(double amplitude, int l, double k, double eta, double radius)
        {
            double wavefunction = amplitude * Math.Pow(k * radius, l) * Math.Exp(-Math.PI * eta / 2) * Math.Pow(2 * eta, l + 1) *
                Math.Sqrt(Math.PI * Math.Exp(2 * Math.PI * eta)) * BesselFunction(l + 1.5, k * radius * eta) / (k * radius);

            return wavefunction;
        }

        #endregion

        #region helper functions
        public double BesselFunction(double nu, double z)
        {
            int terms = 100; // Number of terms in the series expansion

            double result = 0.0;
            double numerator = 1.0;
            double denominator = 1.0;
            double factorial = 1.0;

            for (int n = 0; n < terms; n++)
            {
                double term = numerator / (denominator * factorial);
                result += term;

                numerator *= -0.25 * z * z;
                denominator *= n + 1;
                factorial *= n + 1;
            }

            return result;
        }
        #endregion

        #region physics
        public double KineticEnergy(double mass, double velocity)
        {
            double kineticEnergy = 0.5 * mass * Math.Pow(velocity, 2);
            return kineticEnergy;
        }
        public double GravitationalPotentialEnergy(double mass, double height, double gravitationalAcceleration)
        {
            double gravitationalPotentialEnergy = mass * gravitationalAcceleration * height;
            return gravitationalPotentialEnergy;
        }
        public double GravitationalEnergyMOND(double mass1, double mass2, double distance)
        {
            double gravitationalEnergy = (GravitationalConstant * mass1 * mass2) / (2 * AccelerationConstant) * Math.Log((distance * Math.Sqrt(AccelerationConstant)) / GravitationalConstant + 1);
            return gravitationalEnergy;
        }
        public double GravitationalEnergyNewton(double mass1, double mass2, double distance)
        {
            double gravitationalEnergy = (GravitationalConstant * mass1 * mass2) / distance;
            return gravitationalEnergy;
        }
        public  double GeneralRelativityGravity( double mass1, double mass2, double distance)
        {
            // Schwarzschild radius of the system
            double rs = 2 * GravitationalConstant * (mass1 + mass2) / Math.Pow(SpeedOfLight, 2);

            // Check if the objects are inside their combined Schwarzschild radius
            if (distance <= rs)
            {
                return 0;
            }

            // Calculate the gravitational force using the general relativity equation
            double force = -GravitationalConstant * mass1 * mass2 / (distance * distance) *
                           (1 - rs / distance) / Math.Sqrt(1 - rs / distance);

            return force;
        }
        #endregion

        #region quantum physics
        public double ElectrostaticEnergyC(double charge1, double charge2, double distance)
        {
            double electrostaticEnergy = CoulombConstant * charge1 * charge2 / distance;
            return electrostaticEnergy;
        }
        #endregion
        
        #region vectors
        public Vector3 AddVector(Vector3 vector1, Vector3 vector2)
        {
            Vector3 result = vector1 + vector2;
            return result;
        }
        public Vector3 SubtractVector(Vector3 vector1, Vector3 vector2)
        {
            Vector3 result = vector1 - vector2;
            return result;
        }
        public Vector3 DivideVector(Vector3 vector1, Vector3 vector2)
        {
            Vector3 result = vector1 / vector2;
            return result;
        }
        public Vector3 MultiplyVector(Vector3 vector1, Vector3 vector2)
        {
            Vector3 result = vector1 * vector2;
            return result;
        }
        public double distanceSquaredVector(Vector3 pos1, Vector3 pos2)
        {
            double dx = pos2.X - pos1.X;
            double dy = pos2.Y - pos1.Y;
            double dz = pos2.Z - pos1.Z;

            return dx * dx + dy * dy + dz * dz;
        }
        #endregion

        #region user input
        public static string userInput(string input)
        {
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
                        FieldInfo constantOrVariable = typeof(PhysicsEquations).GetField(inputValue, BindingFlags.Instance | BindingFlags.Public | BindingFlags.NonPublic | BindingFlags.Static);
                        // Check if the input is a valid constant/variable

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
                            return returnValue.ToString();
                        }
                        else if (returnType == typeof(Vector3))
                        {
                            // Cast the returned value to Vector3
                            Vector3 returnVector = (Vector3)result;
                            return returnVector.ToString();
                        }
                        else
                        {
                            // Handle unsupported return type
                            return ("Unsupported return type");
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
                    return "Variable or constant value: " + value;
                    commandFound = true;
                }
            }
            if (!commandFound)
            {
                try
                {
                    if (method != null)
                    {
                        ParameterInfo[] p = method.GetParameters();
                        string parametersNeeded = "\n";
                        foreach (ParameterInfo pi in p)
                        {
                            parametersNeeded += pi.Name + ": " + pi.ParameterType.ToString() + ", ";
                        }
                        return "Invalid command, " + $"variables needed for the function {method.Name}:" + parametersNeeded + $"\n{p.Length} parameters is needed";
                    }

                }
                catch
                {
                    return "Invalid command";

                }

            }
            return "";
        }
        #endregion
    }
}
