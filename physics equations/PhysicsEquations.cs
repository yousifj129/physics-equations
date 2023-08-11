using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;

using static helperFunctions.HelperFunctions;
using System.Reflection;
using System.Numerics.Tensors;
using System.Drawing;

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

        public const double PI = Math.PI;

        public const double goldenRatio = 1.618033988749894848204586834365638117720309179805762862135448622705260462818902449707207204189391;

        public const double EulerNumber = 2.718281828459045235360287471352662497757247093699959574966967627724076630353547594571382178525166427;


        public double random = new Random().NextDouble();

        #endregion

        #region basic math functions
        public double add(double number1, double number2)
        {
            return number1 + number2;
        }
        public double divide(double number1, double number2)
        {
            return number1 / number2;
        }
        public double div(double number1, double number2)
        {
            return number1 / number2;
        }
        public double multiply(double number1, double number2)
        {
            return number1 * number2;
        }
        public double mult(double number1, double number2)
        {
            return number1 * number2;
        }
        public double subtract(double number1, double number2)
        {
            return number1 - number2;
        }
        public double sub(double number1, double number2)
        {
            return number1 - number2;
        }
        public double Random(double minValue, double maxValue)
        {
            Random random = new Random();
            return random.NextDouble() * (maxValue - minValue) + minValue;
        }

        #endregion

        #region Calculus Helpful Functions
        public double[] CalculatePartialDerivatives(double[] variables, double[] metricComponents)
        {
            PartialDerivativesCalculator pdc = new PartialDerivativesCalculator();
            return pdc.CalculatePartialDerivatives(variables, metricComponents);
        }
        public double Sigmoid(double x)
        {
            return 1.0 / (1.0 + Math.Exp(-x));
        }
        public double Sin(double x)
        {
            return Math.Sin(x);
        }
        public double Cos(double x)
        {
            return Math.Cos(x);
        }
        public double Tan(double x)
        {
            return Math.Tan(x);
        }
        public double Sinh(double x)
        {
            return Math.Sinh(x);
        }
        public double Cosh(double x)
        {
            return Math.Cosh(x);
        }
        public double Tanh(double x)
        {
            return Math.Tanh(x);
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

            #region kinetic stuff
        public double KineticEnergy(double mass, double velocity)
        {
            double kineticEnergy = 0.5 * mass * Math.Pow(velocity, 2);
            return kineticEnergy;
        }
        public double KineticVelocity(double mass, double kineticEnergy)
        {
            // Ensure the mass and kinetic energy are positive
            double positiveMass = Math.Abs(mass);
            double positiveKineticEnergy = Math.Abs(kineticEnergy);

            // Calculate the velocity using the formula: velocity = sqrt(2 * kineticEnergy / mass)
            double velocity = Math.Sqrt((2 * positiveKineticEnergy) / positiveMass);

            return velocity;
        }
        public Vector2 KineticVelocityVector2(double mass, double kineticEnergy, double angle)
        {
            // Ensure the mass and kinetic energy are positive
            double positiveMass = Math.Abs(mass);
            double positiveKineticEnergy = Math.Abs(kineticEnergy);

            // Calculate the velocity magnitude using the formula: velocityMagnitude = sqrt(2 * kineticEnergy / mass)
            double velocityMagnitude = Math.Sqrt((2 * positiveKineticEnergy) / positiveMass);

            // Convert the angle from degrees to radians
            double angleInRadians = Math.PI * angle / 180.0;

            // Calculate the horizontal and vertical components of the velocity
            double velocityX = velocityMagnitude * Math.Cos(angleInRadians);
            double velocityY = velocityMagnitude * Math.Sin(angleInRadians);

            // Create and return the velocity as a Vector2
            Vector2 velocity = new Vector2((float)velocityX, (float)velocityY);
            return velocity;
        }
        public Vector3 KineticVelocityVector3(double mass, double kineticEnergy, Vector3 direction)
        {
            // Ensure the mass and kinetic energy are positive
            double positiveMass = Math.Abs(mass);
            double positiveKineticEnergy = Math.Abs(kineticEnergy);

            // Calculate the velocity magnitude using the formula: velocityMagnitude = sqrt(2 * kineticEnergy / mass)
            double velocityMagnitude = Math.Sqrt((2 * positiveKineticEnergy) / positiveMass);

            // Normalize the direction vector
            Vector3 normalizedDirection = Vector3.Normalize(direction);

            // Calculate the velocity components by multiplying the magnitude with the normalized direction
            Vector3 velocity = MultiplyVectorByDouble(normalizedDirection, velocityMagnitude);

            return velocity;
        }
        #endregion

            #region gravity stuff
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
        public double GravitationalForceBergman(double mass1, double mass2, double distance)
        {
            double gravitationalConstant = 6.67430e-11; // Gravitational constant
            double alpha = 0.5; // Constant specific to Bergman theory

            double force = (gravitationalConstant * mass1 * mass2) / (Math.Pow(distance, 2) * (1 + alpha * distance));

            return force;
        }
        //this is wrong, i am just gonna make it later
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
        public Vector3 ModifiedNewtonianDynamicsVelocity(double mass1, Vector3 pos1, double mass2, Vector3 pos2)
        {
            double a0 = 1.2e-10; // MOND acceleration constant
            Vector3 r = pos2 - pos1;
            double R = r.Length();
            double m = mass1 + mass2;
            double u = R / a0;

            double f;
            if (u > 1e-3)
            {
                f = GravitationalConstant * m / (R * R);
            }
            else
            {
                f = GravitationalConstant * m * a0 / (1 + u) * Math.Exp(-u);
            }
             
            return Vector3.Normalize(r) * (float)f * (float)mass2;
        }
        public Vector3 NewtonianGravityVelocity(double mass1, Vector3 pos1, double mass2, Vector3 pos2)
        {
            Vector3 r = pos2 - pos1;
            double R = r.Length();
            double forceMagnitude = GravitationalConstant * mass1 * mass2 / (R * R * R);

            return Vector3.Normalize(r) * (float)forceMagnitude / (float)mass2;
        }
        public Vector3 GeneralRelativityEinsteinSchwarzschild(double mass1, Vector3 pos1, double mass2, Vector3 pos2)
        {
            double m = mass1 + mass2;
            Vector3 r = pos2 - pos1;
            double R = r.Length();

            // Schwarzschild radius of the system
            double rs = 2 * GravitationalConstant * m / (SpeedOfLight * SpeedOfLight);

            // Check if the objects are inside their combined Schwarzschild radius
            if (R <= rs)
            {
                return Vector3.Zero;
            }

            // Schwarzschild metric tensor
            Matrix4x4 g = new Matrix4x4();
            g.M11 = (float)-(1 - rs / R);
            g.M22 = (float)(1 / (1 - rs / R));
            g.M33 = (float)(r.Y * r.Y / R / R);
            g.M44 = (float)(r.Z * r.Z / R / R);
            g.M12 = g.M21 = (float)(-rs * r.X / R / R / R / 2);
            g.M13 = g.M31 = (float)(-rs * r.Y * r.X / R / R / R / 2);
            g.M14 = g.M41 = (float)(-rs * r.Z * r.X / R / R / R / 2);
            g.M23 = g.M32 = (float)(-rs * r.Y * r.Y * r.X / R / R / R / 2);
            g.M24 = g.M42 = (float)(-rs * r.Y * r.Z * r.X / R / R / R / 2);
            g.M34 = g.M43 = (float)(-rs * r.Z * r.Y * r.X / R / R / R / 2);

            // Geodesic equation
            Vector3 d2xdt2 = new Vector3(
                (float)(0.5 * (g.M21 * g.M11 + g.M22 * g.M12 + g.M23 * g.M13 + g.M24 * g.M14)),
                (float)(0.5 * (g.M31 * g.M11 + g.M32 * g.M12 + g.M33 * g.M13 + g.M34 * g.M14)),
                (float)(0.5 * (g.M41 * g.M11 + g.M42 * g.M12 + g.M43 * g.M13 + g.M44 * g.M14))
            ) * (float)m;

            return d2xdt2 * (float)mass2;
        }

        #endregion

            #region elastic energy
        public double CalculateElasticEnergy(double springConstant, double displacement)
        {
            // Ensure the displacement is positive
            double positiveDisplacement = Math.Abs(displacement);

            // Calculate the elastic energy using the formula: E = 0.5 * k * x^2
            double elasticEnergy = 0.5 * springConstant * Math.Pow(positiveDisplacement, 2);

            return elasticEnergy;
        }
        //this is wrong, but i kept it anyway
        public double CalculateElasticEnergyComplex(double springConstantNPerMeter, double displacementMeters, double wireDiameterMeters, double coilCount)
        {
            // Ensure the displacement is positive
            double positiveDisplacement = Math.Abs(displacementMeters);

            // Calculate the effective spring constant based on the spring properties
            double effectiveSpringConstant = (4.0 * springConstantNPerMeter * Math.Pow(wireDiameterMeters, 3)) * (coilCount * coilCount) / (wireDiameterMeters + coilCount * wireDiameterMeters);

            // Calculate the elastic energy using the formula: E = 0.5 * k * x^2
            double elasticEnergyJoules = 0.5 * effectiveSpringConstant * Math.Pow(positiveDisplacement, 2);

            return elasticEnergyJoules;
        }
        #endregion

            #region resistances 
            public double AirResistanceEnergy(double radius, double weight, double velocity, double airDensity, double dragCoefficient)
            {
                // Calculate cross-sectional area
                double crossSectionalArea = Math.PI * Math.Pow(radius, 2);

                // Calculate air resistance energy
                double airResistanceEnergy = (-0.5) * airDensity * crossSectionalArea * Math.Pow(velocity, 2) * dragCoefficient;

                return airResistanceEnergy;
            }
            public double FrictionResistance(double normalForce, double FrictionCoefficient)
            {
                double staticFriction = FrictionCoefficient * normalForce;
                return staticFriction;
            }
            public double FluidResistance(double velocity, double crossSectionalArea, double dragCoefficient, double fluidDensity)
            {
                double fluidResistance = 0.5 * dragCoefficient * fluidDensity * Math.Pow(velocity, 2) * crossSectionalArea;
                return fluidResistance;
            }
            public double ElectricalResistance(double resistivity, double length, double crossSectionalArea)
            {
                double resistance = (resistivity * length) / crossSectionalArea;
                return resistance;
            }
            public double MagneticResistance(double magneticReluctanceCoefficient, double magneticFlux)
            {
                double magneticResistance = magneticReluctanceCoefficient * magneticFlux;
                return magneticResistance;
            }
            public double ThermalResistance(double thermalConductivity, double thickness, double crossSectionalArea)
            {
                double thermalResistance = thickness / (thermalConductivity * crossSectionalArea);
                return thermalResistance;
            }
            public double SoundAbsorptionCoefficient(double totalSurfaceArea, double absorptionArea, double reverberationTime)
            {
                double soundAbsorptionCoefficient = (absorptionArea / totalSurfaceArea) * (1.0 - Math.Exp(-0.161 * reverberationTime));
                return soundAbsorptionCoefficient;
            }
            public double AxialResistance(double crossSectionalArea, double yieldStrength)
            {
                double axialResistance = crossSectionalArea * yieldStrength;
                return axialResistance;
            }
            #endregion

            #region thermodynamics
            public double CalculateHeatTransfer(double temperatureDifference, double thermalConductivity)
            {
                double heatTransfer = temperatureDifference * thermalConductivity;
                return heatTransfer;
            }
            public double CalculateEntropy(double temperature, double heatTransfer)
            {
                double entropyChange = heatTransfer / temperature;
                return entropyChange;
            }
            public double CalculateIdealGasLaw(double pressure, double volume, double temperature)
            {
                double result = (pressure * volume) / temperature;
                return result;
            }
            public double CalculateInternalEnergyThermo(double heatTransfer, double workDone)
            {
                double internalEnergyChange = heatTransfer - workDone;
                return internalEnergyChange;
            }
        #endregion

            #region optics
            public double CalculateRefractionAngle(double angleOfIncidence, double refractiveIndex1, double refractiveIndex2)
            {
                double refractionAngle = Math.Asin(refractiveIndex1 / refractiveIndex2 * Math.Sin(angleOfIncidence));
                return refractionAngle;
            }
            public bool CalculateTotalInternalReflection(double angleOfIncidence, double refractiveIndex1, double refractiveIndex2)
            {
                double criticalAngle = Math.Asin(refractiveIndex2 / refractiveIndex1);
                return angleOfIncidence > criticalAngle;
            }
            public double CalculateLensPower(double focalLength)
            {
                double lensPower = 1.0 / focalLength;
                return lensPower;
            }
            public double CalculateMagnification(double objectHeight, double imageHeight)
            {
                double magnification = imageHeight / objectHeight;
                return magnification;
            }
            public double CalculateImageDistance(double objectDistance, double focalLength)
            {
                double imageDistance = 1 / ((1 / focalLength) - (1 / objectDistance));
                return imageDistance;
            }
            public double CalculateThinLensEquation(double focalLength, double imageDistance = double.NaN, double objectDistance = double.NaN)
            {
                if (double.IsNaN(imageDistance))
                {
                    // Calculate image distance if it is not provided
                    imageDistance = 1 / ((1 / focalLength) - (1 / objectDistance));
                    return imageDistance;
                }
                else if (double.IsNaN(objectDistance))
                {
                    // Calculate object distance if it is not provided
                    objectDistance = 1 / ((1 / focalLength) - (1 / imageDistance));
                    return objectDistance;
                }
                else
                {
                    // Invalid input, both image distance and object distance are provided
                    throw new ArgumentException("Provide either the image distance or the object distance, not both.");
                }
            }
            public double CalculateMirrorEquation(double focalLength, double imageDistance = double.NaN, double objectDistance = double.NaN)
            {
                if (double.IsNaN(imageDistance))
                {
                    // Calculate image distance if it is not provided
                    imageDistance = 1 / ((1 / focalLength) - (1 / objectDistance));
                    return imageDistance;
                }
                else if (double.IsNaN(objectDistance))
                {
                    // Calculate object distance if it is not provided
                    objectDistance = 1 / ((1 / focalLength) - (1 / imageDistance));
                    return objectDistance;
                }
                else
                {
                    // Invalid input, both image distance and object distance are provided
                    throw new ArgumentException("Provide either the image distance or the object distance, not both.");
                }
            }
            public double CalculateSnellLaw(double incidentAngle, double refractiveIndex1, double refractiveIndex2)
            {
                double refractedAngle = Math.Asin(refractiveIndex1 / refractiveIndex2 * Math.Sin(incidentAngle));
                return refractedAngle;
            }
            public double CalculateCriticalAngle(double refractiveIndex1, double refractiveIndex2)
            {
                double criticalAngle = Math.Asin(refractiveIndex2 / refractiveIndex1);
                return criticalAngle;
            }
            public double CalculateDispersion(double refractiveIndexRed, double refractiveIndexBlue, double refractiveIndexYellow)
            {
                double dispersion = (refractiveIndexYellow - 1) / (refractiveIndexBlue - refractiveIndexRed);
                return dispersion;
            }
            #endregion
        
        
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
        public Vector3 MultiplyVectorByDouble(Vector3 v, double d)
        {
            return new Vector3(v.X * ((float)d), v.Y * ((float)d), v.Z * ((float)d));
        }
        #endregion

        #region user input
        public object exampleBasedOnType(Type t)
        {
            if (t == typeof(double))
            {
                return  Math.Round(Random(0,10), 2);
            }
            else if (t == typeof(Vector3))
            {
                return $"({Math.Round(Random(0, 1), 1)},{Math.Round(Random(0, 1), 1)},{Math.Round(Random(0, 1), 1)})";
            }
            else
            {
                return "NAN";
            }
        }
        public string userInput(string input)
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
                            return "("+ deleteSymbolsWithSpace(returnVector.ToString()) + ")";
                        }
                        else if (returnType == typeof(Vector2))
                        {
                            // Cast the returned value to Vector3
                            Vector2 returnVector = (Vector2)result;
                            return "(" + deleteSymbolsWithSpace(returnVector.ToString()) + ")";
                        }
                        else if(returnType == typeof(string))
                        {
                            string returnString = (string)result;
                            return returnString;
                        }
                        else if (returnType == typeof(string[]))
                        {
                            string s = "\n";
                            foreach (string obj in (string[])result)
                            {
                                s+= obj + "\n";
                            }
                            return s;
                        }
                        else if (returnType == typeof(double[]))
                        {
                            string s = "\n";
                            foreach (double obj in (double[])result)
                            {
                                s += obj.ToString() + "\n";
                            }
                            return s;
                        }
                        else if (returnType == typeof(Vector3[]))
                        {
                            string s = "\n";
                            foreach (Vector3 obj in (Vector3[])result)
                            {
                                s += "("+deleteSymbolsWithSpace(obj.ToString()) + ")\n";
                            }
                            return s;
                        }
                        else if (returnType == typeof(Vector2[]))
                        {
                            string s = "\n";
                            foreach (Vector2 obj in (Vector2[])result)
                            {
                                s += "(" + deleteSymbolsWithSpace(obj.ToString()) + ")\n";
                            }
                            return s;
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
                        string s = "";
                        foreach (ParameterInfo pi in p)
                        {
                            parametersNeeded += pi.Name + ": " + pi.ParameterType.ToString() + " \n";
                            
                            s += exampleBasedOnType(pi.ParameterType) + " ";
                        }
                        return "Invalid command, " + $"variables needed for the function {method.Name}:" + 
                            parametersNeeded + $"\n{p.Length} parameters is needed" + 
                            $"\n example: {method.Name} {s}";

                    }

                }
                catch
                {
                    return "Invalid command";

                }

            }
            return "";
        }

        public string getAllVariables()
        {
            Type type = typeof(PhysicsEquations);
            FieldInfo[] fields = type.GetFields(BindingFlags.Public | BindingFlags.NonPublic | BindingFlags.Instance | BindingFlags.Static);
            Dictionary<string,object> variables = new Dictionary<string, object>();
            foreach (FieldInfo field in fields)
            {
                variables.Add(field.Name ,field.GetValue(this));
            }
            string s = "\n";
            foreach (var item in variables)
            {
                s += item.Key + ": " + item.Value + "\n";
            }
            return s;
        }
        public string getAllMethods()
        {
            Type type = typeof(PhysicsEquations);
            MethodInfo[] methods = type.GetMethods(BindingFlags.Public | BindingFlags.NonPublic | BindingFlags.Instance);
            List<string> variables = new List<string>();
            foreach (MethodInfo method in methods)
            {
                variables.Add(method.Name);
            }
            string s = "\n";
            foreach (var item in variables)
            {
                s += item + "\n";
            }
            return s;
        }
        #endregion
        public string[] GetStringArrayWithOnes(int length)
        {
            string[] array = new string[length];
            for (int i = 0; i < length; i++)
            {
                array[i] = "1";
            }
            return array;
        }
        public Vector3[] GetVector3ArrayWithOnes(int length)
        {
            Vector3[] array = new Vector3[length];
            Vector3 onesVector = new Vector3(1, 1, 1);

            for (int i = 0; i < length; i++)
            {
                array[i] = onesVector;
            }

            return array;
        }
    }
}
