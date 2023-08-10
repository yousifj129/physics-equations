using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static helperFunctions.HelperFunctions;

namespace physics_equations
{
    public class PhysicsEquations
    {
        #region constants
        // Gravitational constant (units: m^3 kg^-1 s^-2)
        public const double GravitationalConstant = 6.67430e-11;

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
        public double CalculateComplementaryProbability(double probability)
        {
            double complementaryProbability = 1 - probability;
            return complementaryProbability;
        }
        public double CalculateJointProbability(double probabilityA, double probabilityB)
        {
            double jointProbability = probabilityA * probabilityB;
            return jointProbability;
        }
        public double CalculateProbabilityDensity(double wavefunction)
        {
            double probabilityDensity = Math.Pow(wavefunction, 2);
            return probabilityDensity;
        }
        public double CalculateTransitionProbability(double initialWavefunction, double finalWavefunction)
        {
            double transitionProbability = Math.Pow(Math.Abs(initialWavefunction * finalWavefunction), 2);
            return transitionProbability;
        }
        public double CalculateFermiDiracDistribution(double energy, double temperature, double mu)
        {
            double k = 8.617333262145e-5; // Boltzmann constant
            double fermiDiracDistribution = 1 / (Math.Exp((energy - mu) / (k * temperature)) + 1);
            return fermiDiracDistribution;
        }
        #endregion

        #region wave functions

        public double CalculateBoxWaveFunction(double amplitude, double n, double boxLength, double position)
        {
            double wavefunction = amplitude * Math.Sin((int)n * Math.PI * position / boxLength);
            return wavefunction;
        }
        public  double CalculateGaussianWaveFunction(double amplitude, double sigma, double x, double x0)
        {
            double wavefunction = amplitude * Math.Exp(-(Math.Pow(x - x0, 2)) / (2 * Math.Pow(sigma, 2)));
            return wavefunction;
        }
        public double CalculateHarmonicOscillatorWaveFunction(double amplitude, double n, double x, double x0, double k)
        {
            double wavefunction = amplitude * Math.Pow(k / (Math.PI * (int)n), 0.25) * Math.Exp(-(Math.Pow(x - x0, 2) * k) / (2 * (int)n));
            return wavefunction;
        }
        public double CalculateCoulombWaveFunction(double amplitude, double ls, double k, double eta, double radius)
        {
            int l = (int)ls;
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
    }
}
