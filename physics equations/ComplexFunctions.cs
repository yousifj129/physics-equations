using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace physics_equations
{
    internal class ComplexFunctions
    {
    }
    public class PartialDerivativesCalculator
    {
        public double[] CalculatePartialDerivatives(double[] variables, double[] metricComponents)
        {
            // Define the metric tensor components as functions of the variables
            double g11 = metricComponents[0];
            double g12 = metricComponents[1];
            double g13 = metricComponents[2];
            double g14 = metricComponents[3];
            // Define the variables
            double x = variables[0];
            double y = variables[1];
            double z = variables[2];

            // Calculate the partial derivatives
            double partialDerivativeX = CalculatePartialDerivative(x, y, z, g11, variables, 0);
            double partialDerivativeY = CalculatePartialDerivative(x, y, z, g11, variables, 1);
            double partialDerivativeZ = CalculatePartialDerivative(x, y, z, g11, variables, 2);

            // Create an array to hold the partial derivatives
            double[] partialDerivatives = new double[3];
            partialDerivatives[0] = partialDerivativeX;
            partialDerivatives[1] = partialDerivativeY;
            partialDerivatives[2] = partialDerivativeZ;

            return partialDerivatives;
        }

        private double CalculatePartialDerivative(double x, double y, double z, double g11, double[] variables, int index)
        {
            double epsilon = 1e-6; // Small value for numerical differentiation

            // Perturb the variable at the specified index by epsilon
            double[] perturbedVariables = (double[])variables.Clone();
            perturbedVariables[index] += epsilon;

            // Calculate the function values at the original and perturbed variables
            double originalValue = CalculateFunction(x, y, z, variables, g11);
            double perturbedValue = CalculateFunction(x, y, z, perturbedVariables, g11);

            // Calculate the partial derivative using the central difference formula
            double partialDerivative = (perturbedValue - originalValue) / epsilon;

            return partialDerivative;
        }
        
        
        private double CalculateFunction(double x, double y, double z, double[] variables, double g11)
        {
            // Define the function based on the metric tensor component
            double functionValue = g11 * x * x + y * y + z * z + variables[0] * variables[0];

            return functionValue;
        }
    }
    public class GeneralRelativityCalculator
    {
        // Calculate magnitude of Vector3
        public static double Vector3Magnitude(Vector3 vector)
        {
            return Math.Sqrt(Vector3.Dot(vector, vector));
        }

        // Normalize Vector3 
        public static Vector3 Vector3Normalize(Vector3 vector)
        {
            float magnitude = ((float)Vector3Magnitude(vector));

            // Create a Vector3 with the magnitude value
            Vector3 magnitudeVec = new Vector3(magnitude, magnitude, magnitude);

            return Vector3.Divide(vector, magnitudeVec);
        }
        public static double CalculateSchwarzschildRadius(double mass, double G, double c)
        {
            double rs = (2 * G * mass) / (c * c);
            return rs;
        }

        public static double[] CalculateSchwarzschildMetric(double r, double rs)
        {
            // Metric tensor components
            double g00 = -(1 - rs / r);
            double g11 = 1 / (1 - rs / r);
            double g22 = Math.Pow(r, 2);
            double g33 = Math.Pow(r, 2);

            return new double[] { g00, g11, g22, g33 };
        }
        public static double[,] KerrMetric(double M, double a, double r, double theta, double phi)
        {
            double rho2 = r * r + a * a * Math.Cos(theta) * Math.Cos(theta);
            double delta = r * r - 2 * M * r + a * a;
            double g = r * r - a * a * Math.Cos(theta) * Math.Cos(theta);

            double[,] metric = new double[4, 4];

            metric[0, 0] = -(1 - 2 * M * r / rho2);
            metric[0, 1] = 0;
            metric[0, 2] = 0;
            metric[0, 3] = -a * Math.Sin(theta) * Math.Sin(theta) * 2 * M * r / rho2;

            metric[1, 0] = 0;
            metric[1, 1] = rho2 / delta;
            metric[1, 2] = 0;
            metric[1, 3] = 0;

            metric[2, 0] = 0;
            metric[2, 1] = 0;
            metric[2, 2] = rho2;
            metric[2, 3] = 0;

            metric[3, 0] = -a * Math.Sin(theta) * Math.Sin(theta) * 2 * M * r / rho2;
            metric[3, 1] = 0;
            metric[3, 2] = 0;
            metric[3, 3] = g / rho2;

            return metric;
        }

        public static double CalculateGravityForceSchwarzschild(double mass1, double mass2, double distance, double G, double c)
        {
            // Calculate Schwarzschild radii
            double rs1 = CalculateSchwarzschildRadius(mass1, G, c);
            double rs2 = CalculateSchwarzschildRadius(mass2, G, c);

            // Calculate force
            double force = G * (mass1 * mass2) / (distance * distance);

            // Adjust for Schwarzschild radii
            force *= (1 - rs1 / distance) * (1 - rs2 / distance);

            return force;
        }
        public static double CalculateGravityForceKerr(double mass1, double mass2, double distance,double theta, double phi, double G, double c)
        {
            // Calculate the metric components.
            double[,] metric = KerrMetric(mass1, 0, distance, theta, phi);

            // Calculate the Schwarzschild radii.
            double rs1 = CalculateSchwarzschildRadius(mass1, G, c);
            double rs2 = CalculateSchwarzschildRadius(mass2, G, c);

            // Calculate force
            double force = G * (mass1 * mass2) / (distance * distance);

            // Adjust for Schwarzschild radii and metric.
            force *= (1 - rs1 / distance) * (1 - rs2 / distance) * metric[0,0];

            return force;
        }
        public static double CalculateGravityForceKerrNewman(double mass1, double mass2, double distance, double theta, double phi, double charge1, double charge2, double G, double c)
        {
            // Calculate the metric components.
            double[,] metric = KerrMetric(mass1, 0, distance, theta, phi);

            // Calculate the Schwarzschild radius and electric potential.
            double rs1 = CalculateSchwarzschildRadius(mass1, G, c);
            double rs2 = CalculateSchwarzschildRadius(mass2, G, c);
            double potential1 = charge1 / (4 * Math.PI * c * distance);
            double potential2 = charge2 / (4 * Math.PI * c * distance);

            // Calculate gravitational force.
            double force = G * (mass1 * mass2) / (distance * distance);

            // Adjust for Schwarzschild radii, metric, and electromagnetic effects.
            force *= (1 - rs1 / distance) * (1 - rs2 / distance) * metric[0, 0] + potential1 * potential2;

            return force;
        }

        public static Vector3 CalculateGravityVelocitySchwarzschild(Vector3 pos1, double mass1,
                                               Vector3 pos2, double mass2,
                                               double G, double c)
        {
            // Calculate distance vector between objects
            Vector3 distanceVec = pos2 - pos1;
            double distance = Vector3Magnitude(distanceVec);

            // Calculate Schwarzschild radii
            double rs1 = CalculateSchwarzschildRadius(mass1, G, c);
            double rs2 = CalculateSchwarzschildRadius(mass2, G, c);

            // Calculate force
            double force = G * (mass1 * mass2) / (distance * distance);

            // Adjust for Schwarzschild radii
            force *= (1 - rs1 / distance) * (1 - rs2 / distance);

            // Calculate acceleration from force
            Vector3 accel = Vector3Normalize(distanceVec) * ((float)(force / mass1));

            // Integrate to get velocity
            Vector3 velocity = accel;

            return velocity;
        }

        public static Vector3 CalculateGravityVelocityKerr(Vector3 pos1, double mass1,
                                               Vector3 pos2, double mass2,double theta, double phi,
                                               double G, double c)
        {
            Vector3 distanceVec = pos2 - pos1;
            double distance = Vector3Magnitude(distanceVec);
            // Calculate the metric components.
            double[,] metric = KerrMetric(mass1, 0, distance, theta, phi);

            // Calculate the Schwarzschild radii.
            double rs1 = CalculateSchwarzschildRadius(mass1, G, c);
            double rs2 = CalculateSchwarzschildRadius(mass2, G, c);

            // Calculate force
            double force = G * (mass1 * mass2) / (distance * distance);

            // Adjust for Schwarzschild radii and metric.
            force *= (1 - rs1 / distance) * (1 - rs2 / distance) * metric[0, 0];

            // Calculate acceleration from force
            Vector3 accel = Vector3Normalize(distanceVec) * ((float)(force / mass1));

            // Integrate to get velocity
            Vector3 velocity = accel;

            return velocity;
        }

    }
    public class PartialDerivativesRichardsonExtrapolationCalculator
    {
        public double[] CalculatePartialDerivatives(double[] variables, double[] metricComponents)
        {
            // Define the metric tensor components as functions of the variables
            double g11 = metricComponents[0];
            double g12 = metricComponents[1];
            double g13 = metricComponents[2];
            double g14 = metricComponents[3];
            // Define the variables
            double x = variables[0];
            double y = variables[1];
            double z = variables[2];

            // Calculate the partial derivatives using Richardson extrapolation
            double partialDerivativeX = CalculatePartialDerivativeRichardson(x, y, z, g11, variables, 0);
            double partialDerivativeY = CalculatePartialDerivativeRichardson(x, y, z, g11, variables, 1);
            double partialDerivativeZ = CalculatePartialDerivativeRichardson(x, y, z, g11, variables, 2);

            // Create an array to hold the partial derivatives
            double[] partialDerivatives = new double[3];
            partialDerivatives[0] = partialDerivativeX;
            partialDerivatives[1] = partialDerivativeY;
            partialDerivatives[2] = partialDerivativeZ;

            return partialDerivatives;
        }
        public double RichardsonExtrapolation(double y_n, double y_n_1, int k)
        {
        // Calculate the error estimate
        double error = y_n - y_n_1;

        // Calculate the improved approximation
        double y_new = y_n - (error / (2^k - 1));

        return y_new;
        }
        private double CalculatePartialDerivativeRichardson(double x, double y, double z, double g11, double[] variables, int index)
        {
            double epsilon = 1e-6; // Small value for numerical differentiation

            // Perturb the variable at the specified index by epsilon
            double[] perturbedVariables = (double[])variables.Clone();
            perturbedVariables[index] += epsilon;

            // Calculate the function values at the original and perturbed variables
            double originalValue = CalculateFunction(x, y, z, variables, g11);
            double perturbedValue = CalculateFunction(x, y, z, perturbedVariables, g11);

            // Calculate the partial derivative using the Richardson extrapolation method
            double partialDerivative = RichardsonExtrapolation(originalValue, perturbedValue, 1);

            return partialDerivative;
        }

        private double CalculateFunction(double x, double y, double z, double[] variables, double g11)
        {
            // Define the function based on the metric tensor component
            double functionValue = g11 * x * x + y * y + z * z + variables[0] * variables[0];

            return functionValue;
        }
    }


    public class schrodingerEquations
    {
        public static Vector2[] SolveSchrodingerEquation(Vector2[] potential, Vector2[] positionGrid, double mass, double hbar)
        {
            int n = positionGrid.Length;
            Vector2[] wavefunction = new Vector2[n];
            Vector2[] energyLevels = new Vector2[n];

            // Calculate the Laplacian of the wavefunction using finite differences.
            float deltaX = positionGrid[1].X - positionGrid[0].X;
            float deltaY = positionGrid[1].Y - positionGrid[0].Y;
            Vector2[] laplacian = new Vector2[n];
            for (int i = 1; i < n - 1; i++)
            {
                float d2PsiX = (wavefunction[i + 1].X - 2 * wavefunction[i].X + wavefunction[i - 1].X) / (deltaX * deltaX);
                float d2PsiY = (wavefunction[i + 1].Y - 2 * wavefunction[i].Y + wavefunction[i - 1].Y) / (deltaY * deltaY);
                laplacian[i] = new Vector2(d2PsiX, d2PsiY);
            }

            // Solve the Schrödinger equation for each point in the grid.
            for (int i = 0; i < n; i++)
            {
                float potentialEnergy = potential[i].X + potential[i].Y;
                float kineticEnergy = (float)((hbar * hbar) / (2 * mass) * (laplacian[i].X + laplacian[i].Y));
                energyLevels[i] = new Vector2(potentialEnergy, kineticEnergy);
            }

            return energyLevels;
        }
    }
}
