using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;

namespace testing_stuff_place_2___GR
{
    public class GR
    {
        public const double G = 6.67430e-11;
        public const double c = 299792458;



        public static double[,] stressAndEnergySolve(double[,] metric)
        {
            double[,,] christoffel = CalculateChristoffelSymbolsAdvanced(metric);
            double[,,,] riemann = CalculateRiemannTensor(christoffel);
            double[,] ricci = ricciTensorAdvanced(riemann);
            double R = calculateRscalarAdvanced(metric, ricci);

            double[,] stressAndEnergy = new double[4, 4];
            double dconst = (8 * Math.PI * G) / Math.Pow(c, 4);
            stressAndEnergy = divideMatrix(SubtractMatrix(ricci, multiplyMatrix(multiplyMatrix(metric, -0.5), R)) , dconst);
            
            return stressAndEnergy;
        }

        public static double[,] CalculateGeodesics(double[,] metric, double[] initialPosition, double[] initialVelocity, double stepSize, double finalLambda)
        {
            double[,,] christoffelSymbols = CalculateChristoffelSymbolsAdvanced(metric);
            int dimension = initialPosition.Length;
            double[,] geodesics = new double[dimension, (int)(finalLambda / stepSize) + 1];
            double[] currentPosition = (double[])initialPosition.Clone();
            double[] currentVelocity = (double[])initialVelocity.Clone();
            double currentLambda = 0.0;
            int index = 0;

            while (currentLambda <= finalLambda)
            {
                for (int a = 0; a < dimension; a++)
                {
                    geodesics[a, index] = currentPosition[a];
                }
                currentLambda += stepSize;
                index++;

                double[] acceleration = new double[dimension];

                for (int c = 0; c < dimension; c++)
                {
                    for (int a = 0; a < dimension; a++)
                    {
                        for (int b = 0; b < dimension; b++)
                        {
                            acceleration[c] -= christoffelSymbols[c, a, b] * currentVelocity[a] * currentVelocity[b];
                        }
                    }
                }

                for (int a = 0; a < dimension; a++)
                {
                    currentVelocity[a] += acceleration[a] * stepSize;

                    currentPosition[a] += currentVelocity[a] * stepSize;
                }
            }

            return geodesics;
        }
        #region metric tensors
        public static double[,] schwarzschildMetricTensor(double mass, double r, double theta)
        {
            double[,] metric = new double[4,4];
            double rs = (2 * G * mass)/Math.Pow(c,2);

            
            metric[0,0] = -(1 - (rs/r));
            metric[1,1] = 1/(1 - (rs/r));
            metric[2,2] = Math.Pow(r,2);
            metric[3,3] = Math.Pow(r,2) * Math.Pow(Math.Sin(theta),2);


            return metric;
        }
        #endregion

        #region curvature
        public static double[,,] CalculateChristoffelSymbolsAdvanced(double[,] metricTensor)
        {
            int dimension = metricTensor.GetLength(0);
            double[,] inverseMetricTensor = CalculateInverseMetricTensor(metricTensor);
            double[,,] christoffelSymbols = new double[dimension, dimension, dimension];

            for (int c = 0; c < dimension; c++)
            {
                for (int a = 0; a < dimension; a++)
                {
                    for (int b = 0; b < dimension; b++)
                    {
                        double sum = 0.0;

                        for (int d = 0; d < dimension; d++)
                        {
                            sum += 0.5 * inverseMetricTensor[c, d] * (
                                PartialDerivative(metricTensor, a, b) +
                                PartialDerivative(metricTensor, b, a) -
                                PartialDerivative(metricTensor, a, d)
                            );
                        }

                        christoffelSymbols[c, a, b] = sum;
                    }
                }
            }

            return christoffelSymbols;
        }
        public static double[,,] CalculateChristoffelSymbols(double[,] metricTensor)
        {
            int dimension = metricTensor.GetLength(0);
            double[,,] christoffelSymbols = new double[dimension, dimension, dimension];

            for (int c = 0; c < dimension; c++)
            {
                for (int a = 0; a < dimension; a++)
                {
                    for (int b = 0; b < dimension; b++)
                    {
                        double sum = 0.0;

                        
                            sum += (1 / (2 * metricTensor[c, c])) * (
                                PartialDerivative(metricTensor, c, a) +
                                PartialDerivative(metricTensor, c, b) -
                                PartialDerivative(metricTensor, a, b)
                            );
                        

                        christoffelSymbols[c, a, b] = sum;
                    }
                }
            }

            return christoffelSymbols;
        }
        public static double[,,,] CalculateRiemannTensor(double[,,] christoffelSymbols)
        {
            int dimension = christoffelSymbols.GetLength(0);
            double[,,,] riemannTensor = new double[dimension, dimension, dimension, dimension];

            for (int a = 0; a < dimension; a++)
            {
                for (int b = 0; b < dimension; b++)
                {
                    for (int c = 0; c < dimension; c++)
                    {
                        for (int d = 0; d < dimension; d++)
                        {
                            double sum = 0.0;

                            for (int e = 0; e < dimension; e++)
                            {
                                sum += christoffelSymbols[c, a, e] * christoffelSymbols[e, b, d] -
                                       christoffelSymbols[c, b, e] * christoffelSymbols[e, a, d];
                            }

                            riemannTensor[a, b, c, d] = sum;
                        }
                    }
                }
            }

            return riemannTensor;
        }

        #endregion

        #region ricci and R scalar
        public static double[,] ricciTensor(double[,,,] riemannTensor)
        {
            double[,] ricci = new double[4,4];


            //temp simplify
            ricci[0, 0] = riemannTensor[0, 0, 0, 0] + riemannTensor[1, 0, 1, 0];
            ricci[0, 1] = riemannTensor[0, 0, 0, 1] + riemannTensor[1, 0, 1, 1];
            ricci[1, 0] = riemannTensor[0, 1, 0, 0] + riemannTensor[1, 1, 1, 0];
            ricci[1, 1] = riemannTensor[0, 1, 0, 1] + riemannTensor[1, 1, 1, 1];

            //a more complicated but how its actually done, but it doesnt makea difference, as most of the riemann arent really used
            /* for (int i = 0; i < riemannTensor.GetLength(0); i++)
             {
                 for (int j = 0; j < riemannTensor.GetLength(1); j++)
                 {
                     for (int k = 0; k < riemannTensor.GetLength(2); k++)
                     {
                         for (int l = 0; l < riemannTensor.GetLength(3); l++)
                         {
                             if(i == k)
                             {

                             }
                         }
                     }
                 }
             }*/

            return ricci;
        }
        public static double[,] ricciTensorAdvanced(double[,,,] riemannTensor)
        {
            int dimension = riemannTensor.GetLength(0);
            double[,] ricciTensor = new double[dimension, dimension];

            for (int a = 0; a < dimension; a++)
            {
                for (int b = 0; b < dimension; b++)
                {
                    double sum = 0.0;

                    for (int c = 0; c < dimension; c++)
                    {
                        sum += riemannTensor[a, c, b, c];
                    }

                    ricciTensor[a, b] = sum;
                }
            }

            return ricciTensor;
        }
        public static double calculateRscalar(double[,] metric, double[,] ricci) {
            return (metric[0, 0] * ricci[0, 0] + metric[1, 1] * ricci[1, 1] + metric[2, 2] * ricci[2, 2] + metric[3, 3] * ricci[3, 3]);
        }
        public static double calculateRscalarAdvanced(double[,] metric, double[,] ricci)
        {
            double scalarCurvature = 0.0;

            int dimension = metric.GetLength(0);
            for (int i = 0; i < dimension; i++)
            {
                scalarCurvature += metric[i, i] * ricci[i, i];
            }

            return scalarCurvature;
        }

        #endregion

        #region helper functions
        public static double[,] CalculateInverseMetricTensor(double[,] metricTensor)
        {
            int dimension = metricTensor.GetLength(0);
            Matrix<double> matrix = Matrix<double>.Build.DenseOfArray(metricTensor);
            Matrix<double> inverseMatrix = matrix.Inverse();
            double[,] inverseMetricTensor = inverseMatrix.ToArray();

            return inverseMetricTensor;
        }
        public static double PartialDerivative(double[,] tensor, int a, int b)
        {
            double h = 1;

            double f1 = GetTensorValue(tensor, a + (int)h, b);
            double f2 = GetTensorValue(tensor, a - (int)h, b);

            double der = (f1 - f2) / (2 * h);

            return der;
        }

        // Helper method to get value from tensor
        private static double GetTensorValue(double[,] tensor, int i, int j)
        {
            if (i < 0 || i >= tensor.GetLength(0) ||
               j < 0 || j >= tensor.GetLength(1))
            {
                return 0;
            }

            return tensor[i, j];
        }
        public static double[,] multiplyMatrix(double[,] matrix1, double a)
        {
            double[,] newm = new double[matrix1.GetLength(0), matrix1.GetLength(1)];

            Console.WriteLine(newm.GetLength(0));
            Console.WriteLine(matrix1.Length);
            for (int i = 0; i < newm.GetLength(0) - 1; i++)
            {
                for (int k = 0; k < newm.GetLength(1) - 1; k++)
                {
                    newm[i, k] = matrix1[i, k] * a;
                }
            }
            return newm;
        }
        public static double[,] divideMatrix(double[,] matrix1, double a)
        {
            double[,] newm = new double[matrix1.GetLength(0), matrix1.GetLength(1)];

            for (int i = 0; i < newm.GetLength(0); i++)
            {
                for (int k = 0; k < newm.GetLength(1); k++)
                {
                    newm[i, k] = matrix1[i, k] / a;
                }
            }
            return newm;
        }
        public static double[,] SubtractMatrix(double[,] matrix1, double[,] matrix2)
        {
            double[,] newm = new double[matrix1.GetLength(0), matrix1.GetLength(1)];

            for (int i = 0; i < newm.GetLength(0); i++)
            {
                for (int k = 0; k < newm.GetLength(1); k++)
                {
                    newm[i, k] = matrix1[i, k] - matrix2[i, k];
                }
            }
            return newm;
        }
        #endregion
    }
}
