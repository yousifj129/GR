// See https://aka.ms/new-console-template for more information

using testing_stuff_place_2___GR;
double mass = 10e12;
double r = 1000;
double theta = 70;

double[,] metric = GR.schwarzschildMetricTensor(mass, r, theta);
double[,] testingMetric = new double[,]
{
    { 1, 0 , 0 , 0 },
    { 0, -1 , 0 , 0 },
    { 0, 0 , -1 , 0},
    { 0, 0 , 0 , -1}
};
double[] initialPosition = new double[] { 0, 1000, 70, 30};
double[] initialVelocity = new double[] { -0.25, 0, 0, 0 };

// Set integration parameters
double stepSize = 0.1;
double finalLambda = 10.0;

// Calculate geodesics
double[,] geodesics = GR.CalculateGeodesics(metric, initialPosition, initialVelocity, stepSize, finalLambda);

// Print geodesics
Console.WriteLine("Geodesic Coordinates:");
for (int i = 0; i < geodesics.GetLength(1); i++)
{
    Console.WriteLine($"Step {i}: ({geodesics[0, i]}, {geodesics[1, i]}, {geodesics[2, i]}, {geodesics[3, i]})");
}

Console.ReadLine();