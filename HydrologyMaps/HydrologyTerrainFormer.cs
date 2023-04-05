using System.Numerics;

namespace HydrologyMaps;

public class HydrologyTerrainFormer
{
    private static FastNoiseLite _noiseGenerator = new();

    public HydrologyTerrainFormer()
    {
        _noiseGenerator.SetNoiseType(FastNoiseLite.NoiseType.Perlin);
        _noiseGenerator.SetFrequency(0.003f);
        _noiseGenerator.SetFractalGain(0);
    }

    public static void CarveRivers(float[,] heightMap, List<DirectedNode> riverRoots, float riverDepth,
        float baseRiverWidth, float riverWidthFlowInfluence, float falloffFactor, Random random)
    {
        foreach (var root in riverRoots)
        {
            CarveRiverDFS(heightMap, root, riverDepth, baseRiverWidth, riverWidthFlowInfluence, falloffFactor,  random);
        }
    }


    private static void CarveRiverDFS(float[,] heightMap, DirectedNode currentNode, float riverDepth, float baseRiverWidth, float riverWidthFlowInfluence, float falloffFactor, Random random)
    {
        if (currentNode.Children.Count == 0) return;

        foreach (var childGraphNode in currentNode.Children)
        {
            var childNode = (DirectedNode)childGraphNode;

            // Generate a smoothed path with random twists and turns
            List<Point> smoothedPath = GenerateSmoothedPath(currentNode.Point, childNode.Point, random);

            for (int i = 0; i < smoothedPath.Count - 1; i++)
            {
                Point start = smoothedPath[i];
                Point end = smoothedPath[i + 1];

                double startFlowRate = currentNode.FlowRate / 25;
                double endFlowRate = childNode.FlowRate / 25;

                float startWidth = Lerp(baseRiverWidth, baseRiverWidth * (float)Math.Pow(startFlowRate, riverWidthFlowInfluence), riverWidthFlowInfluence);
                float endWidth = Lerp(baseRiverWidth, baseRiverWidth * (float)Math.Pow(endFlowRate, riverWidthFlowInfluence), riverWidthFlowInfluence);

                float depth = Lerp((float)(riverDepth ), (float)(riverDepth), (float)i / (smoothedPath.Count - 1));
                
                CarveRiverSegment(heightMap, start, end, depth, startWidth, endWidth, falloffFactor, random);
            }

            CarveRiverDFS(heightMap, childNode, riverDepth, baseRiverWidth, riverWidthFlowInfluence, falloffFactor, random);
        }
    }

    private static List<Point> GenerateSmoothedPath(Point start, Point end, Random random)
    {
        List<Point> path = new List<Point> { start, end };
        List<Point> smoothedPath = new List<Point>();

        Vector2 direction = (end - start).ToVector();
        direction /= direction.Length();

        // Add random control points to the path for twists and turns
        int numControlPoints = random.Next(1, 4);
        for (int i = 0; i < numControlPoints; i++)
        {
            float t = (float)i / (numControlPoints + 1);
            Point controlPoint = Point.Lerp(start, end, t);

            // Create a random angle based on the current direction (between -45 and 45 degrees)
            float angle = random.Next(-90, 90) * (float)Math.PI / 180f;

            // Calculate the new direction with the random angle applied
            Vector2 newDirection = new Vector2(
                direction.X * (float)Math.Cos(angle) - direction.Y * (float)Math.Sin(angle),
                direction.X * (float)Math.Sin(angle) + direction.Y * (float)Math.Cos(angle)
            );

            // Calculate the random offset using the new direction
            Vector2 randomOffset = newDirection * random.Next(1, 5);

            controlPoint += Point.FromVector(randomOffset);
            path.Insert(1 + i, controlPoint);
        }

        // Smooth the path using Catmull-Rom splines
        for (float t = 0; t < path.Count - 1; t += 0.5f)
        {
            int index = (int)Math.Floor(t);
            Point p0 = path[Math.Clamp(index - 1, 0, path.Count - 1)];
            Point p1 = path[index];
            Point p2 = path[index + 1];
            Point p3 = path[Math.Clamp(index + 2, 0, path.Count - 1)];

            Vector2 point = CatmullRomSpline(p0.ToVector(), p1.ToVector(), p2.ToVector(), p3.ToVector(), t - index);
            smoothedPath.Add(Point.FromVector(point));
        }

        smoothedPath.Add(path[^1]);
        return smoothedPath;
    }

    private static void CarveRiverSegment(float[,] heightMap, Point start, Point end, float depth, float startWidth,
        float endWidth, float falloffFactor, Random random)
    {
        // Calculate the line equation between start and end points
        Vector2 direction = (end - start).ToVector();
        int numSteps = (int)direction.Length();
        direction /= numSteps;

        // Iterate over the points along the line between start and end points
        for (int i = 0; i <= numSteps; i++)
        {
            Point position = start + Point.FromVector(direction * i);
            int x = Math.Clamp(position.X, 0, heightMap.GetLength(0) - 1);
            int y = Math.Clamp(position.Y, 0, heightMap.GetLength(1) - 1);

            float t = (float)i / numSteps;
            float baseWidth = Lerp(startWidth, endWidth, t);

            // Apply Perlin noise to the river width
            _noiseGenerator.SetSeed(random.Next(10000));
            float perlinValue = 2 * Math.Abs(_noiseGenerator.GetNoise(position.X / 1f, position.Y / 1f, 0));
            float currentWidth = baseWidth * (0.7f + 0.3f * perlinValue);

            // Update the heightMap value by subtracting the depth and considering the width
            for (int offsetX = -((int)currentWidth); offsetX <= (int)currentWidth; offsetX++)
            {
                for (int offsetY = -((int)currentWidth); offsetY <= (int)currentWidth; offsetY++)
                {
                    int nx = Math.Clamp(x + offsetX, 0, heightMap.GetLength(0) - 1);
                    int ny = Math.Clamp(y + offsetY, 0, heightMap.GetLength(1) - 1);
                    float distance = (float)Math.Sqrt(offsetX * offsetX + offsetY * offsetY);


                    if (distance <= currentWidth)
                    {
                        // Calculate depth falloff using the falloffFactor
                        float depthFalloff = 1 - (float)Math.Pow(distance / currentWidth, falloffFactor);
                        float currentDepth = depth * depthFalloff;

                        heightMap[nx, ny] = Math.Max(heightMap[nx, ny] - currentDepth, 0);
                    }
                }
            }
        }
    }
    
    // Calculate the Catmull-Rom spline for a set of control points and a parameter t
    private static Vector2 CatmullRomSpline(Vector2 p0, Vector2 p1, Vector2 p2, Vector2 p3, float t)
    {
        Vector2 a = 2f * p1;
        Vector2 b = p2 - p0;
        Vector2 c = 2f * p0 - 5f * p1 + 4f * p2 - p3;
        Vector2 d = -p0 + 3f * p1 - 3f *
            p2 + p3;

        return 0.5f * (a + (b * t) + (c * t * t) + (d * t * t * t));
    }

    public static float Lerp(float a, float b, float t)
    {
        return a + (b - a) * t;
    }
}