// ReSharper disable CompareOfFloatsByEqualityOperator

namespace RancHydrologyMaps;

public static class HeightmapCleaner
{
    // Since the basemap may have been generated using something like perlin noise, there might be small islands or lakes that we want to remove such that we only one big continuous island 
    // It is assumed that the input heightmap only has values 1 and 0 and this point where 0 is ocean and 1 is land, also assumed that artifacts are round-ish and don't touch the edges of the map
    public static void RemoveArtifacts(float[,] heightmap, int maxDiameterOfArtifacts)
    {
        // This method scans through heightmap and looks for round-ish artifacts. It investigates both land pixels and ocean pixels and checks if the next pixel is the first part of a round shape like:
        /*
               Found first --> ### 
                             ########
                            ####X#####     <-- Find center and estimate size
                             #######
                               ##
         */
        // After a candidate artifact is found we do a flood fill but we limited the flood fill since we might have found something like:
        /* == Not actually an artifact ==
         Found first --> ### 
                       ########
                      ####X#####     <-- Find center and estimate size
                       ###############...
                            ##########...
             ...######################...
        */

        int width = heightmap.GetLength(0);

        for (int y = 1; y < width; y++)
        {
            for (int x = 1; x < width; x++)
            {
                float targetValue = heightmap[x, y];

                // Find width of the top row
                int topRowWidth = 0;
                for (int i = x; i < width && heightmap[i, y] == targetValue; i++)
                {
                    topRowWidth++;
                }

                // Too big to be an artifact, assumed to be the main island shape.
                if (topRowWidth > maxDiameterOfArtifacts)
                {
                    x += topRowWidth;
                    continue;
                }

                // Calculate the top-middle pixel
                int topMiddleX = x + topRowWidth / 2;

                // if there is a pixel above the current pixel at the center of the row then this isn't the top of an artifact
                if (heightmap[topMiddleX, y - 1] == targetValue)
                    continue;

                // Find height of the artifact by scanning down from the top-middle pixel
                int artifactHeight = 0;
                for (int j = y + 1; j < width && heightmap[topMiddleX, j] == targetValue; j++)
                {
                    artifactHeight++;
                }

                // Too big to be an artifact, assumed to be the main island shape.
                if (artifactHeight > maxDiameterOfArtifacts)
                    continue;

                if (topRowWidth <= maxDiameterOfArtifacts && artifactHeight <= maxDiameterOfArtifacts)
                {
                    // Perform a flood-fill to remove the artifact. Flood fill might return an empty list if the flood filled too much, see the ascii figure above
                    List<(int, int)> artifactPixels =
                        FloodFill(heightmap, topMiddleX, y, width, targetValue, maxDiameterOfArtifacts);
                    foreach (var (px, py) in artifactPixels)
                    {
                        heightmap[px, py] = 1 - targetValue; // Change to opposite value
                    }
                }
            }
        }
    }
    
    // Reusable data structures
    private static List<(int, int)> _connectedPixels = new List<(int, int)>();
    private static List<(int, int)> _emptyList = new List<(int, int)>();
    private static Queue<(int, int)> _queue = new Queue<(int, int)>();

    private static List<(int, int)> FloodFill(float[,] heightmap, int x, int y, int width,
        float targetValue, int maxDiameterOfArtifacts)
    {
        int maxArea = (int)(Math.PI * Math.Pow(maxDiameterOfArtifacts / 2.0, 2) * 2f);

        bool[,] visited = new bool[width, width];
        _connectedPixels.Clear();
        _queue.Clear();

        _queue.Enqueue((x, y));

        while (_queue.Count > 0)
        {
            (int cx, int cy) = _queue.Dequeue();

            
            if (cx < 0 || cx >= width || cy < 0 || cy >= width || visited[cx, cy] ||  // Assume that no artifacts touch the edges of the map
                heightmap[cx, cy] != targetValue)
            {
                continue;
            }

            visited[cx, cy] = true;
            _connectedPixels.Add((cx, cy));
            if (_connectedPixels.Count > maxArea)
                return _emptyList;

            _queue.Enqueue((cx - 1, cy));
            _queue.Enqueue((cx + 1, cy));
            _queue.Enqueue((cx, cy - 1));
            _queue.Enqueue((cx, cy + 1));
        }

        return _connectedPixels;
    }
}