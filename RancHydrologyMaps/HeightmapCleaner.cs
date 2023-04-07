namespace RancHydrologyMaps;

public static class HeightmapCleaner
{
    public static void RemoveArtifacts(float[,] heightmap, int maxDiameterOfArtifacts)
    {
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

                if (topRowWidth > maxDiameterOfArtifacts)
                {
                    x += topRowWidth;
                    continue;
                }


                // Calculate the top-middle pixel
                int topMiddleX = x + topRowWidth / 2;

                if (heightmap[topMiddleX, y - 1] == targetValue)
                    continue;
                
                // Find height of the artifact by scanning up and down from the top-middle pixel
                int artifactHeight = 0;
                for (int j = y + 1; j < width && heightmap[topMiddleX, j] == targetValue; j++)
                {
                    artifactHeight++;
                }

                if (artifactHeight > maxDiameterOfArtifacts)
                    continue;
                
                // Decide if the artifact should be removed
                if (topRowWidth <= maxDiameterOfArtifacts && artifactHeight <= maxDiameterOfArtifacts)
                {
                    // Perform a flood-fill to remove the artifact
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
    
            if (cx < 0 || cx >= width || cy < 0 || cy >= width || visited[cx, cy] || heightmap[cx, cy] != targetValue)
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