package pa8;
import java.util.*;
import java.util.AbstractMap.SimpleEntry;
import java.util.ArrayList; 

/** 
 * Graph interface for a directed, connected graph with vertices numbered from 0 to n-1.
 * For weighted graphs, the default weight of each edge is assumed to be 1.
 */
interface Graph {

    /**
     * Add an edge between two vertices.
     * @param v vertex 1 (0-indexed)
     * @param w vertex 2 (0-indexed)
     * 
     * When called for weighted graphs, the default weight of the edge is assumed to be 1.
     */
    public void addEdge(int v, int w);

    /**
     *  Add a weighted edge between two vertices.
     *  @param v vertex 1 (0-indexed)
     *  @param w vertex 2 (0-indexed)
     *  @param weight the weight of the edge
     * 
     * When called for unweighted graphs, it should be equivalent to addEdge(v, w).
     * 
     */
    public void addWeightedEdge(int v, int w, int weight);

    /** 
     * Perform a Breadth-First Search (BFS) starting from the specified vertex.
     * @param start the starting vertex
     * @return a String representing the order of vertices visited, e.g., "0 1 2".
     */ 
    public String bfs(int start);

    /** 
     * Perform a Depth-First Search (DFS) starting from the specified vertex.
     * @param start the starting vertex
     * @return a String representing the order of vertices visited, e.g., "0 2 1".
     */
    public String dfs(int start);

    /**
     * Determine if there is a cycle in the graph.
     * @return true if the graph contains a cycle, false otherwise
     */ 
    public boolean hasCycle();

    /**
     * Find the shortest path between two vertices.
     * Assumes the graph is unweighted (or weights are all 1).
     * 
     * @param v source vertex
     * @param w destination vertex
     * @return a String representing the shortest path, e.g., "0 -> 1 -> 2".
     *         If no path exists, return an empty string or a message indicating no path.
     */ 
    public String shortestPath(int v, int w);
}

    class GraphMatrix implements Graph {
    protected int[][] adjMatrix;
    private int vertices;

    public GraphMatrix(int vertices) {
        this.vertices = vertices;
        this.adjMatrix = new int[vertices][vertices];
    }

    public void addEdge(int v, int w) {
        // add edge between vertices
        this.adjMatrix[v][w] = 1;
    }

    public void addWeightedEdge(int v, int w, int weight) {
        // add weighted edge
        addEdge(v, w);
        this.adjMatrix[v][w] = weight;
    }

    public String bfs(int start) {
        Queue<Integer> queue = new LinkedList<>();
        queue.add(start);
        boolean[] visited = new boolean[this.adjMatrix.length];
        visited[start] = true;
        String result = "";

        while (!queue.isEmpty()) {
            int x = queue.poll();
            result += x + " ";

            for (int i = 0; i < this.adjMatrix.length; i++) {
                if (this.adjMatrix[x][i] != 0 && !visited[i]) {
                    queue.add(i);
                    visited[i] = true;
                }
            }
        }
        return result.trim();
    }

    public String dfs(int start) {
        boolean[] visited = new boolean[this.adjMatrix.length];
        String result = "";
        result = dfsHelper(start, visited, result);
        return result.trim();
    }

    private String dfsHelper(int start, boolean[] visited, String result) {
        if (visited[start]) return result;

        visited[start] = true;
        result += start + " ";

        for (int i = 0; i < this.adjMatrix.length; i++) {
            if (this.adjMatrix[start][i] != 0 && !visited[i]) {
                result = dfsHelper(i, visited, result);
            }
        }
        return result;
    }

    public boolean hasCycle() {
        boolean[] visited = new boolean[this.adjMatrix.length];
        boolean[] recursionStack = new boolean[this.adjMatrix.length];

        for (int i = 0; i < this.adjMatrix.length; i++) {
            if (!visited[i]) {
                if (dfsCyclecheck(i, visited, recursionStack)) {
                    return true;
                }
            }
        }
        return false; 
    }

    private boolean dfsCyclecheck(int node, boolean[] visited, boolean[] recursionStack) {
        visited[node] = true;
        recursionStack[node] = true;

        for (int i = 0; i < this.adjMatrix.length; i++) {
            if (this.adjMatrix[node][i] != 0) {
                if (recursionStack[i]) {
                    return true;
                }
                if (!visited[i] && dfsCyclecheck(i, visited, recursionStack)) {
                    return true;
                }
            }
        }

        recursionStack[node] = false;
        return false;
    }

    public String shortestPath(int v, int w) {
        int[] distances = new int[vertices];
        int[] previous = new int[vertices];

        Arrays.fill(distances, Integer.MAX_VALUE);
        Arrays.fill(previous, -1);

        distances[v] = 0;

        PriorityQueue<int[]> pq = new PriorityQueue<>(Comparator.comparingInt(a -> a[1]));
        pq.add(new int[]{v, 0});

        while (!pq.isEmpty()) {
            int[] current = pq.poll();
            int currentVertex = current[0];
            int currentDist = current[1];

            if (currentVertex == w) {
                break;
            }

            for (int i = 0; i < vertices; i++) {
                if (adjMatrix[currentVertex][i] != 0) {
                    int newDist = currentDist + adjMatrix[currentVertex][i];
                    if (newDist < distances[i]) {
                        distances[i] = newDist;
                        pq.add(new int[]{i, newDist});
                        previous[i] = currentVertex;
                    }
                }
            }
        }

        String result = "";
        int curr = w;
        while (curr != -1) {
            result = curr + " " + result;
            curr = previous[curr];
        }

        return result.trim();
    }
}
    
    class GraphList implements Graph {
    private Map<Integer, ArrayList<Integer>> adjList;
    private int vertices;

    public GraphList(int vertices) {
        this.vertices = vertices;
        adjList = new HashMap<>();
        for (int i = 0; i < vertices; i++) {
            adjList.put(i, new ArrayList<>());
        }
    }

    @Override
    public void addEdge(int v, int w) {
        adjList.get(v).add(w);
        adjList.get(w).add(v);
    }

    @Override
    public void addWeightedEdge(int v, int w, int weight) {
        addEdge(v, w); 
    }

    @Override
    public String bfs(int start) {
        Queue<Integer> queue = new LinkedList<>();
        queue.add(start);

        boolean[] visited = new boolean[this.vertices];
        visited[start] = true;

        String result = "";

        while (!queue.isEmpty()) {
            int vertex = queue.poll();
            result += vertex + " ";

            for (int neighbour : adjList.get(vertex)) {
                if (!visited[neighbour]) {
                    queue.add(neighbour);
                    visited[neighbour] = true;
                }
            }
        }

        return result.trim();
    }

    @Override
    public String dfs(int start) {
        boolean[] visited = new boolean[this.vertices];
        return dfsHelper(start, visited);
    }

    public String dfsHelper(int start, boolean[] visited) {
        String result = "";
        if (visited[start]) return result;

        visited[start] = true;
        result += start + " ";  // Concatenate result

        for (int neighbour : adjList.get(start)) {
            if (!visited[neighbour]) {
                result += dfsHelper(neighbour, visited);
            }
        }
        return result;
    }

    public boolean hasCycle() {
        boolean[] visited = new boolean[this.vertices];
        boolean[] recursionStack = new boolean[this.vertices];

        for (int i = 0; i < this.vertices; i++) {
            if (!visited[i]) {
                if (dfsCyclecheck(i, visited, recursionStack, -1)) {
                    return true;
                }
            }
        }
        return false;
    }

    private boolean dfsCyclecheck(int node, boolean[] visited, boolean[] recursionStack, int parent) {
        visited[node] = true;
        recursionStack[node] = true;

        // Traverse through all the neighbors of the current node
        for (int neighbour : adjList.get(node)) {
            // If the neighbor is not visited and it's not the parent, recurse
            if (!visited[neighbour]) {
                if (dfsCyclecheck(neighbour, visited, recursionStack, node)) {
                    return true;
                }
            }
            // If the neighbor is visited and it's not the parent, a cycle is detected
            else if (recursionStack[neighbour] && neighbour != parent) {
                return true;
            }
        }

        
        recursionStack[node] = false;
        return false;
    }

    @Override
    public String shortestPath(int v, int w) {
        int[] distances = new int[vertices];
        int[] previous = new int[vertices];

        Arrays.fill(distances, Integer.MAX_VALUE);
        Arrays.fill(previous, -1);

        distances[v] = 0;

        PriorityQueue<int[]> pq = new PriorityQueue<>(Comparator.comparingInt(a -> a[1]));
        pq.add(new int[]{v, 0});

        while (!pq.isEmpty()) {
            int[] current = pq.poll();
            int currentVertex = current[0];
            int currentDist = current[1];

            if (currentVertex == w) {
                break;
            }

            for (int neighbour : adjList.get(currentVertex)) {
                int newDist = currentDist + 1;
                if (newDist < distances[neighbour]) {
                    distances[neighbour] = newDist;
                    pq.add(new int[]{neighbour, newDist});
                    previous[neighbour] = currentVertex;
                }
            }
        }

        
        String result = "";
        int curr = w;
        while (curr != -1) {
            result = curr + " " + result;
            curr = previous[curr];
        }

        return result.trim();
    }
}
    
    class GraphMatrixWeighted extends GraphMatrix {
    
        public GraphMatrixWeighted(int vertices) {
            super(vertices);
        }
        public void addWeightedEdge(int v, int w, int weight) {
            this.adjMatrix[v][w] = weight;
            this.adjMatrix[w][v] = weight;
    }

    public String bfs(int start) {
        Queue<Integer> queue = new LinkedList<>();
        boolean[] visited = new boolean[this.adjMatrix.length];
        String result = "";
    
        queue.add(start);
        visited[start] = true;
    
        while (!queue.isEmpty()) {
            int current = queue.poll();
            result += current + " ";
    
            for (int i = 0; i < this.adjMatrix.length; i++) {
                if (this.adjMatrix[current][i] != 0 && !visited[i]) {
                    queue.add(i);
                    visited[i] = true;
                }
            }
        }
    
        return result.trim();
    }
    
    public String dfs(int start) {
        boolean[] visited = new boolean[this.adjMatrix.length];
        return dfsHelper(start, visited);
    }
    
    public String dfsHelper(int current, boolean[] visited) {
            String result = "";
            visited[current] = true;
            result += current + " ";
        
            for (int i = 0; i < this.adjMatrix.length; i++) {
                if (this.adjMatrix[current][i] != 0 && !visited[i]) {
                    result += dfsHelper(i, visited);
                }
            }
        
            return result;
        }
    
    public boolean hasCycle() {
        boolean[] visited = new boolean[this.adjMatrix.length];
        boolean[] recursionStack = new boolean[this.adjMatrix.length];
    
        for (int i = 0; i < this.adjMatrix.length; i++) {
            if (!visited[i]) {
                if (dfsCycleCheck(i, visited, recursionStack)) {
                    return true;
                }
            }
        }
        return false;
    }
    
    private boolean dfsCycleCheck(int node, boolean[] visited, boolean[] recursionStack) {
        visited[node] = true;
        recursionStack[node] = true;
    
        for (int i = 0; i < this.adjMatrix.length; i++) {
            if (this.adjMatrix[node][i] != 0) {
                if (recursionStack[i]) {
                    return true;  // Cycle detected
                }
                if (!visited[i] && dfsCycleCheck(i, visited, recursionStack)) {
                    return true;
                }
            }
        }
    
        recursionStack[node] = false;
        return false;
    }
    
    public String shortestPath(int v, int w) {
        int[] distances = new int[this.adjMatrix.length];
        int[] previous = new int[this.adjMatrix.length];
        Arrays.fill(distances, Integer.MAX_VALUE);
        Arrays.fill(previous, -1);
        distances[v] = 0;
    
        PriorityQueue<int[]> pq = new PriorityQueue<>(Comparator.comparingInt(a -> a[1]));
        pq.add(new int[]{v, 0});
    
        while (!pq.isEmpty()) {
            int[] current = pq.poll();
            int currentVertex = current[0];
            int currentDist = current[1];
    
            if (currentVertex == w) {
                break;
            }
    
            for (int i = 0; i < this.adjMatrix.length; i++) {
                if (this.adjMatrix[currentVertex][i] != 0) {
                    int newDist = currentDist + this.adjMatrix[currentVertex][i];
                    if (newDist < distances[i]) {
                        distances[i] = newDist;
                        pq.add(new int[]{i, newDist});
                        previous[i] = currentVertex;
                    }
                }
            }
        }
    
        
        String result = "";
        int curr = w;
        while (curr != -1) {
            result = curr + " " + result;
            curr = previous[curr];
        }
    
        return result;
    }
    
}

    class GraphListWeighted extends GraphList {
        private Map<Integer, List<AbstractMap.SimpleEntry<Integer, Integer>>> adjListWeighted;
    
        public GraphListWeighted(int vertices) {
            super(vertices);
            adjListWeighted = new HashMap<>();
            for (int i = 0; i < vertices; i++) {
                adjListWeighted.put(i, new ArrayList<>());
            }
        }
    
        
        public void addWeightedEdge(int v, int w, int weight) {
            adjListWeighted.get(v).add(new AbstractMap.SimpleEntry<>(w, weight));
            adjListWeighted.get(w).add(new AbstractMap.SimpleEntry<>(v, weight)); // Assuming the graph is undirected
        }
    
        
        public String bfs(int start) {
            Queue<Integer> queue = new LinkedList<>();
            boolean[] visited = new boolean[adjListWeighted.size()];
            queue.add(start);
            visited[start] = true;
            String result = "";
    
            while (!queue.isEmpty()) {
                int current = queue.poll();
                result += current + " ";
    
                for (AbstractMap.SimpleEntry<Integer, Integer> neighbor : adjListWeighted.get(current)) {
                    int neighborVertex = neighbor.getKey();
                    if (!visited[neighborVertex]) {
                        visited[neighborVertex] = true;
                        queue.add(neighborVertex);
                    }
                }
            }
            return result.trim();
        }
    
        
        public String dfs(int start) {
            boolean[] visited = new boolean[adjListWeighted.size()];
            return dfsHelper(start, visited, "");
        }
    
        private String dfsHelper(int current, boolean[] visited, String result) {
            if (visited[current]) return result;
            visited[current] = true;
            result += current + " ";
    
            for (AbstractMap.SimpleEntry<Integer, Integer> neighbor : adjListWeighted.get(current)) {
                result = dfsHelper(neighbor.getKey(), visited, result);
            }
            return result;
        }
    
        public boolean hasCycle() {
            boolean[] visited = new boolean[adjListWeighted.size()];
            boolean[] recursionStack = new boolean[adjListWeighted.size()]; // To track nodes in the current recursion stack
            
            for (int i = 0; i < adjListWeighted.size(); i++) {
                if (!visited[i]) {
                    if (dfsCycleCheck(i, -1, visited, recursionStack)) {
                        return true;
                    }
                }
            }
            return false;
        }
        
        private boolean dfsCycleCheck(int current, int parent, boolean[] visited, boolean[] recursionStack) {
            visited[current] = true;
            recursionStack[current] = true; // Mark the current node as part of the recursion stack
            
            // Explore all neighbors of the current node
            for (AbstractMap.SimpleEntry<Integer, Integer> neighbor : adjListWeighted.get(current)) {
                int neighborVertex = neighbor.getKey();
                
                // If the neighbor is not visited, recurse on it
                if (!visited[neighborVertex]) {
                    if (dfsCycleCheck(neighborVertex, current, visited, recursionStack)) {
                        return true;
                    }
                }
                // If the neighbor is in the recursion stack and is not the parent, a cycle is found
                else if (recursionStack[neighborVertex] && neighborVertex != parent) {
                    return true;
                }
            }
            
            recursionStack[current] = false; // Backtrack, remove the current node from recursion stack
            return false;
        }
    
        
        public String shortestPath(int v, int w) {
            int[] distances = new int[adjListWeighted.size()];
            int[] previous = new int[adjListWeighted.size()];
            Arrays.fill(distances, Integer.MAX_VALUE);
            Arrays.fill(previous, -1);
            distances[v] = 0;
    
            PriorityQueue<AbstractMap.SimpleEntry<Integer, Integer>> pq =
                new PriorityQueue<>(Comparator.comparingInt(AbstractMap.SimpleEntry::getValue));
            pq.add(new AbstractMap.SimpleEntry<>(v, 0));
    
            while (!pq.isEmpty()) {
                AbstractMap.SimpleEntry<Integer, Integer> current = pq.poll();
                int currentVertex = current.getKey();
                int currentDist = current.getValue();
    
                if (currentVertex == w) break;
    
                for (AbstractMap.SimpleEntry<Integer, Integer> neighbor : adjListWeighted.get(currentVertex)) {
                    int neighborVertex = neighbor.getKey();
                    int weight = neighbor.getValue();
                    int newDist = currentDist + weight;
    
                    if (newDist < distances[neighborVertex]) {
                        distances[neighborVertex] = newDist;
                        previous[neighborVertex] = currentVertex;
                        pq.add(new AbstractMap.SimpleEntry<>(neighborVertex, newDist));
                    }
                }
            }
    
            if (distances[w] == Integer.MAX_VALUE) return "No path found";
    
            // Build path from end node
            String result = "";
            for (int at = w; at != -1; at = previous[at]) {
                result = at + " " + result;
            }
            return result.trim();
        }
    }

public class Main{
    public static void main(String[] args) {
        GraphMatrix graphMatrix = new GraphMatrix(5);
        graphMatrix.addEdge(0, 1);
        graphMatrix.addEdge(0, 2);
        graphMatrix.addEdge(1, 2);
        graphMatrix.addEdge(1, 3);
        graphMatrix.addEdge(2, 4);
        
        System.out.println("BFS (GraphMatrix): " + graphMatrix.bfs(0));
        System.out.println("DFS (GraphMatrix): " + graphMatrix.dfs(0));
        System.out.println("Cycle detected in GraphMatrix: " + graphMatrix.hasCycle());
        System.out.println("Shortest path (GraphMatrix): " + graphMatrix.shortestPath(0, 4));
        
        GraphList graphList = new GraphList(5);
        graphList.addEdge(0, 1);
        graphList.addEdge(0, 2);
        graphList.addEdge(1, 2);
        graphList.addEdge(1, 3);
        graphList.addEdge(2, 4);
        
        System.out.println("BFS (GraphList): " + graphList.bfs(0));
        System.out.println("DFS (GraphList): " + graphList.dfs(0));
        System.out.println("Cycle detected in GraphList: " + graphList.hasCycle());
        System.out.println("Shortest path (GraphList): " + graphList.shortestPath(0, 4));
        
        GraphMatrixWeighted graphMatrixWeighted = new GraphMatrixWeighted(5);
        graphMatrixWeighted.addWeightedEdge(0, 1, 10);
        graphMatrixWeighted.addWeightedEdge(0, 2, 5);
        graphMatrixWeighted.addWeightedEdge(1, 2, 2);
        graphMatrixWeighted.addWeightedEdge(1, 3, 1);
        graphMatrixWeighted.addWeightedEdge(2, 4, 3);
        
        System.out.println("BFS (GraphMatrixWeighted): " + graphMatrixWeighted.bfs(0));
        System.out.println("DFS (GraphMatrixWeighted): " + graphMatrixWeighted.dfs(0));
        System.out.println("Cycle detected in GraphMatrixWeighted: " + graphMatrixWeighted.hasCycle());
        System.out.println("Shortest path (GraphMatrixWeighted): " + graphMatrixWeighted.shortestPath(0, 4));
    }
}










