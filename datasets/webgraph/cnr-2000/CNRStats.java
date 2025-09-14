import it.unimi.dsi.webgraph.*;
import java.io.*;
import java.util.*;

public class CNRStats {
    public static void main(String[] args) throws Exception {
        System.out.println("=== CNR-2000 Graph Statistics ===");
        System.out.println("Loading graph...");
        
        ImmutableGraph graph = BVGraph.loadMapped("cnr-2000");
        
        int numNodes = graph.numNodes();
        long numArcs = graph.numArcs();
        
        System.out.println("\n--- Basic Properties ---");
        System.out.println("Nodes: " + numNodes);
        System.out.println("Arcs: " + numArcs);
        System.out.printf("Average degree: %.3f\n", (double) numArcs / numNodes);
        System.out.printf("Density: %.6f\n", (double) numArcs / (numNodes * (numNodes - 1)));
        
        // Degree analysis
        System.out.println("\n--- Degree Analysis ---");
        int maxOutdegree = 0;
        int maxIndegree = 0;
        int danglingNodes = 0;  // nodes with outdegree 0
        int sinkNodes = 0;      // nodes with indegree 0
        
        int[] indegree = new int[numNodes];
        long totalOutdegree = 0;
        long totalIndegree = 0;
        
        // Calculate outdegrees and indegrees
        for (int node = 0; node < numNodes; node++) {
            if (node % 50000 == 0) {
                System.out.println("  Processing node " + node + "/" + numNodes);
            }
            
            int outdeg = graph.outdegree(node);
            totalOutdegree += outdeg;
            if (outdeg == 0) danglingNodes++;
            if (outdeg > maxOutdegree) maxOutdegree = outdeg;
            
            // Count indegrees
            LazyIntIterator successors = graph.successors(node);
            int successor;
            while ((successor = successors.nextInt()) >= 0) {
                indegree[successor]++;
            }
        }
        
        // Analyze indegrees
        for (int i = 0; i < numNodes; i++) {
            totalIndegree += indegree[i];
            if (indegree[i] == 0) sinkNodes++;
            if (indegree[i] > maxIndegree) maxIndegree = indegree[i];
        }
        
        System.out.println("Maximum outdegree: " + maxOutdegree);
        System.out.println("Maximum indegree: " + maxIndegree);
        System.out.println("Dangling nodes (outdegree=0): " + danglingNodes + " (" + 
                         String.format("%.2f%%", 100.0 * danglingNodes / numNodes) + ")");
        System.out.println("Sink nodes (indegree=0): " + sinkNodes + " (" +
                         String.format("%.2f%%", 100.0 * sinkNodes / numNodes) + ")");
        
        // Degree distribution analysis
        System.out.println("\n--- Degree Distribution ---");
        Map<Integer, Integer> outdegreeHist = new HashMap<>();
        Map<Integer, Integer> indegreeHist = new HashMap<>();
        
        for (int node = 0; node < numNodes; node++) {
            int outdeg = graph.outdegree(node);
            int indeg = indegree[node];
            
            outdegreeHist.put(outdeg, outdegreeHist.getOrDefault(outdeg, 0) + 1);
            indegreeHist.put(indeg, indegreeHist.getOrDefault(indeg, 0) + 1);
        }
        
        // Top outdegrees
        System.out.println("Top 10 outdegree values:");
        outdegreeHist.entrySet().stream()
            .sorted((e1, e2) -> e2.getKey().compareTo(e1.getKey()))
            .limit(10)
            .forEach(e -> System.out.println("  Outdegree " + e.getKey() + ": " + e.getValue() + " nodes"));
        
        System.out.println("Top 10 indegree values:");
        indegreeHist.entrySet().stream()
            .sorted((e1, e2) -> e2.getKey().compareTo(e1.getKey()))
            .limit(10)
            .forEach(e -> System.out.println("  Indegree " + e.getKey() + ": " + e.getValue() + " nodes"));
        
        // Component analysis (simple connected components)
        System.out.println("\n--- Connectivity Analysis ---");
        boolean[] visited = new boolean[numNodes];
        int components = 0;
        int largestComponent = 0;
        
        for (int start = 0; start < numNodes; start++) {
            if (!visited[start]) {
                // BFS to find component size
                Queue<Integer> queue = new LinkedList<>();
                queue.offer(start);
                visited[start] = true;
                int componentSize = 0;
                
                while (!queue.isEmpty()) {
                    int node = queue.poll();
                    componentSize++;
                    
                    LazyIntIterator successors = graph.successors(node);
                    int successor;
                    while ((successor = successors.nextInt()) >= 0) {
                        if (!visited[successor]) {
                            visited[successor] = true;
                            queue.offer(successor);
                        }
                    }
                }
                
                components++;
                if (componentSize > largestComponent) {
                    largestComponent = componentSize;
                }
            }
        }
        
        System.out.println("Weakly connected components: " + components);
        System.out.println("Largest component size: " + largestComponent + " (" +
                         String.format("%.2f%%", 100.0 * largestComponent / numNodes) + ")");
        
        // File size and compression info
        System.out.println("\n--- Storage Information ---");
        File graphFile = new File("cnr-2000.graph");
        File offsetFile = new File("cnr-2000.offsets");
        File oblFile = new File("cnr-2000.obl");
        
        long totalSize = graphFile.length() + offsetFile.length() + oblFile.length();
        
        System.out.println("Graph file size: " + formatBytes(graphFile.length()));
        System.out.println("Offset file size: " + formatBytes(offsetFile.length()));
        System.out.println("OBL file size: " + formatBytes(oblFile.length()));
        System.out.println("Total BVGraph size: " + formatBytes(totalSize));
        System.out.printf("Bits per arc: %.3f\n", (totalSize * 8.0) / numArcs);
        System.out.printf("Bits per node: %.3f\n", (totalSize * 8.0) / numNodes);
        
        System.out.println("\n=== Analysis Complete ===");
    }
    
    private static String formatBytes(long bytes) {
        if (bytes < 1024) return bytes + " B";
        if (bytes < 1024 * 1024) return String.format("%.1f KB", bytes / 1024.0);
        if (bytes < 1024 * 1024 * 1024) return String.format("%.1f MB", bytes / (1024.0 * 1024));
        return String.format("%.1f GB", bytes / (1024.0 * 1024 * 1024));
    }
}