//Inspired from A Faster Algorithm for Betweenness Centrality paper.

//Takes O(nm) and O(nm + n^2log n) time on unweighted and weighted networks.

#include <iostream>
#include <cstdio>
#include <fstream>
#include <sys/time.h>
#include <vector>
#include <queue>
#include <stack>
#include <map>
#include <set>
#include <list>

#include <cstring>
#include <sstream>
#include <algorithm>
#include <limits>

using namespace std;

double clkbegin, clkend, t;

double rtclock(void)
{
  struct timezone Tzp;
  struct timeval Tp;
  int stat;
  stat = gettimeofday (&Tp, &Tzp);
  if (stat != 0) printf("Error return from gettimeofday: %d",stat);
  return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}

// The new adjacency list type.
typedef vector<vector<int> > adjacency_list;

// BFS algorithm is used to calculate all the single source shortest paths in a non weighted graph and the source's closeness.
float bfs_SSSP(int src, int n, stack<int> &visitStack, vector<int> &sigma, list<int> *pred, const adjacency_list &adjList) {
    // Closeness counter.
    float closeness = 0;
    
    // Vector that holds the distances from the source.
    vector<int> dist;
    dist.resize(n, -1);
    dist[src] = 0;

    // Queue used for the Bfs algorithm.
    queue<int> visitQueue;
    visitQueue.push(src);
    //cout<<"debug 30"<<endl;
    // While there are still elements in the queue.
    while (!visitQueue.empty()) {
        // Pop the first.
        int v = visitQueue.front();
        visitQueue.pop();
        visitStack.push(v);
        
        // Closeness part aggregation.
        closeness += dist[v];

        // Check the neighbors w of v.
        //for (vector<int>::iterator it = adjList[v].begin(); it != adjList[v].end(); it++) {
        cout<<"adj size : "<<adjList[v].size()<<"; this is v :"<<v<<endl;
        for(int i=0;i<adjList[v].size();i++){    
            int w = adjList[v].at(i);
            // Node w found for the first time?
            if(v==62)cout<<"debug x"<<v<<" "<<w<<endl;
            if (dist[w] < 0) {
                visitQueue.push(w);
                dist[w] = dist[v] + 1;
            }

            // Is the shortest path to w via u?
            if (dist[w] == dist[v] + 1) {
                pred[w].push_back(v);
                sigma[w] += sigma[v];
            }

        }

    }    
    // Closeness part inversion.
    if (closeness!=0) {
        return 1.0 / closeness;
    } else {
        return 0;
    }
}

// Prompts the user to type an input filename and returns the file pointer.
FILE * readPrompt() {
    FILE * fp;
    char str[100];
    
    // Prompt filename input.
    cout << "Enter filename: ";
    cin >> str;
    
    // Append .net if the user did not provide it.
    if (strstr(str, ".net") == NULL) {
        strcat(str, ".net");
    }
    
    // Open the file and return fp, if it exists, otherwise exit with error.
    fp = fopen(str, "r");
    if (fp == NULL) {
        cout << "File '" << str << "' not found.";
        exit(EXIT_FAILURE);
    } else {
        return fp;
    }
}

// Prints input file statistics just after input has finished.
void printInputStats(bool isWeigthed, int n, int e) {
    ofstream out;
    out.open ("out_graph_stats.txt");
    cout << "\n==================="
            << "\nINPUT GRAPH STATS"
            << "\n>Weighted: " << boolalpha << bool(isWeigthed)
            << "\n>#ofNodes: " << n
            << "\n>#ofEdges: " << e
            << "\n===================\n\n";
    out << "Weighted: " << boolalpha << bool(isWeigthed)
            << "\n>#ofNodes: " << n
            << "\n>#ofEdges: " << e;
    out.close();
}

// Reads an input file and fills up the adjacency list as well as the edges.
void readGraph(int &n, bool &isWeigthed, adjacency_list &adjList) {

    int e = 0; // Total number of edges (for statistics).
    isWeigthed = false;

    char * line = NULL;
    size_t len = 0;
    FILE * fp = readPrompt();

    // Find n, the total number of nodes.
    if (getline(&line, &len, fp) != -1) {
    	        strtok(line, " ");
	        n = atoi(strtok(NULL, " "));

    }
    cout<<"read number of nodes :"<<n<<endl;
    // Reserve n space for adjacency list. If it fails, n was not parsed.
    if (n) {
        adjList.reserve(n);
    } else {
        cout << "Malformed input. Number of nodes undefined.";
        exit(EXIT_FAILURE);
    }
    cout<<"debug 1\n";

    // Read the nodes and the edges, one by one, and fill up adjList and edgeBetweenness.
    int start, end, weight;
    while (getline(&line, &len, fp) != -1) {
        e += 1;
        start = atoi(strtok(line, " "));//-1;
        end = atoi(strtok(NULL, " "));//-1;

        adjList[start].push_back(end);
        adjList[end].push_back(start);

    }

    if (line) {
        free(line);
    }
    
    // Print statistics after reading.
    printInputStats(isWeigthed, n, e);
}

// Clears the variables or re-initializes to 0, so that they are ready for the next loop.
void resetVariables(int src, int n, list<int> *pred, vector<int> &sigma, vector<float> &delta) {
    for (int i = 0; i < n; i++) {
        pred[i].clear();
    }

    sigma.clear();
    sigma.resize(n, 0);
    sigma[src] = 1;

    delta.clear();
    delta.resize(n, 0);
}

// Prints Closeness Centrality.
void printCloseness( int n, vector<float> closeness, bool normalize) {
    float nrml = 1;
    if (normalize) {
        nrml = 1.0/(n - 1);
    } 
    ofstream out;
    out.open ("out_closeness.txt");
    cout << "> Closeness Centrality" << endl;    
    for (int i = 0; i < n; i++) {
        cout << "Node " << i << ": " << closeness[i] / nrml << endl;
        out << "Node " << i << ": " << closeness[i] / nrml << endl;
    }
    out.close();
}

// Prints Node Betweenness Centrality.
void printNodeBetweenness( int n, vector<float> nodeBetweenness, bool normalize) {
    float nrml = 1;
    if (normalize) {
        nrml = (n - 1)*(n - 2);
    }
    ofstream out;
    out.open ("out_node_betweenness.txt");
    cout << endl << "> Node Betweenness Centrality" << endl;
    for (int i = 0; i < n; i++) {
        cout << "Node " << i << ": " << nodeBetweenness[i] / nrml << endl;
        out << "Node " << i << ": " << nodeBetweenness[i] / nrml << endl;
    }
    out.close();
}

// Prints Edge Betweenness Centrality.
void printEdgeBetweenness( int n, map<string, float> edgeBetweenness, bool normalize) {
    float nrml = 1;
    if (normalize) {
        nrml = n * (n - 1);
    }
    ofstream out;
    out.open ("out_edge_tweenness.txt");
    cout << endl << "> Edge Betweenness Centrality" << endl;
    for (map<string, float>::iterator it = edgeBetweenness.begin(); it != edgeBetweenness.end(); it++) {
        cout << "Edge " << it->first << ": " << it->second / nrml << endl;
        out << "Edge " << it->first << ": " << it->second / nrml << endl;
    }
    out.close();
}

void printadjlist(int i, adjacency_list& adjList){
    cout<<"list for node "<<i<<endl;
    for(int j =0;j<adjList[i].size();j++){
        //adjList[i].at(j) = adjList[i].at(j)+1;
        cout<<adjList[i].at(j)<<" ";
    }
    cout<<endl;
}

int main(void) {
	int n; // Number of nodes
    bool isWeigthed; // Weighted graph check.
    adjacency_list adjList; // Adjacency list.

    // Centrality measures.
    map<string, float> edgeBetweenness;
    vector<float> nodeBetweenness;
    nodeBetweenness.resize(n, 0);
    vector<float> closeness;
    closeness.resize(n, 0);

    // Input is read, and values are set to all the arguments.
    readGraph(n, isWeigthed, adjList);
    //copyvector(adjList,adjList2,n);
    list<int> pred[n]; // List of predecessors of node v.
    cout<<"debug 20"<<endl;
    vector<int> sigma;
    cout<<"debug 21"<<endl;
    vector<float> delta;
    stack<int> visitStack; // Stack that holds the inverse order of visited nodes.
    clkbegin = rtclock();
    
    // For each node of the graph.
    
    for (int src = 0; src < n; src++) { 
        // Prepare the variables for the next loop.
        //adjacency_list adjList; // Adjacency list.
        //readGraph(n, isWeigthed, adjList, edgeBetweenness);
        resetVariables(src, n, pred, sigma, delta);
        cout<<"debug 22.1 : "<<src<<endl;
        closeness[src] = bfs_SSSP(src, n, visitStack, sigma, pred, adjList);

        while (!visitStack.empty()) {
            int w = visitStack.top();
            visitStack.pop();
            
            // For each predecessors of node w, do the math!
            for (list<int>::iterator it = pred[w].begin(); it != pred[w].end(); it++) {
                int v = *it;
                float c = ((float) sigma[v] / (float) sigma[w]) * (1.0 + delta[w]);

                delta[v] += c;

            }
            // Node betweenness aggregation part.
            if (w != src) {
                nodeBetweenness[w] += delta[w];
            }
        }
    }

    clkend = rtclock();
    t = clkend-clkbegin;
    cout << "\n" ;
    cout << "Time Taken : " << t;
    cout <<" Approx GFlops : " << 1.0*n*n/t/1e9 << "\n";
    cout << "\n";

    // Printing output.
    printCloseness(n, closeness, true);
    printNodeBetweenness(n, nodeBetweenness, true);

    return 0;
}