#include "bfs.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstddef>
#include <omp.h>

#include "../common/CycleTimer.h"
#include "../common/graph.h"

#define ROOT_NODE_ID 0
#define NOT_VISITED_MARKER -1

void vertex_set_clear(vertex_set* list) {
    list->count = 0;
    
    // Use OpenMP to parallelize the clearing of the `present` array
    #pragma omp parallel for
    for (int i = 0; i < list->max_vertices; i++) {
        list->present[i] = false;
    }
}

void vertex_set_init(vertex_set* list, int count) {
    list->max_vertices = count;
    list->vertices = (int*)malloc(sizeof(int) * list->max_vertices);
    list->present = (bool*)malloc(sizeof(bool) * list->max_vertices);
    vertex_set_clear(list);
}

#define THREADS 8
// Take one step of "top-down" BFS.  For each vertex on the frontier,
// follow all outgoing edges, and add all neighboring vertices to the
// new_frontier.
void top_down_step(
    Graph g,
    vertex_set* frontier,
    vertex_set* new_frontier,
    int* distances)
{

    //GO THROUGH EVERY NODE IN OUR CURRENT FRONTIER
    //#pragma omp parallel for 
    // #define THREADS 4
    #pragma omp parallel for schedule(dynamic, 100)
    //for each node id in the current frontier 
    for (int i=0; i<frontier->count; i++) {
        //get the node id for each node in the frontier
        int node = frontier->vertices[i];

        int start_edge = g->outgoing_starts[node];
        //if node is the last node in the graph (node == g.num_nodes - 1). Getting indices for outgoin, dense array stuff!
        int end_edge = (node == g->num_nodes - 1)
                           ? g->num_edges // The last edge in the entire graph
                           : g->outgoing_starts[node + 1]; // Start of the next node's edges

        // attempt to add all neighbors to the new frontier
        // #pragma omp parallel for 
        for (int neighbor=start_edge; neighbor<end_edge; neighbor++) {
            int outgoing = g->outgoing_edges[neighbor];

            //IF NOT VISITED
            //SOMETHING ABOUT THIS VARIABLE DECLARATION IS SUS 
            int dist = -1; 
            #pragma omp atomic read
            dist = distances[outgoing];
            //test and test and set
            if (dist == NOT_VISITED_MARKER) { //reflect the lowest possible distance 
                if (__sync_bool_compare_and_swap(&distances[outgoing], NOT_VISITED_MARKER, distances[node] + 1)) {  
                    //if we are visiting this neigbout from the first time, do the following
                    int index = 0;
                    #pragma omp atomic capture
                    index = new_frontier->count++;
                    new_frontier->vertices[index] = outgoing;
                }
            }
        }
    }
}

// Implements top-down BFS.
//
// Result of execution is that, for each node in the graph, the
// distance to the root is stored in sol.distances.
void bfs_top_down(Graph graph, solution* sol) {

    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    vertex_set* frontier = &list1;
    vertex_set* new_frontier = &list2;

    // initialize all nodes to NOT_VISITED
    for (int i=0; i<graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER;

    // setup frontier with the root node
    frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    sol->distances[ROOT_NODE_ID] = 0;

    while (frontier->count != 0) {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(new_frontier);

        top_down_step(graph, frontier, new_frontier, sol->distances);

#ifdef VERBOSE
    double end_time = CycleTimer::currentSeconds();
    printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif

        // swap pointers
        vertex_set* tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
    }
}

void bottom_up_step(
    Graph g,
    vertex_set* frontier,
    vertex_set* new_frontier,
    int* distances)
{

    //j is the node we are checking to see if is in the fronteir
    int counter = 0;
    #pragma omp parallel for reduction(+:counter) //unifies thread local counters at the end 
    for (int node=0; node < g->num_nodes; node++) {
        //printf("hi im checking node %d\n", node);
        if (distances[node] == NOT_VISITED_MARKER) {
            //check if it has an incoming edge from a node in the fronteir 
            for (int i=0; i<frontier->count; i++) {
                //collecting all the nodes incoming edges 
                int start_edge = g->incoming_starts[node];
                int end_edge = (node == g->num_nodes - 1)
                           ? g->num_edges
                           : g->incoming_starts[node + 1];
                //going through all the nodes who are incoming to me 
                for (int neighbor=start_edge; neighbor<end_edge; neighbor++) {
                    int incoming = g->incoming_edges[neighbor];
                    if (frontier->present[incoming]) {
                        if (__sync_bool_compare_and_swap(&distances[node], NOT_VISITED_MARKER, distances[incoming] + 1)) { //my distance is the incoming edges already collected distance + 1  
                            //int index = 0;
                            new_frontier->present[node] = true;
                            // frontier->vertices[counter] = node; 
                            counter++; 
                            // #pragma omp atomic capture
                            // counter = new_frontier->count++; //add this node to the new fonteir
                            // new_frontier->vertices[index] = node;
                        }
                    }
                }
            }
            
        }
    }
    //counting how many nodes are there in the next fronteir
    new_frontier->count = counter;

    //add the actual nodes to the fronteir


}


void bfs_bottom_up(Graph graph, solution* sol)
{

    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    vertex_set* frontier = &list1;
    vertex_set* new_frontier = &list2;

    // initialize all nodes to NOT_VISITED
    for (int i=0; i<graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER;

    // setup frontier with the root node
    //frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    frontier->present[ROOT_NODE_ID] = true;
    sol->distances[ROOT_NODE_ID] = 0;

    while (frontier->count != 0) {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(new_frontier);

        bottom_up_step(graph, frontier, new_frontier, sol->distances);

#ifdef VERBOSE
    double end_time = CycleTimer::currentSeconds();
    printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif

        // swap pointers
        vertex_set* tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
    }
}

void bfs_hybrid(Graph graph, solution* sol)
{
    // CS149 students:
    //
    // You will need to implement the "hybrid" BFS here as
    // described in the handout.
}
