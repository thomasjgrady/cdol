#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "cdol.h"

#define VECTOR_IMPL(T, name_prefix, free_elements, free_func) \
    name_prefix##_Vector* name_prefix##_Vector_new(size_t default_capacity) { \
        name_prefix##_Vector* v = malloc(sizeof(name_prefix##_Vector)); \
        v->length = 0; \
        v->capacity = default_capacity; \
        v->data = malloc(default_capacity*sizeof(T)); \
        return v; \
    } \
    name_prefix##_Vector* name_prefix##_Vector_copy(name_prefix##_Vector* v) { \
        name_prefix##_Vector* u = name_prefix##_Vector_new(v->capacity); \
        u->length = v->length; \
        memcpy(u->data, v->data, u->length*sizeof(T)); \
        return u; \
    } \
    void name_prefix##_Vector_free(name_prefix##_Vector* v) { \
        if (free_elements) { \
            for (size_t i = 0; i < v->length; i++) { \
                free_func(v->data[i]); \
            } \
        } \
        free(v->data); \
        free(v); \
    } \
    void name_prefix##_Vector_free_shallow(name_prefix##_Vector* v) { \
        free(v->data); \
        free(v); \
    } \
    void name_prefix##_Vector_push(name_prefix##_Vector* v, T item) { \
        if (v->length >= v->capacity) { \
            v->capacity = v->capacity == 0 ? 1 : 2*v->capacity; \
            v->data = realloc(v->data, v->capacity*sizeof(T)); \
        } \
        v->data[v->length] = item; \
        v->length++; \
    } \
    void name_prefix##_Vector_clear(name_prefix##_Vector* v) { \
        v->length = 0; \
    }

VECTOR_IMPL(AD_Graph_Vertex*, AD_Graph_Vertex_Ptr, true, AD_Graph_Vertex_free);
VECTOR_IMPL(AD_Graph_Edge*, AD_Graph_Edge_Ptr, true, AD_Graph_Edge_free);
VECTOR_IMPL(size_t, Index, false, );

void Index_Vector_printf(Index_Vector* v) {
    printf("Index Vector: {\n");
    printf("\tlength   = %lu\n", v->length);
    printf("\tcapacity = %lu\n", v->capacity);
    printf("\tdata     = [");
    for (size_t i = 0; i < v->length; i++) {
        if (i < v->length-1) {
            printf("%lu, ", v->data[i]);
        }
        else {
            printf("%lu", v->data[i]);
        }
    }
    printf("]\n}\n");
}

AD_Graph_Vertex* AD_Graph_Vertex_new(size_t n_data_bytes, size_t n_grad_bytes) {
    AD_Graph_Vertex* v = malloc(sizeof(AD_Graph_Vertex));
    v->data = malloc(n_data_bytes);
    v->grad = n_grad_bytes > 0 ? malloc(n_grad_bytes) : NULL;
    return v;
}

void AD_Graph_Vertex_free(AD_Graph_Vertex* v) {
    free(v->data);
    if (v->grad != NULL) {
        free(v->grad);
    }
    free(v);
}

Differentiable_Operation* Differentiable_Operation_new(AD_Function forward, AD_Function adjoint) {
    Differentiable_Operation* op = malloc(sizeof(Differentiable_Operation));
    op->forward = forward;
    op->adjoint = adjoint;
    return op;
}

AD_Graph_Edge* AD_Graph_Edge_new(size_t n_in, size_t* in, size_t n_out, size_t* out, Differentiable_Operation* op) {
    AD_Graph_Edge* e = malloc(sizeof(AD_Graph_Edge));
    e->incoming = Index_Vector_new(n_in);
    e->outgoing = Index_Vector_new(n_in);
    memcpy(e->incoming->data, in, n_in*sizeof(size_t));
    memcpy(e->outgoing->data, out, n_out*sizeof(size_t));
    e->incoming->length = n_in;
    e->outgoing->length = n_out;
    e->op = op;
    return e;
}

void AD_Graph_Edge_free(AD_Graph_Edge* e) {
    Index_Vector_free(e->incoming);
    Index_Vector_free(e->outgoing);
    free(e->op);
    free(e);
}

AD_Graph* AD_Graph_new() {
    AD_Graph* graph = malloc(sizeof(AD_Graph));
    graph->vertices = AD_Graph_Vertex_Ptr_Vector_new(64);
    graph->edges = AD_Graph_Edge_Ptr_Vector_new(64);
    return graph;
}

void AD_Graph_free(AD_Graph* graph) {
    AD_Graph_Vertex_Ptr_Vector_free(graph->vertices);
    AD_Graph_Edge_Ptr_Vector_free(graph->edges);
    free(graph);
}

size_t AD_Graph_add_vertex(AD_Graph* graph, AD_Graph_Vertex* v) {
    AD_Graph_Vertex_Ptr_Vector_push(graph->vertices, v);
    return graph->vertices->length-1;
}

size_t AD_Graph_add_edge(AD_Graph* graph, AD_Graph_Edge* e) {
    AD_Graph_Edge_Ptr_Vector_push(graph->edges, e);
    return graph->edges->length-1;
}

Index_Vector* AD_Graph_sources(AD_Graph* graph) {

    // Create a temporary array to store incoming counts for each vertex
    size_t* incoming_counts = calloc(graph->vertices->length, sizeof(size_t));

    // Loop over each edge and and add 1 to each vertex in the outgoing section
    AD_Graph_Edge* e;
    for (size_t i = 0; i < graph->edges->length; i++) {
        e = graph->edges->data[i];
        for (size_t j = 0; j < e->outgoing->length; j++) {
            incoming_counts[e->outgoing->data[j]]++;
        }
    }

    // Collect pointers to all vertices with zero incoming edges
    Index_Vector* v = Index_Vector_new(16);
    for (size_t i = 0; i < graph->vertices->length; i++) {
        if (incoming_counts[i] == 0) {
            Index_Vector_push(v, i);
        }
    }

    // Free temporary array
    free(incoming_counts);

    return v;
}

Index_Vector* AD_Graph_sinks(AD_Graph* graph) {

    // Create a temporary array to store incoming counts for each vertex
    size_t* outgoing_counts = calloc(graph->vertices->length, sizeof(size_t));

    // Loop over each edge and and add 1 to each vertex in the incoming section
    AD_Graph_Edge* e;
    for (size_t i = 0; i < graph->edges->length; i++) {
        e = graph->edges->data[i];
        for (size_t j = 0; j < e->incoming->length; j++) {
            outgoing_counts[e->incoming->data[j]]++;
        }
    }

    // Collect pointers to all vertices with zero outgoing edges
    Index_Vector* v = Index_Vector_new(16);
    for (size_t i = 0; i < graph->vertices->length; i++) {
        if (outgoing_counts[i] == 0) {
            Index_Vector_push(v, i);
        }
    }

    // Free temporary array
    free(outgoing_counts);

    return v;
}

void AD_Graph_forward(AD_Graph* graph, Index_Vector* sources) {

    // TODO: move most of this to a separate function to allow computation of
    // control flow to be amortized as an upfront cost
    
    // Local variables for convenience
    AD_Graph_Edge* edge;

    // Allocate a table which keeps track of the number of outgoing edges from
    // each vertex who's associated operation has been called. The count for
    // each vertex should reach zero when computation is finished
    size_t* outgoing_counts = calloc(graph->vertices->length, sizeof(size_t));
    for (size_t i = 0; i < graph->edges->length; i++) {
        edge = graph->edges->data[i];
        for (size_t j = 0; j < edge->incoming->length; j++) {
            outgoing_counts[edge->incoming->data[j]]++;
        }
    }

    // Allocate a table of booleans to indicate whether a particular vertex
    // is in the frontier set
    bool* in_frontier = calloc(graph->vertices->length, sizeof(size_t));
    for (size_t i = 0; i < sources->length; i++) {
        in_frontier[sources->data[i]] = true;
    }
    
    // Allocate vector to hold currently selected edges
    Index_Vector* selected_edges = Index_Vector_new(16);
    
    // Allocate vectors to hold input/output vertices
    AD_Graph_Vertex_Ptr_Vector* input_vertices  = AD_Graph_Vertex_Ptr_Vector_new(16);
    AD_Graph_Vertex_Ptr_Vector* output_vertices = AD_Graph_Vertex_Ptr_Vector_new(16);
    
    // Run computation
    while (true) {

        // Select edges for which all of their incoming vertices are in the
        // frontier set
        Index_Vector_clear(selected_edges);
        for (size_t i = 0; i < graph->edges->length; i++) {
            edge = graph->edges->data[i];
            bool all_in_frontier = true;
            for (size_t j = 0; j < edge->incoming->length; j++) {
                all_in_frontier &= in_frontier[edge->incoming->data[j]];
            }

            if (all_in_frontier) {
                Index_Vector_push(selected_edges, i);
            }
        }

        // If there were no selected edges, computation is finished
        if (selected_edges->length == 0) {
            break;
        }

        // For each of the selected edges, assemble a vector of input/output
        // vertices and run the computation. After computation, decrement the
        // counter corresponding to the incoming vertices
        for (size_t i = 0; i < selected_edges->length; i++) {
            
            edge = graph->edges->data[selected_edges->data[i]];

            AD_Graph_Vertex_Ptr_Vector_clear(input_vertices);
            for (size_t j = 0; j < edge->incoming->length; j++) {
                AD_Graph_Vertex_Ptr_Vector_push(
                        input_vertices, 
                        graph->vertices->data[edge->incoming->data[j]]
                );
            }

            AD_Graph_Vertex_Ptr_Vector_clear(output_vertices);
            for (size_t j = 0; j < edge->outgoing->length; j++) {
                AD_Graph_Vertex_Ptr_Vector_push(
                        output_vertices, 
                        graph->vertices->data[edge->outgoing->data[j]]
                );
            }

            edge->op->forward(input_vertices->data, output_vertices->data);

            for (size_t j = 0; j < edge->incoming->length; j++) {
                outgoing_counts[edge->incoming->data[j]]--;
            }
        }
        
        // Update the froniter array using selected edges and outgoing counts
        for (size_t i = 0; i < selected_edges->length; i++) {
            edge = graph->edges->data[selected_edges->data[i]];

            for (size_t j = 0; j < edge->incoming->length; j++) {
                if (outgoing_counts[edge->incoming->data[j]] == 0) {
                    in_frontier[edge->incoming->data[j]] = false;
                }
            }

            for (size_t j = 0; j < edge->outgoing->length; j++) {
                in_frontier[edge->outgoing->data[j]] = true;
            }
        }
    }
    
    // Shallow free here because the AD_Graph_Vertex pointers are shared by the
    // underlying graph itself
    AD_Graph_Vertex_Ptr_Vector_free_shallow(output_vertices);
    AD_Graph_Vertex_Ptr_Vector_free_shallow(input_vertices);

    Index_Vector_free(selected_edges);

    free(in_frontier);
    free(outgoing_counts);

}

void AD_Graph_adjoint(AD_Graph* graph, Index_Vector* sinks) {

    // TODO: move most of this to a separate function to allow computation of
    // control flow to be amortized as an upfront cost
    
    // Local variables for convenience
    AD_Graph_Edge* edge;

    // Allocate a table which keeps track of the number of incoming edges from
    // each vertex who's associated adjoint operation has been called. The count
    // for each vertex should reach zero when computation is finished
    size_t* incoming_counts = calloc(graph->vertices->length, sizeof(size_t));
    for (size_t i = 0; i < graph->edges->length; i++) {
        edge = graph->edges->data[i];
        for (size_t j = 0; j < edge->outgoing->length; j++) {
            incoming_counts[edge->outgoing->data[j]]++;
        }
    }

    // Allocate a table of booleans to indicate whether a particular vertex
    // is in the frontier set
    bool* in_frontier = calloc(graph->vertices->length, sizeof(size_t));
    for (size_t i = 0; i < sinks->length; i++) {
        in_frontier[sinks->data[i]] = true;
    }
    
    // Allocate vector to hold currently selected edges
    Index_Vector* selected_edges = Index_Vector_new(16);
    
    // Allocate vectors to hold input/output vertices
    AD_Graph_Vertex_Ptr_Vector* input_vertices  = AD_Graph_Vertex_Ptr_Vector_new(16);
    AD_Graph_Vertex_Ptr_Vector* output_vertices = AD_Graph_Vertex_Ptr_Vector_new(16);
    
    // Run computation
    while (true) {

        // Select edges for which all of their incoming vertices are in the
        // frontier set
        Index_Vector_clear(selected_edges);
        for (size_t i = 0; i < graph->edges->length; i++) {
            edge = graph->edges->data[i];
            bool all_in_frontier = true;
            for (size_t j = 0; j < edge->outgoing->length; j++) {
                all_in_frontier &= in_frontier[edge->outgoing->data[j]];
            }

            if (all_in_frontier) {
                Index_Vector_push(selected_edges, i);
            }
        }

        // If there were no selected edges, computation is finished
        if (selected_edges->length == 0) {
            break;
        }

        // For each of the selected edges, assemble a vector of input/output
        // vertices and run the computation. After computation, decrement the
        // counter corresponding to the outgoing vertices
        for (size_t i = 0; i < selected_edges->length; i++) {
            
            edge = graph->edges->data[selected_edges->data[i]];

            AD_Graph_Vertex_Ptr_Vector_clear(input_vertices);
            for (size_t j = 0; j < edge->incoming->length; j++) {
                AD_Graph_Vertex_Ptr_Vector_push(
                        input_vertices, 
                        graph->vertices->data[edge->incoming->data[j]]
                );
            }

            AD_Graph_Vertex_Ptr_Vector_clear(output_vertices);
            for (size_t j = 0; j < edge->outgoing->length; j++) {
                AD_Graph_Vertex_Ptr_Vector_push(
                        output_vertices, 
                        graph->vertices->data[edge->outgoing->data[j]]
                );
            }

            edge->op->adjoint(input_vertices->data, output_vertices->data);

            for (size_t j = 0; j < edge->outgoing->length; j++) {
                incoming_counts[edge->outgoing->data[j]]--;
            }
        }
        
        // Update the froniter array using selected edges and incoming counts
        for (size_t i = 0; i < selected_edges->length; i++) {
            edge = graph->edges->data[selected_edges->data[i]];

            for (size_t j = 0; j < edge->outgoing->length; j++) {
                if (incoming_counts[edge->outgoing->data[j]] == 0) {
                    in_frontier[edge->outgoing->data[j]] = false;
                }
            }

            for (size_t j = 0; j < edge->incoming->length; j++) {
                in_frontier[edge->incoming->data[j]] = true;
            }
        }
    }
    
    // Shallow free here because the AD_Graph_Vertex pointers are shared by the
    // underlying graph itself
    AD_Graph_Vertex_Ptr_Vector_free_shallow(output_vertices);
    AD_Graph_Vertex_Ptr_Vector_free_shallow(input_vertices);

    Index_Vector_free(selected_edges);

    free(in_frontier);
    free(incoming_counts);

}
