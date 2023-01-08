#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "libcdol.h"

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
VECTOR_IMPL(AD_Graph_Vertex_Ptr_Vector*, AD_Graph_Vertex_Ptr_Vector, true, AD_Graph_Vertex_Ptr_Vector_free);
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

void grad_buf_init_zeros(size_t n, void* ptr) {
    memset(ptr, 0, n);
}

AD_Graph_Vertex* AD_Graph_Vertex_new(size_t n_data_bytes, size_t n_grad_bytes, grad_buf_init_fn f) {
    AD_Graph_Vertex* v = malloc(sizeof(AD_Graph_Vertex));
    v->data = malloc(n_data_bytes);
    v->grad = n_grad_bytes > 0 ? malloc(n_grad_bytes) : NULL;
    if (v->grad != NULL && f != NULL) {
        f(n_grad_bytes, v->grad);
    }
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
    e->outgoing = Index_Vector_new(n_out);
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

void AD_Graph_printf(AD_Graph* graph) {
    printf("AD_Graph {\n");
    printf("\tNumber of vertices = %lu\n", graph->vertices->length);
    printf("\tEdges = [\n");
    AD_Graph_Edge* edge;
    for (size_t i = 0; i < graph->edges->length; i++) {
        edge = graph->edges->data[i];
        printf("\t\t[");
        for (size_t j = 0; j < edge->incoming->length; j++) {
            if (j < edge->incoming->length-1) {
                printf("%lu, ", edge->incoming->data[j]);
            }
            else {
                printf("%lu", edge->incoming->data[j]);
            }
        }
        printf("] ---{ op: %p }---> [", edge->op->forward);
        for (size_t j = 0; j < edge->outgoing->length; j++) {
            if (j < edge->outgoing->length-1) {
                printf("%lu, ", edge->outgoing->data[j]);
            }
            else {
                printf("%lu", edge->outgoing->data[j]);
            }
        }
        printf("]\n");
    }
    printf("\t]\n}\n");
}

size_t AD_Graph_add_vertex(AD_Graph* graph, AD_Graph_Vertex* v) {
    AD_Graph_Vertex_Ptr_Vector_push(graph->vertices, v);
    return graph->vertices->length-1;
}

size_t AD_Graph_add_edge(AD_Graph* graph, AD_Graph_Edge* e) {
    AD_Graph_Edge_Ptr_Vector_push(graph->edges, e);
    return graph->edges->length-1;
}

AD_Graph_Control_Flow* AD_Graph_Control_Flow_new() {
    AD_Graph_Control_Flow* control_flow = malloc(sizeof(AD_Graph_Control_Flow));
    control_flow->sources = Index_Vector_new(16);
    control_flow->sinks   = Index_Vector_new(16);
    control_flow->edge_order = Index_Vector_new(64);
    control_flow->op_inputs  = AD_Graph_Vertex_Ptr_Vector_Vector_new(64);
    control_flow->op_outputs = AD_Graph_Vertex_Ptr_Vector_Vector_new(64);
    return control_flow;
}

void AD_Graph_Control_Flow_free(AD_Graph_Control_Flow* control_flow) {

    Index_Vector_free(control_flow->sources); 
    Index_Vector_free(control_flow->sinks); 
    Index_Vector_free(control_flow->edge_order);

    // Only shallow free here because the pointers are shared with the underlying
    // graph
    for (size_t i = 0; i < control_flow->op_inputs->length; i++) {
        AD_Graph_Vertex_Ptr_Vector_free_shallow(control_flow->op_inputs->data[i]);
    }
    for (size_t i = 0; i < control_flow->op_outputs->length; i++) {
        AD_Graph_Vertex_Ptr_Vector_free_shallow(control_flow->op_outputs->data[i]);
    }
    AD_Graph_Vertex_Ptr_Vector_Vector_free_shallow(control_flow->op_inputs);
    AD_Graph_Vertex_Ptr_Vector_Vector_free_shallow(control_flow->op_outputs);
    free(control_flow);
}

/* 
 * Serializes the graph into a list of edges and their corresponding
 * input/outputs nodes. This list can then be iterated over and the
 * function at each edge applied to perform forward computation, and the
 * reverse done to perform adjoint computation.
 *
 * TODO: Cycle detection. Not all graphs are valid. Should this be up to the
 *  user to validate before calling this function?
 *
 * TODO: Control flow support. E.g. change the control_flow structure to allow
 *  for some kind of simple branching logic when it is used. This is very
 *  complicated, so I will ignore it for now :).
 *
 * TODO: Heuristics on which edges should be added in what order. Maybe some
 *  kind of BFS-type metric would be good? When introducing pipeline parallelism
 *  through graph cuts, it is likely going to be more efficient to have
 *  breadth-first execution patterns (though I have no evidence for this
 *  currently, just speculation).
 */
void AD_Graph_Control_Flow_compute(AD_Graph* graph, AD_Graph_Control_Flow* control_flow) {

    Index_Vector_clear(control_flow->sources);
    Index_Vector_clear(control_flow->sinks);
    Index_Vector_clear(control_flow->edge_order);

    // TODO: Leaks memory if this function is called on the same control_flow
    // input more than once!
    AD_Graph_Vertex_Ptr_Vector_Vector_clear(control_flow->op_inputs);
    AD_Graph_Vertex_Ptr_Vector_Vector_clear(control_flow->op_outputs);

    AD_Graph_Edge* edge;
    
    // Count the in/out degree of each vertex using connectivity information.
    // The out degree buffer is also used to keep track of which functions
    // have been serialized into the control_flow structure, thereby informing
    // whether a particular vertex should remain in the frontier set.
    size_t* vertex_degree_in  = calloc(graph->vertices->length, sizeof(size_t));
    size_t* vertex_degree_out = calloc(graph->vertices->length, sizeof(size_t));

    for (size_t i = 0; i < graph->edges->length; i++) {
        edge = graph->edges->data[i];
        for (size_t j = 0; j < edge->outgoing->length; j++) {
            vertex_degree_in[edge->outgoing->data[j]]++;
        }
        for (size_t j = 0; j < edge->incoming->length; j++) {
            vertex_degree_out[edge->incoming->data[j]]++;
        }
    }
    
    // Set up the sources/sinks in the control_flow structure to be vertices in
    // the graph whose in/out degrees are zero respectively.
    for (size_t i = 0; i < graph->vertices->length; i++) {
        if (vertex_degree_in[i] == 0) {
            Index_Vector_push(control_flow->sources, i);
        }
        if (vertex_degree_out[i] == 0) {
            Index_Vector_push(control_flow->sinks, i);
        }
    }
    
    // Initialize the frontier set with all sources.
    // TODO: Change this from an O(|V|) array to something else. In large graphs
    // this becomes memory intensive.
    bool* in_frontier = calloc(graph->vertices->length, sizeof(size_t));
    for (size_t i = 0; i < control_flow->sources->length; i++) {
        in_frontier[control_flow->sources->data[i]] = true;
    }

    Index_Vector* selected_edges = Index_Vector_new(16);
    
    // Iterate until we have selected all edges
    while (control_flow->edge_order->length < graph->edges->length) {
        
        // Select all edges whose inputs are entirely within the frontier set.
        Index_Vector_clear(selected_edges);
        for (size_t i = 0; i < graph->edges->length; i++) {
            edge = graph->edges->data[i];
            bool all_in_frontier = true;
            for (size_t i = 0; i < edge->incoming->length; i++) {
                all_in_frontier &= in_frontier[edge->incoming->data[i]];
            }

            if (all_in_frontier) {
                Index_Vector_push(selected_edges, i);
            }
        }
        
        for (size_t i = 0; i < selected_edges->length; i++) {

            Index_Vector_push(control_flow->edge_order, selected_edges->data[i]);
            edge = graph->edges->data[selected_edges->data[i]];
            
            // Get the input arguments to the function associated with the
            // current edge, and decrement to the function call counter corresponding
            // to each argument by 1.
            AD_Graph_Vertex_Ptr_Vector* inputs = AD_Graph_Vertex_Ptr_Vector_new(edge->incoming->length);
            for (size_t j = 0; j < edge->incoming->length; j++) {
                size_t inc_idx = edge->incoming->data[j];
                vertex_degree_out[inc_idx]--;
                AD_Graph_Vertex_Ptr_Vector_push(inputs, graph->vertices->data[inc_idx]);
                if (vertex_degree_out[inc_idx] == 0) {
                    in_frontier[inc_idx] = false;
                }
            }
            
            // Get the output arguments to the function associated with the current edge
            AD_Graph_Vertex_Ptr_Vector* outputs = AD_Graph_Vertex_Ptr_Vector_new(edge->outgoing->length);
            for (size_t j = 0; j < edge->outgoing->length; j++) {
                in_frontier[edge->outgoing->data[j]] = true;
                AD_Graph_Vertex_Ptr_Vector_push(outputs, graph->vertices->data[edge->outgoing->data[j]]);
            }

            AD_Graph_Vertex_Ptr_Vector_Vector_push(control_flow->op_inputs, inputs);
            AD_Graph_Vertex_Ptr_Vector_Vector_push(control_flow->op_outputs, outputs);
        }
    }

    Index_Vector_free(selected_edges);
    free(in_frontier);
    free(vertex_degree_out);
    free(vertex_degree_in);
}

void AD_Graph_execute(AD_Graph* graph, AD_Graph_Control_Flow* control_flow, AD_Graph_Execution_Mode mode) {
    if (mode == FORWARD) {
        AD_Graph_Edge* edge;
        AD_Graph_Vertex_Ptr_Vector* input;
        AD_Graph_Vertex_Ptr_Vector* output;
        for (size_t i = 0; i < control_flow->edge_order->length; i++) {
            edge = graph->edges->data[control_flow->edge_order->data[i]];
            input = control_flow->op_inputs->data[i];
            output = control_flow->op_outputs->data[i];
            edge->op->forward(input->data, output->data);
        }
    }
    else {
        AD_Graph_Edge* edge;
        AD_Graph_Vertex_Ptr_Vector* input;
        AD_Graph_Vertex_Ptr_Vector* output;
        for (size_t i = control_flow->edge_order->length; i-- > 0;) {
            edge = graph->edges->data[control_flow->edge_order->data[i]];
            input = control_flow->op_inputs->data[i];
            output = control_flow->op_outputs->data[i];
            edge->op->adjoint(input->data, output->data);
        }
    }
}
