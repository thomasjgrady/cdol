#ifndef CDOL_H
#define CDOL_H

#include <stdlib.h>

#define VECTOR_DECL(T, name_prefix) \
    typedef struct { \
        size_t length; \
        size_t capacity; \
        T* data; \
    } name_prefix##_Vector; \
    name_prefix##_Vector* name_prefix##_Vector_new(size_t default_capacity); \
    name_prefix##_Vector* name_prefix##_Vector_copy(name_prefix##_Vector* v); \
    void name_prefix##_Vector_free(name_prefix##_Vector* v); \
    void name_prefix##_Vector_free_shallow(name_prefix##_Vector* v); \
    void name_prefix##_Vector_push(name_prefix##_Vector* v, T item); \
    void name_prefix##_Vector_clear(name_prefix##_Vector* v)

VECTOR_DECL(size_t, Index);
void Index_Vector_printf(Index_Vector* v);

typedef struct {
    void* data;
    void* grad;
} AD_Graph_Vertex;

typedef void (*grad_buf_init_fn)(size_t, void*); // TODO: better interface for this
void grad_buf_init_zeros(size_t n, void* ptr);

AD_Graph_Vertex* AD_Graph_Vertex_new(size_t n_data_bytes, size_t n_grad_bytes, grad_buf_init_fn f);
void AD_Graph_Vertex_free(AD_Graph_Vertex* v);

typedef void (*AD_Function)(AD_Graph_Vertex**, AD_Graph_Vertex**);

typedef struct {
    AD_Function forward;
    AD_Function adjoint;
} Differentiable_Operation;

Differentiable_Operation* Differentiable_Operation_new(AD_Function forward, AD_Function adjoint);

typedef struct {
    Index_Vector* incoming;
    Index_Vector* outgoing;
    Differentiable_Operation* op;
} AD_Graph_Edge;

AD_Graph_Edge* AD_Graph_Edge_new(size_t n_in, size_t* in, size_t n_out, size_t* out, Differentiable_Operation* op);
void AD_Graph_Edge_free(AD_Graph_Edge* e);

VECTOR_DECL(AD_Graph_Vertex*, AD_Graph_Vertex_Ptr);
VECTOR_DECL(AD_Graph_Vertex_Ptr_Vector*, AD_Graph_Vertex_Ptr_Vector);
VECTOR_DECL(AD_Graph_Edge*, AD_Graph_Edge_Ptr);

typedef struct {
    AD_Graph_Vertex_Ptr_Vector* vertices;
    AD_Graph_Edge_Ptr_Vector* edges;
} AD_Graph;

AD_Graph* AD_Graph_new();
void AD_Graph_free(AD_Graph* graph);
void AD_Graph_printf(AD_Graph* graph);

size_t AD_Graph_add_vertex(AD_Graph* graph, AD_Graph_Vertex* v);
size_t AD_Graph_add_edge(AD_Graph* graph, AD_Graph_Edge* e);

typedef struct {
    Index_Vector* sources;    // Vector of source nodes in the graph (used to input data)
    Index_Vector* sinks;      // Vector of sink nodes in the graph   (used to read out data)
    Index_Vector* edge_order; // Order in which the operations on edges are applied
    AD_Graph_Vertex_Ptr_Vector_Vector* op_inputs;  // Input arguments to the edges in order
    AD_Graph_Vertex_Ptr_Vector_Vector* op_outputs; // Output arguments to the edges in order
} AD_Graph_Control_Flow;

AD_Graph_Control_Flow* AD_Graph_Control_Flow_new();
void AD_Graph_Control_Flow_free(AD_Graph_Control_Flow* control_flow);
void AD_Graph_Control_Flow_compute(AD_Graph* graph, AD_Graph_Control_Flow* control_flow);

typedef enum {
    FORWARD,
    ADJOINT
} AD_Graph_Execution_Mode;

void AD_Graph_execute(AD_Graph* graph, AD_Graph_Control_Flow* control_flow, AD_Graph_Execution_Mode mode);

#endif
