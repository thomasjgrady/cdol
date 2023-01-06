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

AD_Graph_Vertex* AD_Graph_Vertex_new(size_t n_data_bytes, size_t n_grad_bytes);
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
VECTOR_DECL(AD_Graph_Edge*, AD_Graph_Edge_Ptr);

typedef struct {
    AD_Graph_Vertex_Ptr_Vector* vertices;
    AD_Graph_Edge_Ptr_Vector* edges;
} AD_Graph;

AD_Graph* AD_Graph_new();
void AD_Graph_free(AD_Graph* graph);

size_t AD_Graph_add_vertex(AD_Graph* graph, AD_Graph_Vertex* v);
size_t AD_Graph_add_edge(AD_Graph* graph, AD_Graph_Edge* e);

Index_Vector* AD_Graph_sources(AD_Graph* graph);
Index_Vector* AD_Graph_sinks(AD_Graph* graph);

void AD_Graph_forward(AD_Graph* graph, Index_Vector* sources);
void AD_Graph_adjoint(AD_Graph* graph, Index_Vector* sinks);

#endif
