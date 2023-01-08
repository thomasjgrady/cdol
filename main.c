#include <stdio.h>

#include "cdol.h"

void scalar_mul_forward(AD_Graph_Vertex** input, AD_Graph_Vertex** output) {
    double* x0 = (double*) input[0]->data;
    double* x1 = (double*) input[1]->data;
    double* y  = (double*) output[0]->data;
    *y = (*x0) * (*x1);
}

void scalar_mul_adjoint(AD_Graph_Vertex** input, AD_Graph_Vertex** output) {

    double* x0  = (double*) input[0]->data;
    double* x1  = (double*) input[1]->data;
    double* y   = (double*) output[0]->data;
    double* dx0 = (double*) input[0]->grad;
    double* dx1 = (double*) input[1]->grad;
    double* dy  = (double*) output[0]->grad;

    if (dy == NULL) {
        return;
    }
    if (dx0 != NULL) {
        *dx0 += (*dy) * (*x1);
    }
    if (dx1 != NULL) {
        *dx1 += (*dy) * (*x0);
    }
}

int main(void) {
    
    // Create a new graph
    AD_Graph* g = AD_Graph_new();
    
    // Setup operators
    Differentiable_Operation* mul_op = Differentiable_Operation_new(scalar_mul_forward, scalar_mul_adjoint);

    // Setup graph vertices
    size_t x0_idx = AD_Graph_add_vertex(g, AD_Graph_Vertex_new(sizeof(double), sizeof(double), grad_buf_init_zeros));
    size_t x1_idx = AD_Graph_add_vertex(g, AD_Graph_Vertex_new(sizeof(double), sizeof(double), grad_buf_init_zeros));
    size_t x2_idx = AD_Graph_add_vertex(g, AD_Graph_Vertex_new(sizeof(double), sizeof(double), grad_buf_init_zeros));
    size_t y0_idx = AD_Graph_add_vertex(g, AD_Graph_Vertex_new(sizeof(double), sizeof(double), grad_buf_init_zeros));
    size_t y1_idx = AD_Graph_add_vertex(g, AD_Graph_Vertex_new(sizeof(double), sizeof(double), grad_buf_init_zeros));
    size_t z_idx  = AD_Graph_add_vertex(g, AD_Graph_Vertex_new(sizeof(double), sizeof(double), grad_buf_init_zeros));
    
    // Setup graph edges to compute the following sequence of operations:
    //  y0 = x0*x1
    //  y1 = x1*x2
    //  z  = y0*y1
    size_t e0_idx = AD_Graph_add_edge(g, AD_Graph_Edge_new(2, (size_t[]){x0_idx, x1_idx}, 1, (size_t[]){y0_idx}, mul_op));
    size_t e1_idx = AD_Graph_add_edge(g, AD_Graph_Edge_new(2, (size_t[]){x1_idx, x2_idx}, 1, (size_t[]){y1_idx}, mul_op));
    size_t e2_idx = AD_Graph_add_edge(g, AD_Graph_Edge_new(2, (size_t[]){y0_idx, y1_idx}, 1, (size_t[]){z_idx}, mul_op));

    // Set initial state
    ((double*) g->vertices->data[x0_idx]->data)[0] = 2.0;
    ((double*) g->vertices->data[x1_idx]->data)[0] = 3.0;
    ((double*) g->vertices->data[x2_idx]->data)[0] = 4.0;

    AD_Graph_printf(g);

    // Setup control flow
    AD_Graph_Control_Flow* control_flow = AD_Graph_Control_Flow_new();
    AD_Graph_Control_Flow_compute(g, control_flow);

    printf("sources:    "); Index_Vector_printf(control_flow->sources);
    printf("sinks:      "); Index_Vector_printf(control_flow->sinks);
    printf("edge order: "); Index_Vector_printf(control_flow->edge_order);
    
    // Run forward computation
    AD_Graph_execute(g, control_flow, FORWARD);
    printf("z = %1.4f\n", ((double*) g->vertices->data[z_idx]->data)[0]);

    // Set gradient with some data
    ((double*) g->vertices->data[z_idx]->grad)[0] = 0.01;

    // Run adjoint computation
    AD_Graph_execute(g, control_flow, ADJOINT);
    printf("dJ/dz  = %1.4f\n", ((double*) g->vertices->data[z_idx]->grad)[0]);
    printf("dJ/dx0 = %1.4f\n", ((double*) g->vertices->data[x0_idx]->grad)[0]);
    printf("dJ/dx1 = %1.4f\n", ((double*) g->vertices->data[x1_idx]->grad)[0]);
    printf("dJ/dx2 = %1.4f\n", ((double*) g->vertices->data[x2_idx]->grad)[0]);
}
