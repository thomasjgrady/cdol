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

    AD_Graph_printf(g);
    AD_Graph_Control_Flow* control_flow = AD_Graph_Control_Flow_new();
    AD_Graph_Control_Flow_compute(g, control_flow);

    AD_Graph_execute(g, control_flow, FORWARD);

}
