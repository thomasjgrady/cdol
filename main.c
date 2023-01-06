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
        *dx0 = (*dy) * (*x1);
    }
    if (dx1 != NULL) {
        *dx1 = (*dy) * (*x0);
    }
}

int main(void) {
    
    // Create a new graph
    AD_Graph* g = AD_Graph_new();
    
    // Add four vertices to the graph
    size_t x = AD_Graph_add_vertex(g, AD_Graph_Vertex_new(sizeof(double), sizeof(double)));
    size_t y = AD_Graph_add_vertex(g, AD_Graph_Vertex_new(sizeof(double), sizeof(double)));
    size_t z = AD_Graph_add_vertex(g, AD_Graph_Vertex_new(sizeof(double), sizeof(double)));
    size_t w = AD_Graph_add_vertex(g, AD_Graph_Vertex_new(sizeof(double), sizeof(double)));
    
    // Set the input values at the source vertices
    ((double*) g->vertices->data[x]->data)[0] = 2.0;
    ((double*) g->vertices->data[y]->data)[0] = 3.0;
    
    // Set up a new operator representing scalar multiplication
    Differentiable_Operation* mul_op = Differentiable_Operation_new(scalar_mul_forward, scalar_mul_adjoint);
    
    // z = x*y
    size_t e0 = AD_Graph_add_edge(g, AD_Graph_Edge_new(2, (size_t[]){x, y}, 1, (size_t[]){z}, mul_op));

    // w = z*x
    size_t e1 = AD_Graph_add_edge(g, AD_Graph_Edge_new(2, (size_t[]){z, x}, 1, (size_t[]){w}, mul_op));
    
    // Get the list of indices corresponding to source vertices and call forward
    Index_Vector* sources = AD_Graph_sources(g);
    AD_Graph_forward(g, sources);
    
    // Display values
    double x_val = ((double*) g->vertices->data[x]->data)[0];
    double y_val = ((double*) g->vertices->data[y]->data)[0];
    double z_val = ((double*) g->vertices->data[z]->data)[0];
    double w_val = ((double*) g->vertices->data[w]->data)[0];
    
    printf("==== FORWARD ====\n");
    printf("x = %1.4f\n", x_val);
    printf("y = %1.4f\n", y_val);
    printf("z = x*y = %1.4f\n", z_val);
    printf("w = z*x = %1.4f\n", w_val);
    printf("\n");
    
    // Set the value of the derivative of function J w.r.t w at 0.5
    ((double*) g->vertices->data[w]->grad)[0] = 0.5;

    // Get the list of indices corresponding to sink vertices and call adjoint
    Index_Vector* sinks = AD_Graph_sinks(g);
    AD_Graph_adjoint(g, sinks);

    double x_grad = ((double*) g->vertices->data[x]->grad)[0];
    double y_grad = ((double*) g->vertices->data[x]->grad)[0];
    double z_grad = ((double*) g->vertices->data[z]->grad)[0];
    double w_grad = ((double*) g->vertices->data[w]->grad)[0];

    printf("==== ADJOINT ====\n");
    printf("dJ/dw = %1.4f\n", w_grad);
    printf("dJ/dz = %1.4f\n", z_grad);
    printf("dJ/dy = %1.4f\n", y_grad);
    printf("dJ/dx = %1.4f\n", x_grad);

}
