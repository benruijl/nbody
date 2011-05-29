#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int N; // dimension
double h;

int n = 2; // n bodies
double x_bodies[2][2];
double v_bodies[2][2];
double m_bodies[2];
int current_body;
double G = 6.674e-11;

// gravity test: G is 1 
void gravity(double t, double* x, double* dx, double* result) {
    int i, j;
    double dist;

    for (i = 0; i < N; i++) {
        result[i] = 0;
    }

    for (j = 0; j < n; j++) {
        if (j != current_body) {
            dist = 0; 
            for (i = 0; i < N; i++) {
                dist += (x[i] - x_bodies[j][i]) * (x[i] - x_bodies[j][i]);
            }

            dist = sqrt(dist);
            dist = 1.0 / dist * dist * dist;
            
            for (i = 0; i < N; i++) {
                result[i] += G * m_bodies[j] * (x_bodies[j][i] - x[i]) * dist;
            }
        }
    }
}

// test function: x = t^2
void test(double t, double* x, double* dx, double* result) {
    int i;

    for (i = 0; i < N; i++) {
        result[i] = 1;
    }
}

void runge_kutta(void f(double t, double* x, double* dx, double* result), double t, double* x, double* dx, double* new_x, double* new_dx) {
    double* ap[9]; // array pool
    int i;

    for (i = 0; i < 9; i++) {
        ap[i] = (double*) malloc( N * sizeof(double));
    }

    f(t, x, dx, ap[0]); 
    for (i = 0; i < N; i++) {
        ap[0][i] = 0.5 * h * h * ap[0][i]; 
    }

    for (i = 0; i < N; i++) {
        ap[1][i] = x[i] + 0.5 * h * dx[i] + 0.25 * ap[0][i];
        ap[2][i] = dx[i] + 1.0 / h * ap[0][i];
    }

    f(t + 0.5 * h, ap[1], ap[2], ap[3]); 
    for (i = 0; i < N; i++) {
        ap[3][i] = 0.5 * h * h * ap[3][i]; 
    }

    for (i = 0; i < N; i++) {
        ap[4][i] = dx[i] + 1.0 / h * ap[3][i];
    }

    f(t + 0.5 * h, ap[1], ap[4], ap[5]); 
    for (i = 0; i < N; i++) {
        ap[5][i] = 0.5 * h * h * ap[5][i]; 
    }

    for (i = 0; i < N; i++) {
        ap[6][i] = x[i] + h * dx[i] + ap[5][i];
        ap[7][i] = dx[i] + 2.0 / h * ap[5][i];
    }

    f(t + 0.5 * h, ap[6], ap[7], ap[8]); 
    for (i = 0; i < N; i++) {
        ap[8][i] = 0.5 * h * h * ap[8][i]; 
    }

    /* Set the new values */     
    for (i = 0; i < N; i++) {
        new_x[i] = x[i] + h * dx[i] + 1.0/3.0 * (ap[0][i] + ap[3][i] + ap[5][i]);
        new_dx[i] = dx[i] + 1.0/(3.0 * h) * (ap[0][i] + 2 * ap[3][i] + 2 * ap[5][i] + ap[8][i]);
    }

    for (i = 0; i < 9; i++) {
        free(ap[i]);
    }
}

int main(int argc, char** argv) {
    int i, j;

    N = 2;
    h = 0.00001;

    double t_0 = 0;

    x_bodies[0][0] = 0;
    x_bodies[0][1] = 0;
    v_bodies[0][0] = 0;
    v_bodies[0][1] = 0;
    m_bodies[0] = 2e10;//divided by e20
    m_bodies[1] = 6e4;
    
    x_bodies[1][0] = 1.5e2; // divided by e20
    x_bodies[1][1] = 0;
    v_bodies[1][0] = 0;
    v_bodies[1][1] = 25780; 

    FILE* out = fopen("data.txt", "w");
    
    for (i = 0; i < 100000; i++) {
        for (j = 0; j < n; j++) {
            current_body = j;
            runge_kutta(gravity, t_0 + i * h, x_bodies[j], v_bodies[j], x_bodies[j], v_bodies[j]);
        }
        
        fprintf(out, "%f %f\n", x_bodies[1][0], x_bodies[1][1]);
    }

    fclose(out);

    return 0;
}

