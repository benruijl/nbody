#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int N; // dimension
double h;
int n; // n bodies
double** x_bodies;
double** v_bodies;
double* m_bodies;
int current_body;

static double G = 6.674e-11;

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
            dist = 1.0 / (dist * sqrt(dist)); 
            
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
    int i, j, k;
    int iterations;
    double t_0 = 0;
    double** x_new;
    double** v_new;

    if (argc < 2) {
        printf("No configuration file specified.");
        return 1;
    }

    FILE* config = fopen(argv[1], "r");

    if (!config) {
        printf("Could not open file: %s", argv[1]);
        return 2;
    }

    fscanf(config, "%d", &iterations);
    fscanf(config, "%lf", &h);
    fscanf(config, "%lf", &t_0);
    fscanf(config, "%d", &N);
    fscanf(config, "%d", &n);
    
    x_bodies = (double**) malloc(n * sizeof(double*));
    v_bodies = (double**) malloc(n * sizeof(double*));
    m_bodies = (double*) malloc(n * sizeof(double));
    
    x_new = (double**) malloc(n * sizeof(double*));
    v_new = (double**) malloc(n * sizeof(double*));

    for (i = 0; i < n; i++) {
        x_bodies[i] = (double*) malloc(N * sizeof(double));
        v_bodies[i] = (double*) malloc(N * sizeof(double));
        x_new[i] = (double*) malloc(N * sizeof(double));
        v_new[i] = (double*) malloc(N * sizeof(double));
        
        for (j = 0; j < N; j++) {
                fscanf(config, "%lf", &x_bodies[i][j]);
        }
        
        for (j = 0; j < N; j++) {
                fscanf(config, "%lf", &v_bodies[i][j]);
        }

        fscanf(config, "%lf", &m_bodies[i]);
    }

    FILE* out = fopen("data.txt", "w");

    for (i = 0; i < iterations; i++) {
        for (j = 0; j < n; j++) {
            current_body = j;
            runge_kutta(gravity, t_0 + i * h, x_bodies[j], v_bodies[j], x_new[j], v_new[j]);
        }

        for (j = 0; j < n; j++) {
            double* temp = x_bodies[j];
            x_bodies[j] = x_new[j];
            x_new[j] = temp;

            temp = v_bodies[j];
            v_bodies[j] = v_new[j];
            v_new[j] = temp;

            for (k = 0; k < N; k++) {
                fprintf(out, "%lf ", x_bodies[j][k]);
            }
        }
        fprintf(out, "\n");
    }

    fclose(out);

    return 0;
}

