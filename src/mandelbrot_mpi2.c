#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "mpi.h"

#define LEADER 0
#define TAG 10

int block_height;

double c_x_min;
double c_x_max;
double c_y_min;
double c_y_max;

double pixel_width;
double pixel_height;

int iteration_max = 200;

int image_size;
unsigned char **image_buffer;

int i_x_max;
int i_y_max;
int image_buffer_size;

int rgb_size;

int gradient_size = 16;
int colors[17][3] = {
                        {66, 30, 15},
                        {25, 7, 26},
                        {9, 1, 47},
                        {4, 4, 73},
                        {0, 7, 100},
                        {12, 44, 138},
                        {24, 82, 177},
                        {57, 125, 209},
                        {134, 181, 229},
                        {211, 236, 248},
                        {241, 233, 191},
                        {248, 201, 95},
                        {255, 170, 0},
                        {204, 128, 0},
                        {153, 87, 0},
                        {106, 52, 3},
                        {16, 16, 16},
                    };

void allocate_image_buffer(){
    /* Variable image_buffer_size is different at each process. While the leader */
    /* allocates the full buffer, the followers allocates a partial buffer. */
    rgb_size = 3;
    image_buffer = (unsigned char **) malloc(sizeof(unsigned char *) * image_buffer_size);

    for(int i = 0; i < image_buffer_size; i++){
            image_buffer[i] = (unsigned char *) malloc(sizeof(unsigned char) * rgb_size);
        };
};

void init(int argc, char *argv[], int my_id, int num_tasks){
    if(argc < 6){
        printf("usage: ./mandelbrot_seq c_x_min c_x_max c_y_min c_y_max image_size\n");
        printf("examples with image_size = 11500:\n");
        printf("    Full Picture:         ./mandelbrot_seq -2.5 1.5 -2.0 2.0 11500\n");
        printf("    Seahorse Valley:      ./mandelbrot_seq -0.8 -0.7 0.05 0.15 11500\n");
        printf("    Elephant Valley:      ./mandelbrot_seq 0.175 0.375 -0.1 0.1 11500\n");
        printf("    Triple Spiral Valley: ./mandelbrot_seq -0.188 -0.012 0.554 0.754 11500\n");
        exit(0);
    }
    else{
        sscanf(argv[1], "%lf", &c_x_min);
        sscanf(argv[2], "%lf", &c_x_max);
        sscanf(argv[3], "%lf", &c_y_min);
        sscanf(argv[4], "%lf", &c_y_max);
        sscanf(argv[5], "%d", &image_size);

        i_x_max           = image_size;
        i_y_max           = image_size;

        pixel_width       = (c_x_max - c_x_min) / i_x_max;
        pixel_height      = (c_y_max - c_y_min) / i_y_max;

        // printf("TASK %d -- block_height = %d / %d\n", my_id, i_y_max, num_tasks);

        block_height = (int) i_y_max / num_tasks;

        // printf("TASK %d -- block_height = %f\n", my_id, block_height);

        /* The leader is responsible for the full buffer. */
        if(my_id == LEADER){
            image_buffer_size = image_size * image_size;
        }
        /* followers is responsible for the local buffer. */
        else{
            image_buffer_size = image_size * block_height;
        }
    };
};

void update_rgb_buffer(int iteration, int x, int y, int my_id){
    int color, buffer_height;
    
    /* The leader updates directly in the full buffer. */
    if(my_id == LEADER){
        buffer_height = i_y_max;
    }
    /* The followers updates in the partial buffer. */
    else{
        buffer_height = block_height;
    }

    if(iteration == iteration_max){
        image_buffer[(buffer_height * y) + x][0] = colors[gradient_size][0];
        image_buffer[(buffer_height * y) + x][1] = colors[gradient_size][1];
        image_buffer[(buffer_height * y) + x][2] = colors[gradient_size][2];
    }
    else{
        color = iteration % gradient_size;

        image_buffer[(buffer_height * y) + x][0] = colors[color][0];
        image_buffer[(buffer_height * y) + x][1] = colors[color][1];
        image_buffer[(buffer_height * y) + x][2] = colors[color][2];
    };
};

void write_to_file(){
    FILE * file;
    char * filename               = "output.ppm";
    char * comment                = "# ";

    int max_color_component_value = 255;

    file = fopen(filename,"wb");

    fprintf(file, "P6\n %s\n %d\n %d\n %d\n", comment,
            i_x_max, i_y_max, max_color_component_value);

    for(int i = 0; i < image_buffer_size; i++){
        fwrite(image_buffer[i], 1 , 3, file);
    };

    fclose(file);
};

void compute_mandelbrot(int i_y_start, int i_y_end, int my_id){
    double z_x;
    double z_y;
    double z_x_squared;
    double z_y_squared;
    double escape_radius_squared = 4;

    int iteration;
    int i_x;
    int i_y, i_yy;

    double c_x;
    double c_y;

    i_yy = 0;

    for(i_y = i_y_start; i_y <= i_y_end; i_y++){
        c_y = c_y_min + i_y * pixel_height;

        if(fabs(c_y) < pixel_height / 2){
            c_y = 0.0;
        };

        // Transformation to the local buffer:
        // i_yy = i_y % block_height;

        for(i_x = 0; i_x < i_x_max; i_x++){
            c_x         = c_x_min + i_x * pixel_width;

            z_x         = 0.0;
            z_y         = 0.0;

            z_x_squared = 0.0;
            z_y_squared = 0.0;

            for(iteration = 0;
                iteration < iteration_max && \
                ((z_x_squared + z_y_squared) < escape_radius_squared);
                iteration++){
                z_y         = 2 * z_x * z_y + c_y;
                z_x         = z_x_squared - z_y_squared + c_x;

                z_x_squared = z_x * z_x;
                z_y_squared = z_y * z_y;
            };
            update_rgb_buffer(iteration, i_x, i_yy, my_id);
        };
        i_yy++;
        // printf("TASK %d -- i_y = %d and i_y' = %d\n", my_id, i_y, i_yy);
    };
};

int main(int argc, char *argv[]){
    int num_tasks, task_id, height_start, height_end;
    
    MPI_Status mpi_status;
    
    /* Initialization */
    MPI_Init(&argc, &argv);
    
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &task_id);
    
    // printf("MPI task %d [of %d] has started...\n", task_id, num_tasks);

    init(argc, argv, task_id, num_tasks);

    // printf("TASK %d -- image_size = %d and block_height = %d\n", task_id, image_size, block_height);

    allocate_image_buffer();

    // printf("TASK %d -- image_buffer_size = %d\n", task_id, image_buffer_size);

    // printf("TASK %d -- block_height = %.1f\n", task_id, block_height);

    if(task_id == LEADER){
        height_start = 0;
        height_end = block_height - 1;
    }
    else{
        height_start = task_id * block_height;
        height_end = (task_id + 1) * block_height - 1;
    }

    // printf("TASK %d -- height_start = %d and height_end = %d\n", task_id, height_start, height_end);

    compute_mandelbrot(height_start, height_end, task_id);

    if(task_id == LEADER){
        int busy_followers = num_tasks - 1;
        int recv_buffer_size = image_size * block_height * rgb_size;

        unsigned char **aux_buffer = (unsigned char **) malloc(sizeof(unsigned char *) * recv_buffer_size);

        for(int i = 0; i < recv_buffer_size; i++){
                aux_buffer[i] = (unsigned char *) malloc(sizeof(unsigned char) * rgb_size);
            };
        
        // while(busy_followers--){
        //     // MPI_Recv(&(aux_buffer[0][0]), recv_buffer_size, MPI_UNSIGNED_CHAR, MPI_ANY_SOURCE, TAG, MPI_COMM_WORLD, &mpi_status);
        //     // int source_id = mpi_status.MPI_SOURCE;

        //     // for(int i = 0; i < block_height; i++){
        //     //     int ii = block_height*source_id;
        //     //     for(int j = 0; j < i_x_max; j++){
        //     //         image_buffer[ii + j] = aux_buffer[(block_height * i) + j];
        //     //         ii += block_height;
        //     //     }
        //     // }
        //     // memcpy(&image_buffer[0], aux_buffer, sizeof(unsigned char) * recv_buffer_size);
        // }

        // write_to_file();
    }
    else{
        if(task_id == 1){
                write_to_file();
            }
        // MPI_Send(&(image_buffer[0][0]), image_buffer_size*rgb_size, MPI_UNSIGNED_CHAR, LEADER, TAG, MPI_COMM_WORLD);
    }

    // init(argc, argv);

    // allocate_image_buffer();

    // compute_mandelbrot();

    // write_to_file();

    /* End MPI code */
    MPI_Finalize();
};
