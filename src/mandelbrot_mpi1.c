#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define LEADER 0

int step_size;

int rgb_size = 3;

double c_x_min;
double c_x_max;
double c_y_min;
double c_y_max;

double pixel_width;
double pixel_height;

int iteration_max = 200;

int image_size;
unsigned char **image_buffer;

int i_x_min;
int i_y_min;
int i_x_max;
int i_y_max;
int image_buffer_size;

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
    image_buffer = (unsigned char **) malloc(sizeof(unsigned char *) * image_buffer_size);

    for(int i = 0; i < image_buffer_size; i++){
        image_buffer[i] = (unsigned char *) malloc(sizeof(unsigned char) * rgb_size);
    };
};

void init(int argc, char *argv[], int my_id, int numtasks){
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

        i_x_max = image_size;
        i_y_max = image_size;

        pixel_width  = (c_x_max - c_x_min) / i_x_max;
        pixel_height = (c_y_max - c_y_min) / i_y_max;

        step_size = image_size / numtasks;

        /* LEADER task will allocate full buffer. */
        if(my_id == LEADER){
            image_buffer_size = image_size * image_size;
        }
        /* Other tasks will allocate local buffer. */
        else{
            image_buffer_size = image_size * step_size;
        }
    };
};

void update_rgb_buffer(int iteration, int x, int y){
    int color;

    if(iteration == iteration_max){
        image_buffer[(step_size * y) + x][0] = colors[gradient_size][0];
        image_buffer[(step_size * y) + x][1] = colors[gradient_size][1];
        image_buffer[(step_size * y) + x][2] = colors[gradient_size][2];
    }
    else{
        color = iteration % gradient_size;

        image_buffer[(step_size * y) + x][0] = colors[color][0];
        image_buffer[(step_size * y) + x][1] = colors[color][1];
        image_buffer[(step_size * y) + x][2] = colors[color][2];
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

void compute_mandelbrot(int my_id){
    double z_x;
    double z_y;
    double z_x_squared;
    double z_y_squared;
    double escape_radius_squared = 4;

    int iteration;
    int i_x;
    int i_y;

    double c_x;
    double c_y;

    for(i_y = 0; i_y < step_size; i_y++){
        c_y = (c_y_min + my_id*step_size) + i_y * pixel_height;

        if(fabs(c_y) < pixel_height / 2){
            c_y = 0.0;
        };

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

            update_rgb_buffer(iteration, i_x, i_y);
        };
    };
};

int main(int argc, char *argv[]){
    int numtasks, taskid;

    /*** INITIALIZATIONS ***/
    MPI_Init(&argc, &argv);
    MPI_Status *status;
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

    init(argc, argv, taskid, numtasks);

    printf("==> image_buffer_size = %d at task %d\n", image_buffer_size, taskid);

    /* FOLLOWER tasks allocates local buffer. */
    allocate_image_buffer();

    compute_mandelbrot(taskid);

    /* TODO: Send and Recive local buffer and copy it to full buffer. */
    if(taskid == LEADER){
        int running_followers = numtasks - 1;
        int recv_buffer_size = image_size * step_size;
        while(running_followers--){
            for(int i = 0; i < recv_buffer_size; i++){
                unsigned char tmp_rgb[rgb_size];
                MPI_Recv(&tmp_rgb, rgb_size, MPI_UNSIGNED_CHAR, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, status);
                memcpy(&image_buffer[status->MPI_SOURCE * step_size + i], tmp_rgb, sizeof(unsigned char) * rgb_size);
            }
        }
    }
    else{
        for(int i = 0; i < image_buffer_size; i++){
            MPI_Send(&image_buffer[i], rgb_size, MPI_UNSIGNED_CHAR, LEADER, 1, MPI_COMM_WORLD);
        }
    }


    if(taskid == LEADER){
        printf("Writing to file...");
        write_to_file();
    }

    MPI_Finalize();

    /*init(argc, argv);
    allocate_image_buffer();

    compute_mandelbrot();

    write_to_file();

    return 0;*/
};
