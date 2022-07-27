#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <pthread.h>

#define NUM_THREADS 4

double limit_size;

/* Structs to store mandelbrot arguments for thread_create() */

typedef struct {
    int min;
    int max;
} Limits;

typedef struct {
    Limits limits;
    int height_start;
    int height_end;
} MPI_Args;

pthread_t tid[NUM_THREADS];

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
    int rgb_size = 3;
    image_buffer = (unsigned char **) malloc(sizeof(unsigned char *) * image_buffer_size);

    for(int i = 0; i < image_buffer_size; i++){
        image_buffer[i] = (unsigned char *) malloc(sizeof(unsigned char) * rgb_size);
    };
};

void contiguous_image_buffer_allocation(){
    unsigned char *data;
    int data_size;
    int rgb_size = 3;

    /* N pointers of M objects => N*M objects */
    data_size = image_buffer_size * rgb_size;

    /* Allocate spaces in memory for data and buffer */
    data = (unsigned char *) malloc(sizeof(unsigned char) * data_size);
    image_buffer = (unsigned char **) malloc(sizeof(unsigned char *) * image_buffer_size);

    /* Contiguously allocation of data and buffer */
    for(int i = 0; i < image_buffer_size; i++){
        image_buffer[i] = &(data[rgb_size*i]);
    }
}

void init(int argc, char *argv[]){
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
        image_buffer_size = image_size * image_size;

        pixel_width       = (c_x_max - c_x_min) / i_x_max;
        pixel_height      = (c_y_max - c_y_min) / i_y_max;

        limit_size = i_x_max / NUM_THREADS; // sets the size of the boundary to split the figure
    };
};

void update_rgb_buffer(int iteration, int x, int y){
    int color;

    if(iteration == iteration_max){
        image_buffer[(i_y_max * y) + x][0] = colors[gradient_size][0];
        image_buffer[(i_y_max * y) + x][1] = colors[gradient_size][1];
        image_buffer[(i_y_max * y) + x][2] = colors[gradient_size][2];
    }
    else{
        color = iteration % gradient_size;

        image_buffer[(i_y_max * y) + x][0] = colors[color][0];
        image_buffer[(i_y_max * y) + x][1] = colors[color][1];
        image_buffer[(i_y_max * y) + x][2] = colors[color][2];
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

void *compute_mandelbrot(void *inputs){
    /* Convert the inputs. */
    MPI_Args *args = (MPI_Args *)inputs;
    
    /* Capture MPI bounds. */
    int height_start = args->height_start;
    int height_end = args->height_end;

    /* Capture thread bounds. */
    Limits limits = args->limits;
    int x_min = limits.min;
    int x_max = limits.max;

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

    for(i_y = height_start; i_y < height_end; i_y++){
        c_y = c_y_min + i_y * pixel_height;

        if(fabs(c_y) < pixel_height / 2){
            c_y = 0.0;
        };

        /* Compute the set within the defined limits. */
        for(i_x = x_min; i_x < x_max; i_x++){
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
    int my_rank, comm_size;
    int height_start, height_end;
    int local_height;
    MPI_Status status;

    
    /* MPI initialization */
    MPI_Init(&argc, &argv);
    
    /* Read the inputs */
    init(argc, argv);

    /* Get process ID and communication size */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    /* Divide the height into blocks */
    int height_step = i_y_max / comm_size;

    /* Defines the boundaries of each process */
    if(my_rank == 0){
        height_start = 0;
        local_height = i_y_max - (comm_size - 1)*height_step;
        height_end = local_height;
    }
    else{
        local_height = height_step;
        height_start = (i_y_max - (comm_size - 1)*height_step) + (my_rank - 1)*local_height;
        height_end = height_start + height_step;
    }

    /* Contiguous allocation of image buffer */
    contiguous_image_buffer_allocation();

    /* Sets the limits of each thread. */
    Limits limits[NUM_THREADS];
    for(int i=0; i < NUM_THREADS; i++) {
        limits[i].min = i*limit_size;
        limits[i].max = (i+1)*limit_size;
    }

    /* Struct to store compute_mandelbrot() arguments */
    MPI_Args args[NUM_THREADS];
    for(int i=0; i < NUM_THREADS; i++) {
        args[i].limits.min = limits[i].min;
        args[i].limits.max = limits[i].max;
        args[i].height_start = height_start;
        args[i].height_end = height_end;
    }

    /* Each process computes the mandelbrot within its limits */
    for(int i=0; i < NUM_THREADS; i++) {
        /* Create threads to compute Mandelbrot set. */
        pthread_create(&tid[i], NULL, compute_mandelbrot, &args[i]);
    }

    for(int i=0; i < NUM_THREADS; i++) {
        /* Wait for all threads to complete. */
        pthread_join(tid[i], NULL);
    }

    /* Process 0 gathers the calculation of the other */
    /* processes and writes it to the file            */
    if(my_rank == 0){
        int recv_start;
        for(int i = 1; i < comm_size; i++){
            recv_start = (i_y_max - (comm_size - 1)*height_step) + (i - 1)*local_height;
            MPI_Recv(&(image_buffer[(i_y_max * recv_start) + 0][0]), height_step*i_x_max*3, MPI_UNSIGNED_CHAR, i, 1, MPI_COMM_WORLD, &status);
        }
    }
    else{
        MPI_Send(&(image_buffer[(i_y_max * height_start) + 0][0]), height_step*i_x_max*3, MPI_UNSIGNED_CHAR, 0, 1, MPI_COMM_WORLD);
    }

    if(my_rank == 0){
        write_to_file();
    }

    MPI_Finalize();
};
