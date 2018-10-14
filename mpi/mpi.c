#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#define GENERATION 40

// find center position of kernel
#define K_CENTER_X 1
#define K_CENTER_Y 1

#define KWIDTH 3
#define KHEIGHT 3

int div2blocks(int, int, int, int*, int* );
int Argms_handler(int argc,char** argv, int* image_type,int* width,int* height);

int offset(int Cols,int i,int offs );
void swap_arrays(unsigned char** A, unsigned char** B);
int ImageChanged(unsigned char* A,unsigned char* B,int array_size);

void convolute(unsigned char* A, unsigned char* B, float** h, int row_start, int row_end, int col_start, int col_end, int width, int image_type);
void convolute_grey(unsigned char* A, unsigned char* B, int i, int j, float** h, int width /*cols*/);
void convolute_rgb(unsigned char* A, unsigned char* B, int i, int j, float** h, int width /*cols*/);



int main(int argc, char** argv) {



        int width=0;
        int height=0;


        char* image_name = malloc((strlen(argv[1])+1) * sizeof(char));
        strcpy(image_name, argv[1]);
        int image_type;

        int rows_per_block = -1;
        int cols_per_block = -1;

        int flag1=0;
        int flag2=0;
        int flag11=0;
        int flag22=0;

        double start_time, end_time, elapsed_time, elapsed;

        int rooted_num_procs;
        int changed;

        MPI_Datatype LineGreyType;
        MPI_Datatype LineRgbType;
        MPI_Datatype ColGreyType;
        MPI_Datatype ColRgbType;

        MPI_File fh;
        MPI_Status status;


        int process_id, num_processes;
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
        MPI_Comm_rank(MPI_COMM_WORLD, &process_id);

        if (Argms_handler(argc, argv, &image_type, &width, &height)) {
                printf("ERROR:Something went wrong in the arguments!!\n");
                exit(1);
        }


        if (process_id == 0) {

                rooted_num_procs = div2blocks(height, width, num_processes, &rows_per_block, &cols_per_block);
                // printf("rows_per_block = %d | cols_per_block = %d \n",rows_per_block,cols_per_block );
        }

        MPI_Bcast(&rooted_num_procs, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&rows_per_block, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&cols_per_block, 1, MPI_INT, 0, MPI_COMM_WORLD);


        int gaussian_blur[3][3] = {{1, 2, 1}, {2, 4, 2}, {1, 2, 1}};

        float **h = malloc(sizeof(float*));
        for (int i = 0; i < 3; i++) {
                h[i] =malloc(3* sizeof(float));
                for (int j = 0; j < 3; j++) {
                        h[i][j] = gaussian_blur[i][j] / 16.0;
                }
        }




        int start_row_per_proc = (process_id / rooted_num_procs) *rows_per_block;
        int start_col_per_proc = (process_id % rooted_num_procs) * cols_per_block;

        int array_size;
        unsigned char* block_array_bef;
        unsigned char* block_array_after;


        if (image_type == 0) {

                array_size   = (rows_per_block + 2)*(cols_per_block + 2)* sizeof(unsigned char);
                block_array_bef   = malloc(array_size);
                block_array_after = malloc(array_size);
        }
        else{
                array_size   = (rows_per_block + 2)*(3*cols_per_block + 2*3)* sizeof(unsigned char);
                block_array_bef   = malloc(array_size);
                block_array_after = malloc(array_size);
        }

        //create datatypes for line and column FOR RGB AND GREY;
        //**ATTENTION** keep in mind that LineType is of size cols+2
        MPI_Type_contiguous(cols_per_block+2, MPI_BYTE, &LineGreyType);
        MPI_Type_commit(&LineGreyType);

        MPI_Type_contiguous(3*cols_per_block+ 2*3, MPI_BYTE, &LineRgbType);
        MPI_Type_commit(&LineRgbType);


        MPI_Type_vector(rows_per_block, 1, cols_per_block+2, MPI_BYTE, &ColGreyType);
        MPI_Type_commit(&ColGreyType);

        MPI_Type_vector(rows_per_block, 3, 3*cols_per_block + 2*3, MPI_BYTE, &ColRgbType);
        MPI_Type_commit(&ColRgbType);


        int error =MPI_File_open(MPI_COMM_WORLD, image_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
        if (error) {
                printf( "Unable to open the given input file!\n" ); fflush(stdout);

        }


        if (image_type == 0) {
                for (int row = 0; row < rows_per_block; row++) {

                        MPI_File_seek( fh, (start_row_per_proc + row ) * width + start_col_per_proc, MPI_SEEK_SET );
                        MPI_File_read(fh, &block_array_bef[offset(cols_per_block +2, row+1, 1)], cols_per_block, MPI_BYTE, &status);

                }
        }
        else{
                for (int row = 0; row < rows_per_block; row++) {

                        MPI_File_seek( fh, 3*(start_row_per_proc + row ) * width + 3*start_col_per_proc, MPI_SEEK_SET );
                        MPI_File_read(fh, &block_array_bef[offset(cols_per_block*3 +2*3, row+1, 3)], 3*cols_per_block, MPI_BYTE, &status);

                }

        }

        MPI_File_close(&fh);



        MPI_Request send_row_n, recv_row_n;
        MPI_Request send_row_s, recv_row_s;

        int north_proc = -1;
        int west_proc  = -1;
        int east_proc  = -1;
        int south_proc = -1;


        if (start_row_per_proc != 0) {
                north_proc = process_id - rooted_num_procs;
        }
        if (start_row_per_proc + rows_per_block != height) {
                south_proc = process_id + rooted_num_procs;
        }
        if (start_col_per_proc != 0) {
                west_proc = process_id - 1;
        }
        if (start_col_per_proc + cols_per_block != width) {
                east_proc = process_id + 1;
        }


        MPI_Barrier(MPI_COMM_WORLD);
        /* Get time before */
        start_time = MPI_Wtime();

        for (int i = 0; i < GENERATION; i++) {

                if ( image_type == 0) {
                        if (west_proc != -1) {

                                MPI_Send(&block_array_bef[offset(cols_per_block +2, 1, 1)], 1, ColGreyType, west_proc, 0, MPI_COMM_WORLD);
                                MPI_Recv(&block_array_bef[offset(cols_per_block +2, 1, 0)], 1, ColGreyType, west_proc, 0, MPI_COMM_WORLD, &status);
                        }

                        if (east_proc != -1) {

                                MPI_Recv(&block_array_bef[offset(cols_per_block +2, 1, cols_per_block + 1)], 1, ColGreyType, east_proc, 0, MPI_COMM_WORLD, &status);
                                MPI_Send(&block_array_bef[offset(cols_per_block +2, 1, cols_per_block)], 1, ColGreyType, east_proc, 0, MPI_COMM_WORLD);

                        }


                        if (north_proc != -1) {

                                MPI_Isend(&block_array_bef[offset(cols_per_block+2, 1,0)], 1, LineGreyType, north_proc, 0, MPI_COMM_WORLD, &send_row_n);
                                MPI_Irecv(&block_array_bef[offset(cols_per_block+2, 0,0)], 1, LineGreyType, north_proc, 0, MPI_COMM_WORLD, &recv_row_n);
                        }


                        if (south_proc != -1) {

                                MPI_Isend(&block_array_bef[offset(cols_per_block+2, rows_per_block, 0 )], 1, LineGreyType, south_proc, 0, MPI_COMM_WORLD, &send_row_s);
                                MPI_Irecv(&block_array_bef[offset(cols_per_block+2, rows_per_block+1, 0 )], 1, LineGreyType, south_proc, 0, MPI_COMM_WORLD, &recv_row_s);
                        }
                }

                else{

                        if (west_proc != -1) {

                                MPI_Send(&block_array_bef[offset(3*cols_per_block +2*3, 1, 1*3)], 1, ColRgbType, west_proc, 0, MPI_COMM_WORLD);
                                MPI_Recv(&block_array_bef[offset(3*cols_per_block +2*3, 1, 0)], 1, ColRgbType, west_proc, 0, MPI_COMM_WORLD, &status);
                        }

                        if (east_proc != -1) {

                                MPI_Recv(&block_array_bef[offset(3*cols_per_block +2*3, 1, 3*cols_per_block + 1*3)], 1, ColRgbType, east_proc, 0, MPI_COMM_WORLD, &status);
                                MPI_Send(&block_array_bef[offset(3*cols_per_block +2*3, 1, 3*cols_per_block)], 1, ColRgbType, east_proc, 0, MPI_COMM_WORLD);

                        }


                        if (north_proc != -1) {

                                MPI_Isend(&block_array_bef[offset(3*cols_per_block+2*3, 1,0)], 1, LineRgbType, north_proc, 0, MPI_COMM_WORLD, &send_row_n);
                                MPI_Irecv(&block_array_bef[offset(3*cols_per_block+2*3, 0,0)], 1, LineRgbType, north_proc, 0, MPI_COMM_WORLD, &recv_row_n);
                        }


                        if (south_proc != -1) {

                                MPI_Isend(&block_array_bef[offset(3*cols_per_block+2*3, rows_per_block, 0 )], 1, LineRgbType, south_proc, 0, MPI_COMM_WORLD, &send_row_s);
                                MPI_Irecv(&block_array_bef[offset(3*cols_per_block+2*3, rows_per_block+1, 0 )], 1, LineRgbType, south_proc, 0, MPI_COMM_WORLD, &recv_row_s);
                        }

                }


                // convoluteInner()
                convolute(block_array_bef, block_array_after, h, 1, rows_per_block, 1, cols_per_block,cols_per_block,image_type);



                flag1=0;
                flag2=0;
                flag11 =0;
                flag22=0;
                while (1) {
                        if (north_proc != -1) {
                                MPI_Test(&recv_row_n, &flag1, &status);
                                if (flag1 == 1) {
                                        // convolute(up);
                                        convolute(block_array_bef, block_array_after, h, 1, 1, 1, cols_per_block, cols_per_block,image_type);

                                        flag11=1;
                                }
                        }
                        else{
                                flag11 =1;
                        }
                        if (south_proc != -1) {
                                MPI_Test(&recv_row_s, &flag2, &status);
                                if (flag2 == 1) {
                                        // convolute(down)
                                        convolute(block_array_bef, block_array_after, h, rows_per_block, rows_per_block, 1, cols_per_block, cols_per_block,image_type);
                                        flag22=1;
                                }
                        }
                        else{
                                flag22 =1;
                        }

                        if (flag11 == 1 && flag22 == 1) {
                                break;
                        }
                }



                if (west_proc != -1) {
                        convolute(block_array_bef, block_array_after, h, 1, rows_per_block, 1, 1, cols_per_block,image_type);
                }
                if (east_proc != -1) {
                        convolute(block_array_bef, block_array_after, h, 1, rows_per_block, cols_per_block, cols_per_block, cols_per_block,image_type);
                }
                if (north_proc != -1) {
                        MPI_Wait(&recv_row_n, &status);
                        convolute(block_array_bef, block_array_after, h, 1, 1, 1, cols_per_block, cols_per_block,image_type);
                }
                if (south_proc != -1) {
                        MPI_Wait(&recv_row_s, &status);
                        convolute(block_array_bef, block_array_after, h, rows_per_block, rows_per_block, 1, cols_per_block, cols_per_block,image_type);
                }


                /* Wait to have sent all borders */
                if (north_proc != -1)
                        MPI_Wait(&send_row_n, &status);
                if (south_proc != -1)
                        MPI_Wait(&send_row_s, &status);


                swap_arrays(&block_array_bef, &block_array_after);

                /*check if image changed */
                changed =0;
                changed = ImageChanged(block_array_bef,block_array_after,array_size);
                if ( changed ) {
                        // printf("GENERATION : %d Image Changed = %d \n",i,changed );
                }
                else{
                        // printf("GENERATION : %d Image NOT Changed = %d \n",i,changed );
                        break;
                }


        }

        /* Get time elapsed */
        end_time = MPI_Wtime();
        elapsed_time = end_time - start_time;


        /* Parallel write */
        char *outImage = malloc((strlen(image_name) + strlen("output_images/filtered_")) * sizeof(char));
        //char *outImage = malloc((strlen(image_name) + strlen("filtered_")) * sizeof(char));

        strcpy(outImage, "output_images/filtered_");
        //strcpy(outImage, "filtered_");

        strcat(outImage, image_name);
        MPI_File outFile;
        error = MPI_File_open(MPI_COMM_WORLD, outImage, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &outFile);
        if (error) {
                printf( "Unable to open the output file..No space! \n" ); fflush(stdout);
        }


        if (image_type == 0) {
                for (int row = 0; row < rows_per_block; row++) {

                        MPI_File_seek(outFile, (start_row_per_proc + row) * width + start_col_per_proc, MPI_SEEK_SET);

                        MPI_File_write(outFile, &block_array_bef[offset(cols_per_block+2, row+1, 1)], cols_per_block, MPI_BYTE, MPI_STATUS_IGNORE);
                }
        }
        else{
                for (int row = 0; row < rows_per_block; row++) {

                        MPI_File_seek(outFile, 3*(start_row_per_proc + row) * width + 3*start_col_per_proc, MPI_SEEK_SET);

                        MPI_File_write(outFile, &block_array_bef[offset(3*cols_per_block+2*3, row+1, 3)], 3*cols_per_block, MPI_BYTE, MPI_STATUS_IGNORE);
                }
        }

        MPI_File_close(&outFile);


        MPI_Reduce(&elapsed_time, &elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if(process_id == 0)
                printf("Execution time: %f\n",elapsed);



        free(block_array_bef);
        free(block_array_after);
        free(image_name);
        free(outImage);

        MPI_Type_free(&LineGreyType);
        MPI_Type_free(&ColGreyType);
        MPI_Type_free(&LineRgbType);
        MPI_Type_free(&ColRgbType);


        // Finalize the MPI environment. No more MPI calls can be made after this
        MPI_Finalize();

        return 0;
}

/*
   ███████ ██    ██ ███    ██  ██████ ████████ ██  ██████  ███    ██ ███████
   ██      ██    ██ ████   ██ ██         ██    ██ ██    ██ ████   ██ ██
   █████   ██    ██ ██ ██  ██ ██         ██    ██ ██    ██ ██ ██  ██ ███████
   ██      ██    ██ ██  ██ ██ ██         ██    ██ ██    ██ ██  ██ ██      ██
   ██       ██████  ██   ████  ██████    ██    ██  ██████  ██   ████ ███████
 */


inline void convolute(unsigned char* A, unsigned char* B, float** h, int row_start, int row_end, int col_start, int col_end, int width /*cols*/,int image_type) {
        int i, j;

        if (image_type == 0) {

                for (i = row_start; i <= row_end; ++i)
                        for (j = col_start; j <= col_end; ++j)
                                convolute_grey(A, B, i, j, h, width);
        }
        else{

                for (i = row_start; i <= row_end; ++i)
                        for (j = col_start; j <= col_end; j++)
                                convolute_rgb(A, B, i, 3*j, h, 3*width + 6);
        }
}

inline void convolute_rgb(unsigned char* A, unsigned char* B, int i, int j, float** h, int width /*cols*/)
{
        float val_red = 0;
        float val_green = 0;
        float val_blue = 0;
        int ii=0;
        int jj=0;
        for (int m = 0; m < 3; ++m)   // kernel rows
        {
                for (int n = 0; n < 3; ++n)   // kernel columns
                {

                        ii = i + (m - K_CENTER_Y);

                        jj = j + 3*(n - K_CENTER_X);


                        val_red+=A[(width)*ii + jj] * h[m][n];
                        val_green+=A[(width)*ii + jj+1] * h[m][n];
                        val_blue+=A[(width)*ii + jj+2] * h[m][n];


                }
        }

        B[(width)*i + j] = val_red;
        B[(width)*i + j+1] = val_green;
        B[(width)*i + j+2] = val_blue;


}

inline void convolute_grey(unsigned char* A, unsigned char* B, int i, int j, float** h, int width /*cols*/)
{
        int ii=0;
        int jj=0;
        float val=0;

        for (int m = 0; m < 3; ++m)   // kernel rows
        {
                for (int n = 0; n < 3; ++n)   // kernel columns
                {
                        ii = i + (m - K_CENTER_Y);   // i +m -1 ----> curr_row_block + curr_row_kernel -1

                        jj = j + (n - K_CENTER_X);   // j +n -1 ----> curr_col_block + curr_col_kernel -1


                        val+=A[(width+2)*ii + jj] * h[m][n];


                }
        }

        B[(width+2)*i + j] = val;   // P[i][j] += N[ii][jj] * M[m][n];

}



inline void swap_arrays(unsigned char** A, unsigned char** B){
        unsigned char *temp;
        temp = *A;
        *A = *B;
        *B = temp;
        return;
}


inline int Argms_handler(int argc,char** argv, int* image_type,int* width,int* height){
        if (argc != 5) {
                return -1;
        }
        *width = atoi(argv[3]);
        *height = atoi(argv[4]);

        if (strcmp(argv[2], "grey") == 0 || strcmp(argv[2], "GREY") == 0 ) {
                *image_type = 0;
        } else if (strcmp(argv[2], "rgb") == 0 || strcmp(argv[2], "RGB") == 0 ) {
                *image_type = 1;
        }else{
                return -2;
        }

        return 0;
}

inline int offset(int k,int i,int offs ){  //k = cols or rows
        return k*i + offs;
}


inline int ImageChanged(unsigned char * A, unsigned char * B, int array_size)
{
        int i;
        for(i=0; i<array_size; i++) {
                if(A[i]!= B[i])
                        return 1;
        }
        return 0;
}



inline int div2blocks(int height, int width, int num_processes, int* rows_per_block, int* cols_per_block ){

        double double_rows_per_block,double_cols_per_block;

        double root_num_processes = sqrt(num_processes);

        if ((root_num_processes - floor(root_num_processes)) > 0 ) {
                printf("The root of number of processes doesn't give an integer! \n" );
                exit(1);
        }
        int int_root_num_processes = root_num_processes;
        if ((height % int_root_num_processes) && (width % int_root_num_processes)   ) {
                printf("The rows & cols can't be partitioned equally to the processes! \n" );
                exit(1);
        }
        double_rows_per_block =(double) height/root_num_processes;
        double_cols_per_block =(double) width/root_num_processes;

        if (double_rows_per_block - floor(double_rows_per_block ) > 0 ) {
                printf("The rows that are splitted to processes are not integers! \n" );
                exit(1);
        }
        if (double_cols_per_block - floor(double_cols_per_block) > 0) {
                printf("The columns that are splitted to processes are not integers! \n" );
                exit(1);
        }
        *rows_per_block =height/int_root_num_processes;
        *cols_per_block =width/int_root_num_processes;

        return int_root_num_processes;

}
