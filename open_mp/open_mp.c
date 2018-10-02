#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#define GENERATION 20

// find center position of kernel
#define K_CENTER_X 1
#define K_CENTER_Y 1

#define KWIDTH 3
#define KHEIGHT 3

int div2blocks(int, int, int, int*, int* );
int Argms_handler(int argc,char** argv, int* image_type,int* width,int* height);

int offset(int Cols,int i,int offs );
void swap_arrays(unsigned char** A, unsigned char** B);
int ImageChanged(unsigned char* A,unsigned char* B );

void convolute_grey(unsigned char* A, unsigned char* B, float** h, int row_start, int row_end, int col_start, int col_end, int width, int height);


int main(int argc, char** argv) {



        int width=0;
        int height=0;
        unsigned char* temp_array=NULL;
        int number_offs;

        char * image_src_file = "../input_images/";
        char* image_name = malloc((strlen(argv[1]) + strlen(image_src_file) + 1 ) * sizeof(char));
        strcpy(image_name, image_src_file);
        strcat(image_name, argv[1]);

        printf("image:'%s'\n", image_name);

        char pure_image_name[strlen(argv[1]) - 4 + 1]; // Input string
        strncpy(pure_image_name, argv[1], strlen(argv[1]) - 4);
        pure_image_name[strlen(argv[1]) - 4] = '\0'; // null terminate destination
        printf("pure image:'%s'\n", pure_image_name);

        char * resol = malloc((strlen(argv[3]) + strlen(argv[4]) + 1 + 2 + 4) * sizeof(char));
        strcpy(resol, "_");
        strcat(resol, argv[3]);
        strcat(resol, "_");
        strcat(resol, argv[4]);
        strcat(resol, ".raw");
        printf("resolution:'%s'\n", resol);

        char * outImage = malloc((strlen(resol) + strlen(pure_image_name) + 1 + strlen("output_images/filtred_") ) * sizeof(char));
        strcpy(outImage, "output_images/filtred_");
        strcat(outImage, pure_image_name);
        strcat(outImage, resol);

        printf("output image:'%s'\n", outImage);
        FILE * fp = fopen("waterfall_grey.raw", "r");
        if (fp == NULL)
        {
          printf("olo malakies\n");
        }

        exit(1);

        /*char* image_name = malloc((strlen(argv[1])+1) * sizeof(char));
        strcpy(image_name, argv[1]);*/
        int image_type;

        int rows_per_block = -1;
        int cols_per_block = -1;

        int flag1=0;
        int flag2=0;
        int flag11=0;
        int flag22=0;

        double  start_time, end_time, elapsed_time, elapsed;

        int rooted_num_procs;
        int changed=0;

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
        // printf("width = %d , height = %d\n",width,height );


        if (process_id == 0) {

                rooted_num_procs = div2blocks(height, width, num_processes, &rows_per_block, &cols_per_block);
                printf("rows_per_block = %d | cols_per_block = %d \n",rows_per_block,cols_per_block );

        }

        MPI_Bcast(&rooted_num_procs, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&rows_per_block, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&cols_per_block, 1, MPI_INT, 0, MPI_COMM_WORLD);



        // filtes kernel
        // int edge_detection[3][3] = {{1, 4, 1}, {4, 8, 4}, {1, 4, 1}};
        // int box_blur[3][3] = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};
        int gaussian_blur[3][3] = {{1, 2, 1}, {2, 4, 2}, {1, 2, 1}};

        float **h = malloc(sizeof(float*));
        for (int i = 0; i < 3; i++) {
                h[i] =malloc(3* sizeof(float));
                for (int j = 0; j < 3; j++) {

                        h[i][j] = gaussian_blur[i][j] / 16.0;
                        // h[i][j] = edge_detection[i][j] / 28.0;
                        // h[i][j] = box_blur[i][j] / 9.0;

                }
        }




        int start_row_per_proc = (process_id / rooted_num_procs) *rows_per_block;
        int start_col_per_proc = (process_id % rooted_num_procs) * cols_per_block;
        // printf(" start_row = %d per proc : %d \n",start_row_per_proc,process_id );
        // printf(" start_col = %d per proc : %d \n",start_col_per_proc,process_id );


        unsigned char* block_array_bef   = malloc((rows_per_block + 2)*(cols_per_block + 2)* sizeof(unsigned char));
        unsigned char* block_array_after = malloc((rows_per_block + 2)*(cols_per_block + 2)* sizeof(unsigned char));

        //create datatypes for line and column FOR RGB AND GREY;
        //**ATTENTION** keep in mind that LineType is of size cols+2
        MPI_Type_contiguous(cols_per_block+2, MPI_BYTE, &LineGreyType);
        MPI_Type_commit(&LineGreyType);

        MPI_Type_vector(rows_per_block, 1, cols_per_block+2, MPI_BYTE, &ColGreyType);
        MPI_Type_commit(&ColGreyType);


        MPI_Type_contiguous(3*cols_per_block+2, MPI_BYTE, &LineRgbType);
        MPI_Type_commit(&LineRgbType);

        MPI_Type_vector(3*cols_per_block, 1, rows_per_block, MPI_BYTE, &ColRgbType);
        MPI_Type_commit(&ColRgbType);


        int error =MPI_File_open(MPI_COMM_WORLD, image_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
        if (error) {
                printf( "Unable to open the given input file!\n" ); fflush(stdout);

        }


        if (image_type == 0) {
                for (int row = 0; row < rows_per_block; row++) {

                        MPI_File_seek( fh, (start_row_per_proc + row ) * width + start_col_per_proc, MPI_SEEK_SET );
                        // printf("ROW[%d]=FILE READ for proc_id[%d]: offset(cols_per_block +2,row+1,1) = %d\n",row,process_id,offset(cols_per_block +2,row+1,1) );
                        MPI_File_read(fh, &block_array_bef[offset(cols_per_block +2, row+1, 1)], cols_per_block, MPI_BYTE, &status);

                }
        }

        MPI_File_close(&fh);



        MPI_Request send_row_n, recv_row_n;
        // MPI_Request send_col_w, recv_col_w;
        MPI_Request send_row_s, recv_row_s;
        // MPI_Request send_col_e, recv_col_e;

        // MPI_Request send_nw, recv_nw;
        // MPI_Request send_ne, recv_ne;
        // MPI_Request send_sw, recv_sw;
        // MPI_Request send_se, recv_se;

        // MPI_Request *send_req, *receive_req; //send_requestX,receive_requestX
        // send_req = malloc(8 * sizeof(MPI_Request));
        // receive_req = malloc(8 * sizeof(MPI_Request));
        // neighbors
        int north_proc = -1;
        int west_proc  = -1;
        int east_proc  = -1;
        int south_proc = -1;


        if (start_row_per_proc != 0) {
                north_proc = process_id - rooted_num_procs;
        }
        if (start_row_per_proc + rows_per_block != height) {
                // printf("id = %d : start_row_per_proc = %d + rows_per_block = %d != height = %d  \n",process_id,start_row_per_proc,rows_per_block,height );
                south_proc = process_id + rooted_num_procs;
        }
        if (start_col_per_proc != 0) {
                west_proc = process_id - 1;
        }
        if (start_col_per_proc + cols_per_block != width) {
                east_proc = process_id + 1;
        }


        /* Get time before */
        start_time = MPI_Wtime();

        for (int i = 0; i < GENERATION; i++) {

                if (west_proc != -1) {

                        // MPI_Isend(&block_array_bef[1*(cols_per_block+2) + 1], 1, ColGreyType,  west, 0, MPI_COMM_WORLD, &send_west_req);
                        //
                        // MPI_Irecv(offset(src, 1, 0, cols+2), 1, ColGreyType,  west, 0, MPI_COMM_WORLD, &recv_west_req);
                        // MPI_Sendrecv(buffer, 10, MPI_INT, left, 123, buffer2, 10, MPI_INT, right, 123, MPI_COMM_WORLD, &status);

                        // MPI_Isend(offset(src, 1, 1, cols+2), 1, grey_col_type,  west, 0, MPI_COMM_WORLD, &send_west_req);
                        // MPI_Irecv(offset(src, 1, 0, cols+2), 1, grey_col_type,  west, 0, MPI_COMM_WORLD, &recv_west_req);

                        MPI_Send(&block_array_bef[offset(cols_per_block +2, 1, 1)], 1, ColGreyType, west_proc, 0, MPI_COMM_WORLD);
                        MPI_Recv(&block_array_bef[offset(cols_per_block +2, 1, 0)], 1, ColGreyType, west_proc, 0, MPI_COMM_WORLD, &status);
                }

                if (east_proc != -1) {
                        // MPI_Isend(offset(src, 1, cols, cols+2), 1, grey_col_type,  east, 0, MPI_COMM_WORLD, &send_east_req);
                        // MPI_Irecv(offset(src, 1, cols+1, cols+2), 1, grey_col_type,  east, 0, MPI_COMM_WORLD, &recv_east_req)

                        // MPI_Isend(offset(src, 1, cols, cols+2), 1, grey_col_type,  east, 0, MPI_COMM_WORLD, &send_east_req);
                        // MPI_Irecv(offset(src, 1, cols+1, cols+2), 1, grey_col_type,  east, 0, MPI_COMM_WORLD, &recv_east_req);

                        MPI_Recv(&block_array_bef[offset(cols_per_block +2, 1, cols_per_block + 1)], 1, ColGreyType, east_proc, 0, MPI_COMM_WORLD, &status);
                        MPI_Send(&block_array_bef[offset(cols_per_block +2, 1, cols_per_block)], 1, ColGreyType, east_proc, 0, MPI_COMM_WORLD);

                }


                if (north_proc != -1) {

                        // MPI_Isend(offset(src, 1, 1, cols+2), 1, grey_row_type, north, 0, MPI_COMM_WORLD, &send_north_req);
                        // MPI_Irecv(offset(src, 0, 1, cols+2), 1, grey_row_type, north, 0, MPI_COMM_WORLD, &recv_north_req);
                        MPI_Isend(&block_array_bef[offset(cols_per_block+2, 1,0)], 1, LineGreyType, north_proc, 0, MPI_COMM_WORLD, &send_row_n);
                        MPI_Irecv(&block_array_bef[offset(cols_per_block+2, 0,0)], 1, LineGreyType, north_proc, 0, MPI_COMM_WORLD, &recv_row_n);
                }


                if (south_proc != -1) {

                        // MPI_Isend(offset(src, rows, 1, cols+2), 1, grey_row_type, south, 0, MPI_COMM_WORLD, &send_south_req);
                        // MPI_Irecv(offset(src, rows+1, 1, cols+2), 1, grey_row_type, south, 0, MPI_COMM_WORLD, &recv_south_req);
                        MPI_Isend(&block_array_bef[offset(cols_per_block+2, rows_per_block, 0 )], 1, LineGreyType, south_proc, 0, MPI_COMM_WORLD, &send_row_s);
                        MPI_Irecv(&block_array_bef[offset(cols_per_block+2, rows_per_block+1, 0 )], 1, LineGreyType, south_proc, 0, MPI_COMM_WORLD, &recv_row_s);
                }


                // convoluteInner()
                convolute_grey(block_array_bef, block_array_after, h, 1, rows_per_block, 1, cols_per_block,cols_per_block,rows_per_block);



                flag1=0;
                flag2=0;
                flag11 =0;
                flag22=0;
                while (1) {
                        if (north_proc != -1) {
                                // MPI_Wait(&recv_row_n, &status);
                                MPI_Test(&recv_row_n, &flag1, &status);
                                // printf("flag1 = %d\n",flag1 );
                                if (flag1 == 1) {
                                        // convolute_grey(up);
                                        convolute_grey(block_array_bef, block_array_after, h, 1, 1, 1, cols_per_block, cols_per_block,rows_per_block);

                                        flag11=1;
                                }
                        }
                        else{
                                flag11 =1;
                        }
                        if (south_proc != -1) {
                                // MPI_Wait(&recv_row_s, &status);
                                MPI_Test(&recv_row_s, &flag2, &status);
                                // printf("flag2 = %d\n",flag2 );
                                if (flag2 == 1) {
                                        // convolute_grey(down)
                                        convolute_grey(block_array_bef, block_array_after, h, rows_per_block, rows_per_block, 1, cols_per_block, cols_per_block,rows_per_block);
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
                        convolute_grey(block_array_bef, block_array_after, h, 1, rows_per_block, 1, 1, cols_per_block,rows_per_block);
                }
                if (east_proc != -1) {
                        convolute_grey(block_array_bef, block_array_after, h, 1, rows_per_block, cols_per_block, cols_per_block, cols_per_block,rows_per_block);
                }
                // if (north_proc != -1) {
                //         MPI_Wait(&recv_row_n, &status);
                //         convolute_grey(block_array_bef, block_array_after, h, 1, 1, 1, cols_per_block, cols_per_block,rows_per_block);
                // }
                // if (south_proc != -1) {
                //         MPI_Wait(&recv_row_s, &status);
                //         convolute_grey(block_array_bef, block_array_after, h, rows_per_block, rows_per_block, 1, cols_per_block, cols_per_block,rows_per_block);
                // }


                /* Wait to have sent all borders */
                if (north_proc != -1)
                        MPI_Wait(&send_row_n, &status);
                if (south_proc != -1)
                        MPI_Wait(&send_row_s, &status);


                swap_arrays(&block_array_bef, &block_array_after);
                // unsigned char* tmp=NULL;
                // tmp = block_array_bef;
                // block_array_bef = block_array_after;
                // block_array_after = tmp;

                /*check if image changed */

                // changed =0;
                // changed = ImageChanged(block_array_bef,block_array_after);
                // if (ImageChanged) {
                //     printf("GENERATION : %d Image Changed = %d \n",i,changed );
                // }
                // else{
                //     printf("GENERATION : %d Image NOT Changed = %d \n",i,changed );
                // }


        }

        /* Get time elapsed */
        end_time = MPI_Wtime();
        elapsed_time = end_time - start_time;


        /* Parallel write */
      /*  char *outImage = malloc((strlen(image_name) + 8) * sizeof(char));
        strcpy(outImage, "fltred_");*/

      //  strcat(outImage, image_name);
        MPI_File outFile;
        error = MPI_File_open(MPI_COMM_WORLD, outImage, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &outFile);
        if (error) {
                printf( "Unable to open the output file..No space! \n" ); fflush(stdout);
        }
        for (int row = 0; row < rows_per_block; row++) {

                MPI_File_seek(outFile, (start_row_per_proc + row) * width + start_col_per_proc, MPI_SEEK_SET);

                MPI_File_write(outFile, &block_array_bef[offset(cols_per_block+2, row+1, 1)], cols_per_block, MPI_BYTE, MPI_STATUS_IGNORE);
        }

        MPI_File_close(&outFile);


        MPI_Reduce(&elapsed_time, &elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if(process_id == 0)
          printf("Execution time: %f\n",elapsed);



        free(block_array_bef);
        free(block_array_after);
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


inline void convolute_grey(unsigned char* A, unsigned char* B, float** h, int row_start, int row_end, int col_start, int col_end, int width /*cols*/, int height /*rows*/) {

        int ii,jj;
        float val=0;

#pragma omp parallel for shared(A, B) schedule(static) collapse(4)

        for (int i = row_start; i <= row_end; ++i)          // rows [1,rows_per_block] correct
        {
                for (int j = col_start; j <= col_end; ++j)  // columns [1,cols_per_block] correct
                {
                        //val = 0;
                        for (int m = 0; m < 3; ++m) // kernel rows [0-3] correct
                        {

                                for (int n = 0; n < 3; ++n) // kernel columns [0-3] correct
                                {
                                        // index of input signal, used for checking boundary
                                        ii = i + (m - K_CENTER_Y); // i +m -1 ----> curr_row_block + curr_row_kernel -1

                                        jj = j + (n - K_CENTER_X); // j +n -1 ----> curr_col_block + curr_col_kernel -1

                                        // ignore input samples which are out of bound
                                        // if (ii >= 0 && ii <= height && jj >= 0 && jj <= width+2) {
                                        val+=A[(width+2)*ii + jj] * h[m][n];
                                        // }
                                }
                        }
                        B[(width+2)*i + j] = val; // P[i][j] += N[ii][jj] * M[m][n];
                        val=0;
                }
        }
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

inline int offset(int k,int i,int offs ){ //k = cols or rows
        return k*i + offs;
}

inline int ImageChanged(unsigned char* A,unsigned char* B ){
        return strcmp(A, B);
}


inline int div2blocks(int height, int width, int num_processes, int* rows_per_block, int* cols_per_block ){

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
        *rows_per_block = height/int_root_num_processes;
        *cols_per_block = width/int_root_num_processes;
        // printf("division() : rows_per_block = %d cols_per_block = %d\n",*rows_per_block,*cols_per_block );
        return int_root_num_processes;

}
