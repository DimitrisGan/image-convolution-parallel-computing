#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#define GENERATION 2

// find center position of kernel
#define K_CENTER_X 1
#define K_CENTER_Y 1

#define KWIDTH 3
#define KHEIGHT 3

int div2blocks(int, int, int, int*, int* );
int Argms_handler(int argc,char** argv, int* image_type,int* width,int* height);

int offset(int Cols,int i,int offs );
void swap_arrays(char** A, char** B);
int ImageChanged(char* A,char* B );

int convolute_grey(char* A, char* B, float** h, int row_start, int row_end, int col_start, int col_end, int width, int height);

// void convolute(char *, char *, int, int, int, int, int, int, float**);
// void convolute_grey(char*, char*, int, int, int, int, float**);

int main(int argc, char** argv) {



        int width=0;
        int height=0;
        char* temp_array=NULL;
        int number_offs;

        char* image_name = malloc((strlen(argv[1])+1) * sizeof(char));
        strcpy(image_name, argv[1]);
        int image_type;

        int rows_per_block = -1;
        int cols_per_block = -1;

        int rooted_num_procs;

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
                printf("ERROR!!!\n");
                exit(1);
        }


        if (process_id == 0) {

                rooted_num_procs = div2blocks(height, width, num_processes, &rows_per_block, &cols_per_block);
                printf("rows_per_block = %d | cols_per_block = %d \n\n",rows_per_block,cols_per_block );

        }

        MPI_Bcast(&rooted_num_procs, 1, MPI_INT, 0, MPI_COMM_WORLD);

        MPI_Bcast(&rows_per_block, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&cols_per_block, 1, MPI_INT, 0, MPI_COMM_WORLD);




        if (process_id != 0) {

                printf("process_id= %d \n ",process_id );
                //gia arxh to kanw gia grey iamge            printf("div2blocks rows (height = %d , width = %d , num_processes = %d ) \n",height, width, num_processes );
                printf("rows_per_block = %d | cols_per_block = %d \n\n",rows_per_block,cols_per_block );
        }


        // filtes kernel
        int emboss_filter[3][3]= {{-2,-1,0},{-1,1,1},{0,1,2}};
        int sharpen_filter[3][3]= {{0,-1,0},{-1,5,-1},{0,-1,0}};
        int gaussian_blur[3][3] = {{1, 2, 1}, {2, 4, 2}, {1, 2, 1}};


        float **h = malloc(sizeof(float*));
        for (int i = 0; i < 3; i++) {
                h[i] =malloc(3* sizeof(float));
                for (int j = 0; j < 3; j++) {

                        // h[i][j] = emboss_filter[i][j] / 1.0;
                        // h[i][j] = sharpen_filter[i][j] / 1.0;
                        h[i][j] = gaussian_blur[i][j] / 16.0;

                }
        }




        int start_row_per_proc = (process_id / rooted_num_procs) *rows_per_block;
        int start_col_per_proc = (process_id % rooted_num_procs) * cols_per_block;

        // printf(" start_row = %d per proc : %d \n",start_row_per_proc,process_id );
        // printf(" start_col = %d per proc : %d \n",start_col_per_proc,process_id );


        char* block_array_bef   = calloc((rows_per_block + 2)*(cols_per_block + 2), sizeof(char));
        char* block_array_after = calloc((rows_per_block + 2)*(cols_per_block + 2), sizeof(char));

        // char* block_array_bef   = malloc((rows_per_block + 2)*(cols_per_block + 2)* sizeof(char));
        // char* block_array_after = malloc((rows_per_block + 2)*(cols_per_block + 2)* sizeof(char));




        //create datatypes for line and column FOR RGB AND GREY;
        //**ATTENTION** keep in mind that LineType is of size cols+2
        MPI_Type_contiguous(cols_per_block+2, MPI_BYTE, &LineGreyType);
        MPI_Type_commit(&LineGreyType);

        MPI_Type_vector(cols_per_block, 1, rows_per_block, MPI_BYTE, &ColGreyType);
        MPI_Type_commit(&ColGreyType);


        MPI_Type_contiguous(3*cols_per_block+2, MPI_BYTE, &LineRgbType);
        MPI_Type_commit(&LineRgbType);

        MPI_Type_vector(3*cols_per_block, 1, rows_per_block, MPI_BYTE, &ColRgbType);
        MPI_Type_commit(&ColRgbType);



        printf("image_name = %s\n",image_name );

        int error =MPI_File_open(MPI_COMM_WORLD, image_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
        if (error) {
                printf( "Unable to open file \"temp\"\n" ); fflush(stdout);

        }



        if (image_type == 0) {
                for (int row = 0; row < rows_per_block; row++) {
                        // gia kathe grammh pernoume th grammh ths
                        MPI_File_seek( fh, (start_row_per_proc + row ) * width + start_col_per_proc, MPI_SEEK_SET );
                        number_offs = offset(cols_per_block +2,row+1,1);
                        temp_array = &block_array_bef[number_offs];
                        // printf("ROW[%d]=FILE READ for proc_id[%d]: offset(cols_per_block +2,row+1,1) = %d\n",row,process_id,offset(cols_per_block +2,row+1,1) );
                        MPI_File_read(fh, temp_array, cols_per_block, MPI_BYTE, &status);

                        // MPI_File_set_view

                }
        }

        MPI_File_close(&fh);



        MPI_Request send_row_n, recv_row_n;
        MPI_Request send_col_w, recv_col_w;
        MPI_Request send_row_s, recv_row_s;
        MPI_Request send_col_e, recv_col_e;

        MPI_Request send_nw, recv_nw;
        MPI_Request send_ne, recv_ne;
        MPI_Request send_sw, recv_sw;
        MPI_Request send_se, recv_se;

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
                south_proc = process_id + rooted_num_procs;
        }
        if (start_col_per_proc != 0) {
                west_proc = process_id - 1;
        }
        if (start_col_per_proc + cols_per_block != width) {
                east_proc = process_id + 1;
        }

        for (int i = 0; i < GENERATION; i++) {

                if (west_proc != -1) {

                        // MPI_Isend(&block_array_bef[1*(cols_per_block+2) + 1], 1, ColGreyType,  west, 0, MPI_COMM_WORLD, &send_west_req);
                        //
                        // MPI_Irecv(offset(src, 1, 0, cols+2), 1, ColGreyType,  west, 0, MPI_COMM_WORLD, &recv_west_req);
                        // MPI_Sendrecv(buffer, 10, MPI_INT, left, 123, buffer2, 10, MPI_INT, right, 123, MPI_COMM_WORLD, &status);

                        MPI_Send(&block_array_bef[offset(cols_per_block +2, 1, 1)], 1, ColGreyType, west_proc, 0, MPI_COMM_WORLD);
                        MPI_Recv(&block_array_bef[offset(cols_per_block +2, 1, 0)], 1, ColGreyType, west_proc, 0, MPI_COMM_WORLD, &status);
                }


                if (east_proc != -1) {
                        // MPI_Isend(offset(src, 1, cols, cols+2), 1, grey_col_type,  east, 0, MPI_COMM_WORLD, &send_east_req);
                        // MPI_Irecv(offset(src, 1, cols+1, cols+2), 1, grey_col_type,  east, 0, MPI_COMM_WORLD, &recv_east_req)

                        MPI_Recv(&block_array_bef[offset(cols_per_block +2, 1, cols_per_block + 1)], 1, ColGreyType, east_proc, 0, MPI_COMM_WORLD, &status);
                        MPI_Send(&block_array_bef[offset(cols_per_block +2, 1, cols_per_block)], 1, ColGreyType, east_proc, 0, MPI_COMM_WORLD);

                }


                if (north_proc != -1) {

                        MPI_Isend(&block_array_bef[offset(cols_per_block+2, 1,0)], 1, LineGreyType, north_proc, 0, MPI_COMM_WORLD, &send_row_n);
                        MPI_Irecv(&block_array_bef[offset(cols_per_block+2, 0,0)], 1, LineGreyType, north_proc, 0, MPI_COMM_WORLD, &recv_row_n);
                }


                if (south_proc != -1) {

                        MPI_Isend(&block_array_bef[offset(rows_per_block+2, cols_per_block, 0 )], 1, LineGreyType, south_proc, 0, MPI_COMM_WORLD, &send_row_s);
                        MPI_Irecv(&block_array_bef[offset(rows_per_block+2, cols_per_block+1, 0 )], 1, LineGreyType, south_proc, 0, MPI_COMM_WORLD, &recv_row_s);
                }


                // convoluteInner()
                // int convolute_grey(char* A, char* B, float** h, int row_start, int row_end, int col_start, int col_end, int width, int height) {

                convolute_grey(block_array_bef, block_array_after, h, 2, rows_per_block-1, 2, cols_per_block-1,cols_per_block+2,rows_per_block);
                /* Inner Data Convolute */

                // convolute(block_array_bef, block_array_after, 1, rows_per_block, 1, cols_per_block, cols_per_block, rows_per_block, h);

                if (north_proc != -1) {
                        MPI_Wait(&recv_row_n, &status);
                        // convolute(up)
                }
                if (south_proc != -1) {
                        MPI_Wait(&recv_row_s, &status);
                        //convolute(down)
                }



                /* Wait to have sent all borders */
                if (north_proc != -1)
                        MPI_Wait(&send_row_n, &status);
                if (south_proc != -1)
                        MPI_Wait(&send_row_s, &status);



                /* swap arrays */
                // swap_arrays(&block_array_bef, &block_array_after);

                char* tmp = block_array_bef;
                block_array_bef = block_array_after;
                block_array_after = tmp;

                /*check if image changed */
                // int changed=0;
                // changed = ImageChanged(block_array_bef,block_array_after);
                // printf("ImageChanged = %d in process_id = %d  in GENERATION  %d \n",changed,process_id,i );


        }




        /* Parallel write */
        // char *outImage = malloc((strlen(image_name) +9) * sizeof(char));
        // strcpy(outImage, "emboss_");
        char *outImage = malloc((strlen(image_name) + 7) * sizeof(char));
        strcpy(outImage, "blur_");

        strcat(outImage, image_name);
        MPI_File outFile;
        error = MPI_File_open(MPI_COMM_WORLD, outImage, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &outFile);
        if (error) {
                printf( "Unable to open file \"temp\"\n" ); fflush(stdout);
        }
        for (int row = 0; row < rows_per_block; row++) {

                MPI_File_seek(outFile, (start_row_per_proc + row) * width + start_col_per_proc, MPI_SEEK_SET);
                number_offs = offset(cols_per_block+2,row+1,1);
                temp_array = &block_array_bef[number_offs];
                // tmpbuf = offset(src, i, 1, cols+2);
                // MPI_File_write(outFile, tmpbuf, cols, MPI_BYTE, MPI_STATUS_IGNORE);
                MPI_File_write(outFile, temp_array, cols_per_block, MPI_BYTE, MPI_STATUS_IGNORE);
        }


        // for (int row = 0; row < rows_per_block; row++) {
        //     // gia kathe grammh pernoume th grammh ths
        //     MPI_File_seek( fh, (start_row_per_proc + row ) * width + start_col_per_proc, MPI_SEEK_SET );
        //
        //     MPI_File_read(fh, &block_array_bef[offset(cols_per_block +2,row+1,1)], cols_per_block, MPI_BYTE, &status);
        //
        //     // MPI_File_set_view
        //
        // }



        MPI_File_close(&outFile);

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
   ██████   ██████  ██ ███    ██  ██████
   ██   ██ ██    ██ ██ ████   ██ ██
   ██   ██ ██    ██ ██ ██ ██  ██ ██   ███
   ██   ██ ██    ██ ██ ██  ██ ██ ██    ██
   ██████   ██████  ██ ██   ████  ██████
 */




inline int ImageChanged(char* A,char* B ){
        return strcmp(A, B);
}

// void convolute(char *src, char *dst, int row_from, int row_to, int col_from, int col_to, int width, int height, float** h ) {
//  int i, j;
//  // printf("convolute : WIDTH = %d AND HEIGHT = %d \n", width,height );
//
//      for (i = row_from ; i <= row_to ; i++)
//          for (j = col_from ; j <= col_to ; j++)
//              convolute_grey(src, dst, i, j, width+2, height, h);
//
// }
//
// void convolute_grey(char *src, char *dst, int x, int y, int width, int height, float** h) {
//  int i, j, k, l;
//  // printf("convolute_grey : WIDTH = %d AND HEIGHT = %d \n", width,height );
//  float val = 0;
//  for (i = x-1, k = 0 ; i <= x+1 ; i++, k++)
//      for (j = y-1, l = 0 ; j <= y+1 ; j++, l++)
//          val += src[width * i + j] * h[k][l];
//  dst[width * x + y] = val;
// }



// int convolution(char** A ,char** B ,int** h ,int width , int height);

//convolution_rgb()
// void convolute(uint8_t *src, uint8_t *dst, int row_from, int row_to, int col_from, int col_to, int width, int height, float** h, color_t imageType) {

int convolute_grey(char* A, char* B, float** h, int row_start, int row_end, int col_start, int col_end, int width, int height) {

        int ii,jj;


        for (int i = row_start; i <= row_end; ++i)          // rows [1,rows_per_block] correct
        {
                for (int j = col_start; j <= col_end; ++j)  // columns [1,cols_per_block] correct
                {
                        for (int m = 0; m < 3; ++m) // kernel rows [0-3] correct
                        {

                                for (int n = 0; n < 3; ++n) // kernel columns [0-3] correct
                                {
                                        // index of input signal, used for checking boundary
                                        ii = i + (m - K_CENTER_Y); // i +m -1 ----> curr_row_block + curr_row_kernel -1

                                        jj = j + (n - K_CENTER_X); // j +n -1 ----> curr_col_block + curr_col_kernel -1

                                        // ignore input samples which are out of bound
                                        if (ii >= 0 && ii < height && jj >= 0 && jj < width)
                                                B[width*i + j] += A[width*ii + jj] * h[m][n];
                                        // P[i][j] += N[ii][jj] * M[m][n];
                                }
                        }
                        // P[i][j] += N[ii][jj] * M[m][n];

                }
        }

        return 0;
}


void swap_arrays(char** A, char** B){
        char *temp;
        temp = *A;
        *A = *B;
        *B = temp;
        return;
}



/*
   ██████  ███████ ██████  ███████ ███████  ██████ ████████
   ██   ██ ██      ██   ██ ██      ██      ██         ██
   ██████  █████   ██████  █████   █████   ██         ██
   ██      ██      ██   ██ ██      ██      ██         ██
   ██      ███████ ██   ██ ██      ███████  ██████    ██
 */


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
        // for (size_t i = 0; i < argc; i++) {
        //     printf("argv[%ld] = %s\n",i,argv[i] );
        //     printf("argc = %d\n",argc );
        // }
        return 0;
}

int offset(int k,int i,int offs ){ //k = cols or rows
        return k*i + offs;
}




inline int div2blocks(int height, int width, int num_processes, int* rows_per_block, int* cols_per_block ){

        double root_num_processes = sqrt(num_processes);

        if ((root_num_processes - floor(root_num_processes)) > 0 ) {
                printf("1:DEN DIAIREITAI ME RIZA \n" );
                exit(1);
        }
        int int_root_num_processes = root_num_processes;
        if ((height % int_root_num_processes) && (width % int_root_num_processes)   ) {
                printf("2:DEN DIAIROUNTAI OMORFA TA ROWS KAI COLS \n" );
                exit(1);
        }
        printf("height = %d\n",height );
        *rows_per_block = height/int_root_num_processes;
        *cols_per_block = width/int_root_num_processes;
        // printf("division() : rows_per_block = %d cols_per_block = %d\n",*rows_per_block,*cols_per_block );
        return int_root_num_processes;

}
