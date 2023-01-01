
#include "sigfish.h"
#include "cdtw.h"
#include "cudtw.cuh"
#include "assert.h"

#define MALLOC_CHK(ret)                           \
    {                                             \
        if ((ret) == NULL)                        \
        {                                         \
            printf("Could not allocate memory."); \
        }                                         \
    }

#define CHECK(x)                              \
    do                                        \
    {                                         \
        cudaError_t err = (x);                \
        if (err != cudaSuccess)               \
        {                                     \
            fprintf(stderr, "API error %s              \
                    : % d Returned                      \
                    : % d\n",                 \
                    __FILE__, __LINE__, err); \
            exit(1);                          \
        }                                     \
    } while (0)

/* check whether the last CUDA function or CUDA kernel launch is erroneous and if yes an error message will be printed
and then the program will be aborted*/
#define CUDA_CHK()                      \
    {                                   \
        gpu_assert(__FILE__, __LINE__); \
    }

// enum sigfish_log_level_opt _log_level = LOG_VERB;

static inline void gpu_assert(const char *file, uint64_t line)
{
    cudaError_t code = cudaGetLastError();
    if (code != cudaSuccess)
    {
        fprintf(stderr, "[%s::ERROR]\033[1;31m Cuda error: %s \n in file : %s line number : %lu\033[0m\n",
                __func__, cudaGetErrorString(code), file, line);
        if (code == cudaErrorLaunchTimeout)
        {
            fprintf(stderr, "%s", "The kernel timed out. You have to first disable the cuda "
                                  "time out.");
            fprintf(
                stderr,
                "On Ubuntu do the following\nOpen the file /etc/X11/xorg.conf\nYou "
                "will have a section about your NVIDIA device. Add the following "
                "line to it.\nOption \"Interactive\" \"0\"\nIf you do not have a "
                "section about your NVIDIA device in /etc/X11/xorg.conf or you do "
                "not have a file named /etc/X11/xorg.conf, run the command sudo "
                "nvidia-xconfig to generate a xorg.conf file and do as above.\n\n");
        }
        exit(-1);
    }
}

__device__ aln_t *init_aln2(aln_t *cuda_aln){
    float score = INFINITY;
    float score2 = INFINITY;
    int32_t pos = -1;
    int32_t rid = -1;
    char d = 0;
    for (int l = 0; l < SECONDARY_CAP; l++){
        aln_t tmp = {rid, pos, pos, score, score2, d, 0};
        cuda_aln[l] = tmp;
    }
    return cuda_aln;
}

__device__ void update_aln2(aln_t *aln, float score, int32_t rid, int32_t pos, char d, float *cost, int32_t qlen, int32_t rlen){
    int l = 0;
    for (; l < SECONDARY_CAP; l++){
        if (score > aln[l].score)
        {
            break;
        }
        else
        {
            continue;
        }
    }

    if (l != 0){
        for (int m = 0; m < l - 1; m++)
        {
            aln[m] = aln[m + 1];
        }
        aln[l - 1].score = score;
        aln[l - 1].pos_end = pos;
        aln[l - 1].rid = rid;
        aln[l - 1].d = d;

        aln[l - 1].pos_st = pos - qlen + 1;
    }
}

__device__ float min3(float a, float b, float c){
    float min;

    min = a;
    if (b < min)
        min = b;
    if (c < min)
        min = c;
    return min;
}

__device__ void dtw_subsequence(float *x, float *y, int n, int m, float *cost){

    cost[0] = fabs(x[0] - y[0]);


    for (int i = 1; i < n; i++){
        cost[i * m] = fabs(x[i] - y[0]) + cost[(i - 1) * m];
    }
         

    for (int j = 1; j < m; j++){ // subsequence variation: D(0,j) := c(x0, yj)
        cost[j] = fabs(x[0] - y[j]); 
    }
         

    for (int i = 1; i < n; i++){
        for (int j = 1; j < m; j++){
            cost[i * m + j] = fabs(x[i] - y[j]) + min3(cost[(i - 1) * m + j], cost[(i - 1) * m + (j - 1)], cost[i * m + (j - 1)]);
        }  
    }
      
}

__global__ void dtw_single2(core_t *core, uint64_t *cuda_len_raw_signal, size_t *cuda_et_n, 
                    int64_t *cuda_qstart, int64_t *cuda_qend, float *cuda_query, int32_t *cuda_rlen, int max_qlen, 
                    int max_rlen, float *cuda_cost, float *cuda_core_ref_f, float *cuda_db_aln_score, float *cuda_db_aln_score2, 
                    int32_t *cuda_db_aln_pos_st, int32_t *cuda_db_aln_pos_end, int32_t *cuda_db_aln_rid, char *cuda_db_aln_d, 
                    uint8_t *cuda_db_aln_mapq, aln_t *cuda_aln, int32_t SIZE_NUM_REF, int32_t *cuda_ref_st_offset, int8_t core_opt_flag)

                    //                     int32_t *cuda_aln_rid, int32_t *cuda_aln_pos_st, int32_t *cuda_aln_pos_end, float *cuda_aln_score, float *cuda_aln_score2, 
                    // uint8_t *cuda_aln_mapq, char *cuda_aln_d
{
    int i = threadIdx.x;

    if (cuda_len_raw_signal[i] > 0 && cuda_et_n[i] > 0){ // some checks to see if a good read

        // aln_t tempaln = {cuda_aln_rid[i], cuda_aln_pos_st[i], cuda_aln_pos_end[i], cuda_aln_score[i], cuda_aln_score2[i], cuda_aln_mapq[i], cuda_aln_d[i]};
    
        // aln_t *aln = init_aln2(cuda_aln + (SECONDARY_CAP * i)); // initialise a alignment struct

        aln_t *aln = init_aln2(cuda_aln); // initialise a alignment struct

        int64_t start_idx = cuda_qstart[i]; // starting index of the query
        int64_t end_idx = cuda_qend[i];     // ending index of the query

        int32_t qlen = end_idx - start_idx; // query chunk length

        int8_t rna  = core_opt_flag & SIGFISH_RNA; // if data is RNA                                       

        for (int j = 0; j < SIZE_NUM_REF; j++){

            int32_t rlen = cuda_rlen[j];

            dtw_subsequence((cuda_query + (max_qlen * i)) , (&cuda_core_ref_f[j] + (i)), qlen, rlen, (cuda_cost + (max_qlen * max_rlen *i))); /// check correct values for i and j passed

             for (int k = (qlen - 1) * rlen; k < qlen * rlen; k += qlen){
                 float min_score = INFINITY;
                 int32_t min_pos = -1;
                for (int m = 0; m < qlen && k + m < qlen * rlen; m++)
                 {
                     if (cuda_cost[(k + m) + ( max_qlen * max_rlen * j)] < min_score)
                     {
                         min_score = cuda_cost[(k + m) + ( max_qlen * max_rlen * j)];
                         min_pos = m + k;
                     }
                 }
                 update_aln2(aln, min_score, j, min_pos - (qlen - 1) * rlen, '+', (cuda_cost +  max_qlen * max_rlen * j), qlen, rlen);
             }

            //  if (!rna){ // if DNA we must consider the reverse strand as well
            //      dtw_subsequence((cuda_query + ( max_qlen * i)), core->ref->reverse[j], qlen, rlen, (cuda_cost +  max_qlen * max_rlen * j));

            //      for (int k = (qlen - 1) * rlen; k < qlen * rlen; k += qlen)
            //      {
            //          float min_score = INFINITY;
            //          int32_t min_pos = -1;
            //          for (int m = 0; m < qlen && k + m < qlen * rlen; m++)
            //         {
            //              if (cuda_cost[(k + m) + ( max_qlen * max_rlen * j)] < min_score)
            //              {
            //                  min_score = cuda_cost[(k + m) + ( max_qlen * max_rlen * j)];
            //                  min_pos = m + k;
            //              }
            //          }
            //          update_aln2(aln, min_score, j, min_pos - (qlen - 1) * rlen, '-', (cuda_cost +  max_qlen * max_rlen * j), qlen, rlen);
            //      }
            // }
        }

        cuda_db_aln_score[i] = aln[SECONDARY_CAP - 1].score;
        cuda_db_aln_score2[i] = aln[SECONDARY_CAP - 2].score;
        cuda_db_aln_pos_st[i] = aln[SECONDARY_CAP - 1].d == '+' ? aln[SECONDARY_CAP - 1].pos_st : cuda_rlen[aln[SECONDARY_CAP - 1].rid] - aln[SECONDARY_CAP - 1].pos_end;
        cuda_db_aln_pos_end[i] = aln[SECONDARY_CAP - 1].d == '+' ? aln[SECONDARY_CAP - 1].pos_end : cuda_rlen[aln[SECONDARY_CAP - 1].rid] - aln[SECONDARY_CAP - 1].pos_st;

        cuda_db_aln_pos_st[i] += cuda_ref_st_offset[aln[SECONDARY_CAP - 1].rid];
        cuda_db_aln_pos_end[i] += cuda_ref_st_offset[aln[SECONDARY_CAP - 1].rid];
        cuda_db_aln_rid[i] = aln[SECONDARY_CAP - 1].rid;
        cuda_db_aln_d[i] = aln[SECONDARY_CAP - 1].d;

        int mapq = (int)round(500 * (cuda_db_aln_score2[i] - cuda_db_aln_score[i]) / cuda_db_aln_score[i]);
        if (mapq > 60)
        {
            mapq = 60;
        }
        cuda_db_aln_mapq[i] = mapq;

        free(aln);
    }
}

int max(int arr[], int SIZE) {

    int max = arr[0];
  
    for (int i = 1; i < SIZE; i++)
        if (arr[i] > max)
            max = arr[i];
  
    return max;
}

void dtw_cuda_db(core_t *core, db_t *db){
    // For error checking
    cudaError_t code;

    int32_t SIZE = db->n_rec;
    int32_t SIZE_NUM_REF = core->ref->num_ref;

    //  Pointers for elements to save in cuda
    core_t *cuda_core;
    db_t *cuda_db;
    uint64_t *cuda_len_raw_signal;
    size_t *cuda_et_n;
    int64_t *cuda_qstart;
    int64_t *cuda_qend;
    float *cuda_db_aln_score;
    float *cuda_db_aln_score2;
    int32_t *cuda_db_aln_pos_st;
    int32_t *cuda_db_aln_pos_end;
    int32_t *cuda_db_aln_rid;
    char *cuda_db_aln_d;
    uint8_t *cuda_db_aln_mapq;
    float *cuda_query;
    float *cuda_cost;
    int32_t *cuda_rlen;
    float *cuda_core_ref_f;
    int32_t *cuda_core_ref_len;
    int32_t *cuda_ref_st_offset;

    aln_t *cuda_aln;

    // Pointers for aln structure
    int32_t cuda_aln_rid;
    int32_t cuda_aln_pos_st;
    int32_t cuda_aln_pos_end;
    float cuda_aln_score;
    float cuda_aln_score2;
    uint8_t cuda_aln_mapq;
    char cuda_aln_d;

    // Arrays
    uint64_t len_raw_signal[SIZE];
    size_t et_n[SIZE]; //long unsigned int
    int64_t qstart[SIZE];
    int64_t qend[SIZE];
    int32_t qlen[SIZE];
    int32_t rlen[SIZE_NUM_REF];
    float core_ref_f[SIZE_NUM_REF];
    float *core_ref_rev[SIZE_NUM_REF];
    float db_aln_score[SIZE];
    float db_aln_score2[SIZE];
    int32_t db_aln_pos_st[SIZE];
    int32_t db_aln_pos_end[SIZE];
    int32_t db_aln_rid[SIZE];
    char db_aln_d[SIZE];
    uint8_t db_aln_maqq[SIZE];
    int32_t ref_st_offset[SIZE_NUM_REF];

    //variables
    int8_t core_opt_flag = core->opt.flag;

    // Assign values to Arrays
    for(int i=0; i<SIZE; i++){
        len_raw_signal[i] = db->slow5_rec[i]->len_raw_signal;
        et_n[i] = db->et[i].n;
        qstart[i] = db->qstart[i];
        qend[i] = db->qend[i];
        qlen[i] = qend[i] - qstart[i];
    }

    for(int i=0; i<SIZE_NUM_REF; i++){
        rlen[i] = core->ref->ref_lengths[i];
        core_ref_f[i] = *core->ref->forward[i];  // Check again pointer ???
        // core_ref_rev[i] = core->ref->reverse[i]; // Why this gives segmentation fault??????
        ref_st_offset[i] = core->ref->ref_st_offset[i]; // Check this if Error

        // fprintf(stderr, "core_ref_f %i = %f\n",i, core_ref_f[i]);
    }

    int max_qlen = max(qlen, SIZE);
    int max_rlen = max(rlen, SIZE_NUM_REF);

    float query[SIZE][max_qlen] = {};

    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < qlen[i]; j++) {
            if (!(core->opt.flag & SIGFISH_INV) && core_opt_flag){
                query[i][qlen[i] - 1 - j] = db->et[i].event[j + qstart[i]].mean;
            } else {
                query[i][j] = db->et[i].event[j + qstart[i]].mean;
            }
        }
    }

    // Allocate memory in GPU
    cudaMalloc((void **)&cuda_db, sizeof(db_t));
    CUDA_CHK();
    cudaMalloc((void **)&cuda_core, sizeof(core_t));
    CUDA_CHK();
    cudaMalloc((void **)&cuda_len_raw_signal, sizeof(uint64_t)*SIZE);
    CUDA_CHK();
    cudaMalloc((void **)&cuda_et_n, sizeof(size_t)*SIZE);
    CUDA_CHK();
    cudaMalloc((void **)&cuda_qstart, sizeof(int64_t)*SIZE);
    CUDA_CHK();
    cudaMalloc((void **)&cuda_qend, sizeof(int64_t)*SIZE);
    CUDA_CHK();
    cudaMalloc((void **)&cuda_query, sizeof(float)*SIZE*max_qlen);
    CUDA_CHK();
    cudaMalloc((void **)&cuda_cost, sizeof(float)*SIZE*max_qlen*max_rlen);
    CUDA_CHK();
    cudaMalloc((void **)&cuda_rlen, sizeof(int32_t)*SIZE_NUM_REF);
    CUDA_CHK();
    cudaMalloc((void **)&cuda_core_ref_f, sizeof(float)*SIZE_NUM_REF);
    CUDA_CHK();
    cudaMalloc((void **)&cuda_db_aln_pos_st, sizeof(int32_t)*SIZE);
    CUDA_CHK();
    cudaMalloc((void **)&cuda_db_aln_pos_end, sizeof(int32_t)*SIZE);
    CUDA_CHK();
    cudaMalloc((void **)&cuda_db_aln_rid, sizeof(int32_t)*SIZE);
    CUDA_CHK();
    cudaMalloc((void **)&cuda_db_aln_d, sizeof(char)*SIZE);
    CUDA_CHK();
    cudaMalloc((void **)&cuda_db_aln_score, sizeof(float)*SIZE);
    CUDA_CHK();
    cudaMalloc((void **)&cuda_db_aln_score2, sizeof(float)*SIZE);
    CUDA_CHK();
    cudaMalloc((void **)&cuda_db_aln_mapq, sizeof(uint8_t)*SIZE);
    CUDA_CHK();
    cudaMalloc((void **)&cuda_ref_st_offset, sizeof(int32_t)*SIZE);
    CUDA_CHK();
    cudaMalloc((void **)&cuda_aln, sizeof(aln_t)*SECONDARY_CAP*SIZE);
    CUDA_CHK();

    cudaMalloc((void **)&cuda_aln_rid, sizeof(int32_t)*SECONDARY_CAP*SIZE);
    CUDA_CHK();
    cudaMalloc((void **)&cuda_aln_pos_st, sizeof(int32_t)*SECONDARY_CAP*SIZE);
    CUDA_CHK();
    cudaMalloc((void **)&cuda_aln_pos_end, sizeof(int32_t)*SECONDARY_CAP*SIZE);
    CUDA_CHK();
    cudaMalloc((void **)&cuda_aln_score, sizeof(float)*SECONDARY_CAP*SIZE);
    CUDA_CHK();
    cudaMalloc((void **)&cuda_aln_score2, sizeof(float)*SECONDARY_CAP*SIZE);
    CUDA_CHK();
    cudaMalloc((void **)&cuda_aln_mapq, sizeof(uint8_t)*SECONDARY_CAP)*SIZE;
    CUDA_CHK();
    cudaMalloc((void **)&cuda_aln_d, sizeof(char)*SECONDARY_CAP*SIZE);
    CUDA_CHK();

    // copy content from main memory to cuda memory
    cudaMemcpy(cuda_db, db, sizeof(db_t), cudaMemcpyHostToDevice);
    CUDA_CHK();
    // cudaMemcpy(cuda_core, core, sizeof(core_t), cudaMemcpyHostToDevice);
    // CUDA_CHK();
    cudaMemcpy(cuda_len_raw_signal, len_raw_signal, sizeof(uint64_t)*SIZE, cudaMemcpyHostToDevice);
    CUDA_CHK();
    cudaMemcpy(cuda_et_n, et_n, sizeof(size_t)*SIZE, cudaMemcpyHostToDevice);
    CUDA_CHK();
    cudaMemcpy(cuda_qstart, qstart, sizeof(int64_t)*SIZE, cudaMemcpyHostToDevice);
    CUDA_CHK();
    cudaMemcpy(cuda_qend, qend, sizeof(int64_t)*SIZE, cudaMemcpyHostToDevice);
    CUDA_CHK();
    cudaMemcpy(cuda_query, query, sizeof(float)*SIZE*max_qlen, cudaMemcpyHostToDevice);
    CUDA_CHK();
    cudaMemcpy(cuda_rlen, rlen, sizeof(int32_t)*SIZE_NUM_REF, cudaMemcpyHostToDevice);
    CUDA_CHK();
    cudaMemcpy(cuda_core_ref_f, core_ref_f, sizeof(float)*SIZE_NUM_REF, cudaMemcpyHostToDevice);
    CUDA_CHK();
    cudaMemcpy(cuda_ref_st_offset, ref_st_offset, sizeof(int32_t)*SIZE_NUM_REF, cudaMemcpyHostToDevice);
    CUDA_CHK();

    // Execute cuda kernal
    dtw_single2<<<1, SIZE>>>(cuda_core, cuda_len_raw_signal, cuda_et_n, cuda_qstart, 
                            cuda_qend, cuda_query, cuda_rlen, max_qlen, max_rlen, cuda_cost, 
                            cuda_core_ref_f, cuda_db_aln_score, cuda_db_aln_score2, cuda_db_aln_pos_st, 
                            cuda_db_aln_pos_end, cuda_db_aln_rid, cuda_db_aln_d, cuda_db_aln_mapq, cuda_aln, 
                            SIZE_NUM_REF, cuda_ref_st_offset, core_opt_flag);

                            // , 
                            // cuda_aln_rid, cuda_aln_pos_st, cuda_aln_pos_end, cuda_aln_score, cuda_aln_score2, cuda_aln_mapq, cuda_aln_d
    CUDA_CHK();

    cudaDeviceSynchronize();
    CUDA_CHK();

    // copy content from cuda memory to main memory
    cudaMemcpy(db_aln_score, cuda_db_aln_score, sizeof(float)*SIZE, cudaMemcpyDeviceToHost);
    CUDA_CHK();
    cudaMemcpy(db_aln_score2, cuda_db_aln_score2, sizeof(float)*SIZE, cudaMemcpyDeviceToHost);
    CUDA_CHK();
    cudaMemcpy(db_aln_pos_st, cuda_db_aln_pos_st, sizeof(int32_t)*SIZE, cudaMemcpyDeviceToHost);
    CUDA_CHK();
    cudaMemcpy(db_aln_pos_end, cuda_db_aln_pos_end, sizeof(int32_t)*SIZE, cudaMemcpyDeviceToHost);
    CUDA_CHK();
    cudaMemcpy(db_aln_rid, cuda_db_aln_rid, sizeof(float)*SIZE, cudaMemcpyDeviceToHost);
    CUDA_CHK();
    cudaMemcpy(db_aln_d, cuda_db_aln_d, sizeof(char)*SIZE, cudaMemcpyDeviceToHost);
    CUDA_CHK();
    cudaMemcpy(db_aln_maqq, cuda_db_aln_mapq, sizeof(uint8_t)*SIZE, cudaMemcpyDeviceToHost);
    CUDA_CHK();

    // Asign values back to structs
    for(int i=0; i<SIZE; i++){
        // fprintf(stderr, "db_aln_score -%i = %i\n",i, db_aln_maqq[i]);
        db->aln[i].score = db_aln_score[i];
        db->aln[i].score2 = db_aln_score2[i];
        db->aln[i].pos_st = db_aln_pos_st[i];
        db->aln[i].pos_end = db_aln_pos_end[i];
        db->aln[i].rid = db_aln_rid[i];
        db->aln[i].d = db_aln_d[i];
        db->aln[i].mapq = db_aln_maqq[i];
    }

    // cudaFree(cuda_core); CUDA_CHK();
    cudaFree(cuda_db); CUDA_CHK();
    cudaFree(cuda_len_raw_signal); CUDA_CHK();
    cudaFree(cuda_et_n); CUDA_CHK();
    cudaFree(cuda_qstart); CUDA_CHK();
    cudaFree(cuda_qend); CUDA_CHK();
    cudaFree(cuda_query); CUDA_CHK();
    cudaFree(cuda_cost); CUDA_CHK();    
    cudaFree(cuda_core_ref_f); CUDA_CHK();
    cudaFree(cuda_db_aln_pos_st); CUDA_CHK();
    cudaFree(cuda_db_aln_pos_end); CUDA_CHK();
    cudaFree(cuda_db_aln_rid); CUDA_CHK();
    cudaFree(cuda_db_aln_d); CUDA_CHK();
    cudaFree(cuda_db_aln_score); CUDA_CHK();    
    cudaFree(cuda_db_aln_score2); CUDA_CHK();

    return;
}