#include <cstdio>
#include <cstdlib>
#include <vector>

__global__ void zerofill(int *bucket, int N){
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index > N) return;
  bucket[index] = 0;
}

__global__ void bucket_in(int *bucket, int *key, int N){
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index > N) return;
  atomicAdd(&bucket[key[index]], 1);
}

__global__ void key_in(int *bucket, int *key, int i, int j, int N){
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  if (index > N) return;
  if (j <= index && index < j + bucket[i]) key[index] = i;
}

int main() {
  const int BLOCK = 1024;
  int n = 50;
  int range = 5;

  //std::vector<int> key(n);
  int *key;
  cudaMallocManaged(&key, n * sizeof(int));

  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]);
  }
  printf("\n");

  //std::vector<int> bucket(range); 
  int *bucket;
  cudaMallocManaged(&bucket, range * sizeof(int));
  
/*
for (int i=0; i<range; i++) {
    bucket[i] = 0;
  }
*/
  zerofill<<<(range + BLOCK - 1)/BLOCK,BLOCK>>>(bucket, range);
  cudaDeviceSynchronize();

/*
  for (int i=0; i<n; i++) {
    bucket[key[i]]++;
  }
*/
  bucket_in<<<(range + BLOCK - 1)/BLOCK,BLOCK>>>(bucket, key, n);
  cudaDeviceSynchronize();
  
/*
  for (int i=0, j=0; i<range; i++) {
    for (; bucket[i]>0; bucket[i]--) {
      key[j++] = i;
    }
  }
*/
  for (int i=0, j=0; i<range; i++) {
    key_in<<<(range + BLOCK - 1)/BLOCK,BLOCK>>>(bucket, key, i, j, n);
    cudaDeviceSynchronize();
    j += bucket[i];
  }

  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");
  cudafree(key);
  cudafree(bucket);
}
