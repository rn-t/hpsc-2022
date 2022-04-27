#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <immintrin.h>

int main() {
  const int N = 8;
  float x[N], y[N], m[N], fx[N], fy[N], comp[N], temp[N];
  for(int i=0; i<N; i++) {
    x[i] = drand48();
    y[i] = drand48();
    m[i] = drand48();
    fx[i] = fy[i] = 0;
    comp[i] = i;
  }
  for(int i=0; i<N; i++) {
    //for(int j=0; j<N; j++) {
    __m256 xivec = _mm256_set1_ps(x[i]);
    __m256 xjvec = _mm256_load_ps(x);
    __m256 yivec = _mm256_set1_ps(y[i]);
    __m256 yjvec = _mm256_load_ps(y);
      
      //if(i != j) {
      __m256 compvec = _mm256_load_ps(comp);
      __m256 ivec = _mm256_set1_ps(i);
      __m256 mask = _mm256_cmp_ps(compvec, ivec, _CMP_NEQ_OQ);
      __m256 zerovec = _mm256_setzero_ps();
        
        //float rx = x[i] - x[j];
        //float ry = y[i] - y[j];
        __m256 rxvec = _mm256_sub_ps(xivec, xjvec);
        rxvec = _mm256_blendv_ps(zerovec, rxvec, mask);

        __m256 ryvec = _mm256_sub_ps(yivec, yjvec);
        ryvec = _mm256_blendv_ps(zerovec, ryvec, mask);
        
        //float r = std::sqrt(rx * rx + ry * ry);
        __m256 rxvec2 = _mm256_mul_ps(rxvec, rxvec);
        __m256 ryvec2 = _mm256_mul_ps(ryvec, ryvec);
        __m256 rvec2 = _mm256_add_ps(rxvec2, ryvec2);
        __m256 rvec = _mm256_rsqrt_ps(rvec2);

        //fx[i] -= rx * m[j] / (r * r * r);
        //fy[i] -= ry * m[j] / (r * r * r);
        
        __m256 rvec3_2 = _mm256_mul_ps(rvec, rvec);
        rvec3_2 = _mm256_mul_ps(rvec3_2, rvec);
        
        __m256 mjvec = _mm256_load_ps(m);
        
        __m256 fxjvec = _mm256_mul_ps(rxvec, mjvec);
        fxjvec = _mm256_div_ps(fxjvec, rvec3_2);

        __m256 fyjvec = _mm256_mul_ps(ryvec, mjvec);
        fyjvec = _mm256_div_ps(fyjvec, rvec3_2);
        
        __m256 fxjsum = _mm256_permute2f128_ps(fxjvec, fxjvec, 1);
        fxjsum = _mm256_add_ps(fxjvec, fxjsum);
	fxjsum = _mm256_hadd_ps(fxjsum, fxjsum);
	fxjsum = _mm256_hadd_ps(fxjsum, fxjsum);
    
        __m256 fyjsum = _mm256_permute2f128_ps(fyjvec, fyjvec, 1);
        fyjsum = _mm256_add_ps(fyjvec, fyjsum);
	fyjsum = _mm256_hadd_ps(fyjsum, fyjsum);
	fyjsum = _mm256_hadd_ps(fyjsum, fyjsum);

        _mm256_store_ps(temp, fxjsum);
        fx[i] -= temp[0];

        _mm256_store_ps(temp, fyjsum);
        fy[i] -= temp[0];

      //}
    //}
    printf("%d %g %g\n",i,fx[i],fy[i]);
  }
}
