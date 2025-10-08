#pragma once

// SIMD compatibility layer to allow building on platforms without AVX2
// support (for example, Apple Silicon). When AVX is available we defer to
// the native intrinsics, otherwise we provide a scalar fallback that mimics
// the API used by the GLIMPSE models.
//
// ⚠️  PERFORMANCE NOTE (WIP):
// The scalar fallback implementation is SIGNIFICANTLY SLOWER than native AVX2.
// This is a functional implementation to enable compilation on ARM platforms,
// but performance optimization using ARM NEON intrinsics is planned for future
// releases. Users on Apple Silicon should expect reduced performance compared
// to x86_64 systems with AVX2 support.

#if defined(__AVX2__) || defined(__AVX__)
  #include <immintrin.h>
  #define GLIMPSE_SIMD_HAS_AVX 1
#else
  #define GLIMPSE_SIMD_HAS_AVX 0
  #include <algorithm>
  #include <array>
  #include <cstdint>
  #include <cstring>

  namespace glimpse::simd
  {
    struct alignas(32) Vec256f {
      float data[8];
    };

    struct alignas(32) Vec256i {
      std::uint32_t data[8];
    };

    struct alignas(16) Vec128f {
      float data[4];
    };

    template <typename T>
    inline T bit_cast_copy(const void* src)
    {
      T dst;
      std::memcpy(&dst, src, sizeof(T));
      return dst;
    }

    template <typename T>
    inline void bit_copy_to(const T& src, void* dst)
    {
      std::memcpy(dst, &src, sizeof(T));
    }
  }

  using __m256  = glimpse::simd::Vec256f;
  using __m256i = glimpse::simd::Vec256i;
  using __m128  = glimpse::simd::Vec128f;

  // Floating point helpers -------------------------------------------------
  inline __m256 _mm256_set1_ps(float value)
  {
    __m256 result;
    for (int i = 0; i < 8; ++i) result.data[i] = value;
    return result;
  }

  inline __m256 _mm256_load_ps(const float* ptr)
  {
    return glimpse::simd::bit_cast_copy<__m256>(ptr);
  }

  inline void _mm256_store_ps(float* ptr, const __m256& value)
  {
    glimpse::simd::bit_copy_to(value, ptr);
  }

  inline __m256 _mm256_add_ps(const __m256& a, const __m256& b)
  {
    __m256 result;
    for (int i = 0; i < 8; ++i) result.data[i] = a.data[i] + b.data[i];
    return result;
  }

  inline __m256 _mm256_mul_ps(const __m256& a, const __m256& b)
  {
    __m256 result;
    for (int i = 0; i < 8; ++i) result.data[i] = a.data[i] * b.data[i];
    return result;
  }

  inline __m256 _mm256_div_ps(const __m256& a, const __m256& b)
  {
    __m256 result;
    for (int i = 0; i < 8; ++i) result.data[i] = a.data[i] / b.data[i];
    return result;
  }

  inline __m256 _mm256_max_ps(const __m256& a, const __m256& b)
  {
    __m256 result;
    for (int i = 0; i < 8; ++i) result.data[i] = std::max(a.data[i], b.data[i]);
    return result;
  }

  inline __m256 _mm256_min_ps(const __m256& a, const __m256& b)
  {
    __m256 result;
    for (int i = 0; i < 8; ++i) result.data[i] = std::min(a.data[i], b.data[i]);
    return result;
  }

  inline __m256 _mm256_fmadd_ps(const __m256& a, const __m256& b, const __m256& c)
  {
    return _mm256_add_ps(_mm256_mul_ps(a, b), c);
  }

  inline __m256 _mm256_setzero_ps()
  {
    return _mm256_set1_ps(0.0f);
  }

  // Integer helpers --------------------------------------------------------
  inline __m256i _mm256_set1_epi32(int value)
  {
    __m256i result;
    std::uint32_t converted = static_cast<std::uint32_t>(value);
    for (int i = 0; i < 8; ++i) result.data[i] = converted;
    return result;
  }

  inline __m256i _mm256_set_epi32(int e7, int e6, int e5, int e4, int e3, int e2, int e1, int e0)
  {
    __m256i result;
    result.data[0] = static_cast<std::uint32_t>(e0);
    result.data[1] = static_cast<std::uint32_t>(e1);
    result.data[2] = static_cast<std::uint32_t>(e2);
    result.data[3] = static_cast<std::uint32_t>(e3);
    result.data[4] = static_cast<std::uint32_t>(e4);
    result.data[5] = static_cast<std::uint32_t>(e5);
    result.data[6] = static_cast<std::uint32_t>(e6);
    result.data[7] = static_cast<std::uint32_t>(e7);
    return result;
  }

  inline __m256i _mm256_sllv_epi32(const __m256i& values, const __m256i& counts)
  {
    __m256i result;
    for (int i = 0; i < 8; ++i)
    {
      std::uint32_t shift = counts.data[i] & 31u;
      result.data[i] = values.data[i] << shift;
    }
    return result;
  }

  inline __m256 _mm256_castsi256_ps(const __m256i& value)
  {
    return glimpse::simd::bit_cast_copy<__m256>(value.data);
  }

  inline __m256i _mm256_castps_si256(const __m256& value)
  {
    return glimpse::simd::bit_cast_copy<__m256i>(value.data);
  }

  inline __m256 _mm256_blendv_ps(const __m256& a, const __m256& b, const __m256& mask)
  {
    __m256 result;
    for (int i = 0; i < 8; ++i)
    {
      std::uint32_t bits = glimpse::simd::bit_cast_copy<std::uint32_t>(&mask.data[i]);
      result.data[i] = (bits & 0x80000000u) ? b.data[i] : a.data[i];
    }
    return result;
  }

  // Mixed helpers ----------------------------------------------------------
  inline __m128 _mm256_castps256_ps128(const __m256& value)
  {
    __m128 result;
    for (int i = 0; i < 4; ++i) result.data[i] = value.data[i];
    return result;
  }

  inline __m128 _mm256_extractf128_ps(const __m256& value, const int index)
  {
    __m128 result;
    const int offset = (index & 1) ? 4 : 0;
    for (int i = 0; i < 4; ++i) result.data[i] = value.data[offset + i];
    return result;
  }

  // SSE fallbacks ----------------------------------------------------------
  inline __m128 _mm_add_ps(const __m128& a, const __m128& b)
  {
    __m128 result;
    for (int i = 0; i < 4; ++i) result.data[i] = a.data[i] + b.data[i];
    return result;
  }

  inline __m128 _mm_movehdup_ps(const __m128& a)
  {
    __m128 result;
    result.data[0] = a.data[1];
    result.data[1] = a.data[1];
    result.data[2] = a.data[3];
    result.data[3] = a.data[3];
    return result;
  }

  inline __m128 _mm_movehl_ps(const __m128& a, const __m128& b)
  {
    __m128 result;
    result.data[0] = b.data[2];
    result.data[1] = b.data[3];
    result.data[2] = a.data[2];
    result.data[3] = a.data[3];
    return result;
  }

  inline __m128 _mm_add_ss(const __m128& a, const __m128& b)
  {
    __m128 result = a;
    result.data[0] = a.data[0] + b.data[0];
    return result;
  }

  inline float _mm_cvtss_f32(const __m128& a)
  {
    return a.data[0];
  }

  // Convenience helpers ----------------------------------------------------
  #endif  // GLIMPSE_SIMD_HAS_AVX
