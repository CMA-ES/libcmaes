// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2014 Benoit Steiner <benoit.steiner.goog@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef EIGEN_CXX11_TENSOR_TENSOR_FUNCTORS_H
#define EIGEN_CXX11_TENSOR_TENSOR_FUNCTORS_H

namespace Eigen {
namespace internal {


/** \internal
 * \brief Template functor to compute the modulo between an array and a scalar.
 */
template <typename Scalar>
struct scalar_mod_op {
  EIGEN_DEVICE_FUNC scalar_mod_op(const Scalar& divisor) : m_divisor(divisor) {}
  EIGEN_DEVICE_FUNC inline Scalar operator() (const Scalar& a) const { return a % m_divisor; }
  const Scalar m_divisor;
};
template <typename Scalar>
struct functor_traits<scalar_mod_op<Scalar> >
{ enum { Cost = 2 * NumTraits<Scalar>::MulCost, PacketAccess = false }; };


/** \internal
  * \brief Template functor to compute the sigmoid of a scalar
  * \sa class CwiseUnaryOp, ArrayBase::sigmoid()
  */
template <typename T>
struct scalar_sigmoid_op {
  EIGEN_EMPTY_STRUCT_CTOR(scalar_sigmoid_op)
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE T operator()(const T& x) const {
    const T one = T(1);
    return one / (one + std::exp(-x));
  }

  template <typename Packet>
  inline Packet packetOp(const Packet& x) const {
    const Packet one = pset1<Packet>(1);
    return pdiv(one, padd(one, pexp(pnegate(x))));
  }
};

template <typename T>
struct functor_traits<scalar_sigmoid_op<T> > {
  enum {
    Cost = NumTraits<T>::AddCost * 2 + NumTraits<T>::MulCost * 6,
    PacketAccess = packet_traits<T>::HasAdd && packet_traits<T>::HasDiv &&
                   packet_traits<T>::HasNegate && packet_traits<T>::HasExp
  };
};


// Standard reduction functors
template <typename T> struct SumReducer
{
  static const bool PacketAccess = true;
  static const bool IsStateful = false;

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void reduce(const T t, T* accum) const {
    (*accum) += t;
  }
  template <typename Packet>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void reducePacket(const Packet& p, Packet* accum) const {
    (*accum) = padd<Packet>(*accum, p);
  }

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE T initialize() const {
    return static_cast<T>(0);
  }
  template <typename Packet>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Packet initializePacket() const {
    return pset1<Packet>(0);
  }
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE T finalize(const T accum) const {
    return accum;
  }
  template <typename Packet>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Packet finalizePacket(const Packet& vaccum) const {
    return vaccum;
  }
  template <typename Packet>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE T finalizeBoth(const T saccum, const Packet& vaccum) const {
    return saccum + predux(vaccum);
  }
};

template <typename T> struct MeanReducer
{
  static const bool PacketAccess = true;
  static const bool IsStateful = true;

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
  MeanReducer() : scalarCount_(0), packetCount_(0) { }

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void reduce(const T t, T* accum) {
    (*accum) += t;
    scalarCount_++;
  }
  template <typename Packet>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void reducePacket(const Packet& p, Packet* accum) {
    (*accum) = padd<Packet>(*accum, p);
    packetCount_++;
  }

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE T initialize() const {
    return static_cast<T>(0);
  }
  template <typename Packet>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Packet initializePacket() const {
    return pset1<Packet>(0);
  }
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE T finalize(const T accum) const {
    return accum / scalarCount_;
  }
  template <typename Packet>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Packet finalizePacket(const Packet& vaccum) const {
    return pdiv(vaccum, pset1<Packet>(packetCount_));
  }
  template <typename Packet>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE T finalizeBoth(const T saccum, const Packet& vaccum) const {
    return (saccum + predux(vaccum)) / (scalarCount_ + packetCount_ * unpacket_traits<Packet>::size);
  }

  protected:
    int scalarCount_;
    int packetCount_;
};

template <typename T> struct MaxReducer
{
  static const bool PacketAccess = true;
  static const bool IsStateful = false;

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void reduce(const T t, T* accum) const {
    if (t > *accum) { *accum = t; }
  }
  template <typename Packet>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void reducePacket(const Packet& p, Packet* accum) const {
    (*accum) = pmax<Packet>(*accum, p);
  }

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE T initialize() const {
    return -(std::numeric_limits<T>::max)();
  }
  template <typename Packet>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Packet initializePacket() const {
    return pset1<Packet>(-(std::numeric_limits<T>::max)());
  }
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE T finalize(const T accum) const {
    return accum;
  }
  template <typename Packet>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Packet finalizePacket(const Packet& vaccum) const {
    return vaccum;
  }
  template <typename Packet>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE T finalizeBoth(const T saccum, const Packet& vaccum) const {
    return numext::maxi(saccum, predux_max(vaccum));
  }
};

template <typename T> struct MinReducer
{
  static const bool PacketAccess = true;
  static const bool IsStateful = false;

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void reduce(const T t, T* accum) const {
    if (t < *accum) { *accum = t; }
  }
  template <typename Packet>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void reducePacket(const Packet& p, Packet* accum) const {
    (*accum) = pmin<Packet>(*accum, p);
  }

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE T initialize() const {
    return (std::numeric_limits<T>::max)();
  }
  template <typename Packet>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Packet initializePacket() const {
    return pset1<Packet>((std::numeric_limits<T>::max)());
  }
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE T finalize(const T accum) const {
    return accum;
  }
  template <typename Packet>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Packet finalizePacket(const Packet& vaccum) const {
    return vaccum;
  }
  template <typename Packet>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE T finalizeBoth(const T saccum, const Packet& vaccum) const {
    return numext::mini(saccum, predux_min(vaccum));
  }
};


template <typename T> struct ProdReducer
{
  static const bool PacketAccess = true;
  static const bool IsStateful = false;

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void reduce(const T t, T* accum) const {
    (*accum) *= t;
  }
  template <typename Packet>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void reducePacket(const Packet& p, Packet* accum) const {
    (*accum) = pmul<Packet>(*accum, p);
  }

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE T initialize() const {
    return static_cast<T>(1);
  }
  template <typename Packet>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Packet initializePacket() const {
    return pset1<Packet>(1);
  }
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE T finalize(const T accum) const {
    return accum;
  }
  template <typename Packet>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Packet finalizePacket(const Packet& vaccum) const {
    return vaccum;
  }
  template <typename Packet>
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE T finalizeBoth(const T saccum, const Packet& vaccum) const {
    return saccum * predux_mul(vaccum);
  }
};


struct AndReducer
{
  static const bool PacketAccess = false;
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void reduce(bool t, bool* accum) const {
    *accum = *accum && t;
  }
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE bool initialize() const {
    return true;
  }
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE bool finalize(bool accum) const {
    return accum;
  }
};

struct OrReducer {
  static const bool PacketAccess = false;
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void reduce(bool t, bool* accum) const {
    *accum = *accum || t;
  }
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE bool initialize() const {
    return false;
  }
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE bool finalize(bool accum) const {
    return accum;
  }
};

// Argmin/Argmax reducers
template <typename T> struct ArgMaxTupleReducer
{
  static const bool PacketAccess = false;
  static const bool IsStateful = false;

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void reduce(const T t, T* accum) const {
    if (t.second > accum->second) { *accum = t; }
  }
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE T initialize() const {
    return T(0, NumTraits<typename T::second_type>::lowest());
  }
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE T finalize(const T& accum) const {
    return accum;
  }
};

template <typename T> struct ArgMinTupleReducer
{
  static const bool PacketAccess = false;
  static const bool IsStateful = false;

  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE void reduce(const T& t, T* accum) const {
    if (t.second < accum->second) { *accum = t; }
  }
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE T initialize() const {
    return T(0, NumTraits<typename T::second_type>::highest());
  }
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE T finalize(const T& accum) const {
    return accum;
  }
};


// Random number generation
namespace {
#ifdef __CUDA_ARCH__
__device__ int get_random_seed() {
    return clock();
}
#else
int get_random_seed() {
#ifdef _WIN32
    SYSTEMTIME st;
    GetSystemTime(&st);
    return st.wSecond + 1000 * st.wMilliseconds;
#elif defined __APPLE__
    return static_cast<int>(mach_absolute_time());
#else
    timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return static_cast<int>(ts.tv_nsec);
#endif
}
#endif
}

#if !defined (EIGEN_USE_GPU) || !defined(__CUDACC__) || !defined(__CUDA_ARCH__)
// We're not compiling a cuda kernel
template <typename T> class UniformRandomGenerator {

 public:
  static const bool PacketAccess = true;

  UniformRandomGenerator(bool deterministic = true) : m_deterministic(deterministic) {
    if (!deterministic) {
      srand(get_random_seed());
    }
  }
  UniformRandomGenerator(const UniformRandomGenerator& other) {
    m_deterministic = other.m_deterministic;
  }

  template<typename Index>
  T operator()(Index, Index = 0) const {
    return random<T>();
  }
  template<typename Index>
  typename internal::packet_traits<T>::type packetOp(Index, Index = 0) const {
    const int packetSize = internal::packet_traits<T>::size;
    EIGEN_ALIGN_MAX T values[packetSize];
    for (int i = 0; i < packetSize; ++i) {
      values[i] = random<T>();
    }
    return internal::pload<typename internal::packet_traits<T>::type>(values);
  }

 private:
  bool m_deterministic;
};

#if __cplusplus > 199711
template <> class UniformRandomGenerator<float> {
 public:
  static const bool PacketAccess = true;

  UniformRandomGenerator(bool deterministic = true) : m_deterministic(deterministic) {
    if (!deterministic) {
      m_generator.seed(get_random_seed());
    }
  }
  UniformRandomGenerator(const UniformRandomGenerator<float>& other) {
    m_generator.seed(other(0, 0) * UINT_MAX);
    m_deterministic = other.m_deterministic;
  }

  template<typename Index>
  float operator()(Index, Index = 0) const {
    return m_distribution(m_generator);
  }
  template<typename Index>
  typename internal::packet_traits<float>::type packetOp(Index i, Index j = 0) const {
    const int packetSize = internal::packet_traits<float>::size;
    EIGEN_ALIGN_MAX float values[packetSize];
    for (int k = 0; k < packetSize; ++k) {
      values[k] = this->operator()(i, j);
    }
    return internal::pload<typename internal::packet_traits<float>::type>(values);
  }

 private:
  UniformRandomGenerator& operator = (const UniformRandomGenerator&);
  // Make sure m_deterministic comes first to match the layout of the cpu
  // version of the code.
  bool m_deterministic;
  mutable std::mt19937 m_generator;
  mutable std::uniform_real_distribution<float> m_distribution;
};

template <> class UniformRandomGenerator<double> {
 public:
  static const bool PacketAccess = true;

  UniformRandomGenerator(bool deterministic = true) : m_deterministic(deterministic) {
    if (!deterministic) {
      m_generator.seed(get_random_seed());
    }
  }
  UniformRandomGenerator(const UniformRandomGenerator<double>& other) {
    m_generator.seed(other(0, 0) * UINT_MAX);
    m_deterministic = other.m_deterministic;
  }

  template<typename Index>
  double operator()(Index, Index = 0) const {
    return m_distribution(m_generator);
  }
  template<typename Index>
  typename internal::packet_traits<double>::type packetOp(Index i, Index j = 0) const {
    const int packetSize = internal::packet_traits<double>::size;
    EIGEN_ALIGN_MAX double values[packetSize];
    for (int k = 0; k < packetSize; ++k) {
      values[k] = this->operator()(i, j);
    }
    return internal::pload<typename internal::packet_traits<double>::type>(values);
  }

 private:
  UniformRandomGenerator& operator = (const UniformRandomGenerator&);
  // Make sure m_deterministic comes first to match the layout of the cpu
  // version of the code.
  bool m_deterministic;
  mutable std::mt19937 m_generator;
  mutable std::uniform_real_distribution<double> m_distribution;
};
#endif

#else

// We're compiling a cuda kernel
template <typename T> class UniformRandomGenerator;

template <> class UniformRandomGenerator<float> {
 public:
  static const bool PacketAccess = true;

  __device__ UniformRandomGenerator(bool deterministic = true) : m_deterministic(deterministic) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int seed = deterministic ? 0 : get_random_seed();
    curand_init(seed, tid, 0, &m_state);
  }

  __device__ UniformRandomGenerator(const UniformRandomGenerator& other) {
    m_deterministic = other.m_deterministic;
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int seed = m_deterministic ? 0 : get_random_seed();
     curand_init(seed, tid, 0, &m_state);
  }

  template<typename Index>
  __device__ float operator()(Index, Index = 0) const {
    return curand_uniform(&m_state);
  }
  template<typename Index>
  __device__ float4 packetOp(Index, Index = 0) const {
    return curand_uniform4(&m_state);
  }

 private:
  bool m_deterministic;
  mutable curandStatePhilox4_32_10_t m_state;
};

template <> class UniformRandomGenerator<double> {
 public:
  static const bool PacketAccess = true;

  __device__ UniformRandomGenerator(bool deterministic = true) : m_deterministic(deterministic) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int seed = deterministic ? 0 : get_random_seed();
    curand_init(seed, tid, 0, &m_state);
  }
  __device__ UniformRandomGenerator(const UniformRandomGenerator& other) {
    m_deterministic = other.m_deterministic;
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int seed = m_deterministic ? 0 : get_random_seed();
    curand_init(seed, tid, 0, &m_state);
  }
  template<typename Index>
  __device__ double operator()(Index, Index = 0) const {
    return curand_uniform_double(&m_state);
  }
  template<typename Index>
  __device__ double2 packetOp(Index, Index = 0) const {
    return curand_uniform2_double(&m_state);
  }

 private:
  bool m_deterministic;
  mutable curandStatePhilox4_32_10_t m_state;
};

template <> class UniformRandomGenerator<std::complex<float> > {
 public:
  static const bool PacketAccess = false;

  __device__ UniformRandomGenerator(bool deterministic = true) : m_deterministic(deterministic) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int seed = deterministic ? 0 : get_random_seed();
    curand_init(seed, tid, 0, &m_state);
  }
  __device__ UniformRandomGenerator(const UniformRandomGenerator& other) {
    m_deterministic = other.m_deterministic;
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int seed = m_deterministic ? 0 : get_random_seed();
    curand_init(seed, tid, 0, &m_state);
  }
  template<typename Index>
  __device__ std::complex<float> operator()(Index, Index = 0) const {
    float4 vals = curand_uniform4(&m_state);
    return std::complex<float>(vals.x, vals.y);
  }

 private:
  bool m_deterministic;
  mutable curandStatePhilox4_32_10_t m_state;
};

template <> class UniformRandomGenerator<std::complex<double> > {
 public:
  static const bool PacketAccess = false;

  __device__ UniformRandomGenerator(bool deterministic = true) : m_deterministic(deterministic) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int seed = deterministic ? 0 : get_random_seed();
    curand_init(seed, tid, 0, &m_state);
  }
  __device__ UniformRandomGenerator(const UniformRandomGenerator& other) {
    m_deterministic = other.m_deterministic;
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int seed = m_deterministic ? 0 : get_random_seed();
    curand_init(seed, tid, 0, &m_state);
  }
  template<typename Index>
  __device__ std::complex<double> operator()(Index, Index = 0) const {
    double2 vals = curand_uniform2_double(&m_state);
    return std::complex<double>(vals.x, vals.y);
  }

 private:
  bool m_deterministic;
  mutable curandStatePhilox4_32_10_t m_state;
};

#endif


#if (!defined (EIGEN_USE_GPU) || !defined(__CUDACC__) || !defined(__CUDA_ARCH__)) && __cplusplus > 199711
// We're not compiling a cuda kernel
template <typename T> class NormalRandomGenerator {
 public:
  static const bool PacketAccess = true;

  NormalRandomGenerator(bool deterministic = true) : m_deterministic(deterministic), m_distribution(0, 1) {
    if (!deterministic) {
      m_generator.seed(get_random_seed());
    }
  }
  NormalRandomGenerator(const NormalRandomGenerator& other)
      : m_deterministic(other.m_deterministic), m_distribution(other.m_distribution) {
    m_generator.seed(other(0, 0) * UINT_MAX);
  }

  template<typename Index>
  T operator()(Index, Index = 0) const {
    return m_distribution(m_generator);
  }
  template<typename Index>
  typename internal::packet_traits<T>::type packetOp(Index, Index = 0) const {
    const int packetSize = internal::packet_traits<T>::size;
    EIGEN_ALIGN_MAX T values[packetSize];
    for (int i = 0; i < packetSize; ++i) {
      values[i] = m_distribution(m_generator);
    }
    return internal::pload<typename internal::packet_traits<T>::type>(values);
  }

 private:
  bool m_deterministic;
  mutable std::normal_distribution<T> m_distribution;
  mutable std::mt19937 m_generator;
};

#elif defined (EIGEN_USE_GPU) && defined(__CUDACC__) && defined(__CUDA_ARCH__)

// We're compiling a cuda kernel
template <typename T> class NormalRandomGenerator;

template <> class NormalRandomGenerator<float> {
 public:
  static const bool PacketAccess = true;

  __device__ NormalRandomGenerator(bool deterministic = true) : m_deterministic(deterministic) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int seed = deterministic ? 0 : get_random_seed();
    curand_init(seed, tid, 0, &m_state);
  }
  __device__ NormalRandomGenerator(const NormalRandomGenerator<float>& other) {
    m_deterministic = other.m_deterministic;
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int seed = m_deterministic ? 0 : get_random_seed();
    curand_init(seed, tid, 0, &m_state);
  }
  template<typename Index>
   __device__ float operator()(Index, Index = 0) const {
    return curand_normal(&m_state);
  }
  template<typename Index>
   __device__ float4 packetOp(Index, Index = 0) const {
    return curand_normal4(&m_state);
  }

 private:
  bool m_deterministic;
  mutable curandStatePhilox4_32_10_t m_state;
};

template <> class NormalRandomGenerator<double> {
 public:
  static const bool PacketAccess = true;

  __device__ NormalRandomGenerator(bool deterministic = true) : m_deterministic(deterministic) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int seed = deterministic ? 0 : get_random_seed();
    curand_init(seed, tid, 0, &m_state);
  }
  __device__ NormalRandomGenerator(const NormalRandomGenerator<double>& other) {
    m_deterministic = other.m_deterministic;
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int seed = m_deterministic ? 0 : get_random_seed();
    curand_init(seed, tid, 0, &m_state);
  }
  template<typename Index>
  __device__ double operator()(Index, Index = 0) const {
    return curand_normal_double(&m_state);
  }
  template<typename Index>
  __device__ double2 packetOp(Index, Index = 0) const {
    return curand_normal2_double(&m_state);
  }

 private:
  bool m_deterministic;
  mutable curandStatePhilox4_32_10_t m_state;
};

template <> class NormalRandomGenerator<std::complex<float> > {
 public:
  static const bool PacketAccess = false;

  __device__ NormalRandomGenerator(bool deterministic = true) : m_deterministic(deterministic) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int seed = deterministic ? 0 : get_random_seed();
    curand_init(seed, tid, 0, &m_state);
  }
  __device__ NormalRandomGenerator(const NormalRandomGenerator& other) {
    m_deterministic = other.m_deterministic;
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int seed = m_deterministic ? 0 : get_random_seed();
    curand_init(seed, tid, 0, &m_state);
  }
  template<typename Index>
  __device__ std::complex<float> operator()(Index, Index = 0) const {
    float4 vals = curand_normal4(&m_state);
    return std::complex<float>(vals.x, vals.y);
  }

 private:
  bool m_deterministic;
  mutable curandStatePhilox4_32_10_t m_state;
};

template <> class NormalRandomGenerator<std::complex<double> > {
 public:
  static const bool PacketAccess = false;

  __device__ NormalRandomGenerator(bool deterministic = true) : m_deterministic(deterministic) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int seed = deterministic ? 0 : get_random_seed();
    curand_init(seed, tid, 0, &m_state);
  }
  __device__ NormalRandomGenerator(const NormalRandomGenerator& other) {
    m_deterministic = other.m_deterministic;
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int seed = m_deterministic ? 0 : get_random_seed();
    curand_init(seed, tid, 0, &m_state);
  }
  template<typename Index>
  __device__ std::complex<double> operator()(Index, Index = 0) const {
    double2 vals = curand_normal2_double(&m_state);
    return std::complex<double>(vals.x, vals.y);
  }

 private:
  bool m_deterministic;
  mutable curandStatePhilox4_32_10_t m_state;
};

#else

template <typename T> class NormalRandomGenerator {
 public:
  NormalRandomGenerator(bool deterministic = true) : m_deterministic(deterministic) {}

 private:
  bool m_deterministic;
};

#endif


template <typename T, typename Index, size_t NumDims>
class GaussianGenerator {
 public:
  static const bool PacketAccess = false;

  EIGEN_DEVICE_FUNC GaussianGenerator(const array<T, NumDims>& means,
                                      const array<T, NumDims>& std_devs)
      : m_means(means)
  {
    for (size_t i = 0; i < NumDims; ++i) {
      m_two_sigmas[i] = std_devs[i] * std_devs[i] * 2;
    }
  }

  T operator()(const array<Index, NumDims>& coordinates) const {
    T tmp = T(0);
    for (size_t i = 0; i < NumDims; ++i) {
      T offset = coordinates[i] - m_means[i];
      tmp += offset * offset / m_two_sigmas[i];
    }
    return std::exp(-tmp);
  }

 private:
  array<T, NumDims> m_means;
  array<T, NumDims> m_two_sigmas;
};


} // end namespace internal
} // end namespace Eigen

#endif // EIGEN_CXX11_TENSOR_TENSOR_FUNCTORS_H
