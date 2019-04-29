#include <vector>
#include <random>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/sinc.hpp>
#include <boost/math/special_functions/sin_pi.hpp>
#include <benchmark/benchmark.h>

using boost::math::constants::pi;

template<class RandomAccessContainer>
class whittaker_shannon_1 {
public:
    using Real = typename RandomAccessContainer::value_type;
    whittaker_shannon_1(RandomAccessContainer&& y, Real t0, Real h)
      : y_{std::move(y)}, t0_{t0}, h_{h} {}

    Real operator()(Real t) const {
        Real x = pi<Real>()*(t-t0_)/h_;
        Real s = 0;
        for (size_t i = 0; i < y_.size(); ++i) {
            s += y_[i]*boost::math::sinc_pi(x - i*pi<Real>());
        }
        return s;
    }
private:
    RandomAccessContainer y_;
    Real t0_;
    Real h_;
};


template<class Real>
void BM_WhittakerShannon1(benchmark::State& state) {
    std::vector<Real> v(state.range(0));
    std::mt19937 gen(323723);
    std::uniform_real_distribution<Real> dis(-0.95, 0.95);
    for (size_t i = 0; i < v.size(); ++i) {
        v[i] = dis(gen);
    }

    auto ws = whittaker_shannon_1(std::move(v), Real(0), 1/Real(32));
    Real arg = dis(gen);
    for (auto _ : state) {
        benchmark::DoNotOptimize(ws(arg));
    }
    state.SetComplexityN(state.range(0));
}

BENCHMARK_TEMPLATE(BM_WhittakerShannon1, double)->RangeMultiplier(2)
     ->Range(1<<8, 1<<15)->Complexity(benchmark::oN);

template<class RandomAccessContainer>
class whittaker_shannon_2 {
public:
    using Real = typename RandomAccessContainer::value_type;
    whittaker_shannon_2(RandomAccessContainer&& y, Real t0, Real h)
       : y_{std::move(y)}, t0_{t0}, h_{h} {}

    Real operator()(Real t) const {
        Real x = (t-t0_)/h_;
        Real s = 0;
        for (size_t i = 0; i < y_.size(); ++i) {
            Real term = y_[i]/(x-i);
            if(i & 1) {
                s -= term;
            }
            else {
                s += term;
            }

        }
        return s*sin(pi<Real>()*x)/pi<Real>();
    }
private:
    RandomAccessContainer y_;
    Real t0_;
    Real h_;
};

template<class Real>
void BM_WhittakerShannon2(benchmark::State& state) {
    std::vector<Real> v(state.range(0));
    std::mt19937 gen(323723);
    std::uniform_real_distribution<Real> dis(-0.95, 0.95);
    for (size_t i = 0; i < v.size(); ++i) {
        v[i] = dis(gen);
    }

    auto ws = whittaker_shannon_2(std::move(v), Real(0), 1/Real(32));
    Real arg = dis(gen);
    for (auto _ : state) {
        benchmark::DoNotOptimize(ws(arg));
    }
    state.SetComplexityN(state.range(0));
}

BENCHMARK_TEMPLATE(BM_WhittakerShannon2, double)->RangeMultiplier(2)->Range(1<<8, 1<<15)->Complexity(benchmark::oN);

template<class RandomAccessContainer>
class whittaker_shannon_3 {
public:
    using Real = typename RandomAccessContainer::value_type;
    whittaker_shannon_3(RandomAccessContainer&& y, Real t0, Real h)
     : y_{std::move(y)}, t0_{t0}, h_{h}
    {
        for(size_t i = 1; i < y_.size(); i+= 2) {
            y_[i] = -y_[i];
        }
    }

  Real operator()(Real t) const {
      Real x = (t-t0_)/h_;
      Real s = 0;
      Real z0 = x;
      Real z1 = x - 1;
      Real z2 = x - 2;
      Real z3 = x - 3;

      auto it = y_.begin();
      auto end = y_.end();
      while (it != end) {
          Real k0 = 1/z0;
          Real k1 = 1/z1;
          Real k2 = 1/z2;
          Real k3 = 1/z3;

          s += (*it)*k0;
          s += (*it+1)*k1;
          s += (*it+2)*k2;
          s += (*it+3)*k3;

          z0 -= 4;
          z1 -= 4;
          z2 -= 4;
          z3 -= 4;

          it += 4;
      }
      return s*sin(pi<Real>()*x)/pi<Real>();
  }
private:
  RandomAccessContainer y_;
  Real t0_;
  Real h_;
};


template<class Real>
void BM_WhittakerShannon3(benchmark::State& state) {
   std::vector<Real> v(state.range(0));
   std::mt19937 gen(323723);
   std::uniform_real_distribution<Real> dis(-0.95, 0.95);
   for (size_t i = 0; i < v.size(); ++i) {
       v[i] = dis(gen);
   }

   auto ws = whittaker_shannon_3(std::move(v), Real(0), 1/Real(32));
   Real arg = dis(gen);
   for (auto _ : state) {
       benchmark::DoNotOptimize(ws(arg));
   }
   state.SetComplexityN(state.range(0));
}

BENCHMARK_TEMPLATE(BM_WhittakerShannon3, double)->RangeMultiplier(2)->Range(1<<8, 1<<15)->Complexity(benchmark::oN);

BENCHMARK_MAIN();
