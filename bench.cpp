#include <vector>
#include <random>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/sinc.hpp>
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

BENCHMARK_MAIN();
