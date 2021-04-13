#ifndef FTIMER_HPP
#define FTIMER_HPP

#include <chrono>
#include <vector>

class FTimer {

public:
    using clock_t = std::chrono::high_resolution_clock;
    using duration_t = std::chrono::duration<double>;

protected:
    std::vector<duration_t> _measures;

public:
    std::vector<duration_t>& measures() {
        return _measures;
    }
    const std::vector<duration_t>& measures() const {
        return _measures;
    }

    duration_t last() const {
        return _measures.back();
    }

    template<typename F, typename... Args>
    void time(F&& func, Args&&... args) {
        auto start = clock_t::now();
        std::forward<F>(func)(std::forward<Args>(args)...);
        this->_measures.push_back(clock_t::now() - start);
    }

    template<typename F, typename... Args>
    void time(std::size_t i, F&& func, Args&&... args) {
        auto start = clock_t::now();
        std::forward<F>(func)(std::forward<Args>(args)...);
        this->_measures.at(i) += (clock_t::now() - start);
    }

};



#endif /* FTIMER_HPP */
