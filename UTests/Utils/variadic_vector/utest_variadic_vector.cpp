#define VARIADIC_VECTOR_DEBUG

#include "Utils/variadic_container.hpp"

#include "inria/checker.hpp"

class utest_variadic_vector {

    using vector_t = variadic_vector<std::allocator<char>, char, int, float, double>;
    using value_t  = typename vector_t::value_type;

    inria::checker_t check;

    template<std::size_t I, std::size_t J, class Vec>
    static unsigned long get_byte_offset(const Vec& v) {
        return (unsigned long)(((char*)std::get<J>(v._data_tuple))
                               - ((char*)std::get<I>(v._data_tuple)));
    }

    void test_standard_allocation() {
        std::size_t v_cap = 10;

        vector_t v = vector_t();
        v.reserve(v_cap);

        check.equal(v._capacity, v_cap, LOCUS);
        check.equal(get_byte_offset<0,1>(v), v_cap * sizeof(char),  LOCUS);
        check.equal(get_byte_offset<1,2>(v), v_cap * sizeof(int),   LOCUS);
        check.equal(get_byte_offset<2,3>(v), v_cap * sizeof(float), LOCUS);
    };

    void test_constructor() {
        // Test default constructor
        {
            vector_t v = vector_t();

            check.equal(v.size(),     0u, LOCUS);
            check.equal(v.capacity(), 0u, LOCUS);
        }

        // Test constructor, multiple default constructions
        {
            std::size_t v_size = 10;

            vector_t v(v_size);

            check.equal     (v.size(),     v_size, LOCUS);
            check.greater_eq(v.capacity(), v_size, LOCUS);
            for(std::size_t i = 0; i < v_size; ++i) {
                check.equal(v[i], value_t{}, LOCUS);
            }
        }

        // Test constructor from duplicated value
        {
            std::size_t v_size = 5;
            value_t ref{'a', 1, 3.14, 1e-5};

            vector_t v(v_size, ref);

            check.equal     (v.size(),     v_size, LOCUS);
            check.greater_eq(v.capacity(), v_size, LOCUS);
            for(std::size_t i = 0; i < v_size; ++i) {
                check.equal(v[i], ref, LOCUS);
            }
        }

        // Test constructor from range
        {
            std::initializer_list<value_t> l = {{'b',1,3.14,1e-5}, {'c',2,2.71,2e-6}};

            vector_t v(std::begin(l), std::end(l));

            std::size_t i = 0;
            for(auto ref : l) {
                check.equal(v[i], ref, LOCUS);
                ++i;
            }
        }

        // Test copy constructor
        {
            value_t ref{'a', 1, 3.14, 1e-5};
            vector_t v_ref(5, ref);

            vector_t v(v_ref);

            check.equal(v.size(), v_ref.size(), LOCUS);

            for(std::size_t i = 0; i < v.size(); ++i) {
                check.equal(v[i], ref,      LOCUS);
                check.equal(v[i], v_ref[i], LOCUS);
            }

        }
        // Test move constructor
        {
            value_t ref{'a', 1, 3.14, 1e-5};
            vector_t v_source(5, ref);
            auto v_source_data = v_source.data();
            std::size_t v_source_size = v_source.size();
            std::size_t v_source_capacity = v_source.capacity();

            vector_t v(std::move(v_source));

            check.equal(v.data(),     v_source_data, LOCUS);
            check.equal(v.size(),     v_source_size, LOCUS);
            check.equal(v.capacity(), v_source_capacity, LOCUS);

            check.equal(v_source.size(),     0u, LOCUS);
            check.equal(v_source.capacity(), 0u, LOCUS);

            for(std::size_t i = 0; i < v.size(); ++i) {
                check.equal(v[i], ref, LOCUS);
            }
        }
    }


    void test_assignment_operator() {
        // Test copy assignment
        {
            value_t ref{'a', 1, 3.14, 1e-5};
            const vector_t v_source(5, ref);

            vector_t v = vector_t();
            check.equal(v.size(),     0u, LOCUS);
            check.equal(v.capacity(), 0u, LOCUS);

            v = v_source;

            check.different(v.data(), v_source.data(), LOCUS);
            check.equal(v.size(),     v_source.size(), LOCUS);
            check.equal(v.capacity(), v_source.capacity(), LOCUS);

            for(std::size_t i = 0; i < v.size(); ++i) {
                check.equal(v[i], ref, LOCUS);
            }
        }
    }

public:
    bool run() {
        test_standard_allocation();
        test_constructor();

        check.print_summary();
        return check.ok();
    }
};

int main() {
    return utest_variadic_vector{}.run() ? EXIT_SUCCESS : EXIT_FAILURE;
}
