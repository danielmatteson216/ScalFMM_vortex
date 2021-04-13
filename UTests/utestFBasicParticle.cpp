#ifdef NDEBUG // We want the checks to be run even in release
#undef NDEBUG
#endif

#include <cassert>

#include "Components/FBasicParticle.hpp"

#include "Adaptive/FVariadicParticleContainer.hpp"


namespace sd = scalfmm::details;

namespace sfinae {
    template<typename ... Ts> using void_t = void;

}

struct test_FBasicParticle_simple {


    using FReal = double;
    enum {Dim = 3};

    using Particle = FBasicParticle<FReal, Dim, float, scalfmm::pack<4, FReal> >;

    // Check the attribute count
    static_assert(Particle::NbAttributes == 5, "Wrong count of attributes");

    using Container = FVariadicParticleContainer<Particle>;

    void run() {
        bool asserts_are_run = false;
        assert(asserts_are_run = true);
        if(!asserts_are_run) {
            std::cerr << "Asserts are not run";
            exit(-1);
        }

        test_default_constructor();
        test_variadic_constructor();
        test_incomplete_variadic_constructor();
        test_variadic_constructor_overflow_failure();
        test_constructor_from_tuple();
        test_position_getter();
        test_position_setter();
        test_attribute_getter();
        test_attribute_setter();
        test_compare_equal();
        test_pull_push();
    }

    /**
     * \brief Test the default constructor
     *
     * Expected result: particle builds
     */
    void test_default_constructor() {
        Particle p{};
    }

    /**
     * \brief Test the variadic constructor
     *
     * Expected result:
     *   - particle builds
     *   - position is right
     *   - attributes are right
     */
    void test_variadic_constructor() {
        const Particle::position_t pos(1.5,2.5,3.5);
        Particle p(pos,1,2,3,4,5);

        assert(1.5 == pos[0]);
        assert(2.5 == pos[1]);
        assert(3.5 == pos[2]);

        assert(std::get<3>(p) == 1);
        assert(std::get<4>(p) == 2);
        assert(std::get<5>(p) == 3);
        assert(std::get<6>(p) == 4);
        assert(std::get<7>(p) == 5);
    }

    /**
     * \brief Test variadic constructor with incomplete arguments
     *
     * Expected result:
     *   - particle builds
     *   - position is right
     *   - attributes are right
     *   - non given attributes are zeroed
     */
    void test_incomplete_variadic_constructor() {
        const Particle::position_t pos(1.5,2.5,3.5);
        Particle p(pos,1,2);

        assert(1.5 == pos[0]);
        assert(2.5 == pos[1]);
        assert(3.5 == pos[2]);

        assert(std::get<3>(p) == 1);
        assert(std::get<4>(p) == 2);
        assert(std::get<5>(p) == 0);
        assert(std::get<6>(p) == 0);
        assert(std::get<7>(p) == 0);
    }

    /**
     * \brief Check if constructor can be called with specified attributes.
     *
     * Holds a #value that is false if constructor `U({0,0,0}, attributes...)`
     * exists.
     *
     * \tparam T Class to check
     * \tparam attributes int parameter pack of attribute values
     */
    template<class T, int ... attributes>
    struct check_constructor_overflow_failure {
        /// Exists if the constructor `U({0,0,0}, attributes...)` exists
        template<class U, class = decltype(U({0,0,0}, attributes...))>
        static constexpr std::false_type check(U*) {return {};}
        /// SFINAE fallback
        static constexpr std::true_type check(...) {return {};}

        /// value is true if the constructor does not exist
        enum {value = decltype(check((T*)nullptr))::value};
    };

    /**
     * \brief Test variadic constructor when too many arguments are given
     *
     * Expected result:
     *   - the constructor instanciated raises a false static_assert
     */
    void test_variadic_constructor_overflow_failure() {
        using test = check_constructor_overflow_failure<Particle,1,2,3,4,5,6>;
        static_assert(test::value, "Variadic constructor accepts too many attributes");
        assert(test::value);
    }

    /**
     * \brief Test constructor from a std::tuple
     *
     * Expected result:
     *   - the particle is created
     *   - the position and attributes are correct
     */
    void test_constructor_from_tuple() {
        Particle p(std::make_tuple(1.5,2.5,3.5,1,2,3,4,5));
        // position
        assert(std::get<0>(p) == 1.5);
        assert(std::get<1>(p) == 2.5);
        assert(std::get<2>(p) == 3.5);
        // attributes
        assert(std::get<3>(p) == 1);
        assert(std::get<4>(p) == 2);
        assert(std::get<5>(p) == 3);
        assert(std::get<6>(p) == 4);
        assert(std::get<7>(p) == 5);
    }


    /**
     * \brief Test position getter
     *
     * Expected result:
     *   - p.position() returns a copy of the particle position
     */
    void test_position_getter() {
        typename Particle::position_t ref{1.5,2.5,3.5};
        const Particle p(ref, 1,2,3,4,5);
        // check value
        auto pos = p.position();
        assert(pos == ref);
        // check that a copy is returned, not a reference
        pos[0] += 1;
        assert(pos != p.position());
    }

    /**
     * \brief Tests the position setter
     *
     * Expected result:
     *   - setting the position changes the particle position
     */
    void test_position_setter() {
        typename Particle::position_t ref{1.5,2.5,3.5};
        Particle p(ref, 1,2,3,4,5);

        p.setPosition({10.25,11.25,12.25});

        assert(std::get<0>(p) == 10.25);
        assert(std::get<1>(p) == 11.25);
        assert(std::get<2>(p) == 12.25);
    }


    /**
     * \brief Test compile time attribute getter
     *
     * Expected result:
     *   - the attributes are returned using indices starting from 0
     */
    void test_attribute_getter() {
        // const
        const Particle p1(std::make_tuple(1.5,2.5,3.5,1,2,3,4,5));
        assert(p1.attribute<0>() == 1);
        assert(p1.attribute<1>() == 2);
        assert(p1.attribute<2>() == 3);
        assert(p1.attribute<3>() == 4);
        assert(p1.attribute<4>() == 5);
        // non const
        Particle p2(std::make_tuple(1.5,2.5,3.5,1,2,3,4,5));
        assert(p2.attribute<0>() == 1);
        assert(p2.attribute<1>() == 2);
        assert(p2.attribute<2>() == 3);
        assert(p2.attribute<3>() == 4);
        assert(p2.attribute<4>() == 5);
    }

    /**
     * \brief Test compile time attribute setter
     *
     * Expected result:
     *   - the attributes are set using indices starting from 0
     */
    void test_attribute_setter() {
        Particle p1(std::make_tuple(1.5,2.5,3.5,0,0,0,0,0));
        p1.attribute<0>() = 1;
        p1.attribute<1>() = 2;
        p1.attribute<2>() = 3;
        p1.attribute<3>() = 4;
        p1.attribute<4>() = 5;

        assert(p1.attribute<0>() == 1);
        assert(p1.attribute<1>() == 2);
        assert(p1.attribute<2>() == 3);
        assert(p1.attribute<3>() == 4);
        assert(p1.attribute<4>() == 5);
    }

    /**
     * \brief Test particle comparison
     *
     * Expected result:
     *   - Two particles that have the same position and the same attributes
     *     compare equal.
     */
    void test_compare_equal() {
        Particle p1(std::make_tuple(1.5,2.5,3.5,1,2,3,4,5));
        Particle p2(std::make_tuple(1.5,2.5,3.5,1,2,3,4,5));
        Particle p3(std::make_tuple(1.5,2.5,3.5,1,2,0,4,5));
        Particle p4(std::make_tuple(1.5,  0,3.5,1,2,3,4,5));

        assert(p1 == p2);
        assert(p2 == p1);
        assert(p1 != p3);
        assert(p1 != p4);
        assert(p3 != p4);

    }

    void test_pull_push() {
        Particle p({2.3,3.4,4.5}, 5, 6, 7, 8, 9);

        Container container;
        Container other_container;

        container.push(p);
        other_container.push(container[0]);


        assert(Particle(container[0]) == p);
        assert(container.size() == 1);
        assert(Particle(other_container[0]) == p);


    }



};



int main() {


    test_FBasicParticle_simple().run();

}
