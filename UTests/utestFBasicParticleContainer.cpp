#undef NDEBUG
#include <cassert>

#include <vector>

#include "Utils/FOstreamTuple.hpp"

#include "Components/FBasicParticleContainer.hpp"

using FReal = double;

namespace sc = scalfmm;

using sc::pack;
using sc::pack_expand;
using sc::pack_expand_tuple;

class testFBasicParticleContainer {

    using FContainer    = FBasicParticleContainer<FReal, 5, int>;
    using FParticle     = FContainer::FParticle;
    using FParticleData = FParticle::tuple_data_t;

    static_assert(
		  std::is_same<
		  FParticleData,
		  pack_expand_tuple<pack<3,FReal>, pack<5,int> > >::value,
		  "");


    using FPosition  = FPoint<FReal>;

    std::vector<FPosition> pos_vect = {
        {1,5, 9},
        {2,6,10},
        {3,7,11},
        {4,8,12},
    };

    std::vector<FParticle> part_vect = {
        {this->pos_vect[0],  1, 2, 3, 4, 5},
        {this->pos_vect[1],  6, 7, 8, 9,10},
        {this->pos_vect[2], 11,12,13,14,15},
        {this->pos_vect[3], 16,17,18,19,20},
    };


public:
    void run() {
        test_constructor_1();
        test_push();
        test_push_particle();
        test_positions();
        test_attributes();
    }

private:
    void test_constructor_1() {
        FContainer container;

        assert(container.size() == 0);
        assert(container.getNbParticles() == 0);
    }

    void test_push() {
        FContainer container;
        container.push(pos_vect[0], 1,2,3,4,5);

        assert(container.size() == 1);
        assert(container.getNbParticles() == 1);

        assert(container[0] == part_vect[0].as_tuple() );

        container.push(pos_vect[1], 6,7,8,9,10);

        assert(container.size() == 2);
        assert(container.getNbParticles() == 2);

        assert(container[1] == part_vect[1].as_tuple() );
    }

    void test_push_particle() {
        FContainer container;
        container.push(part_vect[0]);

        assert(container.size() == 1);
        assert(container.getNbParticles() == 1);

        assert(container[0] == part_vect[0].as_tuple() );

        container.push(part_vect[1]);

        assert(container.size() == 2);
        assert(container[1] == part_vect[1].as_tuple() );

        container.push(part_vect[2]);
        assert(container.size() == 3);
        assert(container[2] == part_vect[2].as_tuple() );


        container.push(part_vect[3]);
        assert(container.size() == 4);
        assert(container[3] == part_vect[3].as_tuple() );
    }

    void test_positions() {
        FContainer container;
        for(auto p : part_vect) {
            container.push(p);
        }

        FReal * const * positions = container.getPositions();

        assert(pos_vect[0][0] = positions[0][0]);
        assert(pos_vect[0][1] = positions[1][0]);
        assert(pos_vect[0][2] = positions[2][0]);

        assert(pos_vect[1][0] = positions[0][1]);
        assert(pos_vect[1][1] = positions[1][1]);
        assert(pos_vect[1][2] = positions[2][1]);

        assert(pos_vect[2][0] = positions[0][2]);
        assert(pos_vect[2][1] = positions[1][2]);
        assert(pos_vect[2][2] = positions[2][2]);

    }

    void test_attributes() {
        FContainer container;

        int *  dyn_attributes[5];
        int *  attributes[5];

        dyn_attributes[0] = container.getAttribute(0);
        attributes[0]     = container.getAttribute<0>();
        assert(dyn_attributes[0] == attributes[0]);

        dyn_attributes[1] = container.getAttribute(1);
        attributes[1]     = container.getAttribute<1>();
        assert(dyn_attributes[1] == attributes[1]);

        dyn_attributes[2] = container.getAttribute(2);
        attributes[2]     = container.getAttribute<2>();
        assert(dyn_attributes[2] == attributes[2]);

        dyn_attributes[3] = container.getAttribute(3);
        attributes[3]     = container.getAttribute<3>();
        assert(dyn_attributes[3] == attributes[3]);

        dyn_attributes[4] = container.getAttribute(4);
        attributes[4]     = container.getAttribute<4>();
        assert(dyn_attributes[4] == attributes[4]);

        for(std::size_t i = 0; i < 1000; ++i)
        for(auto p : part_vect) {
            container.push(p);
        }

        assert(container.size() == 4000);

        dyn_attributes[0] = container.getAttribute(0);
        attributes[0]     = container.getAttribute<0>();
        assert(dyn_attributes[0] == attributes[0]);

        dyn_attributes[1] = container.getAttribute(1);
        attributes[1]     = container.getAttribute<1>();
        assert(dyn_attributes[1] == attributes[1]);

        dyn_attributes[2] = container.getAttribute(2);
        attributes[2]     = container.getAttribute<2>();
        assert(dyn_attributes[2] == attributes[2]);

        dyn_attributes[3] = container.getAttribute(3);
        attributes[3]     = container.getAttribute<3>();
        assert(dyn_attributes[3] == attributes[3]);

        dyn_attributes[4] = container.getAttribute(4);
        attributes[4]     = container.getAttribute<4>();
        assert(dyn_attributes[4] == attributes[4]);

    }

};


int main() {
    testFBasicParticleContainer().run();
}
