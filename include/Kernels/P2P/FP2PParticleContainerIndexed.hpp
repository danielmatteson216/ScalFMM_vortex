// See LICENCE file at project root
#ifndef FP2PPARTICLECONTAINERINDEXED_HPP
#define FP2PPARTICLECONTAINERINDEXED_HPP

#include "Containers/FVector.hpp"
#include "Components/FParticleType.hpp"

#include "FP2PParticleContainer.hpp"

template<class FReal, int NRHS = 1, int NLHS = 1, int NVALS = 1>
class FP2PParticleContainerIndexed : public FP2PParticleContainer<FReal, NRHS,NLHS,NVALS, FSize> {
    typedef FP2PParticleContainer<FReal, NRHS,NLHS,NVALS, FSize> Parent;

    mutable FVector<FSize> indexes{};

public:

    const FVector<FSize>& getIndexes() const {
        indexes.memocopy(const_cast<FSize*>(std::get<3>(this->data())), this->size());
        return indexes;
    }

};

#endif // FP2PPARTICLECONTAINERINDEXED_HPP
