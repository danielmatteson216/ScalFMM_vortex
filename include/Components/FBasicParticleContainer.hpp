// See LICENCE file at project root

#ifndef FBASIC_PARTICLE_CONTAINER_HPP_
#define FBASIC_PARTICLE_CONTAINER_HPP_

#include <vector>

#include "Utils/FPoint.hpp"
#include "Components/FParticleType.hpp"
#include "Components/FBasicParticle.hpp"

#include "Adaptive/FVariadicParticleContainer.hpp"

#include "Utils/variadic_container.hpp"
#include "Utils/FAlignedAllocator.hpp"
#include "inria/integer_sequence.hpp"

namespace scalfmm {
namespace details {
namespace basic_particle_container {

/** \brief Basic particle container template declaration
 *
 * Used to allow automatic indexing of parameter packs during
 * specialization.
 */
template<class Allocator,
         class PosPack,
         class AttributePack,
         class PosIndexList,
         class AttrIndexList,
         class OtherIndexList,
         class... OtherTypes
         >
class FBasicParticleContainerImpl;

/** \brief Basic particle container implementation
 *
 * A container that can store AttributeCount particle attributes. Each attribute
 * has type Attribute.
 *
 * \tparam Allocator Class used to manage memory allocation
 * \tparam FReal Type of the position members
 * \tparam Attribute Type of the attributes
 * \tparam Dim Space dimension count
 * \tparam AttributeCount Number of attributes for a particle
 * (i.e. same type repeated for each existing particle attribute).
 * \tparam posIndices Index pack to index the position components
 * \tparam attrIndices Index pack to index the attributes
 *
 * \note The template paramters are deduced from the specilization
 * arguments.
 *
 */
template<
    class Allocator,
    class FReal,
    class Attribute,
    std::size_t _Dim,
    std::size_t _AttributeCount,
    std::size_t... posIndices,
    std::size_t... attrIndices,
    class... OtherTypes,
    std::size_t... otherIndices
    >
class FBasicParticleContainerImpl<
    Allocator,
    scalfmm::pack<_Dim, FReal>,
    scalfmm::pack<_AttributeCount, Attribute>,
    inria::index_sequence<posIndices...>,
    inria::index_sequence<attrIndices...>,
    inria::index_sequence<otherIndices...>,
    OtherTypes...
    > :
    public FVariadicParticleContainer<
        FBasicParticle<FReal, _Dim, OtherTypes..., scalfmm::pack<_AttributeCount, Attribute> >,
        Allocator
    >
{
public:
    /// Base type
    using FParticle = FBasicParticle<FReal, _Dim, OtherTypes..., scalfmm::pack<_AttributeCount, Attribute> >;

    enum {
        Dim = _Dim,                              ///< Space dimension count
        AttributeCount = _AttributeCount,        ///< Primary attribute count
        OtherCount     = sizeof...(OtherTypes),  ///< Count of secondary attributes
        AttributeBegin = _Dim + sizeof...(OtherTypes), ///< First index of the primary attributes
    };

private:
    using FBase = FVariadicParticleContainer<FParticle, Allocator>;


    /** \brief Proxy array to return positions
     * \warning updated when getPosition is called.
     */
    mutable FReal* _positions[Dim];

public:
    using value_type = FParticle;

    /// Reuse base constructors
    using FBase::FBase;

    /** \brief Get non const array of position components
     *
     *
     * Returns an array of arrays. Each sub-array stores a component of
     * every particle's position. For instance, `getPositions()[0]` is an
     * array containing all the X components of the particles in the
     * container.
     *
     * ~~~~{.cpp}
     * FReal* Xarray = container.getPositions()[0];
     * FReal* Yarray = container.getPositions()[1];
     * FReal* Zarray = container.getPositions()[2];
     * ~~~~
     *
     * \return The array of position member arrays
     *
     * \note This updates the #_positions proxy
     */
    FReal* const* getPositions() {
        auto l = {(_positions[posIndices] = std::get<posIndices>(this->data())) ...};
        (void)l;
        return _positions;
    }

    /** \brief Get const array of position components
     *
     * See getPositions().
     */
    const FReal* const* getPositions() const {
        auto l = {(_positions[posIndices] = const_cast<FReal*>(std::get<posIndices>(this->data()))) ...};
        (void)l;
        return _positions;
    }



    //DELETE, replace with getPositions
    /** \brief Get non const array of position components
     *
     * See getPositions().
     */
    [[gnu::deprecated]]
    FReal* const* getWPositions() {
        return getPositions();
    }


    /** \brief Attribute getter
     *
     * This gets the underlying array that stores the index_th attribute.
     *
     * ~~~~{.cpp}
     * Attribute* attrArray = container.getAttribute(0);
     * ~~~~
     *
     * \param index The index of the attribute array to get.
     *
     * \return  Index_th attribute array
     *
     * \warning For efficicy purpose, it is highly advised to used the
     * template overload of this member: getAttribute().
     */
    Attribute* getAttribute(const int index) {
        using data_type = decltype(this->data());
        auto my_data = this->data();
        return dynamic_attribute<data_type>::get(index, my_data);
    }

    /** \brief Attribute const getter
     * See getAttribute(const int index).
     */
    const Attribute* getAttribute(const int index) const {
        using data_type = decltype(this->data());
        auto my_data = this->data();
        return dynamic_attribute<data_type>::get(index, my_data);
    }

    /** \brief Compilation time attribute getter
     *
     * This gets the underlying array that stores the index_th attribute.
     *
     * ~~~~{.cpp}
     * Attribute* attrArray = container.getAttribute(0);
     * ~~~~
     *
     * \tparam index The index of the attribute array to get.
     *
     * \return  Index_th attribute array
     */
    template<int index>
    Attribute* getAttribute() {
        return std::get<index+AttributeBegin>(FBase::data());
    }

    /** \brief Compilation time attribute const getter
     * See getAttribute().
     */
    template<int index>
    const Attribute* getAttribute() const {
        return std::get<index+AttributeBegin>(FBase::data());
    }

    //DELETE
    /** \brief Pushes several particles in the container
     *
     * \param arraySize Size of the arrays given as other arguments, equals
     * the number of particles to add to the array
     * \param args Pack of arrays, one per attribute
     * \param positionArray Array of positions, one per particle
     *
     * \tparam Args Parameter pack of types convertible to the Attribute
     * type.
     */
    template<typename... Args>
    void pushArray(const FPoint<FReal,Dim>* positionArray, FSize arraySize, Args*... args) {
        for(FSize idx = 0; idx < arraySize; ++idx) {
            this->push(positionArray[idx], args[idx]...);
        }
    }

    using FBase::push;


    /** \brief Push a particle in the container
     *
     * \tparam Args Parameter pack of types convertible to the Attribute type.
     *
     * \note The pack size must be at most AttributeCount, if it is less,
     * the remaining attributes are zero constructed (i.e. `Attribute(0)`
     * is called to build them).
     *
     * \param position The particle position
     * \param args The particle attributes
     */
    // template<typename... Args>
    // void push(const FPoint<FReal>& position, Args... args) {
    //     this->push_helper(
    //         inria::make_index_sequence<AttributeCount+OtherCount-sizeof...(args)>(),
    //         position, args...
    //         );
    // }

    /** \brief Push a particle in the container (array attributes)
     *
     * \param inParticlePosition The particle position
     * \param values std::array containing the particle's attributes
     */
    void push(const FPoint<FReal>& inParticlePosition, const OtherTypes&... others, const std::array<Attribute , AttributeCount>& values) {
        FBase::push_back(std::get<posIndices>(inParticlePosition)..., others..., values[attrIndices]...);
    }

    /** \brief Push a particle in the container (array attributes, type specifier)
     *
     * This method does nothing different thab its simple overload without
     * the FParticleType argument. The argument is unused.
     *
     * \param inParticlePosition Particle position
     * \param Unused
     * \param values std::array containing the particle's attributes
     */
    void push(const FPoint<FReal>& inParticlePosition, const FParticleType /*type*/,
              const std::array<Attribute , AttributeCount>& values) {
        this->push(inParticlePosition, values);
    }

    /** \brief Push a particle in the container (type specifier)
     *
     * This method does nothing different than its simple overload without
     * the FParticleType argument. The argument is unused.
     *
     * \tparam Args Parameter pack for the pushed attributes, automatically
     * deduced from the arguments
     *
     * \param inParticlePosition Particle position
     * \param unnamed Particle type, does nothing
     * \param args Particle attributes
     */
    template<typename... Args>
    void push(const FPoint<FReal>& inParticlePosition, const FParticleType /*particleType*/, Args... args) {
        FBase::push(inParticlePosition, args...);
    }

    template<typename... Args>
    void push(const FPoint<FReal>& inParticlePosition, Args... args) {
        FBase::push(inParticlePosition, args...);
    }

    //DELETE
    /** \brief Removes particles from the container
     *
     * \param indicesToRemove C array containing the indices in the container of
     * the particles to remove
     * \param nbParticlesToRemove Size of the `indicesToRemove array`
     */
    void removeParticles(const FSize indicesToRemove[], const FSize nbParticlesToRemove) {
        std::vector<FSize> removed;
        for(FSize i = 0; i < nbParticlesToRemove; ++i) {
            auto idx = indicesToRemove[i];
            for(auto&& j : removed) {
                if(i < j)
                    --idx;
            }
            removed.push_back(indicesToRemove[i]);
            FBase::erase(this->cbegin()+idx);
        }
    }

    /** \brief Get a pointer to the beginning of the attributes memory block
     *
     * \return A pointer to the beginning of the attributes memory block
     */
    Attribute* getRawData() {
        return this->getAttribute<0>();
    }

    /** \brief Get a pointer to the beginning of the attributes memory block
     *
     * \return A pointer to the beginning of the attributes memory block
     */
    const Attribute* getRawData() const {
        return this->getAttribute<0>();
    }

    /** \brief Get the leading dimension of the data array
     *
     * The data array is allocated as one memory block that is share among
     * attribute subarrays. Since all attributes have the same type, the
     * distance (aka. leading dimension) between the beginning of each subarray is
     * constant.
     *
     * \note The leading dimension may be bigger than `capacity *
     * sizeof(Attribute)` if the memory is aligned.
     *
     * \return The data array leading dimension
     */
    FSize getLeadingRawData() const {
        return std::get<AttributeBegin+1>(this->data()) - std::get<AttributeBegin>(this->data());
    }


    /** \brief Reset all the particle attributes to 0
     *
     * The particle count and positions are not changed.
     */
    void resetToInitialState() {
        auto l = {(this->resetToInitialState(attrIndices), 0)...};
        (void)l;
    }


    /** \brief Reset an attribute subarray to 0
     *
     * The particle count, positions and other attributes are not changed.
     *
     * \param idxAttribute Index of the subarray to reset
     */
    void resetToInitialState(const unsigned idxAttribute) {
        Attribute* p = this->getAttribute(idxAttribute);
        memset(p, 0, sizeof(Attribute) * this->size());
    }

private:

    template <typename AttributeTuple, std::size_t I = AttributeCount - 1 >
    struct dynamic_attribute : dynamic_attribute<AttributeTuple, I-1>  {
        static_assert(I <= AttributeCount, "Attribute index is greater than attribute count.");
        static auto get(const int& index, AttributeTuple data)
            -> typename std::tuple_element<AttributeBegin, AttributeTuple>::type
        {
            if(I == index) {
                return std::get<AttributeBegin+I>(data);
            } else {
                return dynamic_attribute<AttributeTuple, I-1>::get(index, data);
            }
        }
    };

    template<typename AttributeTuple>
    struct dynamic_attribute<AttributeTuple, 0> {
        static auto get(const int& /*index*/, AttributeTuple data)
            -> typename std::tuple_element<AttributeBegin, AttributeTuple>::type
        {
            return std::get<AttributeBegin>(data);
        }
    };


    void update_proxy_arrays() const {
        noop_t(getAttribute(attrIndices)
                = const_cast<Attribute*>(std::get<AttributeBegin+attrIndices>(FBase::data())) ...);
        noop_t(_positions[posIndices]
               = const_cast<FReal*>(std::get<posIndices>(FBase::data())) ...);

    }


};


/** Helper to get FBasicParticleContainer base
 *
 * A bit of meta-programming to find FBasicParticleContainer base.
 */
template<class FReal,
         class Attribute,
         unsigned AttributeCount,
         std::size_t Dim,
         class Allocator,
         class... OtherTypes
         >
struct FBasicParticleContainerChooser {
private:
    /** \brief Returns the base variadic vector for the FBasicParticleContainer
     *
     * This is a declaration that is only meant to be used at compile
     * time to create the right variadic vector.
     *
     * The return type is:
     * ~~~~
     * variadic_vector<FReal,... , Attribute, ...>
     *              Dim times ^                ^ AttributeCount times
     * ~~~~
     *
     * \tparam T The attribute type for the container
     * \tparam dim_indices integer pack for the position indices
     * \tparam attribute_indices integer pack for the attribute indices
     *
     * \param unnamed The position compounds index_sequence
     * \param unnamed The attributes index_sequence
     *
     * \return The class specialization of the container implementation
     */
    template<std::size_t... dim_indices, std::size_t... attribute_indices>
    static auto getContainer(inria::index_sequence<dim_indices...>, inria::index_sequence<attribute_indices...>)
        -> FBasicParticleContainerImpl<
            Allocator,
            scalfmm::pack<sizeof...(dim_indices), FReal>,
            scalfmm::pack<AttributeCount, Attribute>,
            inria::index_sequence<dim_indices...>,
            inria::index_sequence<attribute_indices...>,
            inria::index_sequence_for<OtherTypes...>,
            OtherTypes...
        >;
public:
    /** \brief FBasicParticleContainer base type
    */
    using type = decltype(getContainer(inria::make_index_sequence<Dim>(),
                                       inria::make_index_sequence<AttributeCount>()));
};


} } } // close namespace scalfmm::details::basic_particle_container

/**
 * @author Quentin Khan (quentin.khan@inria.fr)
 *
 * This container which can hold one type (AttributeClass) for each particle.
 * The memory is allocated in one block and stores the positions and the requested type.
 *
 * For example if one wants to store a struct for each particle:
 * \code
 * struct AStruct{ ... };
 * FBasicParticleContainer<FReal, 1, AStruct> container;
 * \endcode
 * And then access is done using:
 * \code
 * AStruct* strucs = container.getAttributes<0>();
 * \endcode
 * For example if one wants to store 4 doubles for each particles:
 * \code
 * FBasicParticleContainer<FReal, 4, double> container;
 * \endcode
 * And then access is done using:
 * \code
 * double* v1 = container.getAttributes<0>();
 * double* v2 = container.getAttributes<1>();
 * double* v3 = container.getAttributes<2>();
 * double* v4 = container.getAttributes<3>();
 * \endcode
 * The memory is aligned to FP2PDefaultAlignement value.
 */
template<class FReal,
         unsigned AttributeCount,
         class Attribute,
         std::size_t Dim = 3,
         class Allocator = FAlignedAllocator<FP2PDefaultAlignement,char>,
         class... OtherTypes
         >
struct FBasicParticleContainer :
    public scalfmm::details::basic_particle_container::FBasicParticleContainerChooser<FReal, Attribute, AttributeCount, Dim, Allocator, OtherTypes... >::type
{};



#endif
