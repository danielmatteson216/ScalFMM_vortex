#ifndef FVARIADICPARTICLECONTAINER_H
#define FVARIADICPARTICLECONTAINER_H

#include "Utils/FGlobal.hpp"

#include "Utils/FPoint.hpp"
#include "Utils/variadic_container.hpp"
#include "Utils/FAlignedAllocator.hpp"


namespace scalfmm {
namespace details {

/**
 * \brief Helper to generate FVariadicParticleContainer base class
 *
 * \copydetails FVariadicParticleContainer
 */
template<typename Particle, class Allocator>
struct FVariadicParticleContainerBase {

    /// Expand a tuple types and get the corresponding variadic_vector
    template<typename... Ts>
    static auto expand(const std::tuple<Ts...>&)
        -> variadic_vector<Allocator, Ts...>;

    /// Type of the base class
    using type = decltype(expand(std::declval<Particle>().as_tuple()));
};

}} // close namespace scalfmm::details

/**
 * \brief Particle container defined from a structure
 *
 * transparently stores particles in a structure of arrays manner instead of an
 * array of structures.
 *
 * \tparam Particle  The particle descriptor to store
 * \tparam Allocator The allocator to use
 *
 * The Particle type must define an `as_tuple` method that returns a tuple
 * representation of the particle. It is expected to hold position values before
 * attribute values.
 *
 * ~~~
 * posX, posY, poZ, attributes...
 * ~~~
 *
 *
 * *Example*:
 *
 * ~~~{.cpp}
 * struct particle_t {
 *     double x, y;
 *     float phi;
 *
 *     auto as_tuple() const {
 *         return std::make_tuple(x, y, phi);
 *     }
 * };
 *
 * FVariadicparticlecontainer<particle_t> container;
 * ~~~
 *
 */
template<typename Particle, class Allocator = FAlignedAllocator<128,char> >
class FVariadicParticleContainer
    : public scalfmm::details::FVariadicParticleContainerBase<Particle, Allocator>::type
{
    using FBase = typename scalfmm::details::FVariadicParticleContainerBase<Particle, Allocator>::type;

public:

    /// Alias to the particle type
    using particle_t = Particle;

    // Inherit contructors
    using FBase::FBase;

    /**
     * \brief Get particle count in the container
     *
     * \return Particle count
     */
    FSize getNbParticles() const {
        return FSize(FBase::size());
    }

    //DELETE
    /**
     * \brief Clears the container
     */
    [[gnu::deprecated("Use clear() method instead")]]
    void resetNumberOfParticles() {
        FBase::clear();
    }

    /**
     * \brief Push a particle in the container
     *
     * \param particle Particle to push
     */
    void push(const Particle& particle) {
        this->push_back(particle.as_tuple());
    }

    /**
     * \brief Push a particle in the container
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
    template<typename... Args>
    void push(const typename particle_t::position_t& position, Args... args) {
        using vt = typename FBase::value_type;

        constexpr auto Dim = Particle::position_t::Dim;
        constexpr int vt_size = std::tuple_size<vt>::value;
        constexpr int missing_offset = Dim + sizeof...(args);
        constexpr int missing_size = vt_size - missing_offset;

        this->push_impl(
            inria::make_index_sequence<Dim>{},
            inria::make_index_sequence<missing_size>{},
            position, args...);
    }

  template<typename Position, typename... Args, 
	   std::size_t... Is,
	   std::size_t... Js> //, 	   typename... Ts >
  void push_impl(inria::index_sequence<Is...>,
		 inria::index_sequence<Js...>,
		 const Position& position,
		 Args... args)
  {
    using vt = typename FBase::value_type;
    constexpr auto Dim = Particle::position_t::Dim;
    constexpr int missing_offset = Dim + sizeof...(args);

    this->push_back(std::get<Is>(position)..., args...,
		    std::tuple_element_t<missing_offset + Js, vt>{}...);
  }

    /**
     * \brief Get size to send the container
     *
     * \return The size of the serialised particle container in bytes
     *
     * \note Implementation is #getSavedSize_impl
     */
    FSize getSavedSize() const {
        return this->getSavedSize_impl(this->data());
    }
private:
    /**
     * \brief Implementation of #getSavedSize
     *
     * \tparam Types The container types
     *
     * \return The size of the serialised particle container in bytes
     */
    template<class... Types>
    FSize getSavedSize_impl(std::tuple<Types*...>) const {
        FSize tmp = sizeof(FSize);
        // Sum, for T in Types, of (sizeof(T) * this->size())
        auto l = {(tmp += sizeof(Types) * this->size(), 0)..., 0};
        (void)l;
        return tmp;
    }

public:

    /**
     * \brief Serialize particle container to buffer
     *
     * \tparam BufferWriter Buffer class
     *
     * \param buffer Buffer in which to save the container
     *
     * \note Implementation is #save_impl
     */
    template<class BufferWriter>
    void save(BufferWriter& buffer) const {
        this->save_impl(buffer, FBase::indices);
    }
private:
    /**
     * \brief Implementation of #save
     *
     * \tparam BufferWriter Buffer class
     * \tparam Is Index pack for the container types
     *
     * \note Writes #getSavedSize bytes into the buffer
     */
    template<class BufferWriter, std::size_t ... Is>
    void save_impl(BufferWriter& buffer, inria::index_sequence<Is...>) const {
        buffer << FSize(this->size());
        // For each sub-array `s`, write `s` to the buffer
        auto l = {
            (buffer.write(std::get<Is>(this->data()), this->size()),0)...,
            0};
        (void)l;
    }

public:
    /** \brief Restore container content from a buffer
     *
     * \tparam BufferReader Buffer class
     *
     * \param buffer Buffer from which to restore the buffer
     *
     * \note Implementation is #restore_impl
     */
    template<class BufferReader>
    void restore(BufferReader& buffer) {
        this->restore_impl(buffer, FBase::indices);
    }
private:
    /**
     * \brief Implementation of #restore
     *
     * \tparam BufferReader Buffer class
     * \tparam Is Index pack for the container types
     *
     * \note Reads #getSavedSize bytes from the buffer
     */
    template<class BufferReader, std::size_t ... Is>
    void restore_impl(BufferReader& buffer, inria::index_sequence<Is...>) {
        FSize particle_count = 0;
        buffer >> particle_count;

        this->resize(particle_count);

        auto l = {
            (buffer.fillArray(std::get<Is>(this->data()), this->size()),0)...,
            0};
        (void)l;
    }


};

#endif /* FVARIADICPARTICLECONTAINER_H */
