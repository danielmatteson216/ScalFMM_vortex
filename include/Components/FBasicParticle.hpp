/**
 * \brief Implementation file of FBasicParticle class
 * \author Quentin Khan
 *
 * \file
 */

#ifndef FBASIC_PARTICLE_HPP
#define FBASIC_PARTICLE_HPP

#include "inria/integer_sequence.hpp"
#include "Utils/FPoint.hpp"
#include "Utils/FTypePack.hpp"

/**
 * \brief Multi-purpose particle implementation
 *
 * This template implementation of a particle allows simple reuse for several
 * use cases. The aim it to provide an interface that is compatible with the
 * rest of ScalFMM. It is mainly intended to be used as an interface for the
 * particle containers.
 *
 * The Types parameter pack can accept any type that is to be considered as a
 * particle attribute. You can also specify scalfmm::pack type to factorise
 * several types.
 *
 * In the following example, the two specialisations of the class will give the
 * same final structure.
 *
 * ```
 * using FReal = double;
 * static constexpr std::size_t Dim = 3;
 *
 * FBasicParticle<FReal, Dim, int, float, float, float, float>;
 * FBasicParticle<FReal, Dim, int, scalfmm::pack<4, float> >;
 * ```
 *
 * The base of these two classes is
 * ```
 * std::tuple<double, double, double, int, float, float, float, float>;
 * ```
 *
 * \warning Although the classes will have the same final layout, C++ considers
 * these two classes to be different !
 *
 * ##### Example
 *
 * ```
 * // Define a 3D particle with an int attribute
 * using Particle = FBasicParticle<double, 3, int>;
 *
 * Particle p;
 * p.get<>
 * ```
 *
 *
 * \tparam FReal Floating point type
 * \tparam Dim Space dimension count
 * \tparam Types Attributes type list
 *
 */
template<typename _FReal, std::size_t _Dim = 3, typename... Types>
class FBasicParticle : public scalfmm::pack_expand_tuple< scalfmm::pack<_Dim, _FReal>, Types... >{

public:
    /// Storage class : std::tuple<FReal,...(Dim times), Types...>
    using tuple_data_t = scalfmm::
        pack_expand_tuple< scalfmm::pack<_Dim, _FReal>, Types... >;

    /// Expand Types list
    using types_tuple_t = scalfmm::
        pack_expand_tuple< Types... >;

    /// Expand std::enable_if if possible
    template<bool value>
    using sfinae_check = typename std::enable_if<value, bool>::type;

    /**
     * \brief Check parameter pack size vs. attribute + dimension count at
     * compile time
     *
     * \return true if the parameter pack is shorter than the attribute list
     * size + the dimension
     */
    template<typename... Ts>
    constexpr static bool correct_attribute_list_size() {
        return std::tuple_size<tuple_data_t>::value >= sizeof...(Ts) + Dim;
    }

    /// Space dimensions
    constexpr static std::size_t Dim = _Dim;
    /// Size of #tuple_data_t tuple
    constexpr static std::size_t NbAttributes = std::tuple_size<tuple_data_t>::value - Dim;

    /// Floating point type
    using FReal = _FReal;
    /// Position type, required by the FVariadicParticleContainer
    using position_t = FPoint<FReal, Dim>;

    /// Default constructor
    FBasicParticle() = default;
    /// Default copy constructor
    FBasicParticle(const FBasicParticle&) = default;
    /// Default copy operator
    FBasicParticle& operator=(const FBasicParticle&) = default;
    /// Default move constructor
    FBasicParticle(FBasicParticle&&) = default;
    /// Default move operator
    FBasicParticle& operator=(FBasicParticle&&) = default;

    /**
     * \brief Constructor from position and types
     *
     * \tparam Ts Attributes parameter pack; if Ts is too long, this will not
     * compile.
     *
     * \param pos Particle position
     * \param ts Attributes
     *
     * \warning There may be less attributes than defined by the particle, in
     * that case the missing ones are zero constructed.
     */
    template<typename... Ts,
             sfinae_check<correct_attribute_list_size<Ts...>()> = 0>
    FBasicParticle(const position_t& pos, Ts&&... ts) :
        FBasicParticle(inria::make_index_sequence<Dim>(),
                       inria::make_index_sequence<std::tuple_size<tuple_data_t>::value-sizeof...(ts)-Dim>(),
                       pos,
                       std::forward<Ts>(ts)...
            )
    {}

    /// Constructor from tuple equivalent to #tuple_data_t
    template<typename... Ts>
    FBasicParticle(const std::tuple<Ts...>& ts) :
        tuple_data_t(ts)
    {}

    /**
     * \brief Position getter
     *
     * The position is stored in the #_data tuple, to extract it we need to
     * recreate a position_t object. This is done by the position_impl() method.
     *
     * \return A new position_t object
     */
    position_t position() const {
        return position_impl(inria::make_index_sequence<Dim>());
    }

    /**
     * \brief Position setter
     *
     * \parma pos The new position
     */
    void setPosition(const position_t& pos) {
        return setPosition_impl(pos, inria::make_index_sequence<Dim>());
    }

    /**
     * \brief Get a reference to the Ith attribute
     *
     * \tparam I Index of the attribute to get
     */
    template<std::size_t I>
    auto attribute() -> decltype(std::get<Dim+I>(*this)) {
        return std::get<Dim+I>(*this);
    }

    /**
     * \brief Get a const reference to the Ith attribute
     *
     * \tparam I Index of the attribute to get
     */
    template<std::size_t I>
    auto attribute() const -> decltype(std::get<Dim+I>(*this)) {
        return std::get<Dim+I>(*this);
    }


    /**
     * \brief Get a tuple filled with copies of the attributes
     */
    types_tuple_t attributes() const {
        return attributes_impl(
            inria::make_index_sequence<std::tuple_size<types_tuple_t>::value>());
    }

    /**
     * \brief Convert particle to a tuple
     */
    tuple_data_t& as_tuple() {
        return *this;
    }

    /**
     * \brief Convert particle to a tuple
     */
    const tuple_data_t& as_tuple() const {
        return *this;
    }

private:
    /**
     * \brief Contructor implementation
     *
     * Builds the particle, zero constructing missing arguuments.
     */
    template<typename... Ts, std::size_t... Is, std::size_t... Js>
    FBasicParticle(inria::index_sequence<Is...>, inria::index_sequence<Js...>,
                   const FPoint<FReal, Dim>& pos, Ts&&... ts) :
        tuple_data_t(pos[Is]..., std::forward<Ts>(ts)...,
               typename std::tuple_element<Dim+Js+sizeof...(ts), tuple_data_t>::type(0)...)
    {
        //static_assert(sizeof...(Ts) == NbAttributes, "Parameter count is incorrect");
    }

    /**
     * \brief #position method implementation
     *
     * \tparam Is Index sequence of the position elements in the tuple, deduced
     * using arguments
     */
    template<std::size_t... Is>
    position_t position_impl(inria::index_sequence<Is...>) const{
        return position_t(std::get<Is>(this->as_tuple())...);
    }

    /**
     * \brief #setPosition method implementation
     *
     * \tparam Is Index sequence of the position elements in the tuple, deduced
     * using arguments
     *
     * \param pos new position
     */
    template<std::size_t... Is>
    void setPosition_impl(const position_t& pos, inria::index_sequence<Is...>) {
        auto l = {std::get<Is>(this->as_tuple()) = pos[Is] ...};
        (void)l;
    }


    /**
     * \brief #attributes method implementation
     *
     * \tparam Is Index sequence of the attributes in the tuple, deduced using
     * arguments
     *
     * \param pos new position
     */
    template<std::size_t... Is>
    types_tuple_t attributes_impl(inria::index_sequence<Is...>) const {
        return types_tuple_t(std::get<Is+Dim>(*this)...);
    }

};


#endif /* FBASIC_PARTICLE_HPP */
