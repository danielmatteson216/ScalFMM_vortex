#ifndef SCALFMM_BOX_HPP_
#define SCALFMM_BOX_HPP_

#include <ostream>
#include <cmath>

#include "FZCurve.hpp"

/* Implements a N dimensions box
 *
 * \author Quentin Khan <quentin.khan@inria.fr>
 *
 * The box is described by two opposite corners : the maximum and
 * the minimum one. All the class transformations maintain this
 * predicate.
 *
 * \tparam FReal Floating number representation.
 * \tparam Dim Space dimension count.
 * \tparam SpaceFillingCurve A templatize implementation of a space filling curve
 *
 */
template<class _Position, template<std::size_t> class SpaceFillingCurve = FZCurve>
class FBox {
public:
    /// Position type
    using position_t = _Position;
    /// Floating number representation
    using FReal = typename position_t::FReal;
    /// Space dimension
    constexpr static const std::size_t Dim  = position_t::Dim;
    /// Space filling curve type
    using space_filling_curve_t = SpaceFillingCurve<Dim>;

private:
    position_t _c1; ///< Minimum corner
    position_t _c2; ///< Maximum corner
    position_t _center; ///< Center

    /** Rearranges the corners to ensure the maximum-minimum predicate */
    void rearrange_corners () {
        for(std::size_t i = 0; i < Dim; ++i){
            if(_c1[i] > _c2[i]) {
                std::swap(_c1[i], _c2[i]);
            }
        }
        _center = (_c1 + _c2) / 2;
    }
public:
    /// Accessor for the minimum corner
    const position_t& c1() const noexcept {
        return _c1;
    }
    /// Accessor for the maximum corner
    const position_t& c2() const noexcept {
        return _c2;
    }

    /** Builds an empty box at the origin */
    FBox() = default;

    /** Copies an existing box */
    FBox(const FBox&) = default;

    /** Copies an other box */
    FBox& operator=(const FBox& other) = default;

    /** Builds a cube from the lower corner and its side length
     *
     * \param min_corner The lowest corner
     * \param side_length The cube's side length
     **/
    [[gnu::deprecated]] FBox(const position_t& min_corner, FReal side_length) :
        _c1(min_corner),
        _c2(min_corner)
    {
        if(side_length < 0) {
            side_length = - side_length;
        }

        for(auto&& d : _c2) {
            d += side_length;
        }

        _center = (_c1 + _c2) / 2;
    }

    /** Builds a cube using its center and width
     *
     * \param width Cube width
     * \param box_center Cube center
     * \param unused Disambiguisation param from lower corner constructor
     */
    FBox(FReal width, const position_t& box_center)
        : _c1(box_center),
          _c2(box_center),
          _center(box_center)
    {
        if(width < 0) {
            width = - width;
        }
        FReal radius = width / 2;

        for(auto&& d : _c1) {
            d -= radius;
        }

        for(auto&& d : _c2) {
            d += radius;
        }
    }

    /** Builds a box from two corners
     *
     * The maximum and minimum corners are deduced from the given corners.
     *
     * \param corner_1 The first corner.
     * \param corner_2 The second corner.
     */
    FBox(const position_t& corner_1, const position_t& corner_2) : _c1(corner_1), _c2(corner_2) {
        rearrange_corners();
    }

    /** Changes the box corners
     *
     * The maximum and minimum corners are deduced from the given corners.
     *
     * \param corner_1 The first corner.
     * \param corner_2 The second corner.
     */
    void set(const position_t& corner_1, const position_t& corner_2) {
        _c1 = corner_1;
        _c2 = corner_2;

        rearrange_corners();
    }

    /** Checks whether a position is within the box bounds
     *
     * \param p The position to check.
     */
    bool contains(const position_t& p) const {
        for(std::size_t i = 0; i < Dim; i++) {
            if(p[i] < _c1[i] || p[i] > _c2[i]) {
                return false;
            }
        }
        return true;
    }

    /** Checks whether an object's position is within the box bounds
     *
     * The object must implement a `position_t position() const;` method.
     *
     * \tparam T The object type.
     * \param obj The object which position to check.
     */
    template<class T>
    bool contains(const T& obj) const {
        return contains(obj.position());
    }


    /** Accessor for the box center */
    const position_t& center() const {
        return _center;
    }

    /** Accessor for the box corners
     *
     * The corners are numbered using a space filling curve.
     *
     * \param idx The corner index.
     * \return The idx'th corner.
     */
    position_t corner(std::size_t idx) const {
        position_t c;
        std::size_t i = 0;
        for(bool choice : space_filling_curve_t::position(idx)) {
            c[i] = choice ? _c2[i] : _c1[i];
            ++i;
        }
        return c;
    }

    /** Setter for the corners
     *
     * Moves a corner to a new position and modifies the relevant other
     * ones. The corners are numbered using a space filling curve.
     *
     * \param idx The moved corner index.
     * \param pos The new corner position.
     */
    void corner(std::size_t idx, const position_t& pos) {
        std::size_t i = 0;
        for(bool choice : space_filling_curve_t::position(idx)) {
            if(choice) {
                _c2[i] = pos[i];
            } else {
                _c1[i] = pos[i];
            }
            ++i;
        }
        rearrange_corners();
    }

    /** Returns the width for given dimension */
    auto width(std::size_t dim) const noexcept
        -> decltype(std::abs(_c2[dim] - _c1[dim]))
    {
        return std::abs(_c2[dim] - _c1[dim]);
    }

    /** Sums the corners of two boxes */
    FBox operator+(const FBox& other) const {
        return FBox(_c1 + other._c1, _c2 + other._c2);
    }

    /** Tests two boxes equality */
    bool operator==(const FBox& other) const {
        return c1() == other.c1() && c2() == other.c2();
    }
    /** Tests two boxes inequality */
    bool operator!=(const FBox& other) const {
        return ! this->operator==(other);
    }

    friend std::ostream& operator<<(std::ostream& os, const FBox& box) {
        return os << "[" << box.c1() << "," << box.c2() << "]";
    }
};

#endif
