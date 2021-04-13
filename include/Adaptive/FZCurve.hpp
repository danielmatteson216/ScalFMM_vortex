#ifndef SCALFMM_SPACE_CURVE_HPP_
#define SCALFMM_SPACE_CURVE_HPP_

#include "Utils/FPoint.hpp"

/** Provides the corner traversal order of an N dimension hypercube
 *
 * The positions returned are array of booleans. Each boolean tells where
 * to place the element in the binary grid.
 *
 * For instance, in 2D:
 *
 *
 *         __0__ __1__
 *        |     |     |
 *       0|     |  X  |  pos(X) = [true, false]
 *        |_____|_____|
 *        |     |     |
 *       1|     |  Y  |  pos(Y) = [true, true ]
 *        |_____|_____|
 *
 *
 * \tparam Dim The hypercube dimension.
 */
template<std::size_t _Dim>
class FZCurve {
public:
    /// Space dimension count
    constexpr static const std::size_t Dim = _Dim;
    /// Position type used
    using position_t = FPoint<bool, Dim>;
    /// Position count in the grid
    constexpr static std::size_t pos_count = Fpow(2,Dim);

private:
    /// Array of positions type
    using position_array_t = std::array<position_t, pos_count>;
    /// Array to cache the positions corresponding to indexes
    static const position_array_t _positions;

    /** Creates an array of positions to initialize #_positions */
    static position_array_t create_array() noexcept {
        position_array_t positions;
        for(std::size_t i = 0; i < pos_count; ++i) {
            for(std::size_t j = Dim-1, k = i; k != 0; --j, k >>=1) {
                positions[i][j] = k % 2;
            }
        }
        return positions;
    }

public:
    /** The position corresponding to an index
     *
     * \param idx The index of the point in the space filling curve
     * \return The position corresponding to the space filling curve index
     */
    static position_t position(std::size_t idx) noexcept {
        return _positions[idx];
    }

    /** Index in the space filling curve of a boolean position
     *
     * \param p The position
     * \return The space filling curve index corresponding to the position
     */
    static std::size_t index(const position_t& p) noexcept {
        std::size_t idx = 0;
        for(auto i : p) {
            idx <<= 1;
            idx += i;
        }
        return idx;
    }

    /** Index in the space filling curve of a real position relative to the center of the hypercube
     *
     * \param p The position
     * \param center The center of the hypercube
     * \return The space filling curve index corresponding to the position
     */
    template<typename FReal>
    static std::size_t index(const FPoint<FReal, Dim>& p, const FPoint<FReal, Dim>& center) noexcept {
        std::size_t idx = 0;
        for(std::size_t i = 0; i < Dim; ++i) {
            idx <<= 1;
            idx += p[i] > center[i];
        }
        return idx;
    }

};

// Initialization of static variable
template<std::size_t _Dim>
const typename FZCurve<_Dim>::position_array_t FZCurve<_Dim>::_positions(FZCurve<_Dim>::create_array());

#endif
