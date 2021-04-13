// See LICENCE file at project root
#ifndef FTESTCELL_HPP
#define FTESTCELL_HPP

#include <cstddef>
#include <ostream>
#include "FBasicCell.hpp"

#include "Kernels/Generic/FGenericData.hpp"


namespace scalfmm {
namespace detail {
namespace FTestCell_impl {

    template<class Tag>
    struct exp_impl {
        using type = long long int;

        type data;

        operator type() const {
            return data;
        }

        type get() const {
            return this->data;
        }
        void set(const type val) {
            this->data = val;
        }
        void reset() {
            this->data = 0;
        }
        template<class BufferWriter>
        void save(BufferWriter& buffer) {
            buffer << this->data;
        }
        template<class BufferReader>
        void restore(BufferReader& buffer) {
            buffer >> this->data;
        }
        FSize getSavedSize() const {
            return sizeof(this->data);
        }
        template <class BufferWriterClass>
        void serialize(BufferWriterClass& buffer) const {
            buffer << this->data;
        }
        /** Deserialize only up data in a buffer */
        template <class BufferReaderClass>
        void deserialize(BufferReaderClass& buffer){
            buffer >> this->data;
        }

        friend std::ostream& operator<<(std::ostream& os, const exp_impl& d) {
            return (os << d.get());
        }
    };

    using multipole_t = scalfmm::detail::FTestCell_impl::exp_impl<class MultipoleTag>;
    using local_expansion_t = scalfmm::detail::FTestCell_impl::exp_impl<class LocalExpansionTag>;

}
}
}


using FTestCell = FGenericData<
    scalfmm::detail::FTestCell_impl::multipole_t,
    scalfmm::detail::FTestCell_impl::local_expansion_t
    >;

#endif //FTESTCELL_HPP
