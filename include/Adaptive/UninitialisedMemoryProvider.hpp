
#ifndef UNINITIALISEDMEMORYPROVIDER_HPP
#define UNINITIALISEDMEMORYPROVIDER_HPP

#include <cmath>
#include <memory>
#include <sstream>

/**
 * \brief Basic item provider
 *
 * A basic uninitialised memory provider. Provides heap allocation of a constant
 * size array of objects. No reallocation/release can happen apart from deleting
 * the provider.
 *
 * \tparam T Type of the object to provide memory for
 */
template<class T>
class UninitialisedMemoryProvider {
public:
    using type = T;

private:
    using uninit_item = typename std::aligned_storage<sizeof(T)>::type;
    const std::size_t _capacity = 1<<3;
    std::size_t _size = 0;
    std::unique_ptr<uninit_item[]> _data =
        std::unique_ptr<uninit_item[]>{new uninit_item[this->_capacity]};

    std::size_t& size() noexcept {
        return this->_size;
    }

    /**
     * \brief Accessor to the underlying data
     */
    type* data() noexcept {
        return reinterpret_cast<type*>(this->_data.get());
    }

public:
    UninitialisedMemoryProvider() = delete;
    UninitialisedMemoryProvider(const UninitialisedMemoryProvider&) = delete;
    UninitialisedMemoryProvider& operator=(const UninitialisedMemoryProvider&) = delete;
    UninitialisedMemoryProvider(UninitialisedMemoryProvider&&) = default;
    UninitialisedMemoryProvider& operator=(UninitialisedMemoryProvider&&) = default;

    /**
     * \brief Create a provider that may contain count items
     *
     * \param count Bucket capacity
     */
    UninitialisedMemoryProvider(std::size_t count) : _capacity(count) {}

    /**
     * \brief Total capacity of the provider
     */
    std::size_t capacity() const noexcept {
        return this->_capacity;
    }

    /**
     * \brief Already provided object count
     */
    std::size_t size() const noexcept {
        return this->_size;
    }

    /**
     * \brief Read-only accessor to the underlying data
     */
    const type* data() const noexcept {
        return reinterpret_cast<type*>(this->_data.get());
    }

    /**
     * \brief Provide memory for given object count
     *
     * Checks that the provider can issue enough memory for the required object
     * count and returns a pointer to the first one.
     *
     * \param count Required object count
     *
     * \return An array of uninitialised objects
     *
     * \exception std::length_error When the requested count cannot be provided
     * (capacity - size < count)
     */
    type* provide(const std::size_t count) {
        if(! this->can_provide(count)) {
            std::stringstream sstr;
            sstr << "Cannot provide enough items, asked for ";
            sstr << count << ", got " << this->capacity() - this->size();
            throw std::length_error(sstr.str());
        }
        return unsafe_provide(count);
    }


    /**
     * \brief Checks whether the provider has enought memory left or not
     *
     * \param count Object count
     *
     * \return True if at least memory for count items is available, false
     * otherwise.
     */
    bool can_provide(const std::size_t count) const noexcept {
        return this->size() + count <= this->capacity();
    }

    /**
     * \brief Provide memory for given object count without checks
     *
     * This version of the #provide method does not check nor throw in case or
     * error, it is up to the user to ensure that enough memory is available
     * before calling it.
     *
     * The provider size is updated accordinly to the request, therefore it is
     * possible to check after the call is size <= capacity.
     *
     * \param count The required object count
     *
     * \return An array of uninitialised objects
     */
    type* unsafe_provide(const std::size_t count) noexcept {
        auto ret = this->data() + this->size();
        this->_size += count;
        return ret;
    }

};


#endif /* UNINITIALISEDMEMORYPROVIDER_HPP */
