#ifndef VARIADIC_CONTAINER
#define VARIADIC_CONTAINER

#include <cassert>
#include <initializer_list>
#include <limits>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>

#include "FOstreamTuple.hpp"

#include "inria/integer_sequence.hpp"

/// To discard fold expression results
struct noop_t {
    template<typename... Types>
    noop_t(const Types&...) {}
};

/// Check whether allocator can rebind to multiple type at the same time
template<typename T, typename... Types>
class check_multi_allocator_rebind {

    template<typename...Ts, typename T::template rebind<Ts...>::other* = nullptr>
    static typename std::conditional<
        (sizeof...(Types) > 1),
        std::true_type,
        std::false_type>
    ::type
    check(std::nullptr_t);

    template<typename...Ts>
    static std::false_type check(...);

public:
    constexpr static bool value = std::is_same<std::true_type, decltype(check<Types...>(nullptr))>::value;
};

template<typename Alloc, typename... Ts>
using has_multi_rebind = typename std::enable_if<check_multi_allocator_rebind<Alloc, Ts...>::value>::type*;

template<typename Alloc, typename... Ts>
using lacks_multi_rebind = typename std::enable_if<!check_multi_allocator_rebind<Alloc, Ts...>::value>::type*;

#ifdef VARIADIC_VECTOR_DEBUG
class utest_variadic_vector;
#endif

// Forward declaration
template<class Allocator, typename TypeTuple, typename IndexList>
class variadic_vector_impl {};

// Forward declaration
template<typename TypeTuple, typename IndexList>
class variadic_vector_iterator {};

/** \brief A variadic vector
 *
 * \tparam Allocator An allocator class such as `Allocator::template
 * rebind<T>::other` is an allocator for type T.
 * \tparam Types A list of types, a vector is created for each type.
 *
 * TODO description
 *
 * \note
 *     - See http://en.cppreference.com/w/cpp/concept/Allocator for
 * information on allocators.
 *     - See http://en.cppreference.com/w/cpp/language/parameter_pack for
 * information on parameter packs in variadic templates.
 */
template<class Allocator, typename... Types, std::size_t... Indices>
class variadic_vector_impl<Allocator, std::tuple<Types...>, inria::index_sequence<Indices...> > {

    #ifdef VARIADIC_VECTOR_DEBUG
    friend class utest_variadic_vector;
    #endif

public:

    using value_type = std::tuple<Types...>;

    using pointer_tuple = std::tuple<Types*...>;
    using const_pointer_tuple = std::tuple<const Types*...>;

    using reference_tuple  = std::tuple<Types&...>;
    using const_reference_tuple  = std::tuple<const Types&...>;

    using rvalue_tuple = std::tuple<Types&&...>;

    using iterator = variadic_vector_iterator<std::tuple<Types...>, inria::index_sequence<Indices...>>;
    using const_iterator = variadic_vector_iterator<std::tuple<const Types...>, inria::index_sequence<Indices...>>;

    using difference_type = std::ptrdiff_t;
    using size_type = std::size_t;

    using allocator_type_tuple = std::tuple<typename Allocator::template rebind<Types>::other...>;

    constexpr static auto indices = inria::make_index_sequence<sizeof...(Types)>();


private:

    template<typename T>
    using result_tuple = std::tuple<decltype((std::declval<Types*>(), std::declval<T>()))...>;


    /// Allocator for the data block
    typename Allocator::template rebind<char>::other _data_allocator;

    /// Current element count
    size_type _size = 0;
    /// Current element capacity
    size_type _capacity = 0;

    /// Tuple of allocators
    allocator_type_tuple _allocator_tuple;

    /** \brief Tuple of arrays
     * Each element in the tuple is a pointer to an address between _data and
     * (_data + _size)
     */
    pointer_tuple _data_tuple = pointer_tuple{(Types*)nullptr...}; // initialized to nullptr

public:

    // Constructors
    explicit variadic_vector_impl(const allocator_type_tuple& alloc)
        : _allocator_tuple(alloc)
    {}

    /// Default constructor
    explicit variadic_vector_impl() = default;

    /** \brief Constructor, fill with copies
     *
     * Fills the vector with `count` copies of `value` after construction.
     *
     * \param count Number of copies
     * \param value Value to copy into the vector
     * \param alloc Allocator that must be used to manage memory
     *
     * \return The new variadic_vector
     */
    variadic_vector_impl(size_type count,
                         const value_type& value,
                         const allocator_type_tuple& alloc = allocator_type_tuple())
        : _allocator_tuple(alloc)
    {
        this->reserve(count);
        for(size_type i = 0; i < count; ++i) {
            push_back(value);
        }
    }

    /** \brief Constructor, count default construction
     *
     * Fills the vector with `count` default constructed elements.
     *
     * \param count Number of elements to insert
     *
     * \return The new variadic_vector
     */
    explicit variadic_vector_impl(size_type count) {
        this->reserve(count);
        this->_size = count;
        for(size_type i = 0; i < count; ++i) {
            noop_t{(std::get<Indices>(_allocator_tuple).construct(std::get<Indices>(_data_tuple)+i), 0)... };
        }
    }

    /** \brief Constructor from range
     *
     * Fills the vector with content from given range.
     *
     * \tparam InputIt Iterator type that can be converted to #value_type
     *
     * \param first TODO describe
     * \param last TODO describe
     * \param alloc TODO describe
     *
     * \return
     */
    template<class InputIt>
    variadic_vector_impl(InputIt first,
                         InputIt last,
                         const allocator_type_tuple& alloc = allocator_type_tuple())
        : _allocator_tuple(alloc)
    {
        this->assign(first, last);
    }

    variadic_vector_impl( variadic_vector_impl& other) {
        this->assign(other.begin(), other.end());
    }

    variadic_vector_impl(variadic_vector_impl&& other) {
        this->swap(other);
    }

    // No initializer list constructor
    variadic_vector_impl(std::initializer_list<value_type>) = delete;

    ~variadic_vector_impl() {
        this->clear();
        if(nullptr != std::get<0>(this->data())) {
            this->deallocate(this->_data_tuple, this->_capacity);
        }
    }

    /** \brief Assignment operator
     *
     * \param other Container to copy
     *
     * \return The new container
     */
    variadic_vector_impl& operator=(const variadic_vector_impl& other) {
        this->clear();
        this->assign(other.begin(), other.end());
        return *this;
    }

    /** \brief Function description
     *
     * \param other Container to move
     *
     * \return The new container
     */
    variadic_vector_impl& operator=(variadic_vector_impl&& other) {
        this->swap(other);
        return *this;
    }

    // No initializer list
    variadic_vector_impl& operator=(std::initializer_list<value_type>) = delete;


    /** get allocators used by vector */
    /** \brief Get allocators used by the vector
     *
     * The vector handles several different types, each of which needs a
     * different allocator.
     *
     * \return A tuple containing the allocators corresponding to the Types
     */
    allocator_type_tuple get_allocator() const {
        return allocator_type_tuple();
    }


    // ELEMENT ACCESS

    /** \brief Access specified element with bound checking
     *
     * \param pos Element position
     *
     * \return A tuple of references to the elements corresponding to each Type
     * at position `pos` in each sub-array
     *
     * \exception std::out_of_range if `pos >= size()`
     */
    reference_tuple at(size_type pos) {
        if(! (pos < size())) {
            throw std::out_of_range("variadic_vector::at: Position not in range");
        }
        return reference_tuple(std::get<Indices>(_data_tuple)[pos]...);
    }

    /** \brief Access specified element with bound checking
     *
     * \param pos Element position
     *
     * \return A tuple of const references to the elements corresponding to each Type
     * at position `pos` in each sub-array
     *
     * \exception std::out_of_range if `pos >= size()`
     */
    const_reference_tuple at(size_type pos) const {
        if(! (pos < size())) {
            throw std::out_of_range("variadic_vector::at: Position not in range");
        }
        return const_reference_tuple(std::get<Indices>(_data_tuple)[pos]...);
    }


    /** \brief Access specified element
     *
     * \param pos Element position
     *
     * \return A tuple of references to the elements corresponding to each Type
     * at position `pos` in each sub-array
     */
    reference_tuple operator[](size_type pos) {
        return reference_tuple(std::get<Indices>(_data_tuple)[pos]...);
    }

    /** \brief Access specified element
     *
     * \param pos Element position
     *
     * \return A tuple of const references to the elements corresponding to each Type
     * at position `pos` in each sub-array
     */
    const_reference_tuple operator[](size_type pos) const  {
        return const_reference_tuple(std::get<Indices>(_data_tuple)[pos]...);
    }

    /** \brief Acces the first element
     *
     * \return A tuple of references to the first elements in the vector
     *
     * \note Behaviour is undefined if the vector is empty
     */
    reference_tuple front() {
        return operator[](0);
    }

    /** \brief Access the first element
     *
     * \return A tuple of const references to the first elements in the vector
     *
     * \note Behaviour is undefined if the vector is empty
     */
    const_reference_tuple front() const {
        return operator[](0);
    }

    /** \brief Access the last element
     *
     * \return A tuple of references to the last elements in the vector
     *
     * \note Behaviour is undefined if the vector is empty
     */
    reference_tuple back() {
        return operator[](size()-1);
    }

    /** \brief Access the last element
     *
     * \return A tuple of const references to the last elements in the vector
     *
     * \note Behaviour is undefined if the vector is empty
     */
    const_reference_tuple back() const {
        return operator[](size()-1);
    }

    /** \brief Direct access to the underlying arrays
     *
     * The variadic vector maintains several sub-arrays in parallel, this allows
     * a direct acces to the subarrays in a C-like fashion.
     *
     * \return A tuple of pointers to the sub-arrays
     */
    pointer_tuple data() noexcept {
        return _data_tuple;
    }

    /** \brief Direct access to the underlying arrays
     *
     * The variadic vector maintains several sub-arrays in parallel, this allows
     * a direct acces to the subarrays in a C-like fashion.
     *
     * \return A tuple of const pointers to the sub-arrays
     */
    const_pointer_tuple data() const noexcept {
        return const_pointer_tuple(_data_tuple);
    }

    // ITERATORS


    /** \brief Get an iterator to the vector beginning
     *
     * \return An iterator to the vector first element
     */
    iterator begin() {
        return iterator(_data_tuple);
    }

    /** \brief Get a constant iterator to the vector beginning
     *
     * \return A constant iterator to the vector first element
     */
    const_iterator begin() const {
        return const_iterator(_data_tuple);
    }

    /** \brief Get a constant iterator to the vector beginning
     *
     * \return A constant iterator to the vector first element
     */
    const_iterator cbegin() const {
        return begin();
    }

    /** \brief Get an iterator to the vector end
     *
     * \return An iterator past the vector last element
     */
    iterator end() {
        return iterator(_data_tuple) + static_cast<difference_type>(this->size());
    }

    /** \brief Get a constant iterator to the vector end
     *
     * \return A constant iterator past the vector last element
     */
    const_iterator end() const {
        return const_iterator(_data_tuple) + static_cast<difference_type>(this->size());
    }

    /** \brief Get a constant iterator to the vector end
     *
     * \return A constant iterator past the vector last element
     */
    const_iterator cend() const {
        return end();
    }


    // CAPACITY

    /** \brief Check whether the container is empty
     *
     * \return true if the container is empty, false otherwise
     */
    bool empty() const noexcept {
        return (size() == 0);
    }

    /** \brief Get the number of elements in the vector
     *
     * \return The count of elements stored in the vector
     */
    size_type size() const noexcept {
        return _size;
    }


    /** \brief Get the maximum possible number of elements
     *
     * \return The maximum possible size of the vector
     */
    size_type max_size() const noexcept {
        return std::numeric_limits<size_type>::max();
    }


private:
    /** \brief Moves data from origin to destination sets of arrays
     *
     * \warning origin and destination must not overlap !
     *
     * \param origin Tuple of pointers to the data to be moved
     * \param origin_size Size of the origin arrays
     * \param destination Tuple of pointers to the new memory
     * \param destination_size Size of the destination arrays
     */
    void move_data(pointer_tuple origin, size_type origin_size, pointer_tuple destination, size_type destination_size) {
        noop_t{
            (move_data_impl(std::get<Indices>(origin),
                            origin_size,
                            std::get<Indices>(destination),
                            destination_size),
             0) ...};
    }

    /** \brief Move data from origin to destination array
     *
     * \warning origin and destination must not overlap !
     *
     * \param origin Pointer to the data to be moved
     * \param origin_size Size of the origin arrays
     * \param destination Pointer to the new memory
     * \param destination_size Size of the destination arrays
     */
    template<typename T>
    void move_data_impl(T* origin, size_type origin_size, T* destination, size_type destination_size) {
        size_type limit = origin_size < destination_size ? origin_size : destination_size;
        typename Allocator::template rebind<T>::other alloc;
        for(std::size_t i = 0; i < limit; ++i, ++destination, ++origin) {
            alloc.construct(destination, std::move(*origin));
        }
    }


    /** \brief Destroy data from a tuple of arrays
     *
     * \param data_tuple Tuple of pointers to the data to be destroyed
     *
     * \note data_tuple arrays are supposed to have a length of #size.
     */
    void destroy_data(pointer_tuple data_tuple) {
        noop_t{(this->array_destroy_data(std::get<Indices>(data_tuple)), 0) ...};
    }

    /** \brief Destroy data from an array
     *
     * \tparam T Type of the data pointer
     *
     * \param data_ptr Pointer to the data to be destroyed
     *
     * \note data_ptr array is supposed to have a length of #size.
     */
    template<typename T>
    void array_destroy_data(T* data_ptr) {
        typename Allocator::template rebind<T>::other alloc;
        for(std::size_t i = 0; i < this->size(); ++i, ++data_ptr) {
            alloc.destroy(data_ptr);
        }
    }


    /** \brief Allocate memory when a multi-allocator is available
     *
     * The vector contains several sub-arrays of types Types...
     *
     * One may write a custom allocator to manage the way the memory is
     * allocated, for example to guarantee that the memory is allocated as a
     * unique chunk. To do this, a custom allocator must provide a public
     * `rebind` structure that can take a parameter pack as template parameter.
     *
     * The multi-allocator must provide the following interface :
     * ~~~~{.cpp}
     * template<typename... Ts>
     * class allocator {
     *     template<typename... Us>
     *     struct rebind {
     *         using other = <your multi-allocator>
     *     };
     *
     *     /// Allocate space for sizeof...(Ts) sub-arrays of size element_count.
     *     /// \return a tuple of pointers to the sub-arrays for each type
     *     std::tuple<Ts*...> allocate(size_type element_count) };
     *
     *     /// Free data_tuple sub-arrays memory
     *     void deallocate(std::tuple<Ts*...> data_tuple, size_type data_element_count)
     * };
     * ~~~~
     *
     * This template overload will only exist when such an allocator is given as
     * the Allocator template parameter
     *
     * \tparam Ts Types to rebind the Allocator to. Normally equivalent to the
     * class parameter pack Types. This indirect definition allows using SFINAE.
     * \tparam has_multi_rebind SFINAE parameter to disable this method if the
     * Allocator is not a multi allocator.
     *
     * \param new_data Tuple of pointers to set to the newly allocated sub-arrays
     * \param new_cap New sub-arrays capacity
     */
    template<typename... Ts, has_multi_rebind<Allocator, Ts...> = nullptr>
    void allocate(pointer_tuple& new_data, size_type& new_cap) {
        // get multi-allocator
        typename Allocator::template rebind<Types...>::other alloc;
        new_data = alloc.allocate(new_cap);
    }

    /** \brief Deallocate memory when a multi-allocator is available
     *
     * \tparam Ts Types to rebind the allocator to
     * \tparam has_multi_rebind SFINAE parameter to take out this method if the
     * Allocator is not a multi-allocator.
     *
     * \param data_tuple Tuple of pointers to the sub-arrays
     * \param data_capacity Sub-arrays capacity
     */
    template<typename... Ts, has_multi_rebind<Allocator, Ts...> = nullptr>
    void deallocate(pointer_tuple& data_tuple, size_type& data_capacity) {
        // get multi-allocator
        typename Allocator::template rebind<Types...>::other alloc;
        alloc.deallocate(data_tuple, data_capacity);
    }


    /** \brief Variadic sum
     *
     * \tparam Ts Operands types
     *
     * \param ts Operands
     *
     * \return The sum of ts...
     */
    template<typename... Ts>
    auto sum(const Ts&... ts) -> typename std::common_type<Ts...>::type {
        typename std::common_type<Ts...>::type res = {};
        auto l = {(res += ts, 0)...}; (void)l;
        return res;
    }


    /** \brief Set the sub-array tuple from a standard allocation
     *
     * \note Only used when the Allocator is not a multi-allocator.
     *
     * This function works recursively. The first call is assumed t
     *
     *
     * This overload will set the J_th element
     * of `data_tuple` to `data_ptr` and will the advance it to the end of the
     * sub-array for a subsequent call.
     *
     * \tparam I Current sub-array index
     * \tparam J Next sub-array index
     * \tparam Is Next to next sub-array indices
     *
     * \param data_ptr Pointer to available memory
     * \param data_tuple Data tuple to set, this will set the J-th pointer
     * \param data_capacity Sub-array expected capacity
     */
    template<size_type I, size_type J, size_type... Is>
    void set_data_tuple(char*& data_ptr, pointer_tuple& data_tuple, const size_type& data_capacity) {
        assert(I != 0 || ((char*)std::get<0>(data_tuple)) == data_ptr);
        data_ptr += sizeof(typename std::tuple_element<I, value_type>::type) * data_capacity;
        std::get<J>(data_tuple) = reinterpret_cast<
            typename std::tuple_element<J, pointer_tuple>::type>(data_ptr);
        set_data_tuple<J, Is...>(data_ptr, data_tuple, data_capacity);
    }


    /** \brief Set the sub-array recursion end
     *
     * \note Only used when the Allocator is not a multi-allocator.
     *
     * This overload is a no-op, used to end the recursion of the other
     * overloads.
     *
     * \tparam I Unused, should be the index of the last element of the class
     * parameter pack Types
     *
     * \param unnamed unused
     * \param unnamed unused
     * \param unnamed unused
     */
    template<size_type I>
    void set_data_tuple(char*&, pointer_tuple&, const size_type&) {}


    /** \brief Allocate memory when only a standard allocator is available
     *
     * \note Only available when the Allocator is not a multi-allocator.
     *
     * This version of the allocate member uses the Allocator to get enough
     * space to store the sizeof...(Types) sub-arrays.
     *
     * The new allocated memory is then split using set_data_tuple().
     *
     * \tparam Ts Parameter pack of types to allocate memory for. Normally
     * equivalent to the class parameter pack Types. This indirect definition
     * allows using SFINAE.
     * \tparam lacks_multi_rebind SFINAE type used to disable this method if the
     * Allocator is a multi-allocator.
     *
     * \param new_data Pointer tuple to set to the newly allocated memory
     * \param new_cap  Number of elements to allocate memory for
     * (ie. element-wise size of each sub-array)
     */
    template<typename... Ts, lacks_multi_rebind<Allocator, Ts...> = nullptr>
    void allocate(pointer_tuple& new_data, size_type& new_cap) {
        // allocate data
        char* memory = _data_allocator.allocate(new_cap * sum(sizeof(Ts)...));
        // setup pointer tuple to point to the right parts of the vector
        new_data = std::make_tuple((Ts*) memory...);
        this->set_data_tuple<Indices...>(memory, new_data, new_cap);
    }

    /** \brief Deallocate memory when only a standard allocator is available
     *
     * \note Only available when the Allocator is not a multi-allocator.
     *
     * \tparam Ts Types to rebind the allocator to
     * \tparam lacks_multi_rebind SFINAE type used to disable this method if the
     * Allocator is a multi-allocator.
     *
     * \param data_tuple Tuple of pointers to the sub-arrays
     * \param data_capacity Sub-arrays capacity
     */
    template<typename... Ts, lacks_multi_rebind<Allocator, Ts...> = nullptr>
    void deallocate(pointer_tuple& data_tuple, size_type& data_capacity) {
        _data_allocator.deallocate(
            reinterpret_cast<char*>(std::get<0>(data_tuple)),
            data_capacity  * sum(sizeof(Types)...)
            );
    }

    /** \brief Set sub-array memory to 0
     *
     * The memory has to be set to 0 when using vectorial instruction sets,
     * otherwise some calculations may end up false.
     *
     * \param data_tuple Tuple of pointers to the sub-arrays
     * \param cap Sub-arrays capacity
     */
    void set_to_zero(pointer_tuple data_tuple, size_type cap) {
        auto l = { (memset(std::get<Indices>(data_tuple), 0, sizeof(Types) * cap), 0) ...};
        (void)l;
    }


    /** \brief Reallocates memory unconditionnaly to new_cap
     *
     * If new cap is smaller than the current size, the overflowing elements are
     * destroyed.
     *
     * \param new_cap The new capacity
     *
     * \tparam has_multi_rebind SFINAE check to choose the right reallocate
     * version.
     */
    void reallocate(size_type new_cap) {
        #if defined(SCALFMM_USE_AVX)
        {
            std::size_t avx_mult = new_cap % 8;
            if(0 != avx_mult) {
                new_cap += 8 - avx_mult;
            }
        }
        #endif
        // Reserve new storage
        pointer_tuple new_data;
        this->allocate<Types...>(new_data, new_cap);

        #if defined(SCALFMM_USE_AVX)
        set_to_zero(new_data, new_cap);
        #endif

        // Check that we do hold data before attempting to free or copy
        if(nullptr != std::get<0>(this->data())) {
            // Move data to new storage
            this->move_data(this->_data_tuple, this->size(), new_data, new_cap);
            // Free old storage
            this->destroy_data(this->_data_tuple);
            this->deallocate<Types...>(this->_data_tuple, this->_capacity);
        }
        // Set new storage
        this->_data_tuple = new_data;

        if(this->size() > new_cap) {
            this->_size = new_cap;
        }
        this->_capacity = new_cap;
    }

public:

    /** \brief Reserve storage
     *
     * If the capacity() is less than `new_cap`, allocates new memory, otherwise
     * does nothig.
     *
     * \param new_cap New storage capacity
     */
    void reserve(size_type new_cap) {
        if(capacity() < new_cap) {
            this->reallocate(new_cap);
        }
    }


    /** \brief Number of elements that can currently be held in storage
     *
     * \return Current storage capacity
     */
    size_type capacity() const noexcept {
        return _capacity;
    }


    /** \brief Reduces memory usage by freeing memory
     *
     * If capacity() is greater than size(), reallocate memory with size() as
     * the new capacity.
     */
    void shrink_to_fit() {
        this->reallocate(this->size());
    }

    // MODIFIERS

    /** \brief Assigns values to the container by copy
     *
     * The container is cleared and the new values are inserted.
     *
     * \param count  Number of copies to insert
     * \param values Value to copy
     */
    void assign(size_type count, const value_type& values) {
        clear();
        insert(begin(), count, values);
    }

    /** \brief Assigns values to the container from range
     *
     * The container is cleared and the new values inserted
     *
     * \tparam TupleIterator Iterator convertible to a tuple of values that
     * can be converted to elements of #value_type.
     *
     * \param first_it Range beginning
     * \param last_it  Range end
     *
     * \warning Behaviour is undefined if the range belongs to the container.
     */
    template<typename TupleIterator>
    void assign(TupleIterator first_it, TupleIterator last_it) {
        clear();
        this->reserve(last_it - first_it);
        insert(begin(), first_it, last_it);
    }


    /** \brief Assigns values to the container from initializer_list
     *
     * The container is cleared and the new values inserted
     *
     * \param ilist Initializer_list to assign to the container
     */
    void assign(std::initializer_list<value_type> ilist) {
        this->clear();
        this->reserve(ilist.size());
        this->insert(this->begin(), std::begin(ilist), std::end(ilist));
    }

    /** \brief Deletes all elements from the vector */
    void clear() {
        destroy_data(this->data());
        this->_size = 0;
    }

private:

    /** Check the capacity compared to new_size
     *
     * \param new_size Size to check capacity with
     *
     * \return true if capacity was increased, false otherwise.
     */
    bool check_capacity(size_type new_size) {
        if(new_size > this->capacity()) {
            size_type new_capacity =
                 new_size + this->capacity();
            this->reserve(new_capacity);
            return true;
        }
        return false;
    }


    /** \brief Shifts [pos, end) values to the right in a sub-array
     *
     * \warning This subfunction does not check for size / capacity. You must do
     * it before calling it.
     *
     * \tparam T   Sub-array type
     * \tparam Idx Sub-array index in the data() tuple
     *
     * \param pos   Position from which to shift the values
     * \param count Shift count
     */
    template<typename T, std::size_t Idx>
    T* array_move_to_right(const T* pos, size_type count) {
        typename Allocator::template rebind<T>::other alloc;
        // Start at the end of the array to avoid overlapping errors
        T* p = std::get<Idx>(this->data()) + this->size() - 1;
        // New end
        T* new_p = p + count;
        // Shift values
        while(p >= pos) {
            alloc.construct(new_p, *p);
            alloc.destroy(p); // Destroy p so that it can be overwritten
            --p;
            --new_p;
        }
        return p;
    }


    /** \brief Insert values in a sub-array by copy
     *
     * \warning This subfunction does not check for size / capacity. You must do
     * it before calling it.
     *
     * \tparam T Sub-array type
     * \tparam Idx Sub-array index in the data() tuple
     * \tparam U Value to insert type
     *
     * \param pos Position before which to insert the values
     * \param count Number of copies to insert
     * \param value Value to move or copy (universal reference)
     *
     * \warning If value is an rvalue, you shouldn't call this with count other
     * than 1.
     */
    template<typename T, std::size_t Idx, typename U>
    T* array_insert(const T* pos, size_type count, U&& value) {
        typename Allocator::template rebind<T>::other alloc;

        if(pos != std::get<Idx>(this->end())) {
            this->array_move_to_right<T, Idx>(pos, count);
        }

        T* p = const_cast<T*>(pos);
        T const * const end_pos = p + count;
        while(p != end_pos) {
            alloc.construct(p, std::forward<U>(value));
            ++p;
        }

        return p - count;
    }

    /** \brief Implements the insert algorithm for value movement or copy
     *
     * \tparam Values a pack of values or references
     *
     * \param pos Position before which to insert the values
     * \param count Number of copies to insert
     * \param value Value to copy
     *
     * \warning If value is an rvalue, you shoulnd't call this with count other
     * than 1.
     */
    template <typename... Values>
    iterator insert_impl(
        iterator pos,
        size_type count,
        Values&&... values)
    {
        difference_type pos_idx = pos - this->begin();

        if(check_capacity(this->size() + count)) {
            pos = this->begin() + pos_idx;
        }

        auto it = std::make_tuple (
            this->array_insert<Types, Indices>(
                std::get<Indices>(pos),
                count,
                std::forward<Values>(values)
                ) ... );

        this->_size += count;
        return iterator(it);
    }

public:
    /** \brief Insert an element by copy or move
     *
     * \tparam ValueTuple Tuple of types convertible to the Types parameter pack
     *
     * \param pos Position before which to insert the value
     * \param values Tuple of values to insert in each sub-array
     *
     * \return An iterator pointing to the inserted element
     */
    template <typename ValueTuple>
    iterator insert(iterator pos, ValueTuple&& values) {
        static_assert(std::tuple_size<std::decay_t<ValueTuple>>::value == sizeof...(Indices),
                      "Given tuple does not have the right size");
        return insert_impl(pos, 1, std::get<Indices>(std::forward<ValueTuple>(values))...);
    }

    /** \brief Insert copies of an element
     *
     * \tparam ValueTuple Tuple of types convertible to the Types parameter pack
     *
     * \param pos Position before which to insert the value
     * \param count Number of copies to insert&
     * \param values Tuple of values to be copied into the sub-arrays
     *
     * \return An iterator pointing to the first inserted element
     */
    template <typename ValueTuple>
    iterator insert(iterator pos, size_type count, ValueTuple&& values) {
        static_assert(std::tuple_size<std::decay_t<ValueTuple>>::value == sizeof...(Indices),
                      "Given tuple does not have the right size");
        return insert_impl(pos, count, std::get<Indices>(std::forward<ValueTuple>(values))...);
    }

private:
    /** \brief Insert values in a sub-array from a range
     *
     * \warning This subfunction does not check for size / capacity. You must do
     * it before calling it.
     *
     * \param pos Position before which to insert the values
     * \param count Number of copies to insert
     * \param value Value to copy
     *
     * \warning If value is an rvalue, you shoulnd't call this with count other
     * than 1.
     */
    template<typename T, std::size_t Idx, typename Iterator>
    T* array_insert_range(const T* pos, Iterator first, Iterator last) {
        typename Allocator::template rebind<T>::other alloc;
        size_type count = last - first;
        // We need to be sure we have the right pointer if there is some
        // reallocation
        this->array_move_to_right<T, Idx>(pos, count);
        T* p = const_cast<T*>(pos);
        while(first != last) {
            alloc.construct(p, *first);
            ++p;
            ++first;
        }

        return p - count;
    }

    /** Trait that checks whether a type is a tuple of iterators or not */
    template<class T>
    struct is_tuple_get_like {
        template<class U, class = decltype(std::get<0>(std::declval<U>()))>
        static constexpr bool check(U*) {return true;}
        static constexpr bool check(...) {return false;}
        enum {value = check(static_cast<T*>(nullptr))};
    };

public:

    /** \brief Insert elements from a range
     *
     * \tparam TupleIterator Iterator type convertible to a tuple of iterators
     * that are convertible to the elements of #value_type.
     *
     * \param pos Position before which to insert values
     * \param first Begining of range to copy
     * \param last End of range to copy
     *
     * \return An iterator pointing to the first element inserted
     *
     * \warning Behaviour is undefined if the range belongs to the container
     */
    template<class TupleIterator,
             typename std::enable_if<is_tuple_get_like<TupleIterator>::value, char>::type = 0
             >
    iterator insert(iterator pos, TupleIterator first, TupleIterator last) {
        size_type pos_idx = pos - this->begin();
        size_type count = last - first;


        if(check_capacity(this->size() + count)) {
            pos = this->begin() + pos_idx;
        }

        auto it = iterator {
            this->array_insert_range<Types, Indices>( // Magic happens here
                std::get<Indices>(pos),
                std::get<Indices>(first),
                std::get<Indices>(last)
                ) ... };

        this->_size += count;

        return it;
    }


    /** \brief Insert elements from a range
     *
     * \tparam ValueIterator Iterator type convertible to a tuple of elements
     * that are convertible to the elements of #value_type.
     *
     * \param pos Position before which to insert values
     * \param first Begining of range to copy
     * \param last End of range to copy
     *
     * \return An iterator pointing to the first element inserted
     *
     * \warning Behaviour is undefined if the range belongs to the container
     */
    template<class ValueIterator,
             typename std::enable_if<! is_tuple_get_like<ValueIterator>::value, char>::type = 0
             >
    iterator insert(iterator pos, ValueIterator first, ValueIterator last) {
        auto ret_offset = pos - this->cbegin();
        while(first != last) {
            insert(pos, *first);
            ++pos;
            ++first;
        }
        return this->begin() + ret_offset;
    }

    /** Erase an element */
    iterator erase(const_iterator pos) {
        return erase(pos, pos+1);
    }

private:

    /** \brief Moves [pos, end) count steps to the left in a sub-array, erasing
     * existing items
     *
     * \tparam T Sub-array type
     * \tparam Idx Sub-array index in the data() tuple
     *
     * \param pos Position from which to shift to the left
     * \param count Shift steps
     *
     * \return A pointer to the new position cooresponding to `pos`
     *
     * \warning This methods is to be applied to sub-arrays as part of the #erase
     * algorithm, therefore it does not change the vector properties.
     */
    template<typename T, std::size_t Idx>
    T* array_move_to_left(const T* pos, difference_type count) {
        typename Allocator::template rebind<T>::other alloc;
        T* end_ptr = std::get<Idx>(this->data()) + this->size();
        T* p = const_cast<T*>(pos);
        T* new_p = p - count;
        // Copy elements from left to right to avoid overlapping
        while(p < end_ptr) {
            alloc.destroy(new_p);
            alloc.construct(new_p, *p);

            ++p;
            ++new_p;
        }
        // destroy remaining items
        while(new_p < end_ptr) {
            alloc.destroy(new_p);
            ++new_p;
        }

        return const_cast<T*>(pos) - count;
    }

public:

    /** \brief Erase an interval
     *
     * \param first Iterator to the to be erased interval beginning
     * \param last Iterator past the to be erased interval end
     *
     * \return An iterator to the element past the interval, end() if there is
     * none
     *
     * \warning Behaviour is undefined if the interval is not part of the
     * container.
     */
    iterator erase(const_iterator first, const_iterator last) {
        const difference_type distance = last - first;
        auto it = iterator(
            this->array_move_to_left<Types, Indices>(
                std::get<Indices>(last),
                distance
                )...
            );
        this->_size -= static_cast<size_type>(distance);
        return it;
    }


    /** \brief Add an element to the end by copy
     *
     * \param values Tuple of const references to values to be copied at the end
     * of the vector
     */
     void push_back(const const_reference_tuple& values) {
         this->insert(this->end(), values);
     }


    /** \brief Add an elements to the end by copy
     *
     * \tparam Args Types of the values to insert into the sub-arrays
     *
     * \param values Tuple of values to insert
     */
    template<typename... Args>
    void push_back(const std::tuple<Args...>& values) {
        static_assert(sizeof...(Args) == sizeof...(Indices),
                      "Given tuple does not have the right size");

        this->insert(this->end(), values);
    }


    /** \brief Add an element to the end by copy
     *
     * \param values Values to insert into the sub-arrays
     */
    void push_back(Types const&... values) {
        this->insert_impl(this->end(), 1, values...);
    }


    /** \brief Add an element to the end by move
     *
     * \param values Values to move into the sub-arrays
     */
    void push_back(Types&&... values) {
        this->insert_impl(this->end(), 1, values...);
    }


    /** \brief Remove the last element */
    void pop_back() {
        this->erase(this->end()-1);
    }


    /** \brief Change the number of elements stored
     *
     * Resizes the container to given size. If the container holds more than
     * `count` elements, the excess ones are deleted, if it holds less than
     * `count` values, the new ones are default constructed.
     *
     * \param count New container size
     */
    void resize(size_type count) {
        this->resize(count, value_type{} );
    }


    /** \brief Change the number of elements stored
     *
     * Resizes the container to given size. If the container holds more than
     * `count` elements, the excess ones are deleted; if it holds less than
     * `count` elements, the new ones are copy constructed from `values`?
     *
     * \param count New container size
     * \param values Elements to copy if `count > size()`
     */
    void resize(size_type count, const value_type& values) {
        if(count < this->size()) {
            this->erase(this->cbegin() + count, this->cend());
        } else if(count > this->size()) {
            this->insert(this->end(), count - this->size(), values);
        }
    }


    /** \brief Change the number of elements stored
     *
     * Resizes the container to given size. If the container holds more than
     * `count` elements, the excess ones are deleted; if it holds less than
     * `count` elements, the new ones are copy constructed from `values`?
     *
     * \param count New container size
     * \param values Elements to copy if `count > size()`
     */
    void resize(size_type count, const Types&... values) {
        resize_impl(count, values..., indices);
    }

    /** \brief Swaps the contents of two vectors
     *
     * \param other Vector to swap contents with
     */
    void swap(variadic_vector_impl& other) noexcept {
        using std::swap; // uses std::swap if no swap function have been created
        swap(_data_allocator, other._data_allocator);
        swap(_size, other._size);
        swap(_capacity, other._capacity);
        swap(_allocator_tuple, other._allocator_tuple);
        swap(_data_tuple, other._data_tuple);
    }


    /** \brief Checks whether the contents of two vectors are equal
     *
     * \param lhs Left hand vector
     * \param rhs Right hand vector
     */
    friend bool operator==(const variadic_vector_impl& lhs, const variadic_vector_impl& rhs) {
        return (lhs.size() == rhs.size())
            && std::equal(lhs.begin(), lhs.end(), rhs.begin());
    }

    /** \brief Checks whether the contents of two vectors are different
     *
     * \param lhs Left hand vector
     * \param rhs Right hand vector
     */
    friend bool operator!=(const variadic_vector_impl& lhs, const variadic_vector_impl& rhs) {
        return !(lhs == rhs);
    }



    static_assert(sizeof...(Types) >= 1, "The vector must be instanciated with one or more types.");
};






template<typename... Types, std::size_t... Indices>
class variadic_vector_iterator<std::tuple<Types...>, inria::index_sequence<Indices...> >
    : public std::tuple<Types*...> {

    using base_t            = std::tuple<Types*...>;
    using value_tuple       = std::tuple<typename std::remove_const<Types>::type...>;
    using const_value_tuple = std::tuple<const Types...>;

public:
    using difference_type   = std::ptrdiff_t;
    using value_type        = std::tuple<Types...>;
    using reference         = std::tuple<Types&...>;
    using pointer           = void;
    using iterator_category = std::random_access_iterator_tag;
   using base_seq = inria::index_sequence<Indices...> ;
    using pointer_tuple     = std::tuple<Types*...>;


#ifdef __INTEL_COMPILER
    using std::tuple<Types*...>::tuple;
#else
  using typename base_t::tuple;
#endif
    variadic_vector_iterator(base_t tup_in) : base_t(tup_in){}
    variadic_vector_iterator(Types*... val) : base_t(val...){}
  variadic_vector_iterator(base_t tup_in, base_seq t) : base_t(tup_in){}
    /*
     */
    //    variadic_vector_iterator(base_t tup_in):(tup_in){}
    /** \brief Assignment addition operator
     *
     * \param n Shift steps to apply to the iterator
     *
     * \return The moved iterator
     */
    variadic_vector_iterator& operator+=(difference_type n) {
        noop_t{(std::get<Indices>(*this) += n)...};
        return (*this);
    }

    /** \brief Prefix increment operator
     *
     * Moves the iterator one step to the right.
     *
     * \return The incremented iterator
     */
    variadic_vector_iterator& operator++() {
        return (*this) += 1;
    }

    /** \brief Postfix increment operator
     *
     * Moves the operator to the right.
     *
     * \param unnamed Used to differenciate the overloads
     *
     * \return A copy of the iterator before incrementation
     */
    variadic_vector_iterator operator++(int) {
        variadic_vector_iterator tmp(*this);
        *this += 1;
        return tmp;
    }

    /** \brief Addition operator
     *
     * \param lhs Iterator
     * \param n   Shift from iterator
     *
     * \return A new iterator shifted by `n` of `lhs`
     */
    friend variadic_vector_iterator operator+(variadic_vector_iterator lhs, difference_type n) {
        return lhs += n;
    }

    /** \brief Addition operator
     *
     * \param n   Shift from iterator
     * \param lhs Iterator
     *
     * \return A new iterator shifted by `n` of `lhs`
     */
    friend variadic_vector_iterator operator+(difference_type n, variadic_vector_iterator lhs) {
        return lhs += n;
    }

    /** \brief Assignment substraction operator
     *
     * \param n Shift steps to apply to the iterator
     *
     * \return The moved iterator
     */
    variadic_vector_iterator& operator-=(difference_type n) {
        return (*this) += -n;
    }

    /** \brief Prefix decrement operator
     *
     * Moves the iterator one step to the left.
     *
     * \return The decremented iterator
     */
    variadic_vector_iterator& operator--() {
        return (*this) += -1;
    }

    /** \brief Postfix decrement operator
     *
     * Moves the operator to the left.
     *
     * \param unnamed Used to differenciate the overloads
     *
     * \return A copy of the iterator before decrementation
     */
    variadic_vector_iterator operator--(int) {
        variadic_vector_iterator tmp(*this);
        *this += -1;
        return tmp;
    }

    /** \brief Subtraction operator
     *
     * \param lhs Iterator
     * \param n   Opposite shift from iterator
     *
     * \return A new iterator shifted by `-n` of `lhs`
     */
    friend variadic_vector_iterator operator-(variadic_vector_iterator lhs, difference_type n) {
        return lhs += -n;
    }


    /** \brief Distance operator
     *
     * Computes the distance between this iterator and the `other` one.
     *
     * \param other The other iterator
     *
     * \return The distance between the two operators
     */
    difference_type operator-(const variadic_vector_iterator<value_tuple, inria::index_sequence<Indices...> >& other) const {
        return std::get<0>(*this) - std::get<0>(other);
    }

    /** \brief Distance operator
     *
     * Computes the distance between this iterator and the `other` one.
     *
     * \param other The other iterator
     *
     * \return The distance between the two operators
     */
    difference_type operator-(const variadic_vector_iterator<const_value_tuple, inria::index_sequence<Indices...> >& other) const {
        return std::get<0>(*this) - std::get<0>(other);
    }

    /** \brief Dereference operator
     *
     * \return A tuple of references to the values of the pointer element
     */
    reference operator*() {
        return reference(*std::get<Indices>(*this)...);
    }

    /** \brief Shifted dereference operator
     *
     * \param n Shift from this iterator of the iterator to dereference
     *
     * \return A tuple of references to the values of the pointer element
     */
    reference operator[](difference_type n) {
        return *(*this + n);
    }
};


template<typename Allocator, typename... Types>
class variadic_vector_type {
    template<std::size_t... Is>
    constexpr static auto get_type(std::tuple<Types...>, inria::index_sequence<Is...>)
        -> variadic_vector_impl<Allocator, std::tuple<Types...>, inria::index_sequence<Is...>> {
        return variadic_vector_impl<Allocator, std::tuple<Types...>, inria::index_sequence<Is...>>();
    }
public:
    using type = decltype(get_type(std::declval<std::tuple<Types...>>(), inria::make_index_sequence<sizeof...(Types)>()));
};

/** \brief Variadic vector type
 *
 * \tparam Allocator Allocator type to use
 * \tparam Types Parameter pack for type to store
 */
template<typename Allocator, typename... Types>
struct variadic_vector :
#ifndef DOXYGEN_DOC
    public variadic_vector_type<Allocator, Types...>::type
{
    using base_t = typename variadic_vector_type<Allocator, Types...>::type;
    using base_t::base_t;
};
#else
    public variadic_vector_impl<Allocator, std::tuple<Types...>, inria::index_sequence<Indices...> >
{};
#endif



#endif
