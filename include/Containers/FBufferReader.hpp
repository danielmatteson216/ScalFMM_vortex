// See LICENCE file at project root
#ifndef FBUFFERREADER_HPP
#define FBUFFERREADER_HPP

#include <memory>
#include <algorithm>
#include "FAbstractBuffer.hpp"
#include "FBufferWriter.hpp"
#include "Utils/FAssert.hpp"

/**
 * \brief Provides memory management and conversion to basic type
 * \author Cyrille Piacibello, Berenger Bramas, Quentin Khan
 *
 * This class is meant to ease data deserialisation through external libraries
 * such as MPI.
 *
 *   - Reserve memory.
 *   - Pass the data() pointer the the library that will fill the buffer.
 *   - Read using the getValue, fillValue, fillArray methods.
 *
 * An internal index is kept to know what data has been read.
 */
class FBufferReader : public FAbstractBufferReader {
    FSize arrayCapacity;            ///< Allocated space
    std::unique_ptr<char[]> array;  ///< Allocated array
    FSize currentIndex;             ///< First unread byte

public :

    /**
     * \brief Construct the reader
     *
     * \param capacity Buffer capacity in bytes
     */
    explicit FBufferReader(const FSize capacity = 512)
        : arrayCapacity(capacity),
          array(new char[capacity]),
          currentIndex(0)
    {
        FAssertLF(array, "Cannot allocate array");
    }
    /**
     * \brief Construct the reader
     *
     * \param capacity Buffer capacity in bytes
     */
    explicit FBufferReader( FBufferWriter& buf)
        : arrayCapacity(buf.getCapacity() )
    {
      this->cleanAndResize(arrayCapacity) ;
      std::unique_ptr<char[]> arraytmp(new char[arrayCapacity]);
      std::copy(buf.data(),buf.data()+arrayCapacity, arraytmp.get());
      array = std::move(arraytmp);
      currentIndex = buf.getSize() ;

        FAssertLF(array, "Cannot allocate array");
    }


    virtual ~FBufferReader() = default;

    /**
     * \brief Discard current buffer and allocate a new one if needed
     *
     * The read data index is reset.
     *
     * \param capacity New buffer capacity
     */
    void cleanAndResize(const FSize capacity) {
        if(capacity != arrayCapacity){
            arrayCapacity = capacity;
            array.reset(new char[capacity]);
        }
        currentIndex = 0;
    }

    /**
     * \brief Reserve memory if needed, existing data is copied
     *
     * If the current capacity is enough, this operation is a no-op. Otherwise a
     * new buffer is allocated and data is copied from the old buffer before
     * discarding it.
     *
     * \param newCapacity Required minimal capacity
     */
    void reserve(const FSize newCapacity) {
        if(newCapacity > arrayCapacity) {
            std::unique_ptr<char[]> new_array(new char[newCapacity]);
            // The array is modified outside the class, we don't know what must
            // be copied or not : we copy everything
            std::copy(array.get(), array.get()+arrayCapacity, new_array.get());
            array = std::move(new_array);
            arrayCapacity = newCapacity;
        }
    }


    /**
     * \brief Get pointer to allocated memory
     */
    char* data() override {
        return array.get();
    }

    /**
     * \brief Get pointer to allocated memory
     */
    const char* data() const override  {
        return array.get();
    }

    /**
     * \brief Return already deserialised data size in bytes
     *
     * \return The size of the subpart of the buffer that has already been
     * deserialised.
     */
    FSize getSize() const override {
        return currentIndex;
    }

    /**
     * \brief Return the allocated buffer capacity in bytes
     */
    FSize getCapacity() const{
        return arrayCapacity;
    }

    /**
     * \brief Move the read index to a position
     *
     * \param index Position from the beginning of the buffer to move the read
     *              index to.
     */
    void seek(const FSize index) override {
        FAssertLF(index <= arrayCapacity,
                  "FBufferReader :: can't move index because buffer isn't ",
                  "long enough ", index, " ", arrayCapacity);
        currentIndex = index;
    }

    /**
     * \brief Get the read index value
     */
    FSize tell() const override  {
        return currentIndex;
    }

    /**
     * \brief Deserialise an object at read index
     *
     * The object is byte copied into a new instance of T which is returned by
     * value.
     *
     * The read index is incremented by sizeof(T).
     *
     * \tparam T Type of the object to deserialise
     */
    template <class T>
    T getValue(){
        FAssertLF(currentIndex + FSize(sizeof(T)) <= arrayCapacity,
                  "The buffer does not have enough remaining memory to read a ",
                  " value of given type");
        T value;
        fillValue<T>(&value);
        return value;
    }

    /**
     * \brief Deserialise an object at a specified index
     *
     * The object is byte copied into a new instance of T which is returned by
     * value.
     *
     * The read index is changed to `index + sizeof(T)`.
     *
     * \tparam T Type of the object to deserialise
     *
     * \param index Position of the object
     */
    template <class T>
    T getValue(const FSize index){
        seek(index);
        return getValue<T>();
    }

    /**
     * \brief Deserialise an object to given address
     *
     * Equivalent to `fillArray(ptr, 1)`.
     *
     * \tparam T Type of the object to deserialise
     *
     * \param ptr Pointer to the object to copy to
     */
    template <class T>
    void fillValue(T* const ptr){
        fillArray(ptr, 1);
    }

    /**
     * \brief Deserialise contiguous values
     *
     * The objects are byte copied to the array pointed to by inArray.
     *
     * The read index is incremented by `count * sizeof(T)`
     *
     * \tparam T Type of the object to deserialise
     *
     * \param inArray Array of objects to copy the values to
     * \param count Object count in the array
     */
    template <class T>
    void fillArray(T* const inArray, const FSize count){
        FAssertLF(currentIndex + FSize(sizeof(T))*count <= arrayCapacity );
        std::memcpy(reinterpret_cast<char*>(inArray), &array[currentIndex], sizeof(T)*count);
        currentIndex += sizeof(T)*count;
    }

    /**
     * \brief Stream-like deserialisation
     *
     * See fillArray.
     *
     * \tparam T Type of the object to deserialise
     *
     * \param object Object to deserialise to
     */
    template <class T>
    FBufferReader& operator>>(T& object){
        fillValue(&object);
        return *this;
    }

};
#endif
