// See LICENCE file at project root
#ifndef FBUFFERWRITER_HPP
#define FBUFFERWRITER_HPP

#include <memory>
#include "FAbstractBuffer.hpp"
#include "../Utils/FAssert.hpp"

/**
 * \brief Provides memory management and byte serialisation
 * \author Cyrille Piacibello, Berenger Bramas, Quentin Khan
 *
 * This class is meant to ease data byte serialisation to pass through external
 * libraries such as MPI.
 *
 *   - Reserve memory.
 *   - Insert objects.
 *   - Pass data() pointer to the library.
 *
 * An internal index is kept to know how much data has been written.
 */
class FBufferWriter : public FAbstractBufferWriter {
    FSize arrayCapacity;              ///< Allocated space
    std::unique_ptr<char[]> array;    ///< Allocated array
    FSize currentIndex;               ///< Currently filled space

    /**
     * \brief Ensure minimum remaining space in the buffer
     *
     * This methods checks whether there is at least min_rem_cap space left in
     * the buffer. If not, a new buffer is allocated and the old data is copied
     * into the new one.
     *
     * \param min_rem_cap Minimum remaining space in the buffer
     */
    void expandIfNeeded(const FSize min_rem_cap) {
        if( arrayCapacity < currentIndex + min_rem_cap){
            arrayCapacity = FSize(double(currentIndex + min_rem_cap + 1) * 1.5);
            char* arrayTmp = new char[arrayCapacity];
            std::copy(array.get(), array.get()+currentIndex, arrayTmp);
            array.reset(arrayTmp);
        }
    }

public:

    /**
     * \brief Construct the writer
     *
     * \param capacity Buffer capacity in bytes
     */
    explicit FBufferWriter(const FSize capacity = 1024)
        : arrayCapacity(capacity),
          array(new char[capacity]),
          currentIndex(0)
    {
        FAssertLF(array, "Cannot allocate array");
    }

    virtual ~FBufferWriter() = default;

    /**
     * \brief Change buffer capacity and copy existing data to the new buffer
     *
     * If write index is bigger than new_capacity, the extra data is lost and
     * the write index is changed accordingly.
     */
    void resize(const FSize new_capacity){
        if(new_capacity != arrayCapacity){
            arrayCapacity = new_capacity;
            char* arrayTmp = new char[arrayCapacity];
            currentIndex = (currentIndex < arrayCapacity ? currentIndex : arrayCapacity-1);
            std::copy(array.get(), array.get()+currentIndex, arrayTmp);
            array.reset(arrayTmp);
        }
    }

    /**
     * \brief Reset the write index
     *
     * Subsequent writes will overwrite the buffer.
     */
    void reset() override {
        currentIndex = 0;
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

    std::unique_ptr<char[]>* ptrData(){
      return &array;
    }
    /**
     * \brief Return the written memory size
     */
    FSize getSize() const override  {
        return currentIndex;
    }

    /**
     * \brief Return the allocated buffer capacity in bytes
     */
    FSize getCapacity() const {
        return arrayCapacity;
    }

    /**
     * \brief Byte copy an object in the buffer
     *
     * The write index is increased by `sizeof(ClassType)`
     *
     * \tparam ClassType Type of the object to write
     *
     * \param object Object to write to the buffer
     */
    template <class ClassType>
    void write(const ClassType& object) {
        write(&object, 1);
    }

    /**
     * \brief Byte copy an object at given position in the buffer
     *
     * \tparam ClassType Type of the object to write
     *
     * \param object Object to write to the buffer
     */
    template <class ClassType>
    void writeAt(const FSize position, const ClassType& object){
        FAssertLF(position+FSize(sizeof(ClassType)) <= currentIndex);
        memcpy(&array[position], &object, sizeof(ClassType));
    }

    /**
     * \brief Byte copy contiguous objects to the buffer
     *
     * \tparam ClassType Type of the object to write
     *
     * \param objects Array of objects to write
     * \param count Object count in the array
     */
    template <class ClassType>
    void write(const ClassType* const objects, const FSize inSize){
        expandIfNeeded(sizeof(ClassType) * inSize);
        memcpy(&array[currentIndex], objects, sizeof(ClassType)*inSize);
        currentIndex += sizeof(ClassType)*inSize;
    }

    /**
     * \brief Stream-like serialisation
     *
     * \tparam ClassType Type of the object to write
     */
    template <class ClassType>
    FBufferWriter& operator<<(const ClassType& object){
        write(object);
        return *this;
    }
};


#endif // FBUFFERWRITER_HPP
