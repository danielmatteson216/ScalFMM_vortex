// See LICENCE file at project root
#
#include <cstddef>

#include "FUTester.hpp"

#include "Containers/FBufferReader.hpp"
#include "Containers/FBufferWriter.hpp"



/** this class test the buffers container */
class TestBuffer : public FUTester<TestBuffer> {

        // test size
        void TestWriteRead(){
            FBufferWriter writer;

            const int BytesTested = static_cast<int>(sizeof(int)+sizeof(char)+sizeof(double)+sizeof(float));
            const int NbTest = 5;
            for(int idxWrite = 0 ; idxWrite < NbTest ; ++idxWrite){
                writer << idxWrite << char(idxWrite) << double(idxWrite) << float(idxWrite);
            }

            uassert(writer.getCapacity() >= NbTest*BytesTested);
            uassert(writer.getSize() == NbTest*BytesTested);

            FBufferReader reader(writer.getSize());
            uassert(reader.getSize() == 0);
            uassert(reader.getCapacity() == writer.getSize());

            memcpy(reader.data(), writer.data(), writer.getSize());
            for(int idxRead = 0 ; idxRead < NbTest ; ++idxRead){
                int intval;
                char charval;
                double doubleval;
                float floatval;
                reader >> intval >> charval >> doubleval >> floatval;

                uassert(intval == idxRead);
                uassert(charval == char(idxRead));
                uassert(doubleval == double(idxRead));
                uassert(floatval == float(idxRead));

                uassert(reader.tell() == (BytesTested * (idxRead+1)));
            }

            uassert(reader.tell() == reader.getCapacity());
            reader.seek(0);
            uassert(reader.tell() == 0);

            for(int idxRead = 0 ; idxRead < NbTest ; ++idxRead){
                uassert(reader.getValue<int>() == idxRead);
                uassert(reader.getValue<char>() == char(idxRead));
                uassert(reader.getValue<double>() == double(idxRead));
                uassert(reader.getValue<float>() == float(idxRead));

                uassert(reader.tell() == (BytesTested * (idxRead+1)));
            }

            uassert(reader.tell() == reader.getCapacity());
        }


        void TestWriteAt(){
            FBufferWriter writer;

            const int SizeOfInt = int(sizeof(int));
            const int NbTest = 5;
            for(int idxWrite = 0 ; idxWrite < NbTest ; ++idxWrite){
                const FSize position = writer.getSize();

                uassert(position == (NbTest * SizeOfInt * idxWrite) + (idxWrite * SizeOfInt));

                writer.FBufferWriter::write<int>(0);

                uassert(writer.getSize() == (NbTest * SizeOfInt * idxWrite) + (idxWrite * SizeOfInt) + SizeOfInt);

                for(int count = 1 ; count <= NbTest ; ++count) writer << count;
                writer.writeAt(position, idxWrite);
            }

            FBufferReader reader(writer.getSize());
            memcpy(reader.data(), writer.data(), writer.getSize());

            for(int idxWrite = 0 ; idxWrite < NbTest ; ++idxWrite){
                uassert(reader.FBufferReader::getValue<int>() == idxWrite);
                for(int count = 1 ; count <= NbTest ; ++count){
                    uassert(reader.FBufferReader::getValue<int>() == count);
                }
            }
        }


        // set test
        void SetTests(){
            AddTest(&TestBuffer::TestWriteRead,"Test Write then Read");
            AddTest(&TestBuffer::TestWriteAt,"Test Write at then Read");
        }
};

// You must do this
TestClass(TestBuffer)
