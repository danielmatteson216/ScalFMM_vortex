// See LICENCE file at project root

#ifndef FUSERLEAFCONTAINER_HPP
#define FUSERLEAFCONTAINER_HPP

/**
 * @file This file contains a class that is an implementation of
 * FParticleContainer, designed to be used by the User Kernel API.
 */

#include <Components/FBasicParticleContainer.hpp>


/**
 * @author Piacibello
 *
 * @brief This class define another Particle Container, with dynamic
 * (i.e. no template static) storage in order to store the user
 * particles informations
 */
template<class FReal>
class FUserLeafContainer : public FP2PParticleContainerIndexed<FReal>{

    void * userAttributes;

    static Scalfmm_Leaf_Descriptor user_leaf_descriptor;
    using Parent = FP2PParticleContainerIndexed<FReal>;

public:

    static void Init(Scalfmm_Leaf_Descriptor leaf_descriptor){
        user_leaf_descriptor = leaf_descriptor;
    }

    FUserLeafContainer(const FUserLeafContainer&) = delete;
    FUserLeafContainer& operator =(const FUserLeafContainer&) = delete;

    FUserLeafContainer() : userAttributes(nullptr){

    }

    void setContainer(void * inputPtr){
        userAttributes = inputPtr;
    }
    void * getContainer() const {
        return userAttributes;
    }

    static Callback_init_leaf GetInitLeaf(){
        return user_leaf_descriptor.user_init_leaf;
    }

    static Callback_free_leaf GetFreeLeaf(){
        return user_leaf_descriptor.user_free_leaf;
    }

    /**
     * @brief Functions to deal with serialization.
     *
     */
    template <class BufferWriterClass>
    void save(BufferWriterClass& buffer) const{
        Parent::template save<>(buffer);
        if(user_leaf_descriptor.user_copy_leaf){
            FSize nbPart = Parent::getNbParticles();
            FSize sizeToSave = user_leaf_descriptor.user_get_size(nbPart);
            char * temp = new char[sizeToSave];
            user_leaf_descriptor.user_copy_leaf(nbPart,getContainer(),temp);
            buffer.write(temp,sizeToSave);
            delete [] temp;
        }else{
            std::cout<<"No user_copy_leaf function set\nExiting\n";
            exit(0);
        }
    }
    template <class BufferReaderClass>
    void restore(BufferReaderClass& buffer){
        Parent::template restore<>(buffer);
        if(user_leaf_descriptor.user_restore_leaf){
            FSize nbPart = Parent::getNbParticles();
            FSize sizeToSave = user_leaf_descriptor.user_get_size(nbPart);

            char * temp = new char[sizeToSave];
            buffer.fillArray(temp,sizeToSave);
            this->setContainer(user_leaf_descriptor.user_restore_leaf(nbPart,temp));
            delete [] temp;
        }else{
            std::cout<<"No user_restore_leaf function set\nExiting\n";
            exit(0);
        }
    }

    FSize getSavedSize() const {
        //Size of internal datas
        FSize res = Parent::getSavedSize();
        FSize userSize = 0;
        if(!user_leaf_descriptor.user_get_size){
            std::cout<<"No get_save_size function set\nExiting\n";
            exit(0);
        }else{
            FSize nbPart = Parent::getNbParticles();
            userSize = user_leaf_descriptor.user_get_size(nbPart);
            res += userSize;
        }

        return res;
    }
};



#endif // FUSERLEAFCONTAINER_HPP
