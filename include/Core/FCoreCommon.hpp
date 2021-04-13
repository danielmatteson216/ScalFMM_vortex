// See LICENCE file at project root
#ifndef FCORECOMMON_HPP
#define FCORECOMMON_HPP

#include <vector>
#include "Utils/FGlobal.hpp"
#include "Utils/FAssert.hpp"

/**
 * @brief The FFmmOperations enum
 * To chose which operation has to be performed.
 */
///
/// \brief The FFmmOperations enum
///  To chose which operation has to be performed in the methode execute of the Fmm algorithm.
///
enum FFmmOperations {
    FFmmP2P   = (1 << 0),  ///< Particles to Particles operator (Near field)
    FFmmP2M  = (1 << 1),   ///< Particles to Multipole operator (Far field)
    FFmmM2M = (1 << 2),    ///< Multipole to Multipole operator (Far field)
    FFmmM2L  = (1 << 3),   ///< Multipole to Local operator     (Far field)
    FFmmL2L  = (1 << 4),   ///< Local to Local operator         (Far field)
    FFmmL2P  = (1 << 5),   ///< Local to Particles operator     (Far field)
    FFmmP2L  = (1 << 6),   ///< Particles to Local operator     (Far field in Adaptive algorithm)
    FFmmM2P  = (1 << 7),   ///< Multipole to Particles operator (Far field in Adaptive algorithm)
//
    FFmmNearField = FFmmP2P,  ///< Near field operator
    FFmmFarField  = (FFmmP2M|FFmmM2M|FFmmM2L|FFmmL2L|FFmmL2P|FFmmM2P|FFmmP2L),  ///< Only Far Field operators
//
    FFmmNearAndFarFields = (FFmmNearField|FFmmFarField)  ///< Near and far field operators
};
///
/// \brief FFmmOperations_string converts the FFmmOperations (enum) in string
/// \param value the FFmmOperations
/// \return the string corresponding to the FFmmOperations
///
inline std::string FFmmOperations_string(/*enum FFmmOperations*/ const unsigned int & value){

  std::string op("");
  if (value & FFmmP2P)
    op += " FFmmP2P |";
  if (value & FFmmP2M)
    op += " FFmmP2M |";
  if (value & FFmmM2M)
    op += " FFmmM2M |";
  if (value & FFmmM2L)
    op += " FFmmM2L |";
  if (value & FFmmL2L)
    op += " FFmmL2L |";
  if (value & FFmmL2P)
    op += " FFmmL2P |";
  op.erase(op.size()-2,op.size()-1);
  return op;
};

/**
 * \brief Algorithm interface
 *
 * All algorithms should implement this interface.
 */
struct FAlgorithmInterface {

    /** \brief Destructor */
    virtual ~FAlgorithmInterface() {}

    /** \brief Get algorithm short description
     *
     * \return A short string identifying the algorithm implementation
     */
    virtual std::string name() const {
        return "error: unnamed algorithm";
    }

    /** \brief Get algorithm environment description
     *
     * Easily parsable string describing the algorithm settings and environment
     * (thread count, specific settings, etc).
     *
     * \return String describing algorithm settings
     */
    virtual std::string description() const {
        return "error: description not implemented";
    }

    /** \brief Run specific steps of the algorithm
     *
     * \param operations Specifies the algorithm operations to run, see
     * FFmmOperations.
     */
    virtual void execute(const unsigned int operations) = 0;

    /** \brief Run the algorithm */
    virtual void execute() {
        this->execute(FFmmNearAndFarFields);
    }
};





/**
 * \brief Base class of algorithms
 *
 * This class is an abstract algorithm to be able to use the FAlgorithmBuilder
 * and execute from an abstract pointer.
 */
class FAbstractAlgorithm : public FAlgorithmInterface {
protected:

    int upperWorkingLevel; ///< Where to start the work
    int lowerWorkingLevel; ///< Where to end the work (exclusive)
    int nbLevelsInTree;    ///< Height of the tree

    void setNbLevelsInTree(const int inNbLevelsInTree){
        nbLevelsInTree    = inNbLevelsInTree;
        lowerWorkingLevel = nbLevelsInTree;
    }

    void validateLevels() const {
        FAssertLF(FAbstractAlgorithm::upperWorkingLevel <= FAbstractAlgorithm::lowerWorkingLevel);
        FAssertLF(2 <= FAbstractAlgorithm::upperWorkingLevel);
    }
    virtual void executeCore(const unsigned operationsToProceed) = 0;

public:
    FAbstractAlgorithm()
        : upperWorkingLevel(2), lowerWorkingLevel(0), nbLevelsInTree(-1){
    }

    virtual ~FAbstractAlgorithm(){
    }

    /** \brief Execute the whole fmm for given levels. */
    virtual void execute(const int inUpperWorkingLevel, const int inLowerWorkingLevel) final {
        upperWorkingLevel = inUpperWorkingLevel;
        lowerWorkingLevel = inLowerWorkingLevel;
        validateLevels();
        executeCore(FFmmNearAndFarFields);
    }

    /** \brief Execute only some FMM operations for given levels. */
    virtual void execute(const unsigned operationsToProceed, const int inUpperWorkingLevel, const int inLowerWorkingLevel) final {
        upperWorkingLevel = inUpperWorkingLevel;
        lowerWorkingLevel = inLowerWorkingLevel;
        validateLevels();
        executeCore(operationsToProceed);
    }

    using FAlgorithmInterface::execute;

    /** \brief Execute only some steps. */
    virtual void execute(const unsigned operationsToProceed) override final {
        upperWorkingLevel = 2;
        lowerWorkingLevel = nbLevelsInTree;
        validateLevels();
        executeCore(operationsToProceed);
    }

    /// Build and dill vector of the  MortonIndex Distribution at Leaf level
    ///  p = mpi process id then
    ///  [mortonLeafDistribution[2*p], mortonLeafDistribution[2*p+1]  is the morton index shared by process p
    virtual void getMortonLeafDistribution(std::vector<MortonIndex> & mortonLeafDistribution) {
      mortonLeafDistribution.resize(2) ;
      mortonLeafDistribution[0] = 0 ;
      mortonLeafDistribution[1] = 8 << (nbLevelsInTree-1) ;

    };
};




#endif // FCORECOMMON_HPP
