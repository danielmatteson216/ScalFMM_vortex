#ifndef TCLI_EXCEPTIONS_HPP
#define TCLI_EXCEPTIONS_HPP

#include <stdexcept>

namespace inria {
namespace tcli {

struct parameter_conflict : std::logic_error {
    using std::logic_error::logic_error;
};

struct invalid_parameter : std::logic_error {
    using std::logic_error::logic_error;
};

struct unknown_parameter : std::runtime_error {
    using std::runtime_error::runtime_error;
};

struct missing_required_parameter : std::runtime_error {
    using std::runtime_error::runtime_error;
};

struct partial_parse : std::runtime_error {
    using std::runtime_error::runtime_error;
};

struct exclusive_parameters : std::runtime_error {
    using std::runtime_error::runtime_error;
};

struct parse_error : std::runtime_error {
    using std::runtime_error::runtime_error;

    bool missing_required = false;
    bool unknown_parameter = false;
};

}} // close namespace inria::tcli

#endif /* TCLI_EXCEPTIONS_HPP */
