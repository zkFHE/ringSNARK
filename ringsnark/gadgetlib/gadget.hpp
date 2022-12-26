/** @file
 *****************************************************************************
 * @author     This file is part of libsnark, and is developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef GADGET_HPP_
#define GADGET_HPP_

#include <ringsnark/gadgetlib/protoboard.hpp>

namespace ringsnark {

    template<typename RingT>
    class gadget {
    protected:
        protoboard<RingT> &pb;
        const std::string annotation_prefix;
    public:
        gadget(protoboard<RingT> &pb, const std::string &annotation_prefix = "");
    };

} // ringsnark
#include <ringsnark/gadgetlib/gadget.tcc>

#endif // GADGET_HPP_
