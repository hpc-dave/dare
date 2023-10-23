/*
 * MIT License
 *
 * Copyright (c) 2022 David Rieder

 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef MPI_MPITYPECONVERTER_H_
#define MPI_MPITYPECONVERTER_H_

#include <mpi.h>

namespace dare::mpi {

/*!
 * \brief declaration of the generic function
 * Only specialization of this function are usable.
 * If you end up here, that means that the type that you
 * tried to provide to the MPI-routine wrappers have no
 * knowledge of how to convert the specified type to an
 * MPI_Datatype. That means you either have to use MPI_routines
 * directly by yourself or define your custom MPI-datatype and
 * add a template specialization
 */
template <typename T>
inline MPI_Datatype GetMPIType();

/* !
 * \brief Converts basic data type to MPI-datatype
 */
template <>
inline MPI_Datatype GetMPIType<long double>() {
    return MPI_LONG_DOUBLE;
}

/* \brief . */
template <>
inline MPI_Datatype GetMPIType<double>() {
    return MPI_DOUBLE;
}

/* \brief . */
template <>
inline MPI_Datatype GetMPIType<float>() {
    return MPI_FLOAT;
}

/* \brief . */
template <>
inline MPI_Datatype GetMPIType<bool>() {
    return MPI_CXX_BOOL;
}

/* \brief . */
template <>
inline MPI_Datatype GetMPIType<int8_t>() {
    return MPI_INT8_T;
}

/* \brief . */
template <>
inline MPI_Datatype GetMPIType<uint8_t>() {
    return MPI_UINT8_T;
}

/* \brief . */
template <>
inline MPI_Datatype GetMPIType<int16_t>() {
    return MPI_INT16_T;
}

/* \brief . */
template <>
inline MPI_Datatype GetMPIType<uint16_t>() {
    return MPI_UINT16_T;
}

/* \brief . */
template <>
inline MPI_Datatype GetMPIType<int32_t>() {
    return MPI_INT32_T;
}

/* \brief . */
template <>
inline MPI_Datatype GetMPIType<uint32_t>() {
    return MPI_UINT32_T;
}

/* \brief . */
template <>
inline MPI_Datatype GetMPIType<int64_t>() {
    return MPI_INT64_T;
}

/* \brief . */
template <>
inline MPI_Datatype GetMPIType<uint64_t>() {
    return MPI_UINT64_T;
}

}  // namespace dare::mpi

#endif  // MPI_MPITYPECONVERTER_H_
