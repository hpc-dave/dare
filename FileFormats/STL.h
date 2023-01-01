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

#ifndef FILEFORMATS_STL_H_
#define FILEFORMATS_STL_H_

#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <cstring>
#include "../Utilities/Vector.h"

namespace dare::ff{

/*! \class STLfacet
 * \brief represents a single facet of an STL file
 * \tparam VT value type of the points, usually that should be floating point
 */
template<typename VT>
class STLfacet{
public:
    using ValueType = VT;   //!< value type definition
    using PointType = dare::utils::Vector<3, ValueType>;    //!< point type defintion (compatible with VectorType)
    using VectorType = dare::utils::Vector<3, ValueType>;   //!< vector type defintion (compatible with PointType)

    /*
     * \brief constructor with 3 points
     * @param v1 first vector
     * @param v2 second vector
     * @param v3 third vector
     *
     * This constructor will automatically compute the normal, depending
     * on the provided points!
     */
    STLfacet(const PointType& v1, const PointType& v2, const PointType& v3);

    /*!
     * \brief default constructor
     * all points and the normal are set to the default values of the Vector class
     */
    STLfacet();

    /*!
     * \brief copy constructor
     * @param f other facet
     */
    STLfacet(const STLfacet<VT>& f);

    /*!
     * \brief assignment operator
     * @param f other facet
     */
    STLfacet& operator=(const STLfacet<VT>& f);

    /*!
     * \brief applies translational motion to facet
     * @param trajectory trajectory along which the facet moves
    */
    STLfacet& Move(const VectorType& trajectory);

    /*!
     * \brief returns the normal vector of the facet
     */
    const VectorType& GetNormal();
    const VectorType& GetNormal() const;

    /*!
     * \brief returns array of points
     */
    const PointType* GetPoints();
    const PointType* GetPoints() const;

    /*!
     * \brief helper to output the facet
     */
    friend std::ostream& operator<<(std::ostream& os, const STLfacet<VT>& f) {
        os << "P:";
        for (int i{0}; i < 3; ++i)
            os << ' ' << f.points[i];

        os << " N: " << f.normal;
        return os;
    }

private:
    /*!
     * \brief computes the normal of the provided values
     */
    void ComputeNormal();

    PointType points[3];    //!< 3 points of the facet
    VectorType normal;      //!< normal of the facet
};

template<typename VT>
class STL{
public:
    using ValueType = VT;
    using Facet = STLfacet<ValueType>;
    using TriangleContainer = std::vector<Facet>;

    STL();
    STL(const STL& other);
    explicit STL(const TriangleContainer& facets);
    explicit STL(const std::string& file);

    const TriangleContainer& GetFacets();
    const TriangleContainer& GetFacets() const;

    bool ReadFile(const std::string& file);
    bool WriteToFile(const std::string& file, bool write_ASCII = false);

private:

    bool IsASCII(std::ifstream& in);

    bool ReadASCII(std::ifstream& in);
    bool ReadBinary(std::ifstream& in);

    bool WriteToASCII(std::ostream& ofs);
    bool WriteToBinary(std::ostream& ofs);

    TriangleContainer facets;
};

} // end namespace dare::ff

#include "STL.inl"

#endif // FILEFORMATS_STL_H_
