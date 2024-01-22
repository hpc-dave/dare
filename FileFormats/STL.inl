/*
 * MIT License
 *
 * Copyright (c) 2024 David Rieder

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

namespace dare::ff {

template <typename VT>
STLfacet<VT>::STLfacet(const STLfacet<VT>::PointType& v1,
                       const STLfacet<VT>::PointType& v2,
                       const STLfacet<VT>::PointType& v3)
    : points{v1, v2, v3} {
    ComputeNormal();
}

template <typename VT>
STLfacet<VT>::STLfacet() {
}

template <typename VT>
STLfacet<VT>::STLfacet(const STLfacet<VT>& f)
    : points{f.points[0], f.points[1], f.points[2]}, normal(f.normal) {
}

template <typename VT>
STLfacet<VT>& STLfacet<VT>::Move(const STLfacet<VT>::VectorType& trajectory) {
    points[0] += trajectory;
    points[1] += trajectory;
    points[2] += trajectory;
    return *this;
}

template <typename VT>
STLfacet<VT>& STLfacet<VT>::operator=(const STLfacet<VT>& f) {
    if (&f == this)
        return *this;

    points[0] = f.points[0];
    points[1] = f.points[1];
    points[2] = f.points[2];

    normal = f.normal;

    return *this;
}

template <typename VT>
const typename STLfacet<VT>::VectorType& STLfacet<VT>::GetNormal() {
    return normal;
}

template <typename VT>
const typename STLfacet<VT>::VectorType& STLfacet<VT>::GetNormal() const {
    return normal;
}

template <typename VT>
const typename STLfacet<VT>::PointType* STLfacet<VT>::GetPoints() {
    return points;
}

template <typename VT>
const typename STLfacet<VT>::PointType* STLfacet<VT>::GetPoints() const {
    return points;
}

template <typename VT>
void STLfacet<VT>::ComputeNormal() {
    const VectorType v1(points[1] - points[0]);
    const VectorType v2(points[2] - points[0]);
    normal = v1.cross(v2);
#ifndef NDEBUG
    if (normal.length() == static_cast<VT>(0))
        std::cerr << "determined facet with normal length of 0, expect division by 0\n";
#endif
    normal /= normal.length();
}

template <typename VT>
STL<VT>::STL() {
}

template <typename VT>
STL<VT>::STL(const STL<VT>& other) : facets(other.facets()) {
}

template <typename VT>
STL<VT>::STL(const typename STL<VT>::TriangleContainer& c) : facets(c) {
}

template <typename VT>
STL<VT>::STL(const std::string& file) {
    ReadFile(file);
}

template <typename VT>
bool STL<VT>::ReadFile(const std::string& file) {
    std::ifstream in(file, std::ios::in | std::ios::binary);
    if (!in) {
        std::cerr << "Could not open file '" << file << "'! Aborting read operation!" << std::endl;
        return false;
    }
    bool success{false};
    if (IsASCII(in))
        success = ReadASCII(in);
    else
        success = ReadBinary(in);

    if (!success) {
        std::cerr << "An error occured while reading the file '" << file << "!" << std::endl;
    }
    return success;
}

template <typename VT>
bool STL<VT>::WriteToFile(const std::string& file, bool write_ASCII) {
    std::ofstream ofs(file, std::ios::out | std::ios::binary);
    if (!ofs) {
        std::cerr << "Could not open file " << file << "! Aborting write operation!" << std::endl;
        return false;
    }
    bool success{false};
    if (write_ASCII)
        success = WriteToASCII(ofs);
    else
        success = WriteToBinary(ofs);

    if (!success) {
        std::cerr << "An error occured while writing to " << file << "!" << std::endl;
    }
    return success;
}

template <typename VT>
bool STL<VT>::IsASCII(std::ifstream& in) {
    char check_line[5];
    try {
        in.read(check_line, 5);
    } catch (const std::exception& ex) {
        std::cerr << "Follwing error occured during reading: " << ex.what() << ". Aborting!" << std::endl;
        return false;
    }
    in.clear();
    in.seekg(0, std::ios::beg);
    return strcmp(check_line, "solid");
}

template <typename VT>
bool STL<VT>::ReadASCII(std::ifstream& in) {
    std::size_t num_facets{0};
    std::string dummy;

    while (getline(in, dummy))
        ++num_facets;

    if (num_facets < 3) {
        std::cerr << "Cannot find any facets in the file! Aborting!" << std::endl;
        return false;
    }

    num_facets = (num_facets - 2) / 7;
    if (num_facets == 0) {
        std::cerr << "File corrupted! Aborting!" << std::endl;
        return false;
    }

    facets.resize(num_facets);

    /* Reset STLfile */
    in.clear();  // clears EOF flag
    in.seekg(0, std::ios::beg);

    /* Process and store data from STL file */
    getline(in, dummy);  // reads "solid <NAME>"

    typename Facet::PointType v[3];

    try {
        for (std::size_t i = 0; i < num_facets; ++i) {
            in >> dummy >> dummy >> dummy >> dummy >> dummy;  // reads "facet normal nx ny nz"
            in >> dummy >> dummy;                             // reads "outer loop"

            for (std::size_t j = 0; j < 3; ++j) {
                in >> dummy;  // reads "vertex"
                in >> v[j].x() >> v[j].y() >> v[j].z();
            }

            in >> dummy >> dummy;  // reads "endloop" and "endfacet"

            facets[i] = Facet(v[0], v[1], v[2]);
        }
    } catch (const std::exception& ex) {
        std::cerr << "Following exception was triggered while reading the STL file:\n"
                  << ex.what() << std::endl;
        return false;
    }

    return true;
}

template <typename VT>
bool STL<VT>::ReadBinary(std::ifstream& in) {
    using INT32 = int32_t;
    using REAL32 = float;

    auto ReadReal = [](std::ifstream& in) {
        REAL32 buf;
        in.read(reinterpret_cast<char*>(&buf), sizeof(REAL32));
        return buf;
    };
    /* Read header information (80 bytes) */
    char header_info[80];
    in.read(header_info, 80);

    /* Read total number of facets */
    char num_facets[sizeof(INT32)];
    in.read(num_facets, sizeof(INT32));
    INT32* num_facets_real = reinterpret_cast<INT32*>(num_facets);
    facets.resize(*num_facets_real);

    /* Process and store all facets in the file */
    for (std::size_t i = 0; i < facets.size(); ++i) {
        /* Create vector instance to temporarily store data */
        try {
            typename Facet::PointType temp[3];

            /* Read normal vector components and discard*/
            for (int8_t j = 0; j < 3; ++j)
                temp[0].data()[j] = ReadReal(in);

            /* Read vertex points */
            for (int8_t j = 0; j < 3; ++j) {
                for (int8_t k = 0; k < 3; ++k)
                    temp[j].data()[k] = ReadReal(in);

                facets[i] = Facet(temp[0], temp[1], temp[2]);
            }

            /* Read facet attribute (16 bits or 2 bytes) */
            char attribute[2];
            in.read(attribute, 2);
        } catch (std::exception& ex) {
            std::cerr << "An error occured while reading facet " << i
                      << ": " << ex.what() << std::endl;
            break;
        }
    }

    return true;
}

template <typename VT>
bool STL<VT>::WriteToASCII(std::ostream& ofs) {
    ofs << "object\n";

    for (auto& f : facets) {
        ofs << "facet normal "
            << f.GetNormal().x() << ' '
            << f.GetNormal().y() << ' '
            << f.GetNormal().z() << '\n';
        ofs << "outer loop\n";
        for (uint8_t i{0}; i < 3; i++)
            ofs << "vertex " << f.GetPoints()[i].x()
                << ' ' << f.GetPoints()[i].y()
                << ' ' << f.GetPoints()[i].z() << '\n';

        ofs << "endloop\n";
        ofs << "endfacet\n";
    }
    ofs << "endsolid";

    ofs.flush();

    return true;
}

template <typename VT>
bool STL<VT>::WriteToBinary(std::ostream& ofs) {
    using REAL32 = float;
    static_assert(sizeof(REAL32) == 4, "type float is not a 32 bit number!");
    uint32_t n_facets{static_cast<uint32_t>(facets.size())};
    char header[80];
    memset(this->header, 0, 80);
    ofs.write(header, 80);

    ofs.write(reinterpret_cast<char*>(&n_facets), 4);

    const uint16_t add_info{0};

    for (const auto& facet : facets) {
        REAL32 buf[3];  // small buffer for the cast to the 32 bit representation
        for (uint8_t i{0}; i < 3; ++i)
            buf[i] = static_cast<REAL32>(facet.GetNormal().data[i]);

        ofs.write(reinterpret_cast<char*>(buf), 3 * sizeof(REAL32));

        for (uint8_t i{0}; i < 3; ++i) {
            for (uint8_t j{0}; j < 3; ++j)
                buf[j] = static_cast<REAL32>(facet.GetPoints()[i].data()[j]);
            ofs.write(reinterpret_cast<char*>(buf), 3 * sizeof(REAL32));
        }

        ofs.write(reinterpret_cast<const char*>(&add_info), 2);
    }
    ofs.flush();

    return true;
}

template <typename VT>
const typename STL<VT>::TriangleContainer& STL<VT>::GetFacets() {
    return facets;
}

template <typename VT>
const typename STL<VT>::TriangleContainer& STL<VT>::GetFacets() const {
    return facets;
}

}  // end namespace dare::ff
