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

#include <iostream>
#include "Utilities/Vector.h"
#include "FileFormats/STL.h"


int main(int argc, char* argv[]){
    dare::utils::Vector<3, int> vec;
    dare::utils::Vector<3, double> vec_d;

    vec.i() = 1;
    vec.j() = 2;
    vec.k() = 3;
    for(auto v : vec){
        std::cout << v << std::endl;
    }
    vec_d.x() = 4;
    vec_d.y() = 2;
    vec_d.z() = 3;
    for(auto v : vec_d){
        std::cout << v << std::endl;
    }

    using P = dare::ff::STLfacet<double>::PointType;
    P v1(0., 0., 0.);
    P v2(1., 0., 0.);
    P v3(0., 1., 0.);
    dare::ff::STLfacet<double> facet(v1, v2, v3);
    std::cout << facet << std::endl;
    std::cout << "Hello world\n";
    return 0;
}
