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

#include "Cartesian.h"

namespace dare::Grid::details::Cartesian {
std::list<std::string> AllocationManager::reg;
bool AllocationManager::RegisterGrid(const std::string& gname) {
    // test if grid with same name was already allocated and register this one
    auto it = std::find(reg.begin(),
                        reg.end(),
                        gname);
    if (it != reg.end()) {
        ERROR << "The requested grid name '" << gname << "'is already registered, chose an alternative!" << ERROR_CLOSE;
        return false;
    }
    reg.push_back(gname);
    return true;
}

bool AllocationManager::DeregisterGrid(const std::string& gname) {
    auto it = std::find(reg.begin(),
                        reg.end(),
                        gname);
    if (it == reg.end()) {
        ERROR << "The grid " << gname << " could not be found during deallocation in the registry" << ERROR_CLOSE;
        return false;
    }
    reg.erase(it);
    return true;
}

}  // end namespace dare::Grid::details::Cartesian
