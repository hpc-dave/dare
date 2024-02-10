/*
 * MIT License
 *
 * Copyright (c) 2024 David Rieder
 *
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

#include "FileSystemOperations.h"
namespace dare::io {
namespace details {

bool Exists(const Path& path) {
    return boost::filesystem::exists(path);
}
bool IsDirectory(const Path& path) {
    return boost::filesystem::is_directory(path);
}

bool IsDirectoryEmpty(const Path& path) {
    return boost::filesystem::is_empty(path);
}

bool CreateDirectory(const Path& path) {
    return boost::filesystem::create_directory(path);
}

bool RemoveAllEntries(const Path& path) {
    try {
        using Iterator = boost::filesystem::directory_iterator;
        for (Iterator end_it, it(path); it != end_it; it++) {
            boost::filesystem::remove_all(it->path());
        }
    } catch (...) {
        return false;
    }
    return true;
}

bool CreateDirectoryRecursive(const Path& path) {
    try {
        boost::filesystem::path path_create;
        auto it = path.begin();
        if (boost::iequals(it->native(), "/")) {
            path_create += *it;
            ++it;
        }
        for (; it != path.end(); it++) {
            path_create += *it;
            path_create += "/";
            if (!boost::filesystem::exists(path_create))
                boost::filesystem::create_directory(path_create);
        }
    } catch (...) {
        return false;
    }
    return true;
}
}  // end namespace details

}  // end namespace dare::io
