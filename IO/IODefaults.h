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

#ifndef IO_IODEFAULTS_H_
#define IO_IODEFAULTS_H_

#include <string>

namespace dare::io {

struct IODefaultOptions {
    std::string PROCESS_FOLDER_BASE_NAME = "";
    std::string OUTPUT_BASE_PATH = "";
    std::string INPUT_BASE_PATH = "";
} IODefaults;

std::string GetFolderName(int process) {
    return IODefaults.PROCESS_FOLDER_BASE_NAME + std::to_string(process);
}

std::string GetBaseOutputPath() {
    return IODefaults.OUTPUT_BASE_PATH;
}

}  // end namespace dare::io

#endif  // IO_IODEFAULTS_H_
