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

#ifndef IO_FILESYSTEMMANAGER_H_
#define IO_FILESYSTEMMANAGER_H_
#include <string>
#include <memory>
#include <boost/algorithm/string/predicate.hpp>
#include "FileSystemOperations.h"
#include "MPI/ExecutionManager.h"
#include "Utilities/Errors.h"

    namespace dare::io {

namespace details {
class FileSystemManager_helper {
public:
    FileSystemManager_helper(const Path& path_out,
                             const Path& path_in,
                             bool clear_contents,
                             bool overwrite,
                             bool check_with_user);
    FileSystemManager_helper();
    ~FileSystemManager_helper();
    FileSystemManager_helper(const FileSystemManager_helper&) = delete;
    FileSystemManager_helper& operator=(const FileSystemManager_helper&) = delete;

    const Path& GetOutputPath() const;
    const Path& GetInputPath() const;
    bool OverwriteFiles() const;
    bool ClearContents() const;
    bool CheckWithUser() const;

    void SetOutputPath(const Path& path);
    void SetInputPath(const Path& path);
    void ClearContents(bool clear);
    void OverwriteFiles(bool overwrite);
    void CheckWithUser(bool check);

private:
    Path base_path_output;
    Path base_path_input;
    bool clear_contents;
    bool overwrite_files;
    bool check_with_user;
};
}  // end namespace details

class FileSystemManager {
public:
    explicit FileSystemManager(dare::mpi::ExecutionManager* ex_man,
                               const Path& base_path = Path(""),
                               bool clear_contents = false,
                               bool overwrite_files = false,
                               bool check_with_user = true);
    ~FileSystemManager();
    FileSystemManager(FileSystemManager& other);    // NOLINT
    FileSystemManager& operator=(FileSystemManager& other); // NOLINT

    void SetOutputPath(const Path& path);
    void SetInputPath(const Path& path);
    void OverwriteFiles(bool overwrite);
    void ClearContents(bool clear);
    void CheckWithUser(bool check);
    const Path& GetOutputPath() const;
    const Path& GetInputPath() const;

    bool CreateDirectory(const std::string& path, bool ignore_user_check = false) const;
    bool CreateCommonDirectory(const std::string& path) const;
    bool CreateDirectory(const Path& path, bool ignore_user_check = false) const;
    bool CreateCommonDirectory(const Path& path) const;

private:
    std::shared_ptr<details::FileSystemManager_helper> settings;
    dare::mpi::ExecutionManager* ex_man;
};

}  // end namespace dare::io

#endif  // IO_FILESYSTEMMANAGER_H_
