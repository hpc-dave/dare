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

#include "FileSystemManager.h"
namespace dare::io {

namespace details {

FileSystemManager_helper::FileSystemManager_helper(const Path& path_out,
                                                   const Path& path_in,
                                                   bool clear,
                                                   bool owrite,
                                                   bool check_w_user)
                                                   : base_path_output(path_out),
                                                     base_path_input(path_in),
                                                     clear_contents(clear),
                                                     overwrite_files(owrite),
                                                     check_with_user(check_w_user) {
}

FileSystemManager_helper::FileSystemManager_helper()
    : FileSystemManager_helper(Path(""), Path(""), false, false, true) {
}

FileSystemManager_helper::~FileSystemManager_helper() {}

const Path& FileSystemManager_helper::GetOutputPath() const {
    return base_path_output;
}
const Path& FileSystemManager_helper::GetInputPath() const {
    return base_path_input;
}
bool FileSystemManager_helper::ClearContents() const {
    return clear_contents;
}
bool FileSystemManager_helper::OverwriteFiles() const {
    return overwrite_files;
}
bool FileSystemManager_helper::CheckWithUser() const {
    return check_with_user;
}

void FileSystemManager_helper::SetOutputPath(const Path& path) {
    base_path_output = path;
}
void FileSystemManager_helper::SetInputPath(const Path& path) {
    base_path_input = path;
}
void FileSystemManager_helper::OverwriteFiles(bool overwrite) {
    overwrite_files = overwrite;
}
void FileSystemManager_helper::ClearContents(bool clear) {
    clear_contents = clear;
}
void FileSystemManager_helper::CheckWithUser(bool check) {
    check_with_user = check;
}

}  // end namespace details

FileSystemManager::FileSystemManager(dare::mpi::ExecutionManager* ex,
                                     const Path& base_path,
                                     bool clear_contents,
                                     bool overwrite_files,
                                     bool check_with_user)
                                     : settings(std::make_shared<details::FileSystemManager_helper>(
                                        base_path, base_path, clear_contents, overwrite_files, check_with_user)),
                                        ex_man(ex) {
    if (!ex) {
        std::cerr << "Execution Manager provided to FileSystemManager is a nullptr!" << std::endl;
    }
}
FileSystemManager::~FileSystemManager() {
}
FileSystemManager::FileSystemManager(FileSystemManager& other)  // NOLINT
    : settings(other.settings), ex_man(other.ex_man) {
}
FileSystemManager& FileSystemManager::operator=(FileSystemManager& other) {  // NOLINT
    if (&other == this)
        return *this;
    settings = other.settings;
    ex_man = other.ex_man;
    return *this;
}

void FileSystemManager::SetOutputPath(const Path& path) {
    settings->SetOutputPath(path);
}
void FileSystemManager::SetInputPath(const Path& path) {
    settings->SetInputPath(path);
}
void FileSystemManager::OverwriteFiles(bool overwrite) {
    settings->OverwriteFiles(overwrite);
}
void FileSystemManager::CheckWithUser(bool check) {
    settings->CheckWithUser(check);
}
void FileSystemManager::ClearContents(bool clear) {
    settings->ClearContents(clear);
}
const Path& FileSystemManager::GetOutputPath() const {
    return settings->GetOutputPath();
}
const Path& FileSystemManager::GetInputPath() const {
    return settings->GetInputPath();
}

bool FileSystemManager::CreateDirectory(const Path& path, bool ignore_user_check) const {
    if (details::Exists(path)) {
        if (details::IsDirectory(path)) {
            if (!details::IsDirectoryEmpty(path) && settings->ClearContents()) {
                if (path.empty()) {
                    ERROR << "Filepath is empty, I will NOT remove all contents!" << ERROR_CLOSE;
                    return false;
                } else {
                    if (settings->CheckWithUser() && !ignore_user_check) {
                        std::string answer;
                        std::cout << "Warning! The directory '" << path.native() << "' is not empty!\n";
                        for (int tries{0}; tries < 10; tries++) {
                            std::cout << "Do you really want to remove everything?[Y/n]" << std::endl;
                            std::cin >> answer;
                            if (boost::iequals(answer, "n")) {
                                std::cout << "Fine, I will not remove your precious items....\n"
                                          << "Check your settings better next time!" << std::endl;
                                return true;
                            }
                            if (boost::iequals(answer, "Y"))
                                break;
                            std::cout << "I did not understand : " << answer << "  -- try again!" << std::endl;
                        }
                        if (!boost::iequals(answer, "Y")) {
                            std::cout << "I tried multiple times, but only received gibberish.\n"
                                      << "Terminating execution now!" << std::endl;
                            ex_man->Terminate(__func__, "input error");
                        }
                    }
                    return details::RemoveAllEntries(path);
                }
            }
        }
    } else {
        if (!details::CreateDirectoryRecursive(path)) {
            ERROR << "Could not create following path: " << path.native() << ERROR_CLOSE;
            return false;
        }
    }
    return true;
}
bool FileSystemManager::CreateDirectory(const std::string& str_path, bool ignore_user_check) const {
    Path path(str_path);
    return CreateDirectory(path, ignore_user_check);
}

bool FileSystemManager::CreateCommonDirectory(const Path& path) const {
    if (ex_man->AmIRoot()) {
        bool ignore_user_check{false};
        if (!CreateDirectory(path, ignore_user_check)) {
            ex_man->Terminate(__func__, "filesystem error");
        }
    }
    return true;
}
bool FileSystemManager::CreateCommonDirectory(const std::string& str_path) const {
    Path path(str_path);
    return CreateCommonDirectory(path);
}
}  // end namespace dare::io
