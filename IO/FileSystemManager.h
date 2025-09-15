/*
 * MIT License
 *
 * Copyright (c) 2025 David Rieder
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

namespace dare {

namespace details {

/*! \struct FileSystemManager_helper
 * @brief a helper for the file system manager
 * used via the pimpl idiom to minimize the footprint
 */
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
    Path base_path_output;  //!< base path for output files
    Path base_path_input;   //!< base path for input files
    bool clear_contents;    //!< flag to clear contents of output directory
    bool overwrite_files;   //!< flag to overwrite existing files
    bool check_with_user;   //!< flag to check with user before deleting or overwriting files
};
}  // end namespace details

/*! \struct FileSystemManager
 * @brief a manager for handling file system operations
 *
 * This class provides a simple interface for managing file system operations.
 * It can be used to create directories, manage paths, and handle file overwriting
 * and clearing of contents. It is designed to work in parallel computing environments,
 * ensuring that operations are performed consistently across all processes.
 */
class FileSystemManager {
public:
    /*!
     * @brief constructor
     * @param ex_man pointer to execution manager
     * @param base_path basic path for input and output files
     * @param clear_contents flag if contents of output directory should be cleared
     * @param overwrite_files flag if files should be overwritten
     * @param check_with_user flag if user should be asked before deleting or overwriting files
     */
    explicit FileSystemManager(dare::ExecutionManager* ex_man,
                               const Path& base_path = Path(""),
                               bool clear_contents = false,
                               bool overwrite_files = false,
                               bool check_with_user = true);
    /*!
     * @brief destructor
     */
    ~FileSystemManager();

    /*!
     * @brief copy constructor
     * @param other instance to copy from
     */
    FileSystemManager(const FileSystemManager& other);    // NOLINT

    /*!
     * @brief copy assignment operator
     * @param other instance to copy from 
     * @return reference to this
     */
    FileSystemManager& operator=(const FileSystemManager& other); // NOLINT

    /*!
     * @brief set the output path
     */
    void SetOutputPath(const Path& path);

    /*!
     * @brief set the input path
     */
    void SetInputPath(const Path& path);

    /*!
     * @brief set the overwrite flag
     * @param overwrite flag for overwriting
     */
    void OverwriteFiles(bool overwrite);

    /*!
     * @brief set flag for clearing the contents of the output directory
     * @param clear flag to set
     */
    void ClearContents(bool clear);

    /*!
     * @brief set the flag for checking with the user before deleting
     * @param check flag to set
     */
    void CheckWithUser(bool check);

    /*!
     * @brief provides the output path
     * @return internal output path
     */
    const Path& GetOutputPath() const;

    /*!
     * @brief provides the input path
     * @return internal input path
     */
    const Path& GetInputPath() const;

    /*!
     * @brief creates a directory at the specified path
     * @param path path in form of a string
     * @param ignore_user_check flag if use should be asked before deleting contents
     * @return true if successful
     * \note This function is executed on each calling process, if only one common directory should
     * be created, use CreateCommonDirectory()
     */
    bool CreateDirectory(const std::string& path, bool ignore_user_check = false) const;

    /*!
     * @brief creates a directory at the specified path
     * @aram path path in form of a string
     * @return true if successful
     * \note only the root process creates the directory, all others wait for it to be created
     */
    bool CreateCommonDirectory(const std::string& path) const;

    /*!
     * @brief creates directory at the specified path
     * @param path path variable
     * @param ignore_user_check flag if use should be asked before deleting contents
     * @return true if successful
     */
    bool CreateDirectory(const Path& path, bool ignore_user_check = false) const;

    /*!
     * @brief creates a directory at the specified path
     * @param path path variable
     * @return true if successful
     * \note only the root process creates the directory, all others wait for it to be created
     */
    bool CreateCommonDirectory(const Path& path) const;

private:
    std::shared_ptr<details::FileSystemManager_helper> settings;    //!< pimpl idiom for settings
    dare::ExecutionManager* ex_man;     //!< reference to the execution manager
};

}  // end namespace dare

#endif  // IO_FILESYSTEMMANAGER_H_
