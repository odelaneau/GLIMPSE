#ifndef __FS_HPP__
#define __FS_HPP__

#include <libgen.h> // Has dirname() / basename()
#include <string>
#include <cstring> // for strdup
#include <fstream>
#include <unistd.h>
/*
The functions dirname() and basename() break a null-terminated pathname string
into directory and filename components. In the usual case, dirname() returns the
string up to, but not including, the final '/', and basename() returns the
component following the final '/'. Trailing '/' characters are not counted as
part of the pathname.
*/

/// @todo move this to an object file
namespace
{
    class NamedFileStream {
    public:
        NamedFileStream(std::string filename) : stream(filename, stream.binary | stream.out | stream.trunc), filename(filename) {}
        std::fstream stream;
        std::string filename;
    };

    NamedFileStream get_temporary_file(int* file_desc = nullptr /** @todo flags */) {
        char *tmpname = strdup("/tmp/tmpfileXXXXXX");
        int fd = mkstemp(tmpname); /// @todo check return code
        std::string filename(tmpname);
        free(tmpname);
        if (fd) {
            NamedFileStream nfs(filename); // Opens as stream
            if (file_desc) {
                *file_desc = fd;
            } else {
                // Close the file descriptor if not used (file is open as stream)
                close(fd);
            }
            if (!nfs.stream.is_open()) {
                std::cerr << "Could not open temporary file stream to : " << filename << std::endl;
                throw "Could not get temporary file";
            }
            return nfs;
        } else {
            throw "Could not get temporary file";
        }

        return NamedFileStream(filename);
    }
} // namespace

#if __cplusplus >= 201703L
    #include <filesystem>
    namespace fs = std::filesystem;
    using fs::remove;
    // This can be used instead of basename() (more portable, if C++17 is available)
    // fs::path( "/foo/bar.txt" ).filename() => "bar.txt"
#else
    #include <stdio.h> // Has remove()
    #include <sys/stat.h>
    namespace fs {
        inline size_t file_size(const std::string& filename) {
            struct stat st;
            if (stat(filename.c_str(), &st) < 0) {
                std::cerr << "Size of file : " << filename << " could not be determined !" << std::endl;
                return 0;
            } else {
                return st.st_size;
            }
        }

        inline bool exists(const std::string& filename) {
            struct stat st;
            return stat(filename.c_str(), &st) == 0;
        }
    };
#endif

#endif /* __FS_HPP__ */