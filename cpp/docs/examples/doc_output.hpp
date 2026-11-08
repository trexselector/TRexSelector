#pragma once
// Redirects std::cout into one text file per documented section, so that the
// guide pages can show the output of an example next to its code.

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

/**
 * @brief Section-wise capture of std::cout for the documentation examples.
 *
 * @details begin(name) closes the current capture file and sends std::cout into
 *          <dir>/<example>.<name>.txt. The destructor restores std::cout.
 */
class DocSections {
public:
    DocSections(std::string example, std::string dir)
        : example_(std::move(example)), dir_(std::move(dir)), original_(std::cout.rdbuf()) {
        std::filesystem::create_directories(dir_);
    }

    ~DocSections() {
        std::cout.flush();
        std::cout.rdbuf(original_);
    }

    DocSections(const DocSections&) = delete;
    DocSections& operator=(const DocSections&) = delete;

    /** @brief Start capturing the section @p name. */
    void begin(const std::string& name) {
        std::cout.flush();
        file_.close();
        file_.open(dir_ + "/" + example_ + "." + name + ".txt");
        std::cout.rdbuf(file_.rdbuf());
    }

private:
    std::string example_;
    std::string dir_;
    std::streambuf* original_;
    std::ofstream file_;
};
