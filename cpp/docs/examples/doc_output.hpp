#pragma once
// Redirects std::cout into one text file per documented section, so that the
// guide pages can show the output of an example next to its code.

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

/**
 * @brief Break a long message into lines at word boundaries, for the output
 *        boxes of the guide pages.
 *
 * @param text  Message on one line.
 * @param width Largest line length.
 *
 * @return The message with line breaks.
 */
inline std::string wrapped(const std::string& text, std::size_t width = 90) {
    std::istringstream words(text);
    std::string word, line, out;
    while (words >> word) {
        if (!line.empty() && line.size() + 1 + word.size() > width) {
            out += line + '\n';
            line.clear();
        }
        line += (line.empty() ? "" : " ") + word;
    }
    return out + line;
}

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
