#include "chunker_header.h"          // chunk
#include "spliter_header.h"          // split_ref
#include "caller/caller_header.h"    // phase
#include <ligater/ligater_header.h>  // ligate
#include "checker_header.h"          // concordance

#include <boost/program_options.hpp>
#include <iostream>
#include <string>
#include <vector>
#include <functional>
#include <unordered_map>

namespace po = boost::program_options;

int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cerr << "Usage: Glimpse2 <command> [options]\n"
                  << "Commands: chunk, split_ref, phase, ligate, concordance\n";
        return 1;
    }

    std::string cmd = argv[1];

    // Map subcommand → handler
    const std::unordered_map<std::string, std::function<void(const std::vector<std::string>&)>> commands = {
        { "chunk",       [](auto args){ chunker().chunk(args); } },
        { "split_ref",   [](auto args){ spliter().phase(args); } },
        { "phase",       [](auto args){ caller().phase(args); } },
        { "ligate",      [](auto args){ ligater().ligate(args); } },
        { "concordance", [](auto args){ checker().check(args); } }
    };

    auto it = commands.find(cmd);
    if (it == commands.end()) {
        std::cerr << "Unknown command: " << cmd << "\n"
                  << "Available commands: chunk, split_ref, phase, ligate, concordance\n";
        return 1;
    }

    // Build vector from argv[2] to argv[argc-1]
    std::vector<std::string> sub_args(argv + 2, argv + argc);

    // Call the selected subcommand
    it->second(sub_args);

    return 0;
}
