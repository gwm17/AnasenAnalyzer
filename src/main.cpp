#include "Analyzer/Application.h"

#include "fmt/core.h"
#include "fmt/format.h"
#include "fmt/std.h"

int main(int argc, char** argv)
{
    if(argc != 2)
    {
        fmt::print("To run AnasenAnalyzer use the following CLI:\n");
        fmt::print("./bin/AnasenAnalyzer <config_file>\n");
        fmt::print("Replace <config_file> with the path to your configuration file.\n");
        return 1;
    }

    fmt::print("-------------Anasen Analysis Package-------------\n");
    AnasenAnalyzer::Application app;
    app.LoadConfig(argv[1]);
    if (!app.IsValid())
    {
        fmt::print("Configuration from {0} was not parsed correctly, shutting down application.\n", argv[1]);
        return 1;
    }
    app.Run();
    fmt::print("-------------------------------------------------\n");

    return 0;
}