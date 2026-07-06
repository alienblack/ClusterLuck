#include "util/signal_handler.h"

#include <csignal>
#include <string>
#include <unistd.h>

namespace {

std::string g_best_output;
bool g_publish_enabled = true;

void handler(int) {
    if (!g_best_output.empty()) {
        ssize_t unused = write(STDOUT_FILENO, g_best_output.data(), g_best_output.size());
        (void)unused;
    }
    _exit(0);
}

}  // namespace

void install_signal_handler() {
    std::signal(SIGTERM, handler);
    std::signal(SIGINT, handler);
}

void publish_best_output(std::string output) {
    if (!g_publish_enabled) return;
    g_best_output = std::move(output);
}

void set_publish_enabled(bool enabled) {
    g_publish_enabled = enabled;
}
