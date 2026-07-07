#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <unistd.h>

#include <cerrno>
#include <csignal>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

extern const unsigned char _binary_manytree_solver_start[];
extern const unsigned char _binary_manytree_solver_end[];
extern const unsigned char _binary_twotree_solver_start[];
extern const unsigned char _binary_twotree_solver_end[];

static pid_t child_pid = -1;

static void forward_signal(int sig) {
    if (child_pid > 0) kill(child_pid, sig);
}

static bool write_all(int fd, const unsigned char *data, size_t len) {
    while (len > 0) {
        ssize_t wrote = write(fd, data, len);
        if (wrote < 0) {
            if (errno == EINTR) continue;
            return false;
        }
        data += wrote;
        len -= static_cast<size_t>(wrote);
    }
    return true;
}

static std::string make_temp_path(const char *suffix) {
    std::ostringstream out;
    out << "/tmp/pace26_combined_" << static_cast<long long>(getpid()) << "_"
        << suffix;
    return out.str();
}

static bool write_file(const std::string &path, const unsigned char *data,
                       size_t len, mode_t mode) {
    int fd = open(path.c_str(), O_CREAT | O_TRUNC | O_WRONLY, mode);
    if (fd < 0) return false;
    bool ok = write_all(fd, data, len);
    if (close(fd) != 0) ok = false;
    chmod(path.c_str(), mode);
    return ok;
}

static bool write_text_file(const std::string &path, const std::string &text) {
    int fd = open(path.c_str(), O_CREAT | O_TRUNC | O_WRONLY, 0600);
    if (fd < 0) return false;
    bool ok = write_all(fd, reinterpret_cast<const unsigned char *>(text.data()),
                        text.size());
    if (close(fd) != 0) ok = false;
    return ok;
}

static bool parse_header(const std::string &input, int &t, int &n) {
    std::istringstream in(input);
    std::string line;
    while (std::getline(in, line)) {
        if (line.size() >= 2 && line[0] == '#' && line[1] == 'p') {
            std::string tag;
            std::istringstream ss(line);
            ss >> tag >> t >> n;
            return ss && t > 0 && n > 0;
        }
    }
    return false;
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    std::ostringstream buffer;
    buffer << std::cin.rdbuf();
    std::string input = buffer.str();
    if (input.empty()) return 1;

    int t = 0;
    int n = 0;
    if (!parse_header(input, t, n)) return 1;

    if (t == 2 && n > 500) {
        while (true) pause();
    }

    const unsigned char *solver_start = nullptr;
    const unsigned char *solver_end = nullptr;
    const char *solver_name = nullptr;
    if (t == 2) {
        solver_start = _binary_twotree_solver_start;
        solver_end = _binary_twotree_solver_end;
        solver_name = "twotree_solver";
    } else {
        solver_start = _binary_manytree_solver_start;
        solver_end = _binary_manytree_solver_end;
        solver_name = "manytree_solver";
    }

    std::string solver_path = make_temp_path(solver_name);
    std::string input_path = make_temp_path("input.nw");
    if (!write_file(solver_path, solver_start,
                    static_cast<size_t>(solver_end - solver_start), 0700)) {
        return 1;
    }
    if (!write_text_file(input_path, input)) {
        unlink(solver_path.c_str());
        return 1;
    }

    signal(SIGTERM, forward_signal);
    signal(SIGINT, forward_signal);

    child_pid = fork();
    if (child_pid < 0) {
        unlink(solver_path.c_str());
        unlink(input_path.c_str());
        return 1;
    }
    if (child_pid == 0) {
        int fd = open(input_path.c_str(), O_RDONLY);
        if (fd < 0) _exit(127);
        if (dup2(fd, STDIN_FILENO) < 0) _exit(127);
        close(fd);
        char *const argv[] = {const_cast<char *>(solver_path.c_str()), nullptr};
        execv(solver_path.c_str(), argv);
        _exit(127);
    }

    int status = 0;
    while (waitpid(child_pid, &status, 0) < 0) {
        if (errno != EINTR) {
            unlink(solver_path.c_str());
            unlink(input_path.c_str());
            return 1;
        }
    }
    child_pid = -1;
    unlink(solver_path.c_str());
    unlink(input_path.c_str());

    if (WIFEXITED(status)) return WEXITSTATUS(status);
    if (WIFSIGNALED(status)) {
        int sig = WTERMSIG(status);
        signal(sig, SIG_DFL);
        raise(sig);
        return 128 + sig;
    }
    return 1;
}
