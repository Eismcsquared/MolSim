#include "Logger.h"
#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/fmt/ostr.h>

std::shared_ptr<spdlog::logger> test_logger = spdlog::basic_logger_mt("simple_test_logger", "logs/test.txt");
