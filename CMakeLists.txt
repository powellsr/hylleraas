cmake_minimum_required(VERSION 3.17)
project(hylleraas)

set(CMAKE_CXX_STANDARD 17)

include_directories(. /usr/local/include)

add_executable(hylleraas
        hylleraas.cpp)

