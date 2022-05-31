cmake_minimum_required(VERSION 2.8)
project (AVOIDREJECT)

file(GLOB avoidreject_source
	"*.hpp"
	"*.cpp"
)

add_executable(avoidreject ${avoidreject_source})

add_definitions(-std=c++11 -Ofast) 

install(TARGETS avoidreject DESTINATION bin)
