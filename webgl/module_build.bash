emcc module.cpp -o module.js -lembind -s MODULARIZE=1 -sEXPORT_NAME="createModule" -sEXPORTED_RUNTIME_METHODS='["ccall", "cwrap"]' -sALLOW_MEMORY_GROWTH -std=c++11 -O3 --closure 1 $@
