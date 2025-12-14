emcc module.cpp -o module.js -lembind -s MODULARIZE=1 -sEXPORT_NAME="createWASMModule" -sEXPORTED_RUNTIME_METHODS='["ccall", "cwrap"]' -sALLOW_MEMORY_GROWTH -std=c++17 -O3 --closure 1 $@
