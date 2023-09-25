all:
	g++ -std=c++20 -Wall -Wextra -g -fsanitize=address,undefined -D_GLIBCXX_DEBUG -Wl,--copy-dt-needed-entries -lgmpxx -o main main.cpp

fast:
	g++ -std=c++20 -Wall -Wextra -O3 -march=native -DNDEBUG -Wl,--copy-dt-needed-entries -lgmpxx -o main main.cpp