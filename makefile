functions: functions.cpp
	g++ functions.cpp -o zeoran

install:
	mv global.h /usr/include
	mv headers.h /usr/include
	mv libraries.h /usr/include
	rm functions.cpp
	cp -r zeoran_data /usr/include
	rm -r zeoran_data
