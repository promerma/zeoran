
CXXINCLUDES=/usr/include/eigen3

zeoran: zeoran.cpp
	g++ zeoran.cpp -o zeoran -I $(CXXINCLUDES)

install:
	cp zeoran /usr/local/bin/
	cp -r zeoran_data /usr/local/share/zeoran

uninstall:
	rm /usr/local/bin/zeoran
	rm -rf /usr/local/share/zeoran

clean:
	rm -f zeoran
