OS = MacOS
VERSION = Natrium
APP_NAME = DENIS-$(VERSION)-$(OS)
CC = g++


all: app.cpp
	$(CC) -O2 -Imodels tinyxml2.cpp app.cpp -o ../bin/$(APP_NAME)

clean:
	rm ../bin/$(APP_NAME)

