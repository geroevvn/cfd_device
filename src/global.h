/*
 * global.h
 *
 *  Created on: Nov 5, 2019
 *      Author: v1
 */

#ifndef GLOBAL_H_
#define GLOBAL_H_

#include "log4cpp/Category.hh"
#include "log4cpp/Appender.hh"
#include "log4cpp/FileAppender.hh"
#include "log4cpp/Layout.hh"
#include "log4cpp/BasicLayout.hh"
#include "log4cpp/SimpleLayout.hh"
#include "log4cpp/PatternLayout.hh"
#include "log4cpp/Priority.hh"


class Logger {

public:
	static Logger* Instance();
    void open_log_file(const char* logFile);
    log4cpp::Category* logging();
    bool close_log_file();

    void EXIT(int err);

private:
	Logger(){ root = 0; layout = 0; appender = 0; };
	Logger(Logger const&){};
	Logger& operator=(Logger const&){};

	log4cpp::Category* root;
	log4cpp::Layout* layout;
	log4cpp::Appender* appender;

	static Logger* instance;
};


#endif /* GLOBAL_H_ */
