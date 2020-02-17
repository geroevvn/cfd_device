/*
 * global.cpp
 *
 *  Created on: Nov 5, 2019
 *      Author: v1
 */

#include "global.h"



Logger* Logger::instance = 0;


Logger* Logger::Instance()
{
	if( !instance )
	{
		instance = new Logger();
	}

	return instance;
}

void Logger::open_log_file(const char* fileName)
{
	if(appender != 0)
	{
		this->close_log_file();
	}

	appender = new log4cpp::FileAppender("default", fileName, false);
	layout = new log4cpp::BasicLayout();
	appender->setLayout(layout);

	root = &log4cpp::Category::getRoot();
	root->setPriority(log4cpp::Priority::INFO);
	root->addAppender(appender);
}

bool Logger::close_log_file()
{
	if(appender != 0)
	{
		appender->close();
		delete appender;
		delete layout;

		appender = 0;
		layout = 0;

		return true;
	}

	return false;
}


log4cpp::Category* Logger::logging()
{
	return root;
}

void Logger::EXIT(int err)
{
	root->error("Error : %d", err);
	this->close_log_file();
}

