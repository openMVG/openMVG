/**
This file enables to send data to a socket (or shared memory) using
the nanomsg library.
See here for more details:
https://github.com/nanomsg/nanomsg
**/

#ifdef USE_NANOMSG
#ifndef NANOMSG
#define NANOMSG

#include <assert.h>
// see http://stackoverflow.com/questions/21768542/libc-h-no-such-file-or-directory-when-compiling-nanomsg-pipeline-sample
#include <unistd.h>
#include <string>
#include <pthread.h>

#include <stdio.h>
#include <stdlib.h>
#include <nanomsg/nn.h>
#include <nanomsg/pair.h>
#include <iostream>
#include <exception>


class NanoMsgSender
{
public:
	/**
	@brief
	Binds a socket to the given URL. You will be able to use another client to connect to this.
	@param url: The url to bind to.
	@param time_out: Time to give up sending data. in milliseconds.
	**/
	NanoMsgSender(std::string url, int time_out=2):throw_timeout_errors(false)
	{
		this->sock=nn_socket (AF_SP, NN_PAIR);
		assert (this->sock >= 0);
  		assert (nn_bind (this->sock, url.c_str()) >= 0);
  		assert (nn_setsockopt (this->sock, NN_SOL_SOCKET, NN_SNDTIMEO, &time_out, sizeof (time_out)) >= 0);


	}
	/**
	@brief 
	If enables, exceptions are thrown when the timeout is reached
	**/
	void setThrowTimeOutErrors(bool do_throw)
	{
		throw_timeout_errors=do_throw;
	}

	bool sendMsg(std::string msg)
	{
		int sz_n = msg.size()+1;// '\0' too
		int rc=nn_send (this->sock, msg.c_str(), sz_n, 0);
		if (rc<0)
		{
		if (nn_errno()==EAGAIN)
		  if (this->throw_timeout_errors)
		  {
		  	throw std::runtime_error("timeout in NanoMsgSender.sendMsg");
		  } 
		}
	}

protected:
	int sock;
	bool throw_timeout_errors;

};



//nanomsg
#endif 
//if use_nanomsg
#endif 