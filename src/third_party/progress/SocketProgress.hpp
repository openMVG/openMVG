
/**
 * @file progress.h
 * @brief Simple progress bar for console application
 * @author Pierre MOULON
 *
 * Copyright (c) 2011, 2012, 2013 Pierre MOULON
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SOCKETPROGRESS
#define SOCKETPROGRESS

#include <iostream>
#include <string>
#include <algorithm>
#include <openMVG/system/nanomsg.hpp>
#include "progress.hpp"

using namespace std;

/**
* @brief Socket_Progress_send sends the incremental progress updates direclty
to a socket. The class uses NanoMsg system, which is build for distributed message
passing in a heterogenous cluster environment.
* Usage :
* Socket_Progress_send my_progress_bar( url, fileList.size() );
* for (list<string>::const_iterator it = fileList.begin(); it != fileList.end(); ++it, ++my_progress_bar)
* the C_Progress_display::operator++ take in charge the display of the progression (***)
* {
*   const string & filename = *it;
*   ... // do something
* }
* The sent result will be like the following :
* "{progress: {"title": "OpenMvgProcessName", "count":"11", "total": "342" }}";
*
**/

class Socket_Progress_send : public C_Progress
{
  public:
    /** @brief Standard constructor
      * @param expected_count The url and the number of step of the process
      * @param msg_title: The title of the message.
      * @param socket_url, e.g. ipc:///tmp/openmvg.ipc for shared memory (local) messages or tcp://<IP>:<PORT> for distributed messages.
      * @param ulExpected_count maximum number of counts to reach.
      **/
    explicit Socket_Progress_send ( std::string msg_title, std::string socket_url, unsigned long ulExpected_count):url(socket_url),
      title(msg_title)
    { 
      restart ( ulExpected_count ); 
    }

    /** @brief Initializer of the C_Progress_display class
     * @param expected_count The number of step of the process
     **/
    void restart ( unsigned long ulExpected_count )
    {
      C_Progress::restart ( ulExpected_count ); //-- Initialize the base class

      nano_msg_sender=std::shared_ptr<NanoMsgSender>(new NanoMsgSender(url));
    } // restart

  private:
    /** @brief Sends the update directly to the socket.
    **/
    void inc_tic()
    {
      
        std::string msg="{\"progress\": {\"title\": \""+this->title+"\", \"count\":"+std::to_string(_count)+", \"total\": "+std::to_string(_expected_count)+"}}";
        nano_msg_sender->sendMsg(msg);

    } // inc_tic

    std::shared_ptr<NanoMsgSender> nano_msg_sender;
    std::string url, title;
};

/** @} */

#endif
