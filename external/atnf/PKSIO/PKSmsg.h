//#---------------------------------------------------------------------------
//# PKSmsg.h: Message handling for the PKSIO classes.
//#---------------------------------------------------------------------------
//# Copyright (C) 2008
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU Library General Public License as published by
//# the Free Software Foundation; either version 2 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
//# License for more details.
//#
//# You should have received a copy of the GNU Library General Public License
//# along with this library; if not, write to the Free Software Foundation,
//# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//# $Id: PKSmsg.h,v 1.2 2008-11-17 06:42:36 cal103 Exp $
//#---------------------------------------------------------------------------
//# Original: 2008/09/18, Mark Calabretta, ATNF
//#---------------------------------------------------------------------------

#ifndef ATNF_PKSMSG_H
#define ATNF_PKSMSG_H

#include <casa/stdio.h>

using namespace std;

// <summary>
// Message handling for the PKSIO classes.
// </summary>

class PKSmsg
{
  public:
    // Constructor.
    PKSmsg();

    // Destructor.
    virtual ~PKSmsg();

    // Set message disposition.  If fd is non-zero messages will be written
    // to that file descriptor, else stored for retrieval by getMsg().
    virtual int setMsg(
        FILE *fd = 0x0);

    // Log a message.
    virtual void logMsg(const char *msg = 0x0);

    // Get a message string, or 0x0 if there is none.  The null-terminated
    // message string may contain embedded newline characters and will have
    // a trailing newline.
    const char *getMsg();

    // Get the next group of messages by type: ERROR, WARNING, or otherwise.
    // The null-terminated message string may contain embedded newline
    // characters but will NOT have a trailing newline.  Call this repeatedly
    // to unwind the message stack (otherwise messages may be lost).
    enum msgType {NORMAL, WARNING, ERROR};
    const char *getMsg(msgType &type);

    // Clear the message buffer.
    void clearMsg(void);

  protected:
    // Initialize messaging.
    void initMsg();

  private:
    // For messaging.
    char *cMsgBuff, *cMsgIdx;
    int   cMsgLen, cNMsg;
    FILE *cMsgFD;
};

#endif
