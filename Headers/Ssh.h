#ifndef MASTERH
#include "MasterHeader.h"

#else
#include <libssh/libssh.h>
#include <libssh/sftp.h>

class Ssh {
public :
	Ssh() ;
	~Ssh() ;
	void session_init() ;

	char hosts[2][500] ;
	ssh_session session_ssh  ;
	sftp_session sftp;


} ;

#endif
