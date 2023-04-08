#include "Headers/Ssh.h"

Ssh::Ssh()
{
strcpy(hosts[0], "guillard@192.168.1.1") ;
}
//----------------------------------------

Ssh::~Ssh()
{
sftp_free(sftp) ;
ssh_disconnect(session_ssh);
ssh_free(session_ssh);
}

//---------------------------------------------
void Ssh::session_init ()
{
int rc, state, i; char * password ;
session_ssh = ssh_new();
ssh_options_set(session_ssh, SSH_OPTIONS_HOST, hosts[0]);
rc = ssh_connect(session_ssh);
if (rc != SSH_OK)
{
  printf("Error connecting to hosts[0]: %s\n", ssh_get_error(session_ssh));
  return ;
}

state = ssh_is_server_known(session_ssh);

if (state != SSH_SERVER_KNOWN_OK)
   {printf("Erreur : serveur inconnu") ; return ; }

password = getpass("Mot de passe de la session ssh : ");
rc = ssh_userauth_password(session_ssh, NULL, password);
i=0 ; while(password[i]!=0) {password[i]=0 ; i++ ; }
if (rc != SSH_AUTH_SUCCESS)
  {
  fprintf(stderr, "Error authenticating with password: %s\n", ssh_get_error(session_ssh));
  ssh_disconnect(session_ssh);
  ssh_free(session_ssh);
  return ;
  }
// On est connecté normalement ...
sftp = sftp_new(session_ssh);
if (sftp == NULL) {fprintf(stderr, "Error allocating SFTP session: %s\n", ssh_get_error(session_ssh)); return ; }

rc = sftp_init(sftp);
if (rc != SSH_OK) {fprintf(stderr, "Error initializing SFTP session: %d.\n", sftp_get_error(sftp)); sftp_free(sftp); return ; }
printf("AllOK in SSH") ;
}

//------------------------------------------
