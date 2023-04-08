#include <cstdlib>
#include <cstdio>
#include <iostream>

#include <libintl.h>
#define _(String) gettext (String)
#define gettext_noop(String) String
#define N_(String) gettext_noop (String)

int main ()
{
 setlocale (LC_ALL, "");
 bindtextdomain ("hello", "/usr/share/locale/");
 textdomain ("hello");
  
  
 printf(_("Ceci est un test. %d %g\n"), 8, 3.14159) ; 
 std::cout << _("Un autre test ") << 3.145 << "\n" ; 
 
}
