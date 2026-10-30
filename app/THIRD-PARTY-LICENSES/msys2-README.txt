MAFFT for Windows: bundled MSYS2 components
==========================================

The Windows build bundles MAFFT's official Windows package (mafft-win, MAFFT
7.526). MAFFT's Windows driver needs a small Unix environment, so that package
carries a partial MSYS2 tree under mafft/usr/bin/. These programs and libraries
are separate works with their own licences. They ship only in the Windows
build. The macOS and Linux MAFFT packages do not contain them.

Versions were read from the version strings inside each binary. "n/r" means the
binary does not record its version.

Component (MSYS2 package)   Version      Files                                  Licence
--------------------------  -----------  -------------------------------------  ------------------------------------------
bash                        4.4.0        bash.exe, sh.exe                       GPL-3.0-or-later
coreutils                   8.25         basename cat chmod comm cp cut date    GPL-3.0-or-later
                                         dirname echo env expr false fold head
                                         id join ln ls md5sum mkdir mktemp mv
                                         od paste printf pwd rm rmdir sleep
                                         sort split stty tail touch tr true
                                         uname uniq wc (.exe)
gawk                        4.1.4        awk.exe                                GPL-3.0-or-later
grep                        2.26         grep.exe                               GPL-3.0-or-later
sed                         4.2.2        sed.exe                                GPL-3.0-or-later
tar                         1.29         tar.exe                                GPL-3.0-or-later
findutils                   4.6.0        xargs.exe                              GPL-3.0-or-later
which                       2.21         which.exe                              GPL-3.0-or-later
gzip                        n/r          gzip.exe                               GPL-3.0-or-later
texinfo                     n/r          info.exe                               GPL-3.0-or-later
readline                    n/r          msys-readline6.dll                     GPL-3.0-or-later
msys2-runtime (Cygwin)      n/r          msys-2.0.dll                           LGPL-3.0-or-later, with the Cygwin linking
                                                                                exception (msys2-cygwin-runtime.txt)
msys2-runtime utilities     n/r          ps.exe                                 GPL-3.0-or-later
gcc-libs                    GCC 5.3.0    msys-gcc_s-seh-1.dll                   GPL-3.0-or-later with the GCC Runtime
                                                                                Library Exception 3.1
gmp                         6.1.0        msys-gmp-10.dll                        LGPL-3.0-or-later or GPL-2.0-or-later
mpfr                        3.1.4        msys-mpfr-4.dll                        LGPL-3.0-or-later
libiconv                    n/r          msys-iconv-2.dll                       LGPL-2.1-or-later
gettext (libintl)           n/r          msys-intl-8.dll                        LGPL-2.1-or-later
ncurses                     6.0          msys-ncursesw6.dll                     MIT/X11 style
pcre                        8.38         msys-pcre-1.dll                        BSD-3-Clause
zlib                        1.2.8        msys-z.dll                             zlib
file / libmagic             5.28         file.exe, msys-magic-1.dll             BSD-2-Clause
util-linux                  2.26.2       more.exe                               BSD (University of California, Berkeley)

Not MSYS2: usr/bin/mafft, usr/bin/mafft-homologs.rb and usr/lib/mafft/ are
MAFFT itself (BSD, mafft-BSD.txt). usr/share/misc/magic is the data file of
file / libmagic above.

Licence texts
-------------
GPL-3.0          ../LICENSE (the app's own licence, shipped beside this folder)
GPL-2.0          msys2-GPL-2.0.txt
LGPL-3.0         msys2-LGPL-3.0.txt
LGPL-2.1         msys2-LGPL-2.1.txt
GCC exception    msys2-GCC-RUNTIME-LIBRARY-EXCEPTION-3.1.txt
Cygwin runtime   msys2-cygwin-runtime.txt
ncurses          msys2-ncurses-MIT.txt
PCRE             msys2-pcre-BSD-3-Clause.txt
zlib             msys2-zlib.txt
file             msys2-file-BSD-2-Clause.txt
more (BSD)       msys2-util-linux-more-BSD.txt

Source code
-----------
These binaries are unmodified copies from the mafft-win package. Source for
every MSYS2 package is published by the MSYS2 project:
  https://github.com/msys2/MSYS2-packages   (build recipes)
  https://repo.msys2.org/msys/sources/      (source archives)
and by each upstream project (for example https://ftp.gnu.org/gnu/ for the GNU
tools). MAFFT's own source is at https://mafft.cbrc.jp/alignment/software/.
